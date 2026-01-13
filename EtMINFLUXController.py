import os
import glob
import sys
import importlib
import enum
import warnings
import time
import math
import json

import specpy
import mouse

import tifffile as tiff
import scipy.ndimage as ndi
import numpy as np
import pyqtgraph as pg

from collections import deque
from datetime import datetime
from inspect import signature
from qtpy import QtCore
from PyQt5 import QtWidgets
from pynput.keyboard import Key, Controller

warnings.filterwarnings("ignore")


class EtMINFLUXController(QtCore.QObject):
    """ Controller part of the View-Controller-based etMINFLUX smart microscopy software. Interacts with an open Imspector software and controls the etMINFLUX experiments.
    In order to use this, you need to have Imspector running before launching this software. The interaction with Imspector uses the specpy package, a Python wrapper for 
    the Imspector API. If Imspector is not running, a mock Imspector connection will be created, which allows to test the GUI. """

    helpPlotDetectedCoordsSignal = QtCore.Signal(object, object, object)

    def __init__(self, widget,  *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._widget = widget
        
        print('Initializing etMINFLUX controller')
        
        # SYSTEM-SPECIFIC SETTINGS - MAKE SURE TO CHANGE THESE VARIABLES IN THE SETUP JSON TO YOUR SPECIFIC SYSTEM CONFIG
        self.setupInfo = self.loadSetupJson(filepath='etMINFLUX_setup.json')
        # default data dir
        self._dataDir = self.setupInfo.get('save_settings').get('save_directory')
        # list of available detectors for MFX imaging
        self.mfxDetectorList = self.setupInfo.get('hardware_settings').get('mfx_detectors_imspector')
        self.mfxDetectorListThreads = self.setupInfo.get('hardware_settings').get('mfx_detectors_imspector_threads')
        # default screen position of top MINFLUX dataset in Workspace widget in Imspector
        self._set_topmfxdataset_button_pos = self.setupInfo.get('gui_settings').get('top_mfx_dataset_pos')
        # if any activation lasers are present
        self._act_lasers_present = len(self.setupInfo.get('hardware_settings').get('mfx_act_lasers')) > 0
        # set default values of various parameters in the GUI from the setup JSON
        self._widget.setDefaultValues(self.setupInfo, self._act_lasers_present)
        
        self._widget.coordTransformWidget.setSaveFolderField(self._dataDir)

        # open imspector connection
        self._imspector = specpy.get_application()
        print('Imspector connected succesfully')

        # create an instance of keyboard emulator
        self.keyboard = Controller()

        # folders for analysis pipelines and transformations
        self.analysisDir = os.path.join('analysis_pipelines')
        if not os.path.exists(self.analysisDir):
            os.makedirs(self.analysisDir)
        sys.path.append(self.analysisDir)
        self.transformDir = os.path.join('transform_pipelines')
        if not os.path.exists(self.transformDir):
            os.makedirs(self.transformDir)
        sys.path.append(self.transformDir)
        # set lists of analysis pipelines and transformations in the widget
        self._widget.setAnalysisPipelines(self.analysisDir)
        self._widget.setTransformations(self.transformDir)

        # set list of available MFX sequences in the widget
        self._widget.setMfxSequenceList(self.setupInfo.get('acquisition_settings').get('minflux_seqs'), thread=0)
        self._widget.setMfxSequenceList(self.setupInfo.get('acquisition_settings').get('minflux_seqs'), thread=1)

        # set list of available lasers for MFX imaging, get this list manually from Imspector control software
        self._widget.setMinfluxExcLaserList(self.setupInfo.get('hardware_settings').get('mfx_exc_lasers'), thread=0)
        self._widget.setMinfluxExcLaserList(self.setupInfo.get('hardware_settings').get('mfx_exc_lasers'), thread=1)
        self._widget.setMinfluxDetectorList(self.setupInfo.get('hardware_settings').get('mfx_detectors_gui'), thread=0)
        self._widget.setMinfluxDetectorList(self.setupInfo.get('hardware_settings').get('mfx_detectors_gui'), thread=1)
        if self._act_lasers_present:
            self._widget.setMinfluxActLaserList(self.setupInfo.get('hardware_settings').get('mfx_act_lasers'))

        # create a helper controller for the coordinate transform pop-out widget
        self.__coordTransformHelper = EtCoordTransformHelper(self, self._widget.coordTransformWidget)
        self.__binaryMaskHelper = BinaryMaskHelper(self, self._widget.binaryMaskWidget)
        self.__analysisHelper = AnalysisImgHelper(self, self._widget.analysisHelpWidget)
        self.__multiMFXROIHelper = AnalysisImgHelper(self, self._widget.multiMFXROIHelpWidget)
        self.__eventViewHelper = EventWidgetHelper(self, self._widget.eventViewWidget)

        # Connect EtMINFLUXWidget button and check box signals
        self._widget.initiateButton.clicked.connect(self.initiate)
        self._widget.loadPipelineButton.clicked.connect(self.loadPipeline)
        self._widget.softResetButton.clicked.connect(self.softReset)
        self._widget.saveCurrentMeasButton.clicked.connect(self.saveMeasurement)
        self._widget.startMFXPhaseButton.clicked.connect(self.toggleMultipleDetectionMFXInitiation)

        self._widget.presetROISizeCheck.clicked.connect(self.togglePresetROISize)
        self._widget.presetMfxRecTimeCheck.clicked.connect(self.togglePresetRecTime)
        self._widget.autoSaveCheck.clicked.connect(self.togglePresetAutoSave)
        self._widget.autoDeleteMFXDatasetCheck.clicked.connect(self.togglePresetAutoDeleteMFX)
        self._widget.twoThreadsMFXCheck.clicked.connect(self.toggleTwoThreadMFX)

        self._widget.coordListWidget.delROIButton.clicked.connect(self.deleteROI)
        self._widget.followROIModeCheck.clicked.connect(self.toggleFollowingROI)
        self._widget.triggerAllROIsCheck.clicked.connect(self.toggleTriggerAllROIs)
        self._widget.confocalFramePauseCheck.clicked.connect(self.toggleConfocalFramePause)

        self._widget.coordTransformWidget.setSaveDirButton.clicked.connect(self.getSaveFolder)
        self._widget.openGuideButton.clicked.connect(self.openGuide)

        # create timer for fixed recording time syncing
        self.recTimeTimerThread = QtCore.QThread(self)
        self.recTimeTimer = QtCore.QTimer()
        self.recTimeTimer.moveToThread(self.recTimeTimerThread)
        self.recTimeTimer.setSingleShot(True)
        self.recTimeTimer.timeout.connect(self.initiate)
        self.recTimeTimerThread.started.connect(self.recTimeTimer.start)

        # create timer for confocal intermittent pause time
        self.confPauseTimerThread = QtCore.QThread(self)
        self.confPauseTimer = QtCore.QTimer()
        self.confPauseTimer.moveToThread(self.confPauseTimerThread)
        self.confPauseTimer.setSingleShot(True)
        self.confPauseTimer.timeout.connect(self.startConfocalScanning)
        self.confPauseTimerThread.started.connect(self.confPauseTimer.start)

        # create timers for confocal intermittent pause countdown timer showing
        self.confPauseGUICountdownTimerThread = QtCore.QThread(self)
        self.confPauseGUICountdownTimer = QtCore.QTimer()
        self.confPauseGUICountdownTimer.moveToThread(self.confPauseGUICountdownTimerThread)
        self.confPauseGUICountdownTimer.setSingleShot(False)
        self.confPauseGUICountdownTimer.timeout.connect(self.updateIntervalTimerDisp)
        self.confPauseGUICountdownTimerThread.started.connect(self.confPauseGUICountdownTimer.start)
        self.confIntervalDispTimer = QtCore.QElapsedTimer()

        # create timer for multiple event site detection until MINFLUX initiation in MultipleDetect ROI following mode
        self.multipleEventDetectionTimerThread = QtCore.QThread(self)
        self.multipleEventDetectionTimer = QtCore.QTimer()
        self.multipleEventDetectionTimer.moveToThread(self.multipleEventDetectionTimerThread)
        self.multipleEventDetectionTimer.setSingleShot(True)
        self.multipleEventDetectionTimer.timeout.connect(self.toggleMultipleDetectionMFXInitiation)
        self.multipleEventDetectionTimerThread.started.connect(self.multipleEventDetectionTimer.start)

        self.helpPlotDetectedCoordsSignal.connect(self.plotDetectedCoords)

        # initiate log for each detected event
        self.resetDetLog()
        # initiate pipeline parameter values
        self.resetPipelineParamVals()
        # initiate run parameters
        self.resetRunParams()
        # initiate other parameters and flags used during experiments
        self.initiateFlagsParams()
        # initiate plotting parameters
        self.initiatePlottingParams()

    def initiateFlagsParams(self):
        # initiate flags and params
        self.__running = False  # run flag
        self.__runningMFX = False  # run MFX flag
        self.__runMode = RunMode.Experiment  # run mode currently used
        self.__validating = False  # validation flag
        self.__busy = False  # running pipeline busy flag
        self.__presetROISize = True
        self.__presetRecTime = True
        self.__lineWiseAnalysis = False
        self.__run_all_aoi = False  # run all detected events/area flag
        self.__autoSaveMeas = False
        self.__autoDelMFX = False
        self.__followingROI = False
        self.__followingROIMode = ROIFollowMode.Single  # roi following mode currently selected
        self.__followingROIContinue = False
        self.__followingROIEnding = False
        self.__plotROI = False
        self.__confocalFramePause = False
        self.__count_conf_channels = 1
        self.__mfx_twothreaded = False
        self.__multipleDetectionMFXInitiation = False
        self.__followingROIMultipleDetectMINFLUXPhase = False
        self.__img_ana = None
        self.__confoffset = [0.0, 0.0]  # offset of confocal scan, in um
        self.__prev_event_coords_deque = deque(maxlen=100)
        self.__prevFrames = deque(maxlen=50)  # deque for previous fast frames
        self.__prevAnaFrames = deque(maxlen=50)  # deque for previous preprocessed analysis frames
        self.__roi_events = deque(maxlen=0)
        self.__roi_sizes = deque(maxlen=0)
        self.binary_mask = None  # binary mask of regions of interest, used by certain pipelines, leave None to consider the whole image
        self.binary_frames = 2  # number of frames to use for calculating binary mask 
        self.__init_frames = 0  # number of frames after initiating etMINFLUX before a trigger can occur, to allow laser power settling etc
        self.__validation_frames = 2  # number of fast frames to record after detecting an event in validation mode
        self.__params_exclude = ['img_ch1', 'img_ch2', 'img_ch3', 'prev_frames', 'binary_mask', 'exinfo', 'testmode', 'presetROIsize']  # excluded pipeline parameters when loading param fields
        self.__rec_time_deadtime = 5  # deadtime when starting MINFLUX recordings, in s - specifically for m2410 that does not lock everything during this deadtime
        self.imgs = []  # list of confocal images from different channels

    def getSaveFolder(self):
        self._dataDir = self._widget.coordTransformWidget.getSaveFolder()
        self._widget.coordTransformWidget.setSaveFolderField(self._dataDir)

    def openGuide(self):
        # load text from markdown file
        guidetext_url = QtCore.QUrl.fromLocalFile("guidetext.md")
        # show text in subwidget
        self._widget.guideWidget.setText(source=guidetext_url)
        # open subwidget
        self._widget.launchHelpWidget(self._widget.guideWidget, init=True)

    def initiatePlottingParams(self):
        self.__colors = deque(['g','g','g','g','g'])

    def popColor(self):
        col = self.__colors.popleft()
        self.__colors.append(col)
        return col

    def initiate(self):
        """ Initiate or stop an etMINFLUX experiment. """
        if not self.__running:
            # read mfx sequence and lasers and laser powers from GUI
            sequenceIdxth0 = self._widget.mfxth0_seq_par.currentIndex()
            self.mfxth0_seq = self._widget.mfx_seqs[sequenceIdxth0]
            laserIdxth0 = self._widget.mfxth0_exc_laser_par.currentIndex()
            self.mfxth0_exc_laser = self._widget.mfx_exc_lasers[laserIdxth0]
            self.mfxth0_exc_pwr = float(self._widget.mfxth0_exc_pwr_edit.text())
            detectorIdxth0 = self._widget.mfxth0_detector_par.currentIndex()
            self.mfxth0_detector = self.mfxDetectorListThreads[detectorIdxth0]
            self.mfxch0_detector = self.mfxDetectorList[detectorIdxth0]
            self.mfx_act_pwr = float(self._widget.mfx_act_pwr_edit.text())
            if self.__mfx_twothreaded:
                sequenceIdxth1 = self._widget.mfxth1_seq_par.currentIndex()
                self.mfxth1_seq = self._widget.mfx_seqs[sequenceIdxth1]
                laserIdxth1 = self._widget.mfxth1_exc_laser_par.currentIndex()
                self.mfxth1_exc_laser = self._widget.mfx_exc_lasers[laserIdxth1]
                self.mfxth1_exc_pwr = float(self._widget.mfxth1_exc_pwr_edit.text())
                detectorIdxth1 = self._widget.mfxth1_detector_par.currentIndex()
                self.mfxth1_detector = self.mfxDetectorListThreads[detectorIdxth1]
                self.mfxch1_detector = self.mfxDetectorList[detectorIdxth1]

            # read param for triggering random ROI from binary masks
            self.__random_roi_bin = self._widget.triggerRandomROICheck.isChecked()
            # read params for analysis pipeline from GUI
            self.__pipeline_param_vals = self.readPipelineParams()
            # reset general run parameters
            self.resetRunParams()
            self.clearConfocalData()
            # Reset parameter for extra information that pipelines can input and output
            self.__exinfo = None

            # check if roi following mode, and if so which type
            self.checkFollowROIMode()
            # reset analysis view widget
            self.__analysisHelper.resetWidget()
            self.__multiMFXROIHelper.resetWidget()
            # launch help widget, if visualization mode or validation mode
            experimentModeIdx = self._widget.experimentModesPar.currentIndex()
            experimentMode = self._widget.experimentModes[experimentModeIdx]
            if experimentMode == 'TestVisualize':
                self.__runMode = RunMode.TestVisualize
            elif experimentMode == 'TestValidate':
                self.__runMode = RunMode.TestValidate
            else:
                self.__runMode = RunMode.Experiment
                self.__plotROI = self._widget.plotROICheck.isChecked()
            if self.__runMode == RunMode.TestValidate or self.__runMode == RunMode.TestVisualize or self.__plotROI:
                self.launchHelpWidget()
            self.resetEventViewWidget()
            if self.__runMode == RunMode.Experiment:
                self.launchEventViewWidget()
            # ensure activate configuration is conf overview
            self.activateConfocalConfig()
            # load selected coordinate transform
            self.loadTransform()
            # read confocal image shape parameters
            self.getConfocalShape()
            # read param for line-wise analysis pipeline runs and decide analysis period
            self.__lineWiseAnalysis = self._widget.lineWiseAnalysisCheck.isChecked()
            if self.__lineWiseAnalysis:
                self.__confocalLinesAnalysisPeriod = int(self._widget.lines_analysis_edit.text())
            else:
                self.__confocalLinesAnalysisPeriod = self.__confocalLinesFrame
            # if confocal frame pause, read pause value
            if self.__confocalFramePause:
                self.pausetime = int(self._widget.conf_frame_pause_edit.text()) * 1000  # time in ms
            # read number of initial frames without analysis
            self.__init_frames = int(self._widget.init_frames_edit.text()) - 1
            # start confocal imaging loop
            self.startConfocalScanning()
            self._widget.initiateButton.setText('Stop')
            self.__running = True
            self.setDetLogLine('experiment_start', None)
        elif self.__runningMFX:
            if self.recTimeTimer.isActive():
                self.recTimeTimer.stop()
            self.recTimeTimerThread.quit()
            self.stopMFX()
            self.scanEnded()
        else:
            # save log and confocal data
            filename_prefix = datetime.now().strftime('%y%m%d-%H%M%S')
            if self.__runMode == RunMode.Experiment:
                self.saveDetLog(filename_prefix)
                self.saveValidationImages(prev=True, prev_ana=True, path_prefix=filename_prefix)
            # disconnect signals and reset parameters
            self._imspector.disconnect_end(self.imspectorLineEvent, 1)  # disconnect Imspector confocal image frame finished to running pipeline
            self._imspector.pause()
            self._widget.initiateButton.setText('Initiate')
            self.resetPipelineParamVals()
            self.resetRunParams()

    def checkFollowROIMode(self):
        """ Check if ROI following mode, and if so which type. """
        if self.__followingROI:
            roiFollowModeIdx = self._widget.roiFollowingModesPar.currentIndex()
            roiFollowMode = self._widget.roiFollowingModes[roiFollowModeIdx]
            if roiFollowMode == 'Single':
                self.__followingROIMode = ROIFollowMode.Single
            elif roiFollowMode == 'Multiple':
                self.__followingROIMode = ROIFollowMode.Multiple
            elif roiFollowMode == 'SingleRedetect':
                self.__followingROIMode = ROIFollowMode.SingleRedetect
            elif roiFollowMode == 'MultipleDetect':
                self.__followingROIMode = ROIFollowMode.MultipleDetect

    def clearConfocalData(self):
        """ Clear confocal data deques. """
        self.__prevFrames.clear()
        self.__prevAnaFrames.clear()

    def imspectorLineEvent(self):
        """ Called when a line in a confocal image is finished in Imspector. """
        self.__confocalLineAnalysisCurr += 1
        self.__confocalLineFrameCurr += 1
        analysisSuccess = False
        if self.__confocalLineAnalysisCurr == self.__confocalLinesAnalysisPeriod:
            self.__confocalLineAnalysisCurr = 0
            analysisSuccess, stopFollowingExperiment = self.analysisPeriodTrigger()
        if self.__confocalLineFrameCurr >= self.__confocalLinesFrame:
            self.__confocalLineFrameCurr = 0
            self.__fast_frame += 1
            self.bufferLatestImages()
            self.updateConfocalFrameDisp()
            if self.__confocalFramePause and not analysisSuccess and not stopFollowingExperiment:
                self.confocalIntermittentPause()

    def confocalIntermittentPause(self):
        """ Pause confocal scanning for a certain time in a confocal timelapse without continuous acquisition. """
        if self.confPauseTimer.isActive():
            self.confPauseTimer.stop()
            self.confPauseGUICountdownTimer.stop()
        self.confPauseTimerThread.quit()
        self.confPauseGUICountdownTimerThread.quit()
        self._imspector.pause()
        self.startConfPauseTimer()

    def startConfocalScanning(self, reconnect_signal=True):
        """ Start confocal scanning in Imspector. """
        if reconnect_signal:
            self._imspector.disconnect_end(self.imspectorLineEvent, 1)
            self._imspector.connect_end(self.imspectorLineEvent, 1)
        self.__confocalLineAnalysisCurr = 0
        self.__confocalLineFrameCurr = 0
        # ensure activate configuration is conf overview
        self.activateConfocalConfig()
        time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_base'))
        # enable auto-looping and start scan
        self._imspector.enable_loop(True)
        self._imspector.start()

    def scanEnded(self):
        """ End a MINFLUX acquisition. """
        if self.__followingROI and self.__followingROIMode == ROIFollowMode.Multiple:
            self.setDetLogLine(f"mfx_end-id{self.__roi_events[-1][1]}-cycle{self.__roiFollowCurrCycle}", None)
            self.setDetLogLine("recording_mode", "roi_follow_multiple") 
        elif self.__followingROI and self.__followingROIMode == ROIFollowMode.MultipleDetect:
            self.setDetLogLine(f"mfx_end-id{self.__roi_events[-1][1]}-cycle{self.__roiFollowCurrCycle}", None)
            self.setDetLogLine("recording_mode", "roi_follow_multipledetect")
        elif self.__followingROI:
            self.setDetLogLine(f"mfx_end-cycle{self.__roiFollowCurrCycle}", None) 
            self.setDetLogLine("recording_mode", "roi_follow_single")               
        elif self.__run_all_aoi:
            self.setDetLogLine(f"mfx_end-id{self.__roi_events.maxlen-len(self.__roi_events)}", None)
            self.setDetLogLine("recording_mode", "all_roi") 
        else:
            self.setDetLogLine("mfx_end", None)
            self.setDetLogLine("recording_mode", "single_roi") 
        if self.__plotROI and not self.__followingROI:
            self.deleteROIGUI(idx=0)  # delete top ROI from list
        if (not self.__run_all_aoi or not self.__roi_events) and (not self.__followingROIContinue) and (not self.__followingROI):
            # single trigger experiments
            time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_roiswitch'))
            self.endExperiment()
            time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_roiswitch'))
            self.continueFastModality()
            self.__fast_frame = 0
        elif self.__followingROI and (self.__followingROIMode == ROIFollowMode.Multiple or self.__followingROIMode == ROIFollowMode.MultipleDetect):
            # multiple roi following experiments - start a new ROI, or run a confocal if a loop of saved ROIs has been completed
            time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_roiswitch'))  # to make sure Imspector has properly turned off MINFLUX recording
            if self.__roi_events[0][1] == np.min([event[1] for event in self.__roi_events]):
                self.__roiFollowCurrCycle += 1
                self.endFollowingROIStep()
                time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_roiswitch'))
                self.continueFastModality()
                self._widget.initiateButton.setText('Stop')
                return
            event = self.__roi_events.popleft()
            self.__roi_events.append(event)
            self.newROIMFX(event[0], roi_idx=event[1])
            if self.__followingROIMode == ROIFollowMode.MultipleDetect:
                self.plotMFXROIs(self.__roi_events, self.__roi_sizes)
                self._widget.initiateButton.setText('Next ROI cycle')
        elif self.__followingROI:
            # if in a following ROI experiment and we want to go to the next step with confocal - only for single and singleredetect
            time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_roiswitch'))  # to make sure Imspector has properly turned off MINFLUX recording
            self.__roiFollowCurrCycle += 1
            self.endFollowingROIStep()
            time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_roiswitch'))
            self.continueFastModality()
        elif self.__followingROIContinue:
            # only triggers if followingROI check box have been unclicked to end a followingROI experiment,
            # in which case one more confocal frame will be run, all data will be saved, before resetting everytning.
            time.sleep(2*self.setupInfo.get('timing_settings').get('sleep_time_roiswitch'))
            self.continueFastModality()
            self.__fast_frame = 0
        elif self.__run_all_aoi:
            # if we are not following, but running all ROIs detected
            time.sleep(2*self.setupInfo.get('timing_settings').get('sleep_time_roiswitch'))  # to make sure Imspector has properly turned off MINFLUX recording
            event = self.__roi_events.popleft()
            self.newROIMFX(event[0], roi_idx=event[1])

    def deleteROIGUI(self, idx):
        self._widget.coordListWidget.deleteCoord(idx)
        self._widget.analysisHelpWidget.removeROI(idx)

    def setDetLogLine(self, key, val, *args):
        """ Set a line of the event detection log. """
        if val is None:
            val = datetime.now().strftime('%Hh%Mm%Ss%fus')
        if args:
            self.__detLog[f"{key}{args[0]}"] = val
        else:
            key = self.getUniqueKey(key)
            self.__detLog[key] = val

    def getUniqueKey(self, key):
        """ Get a unique key for the event log. If the key already exists, append a number to it until it is unique. """
        if key not in self.__detLog:
            return key
        else:
            new_key_idx = 1
            new_key = key+str(new_key_idx)
            if new_key not in self.__detLog:
                return new_key
            while new_key in self.__detLog:
                new_key_idx += 1
                new_key = new_key+str(new_key_idx)
            return new_key

    def endExperiment(self):
        """ Save all etMINFLUX scan metadata/log file and data. """
        # save log and confocal data leading up to event image
        filename_prefix = datetime.now().strftime('%y%m%d-%H%M%S')
        self.saveDetLog(filename_prefix)
        self.saveValidationImages(prev=True, prev_ana=True, path_prefix=filename_prefix)
        if self.__autoSaveMeas:
            # save all MINFLUX data (including multiple ROI, if so)
            self.saveMINFLUXdata(path_prefix=filename_prefix)

    def saveDetLog(self, filename_prefix):
        """ Save the event detection log to a .txt file"""
        if self.__detLog:
            self.setDetLogLine("pipeline", self.getPipelineName())
            self.logPipelineParamVals()
            self.setDetLogLine("pipeline_runtimes", self.__pipeline_runtimes)
            if self.__confocalFramePause:
                self.setDetLogLine('confocal_frame_interval', self.pausetime / 1000)  # time in s
            # log recording mode and rec time, if existing
            if self.__hasRunMFX:
                if self.__presetRecTime:
                    self.setDetLogLine("MFXRecTime", self.__rec_time_scan)
                    self.setDetLogLine("roifollow_mfx_runtime", self.__follow_roi_confocal_interval)
            self.setDetLogLine('experiment_end', None)
            # save log file with temporal info of trigger event
            name = os.path.join(self._dataDir, filename_prefix) + '_log'
            savename = getUniqueName(name)
            log = [f'{key}: {self.__detLog[key]}' for key in self.__detLog]
            with open(f'{savename}.txt', 'w') as f:
                [f.write(f'{st}\n') for st in log]
            self.resetDetLog()

    def endFollowingROIStep(self):
        filename_prefix = datetime.now().strftime('%y%m%d-%H%M%S')
        # save confocal data leading up to event image
        self.saveValidationImages(prev=True, prev_ana=True, path_prefix=filename_prefix)

    def getTransformName(self):
        """ Get the name of the pipeline currently used. """
        transformidx = self._widget.transformPipelinePar.currentIndex()
        transformname = self._widget.transformPipelines[transformidx]
        return transformname

    def getPipelineName(self):
        """ Get the name of the pipeline currently used. """
        pipelineidx = self._widget.analysisPipelinePar.currentIndex()
        pipelinename = self._widget.analysisPipelines[pipelineidx]
        return pipelinename

    def logPipelineParamVals(self):
        """ Put analysis pipeline parameter values in the log file. """
        param_names = list()
        for pipeline_param_name, _ in self.__pipeline_params.items():
            if pipeline_param_name not in self.__params_exclude:
                param_names.append(pipeline_param_name)
        for key, val in zip(param_names, self.__pipeline_param_vals):
            self.setDetLogLine(key, val)

    def continueFastModality(self):
        """ Continue confocal imaging, after a MINFLUX event acquisition has been performed. """
        if self._widget.endlessScanCheck.isChecked() or self.__followingROIContinue:
            if not self.__followingROIContinue:
                self.__exinfo = None
                self.resetRunParams()
            self.startConfocalScanning()
            if not self.__followingROIContinue:
                self.resetEventViewWidget()
                self.launchEventViewWidget()
            self._widget.initiateButton.setText('Stop')
            self.__running = True
        else:
            # disconnect signals and reset parameters
            self._imspector.disconnect_end(self.imspectorLineEvent, 1)  # disconnect Imspector confocal image frame finished to running pipeline
            self._widget.initiateButton.setText('Initiate')
            self.resetPipelineParamVals()
            self.resetRunParams()
            self.resetPrevEventsDeque()

    def loadTransform(self):
        """ Load a previously saved coordinate transform. """
        transformname = self.getTransformName()
        self.transform = getattr(importlib.import_module(f'{transformname}'), f'{transformname}')

    def loadPipeline(self):
        """ Load the selected analysis pipeline, and its parameters into the GUI. """
        pipelinename = self.getPipelineName()
        self.pipeline = getattr(importlib.import_module(f'{pipelinename}'), f'{pipelinename}')
        self.__pipeline_params = signature(self.pipeline).parameters
        # get parameter for how many confocal channels the loaded pipeline uses
        self.__count_conf_channels = len([channel for channel in self.__pipeline_params.keys() if 'img_ch' in channel])
        # set GUI fields for pipeline parameters
        self._widget.initParamFields(self.__pipeline_params, self.__params_exclude, self.__count_conf_channels)

    def setAnalysisHelpImg(self, img):
        """ Set the preprocessed image in the analysis help widget. """
        if self.__fast_frame < self.__init_frames + 3:
            autolevels = True
        else:
            autolevels = False
        if self.__runMode == RunMode.TestVisualize:
            self._widget.analysisHelpWidget.img.setImage(img, autoLevels=autolevels)
            infotext = f'Min: {np.floor(np.min(img))}, max: {np.floor(np.max(img))}'
            self._widget.analysisHelpWidget.info_label.setText(infotext)
        else:
            self._widget.multiMFXROIHelpWidget.img.setImage(img, autoLevels=autolevels)
            infotext = f'Min: {np.floor(np.min(img))}, max: {np.floor(np.max(img))}'
            self._widget.multiMFXROIHelpWidget.info_label.setText(infotext)

    def softReset(self):
        """ Reset all parameters after a soft lock. """
        self.setBusyFalse()
        self.resetHelpWidget()
        self.resetEventViewWidget()
        self.resetRunParams()
        #self.initiateFlagsParams()  # TODO likely do not want to do this, just resets everything that has been chosen in the GUI - is this really correct?
        self.resetDetLog()
        self.__prevFrames = deque(maxlen=500)
        self.__prevAnaFrames = deque(maxlen=500)
        self.__roi_events = deque(maxlen=0)
        self.__roi_sizes = deque(maxlen=0)

    def setBusyFalse(self):
        """ Set busy flag to false. """
        self.__busy = False

    def resetHelpWidget(self):
        self._widget.resetHelpWidget()
        self.__analysisHelper = AnalysisImgHelper(self, self._widget.analysisHelpWidget)
        self.__multiMFXROIHelper = AnalysisImgHelper(self, self._widget.multiMFXROIHelpWidget)

    def resetEventViewWidget(self):
        self._widget.resetEventViewWidget()
        self.__eventViewHelper = EventWidgetHelper(self, self._widget.eventViewWidget)

    def readPipelineParams(self):
        """ Read user-provided analysis pipeline parameter values. """
        param_vals = list()
        for item in self._widget.param_edits:
            param_vals.append(float(item.text()))
        return param_vals

    def launchHelpWidget(self):
        """ Launch help widget that shows the preprocessed images in real-time. """
        if self.__runMode == RunMode.TestVisualize:
            self._widget.launchHelpWidget(self._widget.analysisHelpWidget, init=True)
        else:
            self._widget.launchHelpWidget(self._widget.multiMFXROIHelpWidget, init=True)

    def launchEventViewWidget(self):
        """ Launch event view widget that shows the confocal stack around the event. """
        self._widget.launchHelpWidget(self._widget.eventViewWidget, init=True)

    def resetDetLog(self):
        """ Reset the event log dictionary. """
        self.__detLog = dict()

    def resetPipelineParamVals(self):
        """ Reset the pipeline parameters. """
        self.__pipeline_param_vals = list()

    def resetRunParams(self):
        """ Reset general pipeline run parameters. """
        self.__running = False
        self.__validating = False
        self.__hasRunMFX = False
        self.__img_ana = None
        self.__fast_frame = 0
        self.__post_event_frames = 0
        self.__pipeline_runtimes = []
        self.__roi_events = deque(maxlen=0)
        self.__roi_sizes = deque(maxlen=0)
        self.__roiFollowCurrCycle = 0
        self.__followingROIContinue = False
        self.__multipleDetectionMFXInitiation = False
        self.__followingROIMultipleDetectMINFLUXPhase = False
        self._widget.setConfGUINullMessages()
        self.resetTimerThreads()
        
    def resetPrevEventsDeque(self):
        """ Reset deque where previous event coords are saved, in endless. """
        self.__prev_event_coords_deque = deque(maxlen=100)

    def resetTimerThreads(self):
        if self.recTimeTimer.isActive():
            self.recTimeTimer.stop()
        self.recTimeTimerThread.quit()
        if self.confPauseTimer.isActive():
            self.confPauseTimer.stop()
        self.confPauseTimerThread.quit()
        if self.confPauseGUICountdownTimer.isActive():
            self.confPauseGUICountdownTimer.stop()
        self.confPauseGUICountdownTimerThread.quit()
        if self.multipleEventDetectionTimer.isActive():
            self.multipleEventDetectionTimer.stop()
        self.multipleEventDetectionTimerThread.quit()

    def generateRandomCoord(self, roi_sizes):
        """ Generate random coord inside the binary mask. """         
        if self.__binary_mask is not None:
            possible_pixels = np.where(self.__binary_mask==True)
            rand_idx = np.random.randint(0, len(possible_pixels[0]))
            random_coords = np.array([[possible_pixels[0][rand_idx], possible_pixels[1][rand_idx]]])
            random_coords = np.flip(random_coords, axis=1)
            coords_scan = random_coords[0]
            if not self.__presetROISize:
                roi_size = roi_sizes[0]
            else:
                roi_size = None
        else:
            random_coords = np.array([])
            roi_size = None
            coords_scan = None
        self.__coords_scan_curr = coords_scan
        return random_coords, coords_scan, roi_size

    def analysisPeriodTrigger(self):
        """ Triggers when an analysis period has finished. Depending on experiment modes set: send either to runPipeline, or act in other ways. """
        analysisSuccess = False
        stopFollowingExperiment = False
        # enter if not still running last analysis period trigger
        if not self.__busy:
            # set busy true
            self.__busy = True
            # get confocal image channel(s) from Imspector
            meas = self._imspector.active_measurement()
            stackidx = 3
            self.imgs = []
            self.imgs.append(np.squeeze(meas.stack(stackidx).data()[0]))
            if self.__count_conf_channels >= 2:
                # if dual color confocal scan - find ch2 in imspector measurement data stacks (same size as ch1)
                ch1_sh = np.shape(self.imgs[0])
                for i in range(stackidx+1,10):
                    try:
                        img = np.squeeze(meas.stack(i).data()[0])
                        if np.shape(img) == ch1_sh:
                            self.imgs.append(img)
                            ch2_idx = i
                            break
                    except:
                        pass
            if self.__count_conf_channels == 3:
                # if triple color confocal scan - find ch3 in imspector measurement data stacks (same size as ch1)
                for i in range(ch2_idx+1,10):
                    try:
                        img = np.squeeze(meas.stack(i).data()[0])
                        if np.shape(img) == ch1_sh:
                            self.imgs.append(img)
                            break
                    except:
                        pass
            if self.__followingROIContinue and (self.__followingROIMode == ROIFollowMode.Single or self.__followingROIMode == ROIFollowMode.SingleRedetect or self.__followingROIMode == ROIFollowMode.Multiple):
                # if recording mode is following ROI and we are currently following a ROI
                analysisSuccess = True
                if not self.__followingROI:
                    # if following ROI mode check box has been unchecked - end experiment
                    self.endFollowingROIExperiment()
                    stopFollowingExperiment = True
                elif self.__followingROIMode == ROIFollowMode.SingleRedetect:
                    # if we are redetecting a single ROI
                    coords_detected, roi_sizes = self.runPipeline()
                    coords_scan, roi_size = self.postPipelineFollowingROISingleRedetect(coords_detected, roi_sizes)
                    if len(coords_scan)>0:
                        self.acquireMINFLUXFull(coords_scan, roi_size)
                        self.updateEventViewWidget(self.imgs[0].copy(), coords_scan)
                    else:
                        # end experiment, as we no longer have anything to follow, and we cannot suddenly "find it back" again in the next image
                        self.endFollowingROIExperiment()
                elif self.__followingROIMode == ROIFollowMode.Single:
                    self.acquireMINFLUXMinimal(pos=self.__roi_center_mfx, roi_size_um=self.__roi_size_um_mfx, pos_conf=self.__pos_conf_mfx)
                    self.updateEventViewWidget(self.imgs[0].copy())
                elif self.__followingROIMode == ROIFollowMode.Multiple:
                    # initiate new cycle of following ROI
                    event = self.__roi_events.popleft()
                    self.__roi_events.append(event)
                    # pause fast imaging
                    self.pauseFastModality()
                    self.newROIMFX(coords_scan=event[0], roi_idx=event[1])
            elif (self.__followingROI or self.__followingROIEnding) and self.__followingROIMode == ROIFollowMode.MultipleDetect:
                # if we are in multiple detect MINFLUX following ROI mode, no matter if we have detected event previously or not, or if we are just ending such an experiment
                # TODO test if the followingROIEnding parameter added works as intended - two checks: check that we do not enter here when multipledetect is selected but followingROI checkbox is unchecked, and check that we enter here when we are ending a followingROI multidetect experiment, after it having run and detected events
                if not self.__followingROIMultipleDetectMINFLUXPhase:
                    # check number of events detected and compare event limit for MFX phase initiation, if we are not already in MFX phase
                    eventlimit = int(self._widget.mfx_phase_eventlim_edit.text())
                    if len(self.__roi_events) >= eventlimit:
                        self.toggleMultipleDetectionMFXInitiation()
                if self.__followingROIEnding:
                    # if following ROI mode check box has been unchecked - end experiment
                    self.endFollowingROIExperiment()
                    stopFollowingExperiment = True
                elif self.__followingROIMultipleDetectMINFLUXPhase:
                    # start MFX on the first position in the deque
                    event = self.__roi_events.popleft()
                    self.__roi_events.append(event)
                    self.pauseFastModality()
                    time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_base'))
                    self.newROIMFX(coords_scan=event[0], roi_idx=event[1])
                    analysisSuccess = True
                    self.plotMFXROIs(self.__roi_events, self.__roi_sizes)
                    self.updateEventViewWidget(self.imgs[0].copy())
                    self._widget.initiateButton.setText('Next ROI cycle')
                elif self.__multipleDetectionMFXInitiation:
                    # if time started from first event has run out, if such a timer exists, or if the button has been pressed that manually toggles MFX initiation
                    self.__followingROIMultipleDetectMINFLUXPhase = True
                    # initiate MINFLUX recording loop, according to above/below by popping and readding an event to the deque.
                    event = self.__roi_events.popleft()
                    self.__roi_events.append(event)
                    self.pauseFastModality()
                    self.newROIMFX(coords_scan=event[0], roi_idx=event[1])
                    analysisSuccess = True
                    self.plotMFXROIs(self.__roi_events, self.__roi_sizes)
                    self.updateEventViewWidget(self.imgs[0].copy())
                    self._widget.initiateButton.setText('Next ROI cycle')
                elif self.__fast_frame > self.__init_frames:
                    # if initial settling frames have passed, and we are still just in analysis phase
                    # then, run pipeline
                    coords_detected, roi_sizes = self.runPipeline()
                    events_list = list()
                    if len(coords_detected) > 0:
                        # if it is the first event detected, initiate event detection windw by calling self.startMultipleEventDetectionTimer()
                        try:
                            self.__roi_events
                            # get base event index of already detected events
                            event_idx = np.max([event[1] for event in self.__roi_events])+1
                        except:
                            self.startMultipleEventDetectionTimer()
                            event_idx = 0
                        # add event to deque and add log entry for event detection time with event idx
                        _, _, _, new_event = self.postPipelineExperiment(coords_detected, roi_sizes)
                        # update event view widget with the new events
                        # prep lists for coords and ROI sizes
                        if new_event:
                            if np.size(coords_detected) > 2:
                                for pair in coords_detected:
                                    events_list.append([pair, event_idx])
                                    event_idx += 1
                            else:
                                events_list.append([coords_detected[0], event_idx])
                            for event in events_list:
                                if self.__prevFrames:
                                    self.addEventEventViewWidget(event[0], event[1], self.imgs[0].copy(), add_prev=True)
                                else:
                                    self.addEventEventViewWidget(event[0], event[1], self.imgs[0].copy())
                        self.helpImageGeneration(coords_detected, roi_sizes, eventids=[event[1] for event in events_list])
                        analysisSuccess = False  # as we want to keep running confocal, unless we stop it after time has run out after the first event
                    # update event view widget events with the new confocal data, except the new events
                    if events_list:
                        event_ids = [event[1] for event in events_list]
                    else:
                        event_ids = []
                    self.updateEventViewWidget(self.imgs[0].copy(), exclude_events=event_ids)
                    if self.__roi_events:
                        self.plotMFXROIs(self.__roi_events, self.__roi_sizes)
            else:
                # if no following ROI mode continuation or multiDetect following ROI mode, i.e. we are looking for first events, and if initial settling frames have passed
                if self.__fast_frame > self.__init_frames:
                    coords_detected, roi_sizes = self.runPipeline()
                    # if visualization mode
                    if self.__runMode == RunMode.TestVisualize:
                        self.postPipelineVisualize(coords_detected, roi_sizes)
                    # if validation mode
                    elif self.__runMode == RunMode.TestValidate:
                        self.postPipelineValidate(coords_detected, roi_sizes)
                    # if experiment mode
                    elif self.__runMode == RunMode.Experiment:
                        coords_detected, coords_scan, roi_size, _ = self.postPipelineExperiment(coords_detected, roi_sizes)
                        # move on with triggering modality switch, if we have detected or randomized coords, and if detected coord is sufficiently far away from previous events
                        if len(coords_scan)>0:
                            # check if event coords are close to coords in previous event deques
                            if len(self.__prev_event_coords_deque) > 0:
                                dists = [eucl_dist(c_old, coords_scan) for c_old in [event[0] for event in self.__roi_events]]
                                #dists = [np.sqrt((c_old[0]-coords_scan[0])**2+(c_old[1]-coords_scan[1])**2) for c_old in self.__prev_event_coords_deque]
                                if np.min(dists) > self.setupInfo.get('analysis_settings').get('min_distance_prev_events'):
                                    new_event = True
                                else:
                                    new_event = False
                            else:
                                new_event = True
                            if new_event or 'fake' in self.getPipelineName():
                                self.__prev_event_coords_deque.append(coords_scan)
                                self.acquireMINFLUXFull(coords_scan, roi_size)
                                self.helpImageGeneration(coords_detected, roi_sizes)
                                if self.__prevFrames:
                                    self.initiateEventViewWidget(coords_scan)  # add all prevframes
                                    self.updateEventViewWidget(self.imgs[0].copy())  # update with detected frame
                                else:
                                    self.initiateEventViewWidget(coords_scan, self.imgs[0].copy())  # update with detected frame
                                analysisSuccess = True
            # unset busy flag
            self.setBusyFalse()
        return analysisSuccess, stopFollowingExperiment

    def initiateEventViewWidget(self, event_coords, frame=None):
        """ Initiate the event view widget with the current confocal frame and event coordinates. """
        if frame is not None:
            self.__eventViewHelper.new_event(frame, event_coords)
        else:
            self.__eventViewHelper.new_event(self.__prevFrames, event_coords)

    def addEventEventViewWidget(self, event_coords, event_id, frame, add_prev=False):
        if add_prev:
            allframes_deque = self.__prevFrames.copy()
            allframes_deque.append(frame)
            self.__eventViewHelper.new_event(allframes_deque, event_coords, event_id)
        else:
            self.__eventViewHelper.new_event(frame, event_coords, event_id)

    def updateEventViewWidget(self, new_frame, event_coords=None, exclude_events=[]):
        """ Update the event view widget with a new confocal frame and event coordinates. """
        self.__eventViewHelper.new_frame(new_frame, event_coords, exclude_events)

    def runPipeline(self):
        """ Run the analyis pipeline on the latest confocal analysis period frame. """
        # start of pipeline
        self.__pipeline_start_file_prefix = datetime.now().strftime('%Y%m%d-%Hh%Mm%Ss')  # use for naming files
        # run pipeline
        self.__pipeline_start_time = time.perf_counter()
        coords_detected, roi_sizes, self.__exinfo, self.__img_ana = self.pipeline(*self.imgs, self.__prevFrames, self.__binary_mask,
                                                                        self.__exinfo, self.__presetROISize,
                                                                        *self.__pipeline_param_vals)
        self.__pipeline_end_time = time.perf_counter()
        self.__pipeline_runtimes.append(round((self.__pipeline_end_time-self.__pipeline_start_time)*1e3,3))  # in ms with µs precision
        return coords_detected, roi_sizes

    def postPipelineVisualize(self, coords_detected, roi_sizes):
        """ Called if in visualization mode, and plots the coordinates and adds them to the coords list in the helper widget. Generates random coords if selected. """
        # set analysis image in help widget
        self.setAnalysisHelpImg(self.__img_ana)
        # if random roi(s) from binary
        if self.__random_roi_bin:
            coords_detected, _, _ = self.generateRandomCoord(roi_sizes)
        # if we have some detected or randomized coords
        if coords_detected.size != 0:
            self.helpPlotDetectedCoordsSignal.emit(coords_detected, roi_sizes, list(range(len(coords_detected))))

    def plotDetectedCoords(self, coords_detected, roi_sizes):
        """ Plot detected coordinates in the help widget, and add them to the coords list in the help widget. """
        # plot detected coords in help widget
        colors = [self.popColor() for _ in range(len(coords_detected))]
        self.__analysisHelper.plotScatter(coords_detected, colors=colors)
        self.__analysisHelper.plotRoiRectangles(coords_detected, roi_sizes, colors=colors, presetROISize=self.__presetROISize)
        self.__analysisHelper.updateNumberEventsDisp(numEvents=len(coords_detected))
        if self.__presetROISize:
            roi_sizes = self.getPresetRoiSize(len(coords_detected))
        self._widget.coordListWidget.addCoords(coords_detected, roi_sizes, colors)

    def plotMFXROIs(self, coords_deque, roi_sizes_deque):
        # plot MFX ROIs in ROI deque coords in corresponding help widget for multi MFX ROIs
        colors = [self.popColor() for _ in range(len(coords_deque))]
        coords_events = [event[0] for event in coords_deque]
        eventids = [event[1] for event in coords_deque]
        roi_sizes_events = [event[0] for event in roi_sizes_deque]
        self.__multiMFXROIHelper.plotScatter(coords_events, colors=colors)
        self.__multiMFXROIHelper.plotRoiRectangles(coords_events, roi_sizes_events, colors=colors, presetROISize=self.__presetROISize)
        self.__multiMFXROIHelper.updateNumberEventsDisp(numEvents=len(coords_events))
        if self.__presetROISize:
            roi_sizes_events = self.getPresetRoiSize(len(coords_deque))
        self._widget.multiMFXROIcoordListWidget.addCoords(coords_events, roi_sizes_events, eventids, colors)

    def postPipelineValidate(self, coords_detected, roi_sizes):
        """ Called if in validation mode, and initiates validation run or takes care of validation run when finishing. """
        # set analysis image in help widget, and start to record validation frames after event
        self.setAnalysisHelpImg(self.__img_ana)
        # if currently running validation
        if self.__validating:
            if self.__post_event_frames > self.__validation_frames:
                # if all validation frames have been recorded, pause fast imaging, end recording, continue fast imaging
                self.saveValidationImages(prev=True, prev_ana=True, path_prefix=self.__pipeline_start_file_prefix)
                self.__fast_frame = 0
                self.pauseFastModality()
                self.endExperiment()
                self.continueFastModality()
                self.__fast_frame = 0
                self.__validating = False
            self.__post_event_frames += 1
        # initiate validation
        elif coords_detected.size != 0:
            # if some events where detected and not validating plot detected coords in help widget
            self.helpPlotDetectedCoordsSignal.emit(coords_detected, roi_sizes, list(range(len(coords_detected))))
            # take first detected coords as event
            if np.size(coords_detected) > 2:
                coords_scan = coords_detected[0,:]
            else:
                coords_scan = coords_detected[0]
            # log detected center coordinate
            self.setDetLogLine("conf_x_center", coords_scan[0])
            self.setDetLogLine("conf_y_center", coords_scan[1])
            # flag for start of validation
            self.__validating = True
            self.__post_event_frames = 0

    def postPipelineExperiment(self, coords_detected, roi_sizes, dist_lim=3):
        """ Called if in experiment mode, and returns MINFLUX scan coordinates according to GUI settings, detected or random. """
        new_event = False
        # if to generate random roi(s) from binary, do that and use these as other detected coordinates
        if self.__random_roi_bin:
            coords_detected, _, roi_size = self.generateRandomCoord(roi_sizes)
        # handle detected coordinates, if some events were detected
        if coords_detected.size != 0:
            pipelinename = self.getPipelineName()
            # if we want to run all detected ROIs in following mode
            if self.__run_all_aoi or (self.__followingROI and self.__followingROIMode == ROIFollowMode.Multiple):
                # prep deques for coords and ROI sizes
                if not self.__presetROISize:
                    roi_sizes_idx = 0
                # prep lists for coords and ROI sizes
                areas_of_interest = list()
                roi_sizes_list = list()
                if np.size(coords_detected) > 2:
                    for coord_id, pair in enumerate(coords_detected):
                        areas_of_interest.append([pair, coord_id])
                        if self.__presetROISize:
                            roi_sizes_list.append([None, coord_id])
                        else:
                            roi_sizes_list.append([roi_sizes[roi_sizes_idx], coord_id])  # TODO Check if all this works as intended: in the other codebase, m2205, it is event_idx instead of coord_id, despite that parameter not existing
                            roi_sizes_idx += 1
                else:
                    areas_of_interest.append([coords_detected[0], 0])
                    if self.__presetROISize:
                        roi_sizes_list.append([None, 0])
                    else:
                        roi_sizes_list.append([roi_sizes[0], 0])
                self.__roi_events = deque(areas_of_interest, maxlen=len(areas_of_interest))
                self.__roi_sizes = deque(roi_sizes_list, maxlen=len(areas_of_interest))
                event_sz = self.__roi_sizes.popleft()
                event_co = self.__roi_events.popleft()
                if self.__followingROI and self.__followingROIMode == ROIFollowMode.Multiple:
                    self.__roi_sizes.append(event_sz)
                    self.__roi_events.append(event_co)
                coords_scan = event_co[0]
                roi_size = event_sz[0]
                self._widget.initiateButton.setText('Next ROI')
            # if we want to follow multiple ROIs but with multiple ROI detection in a time window
            elif self.__followingROI and self.__followingROIMode == ROIFollowMode.MultipleDetect:
                # get base event index of detected events
                try: # if events detected already, check max index in deque
                    event_idx = np.max([event[1] for event in self.__roi_events])+1
                except:  # else start from 0
                    event_idx = 0
                if not self.__presetROISize:
                    roi_sizes_idx = 0
                # prep lists for coords and ROI sizes
                areas_of_interest = list()
                roi_sizes_list = list()
                prev_events_coordpairs = [event[0] for event in self.__roi_events]
                if np.size(coords_detected) > 2:
                    for pair in coords_detected:
                        dists_prev = [eucl_dist(pair, pairprev) for pairprev in prev_events_coordpairs]
                        if len(dists_prev) > 0:
                            mindist = np.min(dists_prev)
                        else:
                            mindist = 1e6
                        if event_idx == 0 or mindist > dist_lim:
                            areas_of_interest.append([pair, event_idx])
                            if self.__presetROISize:
                                roi_sizes_list.append([None, event_idx])
                            else:
                                roi_sizes_list.append([roi_sizes[roi_sizes_idx], event_idx])  # TODO Check if all this works as intended: in the other codebase, m2205, it is event_idx instead of coord_id, despite that parameter not existing
                                roi_sizes_idx += 1
                            self.setDetLogLine(f"event_detection-id{event_idx}-confframe{len(self.__prevFrames)}")    
                            event_idx += 1
                            new_event = True
                else:
                    dists_prev = [eucl_dist(coords_detected[0], pairprev) for pairprev in prev_events_coordpairs]
                    if event_idx==0 or np.min(dists_prev) > dist_lim:
                        areas_of_interest.append([coords_detected[0], event_idx])
                        if self.__presetROISize:
                            roi_sizes_list.append([None, event_idx])
                        else:
                            roi_sizes_list.append([roi_sizes[0], event_idx])
                        self.setDetLogLine(f"event_detection-id{event_idx}-confframe{len(self.__prevFrames)}")    
                        new_event = True
                # check if deque length is not zero,; if so, append new coords and roi sizes; if zero, create new deques
                if len(self.__roi_events) == 0:
                    self.__roi_events = deque(areas_of_interest, maxlen=40)
                    self.__roi_sizes = deque(roi_sizes_list, maxlen=40)
                else:
                    for event in areas_of_interest:
                        self.__roi_events.append(event)
                    for event in roi_sizes_list:
                        self.__roi_sizes.append(event)
                coords_scan = np.array([])
                roi_size = None
            # if we want to run all detected ROIs once
            elif self.__exinfo is not None and ('_def' in pipelinename or '_rand' in pipelinename):
                # take specified detected coord (and roi_size if applicable) as event (idx = 0 if we want brightest)
                idx = int(self.__exinfo)
                if idx < len(coords_detected):
                    if np.size(coords_detected) > np.max([2,idx]):
                        coords_scan = coords_detected[idx,:]
                    else:
                        coords_scan = coords_detected[0]
                    if not self.__presetROISize:
                        if np.size(coords_detected) > np.max([2,idx]):
                            roi_size = roi_sizes[idx]
                        else:
                            roi_size = roi_sizes[0]
                    else:
                        roi_size = None
                else:  # if defined coordinate is higher than the number of detected coords, take the last detected coord
                    if np.size(coords_detected) > 2:
                        coords_scan = coords_detected[-1,:]
                    else:
                        coords_scan = coords_detected[0]
                    if not self.__presetROISize:
                        roi_size = roi_sizes[-1]
                    else:
                        roi_size = None
                if self._widget.endlessScanCheck.isChecked():
                    self._widget.initiateButton.setText('Next ROI')
            else:
                # take random detected coords (and roi_size if applicable) as event (idx = 0 if we want brightest)
                idx = 0
                random_coord = False
                if random_coord == True:
                    if np.size(coords_detected) > 2:
                        idx = np.random.randint(len(coords_detected))
                        coords_scan = coords_detected[idx,:]
                    else:
                        coords_scan = coords_detected[0]
                elif random_coord == False:
                    if np.size(coords_detected) > 2:
                        coords_scan = coords_detected[idx,:]
                    else:
                        coords_scan = coords_detected[idx]
                if not self.__presetROISize:
                    roi_size = roi_sizes[idx]
                else:
                    roi_size = None
                if self._widget.endlessScanCheck.isChecked():
                    self._widget.initiateButton.setText('Next ROI')
        else:
            coords_detected = np.array([])
            coords_scan = np.array([])
            roi_size = np.array([])
        if len(coords_scan)>0:
            self.__coords_scan_curr = coords_scan
        else:
            self.__coords_scan_curr = None
        return coords_detected, coords_scan, roi_size, new_event

    def postPipelineFollowingROISingleRedetect(self, coords_detected, roi_sizes):
        """ Called if in experiment mode and following ROI mode and following ROI redetection mode, and returns MINFLUX scan coordinates according to GUI settings. """
        self.getRedetectDistThreshold()
        # if we have any newly detected coordinates
        if coords_detected.size != 0:
            min_dist = 1e10
            min_idx = -1
            # check distance to old saved coordinate for all newly detected coordinates
            for idx, coord_new in enumerate(coords_detected):
                dist = math.dist(self.__coords_scan_curr, coord_new)
                if dist < min_dist:
                    min_dist = dist
                    min_idx = idx
            # if the nearest new coordinate is below threshold, consider that our coordinate that has moved
            if min_dist < self.__follow_roi_redetectdist_threshold:
                if np.size(coords_detected) > 2:
                    coords_scan = coords_detected[min_idx,:]
                else:
                    coords_scan = coords_detected[min_idx]
                if not self.__presetROISize:
                    roi_size = roi_sizes[min_idx]
                else:
                    roi_size = None
            # else do not consider that we have a valid coordinate
            else:
                coords_detected = np.array([])
                coords_scan = np.array([])
                roi_size = np.array([])
        else:
            coords_detected = np.array([])
            coords_scan = np.array([])
            roi_size = np.array([])
        self.__coords_scan_curr = coords_scan
        return coords_scan, roi_size

    def acquireMINFLUXMinimal(self, pos, roi_size_um, pos_conf):
        """ Minimal MINFLUX acquisition function. """
        # pause fast imaging
        self.pauseFastModality()
        # initiate and run scanning
        self.initiateMFX(position=pos, ROI_size=roi_size_um, pos_conf=pos_conf)
        self.startMFX()

    def acquireMINFLUXFull(self, coords_scan, roi_size, logging=True):
        """ Full MINFLUX acquisition function. """
        # pause fast imaging
        self.pauseFastModality()
        if logging:
            self.setDetLogLine("coord_transf_start", None)
        # transform detected coordinate between from pixels to global sample position in um (for conf --> MFX)
        coords_center_um = self.transform(coords_scan, self.__confoffset, self.__transformCoeffs) 
        # get roi size in um from preset values
        if self.__presetROISize:
            roi_size_um_scan = [float(self._widget.size_x_edit.text()), float(self._widget.size_y_edit.text())]
        # or calculate roi size in um from in confocal image pixels
        else:
            coord_um_top = self.transform([coords_scan[0]+roi_size[0]/2, coords_scan[1]+roi_size[1]/2], self.__confoffset, self.__transformCoeffs)
            coord_um_bot = self.transform([coords_scan[0]-roi_size[0]/2, coords_scan[1]-roi_size[1]/2], self.__confoffset, self.__transformCoeffs)
            roi_size_um_scan = np.subtract(coord_um_top, coord_um_bot)
        # get preset recording time, if set
        if self.__presetRecTime:
            self.__rec_time_scan = float(self._widget.mfx_rectime_edit.text())
        # else run indefinitely
        else:
            self.__rec_time_scan = None
        # initiate and run scanning with transformed center coordinate
        self.initiateMFX(position=coords_center_um, ROI_size=roi_size_um_scan, pos_conf=coords_scan)
        if logging:
            # log detected and scanning center coordinate, with different keys if running all ROIs or only one
            if self.__followingROI and self.__followingROIMode == ROIFollowMode.Multiple:
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, idx=self.__roi_events[-1][1], cycle=self.__roiFollowCurrCycle)
                self.setDetLogLine(f"mfx_initiate-id{self.__roi_events[-1][1]}-cycle{self.__roiFollowCurrCycle}", None)
            elif self.__followingROI:
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, cycle=self.__roiFollowCurrCycle)
                self.setDetLogLine(f"mfx_initiate-cycle{self.__roiFollowCurrCycle}", None)
            elif self.__run_all_aoi:
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, idx=self.__roi_events[-1][1])
                self.setDetLogLine(f"mfx_initiate-id{self.__roi_events[-1][1]}", None)
            else:
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan)
                self.setDetLogLine("mfx_initiate", None)
        self.startMFX()

    def endFollowingROIExperiment(self):
        """ End experiment after a recalculated following ROI was not detected. """
        self.pauseFastModality()
        self.endExperiment()
        self.__followingROIContinue = False
        self.__followingROIEnding = False
        self.continueFastModality()

    def helpImageGeneration(self, coords_detected, roi_sizes, eventids=None):
        """ Called after starting MINFLUX acquisition, if some additional info regarding detected coordinates should be viewed/saved. """
        time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_base'))
        if self.__plotROI:
            # set analysis image in help widget and plot detected coords in help widget
            self.helpPlotDetectedCoordsSignal.emit(coords_detected, roi_sizes, eventids)
            # set analysis image in help widget
            self.setAnalysisHelpImg(self.__img_ana)

    def deleteROI(self):
        roi_listidx = self._widget.coordListWidget.list.currentRow()
        if roi_listidx > -1:
            self._widget.coordListWidget.deleteCoord(roi_listidx)
            self._widget.analysisHelpWidget.removeROI(roi_listidx)

    def activateConfocalConfig(self):
        meas = self._imspector.active_measurement()
        cfg = meas.configuration(self.setupInfo.get('acquisition_settings').get('confocal_config'))
        meas.activate(cfg)
    
    def activateMinfluxConfig(self):
        meas = self._imspector.active_measurement()
        cfg = meas.configuration(self.setupInfo.get('acquisition_settings').get('minflux_config'))
        meas.activate(cfg)

    def bufferLatestImages(self):
        # buffer latest fast frame and (if applicable) validation images
        self.__prevFrames.append(np.copy(self.imgs[0]))
        # buffer previous preprocessed analysis frame
        if self.__img_ana is not None:
            self.__prevAnaFrames.append(np.copy(self.__img_ana))

    def newROIMFX(self, coords_scan, roi_idx):
        """ Initiate a new MINFLUX ROI scan with the given coordinates and ROI index. """
        # switch active Imspector window to conf overview
        self.activateConfocalConfig()
        # transform detected coordinate between from pixels to sample position in um
        coords_center_um = self.transform(coords_scan, self.__confoffset, self.__transformCoeffs) 
        if self.__presetROISize:
            roi_size_um_scan = [float(self._widget.size_x_edit.text()), float(self._widget.size_y_edit.text())]
        else:
            event_sz = self.__roi_sizes.popleft()
            if (self.__followingROI) and (self.__followingROIMode == ROIFollowMode.Multiple or self.__followingROIMultipleDetectMINFLUXPhase):
                self.__roi_sizes.append(event_sz)
            coord_um_top = self.transform([coords_scan[0]+event_sz[0]/2, coords_scan[1]+event_sz[1]/2], self.__confoffset, self.__transformCoeffs)
            coord_um_bot = self.transform([coords_scan[0]-event_sz[0]/2, coords_scan[1]-event_sz[1]/2], self.__confoffset, self.__transformCoeffs)
            roi_size_um_scan = np.subtract(coord_um_top, coord_um_bot)
        time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_base'))
        # get preset rec time, if using that mode
        if self.__presetRecTime:
            rec_time_scan = float(self._widget.mfx_rectime_edit.text())
            self.__rec_time_scan = rec_time_scan
        else:
            self.__rec_time_scan = None
        # initiate and run scanning with transformed center coordinate
        self.initiateMFX(position=coords_center_um, ROI_size=roi_size_um_scan, pos_conf=coords_scan)
        # log detected and scanning center coordinate
        if self.__followingROI:
            self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, idx=roi_idx, cycle=self.__roiFollowCurrCycle)
            self.setDetLogLine(f"mfx_initiate-id{roi_idx}-cycle{self.__roiFollowCurrCycle}", None) 
        else:
            self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, idx=roi_idx)
            self.setDetLogLine(f"mfx_initiate-id{roi_idx}", None)
        self.startMFX()

    def logCoordinates(self, coords_scan, coords_center_scan, roi_size, idx=None, cycle=None):
        """ Log detection and scan coordinates. """
        if cycle is not None and idx is not None:
            self.setDetLogLine(f"x_center_px-id{idx}-cycle{cycle}", coords_scan[0])
            self.setDetLogLine(f"y_center_px-id{idx}-cycle{cycle}", coords_scan[1])
            self.setDetLogLine(f"x_center_um-id{idx}-cycle{cycle}", coords_center_scan[0])
            self.setDetLogLine(f"y_center_um-id{idx}-cycle{cycle}", coords_center_scan[1])
            self.setDetLogLine(f"x_roi_size_um-id{idx}-cycle{cycle}", roi_size[0])
            self.setDetLogLine(f"y_roi_size_um-id{idx}-cycle{cycle}", roi_size[1])
        elif cycle is not None:
            self.setDetLogLine(f"x_center_px-cycle{cycle}", coords_scan[0])
            self.setDetLogLine(f"y_center_px-cycle{cycle}", coords_scan[1])
            self.setDetLogLine(f"x_center_um-cycle{cycle}", coords_center_scan[0])
            self.setDetLogLine(f"y_center_um-cycle{cycle}", coords_center_scan[1])
            self.setDetLogLine(f"x_roi_size_um-cycle{cycle}", roi_size[0])
            self.setDetLogLine(f"y_roi_size_um-cycle{cycle}", roi_size[1])
        elif idx is not None:
            self.setDetLogLine(f"x_center_px-id{idx}", coords_scan[0])
            self.setDetLogLine(f"y_center_px-id{idx}", coords_scan[1])
            self.setDetLogLine(f"x_center_um-id{idx}", coords_center_scan[0])
            self.setDetLogLine(f"y_center_um-id{idx}", coords_center_scan[1])
            self.setDetLogLine(f"x_roi_size_um-id{idx}", roi_size[0])
            self.setDetLogLine(f"y_roi_size_um-id{idx}", roi_size[1])
        else:
            self.setDetLogLine("x_center_px", coords_scan[0])
            self.setDetLogLine("y_center_px", coords_scan[1])
            self.setDetLogLine("x_center_um", coords_center_scan[0])
            self.setDetLogLine("y_center_um", coords_center_scan[1])
            self.setDetLogLine("x_roi_size_um", roi_size[0])
            self.setDetLogLine("y_roi_size_um", roi_size[1])

    def initiateMFX(self, position=[0.0, 0.0], ROI_size=[1.0, 1.0], pos_conf=[0, 0]):
        """ Prepare a MINFLUX scan at the defined global µm position. """
        self.setMFXROI(position, ROI_size)
        # always remove threads down to a single thread
        for _ in range(2):
            self._imspector.value_at('Minflux/threads/rem', specpy.ValueTree.Measurement).trigger()
            time.sleep(0.05)
        if self.__mfx_twothreaded:
            # if two-threaded minflux, prepare the minflux acquisition with an extra thread and the start condition for the second thread
            self._imspector.value_at('Minflux/threads/add', specpy.ValueTree.Measurement).trigger()
            time.sleep(0.05)
            self._imspector.value_at('Minflux/threads/settings/1/sta', specpy.ValueTree.Measurement).set('on_thi0_vfn')
        self.setMFXSequence(self.mfxth0_seq, thread=0)
        if self.__mfx_twothreaded:
            # if two-threaded minflux, set the second MFX sequence
            self.setMFXSequence(self.mfxth1_seq, thread=1)
        if self.__followingROI:
            if self.__followingROIMode == ROIFollowMode.Single or self.__followingROIMode == ROIFollowMode.SingleRedetect:
                # parameter save current parameters, to reuse for later rounds
                self.__roi_center_mfx = position
                self.__roi_size_um_mfx = ROI_size
                self.__pos_conf_mfx = pos_conf
            self.__followingROIContinue = True
            # get conf interval
            self.getConfocalInterval()
            # set mfx rec time to conf interval - 2*rec_deadtimes
            mfx_interval = self.__follow_roi_confocal_interval - self.__rec_time_deadtime * 2
            self.setMFXRecTime(mfx_interval)
            # set rectimetimer interval and start timer thread
            self.startRecTimer(deadtime=False, time=mfx_interval+self.__rec_time_deadtime)
        elif self.__presetRecTime:
            self.setMFXRecTime()
            self.startRecTimer()
        else:
            self.setMFXRecTime(0)
        if self.__followingROI and (self.__followingROIMode == ROIFollowMode.Multiple or self.__followingROIMode == ROIFollowMode.Single or self.__followingROIMode == ROIFollowMode.SingleRedetect):
            self.setMFXDataTag(pos_conf, ROI_size, self.__roi_events[-1][1], self.__roiFollowCurrCycle)
        elif self.__followingROI and self.__followingROIMode == ROIFollowMode.MultipleDetect:
            self.setMFXDataTag(pos_conf, ROI_size, self.__roi_events[-1][1], self.__roiFollowCurrCycle)
        else:
            self.setMFXDataTag(pos_conf, ROI_size, self.__roi_events.maxlen-len(self.__roi_events))
        if self._act_lasers_present:
            self.setMFXLasers(self.mfxth0_exc_laser, self.mfxth0_exc_pwr, 0, self.mfx_act_pwr)
        else:
            self.setMFXLasers(self.mfxth0_exc_laser, self.mfxth0_exc_pwr, thread=0)
        channelDetectors = []
        if len(self.mfxch0_detector) > 4:
            channelDetectors.append(self.mfxch0_detector.split('"')[1])
            channelDetectors.append(self.mfxch0_detector.split('"')[3])
        else:
            channelDetectors.append(self.mfxch0_detector)
        if self.__mfx_twothreaded:
            if len(self.mfxch1_detector) > 4:
                channelDetectors.append(self.mfxch1_detector.split('"')[1])
                channelDetectors.append(self.mfxch1_detector.split('"')[3])
            else:
                channelDetectors.append(self.mfxch1_detector)
            channelDetectors = (np.unique(channelDetectors)).tolist()
        elif len(channelDetectors)==1:
            if self.mfxch0_detector != self.mfxDetectorList[0]:
                channelDetectors.append(self.mfxDetectorList[0])
            else:
                channelDetectors.append(self.mfxDetectorList[1])
        self.setMFXDetectorsChannels(channelDetectors)
        self.setMFXDetectors(self.mfxth0_detector, thread=0)
        if self.__mfx_twothreaded:
            # if two-threaded minflux, set second thread minflux laser
            self.setMFXLasers(self.mfxth1_exc_laser, self.mfxth1_exc_pwr, thread=1)
            self.setMFXDetectors(self.mfxth1_detector, thread=1)
        time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_base'))

    def startRecTimer(self, deadtime=True, time=None):
        if deadtime:
            self.rec_time_scan = int((self.__rec_time_scan + self.__rec_time_deadtime) * 1000)  # time in ms
        else:
            self.rec_time_scan = int(time * 1000)  # time in ms
        self.recTimeTimer.setInterval(self.rec_time_scan)
        self.recTimeTimerThread.start()

    def startConfPauseTimer(self):
        self.confPauseTimer.setInterval(self.pausetime)
        self.confPauseTimerThread.start()
        # for display countdown timer
        self.confPauseGUICountdownTimer.setInterval(300)
        self.confPauseGUICountdownTimerThread.start()
        self.confIntervalDispTimer.start()

    def startMultipleDetectMFXInitiationTimer(self):
        time_window = int(1000 * 60 * float(self._widget.multiple_detect_window_edit.text()))  # min to s to ms
        self.multipleEventDetectionTimer.setInterval(time_window)
        self.multipleEventDetectionTimerThread.start()

    def stopMultipleDetectMFXInitiationTimer(self):
        if self.multipleEventDetectionTimer.isActive():
            self.multipleEventDetectionTimer.stop()
        self.multipleEventDetectionTimerThread.quit()

    def getConfocalInterval(self):
        self.__follow_roi_confocal_interval = int(self._widget.follow_roi_interval_edit.text())

    def getRedetectDistThreshold(self):
        self.__follow_roi_redetectdist_threshold = int(self._widget.follow_roi_redetectthresh_edit.text())

    def updateIntervalTimerDisp(self):
        time_left = (int(self._widget.conf_frame_pause_edit.text()) * 1000 - self.confIntervalDispTimer.elapsed())/1000
        if time_left >= 0:
            self._widget.conf_guipausetimer_edit.setText(f'Time until next confocal: {time_left:.0f} s')

    def updateConfocalFrameDisp(self):
        self._widget.conf_frame_edit.setText(f'Confocal frames acquired: {self.__fast_frame+1}')

    def deleteMFXDataset(self, times=1):
        mouse.move(*self.__coordTransformHelper._set_topmfxdataset_button_pos)
        for _ in range(times):
            time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_base')/3)
            mouse.click()
            self.keyboard.press(Key.delete)
            self.keyboard.release(Key.delete)

    def getConfocalShape(self):
        """ Get current confocal scan shape, offset, and pixel resolution from Imspector. """
        x_roisize = self._imspector.value_at('ExpControl/scan/range/x/len', specpy.ValueTree.Measurement).get()*1e6
        y_roisize = self._imspector.value_at('ExpControl/scan/range/y/len', specpy.ValueTree.Measurement).get()*1e6
        x_pixels = self._imspector.value_at('ExpControl/scan/range/x/res', specpy.ValueTree.Measurement).get()
        y_pixels = self._imspector.value_at('ExpControl/scan/range/y/res', specpy.ValueTree.Measurement).get()
        self.__confoffset = [self._imspector.value_at('ExpControl/scan/range/x/off', specpy.ValueTree.Measurement).get()*1e6,
                             self._imspector.value_at('ExpControl/scan/range/y/off', specpy.ValueTree.Measurement).get()*1e6]
        self.__confocalLinesFrame = y_pixels - 5
        self.__transformCoeffs = [x_roisize, y_roisize, x_pixels, y_pixels]

    def setMFXSequence(self, mfx_seq, thread=0):
        """ Sets MINFLUX sequence, according to the GUI choice of the user. """
        self._imspector.value_at('Minflux/threads/settings/'+str(thread)+'/seq', specpy.ValueTree.Measurement).set(mfx_seq)

    def setMFXDataTag(self, position, roi_size, roi_idx, cycle=None):
        """ Sets MINFLUX data tag, according to the event detection. """
        datatag = 'ROI'+str(roi_idx)+'-Pos['+str(position[0])+','+str(position[1])+']'+'-Size['+f'{roi_size[0]:.2f}'+','+f'{roi_size[1]:.2f}'+']'
        if self.__followingROI:
            rectime = float(self.__follow_roi_confocal_interval - self.__rec_time_deadtime * 2)
            datatag = datatag + '-RecTime['+str(rectime)+'s]'
            if cycle is not None:
                datatag = datatag + '-Cycle['+str(self.__roiFollowCurrCycle)+']'
        elif self.__presetRecTime:
            datatag = datatag + '-RecTime['+str(self.__rec_time_scan)+'s]'
        self._imspector.value_at('Minflux/tag', specpy.ValueTree.Measurement).set(datatag)

    def setMFXRecTime(self, rec_time_int=None):
        if rec_time_int is None:
            # setting rec time
            rec_time_int = int(self.__rec_time_scan)
        self._imspector.value_at('Minflux/flow/stop_time', specpy.ValueTree.Measurement).set(rec_time_int)

    def setMFXLasers(self, exc_laser, exc_pwr, thread=0, act_pwr=None):
        """ Sets MINFLUX lasers and laser powers, according to the GUI choice of the user. """
        self._imspector.value_at('Minflux/threads/settings/'+str(thread)+'/exc', specpy.ValueTree.Measurement).set(exc_laser)
        self._imspector.value_at('Minflux/threads/settings/'+str(thread)+'/exp', specpy.ValueTree.Measurement).set(exc_pwr)
        # set activation laser (405)
        if act_pwr != None:
            laser_status_act = True if act_pwr > 0 else False
            # TODO: CHECK IF THIS WORKS AS INTENDED. SHOULD BE CORRECT FOR M2205 - HOW ABOUT FOR M2410?
            self._imspector.value_at('ExpControl/measurement/channels/0/lasers/0/active', specpy.ValueTree.Measurement).set(laser_status_act)
            self._imspector.value_at('ExpControl/measurement/channels/0/lasers/0/power/calibrated', specpy.ValueTree.Measurement).set(act_pwr)

    def setMFXDetectorsChannels(self, detectors):
        """ Sets available MINFLUX detectors in the channels, according to the GUI choice of the user. """
        self._imspector.value_at('ExpControl/measurement/channels/0/detsel/detector', specpy.ValueTree.Measurement).set(detectors[0])
        self._imspector.value_at('ExpControl/measurement/channels/1/detsel/detector', specpy.ValueTree.Measurement).set(detectors[1])

    def setMFXDetectors(self, detector, thread=0):  #, act_pwr):
        """ Sets MINFLUX detector(s), according to the GUI choice of the user. """
        #print([thread, detector])
        self._imspector.value_at('Minflux/threads/settings/'+str(thread)+'/det', specpy.ValueTree.Measurement).set(detector)

    def setMFXROI(self, position, ROI_size):
        """ Set the MINFLUX ROI by specpy, after activating the MFX configuration. """
        time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_base')/6)
        self.activateMinfluxConfig()
        self._imspector.value_at('ExpControl/scan/range/x/len', specpy.ValueTree.Measurement).set(ROI_size[0]*1e-6)  # input units, m
        self._imspector.value_at('ExpControl/scan/range/y/len', specpy.ValueTree.Measurement).set(ROI_size[1]*1e-6)  # input units, m
        self._imspector.value_at('ExpControl/scan/range/x/off', specpy.ValueTree.Measurement).set(position[0]*1e-6)  # input units, m
        self._imspector.value_at('ExpControl/scan/range/y/off', specpy.ValueTree.Measurement).set(position[1]*1e-6)  # input units, m

    def startMFX(self):
        """ Run event-triggered MINFLUX acquisition in small ROI. """
        self._imspector.start()
        self.__runningMFX = True
        self.__hasRunMFX = True

    def stopMFX(self):
        """ Stop MINFLUX measurement. """
        self._imspector.pause()  # it is not this causing it
        self.__runningMFX = False

    def saveValidationImages(self, prev=True, prev_ana=True, path_prefix='YMD-HMS'):
        """ Save the validation fast images of an event detection, fast images and/or preprocessed analysis images. """
        if prev:
            self._saveImage(self.__prevFrames, path_prefix, 'conf-raw')
        if prev_ana:
            self._saveImage(self.__prevAnaFrames, path_prefix, 'conf-analysisprocessed')
        self.clearConfocalData()

    def saveMINFLUXdata(self, path_prefix='YMD-HMS'):
        meas = self._imspector.active_measurement()
        fileName = path_prefix + '_' + 'minflux' + '.msr'
        filePath = os.path.join(self._dataDir, fileName)
        meas.save_as(filePath)
        if self.__autoDelMFX:
            time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_save'))
            self.deleteMFXDataset(times=10)  # deletes 10 datasets, which should always be enough; no simple way to detect how many datasets there are

    def saveMeasurement(self):
        filename_prefix = datetime.now().strftime('%y%m%d-%H%M%S')
        self.saveMINFLUXdata(filename_prefix)

    def _saveImage(self, img, path_prefix, path_suffix):
        fileName = path_prefix + '_' + path_suffix + '.tif'
        filePath = os.path.join(self._dataDir, fileName)
        tiff.imwrite(filePath, img)

    def pauseFastModality(self):
        """ Pause the fast method, when an event has been detected. """
        self._imspector.enable_loop(False)
        if self.__running:
            self._imspector.disconnect_end(self.imspectorLineEvent, 1)
            self._imspector.pause()

    def togglePresetROISize(self):
        if not self._widget.presetROISizeCheck.isChecked():
            self._widget.size_x_edit.setEditable(False)
            self._widget.size_y_edit.setEditable(False)
            self.__presetROISize = False
        else:
            self._widget.size_x_edit.setEditable(True)
            self._widget.size_y_edit.setEditable(True)
            self.__presetROISize = True

    def toggleTwoThreadMFX(self):
        if not self._widget.twoThreadsMFXCheck.isChecked():
            self._widget.mfxth1_seq_par.setEnabled(False)
            self._widget.mfxth1_exc_laser_par.setEnabled(False)
            self._widget.mfxth1_exc_pwr_edit.setEditable(False)
            self._widget.mfxth1_detector_par.setEnabled(False)
            self._widget.mfxth0_seq_label.setText('MFX sequence')
            self._widget.mfxth0_exc_laser_label.setText('MFX exc laser')
            self._widget.mfxth0_exc_pwr_label.setText('MFX exc power (%)')
            self._widget.mfxth0_detector_label.setText('MFX detector(s)')
            self._widget.mfxth1_seq_label.setStyleSheet('color: rgb(83,83,83);')
            self._widget.mfxth1_exc_laser_label.setStyleSheet('color: rgb(83,83,83);')
            self._widget.mfxth1_exc_pwr_label.setStyleSheet('color: rgb(83,83,83);')
            self._widget.mfxth1_detector_label.setStyleSheet('color: rgb(83,83,83);')
            self._widget.mfxth1_seq_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(83,83,83); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
            self._widget.mfxth1_exc_laser_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(83,83,83); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
            self._widget.mfxth1_detector_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(83,83,83); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
            self.__mfx_twothreaded = False
        else:
            self._widget.mfxth1_seq_par.setEnabled(True)
            self._widget.mfxth1_exc_laser_par.setEnabled(True)
            self._widget.mfxth1_exc_pwr_edit.setEditable(True)
            self._widget.mfxth1_detector_par.setEnabled(True)
            self._widget.mfxth0_seq_label.setText('MFX (th0) sequence')
            self._widget.mfxth0_exc_laser_label.setText('MFX (th0) exc laser')
            self._widget.mfxth0_exc_pwr_label.setText('MFX (th0) exc power (%)')
            self._widget.mfxth0_detector_label.setText('MFX (th0) detector(s)')
            self._widget.mfxth1_seq_label.setStyleSheet('color: rgb(217,83,0);')
            self._widget.mfxth1_exc_laser_label.setStyleSheet('color: rgb(217,83,0);')
            self._widget.mfxth1_exc_pwr_label.setStyleSheet('color: rgb(217,83,0);')
            self._widget.mfxth1_detector_label.setStyleSheet('color: rgb(217,83,0);')
            self._widget.mfxth1_seq_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(170,170,170); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
            self._widget.mfxth1_exc_laser_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(170,170,170); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
            self._widget.mfxth1_detector_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(170,170,170); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
            self.__mfx_twothreaded = True

    def toggleConfocalFramePause(self):
        if not self._widget.confocalFramePauseCheck.isChecked():
            self._widget.conf_frame_pause_edit.setEditable(False)
            self.__confocalFramePause = False
        else:
            self._widget.conf_frame_pause_edit.setEditable(True)
            self.__confocalFramePause = True

    def togglePresetRecTime(self):
        if not self._widget.presetMfxRecTimeCheck.isChecked():
            self._widget.mfx_rectime_edit.setEditable(False)
            self.__presetRecTime = False
        else:
            self._widget.mfx_rectime_edit.setEditable(True)
            self.__presetRecTime = True

    def togglePresetAutoSave(self):
        if not self._widget.autoSaveCheck.isChecked():
            self.__autoSaveMeas = False
        else:
            self.__autoSaveMeas = True

    def togglePresetAutoDeleteMFX(self):
        if not self._widget.autoDeleteMFXDatasetCheck.isChecked():
            self.__autoDelMFX = False
        else:
            self.__autoDelMFX = True

    def toggleFollowingROI(self):
        if not self._widget.followROIModeCheck.isChecked():
            if self.__presetRecTime:
                self._widget.mfx_rectime_edit.setEditable(True)
            self._widget.presetMfxRecTimeCheck.setEnabled(True)
            self.__followingROI = False
            if self.__followingROIContinue:
                self.__followingROIEnding = True
        else:
            self._widget.mfx_rectime_edit.setEditable(False)
            self._widget.presetMfxRecTimeCheck.setEnabled(False)
            self.__followingROI = True
            self.__followingROIEnding = False

    def toggleTriggerAllROIs(self):
        if not self._widget.triggerAllROIsCheck.isChecked():
            self.__run_all_aoi = False
        else:
            self.__run_all_aoi = True

    def toggleMultipleDetectionMFXInitiation(self):
        self.__multipleDetectionMFXInitiation = True
        self.stopMultipleDetectMFXInitiationTimer()

    def getPresetRoiSize(self, length):
        roi_sizes = []
        for _ in range(length):
            roi_sizes.append([float(self._widget.size_x_edit.text()), float(self._widget.size_y_edit.text())])
        return roi_sizes

    def loadSetupJson(self, filepath=None):
        """ Load a setup from a json file. """
        filepath = QtWidgets.QFileDialog.getOpenFileName(caption='Load etMINFLUX setup info', filter='JSON files (*.json)')[0]
        with open(filepath, 'r') as f:
            setup_dict = json.load(f)
        return setup_dict


class AnalysisImgHelper():
    """ Analysis image widget help controller. """
    def __init__(self, etMINFLUXController, analysisWidget, *args, **kwargs):
        self.etMINFLUXController = etMINFLUXController
        self._widget = analysisWidget
        # connect signals from widget
        self._widget.setLevelsButton.clicked.connect(self.setLevels)

    def setLevels(self):
        """ Set analysis image min,max levels from ROI. """
        min_val = float(self._widget.levelMinEdit.text())
        max_val = float(self._widget.levelMaxEdit.text())
        self._widget.img.setLevels([min_val, max_val])

    def plotScatter(self, coords, colors):
        self._widget.scatterPlot.setData(x=[coord[0] for coord in coords], y=[coord[1] for coord in coords], pen=pg.mkPen(None), brush=colors, symbol='x', size=8)

    def plotRoiRectangles(self, coords, roi_sizes, colors, presetROISize):
        # remove previously drawn ROIs
        self._widget.removeROIs()
        # create rectangle items for each ROI
        if presetROISize:
            µm_px_size = 1   # µm to pixels  # TODO: fix this, before reading from gui transformation coordinates
            roi_size_fix = [float(self.etMINFLUXController._widget.size_x_edit.text())*µm_px_size, float(self.etMINFLUXController._widget.size_y_edit.text())*µm_px_size]
            roi_sizes = [roi_size_fix for _ in coords]
        for idx, [coord, roi_size] in enumerate(zip(coords, roi_sizes)):
            color = colors[idx]
            roi_temp = pg.PlotCurveItem(x=[coord[0]-roi_size[0]/2,coord[0]+roi_size[0]/2,coord[0]+int(roi_size[0]/2),coord[0]-int(roi_size[0]/2),coord[0]-int(roi_size[0]/2)], y=[coord[1]-roi_size[1]/2,coord[1]-roi_size[1]/2,coord[1]+roi_size[1]/2,coord[1]+roi_size[1]/2,coord[1]-roi_size[1]/2], pen=pg.mkPen(color, width=2))
            self._widget.rois_draw.append(roi_temp)
        # draw ROIs
        self._widget.drawROIs()

    def updateNumberEventsDisp(self, numEvents):
        self.etMINFLUXController._widget.coordListWidget.numevents_edit.setText(f'Number of detected events: {numEvents}')

    def resetWidget(self):
        self._widget.reset()
        self._widget.scatterPlot.setData(x=[], y=[], pen=pg.mkPen(None))
        self.etMINFLUXController._widget.multiMFXROIcoordListWidget.numevents_edit.setText(f'Number of detected events: -')


class EventWidgetHelper():
    """ Event widget helper, with confocal image viewing and intensity plot; help controller. """
    def __init__(self, etMINFLUXController, eventWidget, *args, **kwargs):
        self.etMINFLUXController = etMINFLUXController
        self._widget = eventWidget
        self.zoom_size = 10
        self.intensity_size = 2
        self.frames = []
        self.intensity_traces = []
        self.event_coords = []
        self.currentDisplayListIdx = -1
        self.dropdownItems = []

        # connect signals from widget
        self._widget.image_viewer.setLevelsButton.clicked.connect(self.setLevels)
        self._widget.eventsPar.currentIndexChanged.connect(self.changeDisplayedEvent)

    def new_event(self, frames, coords, event_id):
        if len(np.shape(frames))==3:
            frames_arr = np.array(frames)
        else:
            frames_arr = np.array([frames])
        zoomstack = frames_arr[:, int(coords[1]-self.zoom_size):int(coords[1]+self.zoom_size+1), int(coords[0]-self.zoom_size):int(coords[0]+self.zoom_size+1)]
        intensity_trace = np.mean(zoomstack[:, int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1), int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1)], axis=(1,2))
        event_frame = len(frames_arr)-1
        self.event_coords.append([coords, event_id, event_frame])
        self.intensity_traces.append([intensity_trace, event_id, event_frame])
        self.frames.append([zoomstack, event_id, event_frame])
        self.addEventDropdown(event_id)

    def new_frame(self, frame, coords, exclude_events=[]):
        if coords != None:
            zoomframe = frame[int(coords[1]-self.zoom_size):int(coords[1]+self.zoom_size+1), int(coords[0]-self.zoom_size):int(coords[0]+self.zoom_size+1)]
            newintensity = np.mean(zoomframe[int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1), int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1)], axis=(0,1))
            self.intensity_traces[0] = np.append(self.intensity_traces[0], newintensity)
            self.frames[0] = np.append(self.frames[0], zoomframe)
        else:
            # add frame to all events already present, excluding the freshly added ones with ids in exclude_events
            for event_frames, event_trace, event_coords in zip(self.frames, self.intensity_traces, self.event_coords):
                if event_frames[1] not in exclude_events:
                    zoomframe = frame[int(event_coords[0][1]-self.zoom_size):int(event_coords[0][1]+self.zoom_size+1), int(event_coords[0][0]-self.zoom_size):int(event_coords[0][0]+self.zoom_size+1)]
                    newintensity = np.mean(zoomframe[int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1), int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1)], axis=(0,1))
                    event_frames[0] = np.stack((*event_frames[0], zoomframe))
                    event_trace[0] = np.append(event_trace[0], newintensity)
        if len(self.frames) > 0:
            self.changeDisplayedEvent(self.currentDisplayListIdx)

    def setLevels(self):
        """ Set zoom image min,max levels. """
        min_val = float(self._widget.image_viewer.levelMinEdit.text())
        max_val = float(self._widget.image_viewer.levelMaxEdit.text())
        self._widget.image_viewer.img.setLevels([min_val, max_val])

    def addEventDropdown(self, event_id):
        event_text = f'Event id {event_id}'
        self._widget.eventsPar.addItem(event_text)
        self.dropdownItems.append(event_text)

    def removeEvent(self, event_id):
        event_text = f'Event id {event_id}'
        event_listindex = self.dropdownItems.index(event_text)
        self._widget.eventsPar.removeItem(event_listindex)
        self.dropdownItems.remove(event_text)
        frames_index = [frames_event[1] for frames_event in self.frames].index(event_id)
        del self.frames[frames_index]
        del self.intensity_traces[frames_index]
        del self.event_coords[frames_index]

    def changeDisplayedEvent(self, event_list_index):
        """ Change displayed event in the event widget. """
        self.currentDisplayListIdx = event_list_index  # keeps track of the list_index of the event currently being displayed
        zoomstack = self.frames[event_list_index][0]
        intensity_trace = self.intensity_traces[event_list_index][0]
        event_frame = self.frames[event_list_index][2]
        self._widget.reset()
        self._widget.add_event(zoomstack, intensity_trace, event_frame)

    def resetEvents(self):
        self.frames = []
        self.intensity_traces = []
        self.event_coords = []
        self._widget.reset()


class EtCoordTransformHelper():
    """ Coordinate transform help widget controller. """
    def __init__(self, etMINFLUXController, coordTransformWidget, *args, **kwargs):

        self.etMINFLUXController = etMINFLUXController
        self._widget = coordTransformWidget

        # connect signals from widget
        self.etMINFLUXController._widget.coordTransfCalibButton.clicked.connect(self.calibrationLaunch)
        self._widget.setDeleteMFXDatasetButton.clicked.connect(self.setDeleteMFXDatasetButtonCall)

        # set parameters for automatic mouse control
        self._widget.setDeleteMFXDatasetButtonText(self.etMINFLUXController._set_topmfxdataset_button_pos)
        self._set_topmfxdataset_button_pos = self.etMINFLUXController._set_topmfxdataset_button_pos

    def setDeleteMFXDatasetButtonCall(self):
        mouse.on_click(self.setTopMFXDatasetPos)
        mouse.wait(button='left')
        time.sleep(self.etMINFLUXController.setupInfo.get('timing_settings').get('sleep_time_base'))
        self._widget.setDeleteMFXDatasetButtonText(self._set_topmfxdataset_button_pos)

    def setTopMFXDatasetPos(self):
        mouse_pos = mouse.get_position()
        self._set_topmfxdataset_button_pos = np.array(mouse_pos)
        mouse.unhook_all()

    def calibrationLaunch(self):
        """ Launch calibration. """
        self.etMINFLUXController._widget.launchHelpWidget(self.etMINFLUXController._widget.coordTransformWidget, init=True)


class BinaryMaskHelper():
    """ Binary mask help widget controller. """
    def __init__(self, etMINFLUXController, binaryMaskWidget, *args, **kwargs):

        self.etMINFLUXController = etMINFLUXController
        self._widget = binaryMaskWidget

        # connect signals from widget
        self.etMINFLUXController._widget.binaryMaskButton.clicked.connect(self.binaryMaskLaunch)
        self._widget.recordBinaryMaskButton.clicked.connect(self.initiateBinaryMask)
        self._widget.resetBinaryMaskButton.clicked.connect(self.resetBinaryMask)

    def binaryMaskLaunch(self):
        """ Launch calibration. """
        self.etMINFLUXController._widget.launchHelpWidget(self.etMINFLUXController._widget.binaryMaskWidget, init=True)

    def initiateBinaryMask(self):
        """ Initiate the process of calculating a binary mask of the region of interest. """
        self.etMINFLUXController.launchHelpWidget()
        self.binary_stack = []
        self.binary_frame = 0
        self.etMINFLUXController.confocalLinesAnalysisPeriod = self.etMINFLUXController._imspector.active_measurement().stack(0).sizes()[0] - 7
        self.etMINFLUXController.confocalLineCurr = 0
        self.etMINFLUXController._imspector.connect_end(self.imspectorLineEventBinaryMask, 1)
        self.startConfocalScanning(reconnect_signal=False)
        self._widget.recordBinaryMaskButton.setText('Recording...')

    def resetBinaryMask(self):
        """ Reset binary mask, back to None. """
        self.binary_stack = []
        self.etMINFLUXController.binary_mask = None  # None to consider the whole image

    def imspectorLineEventBinaryMask(self):
        self.etMINFLUXController.confocalLineCurr += 1
        if self.etMINFLUXController.confocalLineCurr == self.etMINFLUXController.confocalLinesAnalysisPeriod:
            self.etMINFLUXController.confocalLineCurr = 0
            # get image
            meas = self.etMINFLUXController._imspector.active_measurement()
            self.binary_stack.append(np.squeeze(meas.stack(0).data()))
            self.binary_frame += 1
            if self.binary_frame == self.etMINFLUXController.binary_frames:
                self.calculateBinaryMask()

    def calculateBinaryMask(self):
        """ Calculate the binary mask of the region of interest, using both positive and negative masks. """
        self.etMINFLUXController._imspector.pause()
        img_mean = np.mean(self.etMINFLUXController.binary_stack, 0)
        img_bin_pos = ndi.filters.gaussian_filter(img_mean, float(self._widget.bin_smooth_edit.text()))
        img_bin_neg = ndi.filters.gaussian_filter(img_mean, float(self._widget.bin_neg_smooth_edit.text()))
        img_bin_pos = np.array(img_bin_pos > float(self._widget.bin_thresh_edit.text()))
        img_bin_neg = np.array(img_bin_neg < float(self._widget.bin_neg_thresh_edit.text()))
        img_bin_border = np.zeros(np.shape(img_bin_pos))
        border_size = int(self._widget.bin_border_size_edit.text())
        img_bin_border[border_size:-border_size, border_size:-border_size] = 1
        self.etMINFLUXController.binary_mask = np.logical_and(np.logical_and(img_bin_pos, img_bin_neg), img_bin_border)
        self._widget.recordBinaryMaskButton.setText('Record binary mask')
        self.etMINFLUXController.setAnalysisHelpImg(self.etMINFLUXController.binary_mask)
        self.etMINFLUXController._imspector.disconnect_end(self.imspectorLineEventBinaryMask, 1)


class RunMode(enum.Enum):
    Experiment = 1
    TestVisualize = 2
    TestValidate = 3


class ROIFollowMode(enum.Enum):
    Single = 1
    Multiple = 2
    SingleRedetect = 3
    MultipleDetect = 4


def insertSuffix(filename, suffix, newExt=None):
    names = os.path.splitext(filename)
    if newExt is None:
        return names[0] + suffix + names[1]
    else:
        return names[0] + suffix + newExt

def getUniqueName(name):
    name, ext = os.path.splitext(name)
    n = 1
    while glob.glob(name + ".*"):
        if n > 1:
            name = name.replace('_{}'.format(n - 1), '_{}'.format(n))
        else:
            name = insertSuffix(name, '_{}'.format(n))
        n += 1
    return ''.join((name, ext))

def eucl_dist(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)


# Copyright (C) 2023-2026 Jonatan Alvelid
#
# EtMINFLUX is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EtMINFLUX is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
