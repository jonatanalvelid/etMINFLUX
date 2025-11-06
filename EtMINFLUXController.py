import os
import glob
import sys
import importlib
import enum
import warnings
import copy
import time
import math

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
from pynput.keyboard import Key, Controller

warnings.filterwarnings("ignore")


class EtMINFLUXController(QtCore.QObject):
    """ Controller part of the View-Controller-based etMINFLUX smart microscopy software. Interacts with an open Imspector software and controls the etMINFLUX experiments.
    In order to use this, you need to have Imspector running before launching this software. The interaction with Imspector uses the specpy package, a Python wrapper for 
    the Imspector API. If Imspector is not running, a mock Imspector connection will be created, which allows to test the GUI. """

    helpPlotDetectedCoordsSignal = QtCore.Signal(object, object)

    def __init__(self, widget,  *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._widget = widget
        
        print('Initializing etMINFLUX controller')
        
        ############# SYSTEM-SPECIFIC SETTINGS - MAKE SURE TO CHANGE THESE VARIABLES TO YOUR SYSTEM CONFIG #############
        # list of minflux sequences that can be triggered (INPUT THE NAMES/IDS EXACTLY AS THEY ARE IN THE SEQUENCES)
        self.mfxSeqList = ['ja_seqTrk3D_Peroxisomes_HaloJF646-m2410','ja_seqTrk3D_Dec2022-m2410base-pex13halojf646test', 'ja_seqTrk3D_Peroxisomes_HaloJF646-m2410','Tracking_2D_Fast', 'ja_seqTrk3D_Dec2022-2Csearch-m2410base','ja_seqTrk3D_Dec2022-2Cspawn-m2410base', 'ja_seqTrk3D_Dec2022-m2410base', 'ja_Tracking_2D_Fast_Nov2024-m2410base', 'ja_Tracking_2D_Nov2024-2Csearch-m2410base', 'ja_Tracking_2D_Nov2024-2Cspawn-m2410base', 'ja_Tracking_2D_Nov2024-m2410base']  # make sure that these options matches exactly those in Imspector
        # default data and transformations dirs (CREATE THESE FOLDERS IF THEY DO NOT EXIST YET, AND CHANGE THESE PATHS TO THAT OF YOUR SYSTEM)
        self._dataDir = os.path.join('C:\\Users\\Abberior_Admin\\Desktop\\Jonatan\\etMINFLUX-data', 'data')
        self._transformsDir = os.path.join('C:\\Users\\Abberior_Admin\\Desktop\\Jonatan\\etMINFLUX-data', 'transforms')
        # default screen position of Auto Repetition button in Imspector
        self._set_repeat_meas_button_pos = [1489,71]
        # default screen position of top MINFLUX dataset in Workspace widget in Imspector
        self._set_topmfxdataset_button_pos = [2157,1211]
        # list of available lasers for MFX imaging (get this list manually from Imspector control software)
        self.mfxExcLaserList = ['MINFLUX640', 'MINFLUX560']  # names for the lasers (only display names in etMINFLUX widget, the indexes below is key to selection in Imspector)
        self.mfxDetectorListView = ['DAPI', 'GFP', 'Cy3', 'Cy5 near', 'Cy5 far', 'Cy5 near & Cy5 far', 'Cy3 & Cy5 near', 'Cy3 & Cy5 far']
        self.mfxDetectorList = ['DET1', 'DET2', 'DET3', 'DET4', 'DET5', '["DET4","DET5"]', '["DET3","DET4"]', '["DET3","DET5"]']  # DET1=DAPI, DET2=GFP, DET3=Cy3, DET4=Cy5 near, DET5=Cy5 far
        self.mfxDetectorListThreads = ['["DET1"]', '["DET2"]', '["DET3"]', '["DET4"]', '["DET5"]', '["DET4","DET5"]', '["DET3","DET4"]', '["DET3","DET5"]']  # DET1=DAPI, DET2=GFP, DET3=Cy3, DET4=Cy5 near, DET5=Cy5 far
        self.laser_exc_idxs = [3,2]  # corresponding list indexes for the above lasers in the laser list in the Channels widget in Imspector 
        # name of the confocal and minflux configurations in the Imspector measurement to be used for the confocal timelapse imaging on which to run the real-time analysis pipeline, and the minflux measurements
        self.__conf_config = 'confocal'
        self.__mfx_config = 'Minflux'
        ################################################################################################################
        
        self._widget.coordTransformWidget.setSaveFolderField(self._dataDir)

        # open imspector connection
        try:
            self._imspector = specpy.get_application()
            print('Imspector connected succesfully')
        except:
            self._imspector = ImspectorMock()
            print('Load mock Imspector connection')

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
        self._widget.setMfxSequenceList(self.mfxSeqList, thread=0)
        self._widget.setMfxSequenceList(self.mfxSeqList, thread=1)

        # set list of available lasers for MFX imaging, get this list manually from Imspector control software
        self._widget.setMinfluxExcLaserList(self.mfxExcLaserList, thread=0)
        self._widget.setMinfluxExcLaserList(self.mfxExcLaserList, thread=1)
        self._widget.setMinfluxDetectorList(self.mfxDetectorListView, thread=0)
        self._widget.setMinfluxDetectorList(self.mfxDetectorListView, thread=1)

        # create a helper controller for the coordinate transform pop-out widget
        self.__coordTransformHelper = EtCoordTransformHelper(self, self._widget.coordTransformWidget, self._transformsDir)
        self.__analysisHelper = AnalysisImgHelper(self, self._widget.analysisHelpWidget)
        self.__eventViewHelper = EventWidgetHelper(self, self._widget.eventViewWidget)

        # Connect EtMINFLUXWidget button and check box signals
        self._widget.initiateButton.clicked.connect(self.initiate)
        self._widget.loadPipelineButton.clicked.connect(self.loadPipeline)
        self._widget.recordBinaryMaskButton.clicked.connect(self.initiateBinaryMask)
        self._widget.resetBinaryMaskButton.clicked.connect(self.resetBinaryMask)
        self._widget.softResetButton.clicked.connect(self.softReset)
        self._widget.saveCurrentMeasButton.clicked.connect(self.saveMeasurement)

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
        # get preset calibration values from calibration GUI
        self.getCalibrationValues()

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
        self.__plotROI = False
        self.__confocalFramePause = False
        self.__dualColor = False
        self.__mfx_twothreaded = False
        self.__img_ana = None
        self.__prev_event_coords_deque = deque(maxlen=100)
        self.__prevFrames = deque(maxlen=50)  # deque for previous fast frames
        self.__prevAnaFrames = deque(maxlen=50)  # deque for previous preprocessed analysis frames
        self.__confoffset = [0.0, 0.0]  # offset of confocal scan, in um
        self.__binary_mask = None  # binary mask of regions of interest, used by certain pipelines, leave None to consider the whole image
        self.__binary_frames = 2  # number of frames to use for calculating binary mask 
        self.__init_frames = 0  # number of frames after initiating etMINFLUX before a trigger can occur, to allow laser power settling etc
        self.__validation_frames = 2  # number of fast frames to record after detecting an event in validation mode
        self.__params_exclude = ['img', 'img_ch2', 'prev_frames', 'binary_mask', 'exinfo', 'testmode', 'presetROIsize']  # excluded pipeline parameters when loading param fields
        self.__rec_time_deadtime = 5  # deadtime when starting MINFLUX recordings, in s - specifically for m2410 that does not lock everything during this deadtime
        self.__min_dist_prev_event = 7  # minimum accepted distance to a previous event during the same recording (i.e. a run in endless mode etc)

    def getSaveFolder(self):
        self._dataDir = self._widget.coordTransformWidget.getSaveFolder()
        self._widget.coordTransformWidget.setSaveFolderField(self._dataDir)

    def initiatePlottingParams(self):
        self.__colors = deque(['g','g','g','g','g'])

    def popColor(self):
        col = self.__colors.popleft()
        self.__colors.append(col)
        return col

    def getCalibrationValues(self):
        """ Get calibration values from the calibration (coordTransform) GUI. """
        self.getTimings()
        self.__min_dist_prev_event = float(self._widget.coordTransformWidget.min_dist_prev_event_edit.text())

    def getTimings(self):
        self._sleepTime = float(self._widget.coordTransformWidget.time_sleep_edit.text())
        self._sleepTimeROISwitch = float(self._widget.coordTransformWidget.time_sleep_roiswitch_edit.text())
        self._saveTime = float(self._widget.coordTransformWidget.save_time_edit.text())

    def initiate(self):
        """ Initiate or stop an etMINFLUX experiment. """
        if not self.__running:
            # get timings from ROI input
            self.getCalibrationValues()
            # read mfx sequence and lasers and laser powers from GUI
            sequenceIdxth0 = self._widget.mfxth0_seq_par.currentIndex()
            self.mfxth0_seq = self._widget.mfx_seqs[sequenceIdxth0]
            laserIdxth0 = self._widget.mfxth0_exc_laser_par.currentIndex()
            self.mfxth0_exc_laser = self._widget.mfx_exc_lasers[laserIdxth0]
            self.mfxth0_exc_pwr = float(self._widget.mfxth0_exc_pwr_edit.text())
            detectorIdxth0 = self._widget.mfxth0_detector_par.currentIndex()
            self.mfxth0_detector = self.mfxDetectorListThreads[detectorIdxth0]
            self.mfxch0_detector = self.mfxDetectorList[detectorIdxth0]
            #self.mfxth0_act_pwr = float(self._widget.mfx_act_pwr_edit.text())
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
            # launch help widget, if visualization mode or validation mode
            # Check if visualization mode, in case launch help widget
            experimentModeIdx = self._widget.experimentModesPar.currentIndex()
            experimentMode = self._widget.experimentModes[experimentModeIdx]
            if experimentMode == 'TestVisualize':
                self.__runMode = RunMode.TestVisualize
            elif experimentMode == 'TestValidate':
                self.__runMode = RunMode.TestValidate
            else:
                self.__runMode = RunMode.Experiment
                self.__plotROI = self._widget.plotROICheck.isChecked()
            # check if visualization or validation mode
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

    def clearConfocalData(self):
        """ Clear confocal data deques. """
        self.__prevFrames.clear()
        self.__prevAnaFrames.clear()

    def finishMFXROIAuto(self):
        """ Trigger this when a preset-rec-time MFX ROI has finished. """
        self.__runningMFX = False
        self.scanEnded()

    def imspectorLineEvent(self):
        """ Called when a line in a confocal image is finished in Imspector. """
        self.__confocalLineAnalysisCurr += 1
        self.__confocalLineFrameCurr += 1
        analysisSuccess = False
        if self.__confocalLineAnalysisCurr == self.__confocalLinesAnalysisPeriod:
            self.__confocalLineAnalysisCurr = 0
            analysisSuccess = self.analysisPeriodTrigger()
        if self.__confocalLineFrameCurr >= self.__confocalLinesFrame:
            self.__confocalLineFrameCurr = 0
            self.__fast_frame += 1
            self.bufferLatestImages()
            self.updateConfocalFrameDisp()
            if self.__confocalFramePause and not analysisSuccess:
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
        self.__confocalLineCurr = 0
        # ensure activate configuration is conf overview
        self.activateConfocalConfig()
        time.sleep(self._sleepTime)
        # enable auto-looping and start scan
        self._imspector.enable_loop(True)
        self._imspector.start()

    def scanEnded(self):
        """ End a MINFLUX acquisition. """
        if self.__followingROI and self.__followingROIMode == ROIFollowMode.Multiple:
            self.setDetLogLine(f"mfx_end-id{self.__roiFollowMultipleCurrIdx}-cycle{self.__roiFollowCurrCycle}", None)
            self.setDetLogLine("recording_mode", "roi_follow_multiple") 
        elif self.__followingROI:
            self.setDetLogLine(f"mfx_end-cycle{self.__roiFollowCurrCycle}", None) 
            self.setDetLogLine("recording_mode", "roi_follow_single")               
        elif self.__run_all_aoi:
            self.setDetLogLine(f"mfx_end-id{self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque)}", None)
            self.setDetLogLine("recording_mode", "all_roi") 
        else:
            self.setDetLogLine("mfx_end", None)
            self.setDetLogLine("recording_mode", "single_roi") 
        if self.__plotROI and not self.__followingROI:
            self.deleteROIGUI(idx=0)  # delete top ROI from list
        if (not self.__run_all_aoi or not self.__aoi_coords_deque) and (not self.__followingROIContinue):
            time.sleep(self._sleepTimeROISwitch)
            self.endExperiment()
            time.sleep(self._sleepTimeROISwitch)
            self.continueFastModality()
            self.__fast_frame = 0
        elif self.__followingROI and self.__followingROIMode == ROIFollowMode.Multiple:
            time.sleep(self._sleepTimeROISwitch)  # to make sure Imspector has properly turned off MINFLUX recording
            self.__roiFollowMultipleCurrIdx += 1
            if self.__roiFollowMultipleCurrIdx == len(self.__aoi_coords_deque):
                self.__roiFollowCurrCycle += 1
                self.__roiFollowMultipleCurrIdx = 0
                self.endFollowingROIStep()
                time.sleep(self._sleepTimeROISwitch)
                self.continueFastModality()
                return
            coords = self.__aoi_coords_deque.popleft()
            self.__aoi_coords_deque.append(coords)
            self.newROIMFX(coords, roi_idx=self.__roiFollowMultipleCurrIdx)
        elif self.__followingROI:
            time.sleep(self._sleepTimeROISwitch)  # to make sure Imspector has properly turned off MINFLUX recording
            self.__roiFollowCurrCycle += 1
            self.endFollowingROIStep()
            time.sleep(self._sleepTimeROISwitch)
            self.continueFastModality()
        elif self.__followingROIContinue:
            # only triggers if followingROI check box have been unclicked to end a followingROI experiment,
            # in which case one more confocal frame will be run, and all data will afterwards be saved,
            # before resetting everytning.
            time.sleep(2*self._sleepTimeROISwitch)
            self.continueFastModality()
            self.__fast_frame = 0
        elif self.__run_all_aoi:
            time.sleep(2*self._sleepTimeROISwitch)  # to make sure Imspector has properly turned off MINFLUX recording
            self.newROIMFX(self.__aoi_coords_deque.popleft(), roi_idx=self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque))

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
            # log rec time, if exisiting
            if self.__hasRunMFX:
                if self.__presetRecTime:
                    self.setDetLogLine("MFXRecTime", self.__rec_time_scan)
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
        self._widget.initParamFields(self.__pipeline_params, self.__params_exclude)
        if 'dualcolor' in pipelinename:
            self.__dualColor = True
        else:
            self.__dualColor = False

    def initiateBinaryMask(self):
        """ Initiate the process of calculating a binary mask of the region of interest. """
        self.launchHelpWidget()
        self.__binary_stack = []
        self.__binary_frame = 0
        self.__confocalLinesAnalysisPeriod = self._imspector.active_measurement().stack(0).sizes()[0] - 7
        self.__confocalLineCurr = 0
        self._imspector.connect_end(self.imspectorLineEventBinaryMask, 1)
        self.startConfocalScanning(reconnect_signal=False)
        self._widget.recordBinaryMaskButton.setText('Recording...')

    def resetBinaryMask(self):
        """ Reset binary mask, back to None. """
        self.__binary_stack = []
        self.__binary_mask = None  # None to consider the whole image

    def imspectorLineEventBinaryMask(self):
        self.__confocalLineCurr += 1
        if self.__confocalLineCurr == self.__confocalLinesAnalysisPeriod:
            self.__confocalLineCurr = 0
            # get image
            meas = self._imspector.active_measurement()
            self.__binary_stack.append(np.squeeze(meas.stack(0).data()))
            self.__binary_frame += 1
            if self.__binary_frame == self.__binary_frames:
                self.calculateBinaryMask()

    def calculateBinaryMask(self):
        """ Calculate the binary mask of the region of interest, using both positive and negative masks. """
        self._imspector.pause()
        img_mean = np.mean(self.__binary_stack, 0)
        img_bin_pos = ndi.filters.gaussian_filter(img_mean, float(self._widget.bin_smooth_edit.text()))
        img_bin_neg = ndi.filters.gaussian_filter(img_mean, float(self._widget.bin_neg_smooth_edit.text()))
        img_bin_pos = np.array(img_bin_pos > float(self._widget.bin_thresh_edit.text()))
        img_bin_neg = np.array(img_bin_neg < float(self._widget.bin_neg_thresh_edit.text()))
        img_bin_border = np.zeros(np.shape(img_bin_pos))
        border_size = int(self._widget.bin_border_size_edit.text())
        img_bin_border[border_size:-border_size, border_size:-border_size] = 1
        self.__binary_mask = np.logical_and(np.logical_and(img_bin_pos, img_bin_neg), img_bin_border)
        self._widget.recordBinaryMaskButton.setText('Record binary mask')
        self.setAnalysisHelpImg(self.__binary_mask)
        self._imspector.disconnect_end(self.imspectorLineEventBinaryMask, 1)

    def setAnalysisHelpImg(self, img):
        """ Set the preprocessed image in the analysis help widget. """
        if self.__fast_frame < self.__init_frames + 3:
            autolevels = True
        else:
            autolevels = False
        self._widget.analysisHelpWidget.img.setImage(img, autoLevels=autolevels)
        infotext = f'Min: {np.floor(np.min(img))}, max: {np.floor(np.max(img))}'
        self._widget.analysisHelpWidget.info_label.setText(infotext)

    def softReset(self):
        """ Reset all parameters after a soft lock. """
        self.setBusyFalse()
        self.resetHelpWidget()
        self.resetEventViewWidget()
        self.resetRunParams()
        self.initiateFlagsParams()
        self.resetDetLog()

    def setBusyFalse(self):
        """ Set busy flag to false. """
        self.__busy = False

    def resetHelpWidget(self):
        self._widget.resetHelpWidget()
        self.__analysisHelper = AnalysisImgHelper(self, self._widget.analysisHelpWidget)

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
        self._widget.launchHelpWidget(self._widget.analysisHelpWidget, init=True)

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
        self.__aoi_coords_deque = deque(maxlen=0)
        self.__aoi_sizes_deque = deque(maxlen=0)
        self.__roiFollowMultipleCurrIdx = 0
        self.__roiFollowCurrCycle = 0
        self.__followingROIContinue = False
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
        # enter if not still running last analysis period trigger
        if not self.__busy:
            # set busy true
            self.__busy = True
            # get image
            meas = self._imspector.active_measurement()
            stackidx = 3
            self.img = np.squeeze(meas.stack(stackidx).data()[0])
            if self.__dualColor:
                # if dual color confocal scan - find ch2 in imspector measurement data stacks (same size as ch1)
                ch1_sh = np.shape(self.img)
                for i in range(stackidx+1,10):
                    try:
                        img = np.squeeze(meas.stack(i).data()[0])
                        if np.shape(img) == ch1_sh:
                            self.img_ch2 = img
                            break
                    except:
                        pass
            # if recording mode is following ROI and we are currently following a ROI
            if self.__followingROIContinue:
                analysisSuccess = True
                if not self.__followingROI:
                    # if following ROI mode check box has been unchecked - end experiment
                    self.endFollowingROIExperiment()
                elif self.__followingROIMode == ROIFollowMode.SingleRedetect:
                    # if we are redetecting a single ROI
                    coords_detected, roi_sizes = self.runPipeline()
                    coords_scan, roi_size = self.postPipelineFollowingROISingleRedetect(coords_detected, roi_sizes)
                    if len(coords_scan)>0:
                        self.acquireMINFLUXFull(coords_scan, roi_size)
                        self.updateEventViewWidget(self.img.copy(), coords_scan)
                    else:
                        # end experiment, as we no longer have anything to follow, and we cannot suddenly "find it back" again in the next image
                        self.endFollowingROIExperiment()
                elif self.__followingROIMode == ROIFollowMode.Single:
                    self.acquireMINFLUXMinimal(pos=self.__roi_center_mfx, roi_size_um=self.__roi_size_um_mfx, pos_conf=self.__pos_conf_mfx)
                    self.updateEventViewWidget(self.img.copy())
                elif self.__followingROIMode == ROIFollowMode.Multiple:
                    # initiate new cycle of following ROI
                    coords = self.__aoi_coords_deque.popleft()
                    self.__aoi_coords_deque.append(coords)
                    #self.__roiFollowMultipleCurrIdx += 1
                    # pause fast imaging
                    self.pauseFastModality()
                    self.newROIMFX(coords, roi_idx=self.__roiFollowMultipleCurrIdx)
            else:
                # if initial settling frames have passed
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
                        coords_detected, coords_scan, roi_size = self.postPipelineExperiment(coords_detected, roi_sizes)
                        # move on with triggering modality switch, if we have detected or randomized coords, and if detected coord is sufficiently far away from previous events
                        if len(coords_scan)>0:
                            # check if event coords are close to coords in previous event deques
                            if len(self.__prev_event_coords_deque) > 0:
                                dists = [np.sqrt((c_old[0]-coords_scan[0])**2+(c_old[1]-coords_scan[1])**2) for c_old in self.__prev_event_coords_deque]
                                if np.min(dists) > self.__min_dist_prev_event:
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
                                    self.updateEventViewWidget(self.img.copy())  # update with detected frame
                                else:
                                    self.initiateEventViewWidget(coords_scan, self.img.copy())  # update with detected frame
                                analysisSuccess = True
            # unset busy flag
            self.setBusyFalse()
        return analysisSuccess

    def initiateEventViewWidget(self, event_coords, frame=None):
        """ Initiate the event view widget with the current confocal frame and event coordinates. """
        if frame is not None:
            self.__eventViewHelper.new_event(frame, event_coords)
        else:
            self.__eventViewHelper.new_event(self.__prevFrames, event_coords)

    def updateEventViewWidget(self, new_frame, event_coords=None):
        """ Update the event view widget with a new confocal frame and event coordinates. """
        self.__eventViewHelper.new_frame(new_frame, event_coords)

    def runPipeline(self):
        """ Run the analyis pipeline on the latest confocal analysis period frame. """
        # start of pipeline
        self.__pipeline_start_file_prefix = datetime.now().strftime('%Y%m%d-%Hh%Mm%Ss')  # use for naming files
        # run pipeline
        self.__pipeline_start_time = time.perf_counter()
        if self.__dualColor:
            coords_detected, roi_sizes, self.__exinfo, self.__img_ana = self.pipeline(self.img, self.img_ch2, self.__prevFrames, self.__binary_mask,
                                                                        self.__exinfo, self.__presetROISize,
                                                                        *self.__pipeline_param_vals)
        else:
            coords_detected, roi_sizes, self.__exinfo, self.__img_ana = self.pipeline(self.img, self.__prevFrames, self.__binary_mask,
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
            self.helpPlotDetectedCoordsSignal.emit(coords_detected, roi_sizes)

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
            self.helpPlotDetectedCoordsSignal.emit(coords_detected, roi_sizes)
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

    def postPipelineExperiment(self, coords_detected, roi_sizes):
        """ Called if in experiment mode, and returns MINFLUX scan coordinates according to GUI settings, detected or random. """
        # if to generate random roi(s) from binary
        if self.__random_roi_bin:
            coords_detected, coords_scan, roi_size = self.generateRandomCoord(roi_sizes)
        # else handle detected coordinates, if some events were detected
        elif coords_detected.size != 0:
            pipelinename = self.getPipelineName()
            # if we want to run all detected ROIs in following mode
            if self.__followingROI and self.__followingROIMode == ROIFollowMode.Multiple:
                # prep circular deques for coords and ROI sizes
                areas_of_interest = list()
                if np.size(coords_detected) > 2:
                    for pair in coords_detected:
                        areas_of_interest.append(pair)
                else:
                    areas_of_interest.append(coords_detected[0])
                self.__aoi_coords_deque = deque(areas_of_interest, maxlen=len(areas_of_interest))
                coords_scan = self.__aoi_coords_deque.popleft()
                self.__aoi_coords_deque.append(coords_scan)
                # if we do not preset ROI size, also make a deque for the detected roi_sizes
                if not self.__presetROISize:
                    self.__aoi_sizes_deque = deque(roi_sizes, maxlen=len(areas_of_interest))
                    roi_size = self.__aoi_sizes_deque.popleft()
                    self.__aoi_sizes_deque.append(roi_size)
                else:
                    roi_size = None
                self._widget.initiateButton.setText('Next ROI')
            # if we want to run all detected ROIs once
            elif self.__run_all_aoi:
                # prep deques for coords and ROI sizes
                areas_of_interest = list()
                if np.size(coords_detected) > 2:
                    for pair in coords_detected:
                        areas_of_interest.append(pair)
                else:
                    areas_of_interest.append(coords_detected[0])
                self.__aoi_coords_deque = deque(areas_of_interest, maxlen=len(areas_of_interest))
                coords_scan = self.__aoi_coords_deque.popleft()
                # if we do not preset ROI size, also make a deque for the detected roi_sizes
                if not self.__presetROISize:
                    self.__aoi_sizes_deque = deque(roi_sizes, maxlen=len(areas_of_interest))
                    roi_size = self.__aoi_sizes_deque.popleft()
                else:
                    roi_size = None
                self._widget.initiateButton.setText('Next ROI')
            elif self.__exinfo is not None and (pipelinename=='peak_detection_def' or pipelinename=='peak_detection_rand'):
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
                # take brightest detected coord (and roi_size if applicable) as event
                idx = 0
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
        self.__coords_scan_curr = coords_scan
        return coords_detected, coords_scan, roi_size

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
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, idx=self.__roiFollowMultipleCurrIdx, cycle=self.__roiFollowCurrCycle)
                self.setDetLogLine(f"mfx_initiate-id{self.__roiFollowMultipleCurrIdx}-cycle{self.__roiFollowCurrCycle}", None)
            elif self.__followingROI:
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, cycle=self.__roiFollowCurrCycle)
                self.setDetLogLine(f"mfx_initiate-cycle{self.__roiFollowCurrCycle}", None)
            elif self.__run_all_aoi:
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, idx=self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque))
                self.setDetLogLine(f"mfx_initiate-id{self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque)}", None)
            else:
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan)
                self.setDetLogLine("mfx_initiate", None)
        self.startMFX()

    def endFollowingROIExperiment(self):
        """ End experiment after a recalculated following ROI was not detected. """
        self.pauseFastModality()
        self.endExperiment()
        self.__followingROIContinue = False
        self.continueFastModality()

    def helpImageGeneration(self, coords_detected, roi_sizes):
        """ Called after starting MINFLUX acquisition, if some additional info regarding detected coordinates should be viewed/saved. """
        time.sleep(self._sleepTime)
        if self.__plotROI:
            # set analysis image in help widget and plot detected coords in help widget
            self.helpPlotDetectedCoordsSignal.emit(coords_detected, roi_sizes)
            # set analysis image in help widget
            self.setAnalysisHelpImg(self.__img_ana)

    def deleteROI(self):
        roi_idx = self._widget.coordListWidget.list.currentRow()
        if roi_idx > -1:
            self._widget.coordListWidget.list.takeItem(roi_idx)
            if self.__runMode == RunMode.Experiment:
                del self.__aoi_coords_deque[roi_idx-1]  # i-1 as we have already popped the first item that we are currently scanning from this deque
                try:
                    del self.__aoi_sizes_deque[roi_idx-1]  # i-1 as we have already popped the first item that we are currently scanning from this deque
                except:
                    print('roi sizes empty - preset size')
            self._widget.analysisHelpWidget.removeROI(roi_idx)

    def activateConfocalConfig(self):
        meas = self._imspector.active_measurement()
        cfg = meas.configuration(self.__conf_config)
        meas.activate(cfg)
    
    def activateMinfluxConfig(self):
        meas = self._imspector.active_measurement()
        cfg = meas.configuration(self.__mfx_config)
        meas.activate(cfg)

    def bufferLatestImages(self):
        # buffer latest fast frame and (if applicable) validation images
        self.__prevFrames.append(np.copy(self.img))
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
            roi_size = self.__aoi_sizes_deque.popleft()
            if self.__followingROI and self.__followingROIMode == ROIFollowMode.Multiple:
                self.__aoi_sizes_deque.append(roi_size)
            coord_um_top = self.transform([coords_scan[0]+roi_size[0]/2, coords_scan[1]+roi_size[1]/2], self.__confoffset, self.__transformCoeffs)
            coord_um_bot = self.transform([coords_scan[0]-roi_size[0]/2, coords_scan[1]-roi_size[1]/2], self.__confoffset, self.__transformCoeffs)
            roi_size_um_scan = np.subtract(coord_um_top, coord_um_bot)
        time.sleep(self._sleepTime)
        # get preset rec time, if using that mode
        if self.__presetRecTime:
            rec_time_scan = float(self._widget.mfx_rectime_edit.text())
            self.__rec_time_scan = rec_time_scan
        else:
            self.__rec_time_scan = None
        # initiate and run scanning with transformed center coordinate
        self.initiateMFX(position=coords_center_um, ROI_size=roi_size_um_scan, pos_conf=coords_scan)
        # log detected and scanning center coordinate
        if self.__followingROI:  # and self.__followingROIMode == ROIFollowMode.Multiple:
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
        if self.__followingROI: # and self.__followingROIMode == ROIFollowMode.Multiple:
            self.setMFXDataTag(pos_conf, ROI_size, self.__roiFollowMultipleCurrIdx, self.__roiFollowCurrCycle)
        else:
            self.setMFXDataTag(pos_conf, ROI_size, self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque))
        self.setMFXLasers(self.mfxth0_exc_laser, self.mfxth0_exc_pwr, thread=0)  #, self.mfx_act_pwr)
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
        time.sleep(self._sleepTime)

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
        self.confPauseGUICountdownTimer.setInterval(1)
        self.confPauseGUICountdownTimerThread.start()
        self.confIntervalDispTimer.start()

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
            time.sleep(self._sleepTime/3)
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

    def setMFXLasers(self, exc_laser, exc_pwr, thread=0):  #, act_pwr):
        """ Sets MINFLUX lasers and laser powers, according to the GUI choice of the user. """
        self._imspector.value_at('Minflux/threads/settings/'+str(thread)+'/exc', specpy.ValueTree.Measurement).set(exc_laser)
        self._imspector.value_at('Minflux/threads/settings/'+str(thread)+'/exp', specpy.ValueTree.Measurement).set(exc_pwr)
        # set activation laser (405)
        #laser_status_act = True if act_pwr > 0 else False
        #self._imspector.value_at('ExpControl/measurement/channels/0/lasers/0/active', specpy.ValueTree.Measurement).set(laser_status_act)
        #self._imspector.value_at('ExpControl/measurement/channels/0/lasers/0/power/calibrated', specpy.ValueTree.Measurement).set(act_pwr)

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
        time.sleep(self._sleepTime/6)
        self.activateMinfluxConfig()
        self._imspector.value_at('ExpControl/scan/range/x/len', specpy.ValueTree.Measurement).set(ROI_size[0]*1e-6)  # input units, m
        self._imspector.value_at('ExpControl/scan/range/y/len', specpy.ValueTree.Measurement).set(ROI_size[1]*1e-6)  # input units, m
        self._imspector.value_at('ExpControl/scan/range/x/off', specpy.ValueTree.Measurement).set(position[0]*1e-6)  # input units, m
        self._imspector.value_at('ExpControl/scan/range/y/off', specpy.ValueTree.Measurement).set(position[1]*1e-6)  # input units, m

    def setMFXROI_legacy(self, position, ROI_size):
        """ Set the MINFLUX ROI by mouse control: drag ROI, and click "Set as MFX ROI"-button. Change button-click with mouse to shortcut click with keyboard emulation. """
        # get monitor px position
        shift = 1
        if self.__presetROISize:
            positions = (position[0] - ROI_size[0]/self.__px_size_mon[0]/2 + shift,
                        position[1] - ROI_size[1]/self.__px_size_mon[1]/2 + shift,
                        position[0] + ROI_size[0]/self.__px_size_mon[0]/2 + shift,
                        position[1] + ROI_size[1]/self.__px_size_mon[1]/2 + shift)
        else:
            positions = (position[0] - ROI_size[0]/2 + shift,
                        position[1] - ROI_size[1]/2 + shift,
                        position[0] + ROI_size[0]/2 + shift,
                        position[1] + ROI_size[1]/2 + shift)
        # click two positions to remove any currently drawn ROIs
        border_margin = 5
        mock_pos1 = [self.__transformCoeffs[0]+border_margin, self.__transformCoeffs[1]+border_margin]  # top left
        mock_pos2 = [mock_pos1[0]+self.__transformCoeffs[2]-border_margin*2, mock_pos1[1]+self.__transformCoeffs[3]-border_margin*2]  # bottom right
        mouse.move(*mock_pos1)
        mouse.click()
        mouse.move(*mock_pos2)
        mouse.click()
        # drag actual ROI
        mouse.drag(*positions, absolute=True, duration=self._mouse_drag_duration)
        # keyboard shortcut clicks to set as MFX ROI
        self.keyboard.press(Key.ctrl)
        self.keyboard.press(Key.shift)
        self.keyboard.press(Key.alt)
        self.keyboard.press('m')
        self.keyboard.release('m')
        self.keyboard.release(Key.alt)
        self.keyboard.release(Key.shift)
        self.keyboard.release(Key.ctrl)

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
            time.sleep(self._saveTime)
            self.deleteMFXDataset(times=10)  # deletes 10 datasets, which should always be enough, no simple way to detect how many datasets there are

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
        else:
            self._widget.mfx_rectime_edit.setEditable(False)
            self._widget.presetMfxRecTimeCheck.setEnabled(False)
            self.__followingROI = True

    def toggleTriggerAllROIs(self):
        if not self._widget.triggerAllROIsCheck.isChecked():
            self.__run_all_aoi = False
        else:
            self.__run_all_aoi = True

    def getPresetRoiSize(self, length):
        roi_sizes = []
        for _ in range(length):
            roi_sizes.append([float(self._widget.size_x_edit.text()), float(self._widget.size_y_edit.text())])
        return roi_sizes


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


class EventWidgetHelper():
    """ Event widget helper, with confocal image viewing and intensity plot; help controller. """
    def __init__(self, etMINFLUXController, eventWidget, *args, **kwargs):
        self.etMINFLUXController = etMINFLUXController
        self._widget = eventWidget
        self.zoom_size = 10
        self.intensity_size = 2
        self.event_coords = []
        self.intensity_trace = []

        # connect signals from widget
        self._widget.image_viewer.setLevelsButton.clicked.connect(self.setLevels)

    def new_event(self, frames, coords):
        if len(np.shape(frames))==3:
            frames_arr = np.array(frames)
        else:
            frames_arr = np.array([frames])
        zoomstack = frames_arr[:, int(coords[1]-self.zoom_size):int(coords[1]+self.zoom_size+1), int(coords[0]-self.zoom_size):int(coords[0]+self.zoom_size+1)]
        intensity_trace = np.mean(zoomstack[:, int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1), int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1)], axis=(1,2))
        self.intensity_trace = intensity_trace
        self._widget.add_event(zoomstack, self.intensity_trace)
        self.event_coords = coords

    def new_frame(self, frame, coords):
        if coords == None:
            coords = self.event_coords
        zoomframe = frame[int(coords[1]-self.zoom_size):int(coords[1]+self.zoom_size+1), int(coords[0]-self.zoom_size):int(coords[0]+self.zoom_size+1)]
        newintensity = np.mean(zoomframe[int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1), int(self.zoom_size-self.intensity_size):int(self.zoom_size+self.intensity_size+1)], axis=(0,1))
        self.intensity_trace = np.append(self.intensity_trace, newintensity)
        self._widget.add_frame_and_intensity(zoomframe, self.intensity_trace)

    def setLevels(self):
        """ Set zoom image min,max levels. """
        min_val = float(self._widget.image_viewer.levelMinEdit.text())
        max_val = float(self._widget.image_viewer.levelMaxEdit.text())
        self._widget.image_viewer.img.setLevels([min_val, max_val])


class EtCoordTransformHelper():
    """ Coordinate transform help widget controller. """
    def __init__(self, etMINFLUXController, coordTransformWidget, saveFolder, *args, **kwargs):

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
        time.sleep(self.etMINFLUXController._sleepTime)
        self._widget.setDeleteMFXDatasetButtonText(self._set_topmfxdataset_button_pos)

    def setTopMFXDatasetPos(self):
        mouse_pos = mouse.get_position()
        self._set_topmfxdataset_button_pos = np.array(mouse_pos)
        mouse.unhook_all()

    def calibrationLaunch(self):
        """ Launch calibration. """
        self.etMINFLUXController._widget.launchHelpWidget(self.etMINFLUXController._widget.coordTransformWidget, init=True)


class RunMode(enum.Enum):
    Experiment = 1
    TestVisualize = 2
    TestValidate = 3


class ROIFollowMode(enum.Enum):
    Single = 1
    Multiple = 2
    SingleRedetect = 3


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


class ImspectorMock():
    """ Mock Imspector object, that handles the calls made from the EtMINFLUXController object above. """
    def __init__(self, *args, **kwargs):
        self._active_measurement = Measurement()

    def active_measurement(self):
        return self._active_measurement
    
    def value_at(self, path, valTree):
        return ValueAt()

    def start(self):
        pass

    def run(self):
        pass

    def pause(self):
        pass


class Measurement():
    def __init__(self, *args, **kwargs):
        self._stack = Stack()

    def stack(self, val):
        return self._stack.data()
    
    def configuration(self, cfg_name):
        return Configuration(cfg_name)

    def activate(self, cfg):
        self._configuration = cfg


class Configuration():
    def __init__(self, name, *args, **kwargs):
        self._name = name


class Stack():
    def __init__(self, *args, **kwargs):
        self._data = Data(3,100,100)

    def data(self):
        return self._data

class Data():
    def __init__(self, z, x, y, *args, **kwargs):
        self._data = np.zeros((z,x,y))
        self._sizes = (x,y)

    def data(self):
        return self._data
    
    def sizes(self):
        return self._sizes    

class ValueAt():
    def __init__(self, *args, **kwargs):
        self._val = False

    def get(self):
        self._val

    def set(self, val):
        self._val = val
    
    def trigger(self):
        pass


# Copyright (C) 2023-2025 Jonatan Alvelid
#
# EtMINFLUX is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ImSwitch is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
