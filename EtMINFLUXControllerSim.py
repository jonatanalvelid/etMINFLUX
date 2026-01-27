import os
import glob
import sys
import importlib
import enum
import warnings
import time
import math
import json

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
from tkinter.filedialog import askopenfilename

warnings.filterwarnings("ignore")


class EtMINFLUXControllerSim(QtCore.QObject):
    """ Controller part of the View-Controller-based etMINFLUX smart microscopy software. Interacts with an open Imspector software and controls the etMINFLUX experiments.
    In order to use this, you need to have Imspector running before launching this software. The interaction with Imspector uses the specpy package, a Python wrapper for 
    the Imspector API. If Imspector is not running, a mock Imspector connection will be created, which allows to test the GUI. """

    helpPlotDetectedCoordsSignal = QtCore.Signal(object, object)

    def __init__(self, widget,  *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._widget = widget
        
        print('Initializing etMINFLUX controller, simulation mode')
        
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
        self.__eventViewHelper = EventWidgetHelper(self, self._widget.eventViewWidget)

        # Connect EtMINFLUXWidget button and check box signals
        self._widget.initiateButton.clicked.connect(self.initiate)
        self._widget.loadPipelineButton.clicked.connect(self.loadPipeline)
        self._widget.softResetButton.clicked.connect(self.softReset)

        self._widget.coordListWidget.delROIButton.clicked.connect(self.deleteROI)
        self._widget.confocalFramePauseCheck.clicked.connect(self.toggleConfocalFramePause)

        self._widget.selectDataLoadButton.clicked.connect(self.getConfocalData)
        self._widget.coordTransformWidget.setSaveDirButton.clicked.connect(self.getSaveFolder)
        
        self._widget.openGuideButton.clicked.connect(self.openGuide)

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

    def initiateFlagsParams(self):
        # initiate flags and params
        self.__running = False  # run flag
        self.__runMode = RunMode.TestVisualize  # run mode currently used
        self.__validating = False  # validation flag
        self.__busy = False  # running pipeline busy flag
        self.__presetROISize = True
        self.__plotROI = False
        self.__confocalFramePause = False
        self.__count_conf_channels = 1
        self.__img_ana = None
        self._preloaded_confocal_data = deque()
        self.__prev_event_coords_deque = deque(maxlen=100)
        self.__prevFrames = deque(maxlen=50)  # deque for previous fast frames
        self.__prevAnaFrames = deque(maxlen=50)  # deque for previous preprocessed analysis frames
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

    def getConfocalData(self):
        confocal_data_file = askopenfilename(title='Select folder...')
        if confocal_data_file:
            confocal_data = tiff.imread(confocal_data_file)
        else:
            confocal_data = None
        self._preloaded_confocal_data = deque()
        for img in confocal_data:
            self._preloaded_confocal_data.append(img)
        self._preloaded_confocal_data.append(np.zeros_like(confocal_data[0]))  # add an empty image at the end, for resetting analysis
        self._widget.setConfocalDataField(confocal_data_file)

    def initiatePlottingParams(self):
        self.__colors = deque(['g','g','g','g','g'])

    def popColor(self):
        col = self.__colors.popleft()
        self.__colors.append(col)
        return col

    def initiate(self):
        """ Initiate or stop an etMINFLUX experiment. """
        if not self.__running:
            # read param for triggering random ROI from binary masks
            self.__random_roi_bin = self._widget.triggerRandomROICheck.isChecked()
            # read params for analysis pipeline from GUI
            self.__pipeline_param_vals = self.readPipelineParams()
            # reset general run parameters
            self.resetRunParams()
            self.clearConfocalData()
            # Reset parameter for extra information that pipelines can input and output
            self.__exinfo = None

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
            # check if visualization or validation mode
            if self.__runMode == RunMode.TestValidate or self.__runMode == RunMode.TestVisualize or self.__plotROI:
                self.launchHelpWidget()
            self.resetEventViewWidget()
            if self.__runMode == RunMode.Experiment:
                self.launchEventViewWidget()
            # read confocal frame time
            self.confocalFrameTime = int(float(self._widget.preload_confocal_frametime_edit.text()) * 1000) # time in ms
            self.pausetime = self.confocalFrameTime
            # start confocal imaging loop
            self.__running = True
            self.startConfocalScanning()
            self._widget.initiateButton.setText('Stop')
        else:
            # save log and confocal data
            filename_prefix = datetime.now().strftime('%y%m%d-%H%M%S')
            if self.__runMode == RunMode.Experiment:
                self.saveValidationImages(prev=True, prev_ana=True, path_prefix=filename_prefix)
            # disconnect signals and reset parameters
            self._widget.initiateButton.setText('Initiate')
            self.resetPipelineParamVals()
            self.resetRunParams()

    def clearConfocalData(self):
        """ Clear confocal data deques. """
        self.__prevFrames.clear()
        self.__prevAnaFrames.clear()

    def confocalIntermittentPause(self):
        """ Pause confocal scanning for a certain time in a confocal timelapse without continuous acquisition. """
        print('Running status in confocalIntermittentPause:', self.__running)
        if self.__running:
            if self.confPauseTimer.isActive():
                self.confPauseTimer.stop()
                self.confPauseGUICountdownTimer.stop()
            print('mid confocalIntermittentPause')
            print(self.confPauseTimerThread.isRunning())
            self.confPauseTimerThread.quit()
            self.confPauseGUICountdownTimerThread.quit()
            self.startConfPauseTimer()
            print('Finish confocalIntermittentPause')
        else:
            self.resetTimerThreads()

    def startConfocalScanning(self):
        """ Start confocal pulling of images. """
        print('Start startconfocalscanning')
        self.analysisPeriodTrigger()
        self.bufferLatestImages()
        self.confocalIntermittentPause()
        print('Finish startconfocalscanning')
        print(self.__running)

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

    def saveDetLog(self, filename_prefix):
        """ Save the event detection log to a .txt file"""
        if self.__detLog:
            self.setDetLogLine("pipeline", self.getPipelineName())
            self.logPipelineParamVals()
            self.setDetLogLine("pipeline_runtimes", self.__pipeline_runtimes)
            if self.__confocalFramePause:
                self.setDetLogLine('confocal_frame_interval', float(self.pausetime / 1000))  # time in s
            self.setDetLogLine('experiment_end', None)
            # save log file with temporal info of trigger event
            name = os.path.join(self._dataDir, filename_prefix) + '_log'
            savename = getUniqueName(name)
            log = [f'{key}: {self.__detLog[key]}' for key in self.__detLog]
            with open(f'{savename}.txt', 'w') as f:
                [f.write(f'{st}\n') for st in log]
            self.resetDetLog()

    def getTransformName(self):
        """ Get the name of the transform currently used. """
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
        self.__img_ana = None
        self.__fast_frame = 0
        self.__post_event_frames = 0
        self.__pipeline_runtimes = []
        self.__roi_events = deque(maxlen=0)
        self.__roi_sizes = deque(maxlen=0)
        self.__roiFollowMultipleCurrIdx = 0
        self.__roiFollowCurrCycle = 0
        self.__followingROIContinue = False
        self._widget.setConfGUINullMessages()
        self.resetTimerThreads()
        
    def resetPrevEventsDeque(self):
        """ Reset deque where previous event coords are saved, in endless. """
        self.__prev_event_coords_deque = deque(maxlen=100)

    def resetTimerThreads(self):
        if self.confPauseTimer.isActive():
            self.confPauseTimer.stop()
        self.confPauseTimerThread.quit()
        if self.confPauseGUICountdownTimer.isActive():
            self.confPauseGUICountdownTimer.stop()
        self.confPauseGUICountdownTimerThread.quit()

    def generateRandomCoord(self, roi_sizes):
        """ Generate random coord inside the binary mask. """         
        if self.binary_mask is not None:
            possible_pixels = np.where(self.binary_mask==True)
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
            self.imgs = []
            currimg = self._preloaded_confocal_data.popleft()
            self._preloaded_confocal_data.append(currimg)
            if np.sum(currimg) == 0:
                self.setBusyFalse()
                self.initiate()
                return analysisSuccess
            self.imgs.append(currimg)
            
            # TODO: FIX for multiple channels
            #if self.__count_conf_channels >= 2:
            #    # if dual color confocal scan - find ch2 in imspector measurement data stacks (same size as ch1)
            #    ch1_sh = np.shape(self.imgs[0])
            #    for i in range(stackidx+1,10):
            #        try:
            #            img = np.squeeze(meas.stack(i).data()[0])
            #            if np.shape(img) == ch1_sh:
            #                self.imgs.append(img)
            #                ch2_idx = i
            #                break
            #        except:
            #            pass
            #if self.__count_conf_channels == 3:
            #    # if triple color confocal scan - find ch3 in imspector measurement data stacks (same size as ch1)
            #    for i in range(ch2_idx+1,10):
            #        try:
            #            img = np.squeeze(meas.stack(i).data()[0])
            #            if np.shape(img) == ch1_sh:
            #                self.imgs.append(img)
            #                break
            #        except:
            #            pass
            
            # if initial settling frames have passed
            coords_detected, roi_sizes = self.runPipeline()
            # if visualization mode
            if self.__runMode == RunMode.TestVisualize:
                self.postPipelineVisualize(coords_detected, roi_sizes)
            # if validation mode
            elif self.__runMode == RunMode.TestValidate:
                self.postPipelineValidate(coords_detected, roi_sizes)
            # if experiment mode
            elif self.__runMode == RunMode.Experiment:
                coords_detected, coords_scan, roi_sizes = self.postPipelineExperiment(coords_detected, roi_sizes)
                # move on with triggering modality switch, if we have detected or randomized coords, and if detected coord is sufficiently far away from previous events
                if len(coords_scan)>0:
                    # check if event coords are close to coords in previous event deques
                    if len(self.__prev_event_coords_deque) > 0:
                        dists = [np.sqrt((c_old[0]-coords_scan[0])**2+(c_old[1]-coords_scan[1])**2) for c_old in self.__prev_event_coords_deque]
                        if np.min(dists) > self.setupInfo.get('analysis_settings').get('min_distance_prev_events'):
                            new_event = True
                        else:
                            new_event = False
                    else:
                        new_event = True
                    if new_event or 'fake' in self.getPipelineName():
                        self.__prev_event_coords_deque.append(coords_scan)
                        self.helpImageGeneration(coords_detected, roi_sizes)
                        if self.__prevFrames:
                            self.initiateEventViewWidget(coords_scan)  # add all prevframes
                            self.updateEventViewWidget(self.imgs[0].copy())  # update with detected frame
                        else:
                            self.initiateEventViewWidget(coords_scan, self.imgs[0].copy())  # update with detected frame
                        analysisSuccess = True
            # unset busy flag
            self.setBusyFalse()
        else:
            print('Skipped analysis period trigger, still busy from last trigger. Use longer confocal frame time.')
            self.setBusyFalse()
            self.initiate()
            return analysisSuccess
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
        coords_detected, roi_sizes, self.__exinfo, self.__img_ana = self.pipeline(*self.imgs, self.__prevFrames, self.binary_mask,
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
            # take first listed detected coord (and roi_size if applicable), if multiple, as event
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
        return coords_detected, coords_scan, roi_size

    def helpImageGeneration(self, coords_detected, roi_sizes):
        """ Called after starting MINFLUX acquisition, if some additional info regarding detected coordinates should be viewed/saved. """
        time.sleep(self.setupInfo.get('timing_settings').get('sleep_time_base'))
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
                del self.__roi_events[roi_idx-1]  # i-1 as we have already popped the first item that we are currently scanning from this deque
                try:
                    del self.__roi_sizes[roi_idx-1]  # i-1 as we have already popped the first item that we are currently scanning from this deque
                except:
                    pass
            self._widget.analysisHelpWidget.removeROI(roi_idx)

    def bufferLatestImages(self):
        # buffer latest fast frame and (if applicable) validation images
        try:
            self.__prevFrames.append(np.copy(self.imgs[0]))
            # buffer previous preprocessed analysis frame
            if self.__img_ana is not None:
                self.__prevAnaFrames.append(np.copy(self.__img_ana))
        except:
            pass

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

    def startRecTimer(self, deadtime=True, time=None):
        if deadtime:
            self.rec_time_scan = int((self.__rec_time_scan + self.__rec_time_deadtime) * 1000)  # time in ms
        else:
            self.rec_time_scan = int(time * 1000)  # time in ms
        self.recTimeTimer.setInterval(self.rec_time_scan)
        self.recTimeTimerThread.start()

    def startConfPauseTimer(self):
        print('Start confocal pause timer for ', self.pausetime/1000, 's')
        self.confPauseTimer.setInterval(self.pausetime)
        self.confPauseTimerThread.start()
        # for display countdown timer
        print('mid startConfPauseTimer')
        self.confPauseGUICountdownTimer.setInterval(1)
        self.confPauseGUICountdownTimerThread.start()
        self.confIntervalDispTimer.start()
        print('Finish startConfPauseTimer')

    def updateIntervalTimerDisp(self):
        time_left = (self.pausetime - self.confIntervalDispTimer.elapsed())/1000
        if time_left >= 0:
            self._widget.conf_guipausetimer_edit.setText(f'Time until next confocal: {time_left:.0f} s')

    def updateConfocalFrameDisp(self):
        self._widget.conf_frame_edit.setText(f'Confocal frames acquired: {self.__fast_frame+1}')

    def saveValidationImages(self, prev=True, prev_ana=True, path_prefix='YMD-HMS'):
        """ Save the validation fast images of an event detection, fast images and/or preprocessed analysis images. """
        if prev:
            self._saveImage(self.__prevFrames, path_prefix, 'conf-raw')
        if prev_ana:
            self._saveImage(self.__prevAnaFrames, path_prefix, 'conf-analysisprocessed')
        self.clearConfocalData()

    def _saveImage(self, img, path_prefix, path_suffix):
        fileName = path_prefix + '_' + path_suffix + '.tif'
        filePath = os.path.join(self._dataDir, fileName)
        tiff.imwrite(filePath, img)

    def toggleConfocalFramePause(self):
        if not self._widget.confocalFramePauseCheck.isChecked():
            self._widget.conf_frame_pause_edit.setEditable(False)
            self.__confocalFramePause = False
        else:
            self._widget.conf_frame_pause_edit.setEditable(True)
            self.__confocalFramePause = True

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
    def __init__(self, etMINFLUXController, coordTransformWidget, *args, **kwargs):

        self.etMINFLUXController = etMINFLUXController
        self._widget = coordTransformWidget

        # connect signals from widget
        self.etMINFLUXController._widget.coordTransfCalibButton.clicked.connect(self.calibrationLaunch)

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
        self._widget.recordBinaryMaskButton.setText('Recording...')
        self.binary_stack = []
        imgs = []
        for _ in range(self.etMINFLUXController.binary_frames):
            imgs.append(self.etMINFLUXController._preloaded_confocal_data.popleft())
            self.binary_stack.append(imgs[-1])
        for img in reversed(imgs):
            self.etMINFLUXController._preloaded_confocal_data.appendleft(img)
        self.calculateBinaryMask()

    def resetBinaryMask(self):
        """ Reset binary mask, back to None. """
        self.etMINFLUXController.binary_stack = []
        self.etMINFLUXController.binary_mask = None  # None to consider the whole image

    def calculateBinaryMask(self):
        """ Calculate the binary mask of the region of interest, using both positive and negative masks. """
        img_mean = np.mean(self.binary_stack, 0)
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


# Copyright (C) 2023-2025 Jonatan Alvelid
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
