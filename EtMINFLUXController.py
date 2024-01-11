import os
import glob
import sys
import importlib
import enum
import warnings
import time
import specpy
import mouse
import threading
from collections import deque
from datetime import datetime
from inspect import signature
from qtpy import QtCore
#from tkinter import Tk, filedialog

import tifffile as tiff
import scipy.ndimage as ndi
import numpy as np
import pyqtgraph as pg
#from scipy.optimize import least_squares

warnings.filterwarnings("ignore")

#_logsDir = os.path.join('C:\\Users\\Abberior_Admin\\Documents\\Jonatan\\etminflux-files', 'recordings', 'logs')
#_dataDir = os.path.join('C:\\Users\\Abberior_Admin\\Documents\\Jonatan\\etminflux-files', 'recordings', 'data')
#_transformsDir = os.path.join('C:\\Users\\Abberior_Admin\\Documents\\Jonatan\\etminflux-files', 'recordings', 'transforms')
_logsDir = os.path.join('C:\\Users\\alvelidjonatan\\Documents\\Data\\etMINFLUX', 'recordings', 'logs')
_dataDir = os.path.join('C:\\Users\\alvelidjonatan\\Documents\\Data\\etMINFLUX', 'recordings', 'data')
_transformsDir = os.path.join('C:\\Users\\alvelidjonatan\\Documents\\Data\\etMINFLUX', 'recordings', 'transforms')

def thread_info(msg):
    print(msg, int(QtCore.QThread.currentThreadId()), threading.current_thread().name)

class EtMINFLUXController(QtCore.QObject):
    """ Linked to EtMINFLUXWidget."""

    def __init__(self, widget,  *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._widget = widget
        
        print('Initializing etMINFLUX controller')

        # open imspector connection
        try:
            self._imspector = specpy.get_application()
            print('Imspector connected succesfully')
        except:
            self._imspector = ImspectorMock()
            print('Load mock Imspector connection')

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

        ## list of trigger modalities
        #self.triggerModalityList = ['Confocal', 'Confocal fast']
        #self._widget.setTriggerModalityList(self.triggerModalityList)

        # list of minflux sequences that can be triggered
        self.mfxSeqList = ['Imaging_2D', 'Imaging_3D', 'Tracking_2D', 'Tracking_2D_Fast', 'ja_triangle_dmp_lipids']  # make sure that these options matches exactly those in Imspector
        self._widget.setMfxSequenceList(self.mfxSeqList)

        # list of available lasers for MFX imaging, get this list manually from Imspector control software
        self.mfxExcLaserList = ['640', '561', '488']
        self._widget.setMinfluxExcLaserList(self.mfxExcLaserList)

        # create a helper controller for the coordinate transform pop-out widget
        self.__coordTransformHelper = EtCoordTransformHelper(self, self._widget.coordTransformWidget, _transformsDir)
        self.__analysisHelper = AnalysisImgHelper(self, self._widget.analysisHelpWidget)

        # Connect EtMINFLUXWidget button and check box signals
        self._widget.initiateButton.clicked.connect(self.initiate)
        self._widget.loadPipelineButton.clicked.connect(self.loadPipeline)
        self._widget.recordBinaryMaskButton.clicked.connect(self.initiateBinaryMask)
        self._widget.setBusyFalseButton.clicked.connect(self.setBusyFalse)
        self._widget.setMFXROICalibrationButton.clicked.connect(self.setMFXROIButtonPosButtonCall)
        self._widget.setRepeatMeasCalibrationButton.clicked.connect(self.setRepeatMeasButtonPosButtonCall)
        self._widget.saveCurrentMeasButton.clicked.connect(self.saveMeasurement)

        self._widget.presetROISizeCheck.clicked.connect(self.togglePresetROISize)
        self._widget.presetMfxRecTimeCheck.clicked.connect(self.togglePresetRecTime)
        self._widget.autoSaveCheck.clicked.connect(self.togglePresetAutoSave)
        self._widget.followROIModeCheck.clicked.connect(self.toggleFollowingROIMode)
        self._widget.followROIRedetectCheck.clicked.connect(self.toggleRedetectROIMode)

        # create timer for fixed recording time syncing
        self.timerThread = QtCore.QThread(self)
        self.recTimeTimer = QtCore.QTimer()
        self.recTimeTimer.setSingleShot(True)
        self.recTimeTimer.timeout.connect(self.initiate)
        self.recTimeTimer.moveToThread(self.timerThread)
        self.timerThread.started.connect(self.recTimeTimer.start)

        # initiate log for each detected event
        self.resetDetLog()
        # initiate pipeline parameter values
        self.resetPipelineParamVals()
        # initiate run parameters
        self.resetRunParams()
        # initiate other parameters and flags used during experiments
        self.initiateFlagsParams()

    def initiateFlagsParams(self):
        # initiate flags and params
        self.__running = False  # run flag
        self.__runningMFX = False  # run MFX flag
        self.__runMode = RunMode.Experiment  # run mode currently used
        self.__validating = False  # validation flag
        self.__busy = False  # running pipeline busy flag
        self.__presetROISize = True
        self.__presetRecTime = False
        self.__lineWiseAnalysis = False
        self.__autoSaveMeas = False
        self.__followingROIMode = False
        self.__followingROIModeContinue = False
        self.__followingROIRedetect = False
        self.__plotROI = False
        self.__conf_config = 'ov conf'
        #self.__bkg = None  # bkg image
        self.__prevFrames = deque(maxlen=10)  # deque for previous fast frames
        self.__prevAnaFrames = deque(maxlen=10)  # deque for previous preprocessed analysis frames
        self.__binary_mask = None  # binary mask of regions of interest, used by certain pipelines, leave None to consider the whole image
        self.__binary_frames = 2  # number of frames to use for calculating binary mask 
        self.__init_frames = 0  # number of frames after initiating etMINFLUX before a trigger can occur, to allow laser power settling etc
        self.__validation_frames = 2  # number of fast frames to record after detecting an event in validation mode
        self.__params_exclude = ['img', 'prev_frames', 'binary_mask', 'exinfo', 'testmode', 'presetROIsize']  # excluded pipeline parameters when loading param fields
        self.__run_all_aoi = False  # run all detected events/area flag
        self.__rec_time_deadtime = 3  # deadtime when starting MINFLUX recordings, in s

        # initiate and set parameters for automatic mouse control
        # office default
        #self.__set_MFXROI_button_pos = [652,65]
        #self.__set_repeat_meas_button_pos = [407,65]
        # lab default
        self.__set_MFXROI_button_pos = [1613,75]
        self._widget.setMFXROICalibrationButtonText(self.__set_MFXROI_button_pos)
        self.__set_repeat_meas_button_pos = [1378,75]
        self._widget.setRepeatMeasCalibrationButtonText(self.__set_repeat_meas_button_pos)

    def getTimings(self):
        self.__sleepTime = float(self._widget.time_sleep_edit.text())
        self.__sleepTimeROISwitch = float(self._widget.time_sleep_roiswitch_edit.text())
        self.__mouse_drag_duration = float(self._widget.drag_dur_edit.text())

    def initiate(self):
        """ Initiate or stop an etMINFLUX experiment. """
        if not self.__running:
            # get timings from ROI input
            self.getTimings()
            # read mfx sequence and lasers and laser powers from GUI
            sequenceIdx = self._widget.mfx_seq_par.currentIndex()
            self.mfx_seq = self._widget.mfx_seqs[sequenceIdx]
            laserIdx = self._widget.mfx_exc_laser_par.currentIndex()
            self.mfx_exc_laser = self._widget.mfx_exc_lasers[laserIdx]
            self.mfx_exc_pwr = float(self._widget.mfx_exc_pwr_edit.text())
            self.mfx_act_pwr = float(self._widget.mfx_act_pwr_edit.text())

            ## read trigger modality type from GUI
            #modalityIdx = self._widget.triggerModalityPar.currentIndex()
            #self.fast_modality = self._widget.triggerModalities[modalityIdx]
            # read param for using all ROIs
            self.__run_all_aoi = self._widget.triggerAllROIsCheck.isChecked()
            # read param for triggering random ROI from binary mask
            self.__random_roi_bin = self._widget.triggerRandomROICheck.isChecked()
            # read params for analysis pipeline from GUI
            self.__pipeline_param_vals = self.readPipelineParams()
            # reset general run parameters
            self.resetRunParams()
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
                self.__plotROI = self._widget.plotROICheck.isChecked()
            # check if visualization or validation mode
            if self.__runMode == RunMode.TestValidate or self.__runMode == RunMode.TestVisualize or self.__plotROI:
                self.launchHelpWidget()
            # load selected coordinate transform
            self.loadTransform()
            self.__transformCoeffs = self.__coordTransformHelper.getTransformCoeffs()
            # read param for line-wise analysis pipeline runs and decide analysis period
            self.__confocalLinesFrame = self._imspector.active_measurement().stack(0).sizes()[0]
            self.__lineWiseAnalysis = self._widget.lineWiseAnalysisCheck.isChecked()
            if self.__lineWiseAnalysis:
                self.__confocalLinesAnalysisPeriod = int(self._widget.lines_analysis_edit.text())
            else:
                self.__confocalLinesAnalysisPeriod = self.__confocalLinesFrame
            self.__confocalLineAnalysisCurr = 0
            self.__confocalLineFrameCurr = 0
            # read number of initial frames without analysis
            self.__init_frames = int(self._widget.init_frames_edit.text()) - 1
            # start confocal imaging loop (in its own thread, to wait for the .run() function that returns a status when the frame finished (or is it measurement?))
            self._imspector.connect_end(self.imspectorLineEvent, 1) # connect Imspector confocal image frame finished to running pipeline
            self.startConfocalScanning()
            self._widget.initiateButton.setText('Stop')
            self.__running = True
        elif self.__runningMFX:
            self.__run_all_aoi = self._widget.triggerAllROIsCheck.isChecked()
            if self.recTimeTimer.isActive():
                self.recTimeTimer.stop()
            self.timerThread.quit()
            self.stopMFX()
            self.scanEnded()
        else:
            # disconnect signals and reset parameters
            self._imspector.disconnect_end(self.imspectorLineEvent, 1)  # disconnect Imspector confocal image frame finished to running pipeline
            self._imspector.pause()
            self._widget.initiateButton.setText('Initiate')
            self.resetPipelineParamVals()
            self.resetRunParams()

    def finishMFXROIAuto(self):
        """ Trigger this when a preset-rec-time MFX ROI has finished. """
        self.__runningMFX = False
        self.__run_all_aoi = self._widget.triggerAllROIsCheck.isChecked()
        self.scanEnded()

    def imspectorLineEvent(self):
        self.__confocalLineAnalysisCurr += 1
        self.__confocalLineFrameCurr += 1
        if self.__confocalLineAnalysisCurr == self.__confocalLinesAnalysisPeriod:
            self.__confocalLineAnalysisCurr = 0
            self.analysisPeriodTrigger()
        if self.__confocalLineFrameCurr >= self.__confocalLinesFrame:
            self.__confocalLineFrameCurr = 0
            self.__fast_frame += 1
            self.bufferLatestImages()

    def startConfocalScanning(self):
        # ensure activate configuration is conf overview
        meas = self._imspector.active_measurement()
        cfg = meas.configuration(self.__conf_config)
        meas.activate(cfg)
        # set repeat measurement and start scan
        mouse.move(*self.__set_repeat_meas_button_pos)
        mouse.click()

    def runConfocalScan(self):
        # ensure activate configuration is conf overview
        meas = self._imspector.active_measurement()
        cfg = meas.configuration(self.__conf_config)
        meas.activate(cfg)
        # start scan
        self._imspector.start()

    def scanEnded(self):
        """ End a MINFLUX acquisition. """
        idx = self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque)
        self.setDetLogLine(f"mfx_end_{idx}", None)
        print([self.__run_all_aoi, self.__aoi_coords_deque, self.__followingROIMode])
        if (not self.__run_all_aoi or not self.__aoi_coords_deque) and (not self.__followingROIMode):
            time.sleep(self.__sleepTimeROISwitch)  # to make sure Imspector has properly turned off MINFLUX recording
            self.endExperiment()
            self.continueFastModality()
            self.__fast_frame = 0
        elif self.__run_all_aoi:
            time.sleep(self.__sleepTimeROISwitch)  # to make sure Imspector has properly turned off MINFLUX recording
            self.newROIMFX(self.__aoi_coords_deque.popleft(), roi_idx=self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque))
        elif self.__followingROIMode:
            time.sleep(self.__sleepTimeROISwitch)  # to make sure Imspector has properly turned off MINFLUX recording
            self.continueFastModality()

    def setDetLogLine(self, key, val, *args):
        if val is None:
            val = datetime.now().strftime('%Hh%Mm%Ss%fus')
        if args:
            self.__detLog[f"{key}{args[0]}"] = val
        else:
            self.__detLog[key] = val

    def endExperiment(self):
        """ Save all etMINFLUX scan metadata/log file and data. """
        self.setDetLogLine("pipeline", self.getPipelineName())
        self.logPipelineParamVals()
        # save log file with temporal info of trigger event
        filename_prefix = datetime.now().strftime('%y%m%d-%H%M%S')
        name = os.path.join(_logsDir, filename_prefix) + '_log'
        savename = getUniqueName(name)
        log = [f'{key}: {self.__detLog[key]}' for key in self.__detLog]
        with open(f'{savename}.txt', 'w') as f:
            [f.write(f'{st}\n') for st in log]
        self.resetDetLog()
        # save confocal data leading up to event image
        self.saveValidationImages(prev=True, prev_ana=False, path_prefix=filename_prefix)
        if self.__autoSaveMeas:
            # save all MINFLUX data (including multiple ROI, if so)
            self.saveMINFLUXdata(path_prefix=filename_prefix)

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
        #params_ignore = ['img','bkg','binary_mask','testmode','exinfo']
        param_names = list()
        for pipeline_param_name, _ in self.__pipeline_params.items():
            if pipeline_param_name not in self.__params_exclude:
                param_names.append(pipeline_param_name)
        for key, val in zip(param_names, self.__pipeline_param_vals):
            self.setDetLogLine(key, val)

    def continueFastModality(self):
        """ Continue the confocal imaging, after a MINFLUX event acquisition has been performed. """
        if self._widget.endlessScanCheck.isChecked() or self.__followingROIModeContinue:
            time.sleep(self.__sleepTime)
            # connect end of frame signals and start scanning
            self.__confocalLineCurr = 0
            # start confocal imaging loop
            self._imspector.connect_end(self.imspectorLineEvent, 1) # connect Imspector confocal image frame finished to running pipeline
            self.startConfocalScanning()
            self._widget.initiateButton.setText('Stop')
            self.__running = True
        else:
            # disconnect signals and reset parameters
            self._imspector.disconnect_end(self.imspectorLineEvent, 1)  # disconnect Imspector confocal image frame finished to running pipeline
            self._widget.initiateButton.setText('Initiate')
            self.resetPipelineParamVals()
            self.resetRunParams()

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

    def initiateBinaryMask(self):
        """ Initiate the process of calculating a binary mask of the region of interest. """
        self.launchHelpWidget()
        self.__binary_stack = []
        self.__binary_frame = 0
        self.__confocalLinesAnalysisPeriod = self._imspector.active_measurement().stack(0).sizes()[0]
        self.__confocalLineCurr = 0
        self._imspector.connect_end(self.imspectorLineEventBinaryMask, 1)
        self.startConfocalScanning()
        self._widget.recordBinaryMaskButton.setText('Recording...')

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
        """ Calculate the binary mask of the region of interest. """
        self._imspector.pause()
        img_mean = np.mean(self.__binary_stack, 0)
        img_bin = ndi.filters.gaussian_filter(img_mean, float(self._widget.bin_smooth_edit.text()))
        self.__binary_mask = np.array(img_bin > float(self._widget.bin_thresh_edit.text()))
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
        #infotext = f'Min: {np.min(img)}, max: {np.max(img/10000)} (rel. change)'
        infotext = f'Min: {np.min(img)}, max: {np.max(img)}'
        self._widget.analysisHelpWidget.info_label.setText(infotext)

    def setBusyFalse(self):
        """ Set busy flag to false. """
        self.__busy = False

    def readPipelineParams(self):
        """ Read user-provided analysis pipeline parameter values. """
        param_vals = list()
        for item in self._widget.param_edits:
            param_vals.append(float(item.text()))
        return param_vals

    def launchHelpWidget(self):
        """ Launch help widget that shows the preprocessed images in real-time. """
        self._widget.launchHelpWidget(self._widget.analysisHelpWidget, init=True)

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
        self.__fast_frame = 0
        self.__post_event_frames = 0
        self.__aoi_coords_deque = deque(maxlen=0)
        self.__aoi_sizes_deque = deque(maxlen=0)
        self.__followingROIModeContinue = False
        if self.recTimeTimer.isActive():
            self.recTimeTimer.stop()
        self.timerThread.quit()

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
            random_coords = np.array([])
            roi_size = None
            coords_scan = None
        return random_coords, coords_scan, roi_size

    def analysisPeriodTrigger(self):
        """ Triggers when an analysis period has finished. Depending on experiment modes set: send either to runPipeline, or act in other ways. """
        # enter if not still running last analysis period trigger
        if not self.__busy:
            # set busy true
            self.__busy = True
            # get image
            meas = self._imspector.active_measurement()
            self.img = np.squeeze(meas.stack(0).data()[0])
            # if recording mode is following ROI and we are currently following a ROI
            if self.__followingROIModeContinue:
                # if we are redetecting the ROI
                if self.__followingROIRedetect:
                    coords_detected, roi_sizes = self.runPipeline()
                    coords_scan, roi_size = self.postPipelineFollowingROI(coords_detected, roi_sizes)
                    self.acquireMINFLUXFull(coords_scan, roi_size, logging=False)
                else:
                    self.acquireMINFLUXMinimal(pos=self.__roi_center_mfx, roi_size=self.__roi_size_mfx, roi_size_um=self.__roi_size_um_mfx, pos_conf=self.__pos_conf_mfx)
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
                        # move on with triggering modality switch, if we have detected or randomized coords
                        if coords_scan is not None:
                            self.acquireMINFLUXFull(coords_scan, roi_size)
                            self.helpImageGeneration(coords_detected, roi_sizes)
            # unset busy flag
            self.setBusyFalse()

    def runPipeline(self):
        """ Run the analyis pipeline on the latest confocal analysis period frame. """
        # start of pipeline
        self.__pipeline_start_file_prefix = datetime.now().strftime('%Y%m%d-%Hh%Mm%Ss')  # use for naming files
        # run pipeline
        self.setDetLogLine("pipeline_start", None)
        if self.__runMode == RunMode.TestVisualize or self.__runMode == RunMode.TestValidate:
            # if test mode: run pipeline w/ analysis image return
            coords_detected, roi_sizes, self.__exinfo, self.__img_ana = self.pipeline(self.img, self.__prevFrames,self.__binary_mask,
                                                                                    (self.__runMode==RunMode.TestVisualize or
                                                                                    self.__runMode==RunMode.TestValidate),
                                                                                    self.__exinfo, self.__presetROISize,
                                                                                    *self.__pipeline_param_vals)
        else:
            # if experiment mode: run pipeline w/o analysis image return
            coords_detected, roi_sizes, self.__exinfo = self.pipeline(self.img, self.__prevFrames, self.__binary_mask,
                                                                    self.__runMode==RunMode.TestVisualize,
                                                                    self.__exinfo, self.__presetROISize,
                                                                    *self.__pipeline_param_vals)
        self.setDetLogLine("pipeline_end", None)
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
            # plot detected coords in help widget
            self.__analysisHelper.plotScatter(coords_detected, color='g')
            self.__analysisHelper.plotRoiRectangles(coords_detected, roi_sizes, color='g', presetROISize=self.__presetROISize)
            self._widget.coordListWidget.addCoords(coords_detected, roi_sizes)

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
            self.__analysisHelper.plotScatter(coords_detected, color='g')
            self.__analysisHelper.plotRoiRectangles(coords_detected, roi_sizes, color='g', presetROISize=self.__presetROISize)
            self._widget.coordListWidget.addCoords(coords_detected, roi_sizes)
            # take first detected coords as event
            if np.size(coords_detected) > 2:
                coords_scan = coords_detected[0,:]
            else:
                coords_scan = coords_detected[0]
            # log detected center coordinate
            self.setDetLogLine("fastscan_x_center", coords_scan[0])
            self.setDetLogLine("fastscan_y_center", coords_scan[1])
            # flag for start of validation
            self.__validating = True
            self.__post_event_frames = 0

    def postPipelineExperiment(self, coords_detected, roi_sizes):
        """ Called if in experiment mode, and returns MINFLUX scan coordinates according to GUI settings, detected or random. """
        # if to generate random roi(s) from binary (flag: __random_roi_bin)
        if self.__random_roi_bin:
            coords_detected, coords_scan, roi_size = self.generateRandomCoord(roi_sizes)
        # else handle detected coordinates, if some events were detected
        elif coords_detected.size != 0:
            # if we want to run all detected ROIs (flag: __run_all_roi)
            if self.__run_all_aoi:
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
                self._widget.initiateButton.setText('Next ROI')
            else:
                # take first detected coord (and roi_size if applicable) as event
                if np.size(coords_detected) > 2:
                    coords_scan = coords_detected[0,:]
                else:
                    coords_scan = coords_detected[0]
                if not self.__presetROISize:
                    roi_size = roi_sizes[0]
                else:
                    roi_size = None
        else:
            ### TODO: have to check if this is correct null-returning if we do not have any detected or random coordinates. Could be that I can return None, None for coords_scan and roi_size as well.
            coords_detected = np.array([])
            coords_scan = np.array([])
            roi_size = np.array([])
        return coords_detected, coords_scan, roi_size

    def postPipelineFollowingROI(self, coords_detected, roi_sizes):
        """ Called if in experiment mode and following ROI mode and following ROI redetection mode, and returns MINFLUX scan coordinates according to GUI settings. """
        if coords_detected.size != 0:
            # TODO Check if some detected coordinate is close to the previous one, by calculating distance from old to all detected ones
            # TODO If one detected is inside some reasonable range, use this as the new position and roi_size, and return these new values

            # take first detected coord (and roi_size if applicable) as event
            if np.size(coords_detected) > 2:
                coords_scan = coords_detected[0,:]
            else:
                coords_scan = coords_detected[0]
            if not self.__presetROISize:
                roi_size = roi_sizes[0]
        return coords_scan, roi_size

    def acquireMINFLUXMinimal(self, pos, roi_size, roi_size_um, pos_conf):
        """ Minimal MINFLUX acquisition function. """
        # pause fast imaging
        self.pauseFastModality()
        # initiate and run scanning
        self.initiateMFX(position=pos, ROI_size=roi_size, ROI_size_um=roi_size_um, pos_conf=pos_conf)
        self.startMFX()

    def acquireMINFLUXFull(self, coords_scan, roi_size, logging=True):
        """ Full MINFLUX acquisition function. """
        if logging:
            self.setDetLogLine("prepause", None)
        # pause fast imaging
        self.pauseFastModality()
        if logging:
            self.setDetLogLine("coord_transf_start", None)
        # transform detected coordinate between from pixels to sample position in um (for conf --> MFX)
        coords_center_scan, coords_center_um, self.__px_size_mon = self.transform(coords_scan, self.__transformCoeffs) 
        # get roi size in um from preset values
        if self.__presetROISize:
            roi_size_scan = [float(self._widget.size_x_edit.text()), float(self._widget.size_y_edit.text())]
            roi_size_um_scan = roi_size_scan
        # or calculate roi size in um from in confocal image pixels
        else:
            coord_top, coord_um_top, _ = self.transform([coords_scan[0]+roi_size[0]/2, coords_scan[1]+roi_size[1]/2],self.__transformCoeffs)
            coord_bot, coord_um_bot, _ = self.transform([coords_scan[0]-roi_size[0]/2, coords_scan[1]-roi_size[1]/2],self.__transformCoeffs)
            roi_size_scan = np.subtract(coord_top, coord_bot)
            roi_size_um_scan = np.subtract(coord_um_top, coord_um_bot)
        # get preset recording time, if set
        if self.__presetRecTime:
            self.__rec_time_scan = float(self._widget.mfx_rectime_edit.text())
        # else run indefinitely
        else:
            self.__rec_time_scan = None
        # initiate and run scanning with transformed center coordinate
        self.initiateMFX(position=coords_center_scan, ROI_size=roi_size_scan, ROI_size_um=roi_size_um_scan, pos_conf=coords_scan)
        if logging:
            # log detected and scanning center coordinate, with different keys if running all ROIs or only one
            if self.__run_all_aoi:
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, idx=self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque))
                self.setDetLogLine(f"mfx_initiate_{self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque)}", None)
            else:
                self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan)
                self.setDetLogLine("mfx_initiate", None)
        self.startMFX()

    def helpImageGeneration(self, coords_detected, roi_sizes):
        """ Called after starting MINFLUX acquisition, if some additional info regarding detected coordinates should be viewed/saved. """
        self.setEventsImage(coords_detected)
        if self.__plotROI:
            # set analysis image in help widget and plot detected coords in help widget
            self.setAnalysisHelpImg(self.__img_ana)
            self.__analysisHelper.plotScatter(coords_detected, color='g')
            self.__analysisHelper.plotRoiRectangles(coords_detected, roi_sizes, color='g', presetROISize=self.__presetROISize)
            self._widget.coordListWidget.addCoords(coords_detected, roi_sizes)

    def setEventsImage(self, coords):
        """ Create an image with pixel-marked events in the imspector measurement. """
        time.sleep(self.__sleepTime*3)
        meas = self._imspector.active_measurement()
        # activate confocal configuration
        meas.activate(meas.configuration('ov conf'))
        # create new data stack with same size as confocal
        meas.create_stack(int, np.roll(np.shape(meas.stack(0).data()),2))
        # get data stack of newly created stack
        data = meas.stack(meas.number_of_stacks()-1).data()
        # add pixels for the detected event coordinates
        for coord_pair in coords:
            data[0,0,int(coord_pair[1]),int(coord_pair[0])] = 1
        # dilate event positions for visibility
        kernel = [[False, True, False],[True, True, True],[False, True, False]]
        for _ in range(2):
            data[0,0] = ndi.binary_dilation(data[0,0], kernel)
        stack = meas.stack(meas.number_of_stacks()-1)
        stack.set_name('DetectedEvents')

    def bufferLatestImages(self):
        # buffer latest fast frame and (if applicable) validation images
        self.__prevFrames.append(np.copy(self.img))
        if self.__runMode == RunMode.TestValidate:
            # if validation mode: buffer previous preprocessed analysis frame
            self.__prevAnaFrames.append(self.__img_ana)

    def newROIMFX(self, coords_scan, roi_idx):
        # switch active Imspector window to conf overview
        meas = self._imspector.active_measurement()
        cfg = meas.configuration(self.__conf_config)
        meas.activate(cfg)
        # transform detected coordinate between from pixels to sample position in um
        coords_center_scan, coords_center_um, self.__px_size_mon = self.transform(coords_scan, self.__transformCoeffs) 
        if self.__presetROISize:
            roi_size_scan = [float(self._widget.size_x_edit.text()), float(self._widget.size_y_edit.text())]
            roi_size_um_scan = roi_size_scan
        else:
            roi_size = self.__aoi_sizes_deque.popleft()
            coord_top, coord_um_top, _ = self.transform([coords_scan[0]+roi_size[0]/2, coords_scan[1]+roi_size[1]/2],self.__transformCoeffs)
            coord_bot, coord_um_bot, _ = self.transform([coords_scan[0]-roi_size[0]/2, coords_scan[1]-roi_size[1]/2],self.__transformCoeffs)
            roi_size_scan = np.subtract(coord_top, coord_bot)
            roi_size_um_scan = np.subtract(coord_um_top, coord_um_bot)
        time.sleep(self.__sleepTime)
        # log detected and scanning center coordinate
        self.logCoordinates(coords_scan, coords_center_um, roi_size_um_scan, idx=roi_idx)
        # get preset rec time, if using that mode
        if self.__presetRecTime:
            rec_time_scan = float(self._widget.mfx_rectime_edit.text())
            self.__rec_time_scan = rec_time_scan
        else:
            self.__rec_time_scan = None
        # initiate and run scanning with transformed center coordinate
        self.initiateMFX(position=coords_center_scan, ROI_size=roi_size_scan, ROI_size_um=roi_size_um_scan, pos_conf=coords_scan)
        self.setDetLogLine(f"mfx_initiate_{roi_idx}", None)
        self.startMFX()

    def logCoordinates(self, coords_scan, coords_center_scan, roi_size, idx=None):
        """ Log detection and scan coordinates. """
        if idx:
            self.setDetLogLine(f"x_center_px_{idx}", coords_scan[0])
            self.setDetLogLine(f"y_center_px_{idx}", coords_scan[1])
            self.setDetLogLine(f"x_center_um_{idx}", coords_center_scan[0])
            self.setDetLogLine(f"y_center_um_{idx}", coords_center_scan[1])
            self.setDetLogLine(f"x_roi_size_um_{idx}", roi_size[0])
            self.setDetLogLine(f"y_roi_size_um_{idx}", roi_size[1])
        else:
            self.setDetLogLine("x_center_px", coords_scan[0])
            self.setDetLogLine("y_center_px", coords_scan[1])
            self.setDetLogLine("x_center_um", coords_center_scan[0])
            self.setDetLogLine("y_center_um", coords_center_scan[1])
            self.setDetLogLine("x_roi_size_um", roi_size[0])
            self.setDetLogLine("y_roi_size_um", roi_size[1])

    def initiateMFX(self, position=[0.0,0.0], ROI_size=[10,10], ROI_size_um=[1.0,1.0], pos_conf=[0,0]):
        """ Prepare a MINFLUX scan at the defined position. """
        self.setMFXROI(position, ROI_size)
        time.sleep(self.__sleepTime)
        self.setMFXSequence(self.mfx_seq)
        if self.__followingROIMode:
            # parameter save current parameters, to reuse for later rounds
            self.__roi_center_mfx = position
            self.__roi_size_mfx = ROI_size
            self.__roi_size_um_mfx = ROI_size_um
            self.__pos_conf_mfx = pos_conf
            self.__followingROIModeContinue = True
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
        self.setMFXDataTag(pos_conf, ROI_size_um, self.__aoi_coords_deque.maxlen-len(self.__aoi_coords_deque))
        self.setMFXLasers(self.mfx_exc_laser, self.mfx_exc_pwr, self.mfx_act_pwr)
        time.sleep(self.__sleepTime)

    def startRecTimer(self, deadtime=True, time=None):
        if deadtime:
            self.rec_time_scan = int((self.__rec_time_scan + self.__rec_time_deadtime) * 1000)  # time in ms
        else:
            self.rec_time_scan = int(time * 1000)  # time in ms
        self.recTimeTimer.setInterval(self.rec_time_scan)
        self.timerThread.start()

    def getConfocalInterval(self):
        self.__follow_roi_confocal_interval = int(self._widget.follow_roi_interval_edit.text())

    def setMFXSequence(self, mfx_seq):
        """ Sets MINFLUX sequence, according to the GUI choice of the user. """
        self._imspector.value_at('Minflux/sequence_id', specpy.ValueTree.Measurement).set(mfx_seq)

    def setMFXDataTag(self, position, roi_size, roi_idx):
        datatag = 'ROI'+str(roi_idx)+'-Pos['+str(position[0])+','+str(position[1])+']'+'-Size['+f'{roi_size[0]:.2f}'+','+f'{roi_size[1]:.2f}'+']'
        if self.__presetRecTime:
            datatag = datatag + '-RecTime['+str(self.__rec_time_scan)+'s]'
        self._imspector.value_at('Minflux/tag', specpy.ValueTree.Measurement).set(datatag)

    def setMFXRecTime(self, rec_time_int=None):
        if rec_time_int is None:
            # setting rec time
            rec_time_int = int(self.__rec_time_scan)
        self._imspector.value_at('Minflux/flow/stop_time', specpy.ValueTree.Measurement).set(rec_time_int)

    def setMFXLasers(self, exc_laser, exc_pwr, act_pwr):
        """ Sets MINFLUX lasers and laser powers, according to the GUI choice of the user. """
        laser_exc_idxs = [6,3,1]
        # set all excitation to 0 and off
        for laser_exc_idx in laser_exc_idxs:
            self._imspector.value_at('ExpControl/measurement/channels/0/lasers/'+str(laser_exc_idx)+'/active', specpy.ValueTree.Measurement).set(False)
            self._imspector.value_at('ExpControl/measurement/channels/0/lasers/'+str(laser_exc_idx)+'/power/calibrated', specpy.ValueTree.Measurement).set(0)
        # set excitation
        laser_exc_idx = str(laser_exc_idxs[self.mfxExcLaserList.index(exc_laser)])
        laser_status_exc = True
        if exc_pwr == 0:
            exc_pwr = 0.1
        self._imspector.value_at('ExpControl/measurement/channels/0/lasers/'+laser_exc_idx+'/active', specpy.ValueTree.Measurement).set(laser_status_exc)
        self._imspector.value_at('ExpControl/measurement/channels/0/lasers/'+laser_exc_idx+'/power/calibrated', specpy.ValueTree.Measurement).set(exc_pwr)
        # set activation
        laser_status_act = True if act_pwr > 0 else False
        self._imspector.value_at('ExpControl/measurement/channels/0/lasers/0/active', specpy.ValueTree.Measurement).set(laser_status_act)
        self._imspector.value_at('ExpControl/measurement/channels/0/lasers/0/power/calibrated', specpy.ValueTree.Measurement).set(act_pwr)
        
    def setMFXROI(self, position, ROI_size):
        """ Set the MINFLUX ROI by mouse control: drag ROI, and click "Set as MFX ROI"-button"""
        # get monitor px position
        shift = 1
        if self.__presetROISize:
            positions = (position[0] - ROI_size[0]/self.__px_size_mon/2 + shift,
                        position[1] - ROI_size[1]/self.__px_size_mon/2 + shift,
                        position[0] + ROI_size[0]/self.__px_size_mon/2 + shift,
                        position[1] + ROI_size[1]/self.__px_size_mon/2 + shift)
        else:
            positions = (position[0] - ROI_size[0]/2 + shift,
                        position[1] - ROI_size[1]/2 + shift,
                        position[0] + ROI_size[0]/2 + shift,
                        position[1] + ROI_size[1]/2 + shift)
        # click two positions to remove any currently drawn ROIs
        border_margin = 5
        mock_pos1 = [self.__transformCoeffs[0]+border_margin, self.__transformCoeffs[1]+border_margin]  # top left
        mock_pos2 = [mock_pos1[0]+self.__transformCoeffs[2]-border_margin*2, mock_pos1[1]+self.__transformCoeffs[2]-border_margin*2]  # bottom right
        mouse.move(*mock_pos1)
        mouse.click()
        mouse.move(*mock_pos2)
        mouse.click()
        # drag actual ROI
        mouse.drag(*positions, absolute=True, duration=self.__mouse_drag_duration)
        mouse.move(*self.__set_MFXROI_button_pos)
        mouse.click()

    def startMFX(self):
        """ Run event-triggered MINFLUX acquisition in small ROI. """
        self._imspector.start()
        self.__runningMFX = True

    def setMFXROIButtonPosButtonCall(self):
        mouse.on_click(self.setMFXROIButtonPos)
        mouse.wait(button='left')
        time.sleep(0.3)
        self._widget.setMFXROICalibrationButtonText(self.__set_MFXROI_button_pos)

    def setRepeatMeasButtonPosButtonCall(self):
        mouse.on_click(self.setRepeatMeasButtonPos)
        mouse.wait(button='left')
        time.sleep(0.3)
        self._widget.setRepeatMeasCalibrationButtonText(self.__set_repeat_meas_button_pos)

    def setMFXROIButtonPos(self):
        mouse_pos = mouse.get_position()
        self.__set_MFXROI_button_pos = np.array(mouse_pos)
        mouse.unhook_all()

    def setRepeatMeasButtonPos(self):
        mouse_pos = mouse.get_position()
        self.__set_repeat_meas_button_pos = np.array(mouse_pos)
        mouse.unhook_all()

    def stopMFX(self):
        """ Stop MINFLUX measurement. """
        if not self.__presetRecTime:
            self._imspector.pause()
        self.__runningMFX = False

    def saveValidationImages(self, prev=True, prev_ana=True, path_prefix='YMD-HMS'):
        """ Save the validation fast images of an event detection, fast images and/or preprocessed analysis images. """
        if prev:
            self._saveImage(self.__prevFrames, path_prefix, 'conf-raw')
            self.__prevFrames.clear()
        if prev_ana:
            self._saveImage(self.__prevAnaFrames, path_prefix, 'conf-ana')
            self.__prevAnaFrames.clear()

    def saveMINFLUXdata(self, path_prefix='YMD-HMS'):
        meas = self._imspector.active_measurement()
        fileName = path_prefix + '_' + 'minflux' + '.msr'
        filePath = os.path.join(_dataDir, fileName)
        meas.save_as(filePath)

    def saveMeasurement(self):
        filename_prefix = datetime.now().strftime('%y%m%d-%H%M%S')
        self.saveMINFLUXdata(filename_prefix)

    def _saveImage(self, img, path_prefix, path_suffix):
        fileName = path_prefix + '_' + path_suffix + '.tif'
        filePath = os.path.join(_dataDir, fileName)
        tiff.imwrite(filePath, img)

    def pauseFastModality(self):
        """ Pause the fast method, when an event has been detected. """
        if self.__running:
            self._imspector.disconnect_end(self.imspectorLineEvent, 1)
            self._imspector.pause()

    def togglePresetROISize(self):
        if not self._widget.presetROISizeCheck.isChecked():
            self._widget.size_x_edit.setReadOnly(True)
            self._widget.size_y_edit.setReadOnly(True)
            self._widget.size_x_edit.setStyleSheet("color: gray;")
            self._widget.size_y_edit.setStyleSheet("color: gray;")
            self.__presetROISize = False
        else:
            self._widget.size_x_edit.setReadOnly(False)
            self._widget.size_y_edit.setReadOnly(False)
            self._widget.size_x_edit.setStyleSheet("color: black;")
            self._widget.size_y_edit.setStyleSheet("color: black;")
            self.__presetROISize = True

    def togglePresetRecTime(self):
        if not self._widget.presetMfxRecTimeCheck.isChecked():
            self._widget.mfx_rectime_edit.setReadOnly(True)
            self._widget.mfx_rectime_edit.setStyleSheet("color: gray;")
            self.__presetRecTime = False
        else:
            self._widget.mfx_rectime_edit.setReadOnly(False)
            self._widget.mfx_rectime_edit.setStyleSheet("color: black;")
            self.__presetRecTime = True

    def togglePresetAutoSave(self):
        if not self._widget.autoSaveCheck.isChecked():
            self.__autoSaveMeas = False
        else:
            self.__autoSaveMeas = True

    def toggleFollowingROIMode(self):
        if not self._widget.followROIModeCheck.isChecked():
            self.__followingROIMode = False
        else:
            self.__followingROIMode = True

    def toggleRedetectROIMode(self):
        if not self._widget.followROIRedetectCheck.isChecked():
            self.__followingROIRedetect = False
        else:
            self.__followingROIRedetect = True

    def getTransformCoeffs(self):
        return self.__transformCoeffs


class AnalysisImgHelper():
    """ Analysis image widget help controller. """
    def __init__(self, etMINFLUXController, analysisWidget, *args, **kwargs):
        self.etMINFLUXController = etMINFLUXController
        self._widget = analysisWidget
        self.levels = [0,1]
        # connect signals from widget
        self._widget.setLevelsButton.clicked.connect(self.setLevels)

    def setLevels(self):
        """ Set analysis image min,max levels from ROI. """
        min_val = float(self._widget.levelMinEdit.text())
        max_val = float(self._widget.levelMaxEdit.text())
        self._widget.img.setLevels([min_val, max_val])

    def plotScatter(self, coords, color):
        self._widget.scatterPlot.setData(x=[coord[0] for coord in coords], y=[coord[1] for coord in coords], pen=pg.mkPen(None), brush=color, symbol='x', size=8)

    def plotRoiRectangles(self, coords, roi_sizes, color, presetROISize):
        # remove previously drawn ROIs
        self._widget.removeROIs()
        # create rectangle items for each ROI
        if presetROISize:
            µm_px_size = self.etMINFLUXController.getTransformCoeffs()[4]/self.etMINFLUXController.getTransformCoeffs()[3]   # µm to pixels
            roi_size_fix = [float(self.etMINFLUXController._widget.size_x_edit.text())*µm_px_size, float(self.etMINFLUXController._widget.size_y_edit.text())*µm_px_size]
            roi_sizes = [roi_size_fix for _ in coords]
        for coord, roi_size in zip(coords, roi_sizes):
            roi_temp = pg.PlotCurveItem(x=[coord[0]-roi_size[0]/2,coord[0]+roi_size[0]/2,coord[0]+int(roi_size[0]/2),coord[0]-int(roi_size[0]/2),coord[0]-int(roi_size[0]/2)], y=[coord[1]-roi_size[1]/2,coord[1]-roi_size[1]/2,coord[1]+roi_size[1]/2,coord[1]+roi_size[1]/2,coord[1]-roi_size[1]/2], pen=pg.mkPen(color, width=2))
            self._widget.rois_draw.append(roi_temp)
        # draw ROIs
        self._widget.drawROIs()


class EtCoordTransformHelper():
    """ Coordinate transform help widget controller. """
    def __init__(self, etMINFLUXController, coordTransformWidget, saveFolder, *args, **kwargs):

        self.etMINFLUXController = etMINFLUXController
        self._widget = coordTransformWidget
        self.__saveFolder = saveFolder

        # initiate coordinate transform parameters
        self.__transformCoeffs = [0,0,0,0,0]
        self.__calibNameSuffix = '_transformParams.txt'

        # connect signals from widget
        self.etMINFLUXController._widget.coordTransfCalibButton.clicked.connect(self.calibrationLaunch)
        self._widget.saveCalibButton.clicked.connect(self.calibrationFinish)
        self._widget.loadCalibButton.clicked.connect(self.calibrationLoad)
        self._widget.conf_top_left_mon_button.clicked.connect(self.getConfocalTopLeftPixel)
        self._widget.conf_bottom_right_mon_button.clicked.connect(self.getConfocalBottomRightPixel)

        self._widget.setCalibrationList(self.__saveFolder)

        self.calibrationLoad()

    def getCurrentMouseCoordsTopLeft(self):
        self._confTopLeftPosition = mouse.get_position()
        mouse.unhook_all()
        self._widget.conf_top_x_mon_edit.setText(str(self._confTopLeftPosition[0]))
        self._widget.conf_top_y_mon_edit.setText(str(self._confTopLeftPosition[1]))

    def getConfocalTopLeftPixel(self):
        mouse.on_click(self.getCurrentMouseCoordsTopLeft)

    def getCurrentMouseCoordsBottomRight(self):
        confBottomRightPosition = mouse.get_position()
        conf_size_px_mon = int(np.mean([np.abs(confBottomRightPosition[0]-self._confTopLeftPosition[0]), np.abs(confBottomRightPosition[1]-self._confTopLeftPosition[1])]))
        self._widget.conf_size_px_mon_edit.setText(str(conf_size_px_mon))
        mouse.unhook_all()

    def getConfocalBottomRightPixel(self):
        mouse.on_click(self.getCurrentMouseCoordsBottomRight)

    def getTransformCoeffs(self):
        """ Get transformation coefficients. """
        return self.__transformCoeffs

    def calibrationLaunch(self):
        """ Launch calibration. """
        self.etMINFLUXController._widget.launchHelpWidget(self.etMINFLUXController._widget.coordTransformWidget, init=True)

    def calibrationFinish(self):
        """ Finish calibration. """
        # get coordinate transform parameter values from transform widget
        self.__transformCoeffs[0] = int(self._widget.conf_top_x_mon_edit.text())
        self.__transformCoeffs[1] = int(self._widget.conf_top_y_mon_edit.text())
        self.__transformCoeffs[2] = int(self._widget.conf_size_px_mon_edit.text())
        self.__transformCoeffs[3] = float(self._widget.conf_size_um_edit.text())
        self.__transformCoeffs[4] = int(self._widget.conf_size_px_edit.text())
        # save coordinate transform parameters
        name = datetime.now().strftime('%y%m%d-%Hh%Mm')
        filename = os.path.join(self.__saveFolder, name) + self.__calibNameSuffix
        np.savetxt(fname=filename, X=self.__transformCoeffs)

    def calibrationLoad(self):
        """ Load a previously saved transformation calibration. """
        calibNameIdx = self._widget.transformCalibrationsPar.currentIndex()
        calibName = self._widget.transformCalibrations[calibNameIdx]
        filename = os.path.join(self.__saveFolder, calibName) + self.__calibNameSuffix
        params = np.loadtxt(fname=filename, dtype=float, delimiter='\t')
        self._widget.conf_top_x_mon_edit.setText(str(int(params[0])))
        self._widget.conf_top_y_mon_edit.setText(str(int(params[1])))
        self._widget.conf_size_px_mon_edit.setText(str(int(params[2])))
        self._widget.conf_size_um_edit.setText(str(params[3]))
        self._widget.conf_size_px_edit.setText(str(int(params[4])))
        self.__transformCoeffs[0] = int(params[0])
        self.__transformCoeffs[1] = int(params[1])
        self.__transformCoeffs[2] = int(params[2])
        self.__transformCoeffs[3] = float(params[3])
        self.__transformCoeffs[4] = int(params[4])


class RunMode(enum.Enum):
    Experiment = 1
    TestVisualize = 2
    TestValidate = 3


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
        self._data = np.zeros((3,100,100))

    def data(self):
        return self.data


class ValueAt():
    def __init__(self, *args, **kwargs):
        self._val = False

    def get(self):
        self._val

    def set(self, val):
        self._val = val
    
    def trigger(self):
        pass


# Copyright (C) 2023-2023 Jonatan Alvelid
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
