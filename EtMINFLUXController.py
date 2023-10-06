import os
import glob
import sys
import importlib
import enum
import warnings
import time
import specpy
import mouse
from collections import deque
from datetime import datetime
from inspect import signature
#from tkinter import Tk, filedialog

import tifffile as tiff
import scipy.ndimage as ndi
import numpy as np
#from scipy.optimize import least_squares

warnings.filterwarnings("ignore")

_logsDir = os.path.join('C:\\Users\\Abberior_Admin\\Documents\\Jonatan\\etminflux-files', 'recordings', 'logs')
_dataDir = os.path.join('C:\\Users\\Abberior_Admin\\Documents\\Jonatan\\etminflux-files', 'recordings', 'data')
_transformsDir = os.path.join('C:\\Users\\Abberior_Admin\\Documents\\Jonatan\\etminflux-files', 'recordings', 'transforms')
#_logsDir = os.path.join('C:\\Users\\alvelidjonatan\\Documents\\Data\\etMINFLUX', 'recordings', 'logs')
#_transformsDir = os.path.join('C:\\Users\\alvelidjonatan\\Documents\\Data\\etMINFLUX', 'recordings', 'transforms')

class EtMINFLUXController():
    """ Linked to EtMINFLUXWidget."""

    def __init__(self, widget,  *args, **kwargs):
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

        # list of trigger modalities
        self.triggerModalityList = ['Confocal', 'Confocal fast']
        self._widget.setTriggerModalityList(self.triggerModalityList)

        # list of minflux sequences that can be triggered
        self.mfxSeqList = ['Imaging_2D', 'Imaging_3D', 'Tracking_2D', 'Tracking_2D_Fast']  # make sure that these options matches exactly those in Imspector
        self._widget.setMfxSequenceList(self.mfxSeqList)

        # list of available lasers for MFX imaging, get this list manually from Imspector control software
        self.mfxExcLaserList = ['640', '561', '488']
        self._widget.setMinfluxExcLaserList(self.mfxExcLaserList)

        # create a helper controller for the coordinate transform pop-out widget
        self.__coordTransformHelper = EtCoordTransformHelper(self, self._widget.coordTransformWidget, _transformsDir)
        self.__analysisHelper = AnalysisImgHelper(self, self._widget.analysisHelpWidget)

        # Connect EtMINFLUXWidget signals
        self._widget.initiateButton.clicked.connect(self.initiate)
        self._widget.loadPipelineButton.clicked.connect(self.loadPipeline)
        self._widget.recordBinaryMaskButton.clicked.connect(self.initiateBinaryMask)
        self._widget.setBusyFalseButton.clicked.connect(self.setBusyFalse)
        self._widget.setMFXROICalibrationButton.clicked.connect(self.setMFXROIButtonPosButtonCall)
        self._widget.setRepeatMeasCalibrationButton.clicked.connect(self.setRepeatMeasButtonPosButtonCall)
        self._widget.presetROISizeCheck.clicked.connect(self.togglePresetROISize)

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
        #self.__bkg = None  # bkg image
        self.__prevFrames = deque(maxlen=10)  # deque for previous fast frames
        self.__prevAnaFrames = deque(maxlen=10)  # deque for previous preprocessed analysis frames
        self.__binary_mask = None  # binary mask of regions of interest, used by certain pipelines, leave None to consider the whole image
        self.__binary_frames = 2  # number of frames to use for calculating binary mask 
        self.__init_frames = 0  # number of frames after initiating etMINFLUX before a trigger can occur, to allow laser power settling etc
        self.__validation_frames = 2  # number of fast frames to record after detecting an event in validation mode
        self.__params_exclude = ['img', 'prev_frames', 'binary_mask', 'exinfo', 'testmode']  # excluded pipeline parameters when loading param fields
        self.__run_all_aoi = False  # run all detected events/area flag
        self.__prev_frames_len = 0  # number of frames in a stack in the confocal measurement in imspector, i.e. frames saved leading up to event detected

        # initiate and set parameters for automatic mouse control
        self.__mouse_drag_duration = 0.04
        self.__sleepTime = 0.2
        self.__set_MFXROI_button_pos = [1529,76]
        self.__set_repeat_meas_button_pos = [1296,72]

    def getTimings(self):
        self.__sleepTime = float(self._widget.time_sleep_edit.text())
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

            # read trigger modality type from GUI
            modalityIdx = self._widget.triggerModalityPar.currentIndex()
            self.fast_modality = self._widget.triggerModalities[modalityIdx]
            # read param for using all ROIs
            self.__run_all_aoi = self._widget.triggerAllROIsCheck.isChecked()
            # Read params for analysis pipeline from GUI
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
            # check if visualization or validation mode
            if self.__runMode == RunMode.TestValidate or self.__runMode == RunMode.TestVisualize:
                self.launchHelpWidget()
            # load selected coordinate transform
            self.loadTransform()
            self.__transformCoeffs = self.__coordTransformHelper.getTransformCoeffs()
            self.__confocalLinesFrame = self._imspector.active_measurement().stack(0).sizes()[0]
            #self.__prev_frames_len = self._imspector.active_measurement().stack(0).sizes()[1]  #TODO: Check that this is the correct number
            #print(f'Prev frames len: {self.__prev_frames_len}')
            self.__confocalLineCurr = 0
            # start confocal imaging loop (in its own thread, to wait for the .run() function that returns a status when the frame finished (or is it measurement?))
            self._imspector.connect_end(self.imspectorLineEvent, 1) # connect Imspector confocal image frame finished to running pipeline
            # TODO: connect signal from end of MINFLUX measurement to scanEnded() TODO DO THIS LATER, WHEN MINFLUX IMAGING IS ACTUALLY STARTING, OR DO IT DIRECTLY PROGRAMMATICALLY AS THE IMAGING HAS TO BE STOPPED FROM HERE
            self.startConfocalScanning()
            self._widget.initiateButton.setText('Stop')
            self.__running = True
        elif self.__runningMFX:
            self.__run_all_aoi = self._widget.triggerAllROIsCheck.isChecked()
            self.stopMFX()
            self.scanEnded()
        else:
            # disconnect signals and reset parameters
            self._imspector.disconnect_end(self.imspectorLineEvent, 1)  # disconnect Imspector confocal image frame finished to running pipeline
            self._imspector.pause()
            self._widget.initiateButton.setText('Initiate')
            self.resetPipelineParamVals()
            self.resetRunParams()

    def imspectorLineEvent(self):
        self.__confocalLineCurr += 1
        if self.__confocalLineCurr == self.__confocalLinesFrame:
            self.__confocalLineCurr = 0
            self.runPipeline()

    def startConfocalScanning(self):
        # set repeat measurement and start scan
        mouse.move(*self.__set_repeat_meas_button_pos)
        mouse.click()

    def scanEnded(self):
        """ End a MINFLUX acquisition. """
        self.setDetLogLine("scan_end",datetime.now().strftime('%Ss%fus'))
        # TODO: save all acquisition data (confocal frames and MINFLUX data) - is it possible to keep multiple confocal frames always, by using a rotating set of like 5 confocal configs (all the same) in different windows? Or also a stack of 5 frames recorded (assuming we can know when one frame finished), that does just keep overwriting each other, and we could know on which frame it finished.
        if not self.__run_all_aoi or self.__aoi_deque.isempty():
            self.endExperiment()
            self.continueFastModality()
            self.__fast_frame = 0
        elif self.__run_all_aoi:
            time.sleep(3)  # to make sure Imspector is ready for setting a new MFX ROI
            # TODO: need to make sure that Imspector is with confocal image at front again, otherwise drawing of new ROI will be done in the MFX image...?
            self.newROIMFX(self.__aoi_deque.popleft(), roi_idx=self.__aoi_deque.maxlen-len(self.__aoi_deque))

    def setDetLogLine(self, key, val, *args):
        if args:
            self.__detLog[f"{key}{args[0]}"] = val
        else:
            self.__detLog[key] = val

    def endExperiment(self):
        """ Save all etMINFLUX scan metadata/log file. """
        self.setDetLogLine("pipeline", self.getPipelineName())
        self.logPipelineParamVals()
        # save log file with temporal info of trigger event
        filename = datetime.utcnow().strftime('%Hh%Mm%Ss%fus')
        name = os.path.join(_logsDir, filename) + '_log'
        savename = getUniqueName(name)
        log = [f'{key}: {self.__detLog[key]}' for key in self.__detLog]
        with open(f'{savename}.txt', 'w') as f:
            [f.write(f'{st}\n') for st in log]
        self.resetDetLog()

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
        params_ignore = ['img','bkg','binary_mask','testmode','exinfo']
        param_names = list()
        for pipeline_param_name, _ in self.__pipeline_params.items():
            if pipeline_param_name not in params_ignore:
                param_names.append(pipeline_param_name)
        for key, val in zip(param_names, self.__pipeline_param_vals):
            self.setDetLogLine(key, val)

    def continueFastModality(self):
        """ Continue the confocal imaging, after a MINFLUX event acquisition has been performed. """
        if self._widget.endlessScanCheck.isChecked():
            time.sleep(0.05)
            # connect end of frame signals and start scanning
            self.__confocalLineCurr = 0
            # start confocal imaging loop (in its own thread, to wait for the .run() function that returns a status when the frame finished (or is it measurement?))
            self._imspector.connect_end(self.imspectorLineEvent, 1) # connect Imspector confocal image frame finished to running pipeline
            # TODO: connect signal from end of MINFLUX measurement to scanEnded() TODO DO THIS LATER, WHEN MINFLUX IMAGING IS ACTUALLY STARTING, OR DO IT DIRECTLY PROGRAMMATICALLY AS THE IMAGING HAS TO BE STOPPED FROM HERE
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
        self.__confocalLinesFrame = self._imspector.active_measurement().stack(0).sizes()[0]
        self.__confocalLineCurr = 0
        self._imspector.connect_end(self.imspectorLineEventBinaryMask, 1)
        self.startConfocalScanning()
        self._widget.recordBinaryMaskButton.setText('Recording...')

    def imspectorLineEventBinaryMask(self):
        self.__confocalLineCurr += 1
        if self.__confocalLineCurr == self.__confocalLinesFrame:
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
        print(np.max(img_bin))

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

    def runPipeline(self):
        """ Run the analyis pipeline, called after every fast method frame. """
        #print("running pipeline")
        if not self.__busy:
            # if not still running pipeline on last frame
            # get image
            meas = self._imspector.active_measurement()
            #print(np.shape(meas.stack(0).data()))
            img = np.squeeze(meas.stack(0).data()[np.mod(self.__fast_frame, self.__prev_frames_len)])  # TODO NEED TO FIGURE OUT WHICH STACK TO TAKE - ALL STACKS IN A TEMPLATE (I.E. ALL CHANNELS IN ALL WINDOWS) ARE SEEMINGLY RANDOMLY ORDERED
            #print(np.shape(img))
            self.__busy = True
            # log start of pipeline
            pipeline_start_time = datetime.now().strftime('%Y%m%d-%Hh%Mm%Ss')  # use for naming files
            self.setDetLogLine("pipeline_start", datetime.now().strftime('%Ss%fus'))

            # run pipeline
            if self.__runMode == RunMode.TestVisualize or self.__runMode == RunMode.TestValidate:
                # if chosen a test mode: run pipeline with analysis image return
                coords_detected, roi_sizes, self.__exinfo, img_ana = self.pipeline(img, self.__prevFrames, self.__binary_mask,
                                                                        (self.__runMode==RunMode.TestVisualize or
                                                                        self.__runMode==RunMode.TestValidate),
                                                                        self.__exinfo, *self.__pipeline_param_vals)
            else:
                # if chosen experiment mode: run pipeline without analysis image return
                coords_detected, roi_sizes, self.__exinfo = self.pipeline(img, self.__prevFrames, self.__binary_mask,
                                                               self.__runMode==RunMode.TestVisualize,
                                                               self.__exinfo, *self.__pipeline_param_vals)
            self.setDetLogLine("pipeline_end", datetime.now().strftime('%Ss%fus'))

            #print(self.__fast_frame)
            if self.__fast_frame > self.__init_frames:
                #print("init frames passed")
                # if initial settling frames have passed
                if self.__runMode == RunMode.TestVisualize:
                    # if visualization mode: set analysis image in help widget
                    self.setAnalysisHelpImg(img_ana)
                elif self.__runMode == RunMode.TestValidate:
                    # if validation mode: set analysis image in help widget,
                    # and start to record validation frames after event
                    self.setAnalysisHelpImg(img_ana)
                    if self.__validating:
                        # if currently validating
                        if self.__post_event_frames > self.__validation_frames:
                            # if all validation frames have been recorded, pause fast imaging,
                            # end recording, and then continue fast imaging
                            self.saveValidationImages(prev=True, prev_ana=True, path_prefix=pipeline_start_time)
                            self.__fast_frame = 0
                            self.pauseFastModality()
                            self.endExperiment()
                            self.continueFastModality()
                            self.__fast_frame = 0
                            self.__validating = False
                        self.__post_event_frames += 1
                    elif coords_detected.size != 0:
                        # if some events where detected and not validating
                        # take first detected coords as event
                        if np.size(coords_detected) > 2:
                            coords_scan = coords_detected[0,:]
                        else:
                            coords_scan = coords_detected[0]
                        # log detected center coordinate
                        self.setDetLogLine("fastscan_x_center", coords_scan[0])
                        self.setDetLogLine("fastscan_y_center", coords_scan[1])
                        if not self.__presetROISize:
                            roi_size = roi_sizes[0]
                        # flag for start of validation
                        self.__validating = True
                        self.__post_event_frames = 0
                elif self.__runMode == RunMode.Experiment and coords_detected.size != 0:
                    # if experiment mode, and some events were detected
                    if self.__run_all_aoi:
                        # use all detected coords as events
                        areas_of_interest = list()
                        if np.size(coords_detected) > 2:
                            for pair in coords_detected:
                                areas_of_interest.append(pair)
                        else:
                            areas_of_interest.append(coords_detected[0])
                        self.__aoi_deque = deque(areas_of_interest, maxlen=len(areas_of_interest))
                        coords_scan = self.__aoi_deque.popleft()
                        if not self.__presetROISize:
                            self.__aoi_sizes_deque = deque(roi_sizes, maxlen=len(areas_of_interest))
                            roi_size = self.__aoi_sizes_deque.popleft()
                        self._widget.initiateButton.setText('Next ROI')
                    else:
                        # take first detected coords as event
                        if np.size(coords_detected) > 2:
                            coords_scan = coords_detected[0,:]
                        else:
                            coords_scan = coords_detected[0]
                        if not self.__presetROISize:
                            roi_size = roi_sizes[0]
                    self.setDetLogLine("prepause", datetime.now().strftime('%Ss%fus'))
                    # pause fast imaging
                    self.pauseFastModality()
                    self.setDetLogLine("coord_transf_start", datetime.now().strftime('%Ss%fus'))
                    # transform detected coordinate between from pixels to sample position in um (for conf --> MFX)
                    coords_center_scan, self.__px_size_mon = self.transform(coords_scan, self.__transformCoeffs) 
                    if self.__presetROISize:
                        roi_size_scan = [float(self._widget.size_x_edit.text()), float(self._widget.size_y_edit.text())]
                    else:
                        coord_top, _ = self.transform([coords_scan[0]+roi_size[0]/2, coords_scan[1]+roi_size[1]/2],self.__transformCoeffs)
                        coord_bot, _ = self.transform([coords_scan[0]-roi_size[0]/2, coords_scan[1]-roi_size[1]/2],self.__transformCoeffs)
                        roi_size_scan = np.subtract(coord_top, coord_bot)
                        
                    # log detected and scanning center coordinate
                    self.logCoordinates(coords_scan, coords_center_scan)
                    # initiate and run scanning with transformed center coordinate
                    self.initiateMFX(position=coords_center_scan, ROI_size=roi_size_scan)
                    self.runMFX()

                    # buffer latest fast frame and save validation images
                    self.__prevFrames.append(img)
                    self.saveValidationImages(prev=True, prev_ana=False, path_prefix=pipeline_start_time)
                    self.__busy = False
                    return
            # buffer latest fast frame and save validation images
            self.__prevFrames.append(img)
            if self.__runMode == RunMode.TestValidate:
                # if validation mode: buffer previous preprocessed analysis frame
                self.__prevAnaFrames.append(img_ana)
            self.__fast_frame += 1
            # unset busy flag
            self.setBusyFalse()
            #print('finished runPipeline')

    def newROIMFX(self, coords_scan, roi_size, roi_idx):
        # transform detected coordinate between from pixels to sample position in um
        coords_center_scan, self.__px_size_mon = self.transform(coords_scan, self.__transformCoeffs) 
        if self.__presetROISize:
            roi_size_scan = [float(self._widget.size_x_edit.text()), float(self._widget.size_y_edit.text())]
        else:
            roi_size = self.__aoi_sizes_deque.popleft()
            roi_size_scan = np.subtract(self.transform([coords_scan[0]+roi_size[0]/2, coords_scan[1]+roi_size[1]/2],self.__transformCoeffs),self.transform([coords_scan[0]-roi_size[0]/2, coords_scan[1]-roi_size[1]/2],self.__transformCoeffs))
        # log detected and scanning center coordinate
        self.logCoordinates(coords_scan, coords_center_scan, idx=roi_idx)
        # initiate and run scanning with transformed center coordinate
        self.initiateMFX(position=coords_center_scan, ROI_size=roi_size_scan)
        self.runMFX()

    def logCoordinates(self, coords_scan, coords_center_scan, idx=None):
        """ Log detection and scan coordinates. """
        if idx:
            self.setDetLogLine(f"x_center_px_{idx}", coords_scan[0])
            self.setDetLogLine(f"y_center_px_{idx}", coords_scan[1])
            self.setDetLogLine(f"x_center_um_{idx}", coords_center_scan[0])
            self.setDetLogLine(f"y_center_um_{idx}", coords_center_scan[1])
            self.setDetLogLine(f"scan_initiate_{idx}", datetime.now().strftime('%Ss%fus'))
        else:
            self.setDetLogLine("x_center_px", coords_scan[0])
            self.setDetLogLine("y_center_px", coords_scan[1])
            self.setDetLogLine("x_center_um", coords_center_scan[0])
            self.setDetLogLine("y_center_um", coords_center_scan[1])
            self.setDetLogLine("scan_initiate", datetime.now().strftime('%Ss%fus'))

    def initiateMFX(self, position=[0.0,0.0], ROI_size=[1.0,1.0]):
        """ Prepare a MINFLUX scan at the defined position. """
        self.setMFXROI(position, ROI_size)
        self.setMFXSequence(self.mfx_seq)
        self.setMFXLasers(self.mfx_exc_laser, self.mfx_exc_pwr, self.mfx_act_pwr)

    def setMFXSequence(self, mfx_seq):
        """ Sets MINFLUX sequence, according to the GUI choice of the user. """
        self._imspector.value_at('Minflux/sequence_id', specpy.ValueTree.Measurement).set(mfx_seq)

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
        if self.__presetROISize:
            positions = (position[0] - ROI_size[0]/self.__px_size_mon/2,
                        position[1] - ROI_size[1]/self.__px_size_mon/2,
                        position[0] + ROI_size[0]/self.__px_size_mon/2,
                        position[1] + ROI_size[1]/self.__px_size_mon/2)
        else:
            positions = (position[0] - ROI_size[0]/2,
                        position[1] - ROI_size[1]/2,
                        position[0] + ROI_size[0]/2,
                        position[1] + ROI_size[1]/2)
        mouse.drag(*positions, absolute=True, duration=self.__mouse_drag_duration)
        mouse.move(*self.__set_MFXROI_button_pos)
        mouse.click()

    def runMFX(self):
        """ Run event-triggered MINFLUX acquisition in small ROI. """
        self._imspector.start()
        self.__runningMFX = True

    def setMFXROIButtonPosButtonCall(self):
        mouse.on_click(self.setMFXROIButtonPos)

    def setRepeatMeasButtonPosButtonCall(self):
        mouse.on_click(self.setRepeatMeasButtonPos)

    def setMFXROIButtonPos(self):
        mouse_pos = mouse.get_position()
        self.__set_MFXROI_button_pos = np.array(mouse_pos)
        mouse.unhook_all()
        #print(self.__set_MFXROI_button_pos)

    def setRepeatMeasButtonPos(self):
        mouse_pos = mouse.get_position()
        self.__set_repeat_meas_button_pos = np.array(mouse_pos)
        mouse.unhook_all()
        #print(self.__set_repeat_meas_button_pos)

    def stopMFX(self):
        """ Stop MINFLUX measurement. """
        self._imspector.pause()
        self.__runningMFX = False

    def saveValidationImages(self, prev=True, prev_ana=True, path_prefix='YMD-HMS'):
        """ Save the validation fast images of an event detection, fast images and/or preprocessed analysis images. """
        if prev:
            #print(f'Length of prev frames: {len(self.__prevFrames)}')
            self._saveImage(self.__prevFrames, path_prefix, 'conf-raw')
            self.__prevFrames.clear()
        if prev_ana:
            #print(f'Length of prev_ana frames: {len(self.__prevAnaFrames)}')
            self._saveImage(self.__prevAnaFrames, path_prefix, 'conf-ana')
            self.__prevAnaFrames.clear()

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
        name = datetime.utcnow().strftime('%y%m%d-%Hh%Mm')
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
