import os
import glob
import sys
import importlib
import enum
import warnings
import specpy
import mouse
from collections import deque
from datetime import datetime
from inspect import signature
from tkinter import Tk, filedialog

import scipy.ndimage as ndi
import numpy as np
from scipy.optimize import least_squares

warnings.filterwarnings("ignore")

# TODO: change this to preferred path for saved logs
_logsDir = os.path.join('C:\\Users\\alvelidjonatan\\Documents\\Data\\etMINFLUX', 'recordings', 'logs_etminflux')

class EtMINFLUXController():
    """ Linked to EtMINFLUXWidget."""

    def __init__(self, widget,  *args, **kwargs):
        self._widget = widget
        
        print('Initializing etMINFLUX controller')

        # open imspector connection
        try:
            self._imspector = specpy.get_application()
        except:
            self._imspector = ImspectorMock()

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
        self.mfxSeqList = ['Imaging2D', 'Imaging3D', 'Tracking2D', 'Tracking3D']
        self._widget.setMfxSequenceList(self.mfxSeqList)

        # list of available lasers for MFX imaging, get this list manually from Imspector control software
        self.mfxExcLaserList = ['640', '561', '488']
        self._widget.setMinfluxExcLaserList(self.mfxExcLaserList)

        # create a helper controller for the coordinate transform pop-out widget
        self.__coordTransformHelper = EtCoordTransformHelper(self, self._widget.coordTransformWidget, _logsDir)

        # Connect EtMINFLUXWidget and communication channel signals
        self._widget.initiateButton.clicked.connect(self.initiate)
        self._widget.loadPipelineButton.clicked.connect(self.loadPipeline)
        self._widget.recordBinaryMaskButton.clicked.connect(self.initiateBinaryMask)
        self._widget.setBusyFalseButton.clicked.connect(self.setBusyFalse)
        self._widget.setMFXROICalibrationButton.clicked.connect(self.setMFXROIButtonPosButtonCall)

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
        self.__runMode = RunMode.Experiment  # run mode currently used
        self.__validating = False  # validation flag
        self.__busy = False  # running pipeline busy flag
        self.__bkg = None  # bkg image
        self.__prevFrames = deque(maxlen=10)  # deque for previous fast frames
        self.__prevAnaFrames = deque(maxlen=10)  # deque for previous preprocessed analysis frames
        self.__binary_mask = None  # binary mask of regions of interest, used by certain pipelines, leave None to consider the whole image
        self.__binary_frames = 10  # number of frames to use for calculating binary mask 
        self.__init_frames = 5  # number of frames after initiating etMINFLUX before a trigger can occur, to allow laser power settling etc
        self.__validation_frames = 5  # number of fast frames to record after detecting an event in validation mode
        self.__params_exclude = ['img', 'prev_frames', 'binary_mask', 'exinfo', 'testmode']  # excluded pipeline parameters when loading param fields
        self.__run_all_aoi = False  # run all detected events/area flag
        self.__areas_of_interest = list()  # list of detected events/areas center coordinates, to go through after pipeline run
        self.__aoi_idx = 0  # index of current area of interest being scanned
        self.__prev_frames_len = 10  # number of frames in a stack in the confocal measurement in imspector, i.e. frames saved leading up to event detected  TODO: set this when starting an etMINFLUX experiment, based on the size of the stack of the confocal image that can be read

        # initiate and set parameters for automatic mouse control
        self.__mouse_drag_duration = 0.025
        self.__setMFXROI_button_pos = [200,100]

    def initiate(self):
        """ Initiate or stop an etMINFLUX experiment. """
        if not self.__running:
            # read mfx sequence and laser from GUI
            sequenceIdx = self._widget.mfx_seq_par.currentIndex()
            self.mfx_seq = self._widget.mfx_seqs[sequenceIdx]
            laserIdx = self._widget.mfx_exc_laser_par.currentIndex()
            self.mfx_exc_laser = self._widget.mfx_exc_lasers[laserIdx]

            # read trigger modality type from GUI
            modalityIdx = self._widget.triggerModalityPar.currentIndex()
            self.fast_modality = self._widget.triggerModalities[modalityIdx]

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
            # start confocal imaging loop (in it's own thread, to wait for the .run() function that returns a status when the frame finished (or is it measurement?))
            self._imspector.connect_end(self.runPipeline, 1) # connect Imspector confocal image frame finished to running pipeline
            # TODO: connect signal from end of MINFLUX measurement to scanEnded() TODO DO THIS LATER, WHEN MINFLUX IMAGING IS ACTUALLY STARTING, OR DO IT DIRECTLY PROGRAMMATICALLY AS THE IMAGING HAS TO BE STOPPED FROM HERE
            self.startConfocalScanning()
            self._widget.initiateButton.setText('Stop')
            self.__running = True
        else:
            # disconnect signals and turn off fast laser
            self._imspector.disconnect_end(self.runPipeline, 1) # disconnect Imspector confocal image frame finished to running pipeline
            # TODO: disconnect signal from end of MINFLUX measurement to scanEnded()  TODO DO THIS LATER, WHEN MINFLUX IMAGING IS ACTUALLY STARTING, OR DO IT DIRECTLY PROGRAMMATICALLY AS THE IMAGING HAS TO BE STOPPED FROM HERE
            # TODO: stop Imspector confocal imaging
            self._widget.initiateButton.setText('Initiate')
            self.resetPipelineParamVals()
            self.resetRunParams()

    def startConfocalScanning(self):
        # TODO: set repeat measurement - mouse control?
        self._imspector.start()  # TODO: potentially use a.run() instead, to return data at finish, instead of connecting scan finish to runPipeline? But then it has to be in a separate thread.

    def scanEnded(self):
        """ End a MINFLUX acquisition. """
        self.setDetLogLine("scan_end",datetime.now().strftime('%Ss%fus'))
        # TODO: save all acquisition data (confocal frames and MINFLUX data) - is it possible to keep multiple confocal frames always, by using a rotating set of like 5 confocal configs (all the same) in different windows? Or also a stack of 5 frames recorded (assuming we can know when one frame finished), that does just keep overwriting each other, and we could know on which frame it finished.
        if not self.__run_all_aoi or self.__aoi_idx+1 == len(self.__areas_of_interest):
            self.endExperiment()
            self.continueFastModality()
            self.__fast_frame = 0
            self.__aoi_idx = 0
        elif self.__run_all_aoi:
            self.__aoi_idx += 1
            self.initiateSlowScan(self.__areas_of_interest[self.__aoi_idx])
            self.runSlowScan()

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
        if self._widget.endlessScanCheck.isChecked() and not self.__running:
            # connect communication channel signals
            # TODO: connect signal from update of image to running pipeline #xxx.sigUpdateImage.connect(self.runPipeline)
            self._widget.initiateButton.setText('Stop')
            self.__running = True
        elif not self._widget.endlessScanCheck.isChecked():
            # TODO: disconnect signal from end of scan to scanEnded() #xxx.sigScanEnded.disconnect(self.scanEnded)
            self._widget.initiateButton.setText('Initiate')
            self.__running = False
            self.resetPipelineParamVals()

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
        self.__binary_stack = None
        # TODO: activate configuration for stack of N confocal imgs
        self._imspector.connect_end(self.calculateBinaryMask, 0)
        self.startConfocalScanning()
        self._widget.recordBinaryMaskButton.setText('Recording...')

    def calculateBinaryMask(self):
        """ Calculate the binary mask of the region of interest. """
        # get image
        meas = self._imspector.active_measurement()
        self.__binary_stack = meas.stack(0).data()  # TODO: NEED TO FIGURE OUT WHICH STACK TO TAKE - ALL STACKS IN A TEMPLATE (I.E. ALL WINDOWS) ARE JUST IN A STACK LIKE THIS
        img_mean = np.mean(self.__binary_stack, 0)  # TODO: need to check the dims of the image stack, so that it is correct
        img_bin = ndi.filters.gaussian_filter(img_mean, np.float(self._widget.bin_smooth_edit.text()))
        self.__binary_mask = np.array(img_bin > np.float(self._widget.bin_thresh_edit.text()))
        self._widget.recordBinaryMaskButton.setText('Record binary mask')
        self.setAnalysisHelpImg(self.__binary_mask)
        self.launchHelpWidget()
        self._imspector.disconnect_end(self.calculateBinaryMask, 0)

    def setAnalysisHelpImg(self, img):
        """ Set the preprocessed image in the analysis help widget. """
        if self.__fast_frame < self.__init_frames + 3:
            autolevels = True
        else:
            autolevels = False
        self._widget.analysisHelpWidget.img.setImage(img, autoLevels=autolevels)
        infotext = f'Min: {np.min(img)}, max: {np.max(img/10000)} (rel. change)'
        self._widget.analysisHelpWidget.info_label.setText(infotext)

    def setBusyFalse(self):
        """ Set busy flag to false. """
        self.__busy = False

    def readPipelineParams(self):
        """ Read user-provided analysis pipeline parameter values. """
        param_vals = list()
        for item in self._widget.param_edits:
            param_vals.append(np.float(item.text()))
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
        if not self.__busy:
            # get image
            meas = self._imspector.active_measurement()
            img = meas.stack(0).data()[np.mod(self.__fast_frame, self.__prev_frames_len)]  # TODO NEED TO FIGURE OUT WHICH STACK TO TAKE - ALL STACKS IN A TEMPLATE (I.E. ALL WINDOWS) ARE JUST IN A STACK LIKE THIS
            # if not still running pipeline on last frame
            self.__busy = True
            # log start of pipeline
            self.setDetLogLine("pipeline_start", datetime.now().strftime('%Ss%fus'))

            # run pipeline
            if self.__runMode == RunMode.TestVisualize or self.__runMode == RunMode.TestValidate:
                # if chosen a test mode: run pipeline with analysis image return
                coords_detected, self.__exinfo, img_ana = self.pipeline(img, self.__bkg, self.__binary_mask,
                                                                        (self.__runMode==RunMode.TestVisualize or
                                                                        self.__runMode==RunMode.TestValidate),
                                                                        self.__exinfo, *self.__pipeline_param_vals)
            else:
                # if chosen experiment mode: run pipeline without analysis image return
                coords_detected, self.__exinfo = self.pipeline(img, self.__bkg, self.__binary_mask,
                                                               self.__runMode==RunMode.TestVisualize,
                                                               self.__exinfo, *self.__pipeline_param_vals)
            self.setDetLogLine("pipeline_end", datetime.now().strftime('%Ss%fus'))

            if self.__fast_frame > self.__init_frames:
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
                            self.saveValidationImages(prev=True, prev_ana=True)
                            self.pauseFastModality()
                            self.endRecording()
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
                        # flag for start of validation
                        self.__validating = True
                        self.__post_event_frames = 0
                elif self.__runMode == RunMode.Experiment and coords_detected.size != 0:
                    # if experiment mode, and some events were detected
                    # take first detected coords as event
                    if self.__run_all_aoi:
                        self.__aoi_idx = 0
                        self.__areas_of_interest = list()
                        if np.size(coords_detected) > 2:
                            for pair in coords_detected:
                                self.__areas_of_interest.append(pair)
                        else:
                            self.__areas_of_interest.append(coords_detected[0])
                        coords_scan = self.__areas_of_interest[self.__aoi_idx]
                    else:
                        if np.size(coords_detected) > 2:
                            coords_scan = coords_detected[0,:]
                        else:
                            coords_scan = coords_detected[0]
                    self.setDetLogLine("prepause", datetime.now().strftime('%Ss%fus'))
                    # pause fast imaging
                    self.pauseFastModality()
                    self.setDetLogLine("coord_transf_start", datetime.now().strftime('%Ss%fus'))
                    # transform detected coordinate between from pixels to sample position in um (for conf --> MFX)
                    coords_center_scan, self.__px_size_mon = self.transform(coords_scan, self.__transformCoeffs) 
                    # log detected and scanning center coordinate
                    self.setDetLogLine("x_center_px", coords_scan[0])
                    self.setDetLogLine("y_center_px", coords_scan[1])
                    self.setDetLogLine("x_center_um", coords_center_scan[0])
                    self.setDetLogLine("y_center_um", coords_center_scan[1])
                    self.setDetLogLine("scan_initiate", datetime.now().strftime('%Ss%fus'))
                    # initiate and run scanning with transformed center coordinate
                    self.initiateMFX(position=coords_center_scan)
                    self.runMFX()

                    # buffer latest fast frame and save validation images
                    self.__prevFrames.append(img)
                    self.saveValidationImages(prev=True, prev_ana=False)
                    self.__busy = False
                    return
            # use latest fast frame as background for next pipeline run
            self.__bkg = img
            # buffer latest fast frame and save validation images
            self.__prevFrames.append(img)
            if self.__runMode == RunMode.TestValidate:
                # if validation mode: buffer previous preprocessed analysis frame
                self.__prevAnaFrames.append(img_ana)
            self.__fast_frame += 1
            # unset busy flag
            self.setBusyFalse()

    def initiateMFX(self, position=[0.0,0.0], ROI_size=[1.0,1.0]):
        """ Initiate a MINFLUX scan, by setting a MFX ROI. """
        self.setMINFLUXROI(position, ROI_size)

    def setMINFLUXROI(self, position, ROI_size):
        """ Set the MINFLUX ROI by mouse control: drag ROI, and click "Set as MFX ROI"-button"""
        positions = (position[0] - ROI_size[0]/self.__px_size_mon/2,
                    position[1] - ROI_size[1]/self.__px_size_mon/2,
                    position[0] + ROI_size[0]/self.__px_size_mon/2,
                    position[1] + ROI_size[1]/self.__px_size_mon/2)
        mouse.drag(*positions, absolute=True, duration=self.__mouse_drag_duration)
        mouse.move(*self.__setMFXROI_button_pos)
        mouse.click()

    def runMFX(self):
        """ Run event-triggered MINFLUX acquisition in small ROI. """
        self._imspector.run()

    def setMFXROIButtonPosButtonCall(self):
        mouse.on_click(self.setMFXROIButtonPos)

    def setMFXROIButtonPos(self):
        mouse_pos = mouse.get_position()
        self.__setMFXROI_button_pos = np.array(mouse_pos)
        mouse.unhook_all()

    def saveValidationImages(self, prev=True, prev_ana=True):
        """ Save the validation fast images of an event detection, fast images and/or preprocessed analysis images. """
        if prev:
            # TODO: save detectorFast frames leading up to event
            # #xxx.sigSaveImage.emit(self.detectorFast, np.array(list(self.__prevFrames)))  # (detector, imagestack)
            self.__prevFrames.clear()
        if prev_ana:
            # TODO: save preprocessed detectorFast frames leading up to event
            # #xxx.sigSaveImage.emit(self.detectorFast, np.array(list(self.__prevAnaFrames)))  # (detector, imagestack)
            self.__prevAnaFrames.clear()

    def pauseFastModality(self):
        """ Pause the fast method, when an event has been detected. """
        if self.__running:
            # TODO: disconnect signal from update of image to running pipeline 
            #xxx.sigUpdateImage.disconnect(self.runPipeline)
            self.__running = False


class EtCoordTransformHelper():
    """ Coordinate transform help widget controller. """
    def __init__(self, etMINFLUXController, coordTransformWidget, saveFolder, *args, **kwargs):

        self.etMINFLUXController = etMINFLUXController
        self._widget = coordTransformWidget
        self.__saveFolder = saveFolder

        # initiate coordinate transform parameters
        self.__transformCoeffs = [591,149,788,80,800]

        # connect signals from widget
        etMINFLUXController._widget.coordTransfCalibButton.clicked.connect(self.calibrationLaunch)
        self._widget.saveCalibButton.clicked.connect(self.calibrationFinish)
        self._widget.conf_top_left_mon_button.clicked.connect(self.getConfocalTopLeftPixel)
        self._widget.conf_bottom_right_mon_button.clicked.connect(self.getConfocalBottomRightPixel)

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
        self.__transformCoeffs[3] = int(self._widget.conf_size_um_edit.text())
        self.__transformCoeffs[4] = int(self._widget.conf_size_px_edit.text())
        # save coordinate transform parameters
        name = datetime.utcnow().strftime('%Hh%Mm')
        filename = os.path.join(self.__saveFolder, name) + '_transformParams.txt'
        np.savetxt(fname=filename, X=self.__transformCoeffs)


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
