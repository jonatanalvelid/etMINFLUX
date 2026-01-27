import os

import pyqtgraph as pg
import numpy as np
import matplotlib.figure as mplfig

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import QPoint
from tkinter import TkVersion
from tkinter.filedialog import askdirectory
from guielements import *


class EtMINFLUXWidget(QtWidgets.QWidget):
    """ View part of the View-Controller-based etMINFLUX smart microscopy software. A GUI for interacting with the Controller, 
    with parameter fields to tweak acquisition settings, confocal and MINFLUX, analysis pipeline settings, and various calibration options. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        print('Initializing etMINFLUX widget')

        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')
        self.setWindowTitle('etMINFLUX')

        # generate dropdown list for analysis pipelines
        self.analysisPipelines = list()
        self.analysisPipelinePar = ComboBox()
        # add all experiment modes in a dropdown list
        self.experimentModes = ['Experiment','TestVisualize','TestValidate']
        self.experimentModesPar = ComboBox()
        self.experimentModesPar_label = FieldLabel('Running mode')
        self.experimentModesPar.addItems(self.experimentModes)
        self.experimentModesPar.setCurrentIndex(0)
        # generate dropdown list for ROI following modes
        self.roiFollowingModes = ['Single','SingleRedetect','Multiple']
        self.roiFollowingModesPar = ComboBox()
        self.roiFollowingModesPar_label = FieldLabel('ROI following mode')
        self.roiFollowingModesPar.addItems(self.roiFollowingModes)
        self.roiFollowingModesPar.setCurrentIndex(0)
        # create lists for current pipeline parameters: labels and editable text fields
        self.param_names = list()
        self.param_edits = list()
        # create buttons for initiating the event-triggered imaging and loading the pipeline
        self.initiateButton = PushButton('Initiate etMINFLUX')
        self.initiateButton.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        self.loadPipelineButton = PushButton('Load pipeline')
        # create buttons for calibrating coordinate transform, recording binary mask, save current measurement
        self.coordTransfCalibButton = PushButton('Calibration settings...')
        self.saveCurrentMeasButton = PushButton('Save curr. meas.')
        self.openGuideButton = PushButton('Open guide...')
        self.binaryMaskButton = PushButton('Binary mask settings...')
        # creat button for unlocking any softlock happening
        self.softResetButton = PushButton('Soft reset')
        # create check box for endless running mode
        self.endlessScanCheck = CheckBox('Endless') 
        # create check box for scan all detected ROIs
        self.triggerAllROIsCheck = CheckBox('Trigger all ROIs')
        # create check box for pre-setting ROI size
        self.presetROISizeCheck = CheckBox('Pre-set ROI size')
        self.presetROISizeCheck.setChecked(True)
        # create check box for turning on ROI following mode (interleaved conf./MFX in one ROI)
        self.followROIModeCheck = CheckBox('ROI following (intermittent confocal)')
        # create check box for pre-setting MFX acquisition recording time
        self.presetMfxRecTimeCheck = CheckBox('Pre-set ROI rec time')
        self.presetMfxRecTimeCheck.setChecked(True)
        # create check box for linewise analysis pipeline runs
        self.lineWiseAnalysisCheck = CheckBox('Run analysis pipeline linewise')
        # create check box for randomizing ROIs from binary mask
        self.triggerRandomROICheck = CheckBox('Random ROI (bin)')
        # create check box for automatic saving
        self.autoSaveCheck = CheckBox('Auto-save .msr after event')
        self.autoSaveCheck.setChecked(False)
        self.autoDeleteMFXDatasetCheck = CheckBox('Auto-delete mfx data after save')
        self.autoDeleteMFXDatasetCheck.setChecked(False)
        # create check box for plotting ROI even in experiment mode
        self.plotROICheck = CheckBox('Plot ROI (experiment mode)')
        # create check box for confocal monitoring pausing between frames
        self.confocalFramePauseCheck = CheckBox('Confocal frame pause (s)')
        # create check box for using two-threaded MINFLUX recording
        self.twoThreadsMFXCheck = CheckBox('Two-threaded MINFLUX')
        # create editable field for number of initial frames without analysis
        self.init_frames_label = FieldLabel('Initial frames')
        self.init_frames_edit = LineEdit(str(0))
        # create editable fields for MFX acquisition parameters
        self.size_x_label = FieldLabel('MFX ROI size X (µm)')
        self.size_x_edit = LineEdit(str(0))
        self.size_y_label = FieldLabel('MFX ROI size Y (µm)')
        self.size_y_edit = LineEdit(str(0))
        self.mfx_rectime_label = FieldLabel('MFX ROI rec time (s)')
        self.mfx_rectime_edit = LineEdit(str(0))
        self.mfxth0_exc_laser_label = FieldLabel('MFX exc laser')
        self.mfxth0_exc_laser = list()
        self.mfxth0_exc_laser_par = ComboBox()
        self.mfxth0_exc_pwr_label = FieldLabel('MFX exc power (%)')
        self.mfxth0_exc_pwr_edit = LineEdit(str(0))
        self.mfxth1_exc_laser_label = FieldLabel('MFX (th1) exc laser')
        self.mfxth1_exc_laser = list()
        self.mfxth1_exc_laser_par = ComboBox()
        self.mfxth1_exc_laser_par.setEnabled(False)
        self.mfxth1_exc_pwr_label = FieldLabel('MFX (th1) exc power (%)')
        self.mfxth1_exc_pwr_edit = LineEdit(str(0))
        self.mfxth1_exc_pwr_edit.setEditable(False)
        self.mfxth0_detector_label = FieldLabel('MFX detector(s)')
        self.mfxth0_detector = list()
        self.mfxth0_detector_par = ComboBox()
        self.mfxth1_detector_label = FieldLabel('MFX (th1) detector(s)')
        self.mfxth1_detector = list()
        self.mfxth1_detector_par = ComboBox()
        self.mfxth0_seq_label = FieldLabel('MFX sequence')
        self.mfxth0_seq = list()
        self.mfxth0_seq_par = ComboBox()
        self.mfxth1_seq_label = FieldLabel('MFX (th1) sequence')
        self.mfxth1_seq = list()
        self.mfxth1_seq_par = ComboBox()
        self.mfx_act_laser_label = FieldLabel('MFX act laser')
        self.mfx_act_laser = list()
        self.mfx_act_laser_par = ComboBox()
        self.mfx_act_pwr_label = FieldLabel('MFX act power (%)')
        self.mfx_act_pwr_edit = LineEdit(str(0))
        # disable th1 fields by default
        self.mfxth1_seq_par.setEnabled(False)
        self.mfxth1_seq_label.setStyleSheet('color: rgb(83,83,83);')
        self.mfxth1_exc_laser_label.setStyleSheet('color: rgb(83,83,83);')
        self.mfxth1_exc_pwr_label.setStyleSheet('color: rgb(83,83,83);')
        self.mfxth1_detector_label.setStyleSheet('color: rgb(83,83,83);')
        self.mfxth1_seq_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(83,83,83); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
        self.mfxth1_exc_laser_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(83,83,83); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
        self.mfxth1_detector_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(83,83,83); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
        # create editable fields for analysis control parameters
        self.lines_analysis_label = FieldLabel('Analysis period (lines)')
        self.lines_analysis_edit = LineEdit(str(0))
        self.conf_frame_pause_edit = LineEdit(str(0))
        self.conf_frame_pause_edit.setEditable(False)
        # create editable fields for ROI following mode
        self.follow_roi_interval_label = FieldLabel('Confocal interval (s)')
        self.follow_roi_interval_edit = LineEdit(str(0))
        self.follow_roi_redetectthresh_label = FieldLabel('Redetect dist threshold (px)')
        self.follow_roi_redetectthresh_edit = LineEdit(str(0))
        # create GUI group titles
        self.misc_title = TitleLabel('Transform, calibration, saving, guide, binary mask')
        self.minflux_title = TitleLabel('MINFLUX imaging parameters')
        self.pipeline_title = TitleLabel('Analysis pipeline')
        self.follow_ROI_mode_title = TitleLabel('ROI following mode')
        self.analysis_control_title = TitleLabel('Analysis control')
        # create updating info boxes
        self.conf_guipausetimer_edit_nullmessage = 'Time until next confocal: not running'
        self.conf_guipausetimer_edit = LineEdit(self.conf_guipausetimer_edit_nullmessage)
        self.conf_guipausetimer_edit.setEditable(False)
        self.conf_frame_edit_nullmessage = 'Confocal frames acquired: not running'
        self.conf_frame_edit = LineEdit(self.conf_frame_edit_nullmessage)
        self.conf_frame_edit.setEditable(False)
        # initiate pipeline parameter field for number of confocal channels
        self.channels_name = FieldLabel('{}'.format('Req. confocal channels'))
        self.channels_edit = LineEdit('-')
        self.channels_edit.setEditable(False)

        # help widget for coordinate transform
        self.coordTransformWidget = CoordTransformWidget(*args, **kwargs)
        # help widget for binary mask
        self.binaryMaskWidget = BinaryMaskWidget(*args, **kwargs)
        # help widget for showing images from the analysis pipelines, i.e. binary masks or analysed images in live
        self.analysisHelpWidget = AnalysisWidget()
        # help widget for showing coords list
        self.coordListWidget = CoordListWidget(*args, **kwargs)
        self.analysisHelpWidget.addWidgetRight(self.coordListWidget)
        # help widget for event viewing
        self.eventViewWidget = EventWidget()
        # help widget for displaying guide
        self.guideWidget = GuideWidget()

        # create grid and grid layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)

        # initialize widget controls
        currentRow = 0
        self.grid.addWidget(self.initiateButton, currentRow, 0)
        self.grid.addWidget(self.endlessScanCheck, currentRow, 1)
        self.grid.addWidget(self.experimentModesPar_label, currentRow, 2)
        self.grid.addWidget(self.experimentModesPar, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.pipeline_title, currentRow, 0, 1, 2)
        self.grid.addWidget(self.minflux_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.loadPipelineButton, currentRow, 0)
        self.grid.addWidget(self.analysisPipelinePar, currentRow, 1)
        self.grid.addWidget(self.mfxth0_seq_label, currentRow, 2)
        self.grid.addWidget(self.mfxth0_seq_par, currentRow, 3)
        self.grid.addWidget(self.twoThreadsMFXCheck, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.channels_name, currentRow, 0)
        self.grid.addWidget(self.channels_edit, currentRow, 1)
        self.grid.addWidget(self.mfxth0_exc_laser_label, currentRow, 2)
        self.grid.addWidget(self.mfxth0_exc_laser_par, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfxth0_exc_pwr_label, currentRow, 2)
        self.grid.addWidget(self.mfxth0_exc_pwr_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfxth0_detector_label, currentRow, 2)
        self.grid.addWidget(self.mfxth0_detector_par, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfxth1_seq_label, currentRow, 2)
        self.grid.addWidget(self.mfxth1_seq_par, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfxth1_exc_laser_label, currentRow, 2)
        self.grid.addWidget(self.mfxth1_exc_laser_par, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfxth1_exc_pwr_label, currentRow, 2)
        self.grid.addWidget(self.mfxth1_exc_pwr_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfxth1_detector_label, currentRow, 2)
        self.grid.addWidget(self.mfxth1_detector_par, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfx_act_laser_label, currentRow, 2)
        self.grid.addWidget(self.mfx_act_laser_par, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfx_act_pwr_label, currentRow, 2)
        self.grid.addWidget(self.mfx_act_pwr_edit, currentRow, 3)
        self.grid.addWidget(self.triggerRandomROICheck, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.size_x_label, currentRow, 2)
        self.grid.addWidget(self.size_x_edit, currentRow, 3)
        self.grid.addWidget(self.triggerAllROIsCheck, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.size_y_label, currentRow, 2)
        self.grid.addWidget(self.size_y_edit, currentRow, 3)
        self.grid.addWidget(self.presetROISizeCheck, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.mfx_rectime_label, currentRow, 2)
        self.grid.addWidget(self.mfx_rectime_edit, currentRow, 3)
        self.grid.addWidget(self.presetMfxRecTimeCheck, currentRow, 4)
        currentRow += 20  # arbitrary large enough number, to separate sections when loading pipelines with many parameters
        self.grid.addWidget(self.follow_ROI_mode_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.roiFollowingModesPar_label, currentRow, 2)
        self.grid.addWidget(self.roiFollowingModesPar, currentRow, 3)
        self.grid.addWidget(self.followROIModeCheck, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.follow_roi_interval_label, currentRow, 2)
        self.grid.addWidget(self.follow_roi_interval_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.analysis_control_title, currentRow, 0, 1, 2)
        self.grid.addWidget(self.follow_roi_redetectthresh_label, currentRow, 2)
        self.grid.addWidget(self.follow_roi_redetectthresh_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.lines_analysis_label, currentRow, 0)
        self.grid.addWidget(self.lines_analysis_edit, currentRow, 1)
        self.grid.addWidget(self.misc_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.plotROICheck, currentRow, 0)
        self.grid.addWidget(self.lineWiseAnalysisCheck, currentRow, 1)
        self.grid.addWidget(self.saveCurrentMeasButton, currentRow, 2)
        self.grid.addWidget(self.autoSaveCheck, currentRow, 3)
        self.grid.addWidget(self.autoDeleteMFXDatasetCheck, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.init_frames_label, currentRow, 0)
        self.grid.addWidget(self.init_frames_edit, currentRow, 1)
        self.grid.addWidget(self.openGuideButton, currentRow, 2)
        self.grid.addWidget(self.coordTransfCalibButton, currentRow, 3)
        self.grid.addWidget(self.binaryMaskButton, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.confocalFramePauseCheck, currentRow, 0)
        self.grid.addWidget(self.conf_frame_pause_edit, currentRow, 1)
        self.grid.addWidget(self.softResetButton, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.conf_guipausetimer_edit, currentRow, 0, 1, 2)
        currentRow += 1
        self.grid.addWidget(self.conf_frame_edit, currentRow, 0, 1, 2)

        frame_gm = self.frameGeometry()
        topLeftPoint = QtWidgets.QApplication.desktop().availableGeometry().topLeft()
        frame_gm.moveTopLeft(topLeftPoint)
        self.move(frame_gm.topLeft())

    def resetHelpWidget(self):
        self.analysisHelpWidget = AnalysisWidget()
        self.coordListWidget = CoordListWidget()
        self.analysisHelpWidget.addWidgetRight(self.coordListWidget)

    def resetEventViewWidget(self):
        self.eventViewWidget = EventWidget()
        
    def setDefaultValues(self, setupInfo: dict, activation_lasers_present: bool):
        """ Set default values from the setup info dictionary loaded from json file. """
        # set confocal template settings default values
        self.coordTransformWidget.template_conf_stackidx_edit.setText(str(setupInfo.get('acquisition_settings').get('confocal_stack_idx')))
        # set binary mask calculation default values
        self.binaryMaskWidget.bin_thresh_edit.setText(str(setupInfo.get('analysis_settings').get('binary_pos_thresh_def')))
        self.binaryMaskWidget.bin_smooth_edit.setText(str(setupInfo.get('analysis_settings').get('binary_pos_smooth_def')))
        self.binaryMaskWidget.bin_neg_thresh_edit.setText(str(setupInfo.get('analysis_settings').get('binary_neg_thresh_def')))
        self.binaryMaskWidget.bin_neg_smooth_edit.setText(str(setupInfo.get('analysis_settings').get('binary_neg_smooth_def')))
        self.binaryMaskWidget.bin_border_size_edit.setText(str(setupInfo.get('analysis_settings').get('binary_border_size_def')))
        # set MINFLUX laser power default values
        self.mfxth0_exc_pwr_edit.setText(str(setupInfo.get('acquisition_settings').get('minflux_exc_power_def')))
        self.mfxth1_exc_pwr_edit.setText(str(setupInfo.get('acquisition_settings').get('minflux_exc_power_th1_def')))
        # set MINFLUX ROI size and recording time default values
        self.size_x_edit.setText(str(setupInfo.get('acquisition_settings').get('minflux_roisize_x_def')))
        self.size_y_edit.setText(str(setupInfo.get('acquisition_settings').get('minflux_roisize_y_def')))
        self.mfx_rectime_edit.setText(str(setupInfo.get('acquisition_settings').get('minflux_roitime_def')))
        # set etMINFLUX analysis control default values
        self.lines_analysis_edit.setText(str(setupInfo.get('analysis_settings').get('analysis_period_lines_def')))
        self.conf_frame_pause_edit.setText(str(setupInfo.get('analysis_settings').get('analysis_confocal_frame_pause_def')))
        # set ROI following mode default values
        self.follow_roi_interval_edit.setText(str(setupInfo.get('roi_following_settings').get('confocal_interval_def')))
        self.follow_roi_redetectthresh_edit.setText(str(setupInfo.get('roi_following_settings').get('redetect_dist_thresh_def')))
        # activation laser powers default values
        if activation_lasers_present:
            self.mfx_act_pwr_edit.setText(str(setupInfo.get('acquisition_settings').get('minflux_act_power_def')))
            self.mfx_act_pwr_edit.setEditable(True)
        else:
            self.mfx_act_pwr_edit.setEditable(False)
            self.mfx_act_pwr_label.setStyleSheet('color: rgb(83,83,83);')
            self.mfx_act_laser_label.setStyleSheet('color: rgb(83,83,83);')
            self.mfx_act_laser_par.setStyleSheet('background-color: rgb(50,50,50); color: rgb(83,83,83); border: 1px solid rgb(100,100,100); selection-color: rgb(217,83,0); selection-background-color: rgb(30,30,30);')
            self.mfx_act_laser_par.setEnabled(False)

    def initParamFields(self, parameters: dict, params_exclude: list, conf_channels: int):
        """ Initialized event-triggered analysis pipeline parameter fields. """
        # remove previous parameter fields for the previously loaded pipeline
        for param in self.param_names:
            self.grid.removeWidget(param)
        for param in self.param_edits:
            self.grid.removeWidget(param)
        
        currentRow = 4
        self.channels_edit.setText(str(conf_channels))
        # initiate parameter fields for all the parameters in the pipeline chosen
        self.param_names = list()
        self.param_edits = list()
        for pipeline_param_name, pipeline_param_val in parameters.items():
            if pipeline_param_name not in params_exclude:
                # create param for input
                param_name = FieldLabel('{}'.format(pipeline_param_name))
                param_value = pipeline_param_val.default if pipeline_param_val.default is not pipeline_param_val.empty else 0
                param_edit = LineEdit(str(param_value))
                # add param name and param to grid
                self.grid.addWidget(param_name, currentRow, 0)
                self.grid.addWidget(param_edit, currentRow, 1)
                # add param name and param to object list of temp widgets
                self.param_names.append(param_name)
                self.param_edits.append(param_edit)

                currentRow += 1

    def setAnalysisPipelines(self, analysisDir):
        """ Set combobox with available analysis pipelines to use. """
        for pipeline in os.listdir(analysisDir):
            if os.path.isfile(os.path.join(analysisDir, pipeline)):
                pipeline = pipeline.split('.')[0]
                self.analysisPipelines.append(pipeline)
        self.analysisPipelinePar.addItems(self.analysisPipelines)
        self.analysisPipelinePar.setCurrentIndex(0)

    def setMfxSequenceList(self, mfxSeqs, thread=0):
        """ Set combobox with available minflux sequences to use. """
        self.mfx_seqs = mfxSeqs
        if thread == 0:
            self.mfxth0_seq_par.addItems(self.mfx_seqs)
            self.mfxth0_seq_par.setCurrentIndex(0)
        elif thread == 1:
            self.mfxth1_seq_par.addItems(self.mfx_seqs)
            self.mfxth1_seq_par.setCurrentIndex(0)

    def setMinfluxExcLaserList(self, excLasers, thread=0):
        """ Set combobox with available excitation lasers to use. """
        self.mfx_exc_lasers = excLasers
        if thread == 0:
            self.mfxth0_exc_laser_par.addItems(self.mfx_exc_lasers)
            self.mfxth0_exc_laser_par.setCurrentIndex(0)
        elif thread == 1:
            self.mfxth1_exc_laser_par.addItems(self.mfx_exc_lasers)
            self.mfxth1_exc_laser_par.setCurrentIndex(1)

    def setMinfluxActLaserList(self, actLasers):
        """ Set combobox with available excitation lasers to use. """
        self.mfx_act_lasers = ['None']+actLasers
        self.mfx_act_laser_par.addItems(self.mfx_act_lasers)
        self.mfx_act_laser_par.setCurrentIndex(0)

    def setMinfluxDetectorList(self, detectors, thread=0):
        """ Set combobox with available excitation lasers to use. """
        self.mfx_detectors = detectors
        if thread == 0:
            self.mfxth0_detector_par.addItems(self.mfx_detectors)
            self.mfxth0_detector_par.setCurrentIndex(0)
        elif thread == 1:
            self.mfxth1_detector_par.addItems(self.mfx_detectors)
            self.mfxth1_detector_par.setCurrentIndex(1)
    
    def setConfGUINullMessages(self):
        self.conf_guipausetimer_edit.setText(self.conf_guipausetimer_edit_nullmessage)
        self.conf_frame_edit.setText(self.conf_frame_edit_nullmessage)

    def launchHelpWidget(self, widget, init=True):
        """ Launch the help widget. """
        if widget == self.eventViewWidget:
            self.eventViewWidget.showwidget(init)
        elif init:
            widget.show()
        else:
            widget.hide()


class AnalysisWidget(QtWidgets.QWidget):
    """ Pop-up widget for the live analysis images or binary masks. """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')

        self.imgVbWidget = pg.GraphicsLayoutWidget()
        self.imgVb = pg.PlotItem(axisItems={'left': pg.AxisItem('left'), 'bottom': pg.AxisItem('bottom')})
        self.imgVbWidget.addItem(self.imgVb, row=1, col=1)

        self.img = pg.ImageItem(axisOrder = 'row-major')

        self.scatterPlot = pg.ScatterPlotItem()
        self.rois_draw = []

        self.imgVb.addItem(self.img)
        self.imgVb.setAspectLocked(True)
        self.imgVb.invertY(True)
        self.imgVb.addItem(self.scatterPlot)

        self.info_label = QtWidgets.QLabel('<image info>')

        # image min,max levels related widgets
        self.setLevelsButton = PushButton('Set levels')
        self.levelMinLabel = FieldLabel('Level min')
        self.levelMinEdit = LineEdit(str(0))
        self.levelMinEdit.setMaximumWidth(50)
        self.levelMaxLabel = FieldLabel('Level max')
        self.levelMaxEdit = LineEdit(str(1))
        self.levelMaxEdit.setMaximumWidth(50)
        
        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)

        self.grid.addWidget(self.info_label, 0, 0)
        self.grid.addWidget(self.levelMinLabel, 0, 1)
        self.grid.addWidget(self.levelMinEdit, 0, 2)
        self.grid.addWidget(self.levelMaxLabel, 0, 3)
        self.grid.addWidget(self.levelMaxEdit, 0, 4)
        self.grid.addWidget(self.setLevelsButton, 0, 6)
        self.grid.addWidget(self.imgVbWidget, 1, 0, 1, 7)

        frame_gm = self.frameGeometry()
        topLeftPoint = QPoint(0,670)
        frame_gm.moveTopLeft(topLeftPoint)
        self.move(frame_gm.topLeft())

    def removeROIs(self):
        for roi in self.rois_draw:
            self.imgVb.removeItem(roi)
        self.rois_draw = []

    def removeROI(self, idx):
        try:
            self.imgVb.removeItem(self.rois_draw[idx])
            del self.rois_draw[idx]
        except:
            pass

    def drawROIs(self):
        for roi in self.rois_draw:
            self.imgVb.addItem(roi)

    def addWidgetRight(self, widget):
        self.grid.addWidget(widget, 1, 7, 1, 1)


class EventWidget(QtWidgets.QWidget):
    """ Pop-up widget for the live view on a detected event, with confocal image stack viewing as well as a plot over intensity. """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')
        self.setWindowTitle("Event view")
        
        # Create image viewer and intensity graph widgets
        self.image_viewer = StackViewerWidget(self)
        self.intensity_graph = IntensityGraphWidget(self)
        
        # Create a grid layout and add widgets
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(self.image_viewer, 0, 0)
        self.grid.addWidget(self.intensity_graph, 0, 1)

        frame_gm = self.frameGeometry()
        topLeftPoint = QPoint(0,700)
        frame_gm.moveTopLeft(topLeftPoint)
        self.move(frame_gm.topLeft())
        
    def add_frame_and_intensity(self, image, intensity):
        # Add frame to image viewer and update intensity trace
        self.image_viewer.add_frame(image)
        # Update intensity graph with the full trace
        self.intensity_graph.update_plot_postevent(intensity)

    def add_event(self, image_stack, intensity_trace):
        self.image_viewer.add_frames_event(image_stack)
        self.intensity_graph.update_plot_event(intensity_trace)

    def reset(self):
        self.image_viewer.reset()
        self.intensity_graph.reset()

    def showwidget(self, init):
        if init:
            self.show()
        else:
            self.hide()


class IntensityGraphWidget(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')

        # Setup matplotlib figure
        self.figure = mplfig.Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ax = self.figure.add_subplot(111)

        # set layout and add figure
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(self.canvas)

        self.reset()
        
    def update_plot(self):
        self.ax.clear()
        # plot vertical lines for event
        self.ax.axvline(x=self.event_frame, linewidth=4, color='g', alpha=0.6)
        self.ax.plot(self.intensity_trace, 'r-')
        self.ax.set_title("Event area intensity")
        self.ax.set_xlabel("Frame")
        self.ax.set_ylabel("Intensity")
        self.canvas.draw()

    def update_plot_event(self, intensity_trace):
        if len(intensity_trace) == 1:
            self.event_frame = 0
        else:
            self.event_frame = len(intensity_trace)
        self.intensity_trace = intensity_trace
        self.update_plot()

    def update_plot_postevent(self, intensity_trace):
        self.intensity_trace = intensity_trace
        self.update_plot()

    def reset(self):
        self.event_frame = None
        self.intensity_trace = []


class StackViewerWidget(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')

        self.stack = []
        self.current_frame = 0
        self.auto_playing = False
        
        # Slider
        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.valueChanged.connect(self.change_frame)
        
        # Auto play button
        self.auto_play_button = PushButton("Auto play")
        self.auto_play_button.setCheckable(True)
        self.auto_play_button.clicked.connect(self.toggle_auto_play)

        # Timer for auto-play
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.next_frame)

        # Image viewbox and image
        self.imgVbWidget = pg.GraphicsLayoutWidget()
        self.imgVb = pg.PlotItem(axisItems={'left': pg.AxisItem('left'), 'bottom': pg.AxisItem('bottom')})
        self.imgVbWidget.addItem(self.imgVb, row=1, col=1)
        self.img = pg.ImageItem(axisOrder = 'row-major')
        self.scatterPlot = pg.ScatterPlotItem()
        self.rois_draw = []
        self.imgVb.addItem(self.img)
        self.imgVb.setAspectLocked(True)
        self.imgVb.invertY(True)
        self.imgVb.addItem(self.scatterPlot)
        self.img.setLevels([0, 100])

        # image min,max levels related widgets
        self.setLevelsButton = PushButton('Set levels')
        self.levelMinLabel = FieldLabel('Level min')
        self.levelMinEdit = LineEdit(str(0))
        self.levelMinEdit.setMaximumWidth(50)
        self.levelMaxLabel = FieldLabel('Level max')
        self.levelMaxEdit = LineEdit(str(100))
        self.levelMaxEdit.setMaximumWidth(50)

        # Create layour and add widgets
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
        self.layout().addWidget(self.imgVbWidget, 0, 0, 1, 5)
        self.layout().addWidget(self.slider, 1, 0, 1, 4)
        self.layout().addWidget(self.auto_play_button, 1, 4)
        self.grid.addWidget(self.levelMinLabel, 2, 0)
        self.grid.addWidget(self.levelMinEdit, 2, 1)
        self.grid.addWidget(self.levelMaxLabel, 2, 2)
        self.grid.addWidget(self.levelMaxEdit, 2, 3)
        self.grid.addWidget(self.setLevelsButton, 2, 4)

        self.reset()

    def add_frame(self, frame):
        self.stack.append(frame)
        self.slider.setMaximum(len(self.stack) - 1)
        self.update_display()

    def add_frames_event(self, frames):
        self.reset()
        for frame in frames:
            self.stack.append(frame)
        self.slider.setMaximum(len(self.stack) - 1)
        self.update_display()

    def change_frame(self):
        self.current_frame = self.slider.value()
        self.update_display()

    def toggle_auto_play(self):
        if self.auto_play_button.isChecked():
            self.timer.start(250)  # set frame interval to 250 ms
        else:
            self.timer.stop()

    def next_frame(self):
        self.current_frame = (self.current_frame + 1) % len(self.stack)
        self.slider.setValue(self.current_frame)
        self.update_display()

    def update_display(self):
        if self.stack:
            frame = self.stack[self.current_frame]
            self.img.setImage(frame, autoLevels=False)
        
    def reset(self):
        self.stack = []
        self.slider.setMaximum(0)


class CoordListWidget(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')

        self.list = QtWidgets.QListWidget()
        self.list.setVerticalScrollMode(QtWidgets.QListWidget.ScrollPerPixel)

        # create widget fields for the ROI list
        self.delROIButton = PushButton('Delete ROI')
        self.numevents_edit_nullmessage = 'Number of detected events: not running'
        self.numevents_edit = LineEdit(self.numevents_edit_nullmessage)
        self.numevents_edit.setEditable(False)
        
        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(self.delROIButton, 0, 0)
        self.grid.addWidget(self.numevents_edit, 1, 0)
        self.grid.addWidget(self.list, 2, 0)

    def addCoords(self, coord_list, roi_sizes, colors):
        self.clearList()
        for coord, roi_size, color in zip(coord_list, roi_sizes, colors):
            self.addCoord(coord, roi_size, color)

    def addCoord(self, coord, roi_size, color):
        listitem = QtWidgets.QListWidgetItem(f'Pos (px): [{coord[0]},{coord[1]}], ROI size (µm): [{roi_size[0]},{roi_size[1]}]')
        self.list.addItem(listitem)

    def deleteCoord(self, idx):
        self.list.takeItem(0)

    def clearList(self):
        last_item = True
        while last_item is not None:
            last_item = self.list.takeItem(0)


class CoordTransformWidget(QtWidgets.QWidget):
    """ Pop-up widget for the coordinate transform between the two etMINFLUX modalities. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')

        # create buttons
        self.setDeleteMFXDatasetButton = PushButton('Top MFX dataset')
        self.setSaveDirButton = PushButton('Choose save directory')
        self.save_dir_edit = LineEdit('Default data folder')
        self.save_dir_edit.setEditable(False)

        # Create titles
        self.gui_calibration_title = TitleLabel('GUI calibration')
        self.template_title = TitleLabel('Template calibration')
        self.saving_title = TitleLabel('Saving')
        self.transform_title = TitleLabel('Transform')

        # generate dropdown list for coordinate transformations
        self.transformPipelines = list()
        self.transformPipelinePar = ComboBox()
        self.transformPipelinePar.setEnabled(False)
        self.transformPipelineLabel = FieldLabel('Transform pipeline')

        # generate fields for the template settings
        self.template_conf_stackidx_label = FieldLabel('Confocal stack index')
        self.template_conf_stackidx_edit = LineEdit(str(0))

        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
    
        currentRow = 0
        self.grid.addWidget(self.gui_calibration_title, currentRow, 0, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.setDeleteMFXDatasetButton, currentRow, 0, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.template_title, currentRow, 0, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.template_conf_stackidx_label, currentRow, 0)
        self.grid.addWidget(self.template_conf_stackidx_edit, currentRow, 1, 1, 2)
        currentRow += 1
        self.grid.addWidget(self.saving_title, currentRow, 0, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.setSaveDirButton, currentRow, 0)
        self.grid.addWidget(self.save_dir_edit, currentRow, 1, 1, 2)
        currentRow += 1
        self.grid.addWidget(self.transform_title, currentRow, 0, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.transformPipelineLabel, currentRow, 0)
        self.grid.addWidget(self.transformPipelinePar, currentRow, 1, 1, 2)
        
        frame_gm = self.frameGeometry()
        topLeftPoint = QPoint(0,670)
        frame_gm.moveTopLeft(topLeftPoint)
        self.move(frame_gm.topLeft())

    def setDeleteMFXDatasetButtonText(self, coords):
        self.setDeleteMFXDatasetButton.setText(f'Top MFX dataset: [{coords[0]},{coords[1]}]')

    def getSaveFolder(self):
        return askdirectory(title='Select folder...')
    
    def setSaveFolderField(self, folder):
        self.save_dir_edit.setText(folder)

    def setTransformations(self, transformDir):
        """ Set combobox with available coordinate transformations to use. """
        for transform in os.listdir(transformDir):
            if os.path.isfile(os.path.join(transformDir, transform)):
                transform = transform.split('.')[0]
                self.transformPipelines.append(transform)
        self.transformPipelinePar.addItems(self.transformPipelines)
        self.transformPipelinePar.setCurrentIndex(0)


class BinaryMaskWidget(QtWidgets.QWidget):
    """ Pop-up widget for controlling the optional binary mask for analysis restriction. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')

        # create title
        self.binary_title = TitleLabel('Binary mask')
        # create buttons
        self.recordBinaryMaskButton = PushButton('Record binary mask')
        self.resetBinaryMaskButton = PushButton('Reset binary mask')
        # create editable fields for binary mask calculation threshold and smoothing
        self.bin_thresh_label = FieldLabel('Bin. pos. threshold (cnts)')
        self.bin_thresh_edit = LineEdit(str(0))
        self.bin_smooth_label = FieldLabel('Bin. pos. smooth (px)')
        self.bin_smooth_edit = LineEdit(str(0))
        self.bin_neg_thresh_label = FieldLabel('Bin. neg. threshold (cnts)')
        self.bin_neg_thresh_edit = LineEdit(str(0))
        self.bin_neg_smooth_label = FieldLabel('Bin. neg. smooth (px)')
        self.bin_neg_smooth_edit = LineEdit(str(0))
        self.bin_border_size_label = FieldLabel('Bin. border size (px)')
        self.bin_border_size_edit = LineEdit(str(0))

        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
    
        currentRow = 0
        self.grid.addWidget(self.binary_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.bin_smooth_label, currentRow, 2)
        self.grid.addWidget(self.bin_smooth_edit, currentRow, 3)
        self.grid.addWidget(self.recordBinaryMaskButton, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.bin_thresh_label, currentRow, 2)
        self.grid.addWidget(self.bin_thresh_edit, currentRow, 3)
        self.grid.addWidget(self.resetBinaryMaskButton, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.bin_neg_smooth_label, currentRow, 2)
        self.grid.addWidget(self.bin_neg_smooth_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.bin_neg_thresh_label, currentRow, 2)
        self.grid.addWidget(self.bin_neg_thresh_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.bin_border_size_label, currentRow, 2)
        self.grid.addWidget(self.bin_border_size_edit, currentRow, 3)
        
        frame_gm = self.frameGeometry()
        topLeftPoint = QPoint(0,670)
        frame_gm.moveTopLeft(topLeftPoint)
        self.move(frame_gm.topLeft())


class GuideWidget(QtWidgets.QWidget):
    """ Pop-up widget for displaying the etMINFLUX experiment guide. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')

        # create text browser
        self.textBrowser = QtWidgets.QTextBrowser()
        self.textBrowser.setStyleSheet('background-color: rgb(200,200,200)')

        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
        
        self.setWindowTitle('etMINFLUX guide')
    
        currentRow = 0
        self.grid.addWidget(self.textBrowser, currentRow, 0, 1, 3)
        
        frame_gm = self.frameGeometry()
        topLeftPoint = QPoint(0,670)
        frame_gm.moveTopLeft(topLeftPoint)
        self.move(frame_gm.topLeft())

    def setText(self, text=None, source=None):
        if text != None:
            self.textBrowser.setText(text)
        elif source != None:
            self.textBrowser.setSource(source)


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