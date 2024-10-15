import os

import pyqtgraph as pg
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtCore import QPoint

from tkinter import TkVersion
from tkinter.filedialog import askdirectory

from guielements import *


class EtMINFLUXWidget(QtWidgets.QWidget):
    """ Widget for controlling the etMINFLUX implementation. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        print('Initializing etMINFLUX widget')
        
        # set graphic style of widget
        self.setStyleSheet('background-color: rgb(70,70,70);')

        # generate dropdown list for analysis pipelines
        self.analysisPipelines = list()
        self.analysisPipelinePar = ComboBox()
        # generate dropdown list for coordinate transformations
        self.transformPipelines = list()
        self.transformPipelinePar = ComboBox()
        self.transformPipelineLabel = FieldLabel('Transform pipeline')
        # add all experiment modes in a dropdown list
        self.experimentModes = ['Experiment','TestVisualize','TestValidate']
        self.experimentModesPar = ComboBox()
        self.experimentModesPar_label = FieldLabel('Experiment mode')
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
        self.coordTransfCalibButton = PushButton('Transform calibration')
        self.recordBinaryMaskButton = PushButton('Record binary mask')
        self.resetBinaryMaskButton = PushButton('Reset binary mask')
        self.saveCurrentMeasButton = PushButton('Save curr. meas.')
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
        self.autoSaveCheck.setChecked(True)
        self.autoDeleteMFXDatasetCheck = CheckBox('Auto-delete mfx data after save')
        self.autoDeleteMFXDatasetCheck.setChecked(True)
        # create check box for plotting ROI even in experiment mode
        self.plotROICheck = CheckBox('Plot ROI (experiment mode)')
        # create check box for confocal monitoring pausing between frames
        self.confocalFramePauseCheck = CheckBox('Confocal frame pause (s)')
        # create editable fields for binary mask calculation threshold and smoothing
        self.bin_thresh_label = FieldLabel('Bin. pos. threshold (cnts)')
        self.bin_thresh_edit = LineEdit(str(10))
        self.bin_smooth_label = FieldLabel('Bin. pos. smooth (px)')
        self.bin_smooth_edit = LineEdit(str(2))
        self.bin_neg_thresh_label = FieldLabel('Bin. neg. threshold (cnts)')
        self.bin_neg_thresh_edit = LineEdit(str(10))
        self.bin_neg_smooth_label = FieldLabel('Bin. neg. smooth (px)')
        self.bin_neg_smooth_edit = LineEdit(str(2))
        self.bin_border_size_label = FieldLabel('Bin. border size (px)')
        self.bin_border_size_edit = LineEdit(str(15))
        # create editable field for number of initial frames without analysis
        self.init_frames_label = FieldLabel('Initial frames')
        self.init_frames_edit = LineEdit(str(0))
        # create editable fields for MFX acquisition parameters
        self.size_x_label = FieldLabel('MFX ROI size X (µm)')
        self.size_x_edit = LineEdit(str(2))
        self.size_y_label = FieldLabel('MFX ROI size Y (µm)')
        self.size_y_edit = LineEdit(str(2))
        self.mfx_rectime_label = FieldLabel('MFX ROI rec time (s)')
        self.mfx_rectime_edit = LineEdit(str(60))
        self.mfx_exc_laser_label = FieldLabel('MFX exc laser')
        self.mfx_exc_laser = list()
        self.mfx_exc_laser_par = ComboBox()
        self.mfx_exc_pwr_label = FieldLabel('MFX exc power (%)')
        self.mfx_exc_pwr_edit = LineEdit(str(4))
        self.mfx_act_pwr_label = FieldLabel('MFX act power (%)')
        self.mfx_act_pwr_edit = LineEdit(str(0))
        self.mfx_act_pwr_edit.setEditable(False)
        self.mfx_seq_label = FieldLabel('MFX sequence')
        self.mfx_seq = list()
        self.mfx_seq_par = ComboBox()
        # create editable fields for analysis control parameters
        self.lines_analysis_label = FieldLabel('Analysis period (lines)')
        self.lines_analysis_edit = LineEdit(str(100))
        self.conf_frame_pause_edit = LineEdit(str(10))
        self.conf_frame_pause_edit.setEditable(False)
        # create editable fields for ROI following mode
        self.follow_roi_interval_label = FieldLabel('Confocal interval (s)')
        self.follow_roi_interval_edit = LineEdit(str(60))
        self.follow_roi_redetectthresh_label = FieldLabel('Redetect dist threshold (px)')
        self.follow_roi_redetectthresh_edit = LineEdit(str(20))
        # create GUI group titles
        self.saving_title = TitleLabel('Saving')
        self.binary_title = TitleLabel('Binary mask')
        self.minflux_title = TitleLabel('MINFLUX imaging parameters')
        self.pipeline_title = TitleLabel('Analysis pipeline')
        self.transform_title = TitleLabel('Coordinate transform and GUI calibration')
        self.follow_ROI_mode_title = TitleLabel('ROI following mode')
        self.analysis_control_title = TitleLabel('Analysis control')
        # create updating info boxes
        self.conf_guipausetimer_edit_nullmessage = 'Time until confocal: not running'
        self.conf_guipausetimer_edit = LineEdit(self.conf_guipausetimer_edit_nullmessage)
        self.conf_guipausetimer_edit.setEditable(False)

        # help widget for coordinate transform
        self.coordTransformWidget = CoordTransformWidget(*args, **kwargs)
        # help widget for showing images from the analysis pipelines, i.e. binary masks or analysed images in live
        self.analysisHelpWidget = AnalysisWidget()#*args, **kwargs)
        # help widget for showing coords list
        self.coordListWidget = CoordListWidget(*args, **kwargs)
        self.analysisHelpWidget.addWidgetRight(self.coordListWidget)

        # create grid and grid layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)

        # initialize widget controls
        currentRow = 0
        self.grid.addWidget(self.initiateButton, currentRow, 0)
        self.grid.addWidget(self.endlessScanCheck, currentRow, 1)
        self.grid.addWidget(self.experimentModesPar_label, currentRow, 2)
        self.grid.addWidget(self.experimentModesPar, currentRow, 3)
        self.grid.addWidget(self.softResetButton, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.pipeline_title, currentRow, 0, 1, 2)
        self.grid.addWidget(self.transform_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.loadPipelineButton, currentRow, 0)
        self.grid.addWidget(self.analysisPipelinePar, currentRow, 1)
        self.grid.addWidget(self.transformPipelineLabel, currentRow, 2)
        self.grid.addWidget(self.transformPipelinePar, currentRow, 3)
        self.grid.addWidget(self.coordTransfCalibButton, currentRow, 4)
        currentRow += 1
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
        currentRow += 1
        self.grid.addWidget(self.minflux_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.mfx_seq_label, currentRow, 2)
        self.grid.addWidget(self.mfx_seq_par, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfx_exc_laser_label, currentRow, 2)
        self.grid.addWidget(self.mfx_exc_laser_par, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.mfx_exc_pwr_label, currentRow, 2)
        self.grid.addWidget(self.mfx_exc_pwr_edit, currentRow, 3)
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
        currentRow += 1
        self.grid.addWidget(self.follow_ROI_mode_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.analysis_control_title, currentRow, 0, 1, 2)
        self.grid.addWidget(self.roiFollowingModesPar_label, currentRow, 2)
        self.grid.addWidget(self.roiFollowingModesPar, currentRow, 3)
        self.grid.addWidget(self.followROIModeCheck, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.lines_analysis_label, currentRow, 0)
        self.grid.addWidget(self.lines_analysis_edit, currentRow, 1)
        self.grid.addWidget(self.follow_roi_interval_label, currentRow, 2)
        self.grid.addWidget(self.follow_roi_interval_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.plotROICheck, currentRow, 0)
        self.grid.addWidget(self.lineWiseAnalysisCheck, currentRow, 1)
        self.grid.addWidget(self.follow_roi_redetectthresh_label, currentRow, 2)
        self.grid.addWidget(self.follow_roi_redetectthresh_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.init_frames_label, currentRow, 0)
        self.grid.addWidget(self.init_frames_edit, currentRow, 1)
        self.grid.addWidget(self.saving_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.confocalFramePauseCheck, currentRow, 0)
        self.grid.addWidget(self.conf_frame_pause_edit, currentRow, 1)
        self.grid.addWidget(self.saveCurrentMeasButton, currentRow, 3)
        self.grid.addWidget(self.autoSaveCheck, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.conf_guipausetimer_edit, currentRow, 0, 1, 2)
        self.grid.addWidget(self.autoDeleteMFXDatasetCheck, currentRow, 4)

        frame_gm = self.frameGeometry()
        topLeftPoint = QtWidgets.QApplication.desktop().availableGeometry().topLeft()
        frame_gm.moveTopLeft(topLeftPoint)
        self.move(frame_gm.topLeft())

    def resetHelpWidget(self):
        self.analysisHelpWidget = AnalysisWidget()
        self.coordListWidget = CoordListWidget()
        self.analysisHelpWidget.addWidgetRight(self.coordListWidget)

    def initParamFields(self, parameters: dict, params_exclude: list):
        """ Initialized event-triggered analysis pipeline parameter fields. """
        # remove previous parameter fields for the previously loaded pipeline
        for param in self.param_names:
            self.grid.removeWidget(param)
        for param in self.param_edits:
            self.grid.removeWidget(param)

        # initiate parameter fields for all the parameters in the pipeline chosen
        currentRow = 3
        
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

    def setTransformations(self, transformDir):
        """ Set combobox with available coordinate transformations to use. """
        for transform in os.listdir(transformDir):
            if os.path.isfile(os.path.join(transformDir, transform)):
                transform = transform.split('.')[0]
                self.transformPipelines.append(transform)
        self.transformPipelinePar.addItems(self.transformPipelines)
        self.transformPipelinePar.setCurrentIndex(0)

    def setMfxSequenceList(self, mfxSeqs):
        """ Set combobox with available minflux sequences to use. """
        self.mfx_seqs = mfxSeqs
        self.mfx_seq_par.addItems(self.mfx_seqs)
        self.mfx_seq_par.setCurrentIndex(0)

    def setMinfluxExcLaserList(self, excLasers):
        """ Set combobox with available excitation lasers to use. """
        self.mfx_exc_lasers = excLasers
        self.mfx_exc_laser_par.addItems(self.mfx_exc_lasers)
        self.mfx_exc_laser_par.setCurrentIndex(0)
    
    def setConfPauseTimerGUINullMessage(self):
        self.conf_guipausetimer_edit.setText(self.conf_guipausetimer_edit_nullmessage)

    def launchHelpWidget(self, widget, init=True):
        """ Launch the help widget. """
        if init:
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
        #self.img.translate(-0.5, -0.5)

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
        #self.mousePositionLabel = pg.LabelItem(text='Cursor: [nan, nan]', justify="right")
        
        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)

        self.grid.addWidget(self.info_label, 0, 0)
        self.grid.addWidget(self.levelMinLabel, 0, 1)
        self.grid.addWidget(self.levelMinEdit, 0, 2)
        self.grid.addWidget(self.levelMaxLabel, 0, 3)
        self.grid.addWidget(self.levelMaxEdit, 0, 4)
        #self.grid.addWidget(self.mousePositionLabel, 0, 5)
        self.grid.addWidget(self.setLevelsButton, 0, 6)
        self.grid.addWidget(self.imgVbWidget, 1, 0, 1, 7)

        frame_gm = self.frameGeometry()
        topLeftPoint = QPoint(0,670)
        frame_gm.moveTopLeft(topLeftPoint)
        self.move(frame_gm.topLeft())

    #    try:
    #        # set up mouse signal proxy to label current mouse hover position in image
    #        pg.SignalProxy(self.scatterPlot.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)
    #    except:
    #        pass

    #def mouseMoved(self, evt):
    #    try:
    #        mousePoint = self.scatterPlot.vb.mapSceneToView(evt[0])
    #        self.mousePositionLabel.setText(f'Cursor: [{mousePoint.x()}, {mousePoint.y()}]')
    #    except:
    #        pass

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
        self.saveCalibButton = PushButton('Save calibration')
        self.loadCalibButton = PushButton('Load calibration')
        #self.setMFXROICalibrationButton = QtWidgets.QPushButton('Set MFX ROI calib.')
        #self.setRepeatMeasCalibrationButton = QtWidgets.QPushButton('Rep. meas. calib.')
        self.setDeleteMFXDatasetButton = PushButton('Top MFX dataset')
        self.setSaveDirButton = PushButton('Choose save directory')
        self.save_dir_edit = LineEdit('Default data folder')
        self.save_dir_edit.setEditable(False)

        # add all previous transforms in folder to be loaded
        self.transformCalibrations = list()
        self.transformCalibrationsPar = ComboBox()
        self.transformCalibrationsPar_label = FieldLabel('Transform calibration')

        # Create editable fields for calibration parameters
        self.conf_top_x_mon_label = FieldLabel('Confocal top left pixel - x (monitor)')
        self.conf_top_x_mon_edit = LineEdit(str(0))
        self.conf_top_left_mon_button = PushButton('Click detect: confocal top left pixel')
        self.conf_top_y_mon_label = FieldLabel('Confocal top left pixel - y (monitor)')
        self.conf_top_y_mon_edit = LineEdit(str(0))
        self.conf_size_px_mon_label = FieldLabel('Confocal image size, monitor (pixels)')
        self.conf_size_px_mon_edit = LineEdit(str(0))
        self.conf_bottom_right_mon_button = PushButton('Click detect: confocal bottom right pixel')
        self.conf_size_um_label = FieldLabel('Confocal image size, scan (µm)')
        self.conf_size_um_edit = LineEdit(str(0))
        self.conf_size_px_label = FieldLabel('Confocal image size, scan (px)')
        self.conf_size_px_edit = LineEdit(str(0))
        # time sleeps and drag durations
        self.time_sleep_label = FieldLabel('Time sleeps (s)')
        self.time_sleep_edit = LineEdit(str(0.3))
        self.time_sleep_roiswitch_label = FieldLabel('Time sleeps - ROI switch (s)')
        self.time_sleep_roiswitch_edit = LineEdit(str(1))
        self.drag_dur_label = FieldLabel('Drag duration (s)')
        self.drag_dur_edit = LineEdit(str(0.15))
        self.save_time_label = FieldLabel('Save time (s)')
        self.save_time_edit = LineEdit(str(3))

        # Create titles
        self.gui_calibration_title = TitleLabel('GUI calibration')
        self.coord_transform_title = TitleLabel('Coordinate transform')
        self.timing_title = TitleLabel('Timing')
        self.saving_title = TitleLabel('Saving')

        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
    
        currentRow = 0
        self.grid.addWidget(self.coord_transform_title, currentRow, 0, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.transformCalibrationsPar_label, currentRow, 0)
        self.grid.addWidget(self.transformCalibrationsPar, currentRow, 1)
        self.grid.addWidget(self.loadCalibButton, currentRow, 2)
        currentRow += 1
        self.grid.addWidget(self.conf_top_x_mon_label, currentRow, 0)
        self.grid.addWidget(self.conf_top_x_mon_edit, currentRow, 1)
        self.grid.addWidget(self.conf_top_left_mon_button, currentRow, 2)
        currentRow += 1
        self.grid.addWidget(self.conf_top_y_mon_label, currentRow, 0)
        self.grid.addWidget(self.conf_top_y_mon_edit, currentRow, 1)
        currentRow += 1
        self.grid.addWidget(self.conf_size_px_mon_label, currentRow, 0)
        self.grid.addWidget(self.conf_size_px_mon_edit, currentRow, 1)
        self.grid.addWidget(self.conf_bottom_right_mon_button, currentRow, 2)
        currentRow += 1
        self.grid.addWidget(self.conf_size_um_label, currentRow, 0)
        self.grid.addWidget(self.conf_size_um_edit, currentRow, 1)
        currentRow += 1
        self.grid.addWidget(self.conf_size_px_label, currentRow, 0)
        self.grid.addWidget(self.conf_size_px_edit, currentRow, 1)
        self.grid.addWidget(self.saveCalibButton, currentRow, 2)
        currentRow += 1
        self.grid.addWidget(self.gui_calibration_title, currentRow, 0, 1, 3)
        currentRow += 1
        #self.grid.addWidget(self.setMFXROICalibrationButton, currentRow, 0)
        #self.grid.addWidget(self.setRepeatMeasCalibrationButton, currentRow, 2)
        self.grid.addWidget(self.setDeleteMFXDatasetButton, currentRow, 1)
        currentRow += 1
        self.grid.addWidget(self.timing_title, currentRow, 0, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.time_sleep_label, currentRow, 0)
        self.grid.addWidget(self.time_sleep_edit, currentRow, 1)
        currentRow += 1
        self.grid.addWidget(self.time_sleep_roiswitch_label, currentRow, 0)
        self.grid.addWidget(self.time_sleep_roiswitch_edit, currentRow, 1)
        currentRow += 1
        self.grid.addWidget(self.drag_dur_label, currentRow, 0)
        self.grid.addWidget(self.drag_dur_edit, currentRow, 1)
        currentRow += 1
        self.grid.addWidget(self.save_time_label, currentRow, 0)
        self.grid.addWidget(self.save_time_edit, currentRow, 1)
        currentRow += 1
        self.grid.addWidget(self.saving_title, currentRow, 0, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.setSaveDirButton, currentRow, 0)
        self.grid.addWidget(self.save_dir_edit, currentRow, 1, 1, 2)
        
        frame_gm = self.frameGeometry()
        topLeftPoint = QPoint(0,670)
        frame_gm.moveTopLeft(topLeftPoint)
        self.move(frame_gm.topLeft())

    def setCalibrationList(self, calibrationsDir):
        """ Set combobox with available coordinate transformations to use. """
        for transform in os.listdir(calibrationsDir):
            if os.path.isfile(os.path.join(calibrationsDir, transform)):
                transform = transform.split('.')[0].split('_')[0]
                self.transformCalibrations.append(transform)
        self.transformCalibrations = list(reversed(self.transformCalibrations))
        self.transformCalibrationsPar.addItems(self.transformCalibrations)
        self.transformCalibrationsPar.setCurrentIndex(0)

    #def setRepeatMeasCalibrationButtonText(self, coords):
    #    self.setRepeatMeasCalibrationButton.setText(f'Rep. meas. calib.: [{coords[0]},{coords[1]}]')

    #def setMFXROICalibrationButtonText(self, coords):
    #    self.setMFXROICalibrationButton.setText(f'Set MFX ROI calib.: [{coords[0]},{coords[1]}]')

    def setDeleteMFXDatasetButtonText(self, coords):
        self.setDeleteMFXDatasetButton.setText(f'Top MFX dataset: [{coords[0]},{coords[1]}]')

    def getSaveFolder(self):
        return askdirectory(title='Select folder...')
    
    def setSaveFolderField(self, folder):
        self.save_dir_edit.setText(folder)


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