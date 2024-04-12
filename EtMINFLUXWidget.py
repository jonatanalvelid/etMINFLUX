import os

import pyqtgraph as pg
from PyQt5 import QtWidgets, QtCore, QtGui


class EtMINFLUXWidget(QtWidgets.QWidget):
    """ Widget for controlling the etMINFLUX implementation. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        print('Initializing etMINFLUX widget')
        
        # generate dropdown list for analysis pipelines
        self.analysisPipelines = list()
        self.analysisPipelinePar = QtWidgets.QComboBox()
        # generate dropdown list for coordinate transformations
        self.transformPipelines = list()
        self.transformPipelinePar = QtWidgets.QComboBox()
        self.transformPipelineLabel = QtWidgets.QLabel('Transform pipeline')
        self.transformPipelineLabel.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        # add all experiment modes in a dropdown list
        self.experimentModes = ['Experiment','TestVisualize','TestValidate']
        self.experimentModesPar = QtWidgets.QComboBox()
        self.experimentModesPar_label = QtWidgets.QLabel('Experiment mode')
        self.experimentModesPar_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.experimentModesPar.addItems(self.experimentModes)
        self.experimentModesPar.setCurrentIndex(0)
        # generate dropdown list for ROI following modes
        self.roiFollowingModes = ['Single','Multiple','SingleRedetect']
        self.roiFollowingModesPar = QtWidgets.QComboBox()
        self.roiFollowingModesPar_label = QtWidgets.QLabel('ROI following mode')
        self.roiFollowingModesPar_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.roiFollowingModesPar.addItems(self.roiFollowingModes)
        self.roiFollowingModesPar.setCurrentIndex(0)
        # create lists for current pipeline parameters: labels and editable text fields
        self.param_names = list()
        self.param_edits = list()
        # create buttons for initiating the event-triggered imaging and loading the pipeline
        self.initiateButton = QtWidgets.QPushButton('Initiate etMINFLUX')
        self.initiateButton.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        self.loadPipelineButton = QtWidgets.QPushButton('Load pipeline')
        self.setMFXROICalibrationButton = QtWidgets.QPushButton('Set MFX ROI calib.')
        self.setRepeatMeasCalibrationButton = QtWidgets.QPushButton('Rep. meas. calib.')
        # create buttons for calibrating coordinate transform, recording binary mask, save current measurement
        self.coordTransfCalibButton = QtWidgets.QPushButton('Transform calibration')
        self.recordBinaryMaskButton = QtWidgets.QPushButton('Record binary mask')
        self.saveCurrentMeasButton = QtWidgets.QPushButton('Save curr. meas.')
        # creat button for unlocking any softlock happening
        self.softResetButton = QtWidgets.QPushButton('Soft reset')
        # create check box for endless running mode
        self.endlessScanCheck = QtWidgets.QCheckBox('Endless')
        # create check box for scan all detected ROIs
        self.triggerAllROIsCheck = QtWidgets.QCheckBox('Trigger all ROIs')
        # create check box for pre-setting ROI size
        self.presetROISizeCheck = QtWidgets.QCheckBox('Pre-set ROI size')
        self.presetROISizeCheck.setChecked(True)
        # create check box for turning on ROI following mode (interleaved conf./MFX in one ROI)
        self.followROIModeCheck = QtWidgets.QCheckBox('ROI following (intermittent confocal)')
        # create check box for pre-setting MFX acquisition recording time
        self.presetMfxRecTimeCheck = QtWidgets.QCheckBox('Pre-set ROI rec time')
        self.presetMfxRecTimeCheck.setChecked(True)
        # create check box for linewise analysis pipeline runs
        self.lineWiseAnalysisCheck = QtWidgets.QCheckBox('Run analysis pipeline linewise')
        # create check box for randomizing ROIs from binary mask
        self.triggerRandomROICheck = QtWidgets.QCheckBox('Random ROI (bin)')
        # create check box for automatic saving
        self.autoSaveCheck = QtWidgets.QCheckBox('Auto-save .msr after event')
        self.autoSaveCheck.setChecked(True)
        # create check box for plotting ROI even in experiment mode
        self.plotROICheck = QtWidgets.QCheckBox('Plot ROI (experiment mode)')
        # create check box for confocal monitoring pausing between frames
        self.confocalFramePauseCheck = QtWidgets.QCheckBox('Confocal frame pause (s)')
        # create editable fields for binary mask calculation threshold and smoothing
        self.bin_thresh_label = QtWidgets.QLabel('Bin. threshold (cnts)')
        self.bin_thresh_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.bin_thresh_edit = QtWidgets.QLineEdit(str(10))
        self.bin_smooth_label = QtWidgets.QLabel('Bin. smooth (px)')
        self.bin_smooth_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.bin_smooth_edit = QtWidgets.QLineEdit(str(2))
        self.bin_neg_thresh_label = QtWidgets.QLabel('Bin. threshold (cnts)')
        self.bin_neg_thresh_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.bin_neg_thresh_edit = QtWidgets.QLineEdit(str(10))
        self.bin_neg_smooth_label = QtWidgets.QLabel('Bin. smooth (px)')
        self.bin_neg_smooth_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.bin_neg_smooth_edit = QtWidgets.QLineEdit(str(2))
        # create editable field for number of initial frames without analysis
        self.init_frames_label = QtWidgets.QLabel('Initial frames')
        self.init_frames_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.init_frames_edit = QtWidgets.QLineEdit(str(0))
        # create editable fields for MFX acquisition parameters
        self.size_x_label = QtWidgets.QLabel('MFX ROI size X (µm)')
        self.size_x_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.size_x_edit = QtWidgets.QLineEdit(str(2))
        self.size_y_label = QtWidgets.QLabel('MFX ROI size Y (µm)')
        self.size_y_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.size_y_edit = QtWidgets.QLineEdit(str(2))
        self.mfx_rectime_label = QtWidgets.QLabel('MFX ROI rec time (s)')
        self.mfx_rectime_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.mfx_rectime_edit = QtWidgets.QLineEdit(str(60))
        self.mfx_exc_laser_label = QtWidgets.QLabel('MFX exc laser')
        self.mfx_exc_laser_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.mfx_exc_laser = list()
        self.mfx_exc_laser_par = QtWidgets.QComboBox()
        self.mfx_exc_pwr_label = QtWidgets.QLabel('MFX exc power (%)')
        self.mfx_exc_pwr_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.mfx_exc_pwr_edit = QtWidgets.QLineEdit(str(4))
        self.mfx_act_pwr_label = QtWidgets.QLabel('MFX act power (%)')
        self.mfx_act_pwr_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.mfx_act_pwr_edit = QtWidgets.QLineEdit(str(0))
        self.mfx_act_pwr_edit.setReadOnly(True)  # make it non-editable currently
        self.mfx_act_pwr_edit.setStyleSheet("color: gray;")
        self.mfx_seq_label = QtWidgets.QLabel('MFX sequence')
        self.mfx_seq_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.mfx_seq = list()
        self.mfx_seq_par = QtWidgets.QComboBox()
        # create editable fields for analysis control parameters
        self.lines_analysis_label = QtWidgets.QLabel('Analysis period (lines)')
        self.lines_analysis_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.lines_analysis_edit = QtWidgets.QLineEdit(str(100))
        self.conf_frame_pause_edit = QtWidgets.QLineEdit(str(10))
        self.conf_frame_pause_edit.setReadOnly(True)  # make it non-editable currently
        self.conf_frame_pause_edit.setStyleSheet("color: gray;")
        # time sleeps and drag durations
        self.time_sleep_label = QtWidgets.QLabel('Time sleeps (s)')
        self.time_sleep_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.time_sleep_edit = QtWidgets.QLineEdit(str(0.3))
        self.time_sleep_roiswitch_label = QtWidgets.QLabel('Time sleeps - ROI switch (s)')
        self.time_sleep_roiswitch_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.time_sleep_roiswitch_edit = QtWidgets.QLineEdit(str(1))
        self.drag_dur_label = QtWidgets.QLabel('Drag duration (s)')
        self.drag_dur_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.drag_dur_edit = QtWidgets.QLineEdit(str(0.15))
        # create editable fields for ROI following mode
        self.follow_roi_interval_label = QtWidgets.QLabel('Confocal interval (s)')
        self.follow_roi_interval_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.follow_roi_interval_edit = QtWidgets.QLineEdit(str(60))
        self.follow_roi_redetectthresh_label = QtWidgets.QLabel('Redetect dist threshold (px)')
        self.follow_roi_redetectthresh_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.follow_roi_redetectthresh_edit = QtWidgets.QLineEdit(str(20))
        # create GUI group titles
        self.timing_title = QtWidgets.QLabel('Timing')
        self.timing_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.timing_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.saving_title = QtWidgets.QLabel('Saving')
        self.saving_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.saving_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.binary_title = QtWidgets.QLabel('Binary mask')
        self.binary_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.binary_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.binary_neg_title = QtWidgets.QLabel('Binary negative mask')
        self.binary_neg_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.binary_neg_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.minflux_title = QtWidgets.QLabel('MINFLUX imaging parameters')
        self.minflux_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.minflux_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.pipeline_title = QtWidgets.QLabel('Analysis pipeline')
        self.pipeline_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.pipeline_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.transform_title = QtWidgets.QLabel('Coordinate transform')
        self.transform_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.transform_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.follow_ROI_mode_title = QtWidgets.QLabel('ROI following mode')
        self.follow_ROI_mode_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.follow_ROI_mode_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.gui_calibration_title = QtWidgets.QLabel('GUI calibration')
        self.gui_calibration_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.gui_calibration_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))
        self.analysis_control_title = QtWidgets.QLabel('Analysis control')
        self.analysis_control_title.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.analysis_control_title.setFont(QtGui.QFont("Arial", weight=QtGui.QFont.Bold))

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
        self.grid.addWidget(self.analysis_control_title, currentRow, 0, 1, 2)
        self.grid.addWidget(self.timing_title, currentRow, 2, 1, 2)
        self.grid.addWidget(self.gui_calibration_title, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.lines_analysis_label, currentRow, 0)
        self.grid.addWidget(self.lines_analysis_edit, currentRow, 1)
        self.grid.addWidget(self.time_sleep_label, currentRow, 2)
        self.grid.addWidget(self.time_sleep_edit, currentRow, 3)
        self.grid.addWidget(self.setRepeatMeasCalibrationButton, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.plotROICheck, currentRow, 0)
        self.grid.addWidget(self.lineWiseAnalysisCheck, currentRow, 1)
        self.grid.addWidget(self.time_sleep_roiswitch_label, currentRow, 2)
        self.grid.addWidget(self.time_sleep_roiswitch_edit, currentRow, 3)
        self.grid.addWidget(self.setMFXROICalibrationButton, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.init_frames_label, currentRow, 0)
        self.grid.addWidget(self.init_frames_edit, currentRow, 1)
        self.grid.addWidget(self.drag_dur_label, currentRow, 2)
        self.grid.addWidget(self.drag_dur_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.confocalFramePauseCheck, currentRow, 0)
        self.grid.addWidget(self.conf_frame_pause_edit, currentRow, 1)
        self.grid.addWidget(self.follow_ROI_mode_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.saving_title, currentRow, 0, 1, 2)
        self.grid.addWidget(self.roiFollowingModesPar_label, currentRow, 2)
        self.grid.addWidget(self.roiFollowingModesPar, currentRow, 3)
        self.grid.addWidget(self.followROIModeCheck, currentRow, 4)
        currentRow += 1
        self.grid.addWidget(self.saveCurrentMeasButton, currentRow, 0)
        self.grid.addWidget(self.autoSaveCheck, currentRow, 1)
        self.grid.addWidget(self.follow_roi_interval_label, currentRow, 2)
        self.grid.addWidget(self.follow_roi_interval_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.follow_roi_redetectthresh_label, currentRow, 2)
        self.grid.addWidget(self.follow_roi_redetectthresh_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.binary_neg_title, currentRow, 2, 1, 3)
        currentRow += 1
        self.grid.addWidget(self.bin_neg_smooth_label, currentRow, 2)
        self.grid.addWidget(self.bin_neg_smooth_edit, currentRow, 3)
        currentRow += 1
        self.grid.addWidget(self.bin_neg_thresh_label, currentRow, 2)
        self.grid.addWidget(self.bin_neg_thresh_edit, currentRow, 3)

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
                param_name = QtWidgets.QLabel('{}'.format(pipeline_param_name))
                param_value = pipeline_param_val.default if pipeline_param_val.default is not pipeline_param_val.empty else 0
                param_edit = QtWidgets.QLineEdit(str(param_value))
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

    def setRepeatMeasCalibrationButtonText(self, coords):
        self.setRepeatMeasCalibrationButton.setText(f'Rep. meas. calib.: [{coords[0]},{coords[1]}]')

    def setMFXROICalibrationButtonText(self, coords):
        self.setMFXROICalibrationButton.setText(f'Set MFX ROI calib.: [{coords[0]},{coords[1]}]')

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

        self.imgVbWidget = pg.GraphicsLayoutWidget()
        self.imgVb = self.imgVbWidget.addViewBox(row=1, col=1)

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
        self.setLevelsButton = QtWidgets.QPushButton('Set levels')
        self.levelMinLabel = QtWidgets.QLabel('Level min')
        self.levelMinLabel.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.levelMinEdit = QtWidgets.QLineEdit(str(0))
        self.levelMinEdit.setMaximumWidth(50)
        self.levelMaxLabel = QtWidgets.QLabel('Level max')
        self.levelMaxLabel.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.levelMaxEdit = QtWidgets.QLineEdit(str(1))
        self.levelMaxEdit.setMaximumWidth(50)
        
        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)

        self.grid.addWidget(self.info_label, 0, 0)
        self.grid.addWidget(self.levelMinLabel, 0, 1)
        self.grid.addWidget(self.levelMinEdit, 0, 2)
        self.grid.addWidget(self.levelMaxLabel, 0, 3)
        self.grid.addWidget(self.levelMaxEdit, 0, 4)
        self.grid.addWidget(self.setLevelsButton, 0, 5)
        self.grid.addWidget(self.imgVbWidget, 1, 0, 1, 6)

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
        self.grid.addWidget(widget, 1, 6, 1, 1)


class CoordListWidget(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.list = QtWidgets.QListWidget()

        # create buttons for control of the ROI list
        self.delROIButton = QtWidgets.QPushButton('Delete ROI')

        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(self.delROIButton, 0, 0)
        self.grid.addWidget(self.list, 1, 0)

    def addCoords(self, coord_list, roi_sizes, colors):
        self.clearList()
        for coord, roi_size, color in zip(coord_list, roi_sizes, colors):
            self.addCoord(coord, roi_size, color)

    def addCoord(self, coord, roi_size, color):
        listitem = QtWidgets.QListWidgetItem(f'Pos (px): [{coord[0]},{coord[1]}], ROI size (µm): [{roi_size[0]},{roi_size[1]}]')
        #listitem.setForeground(color)  # TODO: this does not work, find out how to change text color
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
        self.saveCalibButton = QtWidgets.QPushButton('Save calibration')
        self.loadCalibButton = QtWidgets.QPushButton('Load calibration')

        # add all previous transforms in folder to be loaded
        self.transformCalibrations = list()
        self.transformCalibrationsPar = QtWidgets.QComboBox()
        self.transformCalibrationsPar_label = QtWidgets.QLabel('Transform calibration')
        self.transformCalibrationsPar_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)

        # Create editable fields for calibration parameters
        self.conf_top_x_mon_label = QtWidgets.QLabel('Confocal top left pixel - x (monitor)')
        self.conf_top_x_mon_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.conf_top_x_mon_edit = QtWidgets.QLineEdit(str(0))
        self.conf_top_left_mon_button = QtWidgets.QPushButton('Click detect: confocal top left pixel')
        self.conf_top_y_mon_label = QtWidgets.QLabel('Confocal top left pixel - y (monitor)')
        self.conf_top_y_mon_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.conf_top_y_mon_edit = QtWidgets.QLineEdit(str(0))
        self.conf_size_px_mon_label = QtWidgets.QLabel('Confocal image size, monitor (pixels)')
        self.conf_size_px_mon_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.conf_size_px_mon_edit = QtWidgets.QLineEdit(str(0))
        self.conf_bottom_right_mon_button = QtWidgets.QPushButton('Click detect: confocal bottom right pixel')
        self.conf_size_um_label = QtWidgets.QLabel('Confocal image size, scan (µm)')
        self.conf_size_um_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.conf_size_um_edit = QtWidgets.QLineEdit(str(0))
        self.conf_size_px_label = QtWidgets.QLabel('Confocal image size, scan (px)')
        self.conf_size_px_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.conf_size_px_edit = QtWidgets.QLineEdit(str(0))

        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
    
        currentRow = 0
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

        currentRow += 1
        self.grid.addWidget(self.saveCalibButton, currentRow, 0)

    def setCalibrationList(self, calibrationsDir):
        """ Set combobox with available coordinate transformations to use. """
        for transform in os.listdir(calibrationsDir):
            if os.path.isfile(os.path.join(calibrationsDir, transform)):
                transform = transform.split('.')[0].split('_')[0]
                self.transformCalibrations.append(transform)
        self.transformCalibrations = list(reversed(self.transformCalibrations))
        self.transformCalibrationsPar.addItems(self.transformCalibrations)
        self.transformCalibrationsPar.setCurrentIndex(0)        


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