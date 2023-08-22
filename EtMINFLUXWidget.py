import os

import napari
import pyqtgraph as pg
from napari.utils.translations import trans
from PyQt5 import QtWidgets, QtCore


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
        # generate dropdown list for trigger modalities
        self.triggerModality = list()
        self.triggerModalityPar = QtWidgets.QComboBox()
        self.triggerModalityPar_label = QtWidgets.QLabel('Trigger modality')
        self.triggerModalityPar_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        # add all experiment modes in a dropdown list
        self.experimentModes = ['Experiment','TestVisualize','TestValidate']
        self.experimentModesPar = QtWidgets.QComboBox()
        self.experimentModesPar_label = QtWidgets.QLabel('Experiment mode')
        self.experimentModesPar_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.experimentModesPar.addItems(self.experimentModes)
        self.experimentModesPar.setCurrentIndex(0)
        # create lists for current pipeline parameters: labels and editable text fields
        self.param_names = list()
        self.param_edits = list()
        # create buttons for initiating the event-triggered imaging and loading the pipeline
        self.initiateButton = QtWidgets.QPushButton('Initiate etMINFLUX')
        self.initiateButton.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        self.loadPipelineButton = QtWidgets.QPushButton('Load pipeline')
        # create buttons for calibrating coordinate transform, recording binary mask, loading scan params
        self.coordTransfCalibButton = QtWidgets.QPushButton('Transform calibration')
        self.recordBinaryMaskButton = QtWidgets.QPushButton('Record binary mask')
        # creat button for unlocking any softlock happening
        self.setBusyFalseButton = QtWidgets.QPushButton('Unlock softlock')
        # create check box for endless running mode
        self.endlessScanCheck = QtWidgets.QCheckBox('Endless')
        # create editable fields for binary mask calculation threshold and smoothing
        self.bin_thresh_label = QtWidgets.QLabel('Bin. threshold')
        self.bin_thresh_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.bin_thresh_edit = QtWidgets.QLineEdit(str(10))
        self.bin_smooth_label = QtWidgets.QLabel('Bin. smooth (px)')
        self.bin_smooth_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.bin_smooth_edit = QtWidgets.QLineEdit(str(2))
        # create editable fields for MFX acquisition parameters
        self.size_x_label = QtWidgets.QLabel('Size X (µm)')
        self.size_x_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.size_x_edit = QtWidgets.QLineEdit(str(2))
        self.size_y_label = QtWidgets.QLabel('Size Y (µm)')
        self.size_y_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.size_y_edit = QtWidgets.QLineEdit(str(2))
        self.mfx_exc_laser_label = QtWidgets.QLabel('MFX exc laser')
        self.mfx_exc_laser_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.mfx_exc_laser = list()
        self.mfx_exc_laser_par = QtWidgets.QComboBox()
        self.mfx_exc_pwr_label = QtWidgets.QLabel('MFX exc power (%)')
        self.mfx_exc_pwr_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.mfx_exc_pwr_edit = QtWidgets.QLineEdit(str(5))
        self.mfx_act_pwr_label = QtWidgets.QLabel('MFX act power (%)')
        self.mfx_act_pwr_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.mfx_act_pwr_edit = QtWidgets.QLineEdit(str(0))
        self.mfx_seq_label = QtWidgets.QLabel('MFX sequence')
        self.mfx_seq_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.mfx_seq = list()
        self.mfx_seq_par = QtWidgets.QComboBox()
        # help widget for coordinate transform
        self.coordTransformWidget = CoordTransformWidget(*args, **kwargs)

        # help widget for showing images from the analysis pipelines, i.e. binary masks or analysed images in live
        self.analysisHelpWidget = AnalysisWidget(*args, **kwargs)

        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)

        # initialize widget controls
        currentRow = 0

        self.grid.addWidget(self.initiateButton, currentRow, 0)
        self.grid.addWidget(self.endlessScanCheck, currentRow, 1)
        self.grid.addWidget(self.experimentModesPar_label, currentRow, 2)
        self.grid.addWidget(self.experimentModesPar, currentRow, 3)

        currentRow += 1

        self.grid.addWidget(self.loadPipelineButton, currentRow, 0)
        self.grid.addWidget(self.analysisPipelinePar, currentRow, 1)
        self.grid.addWidget(self.transformPipelinePar, currentRow, 2)
        self.grid.addWidget(self.coordTransfCalibButton, currentRow, 3)

        currentRow += 1

        self.grid.addWidget(self.bin_smooth_label, currentRow, 2)
        self.grid.addWidget(self.bin_smooth_edit, currentRow, 3)

        currentRow += 1

        self.grid.addWidget(self.bin_thresh_label, currentRow, 2)
        self.grid.addWidget(self.bin_thresh_edit, currentRow, 3)

        currentRow += 1

        self.grid.addWidget(self.recordBinaryMaskButton, currentRow, 3)

        currentRow +=1

        self.grid.addWidget(self.triggerModalityPar_label, currentRow, 2)
        self.grid.addWidget(self.triggerModalityPar, currentRow, 3)

        currentRow += 1

        self.grid.addWidget(self.mfx_seq_label, currentRow, 2)
        self.grid.addWidget(self.mfx_seq_par, currentRow, 3)

        currentRow += 1

        self.grid.addWidget(self.mfx_exc_laser_label, currentRow, 2)
        self.grid.addWidget(self.mfx_exc_laser_par, currentRow, 3)

        currentRow +=1 

        self.grid.addWidget(self.setBusyFalseButton, currentRow, 3)

    def initParamFields(self, parameters: dict, params_exclude: list):
        """ Initialized event-triggered analysis pipeline parameter fields. """
        # remove previous parameter fields for the previously loaded pipeline
        for param in self.param_names:
            self.grid.removeWidget(param)
        for param in self.param_edits:
            self.grid.removeWidget(param)

        # initiate parameter fields for all the parameters in the pipeline chosen
        currentRow = 2
        
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

    def setTriggerModalityList(self, modalityNames):
        """ Set combobox with available trigger modalities to use. """
        self.triggerModalities = modalityNames
        self.triggerModalityPar.addItems(self.triggerModalities)
        self.triggerModalityPar.setCurrentIndex(0)

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

        self.imgVb.addItem(self.img)
        self.imgVb.setAspectLocked(True)

        self.info_label = QtWidgets.QLabel('<image info>')
        
        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)

        self.grid.addWidget(self.info_label, 0, 0)
        self.grid.addWidget(self.imgVbWidget, 1, 0)


class CoordTransformWidget(QtWidgets.QWidget):
    """ Pop-up widget for the coordinate transform between the two etMINFLUX modalities. """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.loadLoResButton = QtWidgets.QPushButton('Load low-res calibration image')
        self.loadHiResButton = QtWidgets.QPushButton('Load high-res calibration image')
        self.saveCalibButton = QtWidgets.QPushButton('Save calibration')
        self.resetCoordsButton = QtWidgets.QPushButton('Reset coordinates')

        self.napariViewerLo = EmbeddedNapari()
        self.napariViewerHi = EmbeddedNapari()

        # add points layers to the viewer
        self.pointsLayerLo = self.napariViewerLo.add_points(name="lo_points", symbol='ring', size=20, face_color='green', edge_color='green')
        self.pointsLayerTransf = self.napariViewerHi.add_points(name="transf_points", symbol='cross', size=20, face_color='red', edge_color='red')
        self.pointsLayerHi = self.napariViewerHi.add_points(name="hi_points", symbol='ring', size=20, face_color='green', edge_color='green')
        

        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
    
        currentRow = 0
        self.grid.addWidget(self.loadLoResButton, currentRow, 0)
        self.grid.addWidget(self.loadHiResButton, currentRow, 1)
        
        currentRow += 1
        self.grid.addWidget(self.napariViewerLo.get_widget(), currentRow, 0)
        self.grid.addWidget(self.napariViewerHi.get_widget(), currentRow, 1)

        currentRow += 1
        self.grid.addWidget(self.saveCalibButton, currentRow, 0)
        self.grid.addWidget(self.resetCoordsButton, currentRow, 1)


class EmbeddedNapari(napari.Viewer):
    """ Napari viewer to be embedded in non-napari windows. Also includes a
    feature to protect certain layers from being removed when added using
    the add_image method. Copied from ImSwitch (github.com/kasasxav/ImSwitch).
    """

    def __init__(self, *args, show=False, **kwargs):
        super().__init__(*args, show=show, **kwargs)

        # Monkeypatch layer removal methods
        oldDelitemIndices = self.layers._delitem_indices

        def newDelitemIndices(key):
            indices = oldDelitemIndices(key)
            for index in indices[:]:
                layer = index[0][index[1]]
                if hasattr(layer, 'protected') and layer.protected:
                    indices.remove(index)
            return indices

        self.layers._delitem_indices = newDelitemIndices

        # Make menu bar not native
        self.window._qt_window.menuBar().setNativeMenuBar(False)

        # Remove unwanted menu bar items
        menuChildren = self.window._qt_window.findChildren(QtWidgets.QAction)
        for menuChild in menuChildren:
            try:
                if menuChild.text() in [trans._('Close Window'), trans._('Exit')]:
                    self.window.file_menu.removeAction(menuChild)
            except Exception:
                pass

    def add_image(self, *args, protected=False, **kwargs):
        result = super().add_image(*args, **kwargs)

        if isinstance(result, list):
            for layer in result:
                layer.protected = protected
        else:
            result.protected = protected
        return result

    def get_widget(self):
        return self.window._qt_window

    def addItem(self, item):
        item.attach(self,
                    canvas=self.window.qt_viewer.canvas,
                    view=self.window.qt_viewer.view,
                    parent=self.window.qt_viewer.view.scene,
                    order=1e6 + 8000)
