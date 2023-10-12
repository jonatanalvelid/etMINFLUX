import pyqtgraph as pg
import numpy as np
import sys
from PyQt5 import QtWidgets, QtGui, QtCore


class AnalysisWidget(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.imgVbWidget = pg.GraphicsLayoutWidget()
        self.imgVb = self.imgVbWidget.addViewBox(row=1, col=1)

        self.img = pg.ImageItem(axisOrder = 'row-major')

        self.scatterPlot = pg.ScatterPlotItem()

        self.imgVb.addItem(self.img)
        self.imgVb.setAspectLocked(True)
        self.imgVb.invertY(True)
        self.imgVb.addItem(self.scatterPlot)

        # generate GUI layout
        self.grid = QtWidgets.QGridLayout()
        self.setLayout(self.grid)
        self.grid.addWidget(self.imgVbWidget)


app = QtWidgets.QApplication(sys.argv)
widget = AnalysisWidget()
widget.show()

# set image
widget.img.setImage(np.random.random((400,400)))

# set scatter
num_points = 10
xdata=400*np.random.random((num_points))
ydata=400*np.random.random((num_points))
widget.scatterPlot.setData(x=xdata, y=ydata, pen=pg.mkPen(None), brush='g', symbol='x', size=15)

# set rectangle
roi_size = 5
rois_draw = []
for x,y in zip(xdata,ydata):
    roi_temp = pg.PlotCurveItem(x=[x-roi_size/2,x+roi_size/2,x+roi_size/2,x-roi_size/2,x-roi_size/2], y=[y-roi_size/2,y-roi_size/2,y+roi_size/2,y+roi_size/2,y-roi_size/2], pen=pg.mkPen('g', width=2))
    rois_draw.append(roi_temp)

for roi in rois_draw:
    widget.imgVb.addItem(roi)

#for roi in rois_draw:
#    widget.imgVb.removeItem(roi)
#rois_draw = []


sys.exit(app.exec_())