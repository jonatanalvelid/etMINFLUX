from EtMINFLUXController import EtMINFLUXController
from EtMINFLUXWidget import EtMINFLUXWidget
from qtpy import QtWidgets
import sys

if __name__ == "__main__":
    etMINFLUXapp = QtWidgets.QApplication(sys.argv)
    widget = EtMINFLUXWidget()
    controller = EtMINFLUXController(widget)

    widget.show()
    sys.exit(etMINFLUXapp.exec_())
    