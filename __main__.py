from EtMINFLUXController import EtMINFLUXController
from EtMINFLUXWidget import EtMINFLUXWidget
from EtMINFLUXControllerSim import EtMINFLUXControllerSim
from EtMINFLUXWidgetSim import EtMINFLUXWidgetSim
from qtpy import QtWidgets
import sys
import json

if __name__ == "__main__":
    etMINFLUXapp = QtWidgets.QApplication(sys.argv)

    with open('etMINFLUX_setup.json', 'r') as f:
        setup_dict = json.load(f)
    if setup_dict.get('system_simulation') == True:
        simulation_mode = True
        print('Simulation mode requested. Continuing in simulation mode.')
    else:
        # test imspector connection
        try:
            import specpy
            imspector = specpy.get_application()
            simulation_mode = False
        except:
            simulation_mode = True
            print('Imspector not loaded. Continuing in simulation mode.')

    if simulation_mode:
        print("load simulation")
        widget = EtMINFLUXWidgetSim()
        controller = EtMINFLUXControllerSim(widget)
    else:
        print("load real system")
        widget = EtMINFLUXWidget()
        controller = EtMINFLUXController(widget)

    widget.show()
    sys.exit(etMINFLUXapp.exec_())
