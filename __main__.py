from qtpy import QtWidgets
import sys
import os
import json
from app_paths import ensure_user_tree

if __name__ == "__main__":
    paths = ensure_user_tree()
    home_dir = str(paths["home"])  # Documents\etMINFLUX
    if home_dir not in sys.path:
        sys.path.insert(0, home_dir)
        
    etMINFLUXapp = QtWidgets.QApplication(sys.argv)

    with open(paths["setup_json"], "r", encoding="utf-8") as f:
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
            print('specpy not loaded or Imspector connection to reached. Continuing in simulation mode.')

    if simulation_mode:
        from EtMINFLUXControllerSim import EtMINFLUXControllerSim as Controller
        from EtMINFLUXWidgetSim import EtMINFLUXWidgetSim as Widget
    else:
        from EtMINFLUXController import EtMINFLUXController as Controller
        from EtMINFLUXWidget import EtMINFLUXWidget as Widget
    widget = Widget(paths=paths)
    controller = Controller(widget, paths=paths)

    widget.show()
    sys.exit(etMINFLUXapp.exec_())
