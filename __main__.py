import os, sys

def app_dir() -> str:
    return os.path.dirname(sys.executable) if getattr(sys, "frozen", False) else os.path.dirname(os.path.abspath(__file__))

APP_DIR = app_dir()
SPECPY_DIR = os.path.join(APP_DIR, "specpy")

USE_EXTERNAL_SPECPY = os.path.isdir(SPECPY_DIR)

if USE_EXTERNAL_SPECPY:
    # external-folder mode: import specpy from APP_DIR\specpy
    if APP_DIR not in sys.path:
        sys.path.insert(0, APP_DIR)

    os.add_dll_directory(SPECPY_DIR)
    os.environ["PATH"] = SPECPY_DIR + os.pathsep + os.environ.get("PATH", "")
    print(f"Using external specpy from {SPECPY_DIR}")
else:
    # env-installed mode: do nothing special, rely on normal import resolution
    print("Trying to use specpy from environment (no external folder found)")
    pass

os.environ["QT_LOGGING_RULES"] = "*.warning=false"
from qtpy import QtWidgets
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
            try:
                imspector = specpy.get_application()
                simulation_mode = False
            except:
                simulation_mode = True
                print('Imspector connection not reached. Continuing in simulation mode.')
        except:
            simulation_mode = True
            print('specpy not loaded. Continuing in simulation mode.')


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
