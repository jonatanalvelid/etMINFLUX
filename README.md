# etMINFLUX
Event-triggered MINFLUX controller, written to interact with the Imspector software, in turn controlling a MINFLUX microscope also capable of running confocal imaging. 

## Installation instructions
No installation of the codebase is required; the repository can be cloned/downloaded/copied, and in order to run the widget __main__.py should be run.
The code requires a Python environment (such as a conda environment), has been tested with Python version 3.10, and has the following dependencies:
- specpy (https://pypi.org/project/specpy/, v1.2.3 (see attached wheel file) used during development (works with Python 3.10), v1.2.1 available via pip works with Python 3.6 and has not been tested)
- qtpy (https://pypi.org/project/QtPy/)
- PyQt5
- pyqtgraph
- matplotlib
- tifffile
- numpy
- scipy
- pynput (https://pypi.org/project/pynput/)
- mouse (https://pypi.org/project/mouse/)
Individual implemented pipelines have among the following dependencies:
- pandas
- scikit-image
- opencv (https://pypi.org/project/opencv-python/)
- trackpy (https://pypi.org/project/trackpy/)

## Setup etMINFLUX
To setup the python evniroment needed to run the etMINFLUX software, follow the following instructions.
- Install a virtual environment using the provided etminflux.yml or requirements.txt files in this repository. Conda has been used during development and was used to export these files. 
- Install specpy v1.2.3 using pip from the provided wheel file in the packages folder.
In order to run event-triggered recordings, some settings need to be adjusted in Imspector.
- A keyboard shortcut for setting a marked region as a MINFLUX ROI has to be created, which the code will use by simulating keyboard input. The shortcut should be set to Ctrl+Shift+Alt+m. This has to be done once and will be saved as long as Imspector settings are not reset. If this shortcut is already in use for something else, select another shortcut and adjust the code accordingly; find the relevant code in EtMINFLUXController, line 1312-1320.
- The MINFLUX sequences to be triggered has to be added (replacing) in the list on line 73 in EtMINFLUXController.py.
- Create empty folders dataDir and transformsDir, for data saving and calibration files, and update these folders on line 42-43 in EtMINFLUXController.py. (Recommended: in User\Documents)
- If the specpy package does not find a connection with an open Imspector version, a mock Imspector class will be run that allows testing of the GUI.

## Run an etMINFLUX experiment
- Run the etMINFLUX software by running __main__.py in this folder.

## Publications
The etMINFLUX codebase is used in the following publications:
- Alvelid, Koerfer, Eggeling 2025, "Event-triggered MINFLUX microscopy: smart microscopy to catch and follow rare events", manuscript (bioRxiv: )

## Contact
If you have any questions about the repository, code, want to use it for your own experiments, or anything else, do not hesitate to contact Jonatan Alvelid on GitHub, BlueSky (@jonatanalvelid.bsky.social), by email (jonatan.alvelid(at)scilifelab.se), or open an issue in the repository.
