# etMINFLUX
Event-triggered MINFLUX controller, written to interact with the abberior Imspector software, in turn controlling an Abberior MINFLUX microscope also capable of running confocal imaging. 

## Installation instructions
No installation of the codebase is required; the repository can be cloned/downloaded/copied, and in order to run the widget \_\_main\_\_.py should be run. 

The widget controls a commercial Abberior MINFLUX microscope, and thus requires such a system to function.

The code requires a Python environment (such as a conda environment), has been developed and tested with Python v3.10 and Imspector v16.3.15645 on Windows 10, and has the following dependencies:
- specpy (https://pypi.org/project/specpy/, v1.2.3 used during development (works with Python 3.10); v1.2.1 available via pip works with Python 3.6 and has not been tested. See below for more info on specpy version.)
You can simply import it or use the wheel to install it into your python env.)
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
- Additionally, in the new environment, install specpy v1.2.3 (or any matched specpy version for your version of Imspector, where >1.2.3 should work). An installation wheel file can be found in the local Imspector folder upon installation of Imspector: C:\Imspector\Versions\16.3.xxxx-wxxx-win64-MINFLUX\python\specpy\Python3.xx.x-NumPy1.x.x\specpy-x.x.x-cpxxx-cpxxx-win_amd64.whl (x replaces version numbers). Copy the wheel file to your cloned repository, change folder to your repository in the command line, and install it using pip: pip install specpy-x.x.x-cpxxx-cpxxx-win_amd64.whl (x replaces version numbers).
- Total installation time on a MINFLUX control computer should be on the order of a few minutes.
  
In order to run event-triggered recordings, some settings need to be adjusted in Imspector (<5 min).
- Start Imspector.
- A keyboard shortcut for setting a marked region as a MINFLUX ROI has to be created, which the code will use by simulating keyboard input. The shortcut should be set to Ctrl+Shift+Alt+m. This has to be done once and will be saved as long as Imspector settings are not reset. If this shortcut is already in use for something else, select another shortcut and adjust the code accordingly; find the relevant code in EtMINFLUXController, line 1312-1320.
- The IDs/names of the MINFLUX sequences to be triggered has to be added (replacing) in the list on line 43 in EtMINFLUXController.py.
- If other MINFLUX lasers than the default ones should be used, or if the laser list indexes are different than the default ones (MINFLUX 488, MINFLUX 560, MINFLUX 642 options, for lasers called MINFLUX 488, MINFLUX 560 and MINFLUX 642 in Imspector; 488 = idx 1 (zero-indexing), 560 = idx 3, 642 = idx 6), the lists of laser names (only for display in etMINFLUX) and corresponding list indexes in the Channels widget in Imspector (important for correct automatic selection and power setting) should be adjusted in lines 52 and 53 in EtMINFLUXController.py.
- Create empty folders dataDir and transformsDir, for data saving and calibration files, and update these folders on line 45-46 in EtMINFLUXController.py (recommended: in User\Documents).
- If the specpy package does not find a connection with an open Imspector version, a mock Imspector class will be run that allows testing of the GUI.
- Crete a template for etMINFLUX, and it is recommended to always use the sample template to avoid missing some settings. In this template, have one configuration (as can be seen in the Imspector Configurations widget) and connected measurement window where your confocal images to detect events in will be recorded. The configuration name has to match a setting in the etMINFLUX code, line 55 in EtMINFLUXController.py. The default name is "ov conf". If your confocal configuration has a different name, either rename it to ov conf, or adjust the name in the etMINFLUX code to match the current configuration name.

Additionally, small adjustments have to be made in the etMINFLUX software to calibrate it to the current running system (<1 min). Some GUI calibration has to be performed in order for the software to be able to simulate mouse click and drags for the functions that are missing in the specpy API. Important: while running an experiment, as simulated mouse movements are used to mark and set a MINFLUX ROI and start a measurement, moving the Imspector window, rearranging the measurement window in Imspector or moving the mouse during these occasions may result in failure of the experiments (clicks and drags will not happen where they are supposed to happen). 
- Run the etMINFLUX software by running \_\_main\_\_.py in this folder.
- In the widget that shows, press the Transform calibration button. In the new pop-out widget that shows up, the last saved calibration is automatically loaded. If calibration has to be performed, do the following:
  - (1) click "Click detect: confocal top left pixel button" and then click the top left pixel in your confocal image window used for event detection in Imspector. Do the same for "Click detect: confocal bottom right pixel". The corresponding values to the left in the pop-out widget will change to the pixels you clicked.
  - (2) Save calibration to be able to use the same next time you open the software.
  - (3) Press "Rep. meas. calib." button, and then click the Auto Repetition button in Imspector (looping measurement scan). The default value can be adjusted on line 48 in EtMINFLUXController.
  - (4) Press "Top MFX dataset" button, and then click the line where the top MINFLUX dataset would show up in the Workspace widget in Imspector. The default value can be adjusted on line 50 in EtMINFLUXController.
- Remaining calibration fields (Timing and Distances) are not necessary to change, but can be adjusted upon user disgression. Timing fields control the time in seconds of different sleep events that take place during etMINFLUX experiments in order to ensure that the simulated mouse and keyboard events are consistent. Lowering the values to save a few seconds during measurements might be possible depending on the specs of the PC used to control the microscope (default values are tested to be valid on Abberior MINFLUX control computers).
- The save directory should be set using the option under the Saving title in the calibration pop-out widget. Here, all files will be saved after an etMINFLUX recording has finished.

## Setup etMINFLUX experiment settings
Before starting an etMINFLUX experiment, all acquisition settings pertaining to the confocal timelapse, analysis pipeline, and MINFLUX acquisitions have to be chosen. Below all options to change parameters are listed categorically according to their different subgroups in the widget GUI.
- Under the Analysis pipeline title you find the settings for which real-time analysis pipeline to use, and which parameter values to use for that specific pipeline.
  - Choose the real-time analysis pipeline to apply to the confocal images, using the dropdown menu under the "Analysis pipeline" heading. Press Load pipeline.
  - Adjust the pipeline parameter values to your liking (after previous testing and optimization (see below)). See each provided pipeline file for information on what each parameter does.
- Under the Analysis control title, settings regarding how the analysis is run can be set.
  - The anlysis defaults to running on every frame that is collected, but can instead be run on the image after N lines have been recorded - set this using the Analysis period (lines) field and check the "Run analysis pipeline linewise" checkbox.
  - While running experiments, clicking the checkbox "Plot ROI (experiment mode)" before starting the experiment will allow the pop-out widget showing the latest confocal frame and detected events to be shown during real experiments.
  - The confocal imaging defaults to no pause between frames. If a pause is required in the confocal timelapse, check the "Confocal frame pause" checkbox and set the wanted pause between the end of one frame and the start of the next in seconds using the field.
  - The Initial frames field controls a value that allow you to not run analysis on M initial confocal frames, by setting a non-zero value in the field.
  - The two lines below shows information during experiments regarding the time until the next confocal frame (if using a confocal frame pause), and the number of confocal frames acquired so far since the start or since the last event (in endless mode).
- Under the MINFLUX imaging parameters title, settings regarding the triggered MINFLUX imaging can be set.
  - The MINFLUX sequence to be used can be chosen from the dropdown menu (see above of how to set available sequences in the code).
  - The MINFLUX excitation laser to be used can be chosen from the corresponding dropdown menu (see above how to set available lasers in the code).
  - The MINFLUX excitation power to be used can be set in the corresponding field.
  - If a 405 activation laser (default laser index 0) is to be used, the corresponding field can be enabled in line 115 in EtMINFLUXWidget.py. If so, the activation power can be set using this field.
  - The ROI size, if pre-set (default) and not decided from the analysis pipeline, check the "Pre-set ROI size" checkbox and set the corresponding X and Y field values to wanted size in µm.
  - The ROI recording time can be set in the MFX ROI rec time editable field, in seconds, if the Pre-set ROI rec time checkbox is checked. If not checked, the measurement will run indefinitely until Stop is pressed in the EtMINFLUX widget (Initiate button while running experiment).
  - If all the detected events/ROIs should be run sbusequently, and not a single one, check the checkbox Trigger all ROIs. In this case, after the set MFX ROI rec time has expired, the EtMINFLUX widget will automatically stop the current ROI recording and start a new one in the next area, until all the ROIs have been recorded. After this, it will return to stop the full experiment, or run a new confocal image/timelapse (if Endless mode is on).
  - If the user wants to record at Random ROI sites, check the corresponding box. In this case, after a single confocal frame, random ROI positions will be chosen and the MINFLUX measurements will be run at these positions. The analysis pipeline used does not matter in this case. The areas to use for randomizing the site can be set using the Binary mask control under the Binary mask title.
- Under the Saving title, two options can be controlled.
  - The Auto-save .msr after event option will automatically save the .msr file in Imspector, with the MINFLUX and latest confocal frame data, after the event has finished recording.
  - The Auto-delete MFX data after save will automatically delete the MINFLUX datasets after saving the .msr file, so that these are not saved again when the next event has finished recording. This deletion is also GUI sensitive, and the position of the Workspace widget in Imspector should not change after calibrating its position (see above), as the option to delete datasets is lacking in the specpy API.
  - If auto-save is disabled, the button to save the current measurement can be pressed instead of manually doing this in Imspector.
- Under the ROI following mode title, settings regarding if to use a ROI following mode can be controlled (see publication for further details on these modes).
  - Check the checkbox if a ROI following mode should be used. These modes record multiple MINFLUXs recordings at a single detect event site, with intermittent confocal images between.
  - Decide which following mode to use:
    - "Single" follows a single detect event site with MINFLUX, confocal, MINFLUX,....
    - "SingleRedetect" follows a single detect event site, but re-runs the pipeline on every intermittent confocal image. It checks the distance of the previous ROI position to all the events detected in the latest pipeline run, and if any distance is shorter than that noted in the "Redetect dist threshold" field, those events are linked, and the next MINFLUX ROI will be recorded at the new event position. This is implemented for following bright peaks that move and are detected with some type of peak detection pipeline. 
    - "Multiple" follows mutliple detect event sites by recording all MINFLUX ROIs and then an intermittent confocal image as MINFLUX, MINFLUX,..., confocal, MINFLUX, MINFLUX,....
  - The Confocal interval field controls the rate at which to perform intermittent confocal images for single ROI following mode. With the multiple ROI following mode it controls the recording time of each MINFLUX ROI cycle (the MFX ROI rec time field described above does not do anything in this case).
- OPTIONAL For setting a binary mask from the confocal images, following the following instructions (mostly useful when recording random ROIs).
  - Set the binary border size, outside which everything will be marked as not being considered.
  - Set the binary smoothing radius for the positive mask, and the binary positive intensity threshold above which to mark the areas as being considered (for excluding background).
  - Set the binary smoothing radius for the negative mask, and the binary negative intensity threshold below which to mark the areas as being considered (for excluding bright spots and areas).
  - The two binary masks will be multiplied, and only the areas covered by both will be considered. The resulting binary mask will be shown in a pop-out widget.
  - If a binary mask has been used, and no longer should be used, use the Reset binary mask button. A new mask or resetting of the current mask is required when changing the confocal image size in pixels.

## Run an etMINFLUX experiment
After all the settings listed above have been set to the liking of the user, the experiment phase can begin. An experiment will run for as long as the user defines the experiment, but can last between <10 s up to hours, depending on the application.
- Setup the type of confocal recording you want to perform in a measurement in Imspector, including scanning parameters, lasers, etc. Make sure to not move the actual measurement window after having performed the calibration described above. If so, recalibrate the position of the measurement window. Also ensure that the confocal measurement configuration is called the same as in the etMINFLUX code (see instructions above for how to set this up). 
- Ensure that the Experiment mode dropdown menu is set to Experiment.
- Choose to check the Endless checkbox if you want the experiment to run indefinitely, continuously finding new events after the previous one has finished recording and been saved. 
- Press the initiate button "Initiate etMINFLUX" to start the experiment. Now, the confocal recording will be started and the confocal imaging will be run according to the settings setup in Imspector and the etMINFLUX widget.
- After every confocal analysis period (default 1 frame), the real-time analysis pipeline will be run. If no events are detected, another confocal frame will be run and the process will be repeated.
- If any events are detected, the confocal imaging will be stopped, a ROI will be drawn on the event coordinates (emulated mouse movement), the ROI will be set to a MINFLUX ROI (emulated mouse movement), MINFLUX acquisition parameters and sequence will be set according to your settings in the etMINFLUX widget, the MINFLUX recording time will be set to the value decided in the etMINFLUX widget (0=indefinite if no pre-set rec time is provided), and the MINFLUX measurement will be started.
- The MINFLUX measurement will be run until the MINFLUX timer hits the pre-set rec time, or until the Intiate button is pressed again (Stop), after which another confocal frame will be run. The full Imspector measurement will be saved in an .msr file, including the MINFLUX datasets and last confocal frame. Earlier confocal frames, up until event detection, will be saved in a tiff stack, and any processed confocal images that are returned from the pipeline will also be saved in a tiff stack. A log .txt file will also be saved, including multiple timestamps during the experiment, MFX ROI position, analysis pipeline and pipeline parameters used, analysis pipeline runtimes, and additional metadata information that can be useful in analysis of the event data afterwards. All files will be saved with current timestamps in the name for easy sorting of the data.
- MINFLUX datsets will automatically be names with a timestamp, ROI index, and ROI position in the name. In ROI Follow mode recordings, an additional stamp for which cycle in the following mode it belongs to will be added. 

## Optimize pipeline parameter values
For pipelines such as peak detection pipelines, parameters can be optimized while running testing runs on the microscope.
- (1) Set the experiment mode to TestVisualize, using the dropdown menu to the top right.
- (2) Run experiment as normal, without endless mode selected.
- (3) In the pop-out widget that shows up, the latest confocal image will show. Additionally, all detected events will show up as green crosses in the image, and will be listed in the list widget on the righthand side of the image.
- (4) To adjust pipeline parameter values, stop pipeline, adjust values, and run it again. Repeat until wanted events are detected.
Other types of analysis pipelines might be adjusted in similar ways, or using offline optimization by running the pipeline on a stack of pre-recorded confocal images with the wanted acquisition settings and frame rate.

## Implement new real-time analysis pipelines
The analysis pipelines provided here are a set of pipelines that were developed during the development of the etMINFLUX method. New experiment-specific pipelines can easily be implemented and added as alternatives by following the following instructions, and by looking at the pipeline_template.py provided in the analysis_pipelines folder with additional instructions.
- Start from the pipeline_template.py file. Copy the file, rename the function and rename the file to the same name. (function in function.py)
- Place the file in the analysis_pipelines folder, after it will automatically show up in the etMINFLUX widget.
- Add any number of additional variables specific to that pipeline, all with numerical values, to the default variables that should not be removed (img, prev_frames, binary_mask, exinfo, presetROIsize). The additional variables will show up in the GUI when loading that analysis pipeline.
- Use any python package, and develop your analysis pipeline in the function.
- Follow the return signature as in pipeline_template. Return any detected event coordinates in the form of a numpy array of lists of (x,y) event coordinate pairs. Even with only one detected event, the structure has to be the same, i.e.: np.array([[x1,y1]]).

## Publications
The etMINFLUX codebase is used in the following publications:
- Alvelid, Koerfer, Eggeling 2025, "Event-triggered MINFLUX microscopy: smart microscopy to catch and follow rare events", preprint (bioRxiv: https://doi.org/10.1101/2025.08.27.672674)

## Contact
If you have any questions about the repository, code, want to use it for your own experiments, or anything else, do not hesitate to contact Jonatan Alvelid on GitHub, BlueSky (@jonatanalvelid.bsky.social), by email (jonatan.alvelid(at)scilifelab.se), or open an issue in the repository.
