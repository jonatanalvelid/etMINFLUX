___
## **Note that this branch works with the latest Imspector version m2410 and its specpy version. For use with earlier Imspector version m2205, see corresponding branch.**
___

# etMINFLUX
Event-triggered MINFLUX controller, written to interact with the abberior Imspector software, in turn controlling an abberior Instruments MINFLUX microscope also capable of running confocal imaging. 

## Versions
etMINFLUX has been developed to work with two separate versions of Imspector controlling a MINFLUX microscope. The main branch in the repository, and the released executables, work with Imspector version m2410 (v16.3.21315), as released by abberior Instruments in 2025. If your MINFLUX microscope runs an older version of Imspector, the branch m2205 contains a version of the etMINFLUX code that was originally developed for that Imspector version. Note that this version requires simulated mouse and keyboard controls, due to those specpy and Imspector versions lacking later-added functionality, and thus requires the user to not adjust the Imspector measurement and widget windows during measurements, and always keeping Imspector as the active window.

## Installation instructions
No installation of the codebase is required; the repository can be cloned/downloaded/copied, and in order to run the widget \_\_main\_\_.py should be run (recommended way). Alternatively, the released executable contains the latest (Feb-2026) update of the etMINFLUX code, and can be run without installing a Python environment. The only requirement is to keep the local specpy package folder (not included in the executable, can be found in your local Imspector installation folder: C:\Imspector\Versions\16.3.xxxx-wxxx-win64-MINFLUX\python\specpy\Python3.xx.x-NumPy1.x.x\specpy) named specpy in the same folder as the executable for proper loading of this functionality (includes specpy and related required DLLs.

The code requires a Python environment (such as a conda environment), has been developed and tested with Python v3.10.4 and Imspector v16.3.15645/15635 (m2205) and v16.3.21315 (m2410), and has the dependencies in the list below. It is recommended to use a conda environment, and install everything from conda-forge (everything available except mouse and opencv) to avoid issues with specpy DLL loadings. 
- specpy (v1.2.3 used during development (works with Python 3.10.4); see below for more info on specpy versions and where to find it. Can be omitted if only requiring a simulated etMINFLUX environment. You can simply import it or use the wheel to install it into your python env.)
- qtpy
- pyqt (required =5)
- pyqtgraph
- matplotlib
- tifffile
- numpy (<v2 for specpy, recommended =1.25)
- scipy
- pynput
- mouse (only available from pip: https://pypi.org/project/mouse/)

Individual implemented pipelines have among the following dependencies:
- pandas
- skimage
- opencv (only available from pip, use opencv-python-headless: https://pypi.org/project/opencv-python-headless/)
- trackpy

## Specpy versions
The publicly available v1.2.1 (2017; https://pypi.org/project/specpy/) does not work with etMINFLUX, as it does not contain the required functionality for controlling a MINFLUX microscope and acquisitions. The publicly available specpy version is not expected to be updated by abberior. Instead, the later versions of Imspector ships with a complementary specpy version in the installation folders. Always use the specpy version shipped with your specific version of Imspector. The installation wheel file for specpy can be found in the local Imspector folder upon installation of Imspector: C:\Imspector\Versions\16.3.xxxx-wxxx-win64-MINFLUX\python\specpy\Python3.xx.x-NumPy1.x.x\specpy-x.x.x-cpxxx-cpxxx-win_amd64.whl (x replaces version numbers). 

The continuously shipped specpy versions with newer versions of Imspector shows abberiors commitment to maintain the compatibility, also since latest versions contain updated MINFLUX functionality through specpy, which means that this framework implementation of etMINFLUX will remain functional as long as Imspector is the main control software for abberior MINFLUX microscopes. Should this no longer be the case, the backbone of the etMINFLUX widget will still be possible to upgrade by changing of the specpy calls in the code to any new Python interface that implements similar functionality. 

## Setup etMINFLUX
Upon launching etMINFLUX for the first time, a .../Documents folder will be created, called etMINFLUX. This contains the subfolders config_files, analysis_pipelines, and transform_pipelines, with files used by the program.

If you run the executable, note that it may take ~20-30 s for it to start.

etMINFLUX settings required to be adjusted by the user can all be found in .../Documents/etMINFLUX/config_files/etMINFLUX_setup.json. The code reads this file, and according to the instructions below uses the field values to adjust the etMINFLUX widget to the local MINFLUX setup. If a simulation version of etMINFLUX is requested, adjust the system_simulation parameter in etMINFLUX_setup.json. In this way, etMINFLUX can be run without access to a local specpy and MINFLUX microscope, and allows testing of analysis pipelines by loading pre-recorded confocal timelapse data in tiff stacks. See below for further information about this mode.

To setup the python environment needed to run the etMINFLUX software as a script (running __ main __.py), follow the following instructions. If using the executable, simply run the executable from a folder in which you also have specpy in a folder called specpy.
- Install a virtual environment using the provided etminflux.yml file in this repository. Conda has been used during development and was used to export these files, as such a conda environment is recommended.
- Additionally, in the new environment, either install specpy v1.2.3 (or any matched specpy version >1.2.3 for your version of Imspector) or keep a folder with specpy in the repository folder (running as script) or in the same folder as the .exe (running as executable). An installation wheel file can be found in the local Imspector folder upon installation of Imspector: C:\Imspector\Versions\16.3.xxxx-wxxx-win64-MINFLUX\python\specpy\Python3.xx.x-NumPy1.x.x\specpy-x.x.x-cpxxx-cpxxx-win_amd64.whl (x replaces version numbers), where also the folder can be found in case that method is preferred: C:\Imspector\Versions\16.3.xxxx-wxxx-win64-MINFLUX\python\specpy\Python3.xx.x-NumPy1.x.x\specpy. Copy the wheel file to your cloned repository, change folder to your repository in the command line, and install it using pip: pip install specpy-x.x.x-cpxxx-cpxxx-win_amd64.whl (x replaces version numbers).
- If a simulation version of etMINFLUX is requested (also for executable), adjust the system_simulation parameter in etMINFLUX_setup.json (true = simulated, false = MINFLUX microscope and specpy required).
- Total installation time on a MINFLUX control computer should be on the order of a few minutes.
  
In order to run event-triggered recordings, some settings need to be adjusted in Imspector and etMINFLUX_setup.json (<5 min).
- Start Imspector.
- The IDs/names of the MINFLUX sequences to be triggered has to be added in (replacing) the list in acquisition_settings/minflux_seqs in etMINFLUX_setup.json.
- Adjust/add the names of MINFLUX excitation lasers, activation lasers, and detectors in the the lists in hardware_settings/mfx_exc_lasers, mfx_act_lasers, mfx_detectors_gui mfx_detectors_imspector, and mfx_detectors_imspector_threads in etMINFLUX_setup.json. The names in mfx_exc_lasers, mfx_act_lasers, and mfx_detectors_gui have to match the corresponding names in the Imspector GUI (for detectors, this is single detectors or combination of detectors). The names/lists in mfx_detectors_imspector are backend IDs in Imspector for the detectors, and are named as DETX where X is an ID starting at 1, for the detector that is on top of the list of the GUI, and increases by one for each subsequent detector in the list of detectors of the GUI. The names/lists in mfx_detectors_imspector_threads are yet another naming used in the acquisition threads functionality in Imspector, and this list should match the mfx_detectors_imspector list with the difference being that in mfx_detectors_imspector_threads even the single entries (such as DET1 in the example .json provided) should be list elements in a string ("DET1" -> "[\"DET\"]"). See the example .json and follow the same pattern, or leave untouched if your system has the same detectors available.
- Set local data saving output directory in save_settings/save_directory in etMINFLUX_setup.json (recommended: in User\Documents).
- Start Imspector.
- If the specpy package does not find a connection with an open Imspector version, a mock Imspector class will be run that allows testing of the GUI.
- Crete a template for etMINFLUX, and it is recommended to always use the sample template to avoid missing some settings. In this template, have one configuration (as can be seen in the Imspector Configurations widget) and connected measurement window where your confocal images to detect events in will be recorded. The configuration name has to match the name you provide in acquisition_settings/confocal_config in etMINFLUX_setup.json. The default name is "confocal". The MINFLUX configuration should by default be named "Minflux", as in the default etMINFLUX settings, but if any discrepancy is present, you can adjust this configuration name in acquisitions_settings/minflux_config in etMINFLUX_setup.json.

Additionally, small adjustments have to be made in the etMINFLUX software to calibrate it to the current running system (<1 min). A small GUI calibration has to be performed in order for the software to be able to simulate mouse clicks for a function that is missing in the specpy API. Note: this is only required if using the automatic MINFLUX dataset deletion upon saving option (unused by default). While running an experiment with this setting turned on, moving the Imspector window, rearranging the widgets in Imspector or moving the mouse during these occasions may result in non-deleted datasets (clicks will not happen where they are supposed to happen). 
- In the etMINFLUX widget, press the "Calibration settings..." button. In the pop-out widget, the calibration value saved in gui_settings/top_mfx_dataset_position in etMINFLUX_setup.json is automatically loaded. If calibration has to be performed, do the following:
  - Press the "Top MFX dataset" button, and then click the line where the top MINFLUX dataset would show up in the Workspace widget in Imspector. The default value can be adjusted in the etMINFLUX_setup.json to the new value that shows on the button after calibration, if always using Imspector and the Workspace widget in the same screen positions.
- Timing settings (timing_settings) in etMINFLUX_setup.json controls the time in seconds of different sleep events that take place during etMINFLUX experiments in order to ensure that the Imspector communication runs consistently. With too small values certain specpy call combinations do not run safely at all times, with occasional functionality skips where the function has been called through specpy but Imspector not performing the requested action, and without feedback back through specpy. The default values have been set conservatively to allow the code to safely run on all potential control PCs. Lowering the values to save a few seconds during measurements might be possible depending on the specifications of the PC used to run Imspector (likely on most PCs). During testing, the three timing values have been used at 0.3, 1.0, and 3.0.
- The data saving directory can be set in the GUI using the option under the calibration pop-out widget.

## Setup etMINFLUX experiment settings
Before starting an etMINFLUX experiment, all acquisition settings pertaining to the confocal timelapse, analysis pipeline, and MINFLUX acquisitions have to be chosen. Below all options to change parameters are listed categorically according to their different subgroups in the widget GUI.
- Under the Analysis pipeline title you find the settings for which real-time analysis pipeline to use, and which parameter values to use for that specific pipeline.
  - Choose the real-time analysis pipeline to apply to the confocal images, using the dropdown menu under the "Analysis pipeline" heading. Press Load pipeline.
  - The number of confocal channels required by the pipeline will show up as the first line below the selection. The system supports an arbitrary number of confocal channels, but the pipeline has to take that many confocal input channels, and the microscope configuration has to be set up to match that number of simultaneous confocal channels.
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
  - If an activation laser is to be used, the laser name must be present in hardware_settings/mfx_act_lasers in etMINFLUX_setup.json. If so, the activation power can be set using the then activated field in the GUI and will be turned on and set to that value upon MINFLUX acquisition start.
  - The ROI size, if pre-set (default) and not decided from the analysis pipeline, check the "Pre-set ROI size" checkbox and set the corresponding X and Y field values to wanted size in µm.
  - The ROI recording time can be set in the MFX ROI rec time editable field, in seconds, if the Pre-set ROI rec time checkbox is checked. If not checked, the measurement will run indefinitely until Stop is pressed in the EtMINFLUX widget (Initiate button while running experiment).
  - If all the detected events/ROIs should be run sbusequently, and not a single one, check the checkbox Trigger all ROIs. In this case, after the set MFX ROI rec time has expired, the EtMINFLUX widget will automatically stop the current ROI recording and start a new one in the next area, until all the ROIs have been recorded. After this, it will return to stop the full experiment, or run a new confocal image/timelapse (if Endless mode is on).
  - If the user wants to record at Random ROI sites, check the corresponding box. In this case, after a single confocal frame, random ROI positions will be chosen and the MINFLUX measurements will be run at these positions. The analysis pipeline used does not matter in this case. The areas to use for randomizing the site can be set using the Binary mask control in the separate pop-out widget (see below).
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
- OPTIONAL For setting a binary mask from the confocal images, to be used to mask the area that the analysis pipeline should consider and/or when recording random ROI positions, follow the following instructions.
  - Open the binary mask settings pop-out widget by pressing the Binary mask settings... button.
  - Set the binary border size, outside which everything will be marked as not being considered.
  - Set the binary smoothing radius for the positive mask, and the binary positive intensity threshold above which to mark the areas as being considered (for excluding background).
  - Set the binary smoothing radius for the negative mask, and the binary negative intensity threshold below which to mark the areas as being considered (for excluding bright spots and areas).
  - The two binary masks (positive and negative) will be multiplied, and only the areas covered by both will be considered. The resulting binary mask will be shown in a pop-out widget.
  - If a binary mask has been calculated, it will automatically be used. If the user wishes to no longer use it, use the Reset binary mask button. Note: a new mask or resetting of the current mask is required when changing the confocal image size in pixels.

## Run an etMINFLUX experiment
After all the settings listed above have been set to the liking of the user, the experiment phase can begin. An experiment will run for as long as the user defines the experiment, but can last between <10 s up to hours, depending on the application. Follow the following instructions to run an experiment.
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

## Simulation running mode
EtMINFLUX version m2410 (main branch) implements a new running mode - simulation mode. Through this mode, launched if no specpy connection can be made, the user can gain familiarity to the GUI even in an office setting without access to a MINFLUX microscope. This running mode can also, most importantly, be used to implement, test, and optimize real-time analysis pipelines. Pre-recorded confocal timelapses can be loaded in the GUI, and then, by running etMINFLUX in visualization mode, according to what is described above, the analysis pipeline choosen will run on the loaded dataset as if it was on a microscope. It will stop when the pre-loaded timelapse finishes. Analysis pipeline parameters can then be optimized to real microscope data without access to a microscope, to prepare for real experiments in a laboratory setting, optimizing the whole experimental pipeline. Pre-recorded confocal timelapses can be single, dual, or triple color datasets, and all pipelines supported in a real experimental setting are supported in the simulation mode. Make sure that pre-loaded confocal datasets are in tif format. In order to pre-load confocal multicolor datasets, make sure that the tif files loaded have been saved as ImageJ hyperstacks with the correct dimensionality in terms of channels (colors) and time. 

## Implement new real-time analysis pipelines
The analysis pipelines provided here are a set of pipelines that were developed during the development of the etMINFLUX method. New experiment-specific pipelines can easily be implemented and added as alternatives by following the following instructions, and by looking at the pipeline_template.py provided in the defaults/analysis_pipelines (or .../Documents/etMINFLUX/analysis_pipelines upon running etMINFLUX for the first time as decribed above) folder with additional instructions.
- Start from the pipeline_template.py file. Copy the file, rename the function and rename the file to the same name. (function in function.py)
- Place the file in the .../Documents/etMINFLUX/analysis_pipelines folder, after it will automatically show up in the etMINFLUX widget (restart needed).
- Add any number of additional variables specific to that pipeline, all with numerical values, to the default variables that should not be removed (img, prev_frames, binary_mask, exinfo, presetROIsize). The additional variables will show up in the GUI when loading that analysis pipeline.
- Use any python package, and develop your analysis pipeline in the function. If you require different python modules than in the default pipelines, install them in the environment in which you run etMINFLUX. 
- Follow the return signature as in pipeline_template. Return any detected event coordinates in the form of a numpy array of lists of (x,y) event coordinate pairs. Even with only one detected event, the structure has to be the same, i.e.: np.array([[x1,y1]]).

## Publications
The etMINFLUX codebase is used in the following publications:
- Alvelid, Koerfer, Eggeling 2025, "Event-triggered MINFLUX microscopy: smart microscopy to catch and follow rare events", manuscript (bioRxiv: https://doi.org/10.1101/2025.08.27.672674)

## Contact
If you have any questions about the repository, code, want to use it for your own experiments, or anything else, do not hesitate to contact Jonatan Alvelid on GitHub, BlueSky (@jonatanalvelid.bsky.social), by email (jonatan.alvelid(at)scilifelab.se), or open an issue in the repository.
