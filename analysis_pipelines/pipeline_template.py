import numpy as np

def pipeline_template(img, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None, var1=0, var2=1.5, var3=100):
    """ Template etMINFLUX real-time analysis pipeline.
    The provided method variables have to be included in the signature, and are the following:
    img - current confocal image (np.array)
    prev_frames - previous confocal image(s)
    binary_mask - binary mask of the regions to consider
    exinfo - any Python object that should be returned from the function and input to the next pipeline call on the following confocal frame;
    an example of info to put here can be a pandas DataFrame with tracks from previously connected peak detections in previous frames
    presetROIsize - boolean to if the ROI sizes are preset or not
    These variables can be left as None and not used in the pipeline, but will always be inputted to the pipeline call. 
    Following this, any variables can be listed, all with default values (that will be loaded when the pipeline is loaded in the GUI),
    and it has to be numeric variables of float or int object type.
    Any packages can be used and imported above this function.
    The function name should be identical to the filename, for example in this case: pipeline_template in pipeline_template.py
    Optional: if a second channel from a two-color confocal image should be used in the pipeline, add img_ch2 in the second position in the call signature,
    and include 'dualcolor' somewhere in the pipeline name. In this case, the second channel will be used with a special function call in the code, and 
    that image can be used in any way fit in this function. The return signature does not change from the default. 
    
    The function has to return:
    (1) coordinates of any detected events in the form of a numpy array of lists, where each list is a set of (x,y) event coordinates.
    For example, for a return of three events: coordinates = np.array([[x1, y1], [x2, y2], [x3, y3]]) where (xi, yi) is the coordinates of the ith event. 
    (2) roi_sizes, if those are determined in the pipeline; if they are not determined, return False in this position to use pre-set ROI sizes.
    (3) exinfo, see description of input value above. It can be edited in the pipeline, and returned again.
    (4) return_img, which is the processed input image, that will be displayed in the pop-out widget while running experiments, or when running the TestVisualize mode.
    
    See below for an example of how to use the parameters and how to return them.
    See also the other implemented pipelines for examples on how to develop them.
    """
    
    roi_sizes = False
    return_img = img
    x1,x2,x3,y1,y2,y3 = 15,30,60,12,48,100
    coordinates = np.array([[x1, y1],[x2, y2],[x3, y3]])
    
    return coordinates, roi_sizes, exinfo, return_img