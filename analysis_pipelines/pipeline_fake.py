import numpy as np

def pipeline_fake(img_ch1, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None, xcoord=135, ycoord=135):
    """ Fake analysis pipeline, that just returns the coordinates of a single point based on the provided x and y coordinates. """
    roi_sizes = False
    coordinates = np.array([[xcoord, ycoord]])
    if not presetROIsize:
        roi_sizes = [[10,10]]
    return coordinates, roi_sizes, exinfo, img_ch1