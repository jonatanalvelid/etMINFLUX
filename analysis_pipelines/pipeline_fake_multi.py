import numpy as np

def pipeline_fake_multi(img_ch1, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None, xcoord1=30, ycoord1=30, xcoord2=60, ycoord2=60, xcoord3=100, ycoord3=100):
    """ Fake analysis pipeline, that just returns the coordinates of three points based on the provided x and y coordinates. """
    roi_sizes = False
    coordinates = np.array([[xcoord1, ycoord1],[xcoord2, ycoord2],[xcoord3, ycoord3]])
    if not presetROIsize:
        roi_sizes = [[15,15],[40,40],[25,25]]
    return coordinates, roi_sizes, exinfo, img_ch1