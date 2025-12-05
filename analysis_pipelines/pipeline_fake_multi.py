import numpy as np

def pipeline_fake_multi(img_ch1, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None, xcoord1=300, ycoord1=300, xcoord2=500, ycoord2=500, xcoord3=400, ycoord3=200):
    """ Fake analysis pipeline, that just returns the coordinates of three points based on the provided x and y coordinates. """
    roi_sizes = False
    coordinates = np.array([[xcoord1, ycoord1],[xcoord2, ycoord2],[xcoord3, ycoord3]])
    if not presetROIsize:
        roi_sizes = [[15,15],[40,40],[25,25]]
    return coordinates, roi_sizes, exinfo, img_ch1