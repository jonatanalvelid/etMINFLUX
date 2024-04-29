import numpy as np

def pipeline_fake(img, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None, xcoord=300, ycoord=300):
    roi_sizes = False
    coordinates = np.array([[xcoord, ycoord]])
    if not presetROIsize:
        roi_sizes = [[10,10]]
    return coordinates, roi_sizes, exinfo, img