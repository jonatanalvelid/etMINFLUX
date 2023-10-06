import numpy as np

def aa_pipeline_fake(img, prev_frames=None, binary_mask=None, testmode=False, exinfo=None, xcoord=300, ycoord=300):
    coordinates = np.array([[xcoord, ycoord]])
    roi_sizes = [[10,10]]
    if testmode:
        return coordinates, roi_sizes, exinfo, img
    else:
        return coordinates, roi_sizes, exinfo