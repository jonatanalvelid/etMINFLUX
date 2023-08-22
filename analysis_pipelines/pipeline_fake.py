import numpy as np

def pipeline_fake(img, prev_frames=None, binary_mask=None, testmode=False, exinfo=None, xcoord=300, ycoord=300):
    coordinates = np.array([[xcoord, ycoord]])
    if testmode:
        return coordinates, exinfo, img
    else:
        return coordinates, exinfo