import numpy as np

def aa_pipeline_fake_multi(img, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None, xcoord1=300, ycoord1=300, xcoord2=500, ycoord2=500, xcoord3=400, ycoord3=200):
    roi_sizes = False
    coordinates = np.array([[xcoord1, ycoord1],[xcoord2, ycoord2],[xcoord3, ycoord3]])
    if not presetROIsize:
        roi_sizes = [[15,15],[40,40],[25,25]]
    return coordinates, roi_sizes, exinfo, img