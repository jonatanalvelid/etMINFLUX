import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
from scipy import ndimage as ndi
import cv2


def gag_signalrise_peakdetcheck(img_ch1, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None,
                     num_peaks=300, thresh_abs_lo=1.7, thresh_abs_hi=10):

    """
    Analysis pipeline to be used for the visualize testing mode, that only performs the peak detection,
    allowing tweaking of peak detection parameters.
    
    Common parameters:
    img_ch1 - current image
    prev_frames - previous image(s)
    binary_mask - binary mask of the region to consider
    exinfo - pandas dataframe of the detected vesicles and their track ids from the previous frames

    Pipeline specific parameters:
    num_peaks - number of peaks to track
    thresh_abs_lo - low intensity threshold in img_ana of the peaks to consider
    thresh_abs_hi - high intensity threshold in img_ana of the peaks to consider
    """
    
    roi_sizes = False

    # define non-adjustable parameters
    smoothing_radius_raw = 1.5  # pixels
    
    img_ch1 = np.array(img_ch1).astype('float32')
    img_ana = ndi.gaussian_filter(img_ch1, smoothing_radius_raw)
    # Peak_local_max as a combo of opencv and numpy
    img_ana = np.clip(img_ana, a_min=0, a_max=None)
    img_ana = img_ana.astype('float32')
    # get filter structuring element
    size = int(2 * smoothing_radius_raw)
    footprint = cv2.getStructuringElement(cv2.MORPH_RECT, ksize=[size,size])
    # maximum filter (dilation + equal)
    image_max = cv2.dilate(img_ana, kernel=footprint)
    mask = np.equal(img_ana, np.array(image_max))
    mask &= np.greater(img_ana, thresh_abs_lo)
    mask &= np.less(img_ana, thresh_abs_hi)

    # get coordinates of peaks
    coordinates = np.nonzero(mask)
    intensities = img_ana[coordinates]
    # highest peak first
    idx_maxsort = np.argsort(-intensities)
    coordinates = tuple(arr for arr in coordinates)
    coordinates = np.transpose(coordinates)[idx_maxsort]

    # remove everyhting down to a certain length
    if len(coordinates) > num_peaks:
        coordinates = coordinates[:int(num_peaks),:]

    coordinates = np.flip(coordinates, axis=1)

    return coordinates, roi_sizes, exinfo, img_ana
