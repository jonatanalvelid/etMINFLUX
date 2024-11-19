import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
from scipy import ndimage as ndi
import cv2


def dyn_signalrise_peakdetcheck(img, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None,
                     min_dist=1, num_peaks=1000, thresh_abs_lo=1.1, thresh_abs_hi=15):

    """
    Common parameters:
    img - current image,
    prev_frames - previous image(s)
    binary_mask - binary mask of the region to consider
    testmode - to return preprocessed image or not
    exinfo - pandas dataframe of the detected vesicles and their track ids from the previous frames

    Pipeline specific parameters:
    min_dist - minimum distance in pixels between two peaks
    num_peaks - number of peaks to track
    thresh_abs_lo - low intensity threshold in img_ana of the peaks to consider
    thresh_abs_hi - high intensity threshold in img_ana of the peaks to consider
    """
    
    roi_sizes = False

    # define non-adjustable parameters
    smoothing_radius = 1.5
    dog_lo = 0.05
    dog_hi = 3
    
    # gaussian filter raw image
    img = np.array(img).astype('float32')
    img_filt = ndi.gaussian_filter(img, smoothing_radius)

    # difference of gaussians to get clear peaks separated from spread-out bkg and noise
    img_dog_lo = ndi.gaussian_filter(img_filt, dog_lo)
    img_dog_hi = ndi.gaussian_filter(img_filt, dog_hi)
    img_dog = img_dog_lo - img_dog_hi
    img_dog = np.clip(img_dog, a_min=0, a_max=30000)

    # further filtering to get a better image for peak detection
    img_ana = img_dog.astype('float32')
    img_ana = ndi.gaussian_filter(img_ana, smoothing_radius)  # Gaussian filter the image, to remove noise and so on, to get a better center estimate
    # Peak_local_max all-in-one as a combo of opencv and cupy
    # get filter structuring element
    size = int(2 * min_dist + 1)
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