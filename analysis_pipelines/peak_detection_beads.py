import numpy as np
from scipy import ndimage as ndi
import cv2

def peak_detection_beads(img, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None,
                       maxfilter_kersize=5, thresh_abs=10, smoothing_radius=1, init_smooth=1, 
                       border_limit=15, coord_num_lim_lo=10, coord_num_lim_hi=25):
    """
    Common parameters:
    img - current image,
    prev_frames - previous image(s)
    binary_mask - binary mask of the region to consider
    testmode - to return preprocessed image or not
    exinfo - pandas dataframe of the detected vesicles and their track ids from the previous frames

    Pipeline specific parameters:
    maxfilter_kersize - size of kernel for maximum filtering
    thresh_abs - low intensity threshold in img_ana of the peaks to consider
    smoothing_radius - diameter of Gaussian smoothing of img_ana, in pixels
    border_limit - how much of the border to remove peaks from in pixels
    init_smooth - if to perform an initial smoothing of the raw image or not (bool 0/1)
    """
    roi_sizes = False

    if binary_mask is None or np.shape(binary_mask) != np.shape(img):
        binary_mask = np.ones(np.shape(img)).astype('uint16')

    img = np.array(img).astype('float32')
    if len(prev_frames) != 0:
        prev_frame = np.array(prev_frames[0]).astype('float32')
    else:
        prev_frame = np.zeros(np.shape(img)).astype('float32')
    if init_smooth==1:
        img = ndi.gaussian_filter(img, smoothing_radius)
        prev_frame = ndi.gaussian_filter(prev_frame, smoothing_radius)

    # multiply with binary mask
    img_ana = img * np.array(binary_mask)

    # Peak_local_max as a combo of opencv and numpy
    size = int(2 * maxfilter_kersize + 1)
    img_ana = np.clip(img_ana, a_min=0, a_max=None)
    img_ana = img_ana.astype('float32')
    # get filter structuring element
    footprint = cv2.getStructuringElement(cv2.MORPH_RECT, ksize=[size,size])
    # maximum filter (dilation + equal)
    image_max = cv2.dilate(img_ana, kernel=footprint)
    mask = np.equal(img_ana, np.array(image_max))
    mask &= np.greater(img_ana, thresh_abs)
    
    # get coordinates of peaks
    coordinates = np.nonzero(mask)
    intensities = img_ana[coordinates]
    # highest peak first
    idx_maxsort = np.argsort(-intensities)
    coordinates = tuple(arr for arr in coordinates)
    coordinates = np.transpose(coordinates)[idx_maxsort]

    # remove everything on the border (takes ~2-3ms if there are a lot of detected coordinates, but usually this is not the case)
    imsize = np.shape(img)[0]
    idxremove = []
    for idx, coordpair in enumerate(coordinates):
        if coordpair[0] < border_limit or coordpair[0] > imsize - border_limit or coordpair[1] < border_limit or coordpair[1] > imsize - border_limit:
            idxremove.append(idx)
    coordinates = np.delete(coordinates,idxremove,axis=0)

    # only keep coords inside provided thresholds, after sorting based on intensity
    if len(coordinates) > coord_num_lim_hi:
        coordinates = coordinates[int(coord_num_lim_lo):int(coord_num_lim_hi),:]
    elif len(coordinates) > coord_num_lim_lo:
        coordinates = coordinates[int(coord_num_lim_lo):,:]

    coordinates = np.flip(coordinates, axis=1)

    return coordinates, roi_sizes, exinfo, img_ana
