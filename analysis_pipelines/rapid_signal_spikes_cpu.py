import numpy as np
import scipy.ndimage as ndi
import cv2
from scipy.spatial import cKDTree, distance

def rapid_signal_spikes_cpu(img, prev_frames=None, binary_mask=None, testmode=False, exinfo=None,
                            min_dist=30, thresh_abs=0.3, num_peaks=10, noise_level=5,
                            smoothing_radius=1, ensure_spacing=0, border_limit=10,
                            init_smooth=1):
    """
    Common parameters:
    img - current image,
    prev_frames - previous image(s)
    binary_mask - binary mask of the region to consider
    testmode - to return preprocessed image or not
    exinfo - pandas dataframe of the detected vesicles and their track ids from the previous frames

    Pipeline specific parameters:
    min_dist - minimum distance in pixels between two peaks
    thresh_abs - low intensity threshold in img_ana of the peaks to consider
    num_peaks - number of peaks to track
    noise_level - noise level of not wanted signal, to avoid peak detection of potentially high ratiometric background noise
    smoothing_radius - diameter of Gaussian smoothing of img_ana, in pixels
    ensure_spacing - to ensure spacing between detected peaks or not (bool 0/1)
    border_limit - how much of the border to remove peaks from in pixels
    init_smooth - if to perform an initial smoothing of the raw image or not (bool 0/1)
    """         

    if len(prev_frames) != 0:                 
        prev_frame = np.array(prev_frames)[-1]
    f_multiply = 1e3
    if binary_mask is None:
        print('Binary mask not provided.')
        binary_mask = np.ones(np.shape(img))
    if len(prev_frames) == 0 or np.shape(img) != np.shape(prev_frame):
        print('You have to provide a background image for this pipeline.')
        img_ana = np.zeros(np.shape(img))
    else:
        if init_smooth==1:
            img = ndi.filters.gaussian_filter(img.astype('float32'), 2*smoothing_radius)
            prev_frame = ndi.filters.gaussian_filter(prev_frame.astype('float32'), 2*smoothing_radius)
        
        # subtract last img (noisier, but quicker)
        img_ana = np.subtract(img,prev_frame)
        
        # divide by last image to get percentual change in img
        img_div = prev_frame
        # replace noise with a very high value to avoid detecting noise
        img_div[img_div < noise_level] = 100000
        img_ana = np.true_divide(img_ana, img_div)
        img_ana = img_ana * binary_mask
        
        img_ana = ndi.filters.gaussian_filter(img_ana, smoothing_radius)  # Gaussian filter the image, to remove noise and so on, to get a better center estimate

    # Peak_local_max all-in-one as a combo of opencv and numpy"
    thresh_abs = thresh_abs * f_multiply
    size = np.int(2 * min_dist + 1)
    img_ana = np.clip(img_ana, a_min=0, a_max=None)
    img_ana = (img_ana * f_multiply).astype('float32')

    # get filter structuring element
    footprint = cv2.getStructuringElement(cv2.MORPH_RECT, ksize=[size,size])
    # maximum filter (dilation + equal)
    image_max = cv2.dilate(img_ana, kernel=footprint)
    mask = np.equal(img_ana, image_max)
    mask &= np.greater(img_ana, thresh_abs)
    
    # get coordinates of peaks
    coordinates = np.nonzero(mask)
    intensities = img_ana[coordinates]
    # highest peak first
    idx_maxsort = np.argsort(-intensities)
    coordinates = tuple(arr for arr in coordinates)
    coordinates = np.transpose(coordinates)[idx_maxsort]

    if ensure_spacing==1:
        output = coordinates
        if len(coordinates):
            if len(coordinates) > 1000:
                coordinates = np.array([[]])
                print('Too many coordinates to ensure spacing quickly. Adjust pipeline parameters.')
            else:
                # Use KDtree to find the peaks that are too close to each other
                tree = cKDTree(coordinates, balanced_tree=False, compact_nodes=False, leafsize=50)

                indices = tree.query_ball_point(coordinates, r=min_dist, p=np.inf, return_sorted=False)
                rejected_peaks_indices = set()
                for idx, candidates in enumerate(indices):
                    if idx not in rejected_peaks_indices:
                        # keep current point and the points at exactly spacing from it
                        candidates.remove(idx)
                        dist = distance.cdist(
                            [coordinates[idx]],
                            coordinates[candidates],
                            distance.minkowski,
                            p=np.inf,
                        ).reshape(-1)
                        candidates = [
                            c for c, d in zip(candidates, dist) if d < min_dist
                        ]

                        # candidates.remove(keep)
                        rejected_peaks_indices.update(candidates)
                
                # Remove the peaks that are too close to each other
                output = np.delete(coordinates, tuple(rejected_peaks_indices), axis=0)

        coordinates = output

    # remove everything on the border
    imsize = np.shape(img)[0]
    idxremove = []
    for idx, coordpair in enumerate(coordinates):
        if coordpair[0] < border_limit or coordpair[0] > imsize - border_limit or coordpair[1] < border_limit or coordpair[1] > imsize - border_limit:
            idxremove.append(idx)
    coordinates = np.delete(coordinates,idxremove,axis=0)

    num_peaks = np.int(num_peaks)
    if len(coordinates) > num_peaks:
        coordinates = coordinates[:num_peaks]
        
    if testmode:
        return coordinates, exinfo, img_ana
    else:
        return coordinates, exinfo