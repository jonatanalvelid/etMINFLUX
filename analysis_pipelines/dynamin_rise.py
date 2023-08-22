import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import cupy as cp
from cupyx.scipy import ndimage as ndi
import cv2
import trackpy as tp
import pandas as pd

tp.quiet()

def eucl_dist(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

def dynamin_rise(img, prev_frames=None, binary_mask=None, testmode=False, exinfo=None, min_dist=4, 
                 num_peaks=100, thresh_abs_lo=500, thresh_abs_hi=8000, border_limit=10,
                 smoothing_radius=0.7, memory_frames=5, track_search_dist=4, frames_appear=10,
                 thresh_stayratio=0.7, thresh_intincratio=1.05, thresh_move_dist=3):
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
    border_limit - how much of the border to remove peaks from
    smoothing_radius - diameter of Gaussian smoothing of img_ana, in pixels
    memory_frames - number of frames for which a vesicle can disappear but still be connected to the same track
    track_search_dist - number of pixels a vesicle is allowed to move from one frame to the next
    frames_appear - number of frames ago peaks of interest appeared (to minimize noisy detections and allowing to track intensity change over time before deicision)
    thresh_stayratio - ratio of frames of the frames_appear that the peak has to be present in
    thresh_intincratio - the threshold ratio of the intensity increase in the area of the peak
    thresh_move_dist - the threshold start-end distance a peak is allowed to move during frames_appear
    """
    
    # define non-adjustable parameters
    f_multiply = 1e0
    smoothing_radius_raw = 0.6
    dog_lo = 1
    dog_hi = 3
    intensity_sum_rad = 3
    meanlen = 3
    track_len = 2 * frames_appear + 1  # number of last frames to keep in the tracking (more frames = slower track linking)
    track_search_dist = int(track_search_dist)
    thresh_stayframes = int(thresh_stayratio*frames_appear)
    memory_frames = int(memory_frames)
    
    if binary_mask is None:
        binary_mask = cp.ones(cp.shape(img)).astype('uint16')
    elif np.shape(img) != np.shape(binary_mask):
        binary_mask = cp.ones(cp.shape(img)).astype('uint16')
    else:
        binary_mask = cp.array(binary_mask).astype('uint16')
    img = cp.array(img).astype('float32')
    img_filt = ndi.filters.gaussian_filter(img, sigma=smoothing_radius_raw)
    
    # difference of gaussians to get clear peaks separated from spread-out bkg and noise
    img_dog_lo = ndi.filters.gaussian_filter(img_filt, dog_lo)
    img_dog_hi = ndi.filters.gaussian_filter(img_filt, dog_hi)
    img_dog = img_dog_lo - img_dog_hi
    
    # further filtering to get a better image for peak detection
    img_dog[img_dog < 0] = 0
    img_dog[img_dog > 30000] = 0
    img_dog = img_dog * f_multiply
    img_ana = img_dog * binary_mask
    img_ana = img_dog * cp.array(binary_mask)
    img_ana = ndi.filters.gaussian_filter(img_ana, smoothing_radius)  # Gaussian filter the image, to remove noise and so on, to get a better center estimate
    img_ana[img_ana > thresh_abs_hi] = thresh_abs_hi
    
    # Peak_local_max all-in-one as a combo of opencv and cupy
    size = int(2 * min_dist + 1)
    # get filter structuring element
    footprint = cv2.getStructuringElement(cv2.MORPH_RECT, ksize=[size,size])
    # maximum filter (dilation + equal)
    image_max = cv2.dilate(img_ana.get(), kernel=footprint)
    #return image, image_max
    mask = cp.equal(img_ana, cp.array(image_max))
    mask &= np.greater(img_ana, thresh_abs_lo)
    mask &= np.less(img_ana, thresh_abs_hi)
    
    # get coordinates of peaks
    coordinates = cp.nonzero(mask)
    intensities = img_ana[coordinates]
    # highest peak first
    idx_maxsort = cp.argsort(-intensities).get()
    coordinates = tuple(arr.get() for arr in coordinates)
    coordinates = np.transpose(coordinates)[idx_maxsort]

    # remove everything on the border
    imsize = cp.shape(img)[0]
    idxremove = []
    for idx, coordpair in enumerate(coordinates):
        if coordpair[0] < border_limit or coordpair[0] > imsize - border_limit or coordpair[1] < border_limit or coordpair[1] > imsize - border_limit:
            idxremove.append(idx)
    coordinates = np.delete(coordinates,idxremove,axis=0)

    # remove everyhting down to a certain length
    if len(coordinates) > num_peaks:
        coordinates = coordinates[:int(num_peaks),:]
    
    # add to old list of coordinates
    if exinfo is None:
        exinfo = pd.DataFrame(columns=['particle','t','x','y','intensity'])
    coordinates = coordinates[coordinates[:, 0].argsort()]

    # extract intensities summed around each coordinate
    intensities = []
    for coord in coordinates:
        intensity = cp.sum(img[coord[0]-intensity_sum_rad:coord[0]+intensity_sum_rad+1,coord[1]-intensity_sum_rad:coord[1]+intensity_sum_rad+1])/(2*intensity_sum_rad+1)**2
        intensities.append(intensity)

    # add to old list of coordinates
    if len(exinfo) > 0:
        timepoint = max(exinfo['t'])+1
    else:
        timepoint = 0
    if len(coordinates)>0:
        coords_df = pd.DataFrame(np.hstack((np.array(range(len(coordinates))).reshape(-1,1),timepoint*np.ones(len(coordinates)).reshape(-1,1),coordinates,np.reshape(cp.array(intensities).get(),(-1,1)))),columns=['particle','t','x','y','intensity'])
        tracks_all = exinfo.append(coords_df)
    else:
        tracks_all = exinfo
    
    # event detection
    coords_event = np.empty((0,3))
    if len(tracks_all) > 0:
        # link coordinate traces (only last track_len frames)
        tracks_all = tracks_all[tracks_all['t']>max(tracks_all['t'])-track_len]
        tracks_all = tp.link(tracks_all, search_range=track_search_dist, memory=memory_frames, t_column='t')
        
        # event detection of appearing vesicles
        # conditions:
        # 1. one track appears frames_appear ago
        # 2. track stays for at least thresh_stayframes (7?) frames
        # 3. intensity of track spot increases over thresh_stayframes frames with at least thresh_intincratio (3x?)
        # (4. check that track has not moved too much in the last frames?)
        
        if timepoint >= 2*frames_appear:
            tracks_timepoint = tracks_all[tracks_all['t']==timepoint-frames_appear]
            tracks_before = tracks_all[tracks_all['t']<timepoint-frames_appear]
            tracks_after = tracks_all[tracks_all['t']>timepoint-frames_appear]
            #particle_ids_after = np.unique(tracks_after['particle'])
            particle_ids_before = np.unique(tracks_before['particle'])
            for _,track in tracks_timepoint.iterrows():
                # check for appearing tracks
                particle_id = int(track['particle'])
                if particle_id not in particle_ids_before:
                    # check that it stays for at least thresh_stayframes frames
                    track_self_after = tracks_after[tracks_after['particle']==particle_id]
                    if len(track_self_after) > thresh_stayframes:
                        # check that intensity of spot increases over the thresh_stay frames with at least thresh_intincratio
                        track_self = track_self_after.tail(1)
                        prev_frames = np.array(prev_frames)
                        track_intensity_before = np.sum(np.sum(prev_frames[:, int(track_self['x'])-intensity_sum_rad:int(track_self['x'])+intensity_sum_rad+1,
                                                                         int(track_self['y'])-intensity_sum_rad:int(track_self['y'])+intensity_sum_rad+1],
                                                               axis=1),axis=1)/(2*intensity_sum_rad+1)**2
                        track_intensity_after = track_self_after['intensity']
                        int_before = np.mean(track_intensity_before[0:meanlen])
                        int_init = np.mean(track_intensity_after.iloc[0:meanlen])
                        int_last = np.mean(track_intensity_after.iloc[-(meanlen+1):-1])
                        intincratio_before = int_init/int_before
                        intincratio = int_last/int_init
                        intincrratio_tot = int_last/int_before
                        if intincratio_before > thresh_intincratio or intincratio > thresh_intincratio or intincrratio_tot > thresh_intincratio:
                            # check that track has not moved too much since it appeared
                            track_self_after_start = track_self_after.head(1)
                            track_self_after_end = track_self_after.tail(1)
                            d_start_end = eucl_dist((int(track_self_after_start['x']),int(track_self_after_start['y'])),(int(track_self_after_end['x']),int(track_self_after_end['y'])))
                            if d_start_end < thresh_move_dist:
                                # if all conditions are true: potential appearence event frames_appear ago, save coord of curr position
                                print([intincratio_before, intincratio, intincrratio_tot])
                                coords_event = np.array([[int(track_self['y']), int(track_self['x'])]])
                                break
    if testmode:
        return coords_event, tracks_all, img_ana.get()
    else:
        return coords_event, tracks_all