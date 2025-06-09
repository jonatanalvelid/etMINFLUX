import numpy as np
from scipy import ndimage as ndi
import cv2
import trackpy as tp
import pandas as pd

tp.quiet()

def peak_detection_stationary_dualcolor(img, img_ch2, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None,
                       maxfilter_kersize=5, thresh_abs=10, smoothing_radius=1, border_limit=15, init_smooth=1,
                       num_prev=4, msm_thresh=0.7, ch2sig_thresh_lo=0.2, ch2sig_thresh_hi=10):
    
    """
    Analysis pipeline to detect bright peaks that are stationary in an image, using a maximum intensity detection filter,
    and then checking if they are stationary inside the previous frames. Only return a random stationary peak. This version
    also considers a second channel, from img_ch2, and checks that the mean signal in that channel is within a certain range. 
    
    Common parameters:
    img - current image
    img_ch2 - current image in the second channel
    prev_frames - previous image(s)
    binary_mask - binary mask of the region to consider
    exinfo - pandas dataframe of the detected vesicles and their track ids from the previous frames
    presetROIsize - preset ROI sizes boolean, if used or not (not applicable for this pipeline)

    Pipeline specific parameters:
    maxfilter_kersize - size of kernel for maximum filtering
    thresh_abs - low intensity threshold in img_ana of the peaks to consider
    smoothing_radius - diameter of Gaussian smoothing of img_ana, in pixels
    border_limit - how much of the border to remove peaks from in pixels
    init_smooth - if to perform an initial smoothing of the raw image or not (bool 0/1)
    num_prev - number of frames (num_prev+1) to consider for checking if it is stationary
    msm_thresh - maximum threshold for msm, mean squared movement, between individual frames
    ch2sig_thresh_lo - low threshold for the mean channel 2 signal around a stationary peak
    ch2sig_thresh_hi - high threshold for the mean channel 2 signal around a stationary peak
    """
    roi_sizes = False

    if binary_mask is None or np.shape(binary_mask) != np.shape(img):
        binary_mask = np.ones(np.shape(img)).astype('uint16')

    img = np.array(img).astype('float32')
    img_ch2 = np.array(img_ch2).astype('float32')
    if init_smooth==1:
        img = ndi.gaussian_filter(img, smoothing_radius)
    
    # multiply with binary mask
    img_ana = img * np.array(binary_mask)

    # Peak_local_max in current image as a combo of opencv and numpy
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
    
    coordinates = np.flip(coordinates, axis=1)

    # create initial exinfo, if first analysis frame
    if exinfo is None:
        exinfo = pd.DataFrame(columns=['particle','frame','x','y','intensity','intensitych2'])
    coordinates = coordinates[coordinates[:, 0].argsort()]
    
    # extract intensities summed around each coordinate
    intensity_sum_rad = 3
    intensity_sum_rad_ch2 = 3
    intensities = []
    intensitiesch2 = []
    for coord in coordinates:
        intensity = np.sum(img[coord[1]-intensity_sum_rad:coord[1]+intensity_sum_rad+1,coord[0]-intensity_sum_rad:coord[0]+intensity_sum_rad+1])/(2*intensity_sum_rad+1)**2
        intensitych2 = np.sum(img_ch2[coord[1]-intensity_sum_rad_ch2:coord[1]+intensity_sum_rad_ch2+1,coord[0]-intensity_sum_rad_ch2:coord[0]+intensity_sum_rad_ch2+1])/(2*intensity_sum_rad_ch2+1)**2
        intensities.append(intensity)
        intensitiesch2.append(intensitych2)
    
    # add to old list of coordinates
    if len(exinfo) > 0:
        timepoint = int(max(exinfo['frame'])+1)
    else:
        timepoint = 0
    if len(coordinates)>0:
        coords_df = pd.DataFrame(np.hstack((np.array(range(len(coordinates))).reshape(-1,1),timepoint*np.ones(len(coordinates)).reshape(-1,1),coordinates,np.array(intensities).reshape(-1,1),np.array(intensitiesch2).reshape(-1,1))),columns=['particle','frame','x','y','intensity','intensitych2'])
        tracks_all = pd.concat([exinfo, coords_df])
    else:
        tracks_all = exinfo
    
    memory_frames = 0
    track_search_dist = 7
    coords_events = []
    if timepoint >= num_prev:
        tracks_all = tp.link(tracks_all, search_range=track_search_dist, memory=memory_frames, t_column='frame');
        tracks_use = tracks_all[tracks_all['frame']>timepoint-(num_prev+1)].copy()
        for particle in np.unique(tracks_use['particle']):
            particle_track = tracks_use[tracks_use['particle']==particle]
            if len(particle_track) > num_prev:
                movement = tp.motion.compute_drift(particle_track)
                # get the frame-to-frame square movement array in the last num_prev+1 frames
                sm = [np.sqrt((line1[1]['y']-line2[1]['y'])**2+(line1[1]['x']-line2[1]['x'])**2) for line1, line2 in zip(movement.iterrows(),movement.iloc[1:].iterrows())]
                if np.mean(sm) < msm_thresh:
                    ch2sig = np.mean(particle_track['intensitych2'])
                    if ch2sig > ch2sig_thresh_lo and ch2sig < ch2sig_thresh_hi:
                        coords_events.append([[int(particle_track.iloc[-1]['x']), int(particle_track.iloc[-1]['y'])]])
    
    coords_events = np.array(coords_events)
    
    # generate random coordinate number to use from the user-provided inputs
    if len(coords_events) > 0:
        coord_num = np.random.randint(0, len(coords_events))
        coords_event_return = coords_events[coord_num]
    else:
        coords_event_return = np.empty((0,3))
    
    return coords_event_return, roi_sizes, tracks_all, img_ch2