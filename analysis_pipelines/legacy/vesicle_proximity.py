import numpy as np
import cupy as cp
from cupyx.scipy import ndimage as ndi
import cv2
import trackpy as tp
import pandas as pd

tp.quiet()

def eucl_dist(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

def vesicle_proximity(img, prev_frames=None, binary_mask=None, testmode=False, exinfo=None, min_dist=4, 
                      num_peaks=100, thresh_abs=500, border_limit=10, smoothing_radius=0.7, 
                      ves_dist=7, stat_frames=5, track_search_dist=4, track_mov_thresh=2.5):
    """
    Common parameters:
    img - current image
    prev_frames - previous image(s)
    binary_mask - binary mask of the region to consider
    testmode - to return preprocessed image or not
    exinfo - pandas dataframe of the detected vesicles and their track id from the previous frames

    Pipeline specfic parameters:
    min_dist - minimum distance in pixels between two peaks
    num_peaks - number of peaks to detect
    thresh_abs - intensity threshold in img_ana of the peaks to consider
    border_limit - how much of the border to remove peaks from 
    smoothing_radius - diameter of Gaussian smoothing of img_ana, in pixels
    ves_dist - max distance between two vesicles at the event detection
    stat_frames - number of frames the connecting vesicles have to be stationary for
    track_search_dist - number of pixels a vesicle is allowed to move from one frame to the next
    track_mov_thresh - distance one of the vesicles in a potential event detection has to have travelled inside the last 3*stat_frame frames

    Derived parameters:
    track_len - number of previous frames to keep in the tracking (more frames = slower track linking), has to be large enough to take an informed decision on event
    track_len_thresh - number of frames the vesicles of a potential event has to have been visible, to avoid noisy detections
    memory_frames - number of frames for which a vesicle can disappear but still be connected to the same track
    """
    
    # define non-adjustable parameters
    prev_tracks = exinfo
    f_multiply = 1e1
    stat_frames = int(stat_frames)
    track_search_dist = int(track_search_dist)
    track_len = 4*stat_frames+1
    track_len_thresh = 3/2*stat_frames
    memory_frames = stat_frames
    dog_lo = 1
    dog_hi = 3
    
    if (binary_mask is None) or (np.shape(img) != np.shape(binary_mask)):
        binary_mask = cp.ones(cp.shape(img)).astype('int16')
    else:
        binary_mask = cp.array(binary_mask).astype('int16')
    img = cp.array(img).astype('int16')

    # gaussian filter raw image
    img = ndi.filters.gaussian_filter(img, sigma=smoothing_radius)

    # difference of gaussians to get clear peaks separated from spread-out bkg and noise
    img_dog_lo = ndi.filters.gaussian_filter(img, dog_lo)
    img_dog_hi = ndi.filters.gaussian_filter(img, dog_hi)
    img_dog = img_dog_lo - img_dog_hi

    # further filtering to get a better image for peak detection
    img_dog[img_dog<0] = 0
    img_dog = img_dog*f_multiply
    img_ana = img_dog * cp.array(binary_mask)
    img_ana = ndi.filters.gaussian_filter(img_ana, smoothing_radius)  # Gaussian filter img_ana, to remove noise and so on, to get a better center estimate
    
    # Peak_local_max all-in-one as a combo of opencv and cupy
    size = int(2 * min_dist + 1)
    # get filter structuring element
    footprint = cv2.getStructuringElement(cv2.MORPH_RECT, ksize=[size,size])
    # maximum filter (dilation + equal)
    image_max = cv2.dilate(img_ana.get(), kernel=footprint)
    #return image, image_max
    mask = cp.equal(img_ana,cp.array(image_max))
    mask &= cp.greater(img_ana, thresh_abs)
    
    # get coordinates of peaks
    coordinates = cp.nonzero(mask)
    intensities = img_ana[coordinates]
    # highest peak first
    idx_maxsort = cp.argsort(-intensities).get()
    coordinates = tuple(arr.get() for arr in coordinates)
    coordinates = np.transpose(coordinates)[idx_maxsort]
    coordinates = cp.fliplr(coordinates)

    # remove everything on the border
    imsize = cp.shape(img)[0]
    idxremove = []
    for idx, coordpair in enumerate(coordinates):
        if coordpair[0] < border_limit or coordpair[0] > imsize - border_limit or coordpair[1] < border_limit or coordpair[1] > imsize - border_limit:
            idxremove.append(idx)
    coordinates = np.delete(coordinates,idxremove,axis=0)

    # remove peaks down to a certain number
    if len(coordinates) > num_peaks:
        coordinates = coordinates[:int(num_peaks),:]
    
    # add to old list of coordinates
    if prev_tracks is None:
        prev_tracks = pd.DataFrame(columns=['particle','t','x','y'])
    coordinates = np.flip(coordinates, axis=1)
    coordinates = coordinates[coordinates[:, 0].argsort()]
    if len(prev_tracks) > 0:
        timepoint = max(prev_tracks['t'])+1
    else:
        timepoint = 0
    if len(coordinates)>0:
        coords_df = pd.DataFrame(np.hstack((np.array(range(len(coordinates))).reshape(-1,1),timepoint*np.ones(len(coordinates)).reshape(-1,1),coordinates)),columns=['particle','t','x','y'])
        tracks_all = prev_tracks.append(coords_df)
    else:
        tracks_all = prev_tracks
    
    # link coordinate traces (only last track_len frames)
    coords_events = np.empty((0,2))
    if len(tracks_all) > 0:
        tracks_all = tracks_all[tracks_all['t']>max(tracks_all['t'])-track_len]
        tracks_all = tp.link(tracks_all, search_range=track_search_dist, memory=memory_frames, t_column='t')
        # event detection of fusing vesicles (two detected tracks becoming one)
        # conditions:
        # 1. one track (#1) disappears
        # 2. another track (#2) close by at the time of disappearence
        # 3. both tracks (#1 and #2) have tracked points in 50% of time points leading up to disappearance
        # 4. at least one track has moved an accumulated vectorial distance above a threshold
        # 5. at least one track has moved an accumulated absolute distance above a threshold (TODO: is this not always true if the above is true?)
        d_self = 0
        if timepoint >= stat_frames:
            tracks_timepoint = tracks_all[tracks_all['t']==timepoint-stat_frames]
            tracks_timepoint_around = tracks_all.loc[(tracks_all['t']>timepoint-stat_frames-3) & (tracks_all['t']<=timepoint-stat_frames)]
            tracks_after = tracks_all[tracks_all['t']>timepoint-stat_frames]
            particle_ids_after = np.unique(tracks_after['particle'])
            for _,track_old1 in tracks_timepoint.iterrows():
                event_found = False
                # check for disappearing tracks (1), that have stayed disappeared for more than x number of frames (i.e. check which tracks
                # disappeared at time t=timepoint-x+1 compare to timpoint-x, and check that they have not appeared again after that)
                particle_id_old1 = int(track_old1['particle'])
                if particle_id_old1 not in particle_ids_after:
                    # if disappearing track (1):
                    # (2) check if there were two tracks close to each other at the moment of disappearing, inside ves_dist
                    coord_old1 = (track_old1['x'],track_old1['y'])
                    for _,track_old2 in tracks_timepoint_around.iterrows():
                        particle_id_old2 = int(track_old2['particle'])
                        if particle_id_old1 != particle_id_old2:
                            coord_old2 = (track_old2['x'],track_old2['y'])
                            d = eucl_dist(coord_old1,coord_old2)
                            if d < ves_dist:
                                # if close-by tracks (2):
                                # (3) check if both vesicles have temporally longer tracks than track_len_thresh points in the last 3*stat_frames frames.
                                tracks_timepoint_before = tracks_all.loc[(tracks_all['t']>max(0,timepoint-4*stat_frames)) & (tracks_all['t']<timepoint-stat_frames)]
                                tracks_self1 = tracks_timepoint_before[tracks_timepoint_before['particle']==particle_id_old1]
                                tracks_self1.reset_index(drop=True, inplace=True)
                                tracks_self2 = tracks_timepoint_before[tracks_timepoint_before['particle']==particle_id_old2]
                                tracks_self2.reset_index(drop=True, inplace=True)
                                if (len(tracks_self1) > track_len_thresh) and (len(tracks_self2) > track_len_thresh):
                                    # if temporally long enough tracks (3):
                                    # (4) check that at least one of the vesciles has moved an accumulated distance
                                    # (vectorial) d from its starting position in the last 3*stat_frames before disappearance.
                                    track_self1_start = tracks_self1.head(1)
                                    track_self1_end = tracks_self1.tail(1)
                                    track_self2_start = tracks_self2.head(1)
                                    track_self2_end = tracks_self2.tail(1)
                                    d_start_end1 = eucl_dist((int(track_self1_start['x']),int(track_self1_start['y'])),(int(track_self1_end['x']),int(track_self1_end['y'])))
                                    d_start_end2 = eucl_dist((int(track_self2_start['x']),int(track_self2_start['y'])),(int(track_self2_end['x']),int(track_self2_end['y'])))
                                    if d_start_end1 > track_mov_thresh or d_start_end2 > track_mov_thresh:
                                        # if long enough accumulated vectorial distance for one of the vesicles (4):
                                        # (5) check that atleast one of the vesicles has moved longer than an accumulated
                                        # length (absolute values) of d pixels in the last 3*stat_frames frames before disappearance.
                                        d_self = 0
                                        for idx3,track_self in tracks_self1.iterrows():
                                            if idx3==0:
                                                coord_self_prev = (track_self['x'],track_self['y'])
                                            else:
                                                coord_self_curr = (track_self['x'],track_self['y'])
                                                d_self += eucl_dist(coord_self_prev, coord_self_curr)
                                                if d_self > track_mov_thresh:
                                                    # if all condiitions are true (1-5):
                                                    # potential fusion event stat_frames frames ago, save coordinates of current vesicle position
                                                    track_current = tracks_after[tracks_after['particle']==particle_id_old2]
                                                    if len(track_current) > 0:
                                                        track_current = track_current.tail(1)
                                                        coord_event_current = np.array([[int(track_current['y']), int(track_current['x'])]])
                                                        coords_events = np.append(coords_events, coord_event_current, axis=0)
                                                        event_found = True
                                                        break
                                                else:
                                                    coord_self_prev = coord_self_curr
                                        if not event_found:
                                            d_self = 0
                                            for idx3,track_self in tracks_self2.iterrows():
                                                if idx3==0:
                                                    coord_self_prev = (track_self['x'],track_self['y'])
                                                else:
                                                    coord_self_curr = (track_self['x'],track_self['y'])
                                                    d_self += eucl_dist(coord_self_prev, coord_self_curr)
                                                    if d_self > track_mov_thresh:
                                                        # if all condiitions are true (1-5):
                                                        # potential fusion event stat_frames frames ago, save coordinates of current vesicle position
                                                        track_current = tracks_after[tracks_after['particle']==particle_id_old2]
                                                        if len(track_current) > 0:
                                                            track_current = track_current.tail(1)
                                                            coord_event_current = np.array([[int(track_current['y']), int(track_current['x'])]])
                                                            coords_events = np.append(coords_events, coord_event_current, axis=0)
                                                            break
                                                    else:
                                                        coord_self_prev = coord_self_curr
                                break

    if testmode:
        img_ana = img_ana.get()
        return coords_events, tracks_all, img_ana
    else:
        return coords_events, tracks_all