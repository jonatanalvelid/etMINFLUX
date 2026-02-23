import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
from scipy import ndimage as ndi
import cv2
import trackpy as tp
import pandas as pd
from scipy.optimize import curve_fit

tp.quiet()

def eucl_dist(a,b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

def f(x, A, B):
    return A*x + B

def gag_signalrise(img_ch1, prev_frames=None, binary_mask=None, exinfo=None, presetROIsize=None,
                     min_dist_appear=5, num_peaks=300, thresh_abs_lo=1.7, thresh_abs_hi=10, finalintlo=0.75, 
                     finalinthi=5, border_limit=10, memory_frames=6, track_search_dist=10, frames_appear=4, 
                     thresh_intincratio=1.3, thresh_intincratio_max=15, intincslope=0.1, thresh_move_dist=1.3):
    
    """ 
    Analysis pipeline that detects a slowly increasing spot-like signal (0-maximum intensity inside ~10-20 frames),
    used for detecting slow accumulation of gag signal, with a low confocal frame rate, at potential virus budding sites.
    
    Common parameters:
    img_ch1 - current image
    prev_frames - previous image(s)
    binary_mask - binary mask of the region to consider
    exinfo - pandas dataframe of the detected vesicles and their track ids from the previous frames

    Pipeline specific parameters:
    min_dist_appear - minimum distance in pixels between two peaks
    num_peaks - number of peaks to track
    thresh_abs_lo - low intensity threshold in img_ana of the peaks to consider
    thresh_abs_hi - high intensity threshold in img_ana of the peaks to consider
    finalintlo = low threshold on the final intensity of the peak
    finalinthi = high threshold on the final intensity of the peak
    border_limit - how much of the border to remove peaks from
    memory_frames - number of frames for which a vesicle can disappear but still be connected to the same track
    track_search_dist - number of pixels a vesicle is allowed to move from one frame to the next
    frames_appear - number of frames ago peaks of interest appeared (to minimize noisy detections and allowing to track intensity change over time before deicision)
    thresh_stayratio - ratio of frames of the frames_appear that the peak has to be present in
    thresh_intincratio - the threshold ratio of the intensity increase in the area of the peak
    thresh_intincratio_max - the maximum threshold ratio of the intensity increase in the area of the peak
    intincslope = threshold on the slope of the intensity increase over frames_appear frames, used insted of threshold_intincratio when the data allows (long enough)
    thresh_move_dist - the threshold start-end distance a peak is allowed to move during frames_appear
    """
    
    roi_sizes = False

    # define non-adjustable parameters
    smoothing_radius_raw = 1.5  # pixels
    intensity_sum_rad = 2  # pixels
    frames_appear = int(frames_appear)
    memory_frames = int(memory_frames)
    meanlen = int(np.min([2, frames_appear]))  # for such short frames_appear (3), just use frames_appear instead
    track_search_dist = int(track_search_dist)
    thresh_stayframes = int(frames_appear*0.7)  # can be gone after it appears, for 30% of frames that comes
    
    if binary_mask is None:
        img_ch1 = np.array(img_ch1).astype('float32')
    else:
        img_ch1 = img_ch1 * binary_mask
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
    
    # create initial exinfo, if first analysis frame
    if exinfo is None:
        exinfo = pd.DataFrame(columns=['particle','t','x','y','intensity'])
    coordinates = coordinates[coordinates[:, 0].argsort()]
    
    # extract intensities summed around each coordinate
    intensities = []
    for coord in coordinates:
        intensity = np.sum(img_ch1[coord[0]-intensity_sum_rad:coord[0]+intensity_sum_rad+1,coord[1]-intensity_sum_rad:coord[1]+intensity_sum_rad+1])/(2*intensity_sum_rad+1)**2
        intensities.append(intensity)
    
    # add to old list of coordinates
    if len(exinfo) > 0:
        timepoint = max(exinfo['t'])+1
    else:
        timepoint = 0
    if len(coordinates)>0:
        coords_df = pd.DataFrame(np.hstack((np.array(range(len(coordinates))).reshape(-1,1),timepoint*np.ones(len(coordinates)).reshape(-1,1),coordinates,np.array(intensities).reshape(-1,1))),columns=['particle','t','x','y','intensity'])
        tracks_all = pd.concat([exinfo, coords_df])
    else:
        tracks_all = exinfo
    
    # event detection
    imgsize = np.shape(img_ch1)[0]
    coords_event = np.empty((0,3))
    if len(tracks_all) > 0:
        # link coordinate traces (only last memory_frames+frames_appear frames, in order to be able to link tracks memory_frames ago for when a potential event appeared)
        tracks_all = tracks_all[tracks_all['t']>max(tracks_all['t'])-3*frames_appear]
        tracks_all = tp.link(tracks_all, search_range=track_search_dist, memory=memory_frames, t_column='t')
        
        # event detection of appearing vesicles
        # conditions:
        # 1. one track appears inside frames_appear to 2*frames_appear ago
        # 2. track stays for at least thresh_stayframes frames from time_timepoint
        # 3. there is no other peak next to it when appearing (too uncertain)
        # 4. intensity of track spot increases over frames_appear frames, before and after track appeared, with at least thresh_intincratio
        # 5. check that final intensity is above a certain threshold
        # 6. check that track has not moved too much in the last frames
        # 7. check that final position is not outside the border_limit
        
        if timepoint >= 2*frames_appear:
            tracks_timepoint = tracks_all[tracks_all['t']==timepoint-frames_appear]
            tracks_after = tracks_all[tracks_all['t']>timepoint-frames_appear]
            if timepoint >= 3*frames_appear:
                tracks_prebefore = tracks_all[tracks_all['t']<timepoint-2*frames_appear]
            tracks_before = tracks_all[(tracks_all['t']<timepoint-frames_appear) & (tracks_all['t']>=timepoint-2*frames_appear)]
            if timepoint >= 3*frames_appear:
                particle_ids_before = np.unique(tracks_prebefore['particle'])  # take older frames as check
            else:
                particle_ids_before = np.unique(tracks_before['particle'])  # take last frames_appear as check
            for _, track in tracks_timepoint.iterrows():
                # check for appearing tracks, at 1x frames_appear to 3x frames_appear ago
                particle_id = int(track['particle'])
                if particle_id not in particle_ids_before:
                    # check that it stays for at least thresh_stayframes frames
                    track_self_after = tracks_after[tracks_after['particle']==particle_id]
                    if len(track_self_after) >= thresh_stayframes:
                        track_full = tracks_all[tracks_all['particle']==particle_id]
                        appearance_frame = np.min(track_full['t'])
                        track_appearance = track_full[track_full['t']==appearance_frame]
                        other_tracks_appearance = tracks_all[(tracks_all['t']==appearance_frame) & (tracks_all['particle']!=particle_id)]
                        particle_dists = [eucl_dist((int(track_appearance['x']),int(track_appearance['y'])),(int(x2),int(y2))) for x2,y2 in zip(other_tracks_appearance['x'],other_tracks_appearance['y'])]
                        if np.min(particle_dists) > min_dist_appear:
                            track_self_before = tracks_before[tracks_before['particle']==particle_id]
                            # check that intensity of spot increases over the thresh_stay frames with at least thresh_intincratio
                            track_self = track_self_after.tail(1)
                            prev_frames = np.array(prev_frames).astype('float32')
                            if len(track_self_before) > 0:
                                track_intensity_before = np.array(track_self_before['intensity'])
                            else:
                                track_intensity_before = np.sum(prev_frames[-2*frames_appear:-frames_appear, int(track_self['x'])-intensity_sum_rad:int(track_self['x'])+intensity_sum_rad+1,
                                                                        int(track_self['y'])-intensity_sum_rad:int(track_self['y'])+intensity_sum_rad+1],
                                                                    axis=(1,2))/(2*intensity_sum_rad+1)**2
                            track_intensity_after = np.array(track_self_after['intensity'])
                            if len(track_self_before) == 0:
                                track_intensity_arounddetect = np.array([track_intensity_before[-1], tracks_timepoint[tracks_timepoint['particle']==particle_id]['intensity'].iloc[0], track_intensity_after[0]])
                                int_detect = np.mean(track_intensity_arounddetect)
                                int_before = np.mean(track_intensity_before[:meanlen])
                                int_after = np.mean(track_intensity_after[-meanlen:])
                                intincrratio_before = int_detect/int_before
                                intincrratio_after = int_after/int_detect
                                if (intincrratio_before > thresh_intincratio and intincrratio_before < thresh_intincratio_max) and (intincrratio_after > thresh_intincratio and intincrratio_after < thresh_intincratio_max):
                                    # check that final intensity of track is at least above finalint
                                    if int_after > finalintlo and int_after < finalinthi:
                                        # check that track has not moved too much since it appeared
                                        d_vects = [eucl_dist((int(x1),int(y1)),(int(x2),int(y2))) for x1,y1,x2,y2 in zip(track_self_after['x'].tail(-1),track_self_after['y'].tail(-1),track_self_after['x'],track_self_after['y'])]
                                        if np.mean(d_vects) < thresh_move_dist:
                                            # if all conditions are true: potential appearence event frames_appear ago, save coord of curr position
                                            if int(track_self['x']) > border_limit and int(track_self['x']) < imgsize - border_limit and int(track_self['y']) > border_limit and int(track_self['y']) < imgsize - border_limit:
                                                # last check that event is not inside the border, if it is just continue looking at the next track
                                                coords_event = np.array([[int(track_self['x']), int(track_self['y'])]])
                                                break
                            else:
                                track_intensity_all = np.concatenate((track_intensity_before, np.array([tracks_timepoint[tracks_timepoint['particle']==particle_id]['intensity'].iloc[0]]), track_intensity_after)).tolist()
                                x = np.linspace(0,len(track_intensity_all)-1,len(track_intensity_all))
                                sigma = np.ones(len(x))
                                sigma[[0]] = 0.01
                                popt, _ = curve_fit(f, x, track_intensity_all, sigma=sigma)
                                if popt[0] > intincslope:
                                    # check that final intensity of track is at least above finalint
                                    int_after = np.mean(track_intensity_after[-meanlen:])
                                    if int_after > finalintlo and int_after < finalinthi:
                                        # check that track has not moved too much since it appeared
                                        d_vects = [eucl_dist((int(x1),int(y1)),(int(x2),int(y2))) for x1,y1,x2,y2 in zip(track_self_after['x'].tail(-1),track_self_after['y'].tail(-1),track_self_after['x'],track_self_after['y'])]
                                        if np.mean(d_vects) < thresh_move_dist:
                                            # if all conditions are true: potential appearence event frames_appear ago, save coord of curr position
                                            if int(track_self['x']) > border_limit and int(track_self['x']) < imgsize - border_limit and int(track_self['y']) > border_limit and int(track_self['y']) < imgsize - border_limit:
                                                # last check that event is not inside the border, if it is just continue looking at the next track
                                                coords_event = np.array([[int(track_self['x']), int(track_self['y'])]])
                                                break

    coords_event = np.flip(coords_event, axis=1)

    return coords_event, roi_sizes, tracks_all, img_ana
