# -*- coding: utf-8 -*-

__version__ = '0.92'
# $Source$

import numpy as np

from ..mrrcore import empty
from ..read import read_dicom_set
from ..unwrapping import DerivativeVariance



def read_times(files, mask=None, unwrap=True, verbose=False):
    '''
    Read-in a bunch of image-sets belonging to one timeline.
    Returns a list of 3d-MRRArrays.
    files contains a list of filenames, one for each timepoint.
    
    files : list
        list conaining one filename per timepoint. to be read-in.
    mask : array
        2d-array containing the mask. Needed if unwrapping is used.
    unwrap : [True | False]
        wether to unwrap the images
    '''
    times = []
    for f in files:
        times.append(read_dicom_set(f, unwrap=unwrap, mask=mask, 
                                    verbose=verbose)
                    )
    return times


def create_timeline(times, normalize=True, verbose=False):
    '''
    Create a timeline (3d-MRRArray with z-axis as time) of a list of 
    image-sets as created with mrr.read_times().
    
    If normalize is set to True (default), every image-set is attempted to be 
    normalized, i.e. the images are corrected to represent the same multiple
    of 2*pi.
    '''
    #Normalize
    if normalize:
        for t, img_set in enumerate(times):
            if verbose:
                print '\nNormalizing imageset %i...' % t
            normalize_image_set(img_set, threshold=0.7, verbose=verbose)
    #Create Timeline
    shape = times[0].shape
    unwrapped = times[0].unwrapped
    Timeline = empty((len(times),shape[1],shape[2]), 
                         orig_file='Zeitreihe', 
                         unwrapped=unwrapped
                         )
    for t, timepoint in enumerate(times):
        Timeline[t] = timepoint.mean(axis=0)
    return Timeline


def read_timeline(files, mask=None, unwrap=True,
                  normalize=True, norm_thresh=0.7,
                  return_times=False,
                  verbose=False):
    '''
    Read-in a set of files and create a 3D-MRRArray representing the timeline
    (z-axis corresponing to the time-axis).
    
    mask : array
        2d-array containing a mask, mandatory if unwrapping is used
    unwrap : [True | False]
        wether to unwrap the images individually
    normalize : [True | False]
        wether to normalize the individual image-sets using 
        mrr.normalize_image_set
    norm_thresh : float, optional
        threshold to be used in nomralization
    return_times : [True | False]
        If True, the list of images is returned as well.
    '''
    #Read-In
    times = read_times(files, mask, unwrap, verbose)
    #Normalize
    if normalize:
        for t, img_set in enumerate(times):
            if verbose:
                print '\nNormalizing imageset %i...' % t
            normalize_image_set(img_set, threshold=norm_thresh, 
                                verbose=verbose)
    #Create Timeline
    Timeline = create_timeline(times)
    if return_times:
        return Timeline, times
    else:
        return Timeline


def calculate_jump(img1, img2):
    '''
    Calculates the mean difference of all pixel not masked between two images.
    '''
    # only use pixels not masked out in any image
    # (pixels are set to True if they are valid pixel)
    mask = img1.mask & img2.mask
    return (img2 - img1)[mask == True].mean().phase


def normalize_image_set(img_set, threshold=0.7, continuous=False,
                        verbose=True, test=False):
    '''
    Normalizes an image_set by recursively adding or substracting 1.0 to all
    pixels of individual images.
    '''
    changed = False
    for i, img in enumerate(img_set[1:]):
        prior = img_set[i] if continuous else img_set[0]
        if verbose:
            print 'Normalizing Image %i...' % (i+1)
        jump = calculate_jump(prior, img)
        if test and verbose:
            print 'Durchschnittliche Abweichung: %f' % jump
        if jump > threshold:
            img_set[i+1] -= 1.0
            changed = True
            if verbose:
                print '\tSubtracted 1.0'
        elif jump <= -threshold:
            img_set[i+1] += 1.0
            changed = True
            if verbose:
                print '\tAdded 1.0'
    if changed:
        normalize_image_set(img_set, threshold, verbose, test)


def broaden_mask(img, threshold=0.05, qual=None):
    """
    Mask pixels based on a quality map.
    
    All pixels, whose quality is below the maximum quality of the image times 
    the given threshold are masked. Per default, this means the lowest 5% of 
    quality are masked out.
    If no qualitymap is provided, the Phase Derivative Variance is calculated.
    
    Parameters
    ----------
    img : MRRarray
        the data to be masked
    threshold : float
        Determines how many pixels to mask out
    qual : array (optional)
        Optional qualitymap to use.
    """
    if not np.any(qual):
        qual = DerivativeVariance(img.phase)
        qual = qual[img.mask==True].max()*1.1 - qual
    max_value = qual[img.mask==True].max()
    img['mask'][qual<max_value*threshold] = False
               
               
def unwrap_timeline(timeline, threshold=0.6, verbose=False):
    '''
    Unwraps a timeline as created by create_timeline or read_timeline along
    the time-axis (z-direction).
    The array is changed in place, i.e. the original array will be lost.
    '''
    for t,img in enumerate(timeline[1:]):
        jump = calculate_jump(timeline[t], img)
        if jump < -threshold:
            timeline[t+1:] -= 1.0
            if verbose:
                print 'Subtracted 1.0 at time %i' % t+1
        elif jump > threshold:
            timeline[t+1:] += 1.0
            if verbose:
                print 'Added 1.0 at time %i' % t+1   