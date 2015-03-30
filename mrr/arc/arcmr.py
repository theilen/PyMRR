# -*- coding: utf-8 -*-
"""
Custom functions to evalue phase images with ARF-contrast.

Created on Fri Aug 22 13:02:06 2014

@author: Sebastian Theilenberg
"""

__version__ = "1.2.1"
# $Source$

import numpy as np
import math
import matplotlib.pyplot as plt
import copy

from ..plotting import display
from ..mrrcore import empty_like
from .. import bvalue


def roi_circle(center, radius, mask=None):
    '''
    Creates the coordinates of a circular region-of-interest in an 2D-array
    around the provided center-point. If a mask is present, masked values are
    ignored.
    Returns two lists containing the sorted x- and y-values of that ROI.

    Parameters
    ----------
    center : (x,y)
        tuple containing x- and y-value of the center-point of the desired ROI.
    radius : int
        Integer specifying the desired radius of the circle.
    mask : 2darray (optional)

    Returns
    -------
    coords : ([x],[y])
        2-tuple containing one list of x-coordinates and one of y-coordinates.
    '''
    x,y = center
    roi_x = []
    roi_y = []
    for _y in range(y-radius,y+radius+1):
        for _x in range(x-radius,x+radius+1):
            norm = np.sqrt(np.power(_x-x,2) + np.power(_y-y,2))
            if np.round(norm) <= radius:
                #ignore masked pixels
                if np.any(mask):
                    if mask[_y,_x] == False:
                        continue
                roi_x.append(_x)
                roi_y.append(_y)

    return roi_x, roi_y



def create_overlay(shape, x, y):
    '''
    Returns an array of given shape containing NaN in every pixel except these
    specified by x and y.
    '''
    if not len(shape) == 2:
        raise ValueError('shape has to be 2-dimensional!')
    result = np.empty(shape)
    result[:,:] = np.nan
    result[y,x] = 1
    return result

#transparent colormap
over_red = copy.copy(plt.cm.get_cmap(plt.cm.Reds_r))
over_red.set_bad(alpha=0)


def display_roi(img, roi):
    'display an image with roi marked as an overlay'
    display(img, hold=True)
    plt.imshow(roi, cmap=over_red, alpha=0.4)


def normalize(img, norm_center, norm_radius=10):
    '''
    Returns the mean value and standard deviation of a circular roi in img
    ignoring masked pixels.
    
    Parameters
    ----------
    img : array
        image containing the pixel data
    norm_center : (x,y)
        center of the circular region of interest
    norm_radius : int (optional)
        radius of the circular ROI. Default: 10
        
    Returns
    -------
    m : float
        mean value of the ROI
    s : float
        standard deviation of the ROI   
    '''
    
    #ToDo: keine maske eingelesen?
    x,y = roi_circle(norm_center, norm_radius, img.mask)
    
    #normierung
    m = img.phase[y,x].mean()
    s = img.phase[y,x].std()
    
    #darstellung
    overlay = create_overlay(img.shape, x, y)
    display_roi(img, overlay)
    
    print 'mean phase is {:.3} +- {:.3}'.format(m,s)
    
    return m, s
   

def calc_displacement(phase, norm=0, deviation=0, G=20e-3, delta=20e-3):
    '''
    Calculates the displacement per pixel out of the given phase image.
    
    Parameters
    ----------
    phase : array
        image containing the phase data
    norm : float (optional)
        value to normalize the phase with. Default: 0
    deviation : float (optional)
        standard deviation of norm. Default: 0
    G : float (optional)
        maximal gradient's strength in T/m
    delta : float (optional)
        gradient's length in s
        
    Returns
    -------
    displacement : 2darray
        2-dimensional array containing the displacement per pixel
    '''
    if len(phase.shape) != 2:
        raise IndexError, 'phase has shape {}, expected 2-dimensional data'.format(phase.shape)
        
    const = G * delta * bvalue.gamma
    
    displacement = empty_like(phase)
    displacement['phase'] = - 2.0 * math.pi * (phase.phase - norm)/const * 1e6
    displacement['dev'] = - 2.0 * math.pi * deviation/const * 1e6
    displacement['mask'] = phase['mask'].copy()
    
    return displacement
    
    
        