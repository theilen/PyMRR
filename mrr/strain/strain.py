# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:28:53 2016

@author: theilenberg
"""

import numpy as np
from scipy.ndimage.filters import gaussian_filter
from scipy.integrate import trapz

from ..mrrcore import shrink_mask


def calc_strain(array, gauss=0.6, axis=1, n=2, dist=None):
    """
    Return the calculated phase strain.

    Strain is calculated as derivative (using np.gradient) of the phase of
    *array* in direction of *axis*. Prior to calculation, the phase is
    smoothed using a gaussian filter.

    The resulting strain will likely be overestimated at the borders. To
    account for this, the mask may be adjusted afterwards.

    Parameters:
    -----------
    array : MRRArray, input data
    gauss : float, Size of the gaussian smoothing kernel. Set to 0 to
        deactivate smoothing
    axis : int, axis along to calculate the strain
    n : int, number of pixels the mask is reduced

    Return:
    -------
    strain : MRRArray
    """
    strain = array.phase.copy()
    # strain[array.mask == False] = np.nan
    if array.ndim == 2:
        shrink_indices = [0, 1]
        sigma = [gauss, ]*2
    if array.ndim == 3:
        shrink_indices = [1, 2]
        sigma = [0, gauss, gauss]
    if dist is None:
        dist = [1, ]*array.ndim
    else:
        try:
            match = (len(dist) == array.ndim)
        except TypeError:
            dist = [1, ]*array.ndim
        else:
            if not match:
                raise TypeError("dist does not match array's dimensions!")

    strain = gaussian_filter(strain, sigma)
    strain = np.gradient(strain, edge_order=2, *dist)[axis]

    mask = shrink_mask(array.mask, shrink_indices[0], n=n)
    mask &= shrink_mask(array.mask, shrink_indices[1], n=n)
    res = array.copy()
    res['phase'] = strain
    res['mask'] = mask
    res['dev'] = 0.
    res['phase'][res.mask == False] = np.nan
    return res


def total_strain(strain, times, abs=True):
    """
    Integrate strain along first axis
    """
    assert strain.ndim == 3
    s = strain.phase.copy()
    if abs:
        s = np.abs(s)
    total = strain[0].copy()
    total['phase'] = trapz(s, times, axis=0)
    total['dev'] = 0.
    return total
