# -*- coding: utf-8 -*-
"""
Helper-functions based on bitflags.
"""

___version__ = "1.1"
# $Source$

# Version history
# ===============
#
# version 1.1
# - removed get_charge (moved to ..algorithms.residues)
# - added __all__


import numpy as np

from .. import constants as cnst
from .pixel import create_box, distance, get_neighbours


__all__ = ['set_border',
           'find_nearest_border',
           'find_flag_neighbour',
           'thin_mask']


def set_border(flags):
    '''Set the border pixels of a bitflag-array to BORDER.'''
    flags[0,:] |= cnst.BORDER
    flags[-1,:] |= cnst.BORDER
    flags[:,0] |= cnst.BORDER
    flags[:,-1] |= cnst.BORDER


def thin_mask(flagarray, n=1, flag=cnst.BORDER):
    i = 0
    while i < n:
        i += 1
        remove = []
        changed = False
        border = zip(*np.where(flagarray & flag))
        for pixel in border:
            if not find_flag_neighbour(pixel, flagarray, flag, unset=True):
                continue  # only neighbours that are marked as border
            else:
                remove.append(pixel)
        if remove:
            changed = True
            flagarray[zip(*remove)] &= ~flag
        if not changed:
            break


def find_nearest_border(flag_array, pixel):
    """Return the the nearest border-pixel to pixel."""
    if flag_array[pixel[0],pixel[1]] & cnst.BORDER:
        return pixel
    n = 3
    while True:
        indices = zip(*create_box(*pixel, n=n))  # zip to create index-arrays
        borders = np.where(flag_array[indices] & cnst.BORDER)[0]
        if np.any(borders):
            # all border pixels as tuple
            bps = [zip(*indices)[i] for i in borders]
            nearest = min(bps, key=lambda x: distance(x, pixel))
            return nearest
        n += 2


def find_flag_neighbour(pixel, flagarray, flag, unset=False):
    '''
    Return one pixel of the four adjoining neighbours that has the specified
    flag set.

    If unset is True, return one neighbouring pixel that has the
    specified flag not set.
    Return False if no corresponding pixel is found.
    '''
    pixel = get_neighbours(*pixel, shape=flagarray.shape)
    if not unset:
        for p in pixel:
            if flagarray[p] & flag:
                return p
    else:
        for p in pixel:
            if not flagarray[p] & flag:
                return p
    return False
