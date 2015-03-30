# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 11:25:26 2014

@author: Sebastian Theilenberg
"""

__version__ = "1.1"
# $Source$

# Version history
# ===============
#
# version 1.1
# - added keyword shape to find_unwrapped_neighbour
# - adjusted unwrap_pixel to new find_unwrapped_neighbour


from .. import constants as cnst
from ..array_manipulation import wrap, get_neighbours


def find_unwrapped_neighbour(pixel, flag_array, shape=None):
    '''Find one adjoining unwrapped pixel.'''
    if not shape:
        shape = flag_array.shape
    for by, bx in get_neighbours(*pixel, shape=shape):
        if flag_array[by,bx] & cnst.UNWRAPPED:
            return (by, bx)


def unwrap_pixel(p, flag_array, array, result_array):
    '''Unwraps pixel p and sets the flag_array to UNWRAPPED.'''
    nb = find_unwrapped_neighbour(p, flag_array, shape=flag_array.shape)
    if not nb:
        return False
    diff = wrap(array[p[0],p[1]] - array[nb[0],nb[1]])
    result_array[p[0],p[1]] = result_array[nb[0],nb[1]] + diff
    flag_array[p[0],p[1]] |= cnst.UNWRAPPED
    return True
