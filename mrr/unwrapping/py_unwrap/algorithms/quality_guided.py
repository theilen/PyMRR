# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 10:20:12 2014

@author: Sebastian Theilenberg
"""

__version__ = "1.1"
# $Source$

# Version history
# ===============
#
# changes in version 1.1
# - reimplemented unwrap_quality_guided to use adjoin.SortedList for better
#   performance


import numpy as np

from .adjoin import update_adjoin_list, refill_adjoin, SortedList
from .. import constants as cnst
from ..array_manipulation import unwrap_pixel


def unwrap_quality_guided(array, qualitymap, mask=None,
                          test=False):

    flags = cnst.init_flag(array)
    if np.any(mask):
        if not mask.shape == array.shape:
            raise ValueError("mask does not match array")
        else:
            flags[mask==False] |= cnst.BORDER

    max_list_size = array.shape[0] + array.shape[1]

    result = np.zeros_like(array)
    # sorted list with lowest quality as first element
    adjoin = SortedList(key=lambda x: qualitymap[x])

    # start with highest quality pixel
    sp = find_highest_quality(qualitymap, flags, mask)
    result[sp] = array[sp]
    flags[sp] |= cnst.UNWRAPPED
    update_adjoin_list(adjoin, sp, avoid=cnst.AVOID|cnst.UNWRAPPED,
                       flagarray=flags, qualitymap=qualitymap)

    if test: 
        counter = 0
    # unwrap using flood-fill
    while len(adjoin) > 0:
        if test: 
            counter += 1
        # unwrap highest quality pixel
        p = adjoin.pop()
        unwrap_pixel(p, flags, array, result)
        update_adjoin_list(adjoin, p, avoid=cnst.AVOID|cnst.UNWRAPPED,
                           flagarray=flags, qualitymap=qualitymap,
                           trim=True, max_size=max_list_size)

        if len(adjoin) == 0:
            # refill adjoin-list with POSTPONDED pixels
            refill_adjoin(adjoin, flags, max_list_size)

        if test:
            if counter % 50 == 0:
                np.save('./qg_unwrap{:0>5}.npy'.format(counter), flags)
                np.save('./qg_unwrap_adjoin{:0>5}.npy'.format(counter), adjoin)
                np.save('./qg_unwrap_adjoin_keys{:0>5}.npy'.format(counter),
                        adjoin.get_keys())

    return result


def find_highest_quality(qualitymap, flagarray=None, mask=None):
    if np.any(mask):
        highest = qualitymap[mask==True].max()
    else:
        highest = qualitymap.max()
    # find highest quality
    indices = np.where(np.isclose(qualitymap, highest))
    # find useful pixel
    for y, x in zip(*indices):
        if np.any(flagarray):
            if not flagarray[y,x] & cnst.AVOID:
                return (y, x)
        else:
            return (y, x)


       
    
    
    
    