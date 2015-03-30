# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 13:26:42 2014

@author: Sebastian Theilenberg
"""

__version__ = "1.1"
# $Source$

#Version history
#===============
#
#version 1.1
#- reimplemented to use .adjoin.SortedList

import numpy as np

from .. import constants as cnst
from .residues import find_residues, get_charge
from ..array_manipulation import find_flag_neighbour#, thin_mask
from .adjoin import update_adjoin_list, refill_adjoin, SortedList
from .py_goldstein import unwrap_around_branch_cuts


#import matplotlib.pyplot as plt

__all__ = ['unwrap_mask_cut', 'grow_mask_cuts', 'thin_mask_cuts']


def unwrap_mask_cut(array, quality_map, mask=None):
    """
    Unwrap two-dimensional phase-data using the mask-cut-algorithm.
    
    Parameters
    ----------
    array : numpy.ndarray
        phase data to be unwrapped. The data is supposed to be two-dimensional
        and wrapped to 0.0 to 1.0.
    mask : numpy.ndarray, optional
        mask to mark regions to avoid during the unwrapping process. Optional, 
        but not providing a mask in the presens of noisy regions may 
        significantly slow down performance. The shape of mask has to match 
        the shape of array!
        
    Returns
    -------
    result : numpy.ndarray
        The unwrapped phase-data
    flags : bitflag-array
        Bitflags used during th            #refill adjoin-list with POSTPONDED pixels
            qt_ = refill_adjoin(adjoin, flags, qualitymap, max_list_size, 
                                sort_reverse=False)
            if qt_: 
                qual_thresh = qt_e unwrapping process. The bitflags used are 
        defined in x.constants.
    num_pieces : int
        Number of seperated pieces unwrapped. If this does not equal to one,
        some regions in the unwrapped data were disconnected and the unwrapped
        phase data may be ambiguous.
    """
    flags = cnst.init_flag(array)
    if np.any(mask):
        if not mask.shape == array.shape:
            raise ValueError("mask does not match array")
        else:
            flags[mask==False] |= cnst.BORDER
            
    find_residues(array, flags)

    grow_mask_cuts(flags, quality_map)
    thin_mask_cuts(flags)

    result, num_pieces = unwrap_around_branch_cuts(array, flags)
    
    return result, flags, num_pieces
    

def grow_mask_cuts(flagarray, qual):
    """
    Creates mask cuts based on a growing approach growing along lowest quality
    pixels as specified in qual.
    """
    maxsize = flagarray.shape[0] + flagarray.shape[1]
    residues = zip(*np.where(flagarray & cnst.RESIDUE))
    for res in residues:
        #print 'starting with residue {} (qt = {:.4f})'.format(res, qual_thresh)
        if flagarray[res] & cnst.BALANCED:
            continue
        charge = get_charge(flagarray[res])
        flagarray[res] |= cnst.CUT
        adjoin = SortedList(key=lambda x: 1./qual[x]) #descending quality
        update_adjoin_list(adjoin, res, avoid=cnst.CUT, flagarray=flagarray, 
                           qualitymap=qual, trim=True, max_size=maxsize)
        #print 'adjoin_list: {}, (qt = {:.4f})'.format(adjoin[-3:], qual_thresh)
        while len(adjoin) > 0 and charge != 0:
            p = adjoin.pop()
            #print 'popped pixel {} [{:.4f}] (qt = {})'.format(p, qual[p], qual_thresh)
            flagarray[p] |= cnst.CUT
            if (flagarray[p] & cnst.RESIDUE and 
                not flagarray[p] & cnst.BALANCED):
                #print 'unbalanced resiude!'
                charge += get_charge(flagarray[p])
                flagarray[p] |= cnst.BALANCED
            if flagarray[p] & cnst.BORDER:
                #print 'border!'
                charge = 0
            update_adjoin_list(adjoin, p, avoid=cnst.CUT, flagarray=flagarray, 
                               qualitymap=qual, trim=True, max_size=maxsize)
            #print 'adjoin_list: {}, [{}], (qt = {:.4f})'.format(adjoin[-3:], qual[zip(*adjoin[-3:])], qual_thresh)
            if len(adjoin) == 0:
                #refill adjoin-list with POSTPONDED pixels
                refill_adjoin(adjoin, flagarray, maxsize)
                        

def thin_mask_cuts(flagarray):
    """
    Thins mask cuts without changing the connectivity of the mask.
    """
    while True:
        changed = False
        cuts = zip(*np.where(flagarray & cnst.CUT))
        for pixel in cuts:
            if not find_flag_neighbour(pixel, flagarray, cnst.CUT, unset=True):
                continue #only neighbours that are marked CUT
            if find_flag_neighbour(pixel, flagarray, cnst.RESIDUE|cnst.BORDER):
                continue #neighbour that is RESIDUE
            if mask_connectivity(pixel, flagarray):
                flagarray[pixel] &= ~cnst.CUT
                changed = True
        if not changed:
            return
    

def mask_connectivity(pixel, flagarray):
    '''Returns False, if removing the pixel would change the mask-connectivity'''
    regions = [[[-1,-1,0],[0,1,1]],
               [[0,1,1],[1,1,0]],
               [[1,1,0],[0,-1,-1]],
               [[0,-1,-1],[-1,-1,0]]]
    counter = 0
    for reg in regions:
        indices = [ [pixel[0]+i for i in reg[0]], [pixel[1]+i for i in reg[1]] ]
        try:
            cuts = flagarray[indices]# a,b,c = zip(*indices)
        except IndexError:
            return False
        if cuts[0]:
            if not cuts[2]: counter += 1
        else:
            if cuts[1]:
                counter += 1
                if not cuts[2]: counter += 1
            elif cuts[2]: counter += 1
    return False if counter>2 else True
    
                
            

        
        
