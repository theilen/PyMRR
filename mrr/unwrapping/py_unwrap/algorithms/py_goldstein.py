# -*- coding: utf-8 -*-
"""
Implementation of Goldstein's, Zebker's and Werner's algorithm for unwrapping
a two-dimensional phase array.

@author: Sebastian Theilenberg
"""

__version__ = "1.1"
# $Source$

#Version history
#===============
#
#version 1.1
#- added __all__
#- removed update_adjoin
#- switched to use .adjoin.update_adjoin_list


import numpy as np

from .residues import find_residues, remove_dipoles, get_charge, \
                        place_branch_cut
from .adjoin import update_adjoin_list, SortedList
from ..array_manipulation import find_nearest_border, set_border, \
                                    create_box, unwrap_pixel, distance
from .. import constants as cnst


__all__ = ['unwrap_py_gold', 'unwrap_around_branch_cuts']


def unwrap_py_gold(array, mask=None, dipoles=True, unwrap_by_distance=False,
                   test=False):
    """
    Unwrap two-dimensional phase-data using Goldstein's algorithm.
    
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
    diploles : bool, optional
        If set to True, dipoles in the residues are removed as a preprocessing
        to goldstein's algorithm. (Default: True)
        
    Returns
    -------
    result : numpy.ndarray
        The unwrapped phase-data
    flags : bitflag-array
        Bitflags used during the unwrapping process. The bitflags used are 
        defined in x.constants.
    num_pieces : int
        Number of seperated pieces unwrapped. If this does not equal to one,
        some regions in the unwrapped data were disconnected and the unwrapped
        phase data may be ambiguous.
    """
    #prepare flags
    flags = cnst.init_flag(array)
    if np.any(mask):
        if not mask.shape == array.shape:
            raise ValueError("mask does not match array")
        else:
            flags[mask==False] |= cnst.BORDER
    set_border(flags)

    #residues
    find_residues(array, flags)    
    #branch_cuts
    if dipoles: remove_dipoles(flags)
    generate_branch_cuts(flags)
    #unwrap    
    result, num_pieces = unwrap_around_branch_cuts(
                            array, flags, test=test, 
                            unwrap_by_distance=unwrap_by_distance
                            )
    
    if test:
        np.save('./pg_flags.npy', flags)
    
    return result, flags, num_pieces
    
    
    
def unwrap_around_branch_cuts(array, flag_array, unwrap_by_distance=False,
                              test=False):
    '''Unwrapp the array avoiding crossing branch cuts.'''
    result = np.zeros_like(array)
    if unwrap_by_distance:
        adjoin = SortedList()
    else: 
        adjoin = []
    
    num_pieces = 0
    maxl = 0
        
    while True:
        if test: counter = 0
        num_pieces += 1 #amount of seperated pieces
        #Find starting pixel and unwrap it
        sp = zip(*np.where(flag_array & (cnst.BORDER | cnst.UNWRAPPED) == False))[0]
        #print 'starting region {} with pixel {}'.format(num_pieces, sp)
        result[sp[0],sp[1]] = array[sp[0],sp[1]]
        flag_array[sp[0],sp[1]] |= cnst.UNWRAPPED
        
        if unwrap_by_distance:
            adjoin.key_function = lambda x: 1./distance(x, sp)#reverse-sorted by distance to sp
        update_adjoin_list(adjoin, sp, avoid=cnst.AVOID|cnst.UNWRAPPED, 
                           flagarray=flag_array)
        while len(adjoin) > 0:
            if test: counter += 1
            #unwrap pixel in adjoin_list
            unwrap_next_pixel(adjoin, flag_array, array, result)
            maxl = max(maxl, len(adjoin))
            if test:
                if counter%50 == 0:
                    np.save('./pg_unwrap{:0>5}.npy'.format(counter), flag_array)
                    np.save('./pg_unwrap_adjoin{:0>5}.npy'.format(counter), adjoin)
                    np.save('./pg_unwrap_adjoin_keys{:0>5}.npy'.format(counter),
                            adjoin.get_keys())
        #check for additional separated pieces
        if not np.any(np.where(flag_array & (cnst.BORDER | cnst.CUT | cnst.UNWRAPPED) == False)):
            break
    #Unwrap branch cut pixel
    for sp in zip(*np.where(flag_array & cnst.CUT)):
        if flag_array[sp[0],sp[1]] & (cnst.UNWRAPPED | cnst.BORDER):
            continue
        if unwrap_pixel(sp, flag_array, array, result):
            update_adjoin_list(adjoin, sp, avoid=cnst.AVOID|cnst.UNWRAPPED,
                               flagarray=flag_array)
            while len(adjoin) > 0:
                unwrap_next_pixel(adjoin, flag_array, array, result)
                maxl = max(maxl, len(adjoin))
    
    return result, num_pieces
    

def unwrap_next_pixel(pixel_list, flag_array, array, result_array):
    p = pixel_list.pop()
    if not unwrap_pixel(p, flag_array, array, result_array):
        return False
    update_adjoin_list(pixel_list, p, avoid=cnst.AVOID|cnst.UNWRAPPED,
                       flagarray=flag_array)
    return True
    

def generate_branch_cuts(flag_array, maxboxsize=None):
    '''Generate branch cuts to balance all residues.'''
    #Set maximum size of box to search for neighbouring residues
    if maxboxsize == None:
        maxboxsize = min(*flag_array.shape)//3
        
    residues = zip(*np.where(flag_array & cnst.RESIDUE))
    #Iteration over all unbalanced residues
    for ry, rx in residues:
        if flag_array[ry,rx] & cnst.BALANCED: 
            continue
        charge = get_charge(flag_array[ry,rx])
        flag_array[ry,rx] |= (cnst.ACTIVE | cnst.BALANCED)
        
        for n in range(3,maxboxsize+1,2):
            active_pixels = zip(*np.where(flag_array & cnst.ACTIVE))

            for ay, ax in active_pixels:
                #Check all pixels in the current box. If boxsize is 
                #increased, only new pixels need to be considered, 
                #except if the active pixel was not visited before.
                #print 'searching around pixel {},{} with boxsize {}'.format(ay,ax,n)
                full_box = False if flag_array[ay,ax] & cnst.BALANCED else True

                #print 'box: {}'.format(create_box(ay, ax, n=n, full=full_box))
                for by, bx in create_box(ay, ax, n=n, full=full_box):
                    #found border
                    if flag_array[by,bx] & cnst.BORDER: #expects border-pixels to be marked BORDER
                        place_branch_cut((ay,ax), (by,bx), flag_array)
                        charge = 0
                    elif (flag_array[by,bx] & cnst.RESIDUE) and not (flag_array[by,bx] & cnst.ACTIVE):
                        #non-active residue
                        if not (flag_array[by,bx] & cnst.BALANCED):
                            #unbalanced residue
                            charge += get_charge(flag_array[by,bx])
                            flag_array[by,bx] |= cnst.BALANCED
                        flag_array[by,bx] |= cnst.ACTIVE
                        place_branch_cut((ay,ax), (by,bx), flag_array)
                    if charge == 0: break
                if charge == 0: break
            if charge == 0: break

        if charge != 0:
            #Did not find residues to balance with -> place branch-cut to border
            nearest_border = find_nearest_border(flag_array, (ry,rx))
            place_branch_cut((ry,rx), nearest_border, flag_array)
        active_pixels = np.where(flag_array & cnst.ACTIVE)
        flag_array[active_pixels] &= ~cnst.ACTIVE
        
        


        

    
                
