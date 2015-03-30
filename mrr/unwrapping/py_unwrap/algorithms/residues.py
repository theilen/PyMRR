# -*- coding: utf-8 -*-
"""
Created on Sun Nov 09 12:38:09 2014

@author: Sebastian Theilenberg
"""

__version__ = "1.1"
# $Source$

# Version history
# ===============
#
# version 1.1
# - added get_charge


import numpy as np

from .. import constants as cnst
from ..array_manipulation import wrap, create_box


def get_charge(flag):
    'Returns the charge of a pixel by a given flag'
    if not flag & cnst.RESIDUE:
        return 0
    else:
        return -1 if (flag & cnst.NEG_RES) else 1


def find_residues(array, flag_array):
    # calculate wrapped differences
    dx = wrap(np.diff(array, axis=1))
    dy = wrap(np.diff(array, axis=0))
    # identify residues
    it = np.nditer(flag_array, op_flags=['readwrite'], flags=['multi_index'])
    for flag in it:
        if flag & cnst.BORDER:
            continue
        y, x = it.multi_index
        try:
            # sum of wrapped differences
            s = dy[y,x] + dx[y+1,x] - dy[y,x+1] - dx[y,x]
        except IndexError:
            continue
        else:
            # set flags
            if np.sign(s) == 1:
                flag[...] |= cnst.POS_RES
            elif np.sign(s) == -1:
                flag[...] |= cnst.NEG_RES


def remove_dipoles(flag_array):
    residues = zip(*np.where(flag_array & cnst.RESIDUE))
    for ry, rx in residues:
        if flag_array[ry, rx] & cnst.BALANCED:
            continue
        charge = get_charge(flag_array[ry, rx])
        for p in create_box(ry, rx, 3):  # get_neighbours(ry,rx):
            if ((flag_array[p[0], p[1]] & cnst.RESIDUE)
                    and not (flag_array[p[0], p[1]] & cnst.BALANCED)):
                if get_charge(flag_array[p[0], p[1]]) + charge == 0:
                    place_branch_cut(p, (ry, rx), flag_array)
                    flag_array[p[0], p[1]] |= cnst.BALANCED
                    flag_array[p[0], p[1]] &= ~cnst.RESIDUE
                    flag_array[ry, rx] |= cnst.BALANCED
                    flag_array[ry, rx] &= ~cnst.RESIDUE


def place_branch_cut(pixel1, pixel2, flag_array):
    '''Place a branch cut between two pixels.'''
    # Ensure working from left to right
    if pixel1[1] > pixel2[1]:
        pixel1, pixel2 = pixel2, pixel1
    # determine wether working up->down or vice versa
    if pixel1[0] > pixel2[0]:
        # upward
        pixel = [pixel1[0]+1, pixel1[1]]  # adjust y for first increment
        finish = [pixel2[0]+1, pixel2[1]]
        dy = -1
    else:
        # downward
        pixel = list(pixel1)  # no adjustment needed
        finish = list(pixel2)
        dy = 1
    # create branch-cut
    while True:
        if pixel[0] != finish[0]:
            pixel[0] += dy
        if pixel[1] != finish[1]:
            pixel[1] += 1
        flag_array[pixel[0], pixel[1]] |= cnst.CUT
        if pixel == finish:
            return
