# -*- coding: utf-8 -*-
"""
Created on Sun Nov 09 11:56:02 2014

@author: Sebastian Theilenberg
"""

__version__ = "1.0"
# $Source$


__all__ = ['POS_RES',
           'NEG_RES',
           'CUT',
           'BORDER',
           'BALANCED',
           'UNWRAPPED',
           'ACTIVE',
           'POSTPONED',
           'RESIDUE',
           'AVOID',
           'init_flag'
           ]


import numpy as np


POS_RES = 0x01      # 1
NEG_RES = 0x02      # 2
CUT = 0x04          # 4
BORDER = 0x08       # 8
BALANCED = 0x10     # 16
UNWRAPPED = 0x20    # 32
ACTIVE = 0x40       # 64
POSTPONED = 0x80    # 128

RESIDUE = POS_RES | NEG_RES
AVOID = BORDER | CUT | POSTPONED


def init_flag(array):
    return np.zeros_like(array, dtype=np.ubyte)
