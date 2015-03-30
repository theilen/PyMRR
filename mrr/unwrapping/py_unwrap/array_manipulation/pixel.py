# -*- coding: utf-8 -*-
"""
Helper-functions working with pixel-indices.

@author: Sebastian Theilenberg
"""

__version__ = "1.1"
# $Source$


# version history
# =================
# version 1.1
# -------------
# - added keyword shape to get_neighbours


import numpy as np


def wrap(phase):
    '''Wrap a scalar value or an entire array to -0.5 <= a < 0.5.'''
    if hasattr(phase, '__len__'):
        phase = phase.copy()
        phase[phase > 0.5] -= 1.
        phase[phase <= -0.5] += 1.
        return phase
    else:
        if phase > 0.5:
            phase -= 1.
        elif phase <= -0.5:
            phase += 1.
    return phase


def distance(pixel1, pixel2):
    """Return the euklidean distance of two pixels."""
    result = 0
    for i in range(len(pixel1)):
        result += (pixel1[i] - pixel2[i])**2.
    return np.sqrt(result)


def get_neighbours(y, x, shape=None):
    '''Return a list of the four direct neighbours of (y,x), if existing.'''
    def valid_pixel(p, y_max, x_max):
        return (0 <= p[0] < y_max) and (0 <= p[1] < x_max)

    if not shape:
        shape = (999, 999)

    pixel = [(y-1,x), (y+1,x), (y,x-1), (y,x+1)]
    return [p for p in pixel if valid_pixel(p, *shape)]


def create_box(y, x, n, full=False):
    """
    Return the indices (y,x) of the border-pixels of a squared box of size n
    centered to (y,x). If full==True, returns all pixels in that box except
    the center-one.
    """
    if full:
        indices = []
        for i in range(3, n+1, 2):
            indices.extend(create_box(y, x, i))
    else:
        indices = zip([y-n//2]*n, range(x-n//2, x+n//2+1)) \
            + zip([y+n//2]*n, range(x-n//2, x+n//2+1)) \
            + zip(range(y-n//2+1, y+n//2), [x-n//2]*(n-2)) \
            + zip(range(y-n//2+1, y+n//2), [x+n//2]*(n-2))
    # Remove pixel with negative coordinates
    indices = [i for i in indices if not (i[0] < 0 or i[1] < 0)]
    indices.sort(key=lambda p: distance(p, (y, x)))
    return indices
