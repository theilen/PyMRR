# -*- coding: utf-8 -*-
"""
No docstring yet.

@author: Sebastian Theilenberg
"""

# version history
# ===============
#
# changes in version 1.1
# - reimplemented DerivativeVariance and _phase_derivative_window for better
#   performance
#

__version__ = "1.1"
# $Source$


import numpy as np

from ..array_manipulation import wrap


def derivative_variance(array, k=3):
    '''
    Calculates the Phase Derivative Variance of an array using a k x k window.

    Parameters
    ----------
    array : np.ndarray
        Array containing the phase data
    k : int
        window size, should be uneven
    '''
    result = np.empty(array.shape)
    diff = [wrap(np.diff(array, axis=a)) for a in [0, 1]]

    it = np.nditer(result, flags=['multi_index'])
    while not it.finished:
        index = it.multi_index
        result[index] = _phase_derivative_window(diff, index, k//2)/k**2.
        it.iternext()

    return result


def _phase_derivative_window(diff, index, k):
    '''
    Calculates the phase derivative variance at one pixel.

    Uses the supplied derivatives <diff> in a size <k> neighbourhood of index.
    Since the shape of this neighbourhood is dependent on the axis of the
    derivatives, <axis> has to be supplied as well.
    '''
    y, x = index
    result = 0
    for axis in [0, 1]:
        try:
            diff_win = diff[axis][y-k:y+k+axis, x-k:x+k+(axis % 1)]
        except IndexError:
            y_max, x_max = diff[axis].shape[1:]
            diff_win = diff[axis][
                slice(max(0, y-k), min(y_max, y+k+axis)),
                slice(max(0, x-k), min(x_max, x+k+(axis % 1)))
                ]
        mean = diff_win.mean()
        result += np.power(diff_win - mean, 2.).sum() / k**2.
    return np.sqrt(result)
