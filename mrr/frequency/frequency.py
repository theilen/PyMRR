# -*- coding: utf-8 -*-

"""
Created on Fri Oct 30 15:29:20 2015

@author: Sebastian Theilenberg
"""

import numpy as np
from numpy.polynomial.polynomial import polyfit
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt


__version__ = '0.1'
# $Source$


def _find_extrema(t, data, w=None, s=None, t0=None, t1=None):
    "Find extrema by using a smoothing spline."
    if s is None and w is not None:
            s = t.size
    Sp = UnivariateSpline(t, data, w=w, k=4, s=s)
    roots = Sp.derivative(n=1).roots()
    if t0:
        roots = roots[roots > t0]
    if t1:
        roots = roots[roots < t1]
    return roots


def _calc_frequencies(extrema, weights=None):
    """
    Calculate the frequency from a number of extrema.

    Also works with multiple sets of extrema that are fitted individually
    with the same weights if extrema has shape (m, k) with k the number of
    sets. weights must have shape (m,).
    """
    extrema = np.asarray(extrema)
    n = extrema.shape[0]
    ret = polyfit(np.arange(n), extrema, deg=1, w=weights)
    return 0.5/ret[1]


def _match_extremes(extremes, match):
    if extremes.size < match.size:
        res = np.empty_like(match)
        res[:] = np.nan
        return res
    res = []
    for e_ in match:
        res.append(extremes[(np.abs(extremes-e_)).argmin()])
    return np.asarray(res)


def calculate_frequencies(t, data, N=1000, s=None, weights=None,
                          t0=None, t1=None):

    t = np.asarray(t)
    if t0 is None:
        t0 = t[0]
    if t1 is None:
        t1 = t[-1]
    # spline smoothing factor with weighted data
    if s is None:
        s = t.size
    if weights is None:
        weights = 1./data.dev**2.

    # find mean frequency
    extrema = _find_extrema(t, data, w=weights, t0=t0, t1=t1)
    n = extrema.size
    frequency = _calc_frequencies(extrema)

    # bootstrap
    bootsamples = np.empty((N, t.size))
    bootextremes = np.empty((N, n))
    for i in range(t.size):
        bootsamples[:, i] = np.random.normal(data.phase[i], data.dev[i], N)
    for i, sample in enumerate(bootsamples):
        extr_ = _find_extrema(t, sample, w=weights, s=s, t0=t0, t1=t1)
        if extr_.size != n:
#            plt.plot(t, sample, 'x')
#            plt.plot(t, UnivariateSpline(t, sample, weights, k=4, s=s)(t))
#            print extr_.size
            # TODO: will not work for n(extr_) < n(extrema)
            extr_ = _match_extremes(extr_, extrema)
        bootextremes[i] = extr_
    # handle nans
    if np.any(np.isnan(bootextremes)):
        indices = np.where(np.all(np.isfinite(bootextremes), axis=1))[0]
        print indices.shape
        bootextremes = bootextremes[indices]
        print "found nans", bootextremes.shape
        #return bootsamples
    # find error of extrema
    dextr = bootextremes.std(axis=0)
    # find error of frequencies
    freqs = _calc_frequencies(bootextremes.swapaxes(0, 1), weights=dextr)
    dfreq = freqs.std()

    # recalculate frequency using error of t0
    frequency = _calc_frequencies(extrema, weights=dextr)

    return np.array([frequency, dfreq])


# testing
calc_freq = _calc_frequencies
match_ = _match_extremes
find_extr = _find_extrema
