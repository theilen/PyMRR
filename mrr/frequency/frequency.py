# -*- coding: utf-8 -*-

"""
Created on Fri Oct 30 15:29:20 2015

@author: Sebastian Theilenberg
"""

import numpy as np
from numpy.polynomial.polynomial import polyfit
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import os


__version__ = '0.2'
# $Source$

__all__ = ["calculate_frequency", "average_frequencies"]


def _find_extrema(t, data, w=None, s=None, t0=None, t1=None,
                  show_plot=False, ploterr=None):
    "Find extrema by using a smoothing spline."
    if s is None and w is not None:
            s = t.size
    Sp = UnivariateSpline(t, data, w=w, k=4, s=s)
    roots = Sp.derivative(n=1).roots()
    if t0:
        roots = roots[roots > t0]
    if t1:
        roots = roots[roots < t1]

    if show_plot:
        fig = plt.figure()
        t_ = np.linspace(t[0], t[-1], t.size*10)
        plt.errorbar(t, data, yerr=ploterr, fmt='x')
        plt.plot(t_, Sp(t_))
        for r in roots:
            plt.axvline(r, color='black', ls='--')

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


def calculate_frequency(t, data, N=1000, s=None, weights=None,
                          t0=None, t1=None, show_plot=False):
    """
    Estimate the frequency of data and its error by bootstrapping.

    Parameters
    ----------
    t : 1Darray
        x-values of data (time)
    data : 1d-MRRArray
        the phase data.
    N : int
        number of bootsamples to draw to estimate the errors.
    s : float
        smoothing factor to be used for the spline representation of data.
    weights : 1darray
        weights to be used for the spline representation of data.
    t0 : float
        Extrema before this x-value will not be considered.
    t1 : float
        Extrema after this x-value will not be considered.
    show_plot : bool
        If True, generates a plot showing the estimated extrema and the spline
        representation.
    """
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
    extrema = _find_extrema(t, data.phase, w=weights, s=s, t0=t0, t1=t1,
                            show_plot=show_plot, ploterr=data.dev)
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
            # TODO: will not work for n(extr_) < n(extrema)
            extr_ = _match_extremes(extr_, extrema)
        bootextremes[i] = extr_
    # handle nans
    if np.any(np.isnan(bootextremes)):
        indices = np.where(np.all(np.isfinite(bootextremes), axis=1))[0]
        bootextremes = bootextremes[indices]
        print "found {} bad samples".format(N - indices.size)
    # find error of extrema
    dextr = bootextremes.std(axis=0)
    # find error of frequencies
    freqs = _calc_frequencies(bootextremes.swapaxes(0, 1),
                              weights=1./dextr**2.)
    dfreq = freqs.std()

    # recalculate frequency using errors of extrema
    frequency = _calc_frequencies(extrema, weights=1./dextr**2.)

    return np.array([frequency, dfreq])


def average_frequencies(t, data, error, N=1500, s=None, sweights=None,
                        limits=None, show_plot=False, print_infos=False):
    """
    Estimate the main frequency of every set in data and its error by
    bootstrapping.

    Parameters
    ----------
    t : list
        x-values for every data set
    data : list
        the phase data for every data set
    N : int
        number of bootsamples to draw to estimate the errors.
    s : float or list
        smoothing factor to be used for the spline representation of data.
        If s is a float, it is used for every data set. Otherwise provide a
        list of values with one value for every data set.
    sweights : list
        weights to be used for the spline representation of data.
    limits : list
        list of (t0,t1) pairs used to limit the x-range of the data.
        If only one pair is provided, it is used for every data set. If None,
        no limits are used.
    show_plot : bool
        If True, generates a plot showing the estimated extrema and the spline
        representation.
    print_infos : bool
        If True, additional infos are printed during excecution.
    """
    noofsets = len(data)
    # test data for consistency
    if len(t) != noofsets:
        raise ValueError("Expected x-values for every dataset")
    if weights is not None and len(weights) != noofsets:
        raise ValueError("Weights do not match data")
    if len(error) != noofsets:
        raise ValueError("Error does not match data")

    # set defaults
    if limits is None:
        use_limits = False
    else:
        use_limits = True
    if use_limits and len(limits) == 1:
        limits = limits*noofsets

    if s is not None and not hasattr(s, "__getitem__"):
        s = [s]*noofsets

    bootfrequencies = []
    samplefrequencies = []
    # find frequency of data sets and bootstrap
    for nset in range(noofsets):
        # get data
        t_ = t[nset]
        d_ = data[nset]
        e_ = error[nset]
        if weights is not None:
            w_ = weights[nset]
        else:
            w_ = None
        # set spline smoothing factor
        if s is None:
            # TODO: estimate s for weighted spline
            s_ = t_.size
        else:
            s_ = s[nset]
        # set limits
        if use_limits:
            t0, t1 = limits[nset]
        else:
            t0, t1 = t_[0], t_[-1]
        # find mean frequency
        extrema = _find_extrema(t_, d_, w=w_, s=s_, t0=t0, t1=t1,
                                show_plot=show_plot)
        n = extrema.size
        bootsamples = np.empty((N, t_.size))
        bootextremes = np.empty((N, n))
        for i in range(t_.size):
            bootsamples[:, i] = np.random.normal(d_[i], e_[i], N)
        for i, sample in enumerate(bootsamples):
            extr_ = _find_extrema(t_, sample, w=w_, s=s_, t0=t0, t1=t1,
                                  show_plot=False)
            if extr_.size != n:
                # TODO: will not work for n(extr_) < n(extrema)
                extr_ = _match_extremes(extr_, extrema)
            bootextremes[i] = extr_
        # handle nans
        if np.any(np.isnan(bootextremes)):
            indices = np.where(np.all(np.isfinite(bootextremes), axis=1))[0]
            bootextremes = bootextremes[indices]
            if print_infos:
                print "found {} bad samples in run {}/{}".format(
                    N - indices.size, nset+1, noofsets)

        dextr = bootextremes.std(axis=0)
        f_ = _calc_frequencies(extrema, weights=1./dextr**2.)
        bootfreqs = _calc_frequencies(bootextremes.swapaxes(0, 1),
                                      weights=1./dextr**2.)
        # store for later use
        bootfrequencies.append(bootfreqs)
        samplefrequencies.append(f_)

    # match number of bootsamples
    nvalid = min([b_.size for b_ in bootfrequencies])
    bootfrequencies = [b_[:nvalid] for b_ in bootfrequencies]
    bootfrequencies = np.vstack(bootfrequencies)  # (noofsets, nvalid)
    samplefrequencies = np.array(samplefrequencies)  # (noofsets,)
    # error of frequency per set
    dfreqs = np.std(bootfrequencies, axis=1)  # (noofsets,)
    fw = 1./dfreqs**2.
    # mean frequency per bootsample. shape: (nvalid,)
    mean_set_freqs = np.average(bootfrequencies, weights=fw, axis=0)
    # calculate mean frequency
    f = np.average(samplefrequencies, weights=fw)
    df = mean_set_freqs.std()

    return f, df


def testing(N=500, show_plot=True, print_infos=False):
    files = [f for f in os.listdir('.') if f.endswith('.npy')]
    rawdata = [np.load(f) for f in files]
    t = [d_[0] for d_ in rawdata]
    data = [d_[1] for d_ in rawdata]
    errs = [d_[2] for d_ in rawdata]
    weights = [1./e**2. for e in errs]
    return average_frequencies(t, data, errs, N=500, weights=weights,
                               show_plot=show_plot, print_infos=print_infos)
