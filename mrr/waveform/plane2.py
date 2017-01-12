# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 16:17:32 2016

@author: theilenberg
"""

import numpy as np
import math


try:
    from uncertainties import unumpy, ufloat
    from uncertainties.umath import acos, atan2, degrees
except ImportError:
    UNCERT = False
    from math import acos, atan2, degrees
else:
    UNCERT = True


def _uncert_array(x, sx, syst):
    """
    Returns a uarray if possible, else a np.array ignoring the uncertainties.
    """
    if UNCERT is True:
        if np.any(sx):
            assert x.shape == sx.shape
            s = np.abs(sx)
        else:
            s = 0
        arr = unumpy.uarray(x, s)
        # bulid array of systematic uncertainties
        if np.any(syst):
            systarr = unumpy.uarray([ufloat(0., s_, 'sytematic')
                                    for s_ in np.abs(syst).tolist()])
#            systarr = []
#            for s in syst:
#                systarr.append(ufloat(0., abs(s), tag='systematic'))
            arr += np.array(systarr)
    else:
        arr = np.array(x)
    return arr


def _uncert_nominal(x):
    "Returns the nominal value."
    if UNCERT is True:
        arr = unumpy.nominal_values(x)
    else:
        arr = np.asarray(x)
    if arr.size == 1:
        arr = arr.item()
    return arr


def _uncert_deviation(x):
    "Returns the standard deviation."
    if UNCERT is True:
        arr = unumpy.std_devs(x)
    else:
        arr = np.zeros_like(x)
    if arr.size == 1:
        arr = arr.item()
    return arr


def _uncert_systematic(x):
    "Returns the systematic error component."
    if not hasattr(x, '__len__'):
        x = [x]
    res = []
    for item in x:
        syst_err = math.sqrt(sum(
            e**2.
            for (v, e) in item.error_components().items()
            if v.tag == "systematic"))
        res.append(syst_err)
    if len(res) == 1:
        return res[0]
    else:
        return np.array(res)


class Plane(object):

    def __init__(self, t, p, sd=None, syst=None):
        """
        t : array
            timepoitns
        p : array
            (len(t), 3, 3)-shaped array of coordinates
        sd : array
            standard deviation of p
        syst : array
            systematic error on p
        """
        assert t.ndim == 1
        assert p.shape == (3, 3, t.shape[0])
        if np.any(sd):
            assert p.shape == sd.shape
        if np.any(syst):
            assert p.shape == syst.shape

        self._times = t
        if UNCERT:
            self._points = _uncert_array(p, sd, syst)
        else:
            raise ValueError
        self._has_run = False

    def run(self):
        'Calculate plane coefficients for all timepoints.'
        # calculate normal vectors
        p1, p2, p3 = self._points
        for p in ['p1', 'p2', 'p3']:
            exec("{0} = np.rollaxis({0}, -1, 0)".format(p))
        assert p1.shape[0] == self._times.shape[0]
        n = np.cross(p1-p3, p2-p3)
        norm = np.sqrt(np.sum(_uncert_nominal(n)**2., axis=-1))
        n = (n.T/norm).T
        # correct sign of normal vectors to n*p3 > 0
        n[np.array(
            [np.dot(n[i], p3[i]) for i in xrange(n.shape[0])]) < 0] *= -1
        d = np.array([np.dot(n[i], p3[i]) for i in xrange(n.shape[0])])
        self._n = n
        self._d = d
        self._has_run = True

    def x(self, y, z, return_time=True):
        "Return the x coordinate of pixel (x, y, z) at at all times."
        res = self._calc_coordinates(0, y, z)
        if return_time:
            return self._times, res
        else:
            return res

    def y(self, x, z, return_time=True):
        "Return the y coordinate of pixel (x, y, z) at at all times."
        res = self._calc_coordinates(1, x, z)
        if return_time:
            return self._times, res
        else:
            return res

    def z(self, x, y, return_time=True):
        "Return the z coordinate of pixel (x, y, z) at at all times."
        res = self._calc_coordinates(2, x, y)
        if return_time:
            return self._times, res
        else:
            return res

    def _calc_coordinates(self, index, a, b):
        "Calculate the coordinate index with fixed coordinates a and b."
        indices = range(3)
        assert index in indices
        indices.remove(index)
        # ToDo sort list?
        x = (self._d - a*self._n[:, indices[0]] - b*self._n[:, indices[1]])
        x /= self._n[:, index]
        # TODO: implement single t
        return x

    def _orient(self, axis, deg=True):
        axis = axis/np.linalg.norm(axis, ord=2)
        angle = np.array(
            [acos(np.dot(self._n[i], axis)) for i in xrange(self._n.shape[0])]
            )
        if deg:
            angle = np.array(
                [degrees(angle[i]) for i in xrange(angle.shape[0])]
                )
        return angle

    def orientation(self, axis, deg=True):
        "Return the angle between the normal vector and axis."
        axis = np.array(axis)
        assert axis.size == 3
        return self._orient(axis, deg=deg)

    def inclination(self, deg=True):
        "Return the inclination with respect to the negative y axis."
        return self._orient(axis=[0, 1, 0], deg=deg)

    def azimuth(self, deg=True):
        "Return the azimuth of the plane with respect to the positive x-axis."
        az = np.array(
            [atan2(self._n[i, 2], self._n[i, 0])
                for i in xrange(self._n.shape[0])])
        az = (az + 2.*np.pi) % (2.*np.pi)
        if deg:
            az = np.array([degrees(az[i]) for i in xrange(az.shape[0])])
        return az
