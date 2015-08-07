# -*- coding: utf-8 -*-
"""
Created on Fri Mar 06 20:27:30 2015

@author: Sebastian Theilenberg
"""

import numpy as np
from scipy.interpolate import interp1d
import math
from math import pi

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
        arr = unumpy.uarray(x, sx)
        # bulid array of systematic uncertainties
        systarr = []
        for s in syst:
            systarr.append(ufloat(0., s, tag='systematic'))
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


class moving_plane(object):

    def __init__(self, p1, p2, p3, sd1=None, sd2=None, sd3=None,
                 syst1=None, syst2=None, syst3=None):
        """
        Create a moving plane based on the trajectories of three basepoints.

        If the uncertainties package is present, moving_plane can calculate
        the standard deviations of the calculated values, provided some
        standard deviations for the trajectories.

        p1, p2, p3 : array_like
            trajectories of the three basepoints (t, x, y, z)
        sd1, sd2, sd3 : array_like
            standard deviations of the points p1...p3. (dx, dy, dz)
        syst1, syst2, syst3 : array_like
            systematic deviations of the points p1...p3. (dx, dy, dz)
        """
        points = [p1, p2, p3]
        self._points = []
        for i, p in enumerate(points):
            if not len(p) == 4:
                raise ValueError("expected size 4 array")
            interp = interp1d(p[0], p[1:], copy=True,
                              bounds_error=True, assume_sorted=True)
            self._points.append(interp)
        self._sd = []
        for i, sd in enumerate([sd1, sd2, sd3]):
            # set sd to zero if not present
            if not np.any(sd):
                sd = np.zeros((3, len(points[i])))
            if not len(sd) == 3:
                raise ValueError("expected size 3 array for sd")
            interp = interp1d(points[i][0], sd, copy=True,
                              bounds_error=True, assume_sorted=True)
            self._sd.append(interp)
        self._syst = []
        for i, sd in enumerate([syst1, syst2, syst3]):
            # set sd to zero if not present
            if not np.any(sd):
                sd = np.zeros((3, len(points[i])))
            if not len(sd) == 3:
                raise ValueError("expected size 3 array for syst")
            interp = interp1d(points[i][0], sd, copy=True,
                              bounds_error=True, assume_sorted=True)
            self._syst.append(interp)

    def __call__(self, t, std=False):
        'Returns the planes coefficients (n1, n2, n3, d) at time t.'
        return self._return_result(self._create_plane(t))

    def _return_result(self, result, std=False, seperate_syst=False):
        n = _uncert_nominal(result)
        if std:
            s = _uncert_deviation(result)
            if seperate_syst:
                return n, s, _uncert_systematic(result)
            else:
                return n, s
        else:
            return n

    def _get_points(self, t, index=None):
        points = []
        for p in range(3):
            if index:
                if not p == index:
                    continue
            data = _uncert_array(
                self._points[p](t), self._sd[p](t), self._syst[p](t)
                )
            points.append(data)
        if len(points) == 1:
            points = points[0]
        return points

    def _create_plane(self, t):
        "Creates the plane for a single time t."
        n = self._construct_n(t)
        d = np.dot(n, self._get_points(t, 3))
        plane = np.array([n.tolist()] + [d])
        return plane

    def _construct_n(self, t):
        "Construct a normalized normal vector to p1-p3 and p2-p3 at time t."
        p1, p2, p3 = self._get_points(t)
        n = np.cross(p1-p3, p2-p3)
        norm = np.dot(n, n)**0.5
        return n/norm

    def _calc_param(self, t, index, coords):
        """
        Return the coordinate index at time t with given two coordinates
        coords. coords in the form [x, y, z] with only two of them given.
        """
        try:
            len(t)
        except TypeError:
            t = [t]
        finally:
            t = np.array(t)

        result = []
        for i in range(t.size):
            coords_ = coords[:]
            t_ = t[i]
            n = self._construct_n(t_)
            res = np.dot(self._get_points(t_, 2), n)  # constant of normal form
            for c in range(3):
                if c == index:
                    continue
                else:
                    res -= n[c] * coords_.pop(0)
            result.append(res/n[index])
        return np.array(result)

    def x(self, t, y, z, std=False, syst=False):
        "Returns the x-component of a point in the yz-plane at time t."
        return self._return_result(self._calc_param(t, 0, [y, z]), std, syst)

    def y(self, t, x, z, std=False, syst=False):
        "Returns the y-component of a point in the xz-plane at time t."
        y = self._calc_param(t, 1, [x, z])
        return self._return_result(y, std, syst)

    def z(self, t, x, y, std=False, syst=False):
        "Returns the z-component of a point in the xy-plane at time t."
        return self._return_result(self._calc_param(t, 2, [x, y]), std, syst)

    def inclination(self, t, deg=True, std=False, syst=False):
        "Return the inclination of n with respect to the negative y-axis."
        orient = self._orient(t, axis=np.array([0, -1, 0]), deg=deg)
        return self._return_result(orient, std, syst)

    def orientation(self, t, axis, deg=True, std=False, syst=False):
        "Return the angle of the normal vector n to axis at time t."
        orient = self._orient(t, axis, deg)
        return self._return_result(orient, std, syst)

    def _orient(self, t, axis, deg=True):
        "Return the angle of the nromal vector n to axis at time t."
        n = self._construct_n(t)
        axis = axis/np.linalg.norm(axis, ord=2)
        angle = acos(np.dot(n, axis))
        if deg:
            angle = degrees(angle)
        return angle

    def azimuth(self, t, deg=True, std=False, syst=False):
        "Return the azimuth with respect to the positive x-axis."
        x, y, z = self._construct_n(t)
        azimuth = (atan2(z, x) + 2.*pi) % (2.*pi)
        if deg:
            azimuth = degrees(azimuth)
        return self._return_result(azimuth, std, syst)
