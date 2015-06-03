# -*- coding: utf-8 -*-
"""
Created on Fri Mar 06 20:27:30 2015

@author: theilenberg
"""

import numpy as np
from scipy.interpolate import interp1d
import math


class moving_plane():

    def __init__(self, p1, p2, p3):
        """
        Create a moving plane based on the trajectories of three basepoints.

        p1, p2, p3 : array_like
            trajectories of the three basepoints (t, x, y, z)
        """
        self._points = []
        for p in [p1, p2, p3]:
            if not len(p) == 4:
                raise ValueError("expected size 4 array")
            interp = interp1d(p[0], p[1:], copy=True,
                              bounds_error=True, assume_sorted=True)
            self._points.append(interp)
        self._basepoint = self._points.pop()

    def __call__(self, t):
        'Returns the planes coefficients (n1, n2, n3, d) at time t.'
        return self._create_plane(t)

    def _get_points(self, t):
        points = []
        for i in range(2):
            points.append(self._points[i](t))
        return points

    def _create_plane(self, t):
        "Creates the plane for a single time t."
        n = self._construct_n(t)
        basepoint =self._basepoint(t)
        d = np.dot(n, basepoint)
        plane = np.array([n.tolist()] + [d])
        return plane

    def _construct_n(self, t):
        """Construct a normalized normal vector to p1-p3 and p2-p3."""
        p1, p2 = self._get_points(t)
        p3 = self._basepoint(t)
        n = np.cross(p1-p3, p2-p3)
        norm = np.linalg.norm(n, ord=2)
        return n/norm
        
    def _calc_param(self, t, index, coords):
        """
        Return the coordinate index at time t with given two coordinates 
        coords. coords in the form [x, y, z] with only two of them given.
        """
        n = self._construct_n(t)
        res = np.dot(self._basepoint, n)  # constant of normal form
        for i in range(3):
            if i == index:
                continue
            res -= n[i] * coords.pop(0)
        return res/n[index]

    @property
    def x(self, t, y, z):
        return self._calc_param(t, 0, [y, z])

    @property
    def y(self, t, x, z):
        return self._calc_param(t, 1, [x, z])

    @property
    def z(self, t, x, y):
        return self._calc_param(t, 2, [x, y])

#    def _calc_orientation(self, t, vect):
#        n = self._construct_n(t)
#        x, y, z = n
#        incl = math.acos(z)
#        azimuth = (math.atan2(y, x) + 2.*math.pi) % (2.*math.pi)
#        return incl, azimuth

    def inclination(self, t):
        "Return the inclination of n with respect to the negative y-axis."
        return self.orientation(t, axis=np.array([0, -1 ,0]))

    def orientation(self, t, axis, deg=True):
        n = self._construct_n(t)
        axis = axis/np.linalg.norm(axis, ord=2)
        angle = math.acos(np.dot(n, axis))
        if deg:
            angle = math.degrees(angle)
        return angle

    def azimuth(self, t, deg=True):
        "Return the azimuth with respect to the negative y-axis."
        x, y, z = self._construct_n(t)
        azimuth = (math.atan2(z, x) + 2.*math.pi) % (2.*math.pi)
        return math.degrees(azimuth)

    def get_orientation(self, t):
        """
        Calculates the orientation of the normal vector n at time t.

        Returns
        -------
        inclination : float
            Angle (in radian) between the positive z-axis and n.
        azimuth : float
            Angle (in radian) between the positive x-axis and the projection
            of n onto the xy-plane. Measured from 0 to 2pi.
        """
        if hasattr(t, "__len__"):
            result = np.empty((len(t), 2))
            for i, t_ in enumerate(t):
                result[i] = self._calc_orientation(t_)
        else:
            result = self._calc_orientation(t)
        return np.array(result)
