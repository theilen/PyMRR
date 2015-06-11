# -*- coding: utf-8 -*-
"""
Created on Fri Mar 06 20:27:30 2015

@author: Sebastian Theilenberg
"""

import numpy as np
from scipy.interpolate import interp1d
import math


class moving_plane(object):

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
        "Construct a normalized normal vector to p1-p3 and p2-p3 at time t."
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
        try:
            len(t)
        except TypeError:
            t = [t]
        finally:
            t = np.array(t)

        result = np.empty(t.shape)
        for i in range(t.size):
            coords_ = coords[:]
            t_ = t[i]
            n = self._construct_n(t_)
            res = np.dot(self._basepoint(t_), n)  # constant of normal form
            for c in range(3):
                if c == index:
                    continue
                else:
                    res -= n[c] * coords_.pop(0)
            result[i] = res/n[index]
        return result

    def x(self, t, y, z):
        "Returns the x-component of a point in the yz-plane at time t."
        return self._calc_param(t, 0, [y, z])

    def y(self, t, x, z):
        "Returns the y-component of a point in the xz-plane at time t."
        return self._calc_param(t, 1, [x, z])

    def z(self, t, x, y):
        "Returns the z-component of a point in the xy-plane at time t."
        return self._calc_param(t, 2, [x, y])

    def inclination(self, t):
        "Return the inclination of n with respect to the negative y-axis."
        return self.orientation(t, axis=np.array([0, -1 ,0]))

    def orientation(self, t, axis, deg=True):
        "Return the angle of the nromal vector n to axis at time t."
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
