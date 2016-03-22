# -*- coding: utf-8 -*-
"""
Functions to find pixel coordinates from dicom files.

Based on Digital Imaging and Communications in Medicine (DICOM) Part 3,
2011, section C.7.6.2.1.1, p. 409

Created on Tue Mar 03 17:22:56 2015

@author: Sebastian Theilenberg
"""


import numpy as np
from warnings import warn
import dicom

#from ..read.read import read_dicom_headers


#def read_dicom_headers(dcm):
#    """
#    Read headers of a set of dicom files.
#    """
#    files = nameparser(dcm)
#    headers = [dicom.read_file(f_, stop_before_pixels=True) for f_ in files]
#    if len(headers) == 1:
#        return headers[0]
#    else:
#        return headers


def _read_image_plane_module(dcm):
    """
    Return the image plane module of the dicom header.

    Returns:
    ipp : array
        Image Position Patient (x, y, z)
    iop : array
        Image Orientation Patient (x, y, z, x, y, z),
        first three values in r-direction, second three in c-direction
    ps : array
        Pixel Spacing in direction (r, c)
    st : float
        Slice Thickness
    sl : float
        Slice Location
    """
    ipp = np.array([float(i) for i in dcm.ImagePositionPatient])
    iop = np.array([float(i) for i in dcm.ImageOrientationPatient])
    ps = np.array([float(i) for i in dcm.PixelSpacing])
    st = float(dcm.SliceThickness)
    sl = float(dcm.SliceLocation)
    return ipp, iop, ps, st, sl


def _get_matrix_3d(dcms):
    """
    Defines the 3D affine matrix to map voxel coordinates to mm in the
    Dicom Patient Coordinate System.
    """
    t1, iop, ps, st, sl = _read_image_plane_module(dcms[0])
    X, Y = iop[:3], iop[3:]
    dy, dx = ps
    tn = _read_image_plane_module(dcms[-1])[0]
    n = len(dcms)
    k = [(t1[i]-tn[i])/(1-n) for i in range(3)]
    M = np.array([[Y[0]*dy, X[0]*dx, k[0], t1[0]],
                  [Y[1]*dy, X[1]*dx, k[1], t1[1]],
                  [Y[2]*dy, X[2]*dx, k[2], t1[2]],
                  [0, 0, 0, 1]]
                 )
    return M


def _get_matrix_2d(dcm):
    """
    Defines the 2D affine matrix to map pixel coordinates to mm in the
    Dicom Patient Coordinate System.
    """
    ipp, iop, ps, _, __ = _read_image_plane_module(dcm)
    X, Y = iop[:3], iop[3:]
    dy, dx = ps
    M = np.array([[Y[0]*dy, X[0]*dx, 0, ipp[0]],
                  [Y[1]*dy, X[1]*dx, 0, ipp[1]],
                  [Y[2]*dy, X[2]*dx, 0, ipp[2]],
                  [0, 0, 0, 1]]
                 )
    return M


def get_voxel_size(dcm):
    """
    Returns the voxel size of pixels in dcm in mm as (r, c, s)
    """
    ipp, iop, ps, st, sl = _read_image_plane_module(dcm)
    return np.array(ps.tolist() + [st])


def get_matrix(dcm):
    if type(dcm) == dicom.dataset.FileDataset:
        M = _get_matrix_2d(dcm)
    else:
        M = _get_matrix_3d(dcm)
    return M


def get_position(M, r, c, s=0):
    """
    Returns the voxel coordinates (r, c, s) in mm (x, y, z) in the Dicom
    Patient Coordinate System.
    """
    if s != 0 and not np.any(M[:, 2]):
        warn("z location provided, but matrix seems to be in-plane!",
             UserWarning)
    pos = np.dot(M, np.array([r, c, s, 1]).T)
    return pos[:3]
