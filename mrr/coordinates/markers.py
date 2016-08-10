# -*- coding: utf-8 -*-
"""
This module provides methods to (semi)automatically extract the positions of
markers in a DICOM image as used in MRR.

@author: Sebastian Theilenberg
"""

import numpy as np
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from sklearn.cluster import DBSCAN
from skimage import draw

from .dicom_coordinates import get_matrix, get_position
from ..plotting import display

__version__ = "0.11"


LEFT_CORR = (39.1, 0.0, -29.55)
RIGHT_CORR = (-37.9, 0.0, -29.55)
CRANIAL_CORR = (0.0, 0.0, 26.0)

CORRECTIONS = (LEFT_CORR,
               RIGHT_CORR,
               CRANIAL_CORR)

THREE_MARKERS_BONN = (
    # cranial marker
    {'image': 0,
     'xlimits': (150, 300),
     'ylimits': (100, 250),
     'msize': (500, 1500),
     'bthresh': 3.0},
    # patient right marker
    {'image': 1,
     'xlimits': (250, 350),
     'ylimits': (0, 100),
     'msize': (200, 500),
     'bthresh': 0.8},
    # patient left marker
    {'image': 1,
     'xlimits': (200, 400),
     'ylimits': (400, 512),
     'msize': (200, 500),
     'bthresh': 4.0}
    )

READOUT_SETS = {'Bonn': THREE_MARKERS_BONN}


def read_out_markers(imgset, header, indices=[], zoom=True, average=True):
    """
    indices: [r, c, s, a, b]
    """
    coordinates = []
    M = get_matrix(header)
    print "{:25}{:>9}{:>9}{:>9}".format("", "x", "y", "z")
    for i, index in enumerate(indices):
        try:
            r, c, s, a, b = index
        except ValueError:
            raise ValueError("""
                Wrong number of values, need 3 coordinates and 2 lenghts
                """)

        # calculate coordinates
        pos = get_position(M, r, c, s)
        coordinates.append(pos)

        # display ellipse
        display(imgset[s], hold=True)
        plt.gca().add_artist(Ellipse((c, r), a, b, fill=False, color='red'))
        plt.title("Slice {}".format(s))
        if zoom:
            ymin, ymax = r+b, max(0, r-b)
            xmin, xmax = max(0, c-a), c+a
            plt.ylim(ymin, ymax)
            plt.xlim(xmin, xmax)
            # yaxis is counting fromhigher to lower values, therefore switch
            # ymin and ymax in array slicing
            visible = imgset[s][ymax:ymin, xmin:xmax]
            vmin, vmax = visible.min(), visible.max()
            plt.clim(vmin, vmax)

        pstring = ("Marker {:>3} in slice {:>3}: "
                   "{:>9.3f}{:>9.3f}{:>9.3f}").format(i, s, *pos)
        print pstring

    coordinates = np.asarray(coordinates)
    if average:
        coordinates = coordinates.mean(axis=0)
    return coordinates


def correct_positions(left=None, right=None, cranial=None):
    m_counter = 0
    for marker in [left, right, cranial]:
        if marker is not None:
            if not marker.size == 3:
                raise ValueError(
                    "coordinates should be 3D (got {})".format(marker)
                    )
            marker += CORRECTIONS[m_counter]
        m_counter += 1


def find_marker_regions(img, xlimits=None, ylimits=None,
                        bthresh=1., msize=(150, 1000), dbkwargs=None,
                        verbose=False):
    """
    Find pixel regions in an image belonging to a marker.

    Args:
        img: twodimensional greyscale image data.
        xlimits: tuple of minimal and maximal value of a slice of the image
            to reduce computing time (and memory usage)
        ylimits: see xlimits
        bthresh: float used to create a binary version of img. All pixels are
            set to 1 whose value is greater than img.mean()*bt
        msize: Determines the possible size of a marker region in pixels.
            Expects a tuple (minimum size, maximum size)
        dbkwargs: dictionary containing arguments to sklearn.cluster.DBSCAN
        return_regions: (optional) Whether to return all regions found
        verbose: (optional) Wheter to print out informations and plot results

    Returns:
        A list of regions that are considered a marker based on the provided
        limits. Every region is a list of pixel coordinates in the original
        image.
    """
    dbk = {'eps': 1.5, 'min_samples': 7, 'algorithm': 'auto'}
    if dbkwargs is not None:
        dbk.update(dbkwargs)

    if xlimits is None:
        xlimits = (0, img.shape[1])
    if ylimits is None:
        ylimits = (0, img.shape[0])

    r_corr = np.array([ylimits[0], xlimits[0]])

    # cut image and create binary data
    orig = img[ylimits[0]:ylimits[1], xlimits[0]:xlimits[1]]
    data = orig > orig.mean() * bthresh
    assert np.any(data > 0)

    # fit DBSCAN to binary data to find clusters
    X = np.vstack(np.where(data))
    db = DBSCAN(**dbk).fit(np.asarray(X, dtype=np.float).T)

    # extract regions and find markers candidates
    labels = db.labels_
    u_labels = set(labels)

    regions = {k: labels == k for k in u_labels if k != -1}  # -1: background
    rsizes = {k: r.sum() for k, r in regions.iteritems()}
    markers = [
        k
        for k, s in rsizes.iteritems()
        if (s >= msize[0] and s <= msize[1])
        ]

    if not regions:
        raise ValueError("No regions found! Adjust settings")

    if verbose:
        # create integer image marking different regions
        pltimg = np.zeros_like(orig)
        for k, r in regions.iteritems():
            if k in markers:
                label_ = 2
            else:
                label_ = 1
            ind = X.T[r]
            pltimg[ind[:, 0], ind[:, 1]] = label_
        # mask background regions
        invalid = X.T[labels == -1]
        assert invalid.ndim == 2
        pltimg[invalid[:, 0], invalid[:, 1]] = np.ma.masked

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(orig, cmap='Greys_r', interpolation='none')
        ax.contourf(pltimg, levels=[0.5, 1.5, 2.5],
                    alpha=0.2, colors=['b', 'r'])
        ax.contour(pltimg, levels=(1.95,), colors='r')

        print 'I guess regions {} are markers'.format([
            '{} ({})'.format(k, s)
            for k, s in rsizes.iteritems()
            if k in markers])

    if not markers:
        raise ValueError(
            "No marker region found. Check image limits or markers size")

    # correct marker coordinates to original image data
    marker_regions = [
        X.T[r] + r_corr
        for k, r in regions.iteritems()
        if k in markers
        ]

    return marker_regions


def fit_marker(region, img, verbose=False, log=False):
    """
    Fit an ellipse to a marker region.

    Args:
        region: a list of coordinates specifying the region of the marker.
        img: Twodimensional image data.
        verbose: (optional) Prints out information and plots the result.
        log: (optional) Return the intermediate states of the minimizer for
            debugging purposes.

    Returns:
        The coefficients of the fitted ellpise

            np.array([cy, cx, ry, rx])

        loglist: Only if log == True. Contains for every step of the minimizer
            the current parameters and the residual.
    """

    def cost(params, bimg, loglist=None):
        cy, cx, yr, xr = params
        coords = draw.ellipse(cy, cx, yr, xr, shape=bimg.shape)
        template = np.zeros_like(bimg)
        template[coords] = 1
        res = ((template - bimg)**2).sum()
        if loglist is not None:
            loglist.append((params, res))
        return res

    region = np.asarray(region)
    assert region.shape[1] == 2

    # create binary image
    bimg = np.zeros_like(img)
    bimg[region[:, 0], region[:, 1]] = 1

    initial = (
        region[:, 0].mean(),
        region[:, 1].mean(),
        (region[:, 0].max() - region[:, 0].min())/2.,
        (region[:, 1].max() - region[:, 1].min())/2.
        )

    args = [bimg, ]
    if log:
        loglist = []
        args.append(loglist)
    res = fmin(cost, initial, args=tuple(args), disp=verbose)
    cy, cx, yr, xr = res

    if verbose:
        print 'initial:', initial
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(img, interpolation='none', cmap='Greys_r')
        ax.add_artist(Ellipse((cx, cy), 2.*xr, 2.*yr,
                              color='red', fill=False, linewidth=2))
        ax.plot(cx, cy, 'rx', mew=2)
        ax.set_xlim(cx - 2.*xr, cx + 2.*xr)
        ax.set_ylim(cy - 2.*yr, cy + 2.*yr)

    if log:
        return res, loglist
    else:
        return res
