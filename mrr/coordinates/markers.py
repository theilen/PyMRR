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
from ..plotting.mrrplot import display, create_axes_grid

__version__ = "0.11"


LEFT_CORR = (39.1, 0.0, -29.6)
RIGHT_CORR = (-37.9, 0.0, -29.6)
CRANIAL_CORR = (0.0, 0.0, 26.0)

CORRECTIONS = (LEFT_CORR,
               RIGHT_CORR,
               CRANIAL_CORR)

THREE_MARKERS_BONN = {
    'cranial': {'yarea': (80, 200),
                'xarea': (150, 300)},
    'dorsal': {'yarea': (200, 350)}
    }


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
                        bthresh=40, msize=(200, 1000),
                        verbose=False, plot=False):
    """
    Find pixel regions in an image belonging to a marker.

    Args:
        img: twodimensional greyscale image data.
        xlimits: tuple of minimal and maximal value of a slice of the image
            to reduce computing time (and memory usage)
        ylimits: see xlimits
        bthresh: int used to create a binary version of img. All pixels are
            set to 1 whose value is greater than bthresh
        msize: Determines the possible size of a marker region in pixels.
            Expects a tuple (minimum size, maximum size)
        verbose: (optional) Wheter to print out informations and plot results

    Returns:
        A list of regions that are considered a marker based on the provided
        limits. Every region is a list of pixel coordinates in the original
        image.
    """
    dbk = {'eps': 1.5, 'min_samples': 9, 'algorithm': 'auto'}

    if xlimits is None:
        xlimits = (0, img.shape[1])
    if ylimits is None:
        ylimits = (0, img.shape[0])

    r_corr = np.array([ylimits[0], xlimits[0]])

    # cut image and create binary data
    orig = img[ylimits[0]:ylimits[1], xlimits[0]:xlimits[1]]
    # TODO: Make sure there are no holes in the binary area
    data = orig > bthresh
    assert np.any(data > 0)

    # fit DBSCAN to binary data to find clusters
    X = np.vstack(np.where(data))
    db = DBSCAN(**dbk).fit(np.asarray(X, dtype=np.float).T)

    # extract regions and find markers candidates
    labels = db.labels_
    u_labels = set(labels)

    regions = {k: labels == k for k in u_labels if k != -1}  # -1: background
    if not regions:
        raise ValueError("No regions found! Try adjusting the threshold.")

    rsizes = {k: r.sum() for k, r in regions.iteritems()}
    markers = [
        k
        for k, s in rsizes.iteritems()
        if (s >= msize[0] and s <= msize[1])
        ]

    if plot:
        _plot_markers(img, X.T + r_corr, regions, labels, markers)

    if verbose:
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


def _plot_markers(img, X, regions, labels, markers):
    pltimg = np.zeros_like(img, dtype=np.int)
    for k, r in regions.iteritems():
        if k in markers:
            label_ = 2
        else:
            label_ = 1
        ind = X[r]
        pltimg[ind[:, 0], ind[:, 1]] = label_
    # mask background regions
    invalid = X[labels == -1]
    assert invalid.ndim == 2
    pltimg[invalid[:, 0], invalid[:, 1]] = np.ma.masked

    ax = plt.gca()
    ax.imshow(img, cmap='Greys_r', interpolation='none')
    ax.contourf(pltimg, levels=[0.5, 1.5],
                alpha=0.2, colors=['b'])
    ax.contour(pltimg, levels=(1.95,), colors='r')


def broaden_binary(data, n=1):
    """Grow a binary mask by n pixel."""
    for _ in xrange(n):
        grady, gradx = np.gradient(data)
        data[(np.abs(grady) + np.abs(gradx)) != 0] = True


def fit_marker(region, img, broaden=10, verbose=False, plot=False):
    """
    Fit an ellipse to a marker region.

    Args:
        region: a list of coordinates specifying the region of the marker.
        img: Twodimensional image data.
        verbose: (optional) Prints out information.
        plot: (optional) Plot the result.

    Returns:
        The coefficients of the fitted ellpise

            np.array([cy, cx, ry, rx])

    """
    try:
        data = img.phase
    except AttributeError:
        data = np.array(img)
    assert data.ndim == 2

    region = np.asarray(region)
    assert region.shape[1] == 2

    fitmap = _create_fitmap(data, region, broaden, plot)

    initial = (
        region[:, 0].mean(),
        region[:, 1].mean(),
        (region[:, 0].max() - region[:, 0].min())/2.,
        (region[:, 1].max() - region[:, 1].min())/2.
        )

    res = _fit_marker_ellipse(fitmap, initial)

    if verbose:
        print 'initial values:', initial

    if plot:
        ax = plt.gca()
        ax.imshow(data, interpolation='none', cmap='Greys_r', vmax=0.1)
        cy, cx, ry, rx = res
        ax.add_artist(Ellipse((cx, cy), 2.*rx, 2.*ry,
                              color='red', fill=False, linewidth=2))
        ax.plot(cx, cy, 'rx', mew=2)
        ax.set_xlim(cx - 3.*rx, cx + 3.*rx)
        ax.set_ylim(cy - 3.*ry, cy + 3.*ry)

    return res


def _create_fitmap(data, region, broaden=0, plot=False):
    """Create an image based on data with void data outside of region."""
    ry, rx = region.T
    bimg = np.zeros_like(data, dtype=np.bool)
    bimg[ry, rx] = True
    broaden_binary(bimg, n=broaden)

    if plot:
        ax = plt.gca()
        ax.contour(bimg, levels=(0.95,), colors='y')

    fitmap = np.array(data, dtype=np.float, copy=True)
    fitmap[~bimg] = 0.
    norm = data[ry, rx].max()
    fitmap /= norm

    return fitmap


def _fit_marker_ellipse(fitmap, initial):
    """Fit an ellipse to a fitmap."""

    def cost(params, bimg):
        cy, cx, yr, xr = params
        coords = draw.ellipse(cy, cx, yr, xr, shape=bimg.shape)
        template = np.zeros_like(bimg)
        template[coords] = 1
        res = ((template - bimg)**2).sum()
        return res

    return fmin(cost, initial, args=(fitmap,), disp=False)


def analyze_image_set(images, bthresh=40, yarea=None, xarea=None, plot=False):
    """
    Find markers in a set of images.

    Args:
        images: 3D set of 2D images
        bthresh: Threshold used to find regions in the images. Absolute value.
        xarea, yarea: min, max values for markers. regions outside these values
            will be ignored.
        plot: (optional) Plots the analyzed images
    """
    try:
        data = images.phase
    except AttributeError:
        data = np.array(images)
    else:
        bthresh /= 4096.

    regions = _find_regions_set(data, bthresh, yarea, xarea, plot)
    markers = _fit_markers_set(data, regions, plot)
    if hasattr(images, 'matrix'):
        M = images.matrix
        markers = [get_position(M, s=sl, *m[0:2]) for (sl, m) in markers]
        markers.sort(key=lambda x: x[0])
        markers = np.array(markers)
    return markers


def _find_regions_set(images, bthresh, yarea=None, xarea=None, plot=False):
    """Find regions in a set of images."""
    # test for phase data
    n = images.shape[0]

    if plot:
        _, axes = create_axes_grid(n)

    regions = []
    for i in range(n):
        if plot:
            plt.sca(axes[i])
        try:
            regs_ = find_marker_regions(images[i], bthresh=bthresh,
                                        verbose=False, plot=plot)
        except ValueError:
            pass
        else:
            regions.extend([(i, r) for r in regs_])
        finally:
            if plot:
                axes[i].set_title('Image {}'.format(i))
    if yarea:
        regions = [
            (i, r)
            for (i, r) in regions
            if np.all(r.T[0] > yarea[0]) and np.all(r.T[0] < yarea[1])
            ]
        if plot:
            for axis in axes:
                for y_ in yarea:
                    axis.axhline(y_, color='g')
    if xarea:
        regions = [
            (i, r)
            for (i, r) in regions
            if np.all(r.T[1] > xarea[0]) and np.all(r.T[1] < xarea[1])
            ]
        if plot:
            for axis in axes:
                for x_ in xarea:
                    axis.axvline(x_, color='g')
    regions.sort(key=lambda r: (r[1][:, 1].min(), r[0]))
    return regions


def _fit_markers_set(images, regions, plot=False):
    """Fit markers to regions in a set of images."""
    markers = []
    if plot:
        _, axes = create_axes_grid(len(regions))
    for i, (sl, reg) in enumerate(regions):
        if plot:
            plt.sca(axes[i])
            axes[i].set_title("Image {}".format(sl))
        try:
            marker = fit_marker(reg, images[sl], plot=plot)
        except Exception:
            pass
        else:
            markers.append((sl, np.array(marker)))
    return markers
