# -*- coding: utf-8 -*-
"""
Created on Tue May 19 14:34:50 2015

@author: Sebastian Theilenberg
"""

__version__ = "0.1"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

from .dicom_coordinates import get_matrix, get_position
from ..plotting import display


def read_out_markers(imgset, header, indices=[], zoom=True, average=True):
    """
    
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

        print "Marker {:>3} in slice {:>3}: {:>9.3f}{:>9.3f}{:>9.3f}".format(i, s, *pos)

    coordinates = np.asarray(coordinates)
    if average:
        coordinates = coordinates.mean(axis=0)
    return coordinates
