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
            plt.ylim(r+b, r-b)
            plt.xlim(c-a, c+a)

        pstring = ("Marker {:>3} in slice {:>3}: "
                   "{:>9.3f}{:>9.3f}{:>9.3f}").format(i, s, *pos)
        print pstring

    coordinates = np.asarray(coordinates)
    if average:
        coordinates = coordinates.mean(axis=0)
    return coordinates


left_correction = (39.1, 0.0, -29.55)
right_correction = (-37.9, 0.0, -29.55)
cranial_correction = (0.0, 0.0, 26.0)

corrections = (left_correction,
               right_correction,
               cranial_correction)


def correct_positions(left=None, right=None, cranial=None):
    m_counter = 0
    for marker in [left, right, cranial]:
        if marker is not None:
            if not marker.size == 3:
                raise ValueError(
                    "coordinates should be 3D (got {})".format(marker)
                    )
            marker += corrections[m_counter]
        m_counter += 1
    
