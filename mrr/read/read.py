# -*- coding: utf-8 -*-
"""
last change: Wed Oct 22 16:17 2014

@author: Sebastian Theilenberg
"""

import dicom
import numpy as np
import os
import re
from PIL import Image
from warnings import warn

from ..mrrcore import MRRArray, cond_print, empty, copy_attributes
from ..unwrapping import unwrap_array
from .parse_dicom import parse_parameters, check_sequence
from ..coordinates.dicom_coordinates import get_matrix

__version__ = '1.32'

# $Source$


# version history
# =================
# version 1.3
# -----------
# - read_dicom:
#        - added parsing of parameters
# - read_dicom_set:
#       - added parsing of parameters
#
# version 1.2
# -----------
# - read_dicom:
#     - added test for unwrapper
# - read_data_set:
#   - added support for unwrapper
#   - added test for unwrapper
#   - changed exception handling
#   - edited docstring
# - removed *import from core
#
# version 1.1.1
# -------------
# - switched to relative imports
# - added support for different unwrap-algorithms to read_dicom and
#   read-dicom_set to support unwrapping version 1.2
#
# version 1.1
# -----------
# - read_dicom_set: added IOError if nameparser returns empty list
# - read_mask: added a test for invertet masks
# - nameparser: added ValueError if filename does not match pattern
#   (len(items)<2)


__metaclass__ = type


class MissingFileError(Exception):
    '''Raised, whenever a file to be loaded was not found.'''
    def __init__(self, *args):
        if not args:
            args = ('File not found',)
        Exception.__init__(self, *args)


class UnwrapperError(Exception):
    '''Raised, whenever unwrapping fails.'''
    def __init__(self, *args):
        if not args:
            args = ('Could not unwrap!',)
        Exception.__init__(self, *args)


def nameparser(filename):
    '''
    Parses a given directory and returns all files matching filename_<i>.
    Returns a sorted list of filenames.

    Tested with:
    '188_13-12-10_82'
    '188_13-12-10_82_1'
    '13-12-10_82_1'
    '''
    path = os.path.abspath(filename)
    directory, name = os.path.split(path)
    # create searchpattern
    items = name.split('_')
    # find image-number
    if len(items) < 2:
        raise ValueError(
            "Given filename '{}' does not match parsing pattern!".format(
                filename)
            )
    for i in range(2):
        if re.match(r'\d\d-\d\d-\d\d', items[i]):
            index = i
            break
    try:
        items[index+2] = '(\d+)'
    except IndexError:
        items.append('(\d+)')
    search_pattern = re.compile('_'.join(items)+'$')
    # find files
    files = [f for f in os.listdir(directory) if re.match(search_pattern, f)]
    # sort numerically
    files.sort(key=lambda f: int(re.match(search_pattern, f).groups()[0]))
    # return list of absolute paths
    return [os.path.join(directory, f) for f in files]


def _get_pixel_data(dc):
    """Return the processed pixel data of a dicom object"""
    return np.asarray(dc.pixel_array, dtype=np.float32)/4096.


def _unwrap_pixel_data(data, mask, algorithm, filename=None, **kwargs):
    pixel_data, add = unwrap_array(data, mask, additional=True,
                                   algorithm=algorithm, **kwargs)
    if algorithm in ['py_gold', ]:
        if add[-1] != 1:
            warn(("Found disconnected pieces while unwrapping"
                  "{}").format(filename))
    return pixel_data, add


def _mrr_from_dicom(dc, unwrap, mask, unwrapper, **ukwargs):
    pixel_data = _get_pixel_data(dc)
    if unwrap:
        pixel_data, add = _unwrap_pixel_data(pixel_data, mask, unwrapper,
                                             *ukwargs)
    else:
        add = None

    seq_data = parse_parameters(dc)
    seq_data.update({'unwrapped': unwrap})
    data = MRRArray(pixel_data, mask=mask, **seq_data)

    return data, add


def read_dicom(dicom_file, unwrap=False, mask=None, verbose=True,
               unwrap_data=False, unwrapper='py_gold', **ukwargs):
    '''
    Reads in one dicom-file and returns the image data as 2d-MRRArray.
    If unwrap is set (default: False) the data is automatically unwrapped. If
    so, an array <mask> containing the mask has to be provided as well!
    '''
    dc = dicom.read_file(dicom_file)
    cond_print('Read in file %s' % dicom_file, verbose)
    data, add = _mrr_from_dicom(dc, unwrap, mask, unwrapper, **ukwargs)

    data.orig_file = os.path.basename(dicom_file)

    if unwrap_data:
        return data, add
    else:
        return data


def read_dicom_set(dicom_file, unwrap=False, mask=None, verbose=False,
                   unwrapper='py_gold'):
    '''
    Reads-in the whole set of dicom-files belonging to that series and returns
    it as a multidimensional array with the first dimension being the image
    number.

    Parameters
    ----------
    dicom_file : str
        Path of one arbitrary file of the set to be read-in.
    unwrap : bool (optional)
        whether to unwrap the data while reading it. (Default: False)
    mask : 2darray (optional)
        mask to set in the array. Mandatory if unwrap==True!
    verbose : bool (optional)
        writes information in stdout. (Default: False)
    unwrapper : str
        which unwrapper to use if unwrap==True. (Default: 'py_gold')

    Returns
    -------
    data : 3darray | list of 3darrays
        Read-in data. A MRRArray if only one PTFT was present, a list of
        MRRArrays otherwise, sorted by PTFT. For every MRRArray, the first
        axis corresponds to the image number.
    '''
    files = nameparser(dicom_file)
    if len(files) == 0:
        raise IOError("Did not find any dicom files! Wrong path?")
    nofimages = len(files)
    basename = "_".join(dicom_file.split('_')[:-1])
    cond_print('Found {} file(s) in total'.format(nofimages), verbose)

    # Read files
    dicom_data = [dicom.read_file(f) for f in files]
    # check for nin_ep2d_diff first
    is_epi = check_sequence(dicom_data[0])
    result = []
    images = []
    index = 0
    while index < nofimages:
        dc, _ = _mrr_from_dicom(dicom_data[index], unwrap, mask, unwrapper)

        if is_epi:
            # sort and divide epi files per PTFT
            try:
                ptft = images[-1].PTFT
            except IndexError:
                ptft = dc.PTFT

            if dc.PTFT == ptft:  # Collect data with same PTFT
                images.append(dc)
            else:
                cond_print("new PTFT: {}".format(dc.PTFT), verbose)
                # Write all data in images into one MRRArray
                result.append(_collect_data(images, basename))
                # restart images with new PTFT
                images = []
                images.append(dc)
        else:
            images.append(dc)
        # Increase index before next file
        index += 1
    # collect remaining data
    if images:
        result.append(_collect_data(images))

    # update matrix to 3D
    if nofimages > 1:
        M = get_matrix(dicom_data)
        for d in result:
            d.matrix = M.copy()

    if len(result) == 1:
        result = result[0]
    elif is_epi:
        result.sort(key=lambda x: x.PTFT)

    return result


def _collect_data(images, orig_file=""):
    shape = images[0].shape
    data = empty((len(images), shape[0], shape[1]))
    for i, item in enumerate(images):
        data[i] = item
    copy_attributes(data, images[0])
    data.orig_file = orig_file
    return data


def read_mask(filename, check_invert=True):
    '''
    Reads in a raster graphics file and converts it to a boolean mask.

    Valid pixels are set to True, invalid pixels are set to False. If the mask
    seems to be inverted (e.g. outer parts are valid and inner parts are
    invalid), the mask is corrected. This behavior may be switched off using
    check_invert=False.
    '''
    mask = np.asarray(Image.open(filename), dtype=np.bool)
    if check_invert:
        if mask[0, 0]:
            return np.invert(mask)
    return mask


def read_bitmap(bitmap_file):
    '''
    Reads-In a bitmap file and returns the image-data as numpy-array.
    Do NOT use to image MRI-data!
    '''
    return np.asarray(Image.open(bitmap_file))


def read_dicom_headers(dcm):
    """
    Read headers of a set of dicom files.
    """
    files = nameparser(dcm)
    headers = [dicom.read_file(f_, stop_before_pixels=True) for f_ in files]
    if len(headers) == 1:
        return headers[0]
    else:
        return headers
