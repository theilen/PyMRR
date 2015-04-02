# -*- coding: utf-8 -*-
"""
last change: Wed Oct 22 16:17 2014

@author: Sebastian Theilenberg
"""

__version__ = '1.31'
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

import dicom
import numpy as np
import os
import re
from PIL import Image


from ..mrrcore import MRRArray, cond_print, empty
from ..unwrapping import unwrap_array, valid_unwrapper
from .parse_dicom import read_parameters, parse_parameters


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

    Getestet mit:
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


def read_dicom(dicom_file, unwrap=False, mask=None, verbose=True,
               unwrap_data=False, unwrapper='py_gold', **ukwargs):
    '''
    Reads in one dicom-file and returns the image data as 2d-MRRArray.
    If unwrap is set (default: False) the data is automatically unwrapped. If
    so, an array <mask> containing the mask has to be provided as well!
    '''
    if unwrap is True and not valid_unwrapper(unwrapper):
        raise AttributeError('No algorithm named {}'.format(unwrapper))

    dc = dicom.read_file(dicom_file)
    cond_print('Read in file %s' % dicom_file, verbose)
    pixel_data = np.asarray(dc.pixel_array, dtype=np.float32)/4096.
    seq_data = parse_parameters(dc)

    # Unwrap if specified
    if unwrap:
        if not np.any(mask):
            raise UnwrapperError('Could not unwrap file %s. \
                                  No mask provided!' % dicom_file)
        # try to unwrap
        try:
            pixel_data, add = unwrap_array(pixel_data, mask, additional=True,
                                           algorithm=unwrapper, **ukwargs)
        except:
            raise
            # raise UnwrapperError('Could not unwrap file %s' % dicom_file)
        if unwrapper in ['py_gold', 'c_gold']:
            if add[-1] != 1:
                print ("WARNING: found disconnected "
                       "pieces while unwrapping {}").format(dicom_file)

    # create MRRArray
    seq_data.update({"orig_file": os.path.basename(dicom_file),
                     "unwrapped": unwrap})
    data = MRRArray(pixel_data, mask=mask, **seq_data)

    # return data
    if unwrap_data:
        return data, add
    return data


def read_dicom_set(dicom_file, unwrap=False, mask=None, verbose=True,
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
        wether to unwrap the data while reading it. (Default: False)
    mask : 2darray (optional)
        mask to set in the array. Mandatory if unwrap==True!
    verbose : bool (optional)
        writes information in stdout. (Default: True)
    unwrapper : str
        which unwrapper to use if unwrap==True. (Default: 'py_gold')

    Returns
    -------
    data : 3darray
        MRRArray containing the read-in data. The first axis corresponds to the
        image number.
    '''
    if unwrap is True and not valid_unwrapper(unwrapper):
        raise AttributeError('No algorithm named {}'.format(unwrapper))

    # read-in files
    files = nameparser(dicom_file)
    if len(files) == 0:
        raise IOError("Did not find any dicom files! Wrong path?")
    # read first file
    nofimages = len(files)
    cond_print('Found {} file(s) in total'.format(nofimages), verbose)
    dc = read_dicom(files[0], unwrap, mask, verbose, unwrapper=unwrapper)

    # write MRRArray
    seq_data = read_parameters(files[0])
    seq_data.update({"orig_file": os.path.basename(dicom_file),
                     "unwrapped": unwrap})
    shape = dc.shape
    data = empty((nofimages, shape[0], shape[1]), **seq_data)
    data[0] = dc
    # fill remaining array
    i = j = 1
    while i in xrange(1, nofimages):
        try:
            data[j] = read_dicom(files[i], unwrap, mask, verbose,
                                 unwrapper=unwrapper)
        except:
            cond_print('File "%s" could not be read...' % files[i], verbose)
            raise
        else:
            j += 1
        i += 1
    # Clean size of array
    if j != i:
        data = np.delete(data, np.s_[j:], axis=0)
        cond_print('WARNING: %i file(s) could not be read in! ' % (i-j) +
                   'Read %i files in total.' % j)

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
        if mask[0, 0] is True:
            return np.invert(mask)
    return mask


def read_bitmap(bitmap_file):
    '''
    Reads-In a bitmap file and returns the image-data as numpy-array.
    Do NOT use to image MRI-data!
    '''
    return np.asarray(Image.open(bitmap_file))
