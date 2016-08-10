"""
Provides front-end function unwrap_array to use the unwrapping-algorithms
present in unwrapping.

@author: Sebastian Theilenberg
"""

import numpy as np

from .py_unwrap import unwrap_py_gold, unwrap_quality_guided, \
    unwrap_mask_cut, DerivativeVariance
from ..mrrcore import MRRArray


__version__ = '1.4.1'

# $Source$


# version history
# ===============
# version 1.4.2
# - added unwrap_image as wrapper of unwrap_array to MRRArray

# version 1.4.1
# - added handling of zero-values in qualitymap
# - removed unwrapper c_gold
#
# version 1.4
# - added __all__
# - unwrap_array:
#   - added mask_cut
#   - added test for MRRArray
# - added valid_unwrapper
#
# version 1.3
# - added quality_guided to unwrap_array
#
# version 1.2
# - added py_gold to unwrap_array, major changes in implementation


__all__ = ['unwrap_array', 'unwrap_image', 'valid_unwrapper']


def unwrap_array(array, mask, algorithm='py_gold', additional=False,
                 **keyargs):
    """
    Unwrap two-dimensional using various algorithms.

    Parameters
    ----------
    array : 2darray
        two-dimensional numpy array containing the wrapped data. Note that the
        data should be normalized, i.e. the values should be in the range
        of 0.0 to 1.0. If the data is not provided in this way, an attemp to
        normalize it is made.
    mask : 2darray
        two-dimensional numpy array containing the mask to be used during
        unwrapping. Make sure the masks dimensions match the ones of
        wrapped_array!
    algorithm : string
        One of the following:
            'c_gold' : Wrapper for Ghiglia's C-implementation of Goldstein's
                        algorithm.
            'quality_guided' : pure-python implementation of a quality guided
                                phase-unwrapping algorithm using the inverse
                                phase-derivative as a quality-map.
                                see Ghiglia, Pritt "Two-Dimensional Phase
                                Unwrapping Theory, Algorithms, and Software"
                                John Wiley & Sons 1998, p. 122-126
            'mask_cut' : pure-python implementation of a quality-guided
                            mask-cut-algorithm, using the inverse
                            phase-derivative as quality-map.
                            See Ghiglia, Pritt "Two-Dimensional Phase
                            Unwrapping Theory, Algorithms, and Software"
                            John Wiley & Sons 1998, p. 137-141
            'py_gold' : pure-python implementation of Goldstein's algorithm.
                        (Default)
    additional : bool (optional)
        Whether to return the additional data created during unwrapping.

    Returns
    -------
    unwrapped : 2darray
        The unwrapped array
    data : tuple
        Tuple containing additional data used or created during the unwrapping
        process. The content depends on the unwrapper used:
            quality_guided: quality_map
            mask_cut:       (quality_map, flag_array, num_pieces)
            py_gold:        (flag_array, num_pieces)
    """
    if algorithm not in ['py_gold', 'quality_guided', 'mask_cut']:
        raise AttributeError('No algorithm named {}'.format(algorithm))

    if not np.any(mask):
        raise ValueError("No mask provided!")

    if array.shape != mask.shape:
        raise ValueError("Mask does not fit to array!")

    if isinstance(array, MRRArray):
        array = array['phase'].view(np.ndarray)

    # test if array wrapped to 0..1
    if array.min() < 0.0 or array.max() > 1.0:
        print '''Warning: Array should be normalized to contain data ranging
        from 0.0 to 1.0. Data renormalized, check results!'''
        normalize_array(array)

#    if algorithm == 'c_gold':
#        if not goldstein_available:
#            print 'c_gold is not available, using py_gold instead'
#            algorithm = 'py_gold'
#        else:
#            unwr, data = c_gold(array.astype(np.float32), mask, **keyargs)

    if algorithm == 'py_gold':
        unwr, flags, num = unwrap_py_gold(array, mask, **keyargs)
        data = (flags, num)

    if algorithm == 'quality_guided':
        qual = 1./(DerivativeVariance(array) + 1e-3)  # avoid division by zero
        unwr = unwrap_quality_guided(array, qual, mask, **keyargs)
        data = qual

    if algorithm == 'mask_cut':
        qual = 1./(DerivativeVariance(array) + 1e-3)  # avoid division by zero
        unwr, flags, num = unwrap_mask_cut(array, qual, mask=mask, **keyargs)
        data = (qual, flags, num)

    if additional:
        return unwr, data
    else:
        return unwr


def normalize_array(array):
    '''Wraps an array to 0..1'''
    if array.min() < 0.0:
        array -= array.min()
    if array.max() > 1.0:
        array /= array.max()


def valid_unwrapper(unwrapper):
    '''Returns True if unwrapper can be used with unwrap_array.'''
    alg_ = ['py_gold', 'quality_guided', 'mask_cut']
    if unwrapper not in alg_:
        return False
    else:
        return True


def unwrap_image(image, algorithm='py_gold', additional=False, **keyargs):

    image['phase'], data = unwrap_array(
        image.phase, image.mask, algorithm=algorithm, additional=True,
        **keyargs)
    if additional:
        return data
