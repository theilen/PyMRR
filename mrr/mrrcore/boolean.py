import numpy as np

__version__ = '1.0'
# $Source$

'''Provides boolean operations for MRRArray.'''

#Tested
def all_equal(a, b, field='phase', axis=None):
    '''
    Test wether all elements in a specific field of a and b are equal along a
    given axis.

    Evaluates np.all(a[field]==b[field], axis)
    If field==None the above is evaluated over all fields!
    
    See np.all for more detailled information.
    '''
    if field==None:
        fields = a.dtype.names
        result = np.all(a[fields[0]]==b[fields[0]])
        for f in fields:
            result = result and np.all(a[f]==b[f])
    else:
        result = np.all(a[field]==b[field], axis=axis)
    try:
        return result.item()
    except ValueError:
        return result.view(np.ndarray)


# Tested
def any_equal(a, b, field='phase', axis=None):
    '''
    Test wether at least one element in a specific field along a given axis is
    equal.

    Evaluates np.any(a[field]==b[field], axis)
    If field==None the above is evaluated over all fields!

    See np.any for more detailled information.
    '''
    if field is None:
        fields = a.dtype.names
        result = np.any(a[fields[0]] == b[fields[0]])
        for f in fields:
            result = result or np.all(a[f] == b[f])
    else:
        result = np.any(a[field] == b[field], axis=axis)
    try:
        return result.item()
    except ValueError:
        return result.view(np.ndarray)


def _shrink_mask_1D(mask, n):
    mask = mask.copy()
    for n_ in range(n):
        indices = np.where(mask == True)[0]
        try:
            indices = indices[[0, -1]]
        except IndexError:  # indices either empty or len == 1
            pass
        mask[indices] = False
    return mask


def shrink_mask(mask, axis, n=1):
    """
    Shrink a boolean mask along axis by n pixels.
    """
    return np.apply_along_axis(_shrink_mask_1D, axis, mask, n=n)
