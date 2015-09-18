'''Provides element-wise arithmetics for MRRArray.'''

#Version history
#==================
#
#changes in 1.1
#----------------
#- reimplemented add, subtract, multiply, divide
#   - useable with skalars, arrays and MRRArrays
#- remapped sub, mult and div for backwards compatibility
#- added _combine_mask()
#- added negate()
#- added absolute()

__version__ = '1.1'
# $Source$

import numpy as np

#from .core import MRRArray


def _combine_mask(a, b):
    """Returns the combined mask of a and b."""
    return (a['mask'] & b['mask']).view(np.ndarray)


def negate(a): #tested
    "Negate a, and return a new MRRArray."
    res = a.copy()
    res['phase'] = res.phase * -1.0
    return res


def absolute(a): #tested
    res = a.copy()
    res['phase'] = np.abs(res.phase)
    return res


def add(a, b, out=None): #tested
    '''
    Add a and b elementwise with respect to MRRArrays.

    'phase' is calculated as sum. 'dev' is computed with gaussian propagation
    of uncertainty, assuming only statistical uncertainties present and
    unequal input-arrays. 'mask' is set according to a logical AND.

    Parameters
    ----------
    a : MRRArray
    b : array-like or skalar
    out : MRRArray (optional)
        if provided, the result will be placed in out. Out must have the same
        shape as the expected result.

    Returns
    -------
    out : MRRAarray
    '''
    if isinstance(b, basestring):
        return NotImplemented
    if out is None:
        out = a.copy()

    try:
        b_ = b.phase
    except AttributeError:
        # b is scalar
        b_ = np.asarray(b)
    else:
        #b is MRRArray
        #set standard deviation
        d_ = b.dev
        np.power(out.dev, 2., out=out['dev'])
        np.add(np.power(d_, 2.), out.dev, out=out['dev'])
        np.sqrt(out.dev, out=out['dev'])
        #set mask
        out['mask'] = _combine_mask(a, b)
    #set phase for both cases
    finally:
        not_implemented = np.add(b_, out.phase, out=out['phase'])

    if np.any(not_implemented) == NotImplemented:
        return NotImplemented
    else:
        return out


def subtract(a, b, out=None): #tested
    '''
    Subtract b from a elementwise with respect to MRRArrays.

    'phase' is calculated as sum. 'dev' is computed with gaussian propagation
    of uncertainty, assuming only statistical uncertainties present and
    unequal input-arrays. 'mask' is set according to a logical AND.

    Parameters
    ----------
    a : MRRArray
    b : array or skalar
    out : MRRArray (optional)
        if provided, the result will be placed in out. Out must have the same
        shape as the expected result.

    Returns
    -------
    out : MRRAarray
    '''
    if isinstance(b, basestring):
        return NotImplemented
    return add(a, -b, out)
#backward-compability
sub = subtract


def multiply(a, b, out=None): #tested
    '''
    Multiply a and b elementwise with respect to MRRArrays.

    'phase' is calculated as a*b. 'dev' is computed with gaussian propagation
    of uncertainty, assuming only statistical uncertainties present and
    unequal input-arrays. 'mask' is set according to a logical AND.

    Parameters
    ----------
    a : MRRArray
    b : array or skalar
    out : MRRArray (optional)
        if provided, the result will be placed in out. Out must have the same
        shape as the expected result.

    Returns
    -------
    out : MRRAarray
    '''
    if isinstance(b, basestring):
        return NotImplemented
    if out is None:
        out = a.copy()

    pa = a.phase
    da = a.dev
    try:
        pb = b.phase
    except AttributeError:
        # b is scalar
        pb = np.asarray(b)
        # standard deviation
        np.multiply(out.dev, pb, out=out['dev'])
    else:
        # b is MRRArray
        db = b.dev
        # standard deviation
        np.sqrt(np.power(pb*da, 2.) + np.power(pa*db, 2.), out=out['dev'])
        # new mask
        out['mask'] = _combine_mask(a, b)
    finally:
        # multiplication of the phase
        ni = np.multiply(out.phase, pb, out=out['phase'])

    if np.any(ni) == NotImplemented:
        return NotImplemented
    else:
        return out
#backward compatibility
mult = multiply


def divide(a, b, out=None): #tested
    '''
    Divide a by b elementwise with respect to MRRArrays.

    'phase' is calculated as a/b. 'dev' is computed with gaussian propagation
    of uncertainty, assuming only statistical uncertainties present and
    unequal input-arrays. 'mask' is set according to a logical AND.

    Parameters
    ----------
    a : MRRArray
    b : array or skalar
    out : MRRArray (optional)
        if provided, the result will be placed in out. Out must have the same
        shape as the expected result.

    Returns
    -------
    out : MRRAarray
    '''
    if isinstance(b, basestring):
        return NotImplemented
    if out is None:
        out = a.copy()

    pa = a.phase
    da = a.dev
    try:
        pb = b.phase
    except AttributeError:
        # b is scalar
        pb = np.asarray(b)
        # standard deviation
        np.true_divide(out.dev, pb, out=out['dev'])
    else:
        # b is MRRArray
        db = b.dev
        # standard deviation
        np.sqrt(np.power(da/pb, 2.) + np.power(pa*db, 2)/np.power(pb, 4.),
                out=out['dev'])
        # new mask
        out['mask'] = _combine_mask(a, b)
    finally:
        # divide phase
        ni = np.true_divide(out.phase, pb, out=out['phase'])

    if np.any(ni) == NotImplemented:
        return NotImplemented
    else: return out
#backward compatibility
div = divide


def power(a, b, out=None): #tested
    '''
    Compute a raised to the power of b elementwise with respect to MRRArrays.
    b has to be a skalar value!

    'phase' is calculated as a**b. 'dev' is computed with gaussian propagation
    of uncertainty, assuming only statistical uncertainties present and
    unequal input-arrays. 'mask' remains unaltered.

    Parameters
    ----------
    a : MRRArray
    b : skalar
    out : MRRArray (optional)
        if provided, the result will be placed in out. Out must have the same
        shape as the expected result.

    Returns
    -------
    out : MRRAarray
    '''
    if hasattr(b, '__len__'):
        return NotImplemented

    if out is None:
        out = a.copy()

    np.power(a.phase, b, out=out['phase'])
    out['dev'] = b * a.dev * np.power(a.phase, b-1.)

    return out
