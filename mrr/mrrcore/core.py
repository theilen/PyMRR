# -*- coding: cp1252 -*-


# Version history
# ==================
# changes in 1.2.2
# - added additional attributes to MRRArray
#
# changes in 1.2.1
# - reimplemented load to return MRRArray
# - changes in mrr_mean to handle weights in case of missing deviation-data
# - bugfix in mrr_mean to handle broadcasting in case of axis!=0
# - bugfix in mrr_mean to handle zeros in variance
#
# changes in 1.2
# ----------------
# - major changes in MRRArray:
#   - added numeric magic functions to MRRArray
#      - __add__, __sub__. __mult__, __div__ and the respective i- and r-forms
#      - __pow__, __neg__, __abs__
#   - added min, max and mean as methods to MRRArray
#   - added method copy() to MRRArray
# - reimplemented statistics.mean and added as mrr_mean()
#
# changes in 1.1
# ----------------
# - added mrr_min() and mrr_max()
# - added zeros() and zeros_like()
# - added load() and save()
# - minor change to MRRArray:
#       - changed default-value of mask to None
#       - set mask to True if mask==None
#

__version__ = '1.2.3'
# $Source$

import numpy as np
import os

from .arithmetics import absolute, negate, add, subtract, multiply, divide, \
    power


__all__ = ["MRRArray", "empty", "empty_like",
           "zeros", "zeros_like",
           "save", "load",
           "mrr_min", "mrr_max", "mrr_mean",
           "cond_print"
           ]


class MRRArray(np.ndarray):
    '''
    Class to represent an image or a set of images including associated
    (standard)deviation per pixel and mask.

    The data per pixel is in the format ('phase','dev','mask')
    'phase' : float32
    'dev' : float32
    'mask' : boolean

    For consistency the phasedata should always be provided in multiples
    of 2*pi!

    'mask' is a boolean value per pixel, where False corresponds to an invalid
    pixel, and True corresponds to a valid pixel.

    Additional fields:
    original_file : string, The filename the data came from
    unwrapped : boolean, wether the data was unwrapped
    '''
    _attributes = {"orig_file": None,  # original filename
                   "unwrapped": False,  # flag whether the data was unwrapped
                   "PTFT": None,  # Trigger Delay, if present
                   "delta": None,
                   "Delta": None,
                   "maxGrad": None,
                   "gradstart": None,
                   "bvalue": None
                   }
    _datatype = {'names':   ['phase', 'dev', 'mask'],
                 'formats':  [np.float32, np.float32, np.bool_]
                 }

    def __new__(cls, input_array, dev=None, mask=None, **keypar):
        # cast to custom class
        obj = np.asarray(input_array, dtype=cls._datatype).view(cls)
        # set default values for dev and mask
        if mask is None:
            mask = True
        if input_array.dtype != cls._datatype:
            obj['dev'] = dev
            obj['mask'] = mask
        # add attributes
        for attr in keypar:
            if attr in cls._attributes:
                setattr(obj, attr, keypar[attr])
        for attr in cls._attributes:
            try:
                getattr(cls, attr)
            except AttributeError:
                setattr(cls, attr, cls._attributes[attr])
        # return new object
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        for attr in self._attributes:
            setattr(self, attr, getattr(obj, attr, self._attributes[attr]))
        # self.orig_file = getattr(obj, 'orig_file', None)

    def __array_wrap__(self, obj, context=None):
        return obj.view(type(self))

    def __str__(self):
        s = "MRRArray of shape {}".format(self.shape)
        if self.orig_file:
            s += ", read from file {}".format(self.orig_file)
        s += " ("
        if not self.unwrapped:
            s += "not "
        s += "unwrapped)\n"
        for name, unit in zip(
                ["delta", "bvalue", "Delta", "maxGrad", "gradstart", "POST"],
                ["ms"]*2 + ["s/mm^2"] + ["mT/m"] + ["ms"]*2):
            s += "\t{:<9}:{:>10.3f} {}\n".format(
                name, getattr(self, name), unit
                )
        s += super(MRRArray, self).__str__()
        return s

    # return fields as view of numpy.ndarray
    @property
    def mask(self):
        "Return a view of self's mask."
        return self['mask'].view(np.ndarray)

    @property
    def dev(self):
        "Return a view of self's standard deviation."
        return self['dev'].view(np.ndarray)

    @property
    def phase(self):
        "Return a view of self's phase."
        return self['phase'].view(np.ndarray)

    @property
    def variance(self):
        "Return a view of self's variance as square of the standard deviation."
        return self['dev'].view(np.ndarray)**2.

    # Arithmetics
    def __neg__(self):
        "Return self with negated phase values."
        return negate(self)

    def __abs__(self):
        "Return self with absolute phase values."
        return absolute(self)

    def __add__(self, other):
        "Add other to self, and return a new MRRArray."
        return add(self, other)

    def __radd__(self, other):
        "Add other to self, and return a new MRRArray."
        return add(self, other)

    def __sub__(self, other):
        "Subtract other from self, and return a new MRRArray."
        return subtract(self, other)

    def __rsub__(self, other):
        "Subtract other from self, and return a new MRRArray."
        return add(-self, other)

    def __mul__(self, other):
        "Multiply self with other, and return a new MRRArray."
        return multiply(self, other)

    def __rmul__(self, other):
        "Multiply self with other, and return a new MRRArray."
        return multiply(self, other)

    def __div__(self, other):
        "Divide self by other, and return a new MRRArray."
        return divide(self, other)

    def __rdiv__(self, other):
        "Divide self by other, and return a new MRRArray."
        return multiply(self**-1., other)

    def __truediv__(self, other):
        "Divide self by other, and return a new MRRArray."
        return divide(self, other)

    def __rtruediv__(self, other):
        "Divide self by other, and return a new MRRArray."
        return multiply(self**-1., other)

    def __pow__(self, other):
        "Raise self to the power of other, and return a new MRRArray."
        return power(self, other)

    def __iadd__(self, other):
        "Add other to self in place."
        add(self, other, out=self)
        return self

    def __isub__(self, other):
        "Substract other from self in place."
        subtract(self, other, out=self)
        return self

    def __imul__(self, other):
        "Multiply self with other in place."
        multiply(self, other, out=self)
        return self

    def __idiv__(self, other):
        "Divide self by other in place."
        divide(self, other, out=self)
        return self

    def __ipow__(self, other):
        "Raise self to the power of other in place."
        power(self, other, out=self)
        return self

    def min(self, axis=None, **kwargs):
        "Return the minimum phase value of self with a valid mask."
        return mrr_min(self, axis, **kwargs)

    def max(self, axis=None, **kwargs):
        "Return the maximum phase value of self with a valid mask."
        return mrr_max(self, axis, **kwargs)

    def mean(self, axis=None, weighted=True, unbias=True):
        "Return the mean of self's phase values along axis. See mrr.mean"
        return mrr_mean(self, axis, weighted, unbias)


# Array-Creation
# ---------------
def empty(shape, **keyargs):
    '''Creates an empty MRRArray with given shape.'''
    return MRRArray(np.empty(shape), **keyargs)


def empty_like(a):
    '''Creates an empty MRRArray of shape a.shape.'''
    return MRRArray(np.empty(a.shape))


def zeros(shape, **keyargs):
    '''Creates an MRRArray with given shape filled with zeros.'''
    return MRRArray(np.zeros(shape), **keyargs)


def zeros_like(a):
    '''Creates an MRRArray of shape a.shape filled with zeros.'''
    return MRRArray(np.zeros(a.shape))


def save(filename, a):
    np.save(filename, a)


def load(filename):
    array = MRRArray(np.load(filename))
    array.orig_file = os.path.abspath(filename)
    return array


# statistics
# ----------

def mrr_mean(a, axis=None, weighted=True, unbias=True):
    '''
    Creates a mean-image of an MRRArray ignoring masked values.

    result['phase'] contains the mean value of phase-data over the specified
    axis. result['dev'] is set to the (biased) standard deviation of phase data
    along the specified axis. result['mask'] will be False if any pixel along
    the given axis is False.

    Parameters
    ----------
    array : MRRArray
        values to average
    axis : int (optional)
        The axis along the mean is calculated. If axis = None, the mean is
        calculated over all pixels.
    weighted : bool (optional)
        Whether to use the variance as weights
    unbiased : bool (optional)
        If True, the weighted variance is being unbiased.

    Returns
    -------
    result : MRRArray
        The mean of array along the given axis
    '''
    if weighted:
        if not np.all(np.isfinite(np.sum(a.dev, axis=axis))):
            weighted = False
    weights = 1.0/a.variance if weighted else None

    # handle zeros in variance
    try:
        weights[np.where(~np.isfinite(weights))] = a.variance.mean()  # ok?
    except TypeError:
        pass

    p = np.ma.masked_where(a.mask is False, a.phase)

    av, v1 = np.ma.average(p, axis, weights, returned=True)
    res = MRRArray(av, orig_file=a.orig_file, unwrapped=a.unwrapped)

    if axis:
        if axis != 0:
            av = np.expand_dims(av, axis)
    var = np.ma.average((p-av)**2., axis=axis, weights=weights,
                        returned=False)
    if weighted and unbias:
        v2 = np.ma.sum(np.ma.masked_where(a.mask==False, weights**2.),
                       axis=axis)
        var /= (1 - v2/v1**2.)
    np.sqrt(var, out=res['dev'])

    res['mask'] = np.any(a.mask, axis=axis)

    return res


# miscellaneous
# -------------

def mrr_min(array, axis=None, **kwargs):
    "Return the minimal phase value of array that has a valid mask."
#    return array['phase'][array['mask']==True].min().item()
    return array.phase[array['mask'] is True].min(axis=axis, **kwargs)
    # ... .item()?


def mrr_max(array, axis=None, **kwargs):
    "Return the maximum phase value of array that has a valid mask."
#    return array['phase'][array['mask']==True].max().item()
    return array.phase[array['mask'] is True].max(axis=axis, **kwargs)


def cond_print(string, verbose=True):
    '''Prints string if verbose is valid.'''
    if verbose:
        print string
        return True
    else:
        return False