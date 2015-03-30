# -*- coding: utf-8 -*-
"""
Created on Fri Dec 05 18:05:48 2014

@author: theilenberg
"""

__version__ = '1.0'

# Version history
# ===============
#
# version 1.0
# - added SortedList
# - removed trim_adjoin and sort_adjoin
# - reimplemented all functions to use SortedList

import numpy as np
from bisect import bisect_right, bisect_left

from ..array_manipulation import get_neighbours
from .. import constants as cnst


def update_adjoin_list(adjoin_list, pixel,
                       avoid=None, flagarray=None, qualitymap=None,
                       trim=False, max_size=None):
    '''
    Update <adjoin_list> with the neighbours of <pixel>.

    Parameters
    ----------
    adjoin_list : list
        List to add pixels to.
    pixel : (y,x)
        Pixel to find neighbours of.
    avoid : bitflag (optional)
        pixel with this bitflag are not appended.
    flagarray : 2darray (optional)
        array containing bitflags.
        Mandatory if avoid is not None or trim is True
    qualitymap : 2darray (optional)
        array containing the qualitymap.
        Mandatory if sort_list or trim is True
    trim : bool (optional)
        whether to trim the list. Trimming expects adjoin_list to be of type
        SortedList!
    max_size : int (optional)
        maximum list length if trim is True
    '''
    if np.any(flagarray):
        shape = flagarray.shape
    else:
        shape = None

    for p in get_neighbours(*pixel, shape=shape):
        if avoid:
            if flagarray[p] & avoid:
                continue
        if p in adjoin_list:
            continue
        if trim:
            if len(adjoin_list) > max_size:
                if adjoin_list.key_function(p) < adjoin_list.get_key(0):
                    flagarray[p] |= cnst.POSTPONED
                    continue
        adjoin_list.append(p)

    if trim:
        trim_sorted_list(adjoin_list, max_size, flagarray)


def trim_sorted_list(list_, max_len, flagarray):
    length = len(list_)
    if length > max_len:
        flagarray[zip(*list_[:length//2])] |= cnst.POSTPONED
        del list_[:length//2]
        return True
    return False


def refill_adjoin(adjoin_list, flagarray, max_size):
    post = zip(*np.where(flagarray & cnst.POSTPONED))
    if not np.any(post):
        return False
    adjoin_list.extend(post)
    # unmark pixel
    flagarray[zip(*post)] &= ~cnst.POSTPONED
    # trim
    trim_sorted_list(adjoin_list, max_size, flagarray)


class SortedList(object):
    '''
    Sorted List after
    http://code.activestate.com/recipes/577197-sortedcollection/
    '''

    def __init__(self, iterable=(), key=None):
        self._given_key = key
        key = (lambda x: x) if key is None else key
        items = sorted((key(item), item) for item in iterable)
        self._keys = [k for k, item in items]
        self._items = [item for k, item in items]
        self._key = key

    @property
    def key_function(self):
        '''key function'''
        return self._key

    @key_function.setter
    def key_function(self, key):
        if key is not self._key:
            self.__init__(self._items, key=key)

    @key_function.deleter
    def key_function(self):
        self._key = None

    def __len__(self):
        return len(self._items)

    def __getitem__(self, key):
        return self._items[key]

    def __delitem__(self, key):
        del self._items[key]
        del self._keys[key]

    def __iter__(self):
        return iter(self._items)

    def __reversed__(self):
        return reversed(self._items)

    def __contains__(self, item):
        k = self._key(item)
        l = bisect_left(self._keys, k)
        r = bisect_right(self._keys, k)
        return item in self._items[l:r]

    def __repr__(self):
        return '{}({} key:{})'.format(
            self.__class__.__name__,
            self._items,
            getattr(self._given_key, '__name__', repr(self._given_key))
            )

    def count(self, item):
        k = self._key(item)
        l = bisect_left(self._keys, k)
        r = bisect_right(self._keys, k)
        return self._items[l:r].count(item)

    def index(self, item):
        'Return the position of item.'
        k = self._key(item)
        l = bisect_left(self._keys, k)
        r = bisect_right(self._keys, k)
        return self._items[l:r].index(item) + l

    def insert(self, item):
        'insert item, if equal key is found insert to the right'
        k = self._key(item)
        index = bisect_right(self._keys, k)
        self._keys.insert(index, k)
        self._items.insert(index, item)

    def pop(self):
        'Return and remove last item.'
        del self._keys[-1]
        return self._items.pop()

    def remove(self, item):
        'Remove first occurence of item.'
        index = self.index(item)
        del self._keys[index]
        del self._items[index]

    def append(self, item):
        'Equal to SortedList.insert(item).'
        self.insert(item)

    def extend(self, iterable):
        'Add a sequence.'
        for item in iterable:
            self.insert(item)

    def get_key(self, i):
        'Return key of the item at index i.'
        return self._keys[i]

    def get_keys(self):
        'Return keys of self.'
        return self._keys[:]
