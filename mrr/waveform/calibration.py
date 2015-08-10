# -*- coding: utf-8 -*-
"""
Created on Thu May 21 15:18:56 2015

@author: Sebastian Theilenberg
"""


import math
from scipy import odr
import numpy as np
import matplotlib.pyplot as plt
from warnings import warn
import pickle


def cal_from_file(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)


class Calibrate(object):

    def __init__(self, x=None, y=None, dx=None, dy=None):
        if np.any(x):
            if not np.any(y):
                raise ValueError("no y-values provided!")
            self.set_values(x, y, dx, dy)
        self.result = None

    def __call__(self, x, std=False):
        self._has_run()
        y = self._func_values(x)
        if std:
            err = self._func_err(x)
            return y, err
        return y

    @property
    def values(self):
        return self._x, self._y, self._dx, self._dy

    def set_values(self, x, y, dx=None, dy=None):
        assert len(y) == len(x)
        self._x = np.asarray(x)
        self._xm = self._x.mean()
        self._y = np.asarray(y)
        self._ym = self._y.mean()
        if dx is not None:
            assert len(dx) == len(x)
            self._dx = np.asarray(dx)
        if dy is not None:
            assert len(dy) == len(x)
            self._dy = np.asarray(dy)
        self._correct_to_cg()

    @property
    def slope(self):
        self._has_run()
        return np.hstack((self.result.beta, self.result.sd_beta))

    def _correct_to_cg(self):
        self._x -= self._xm
        self._y -= self._ym

    def save(self, filename):
        self._has_run()
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    def _has_run(self):
        if self.result is None:
            raise AttributeError("No calibration present. Run Calibrate.run()"
                                 "first")

    def run(self):
        Model = odr.Model(self._linear)
        Data = odr.RealData(self._x, self._y,
                            sx=self._dx, sy=self._dy)
        initial = self._initial_guess()
        solver = odr.ODR(Data, Model, initial)
        self.result = solver.run()
        return self.result.beta, self.result.sd_beta

    def _linear(self, P, x):
        "Linear function of regression."
        return P[0]*x  # + P[1]

    def _initial_guess(self):
        dy = self._y.max() - self._y.min()
        dx = self._x.max() - self._x.min()
        return (dy/dx,)  # self._y[0] - self._x[0]*dy/dx)

    def _func_values(self, x):
        "Returns the calibrated values f(x) as np.array."
        x = np.array(x) - self._xm
        # m, n = self.result.beta
        # return x*m + n
        return self._linear(self.result.beta, x) + self._ym

    def _func_err(self, x):
        x = np.asarray(x)
        dm = self.result.sd_beta
        return np.abs(dm * (x - self._xm))

    def plot(self, title=None, annotation=True):
        self._has_run()
        # find values in original system
        t_range = (self._x.max() - self._x.min()) * 0.05
        t_ = np.linspace(self._x.min() + self._xm - t_range,
                         self._x.max() + self._xm + t_range,
                         100)
        err = self._func_err(t_)
        y = self._func_values(t_)

        plt.figure()
        if title is not None:
            plt.title(title)
        plt.minorticks_on()
        # plot margin
        plt.fill_between(t_, y - err, y + err, alpha=0.2, color='teal')
        # plot function
        plt.plot(t_, y, color='teal')
        # plot measurement values
        plt.errorbar(self._x + self._xm, self._y + self._ym,
                     yerr=self._dy, xerr=self._dx,
                     ls='none', color='blue')
        # adjust xlim
#        if plt.xlim()[0] < 0:
#            plt.xlim(xmin=0)
        # labels
        plt.xlabel(r"$U\,\,\left[\si{\milli\volt}\right]$")
        plt.ylabel(r"$x\,\,\left[\si{\milli\meter}\right]$")
        # add function
        if annotation:
            # check latex settings
            if plt.rcParams["text.usetex"] is True:
                if ("\usepackage{siunitx}"
                        not in plt.rcParams["text.latex.preamble"]):
                    warn("Matplotlib.pyplot not fully configured to use"
                         "tex with siunitx. Add '\usepackage{siunitx}'"
                         "to rcParams['text.latex.preambel']",
                         UserWarning)
                else:
                    annotation = self._create_annotation()
                    plt.annotate(annotation, xy=(0.08, 0.82),
                                 xycoords="axes fraction")

    def _create_annotation(self):
        m = self.result.beta[0]
        dm = self.result.sd_beta[0]
        # find precisions of last relevant decimal
        precision = int(math.floor(-math.log10(dm)) + 1)
        # calculate the decimals of the errors at these positions
        err = int(round(dm*10**precision))
        # correct for rounding 0.09 to 0.1
        if err > 9:
            err /= 10
            precision -= 1
        # create annotation string
        annotation = (r"$x = \SI{{{0:.{2}f}\pm{1:.{2}f}}}"
                      r"{{\milli\meter\per\milli\volt}}\cdot U$"
                      ).format(m, dm, precision)
        return annotation
