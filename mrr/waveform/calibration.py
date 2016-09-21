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
        Cal = pickle.load(f)
        if not isinstance(Cal, Calibrate):
            raise AttributeError(
                "File {} did not return a Calibrate object!".format(f)
            )
        return Cal


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
            err = self._func_std(x)
            return y, err
        return y

    @property
    def values(self):
        return self._x, self._y, self._dx, self._dy

    def set_values(self, x, y, dx=None, dy=None):
        assert len(y) == len(x)
        self._x = np.array(x, copy=True)
        self._y = np.array(y, copy=True)
        if dx is not None:
            assert len(dx) == len(x)
            self._dx = np.asarray(dx)
        if dy is not None:
            assert len(dy) == len(y)
            self._dy = np.asarray(dy)

    @property
    def slope(self):
        "Return the slope of the calibration"
        self._has_run()
        return np.hstack((self.result.beta[0],
                          # correct std error to std dev. by residual variance
                          self.result.sd_beta[0]/np.sqrt(self.result.res_var))
                         )

    @property
    def intercept(self):
        "Return the intercept of the calibration"
        return np.hstack((self.result.beta[1],
                          # correct std error to std dev. by residual variance
                          self.result.sd_beta[1]/np.sqrt(self.result.res_var))
                         )

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
        # return self.result.beta, self.result.sd_beta

    def _linear(self, P, x):
        "Linear function of regression."
        return P[0]*x + P[1]

    def _initial_guess(self):
        dy = self._y.max() - self._y.min()
        dx = self._x.max() - self._x.min()
        m = dy/dx
        return (m, self._y[0] - self._x[0]*m)

    def _func_values(self, x):
        "Returns the calibrated values f(x) as np.array."
        x = np.array(x)
        return self._linear(self.result.beta, x)

    def _func_std(self, x):
        "Return the standard deviation on the calibrated values y(x)"
        g = np.array([x, 1])
        return np.sqrt(np.dot(g, np.dot(self.result.cov_beta, g.T)))

    def _func_err(self, x):
        "Return the standard error on the calibrated values y(x)"
        raise NotImplementedError("How to do it? How??? And why?")

    def plot(self, title=None, annotation=False):
        self._has_run()
        # find values in original system
        t_range = (self._x.max() - self._x.min()) * 0.05
        t_ = np.linspace(self._x.min() - t_range,
                         self._x.max() + t_range,
                         100)
        std = self._func_std(t_)
        y = self._func_values(t_)

        plt.figure()
        if title is not None:
            plt.title(title)
        plt.minorticks_on()
        # plot margin
        plt.fill_between(t_, y - std, y + std, alpha=0.2, color='teal')
        # plot function
        plt.plot(t_, y, color='teal')
        # plot measurement values
        plt.errorbar(self._x, self._y,
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
        warn("Annotation not implemented!")
        return ""
        # TODO: create with intercept
        m, dm = self.slope
        c, dc = self.intercept
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
