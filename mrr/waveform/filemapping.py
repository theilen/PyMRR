# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:15:19 2015

@author: theilenberg
"""


import os
import time
import shutil


def copy_wavefiles(source, destination, timezoneadjust=10000, test=False):
    source = os.path.abspath(source)
    destination = os.path.abspath(destination)
    
    extensions = set(['.isf', '.csv'])

    for f in os.listdir(source):
        sourcepath = os.path.join(source, f)
        if not os.path.isfile(sourcepath):
            if test: print '{} not a file'.format(f)
            continue
        base, ext = os.path.splitext(f)
        if ext.lower() not in extensions:
            continue
        # find creation time
        t_epoch = os.path.getctime(sourcepath)
        tstring = time.strftime("%H%M%S", time.gmtime(t_epoch))
        tstring += "{:4.4f}".format(t_epoch % 1.)[1:]
        tfloat = float(tstring)
        tfloat += timezoneadjust
        # copy file
        copypath = os.path.join(destination,
                                '{}_{:.4f}{}'.format(base, tfloat, ext)
                                )
        if test:
            print sourcepath, copypath
            continue
        if not os.path.exists(copypath):
            shutil.copy2(sourcepath, copypath)
        
