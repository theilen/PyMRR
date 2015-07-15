# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:15:19 2015

@author: theilenberg
"""


import os
import time
import shutil
import dicom
import numpy as np
import pickle

from ..read import parse_parameters


_extensions = set(['.isf', '.csv'])

def copy_wavefiles(source, destination, timezoneadjust=10000, test=False):
    source = os.path.abspath(source)
    destination = os.path.abspath(destination)
    
    for f in os.listdir(source):
        sourcepath = os.path.join(source, f)
        if not os.path.isfile(sourcepath):
            if test: print '{} not a file'.format(f)
            continue
        base, ext = os.path.splitext(f)
        if ext.lower() not in _extensions:
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


def create_dicom_times(dirpath, series=[], skip_study=[], toff=0.):
    """
    Parameters:
    -----------
    dirpath : str
        path to the directory of the dicom-files
    series : list
        list of integers indicating series numbers to use
    skip_study : list
        list of integers indiciating study-IDs to ignore
    toff : float
        time in seconds to adjust the tomograph time to the oscilloscope time.
        t_mri = t_osci + toff
    """
    dcmtimes = []
    if not series:
        series = range(99)
    series_list = set([str(i) for i in series])
    study_exceptions = set([str(i) for i in skip_study])

    files = os.listdir(dirpath)
    for i in reversed(range(len(files))):  # TODO reversed needed?
        f = files[i]
        t_ = f.split('_')
        if len(t_) < 4:
            files.remove(f)
        elif (t_[0] in study_exceptions) or (t_[-2] not in series_list):
            files.remove(f)
            
    for f in files:
        fabs = os.path.join(dirpath, f)
        try:
            dcm = dicom.read_file(fabs, stop_before_pixels=True)
        except IOError:
            continue
        try:
            tau = parse_parameters(dcm)['PTFT']
        except:
            continue
        # read time in seconds since epoch
        act, subsec = dcm.AcquisitionTime.split('.')
        date = dcm.AcquisitionDate
        t_ = time.strptime(date + act, "%Y%m%d%H%M%S")
        t_epoch = time.mktime(t_)
        # correct time for oszitime (t_mri = t_oszi + 147s)
        t_epoch -= toff
        # add subseconds
        t_epoch += float("0." + subsec)
        dcmtimes.append((tau, f, t_epoch))
        
    dcmtimes.sort(key=lambda x: x[-1])
    return dcmtimes
        

def create_osci_times(dirpath, acquisitiondate):
    """
    Parameters:
    -----------
    dirpath : str
        path to the directory of the waveform-files
    acquisitiondate : str
        str indiciating the date of the waveform-acquisition. "YYYYMMDD"
    """
    osclist = []

    for f in os.listdir(dirpath):
        fabs = os.path.join(dirpath, f)
        if not os.path.isfile(fabs):
            continue
        elif not os.path.splitext(f)[-1].lower() in _extensions:
            continue
        elif not "CH1" in f:
            continue
        base, ext = os.path.splitext(f)
        traw = base.split('_')[-1]
        tstring, subsec = traw.split('.')
        t_ = time.strptime(acquisitiondate + tstring, "%Y%m%d%H%M%S")
        t_epoch = time.mktime(t_)
        # add subseconds
        t_epoch += float("0." + subsec)
        osclist.append((f, t_epoch))
    
    osclist.sort(key=lambda x: x[-1])
    return osclist


def _find_mean_difference(datalist, mean=3.0):
    data = np.diff([item[-1] for item in datalist])
    indices = np.where(np.isclose(data, mean, atol=1.0))
    return data[indices].mean()


def map_files(dcmlist, osclist, tr=3.0, verbose=False):
    
    if verbose:
        print "{} dicom files, {} waveforms".format(len(dcmlist), len(osclist))
    
    d_dcm = _find_mean_difference(dcmlist)
    d_osc = _find_mean_difference(osclist)
    delta = d_dcm - d_osc
    filemappings = []
    filemapping = {}
    series = 0
    
    for i in range(len(dcmlist)):
        dcm = dcmlist[i]
        # new series?
        if int(dcm[1].split('_')[-2]) != series:
            if filemapping:
                filemappings.append(filemapping)
            filemapping = {}
            series = int(dcm[1].split('_')[-2])
            counter = 0
            if verbose:
                print "New Series", series
        
        t_epoch = dcm[-1]
        match_item = [
            dcm[0], dcm[1],
            time.strftime("%H:%M:%S", time.localtime(t_epoch))
            ]
        match_item[-1] += "{:.3f}".format(t_epoch % 1)[1:]
        try:
            match = filter(lambda x: (
                (t_epoch - x[-1] - delta*counter) < 0. and 
                (t_epoch - x[-1] - delta*counter) > -tr), 
                osclist)[0]
        except IndexError:
            match_item.extend(["-", "-"])
        else:
            match_item.extend([
                match[0], time.strftime("%H:%M:%S", time.localtime(match[-1]))
                ])
            match_item[-1] += "{:.3f}".format(match[-1] % 1)[1:]
        finally:
            filemapping[dcm[1]] = tuple(match_item)
    
    filemappings.append(filemapping)
    return filemappings


def _create_latex_header(title=""):
    head = latex_head['preamble']
    head += r"\author{{{}}}".format(latex_head['title'])
    head += latex_head['header']
    return head

def _create_latex_foot():
    return latex_head['foot']
    
latex_head = {'preamble': r"""\documentclass[8pt]{{article}}
\usepackage{siunitx}
\usepackage{longtable}
\usepackage[top=1cm, bottom=1.5cm]{geometry}
\title{Mapping Dicom and Oszilloscope Files}""",
'title': "",
'header': r"""\date{}
\begin{document}
\maketitle
\begin{longtable}{c||l|l||c|c|l}
$\tau$ [\si{\milli\second}] 
& dicom
& waveform
& dicom time
& waveform time
& Comments \\
\hline\hline
\endhead
""",
'foot': r"""\end{longtable}
\end{document}
"""
}


def print_mapping(filename, filemapping, title=""):
    with open(filename + '.tex', 'w') as f:
        f.write(_create_latex_header(title))

        keys = filemapping.keys()
        keys.sort(key=lambda x: int(x.split('_')[-1]))
        keys.sort(key=lambda x: int(x.split('_')[-2]))
        keys.sort(key=lambda x: float(filemapping[x][0]))

        for i in range(len(keys)):
            item = filemapping[keys[i]]
            fstring = "{:>5.1f} & {:>18} & {:>7} & {} & {} & \\\\\n".format(
                item[0], item[1].replace('_', '\_'),
                item[3][:7], item[2], item[4]
                )
            f.write(fstring.replace("\\", "\\"))
        f.write(_create_latex_foot())


def save_mapping(filename, filemapping):
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(filemapping, f, protocol=pickle.HIGHEST_PROTOCOL)
    return True

def load_mapping(filename):
    return pickle.load(filename)