# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 14:28:39 2014

@author: Sebastian Theilenberg
"""

#version-history
#
#changes in version 1.2.1
#- bugfixes in calculate_start und _plot_calculation
#
#changes in version 1.2
#- removed get_optic, scale_waveform & calculate_scaling
#- added calibrate_linear & scale_linear
#- added calculate_start
#- added __all__
#
#changes in version 1.1
#- renamed scale_waveform() to calculate_scaling()
#- added new scale_waveform()
#- renamed get_waveform() to get_optic()
#- added new get_waveform()
#- added read-in support for old oscilloscope

__version__ = '1.2.1b'
# $Source$


import numpy as np
from scipy import optimize
from scipy.stats import linregress
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os.path

from .tektronix import read_isf_files


__all__ = ['get_waveform',      # read-in waveform
           'calibrate_linear',
           'scale_linear',      # scale a waveform by a linear function
           'smooth',            # smooth wavefroms
           'calculate_start'    # estimate start of falling motion
           ]


_read_sets = {'new': {'delimiter': ',',  # Spalten-Trennung
                      'skip_header': 21,  # Header-Zeilen
                      'unpack': True},
              'old': {'delimiter': ',',  # Spalten-Trennung
                      'skip_header': 6,  # Header-Zeilen
                      'unpack': True,
                      'start_col': 3}
              }


def get_waveform(filename, skip=5, read_in='new', **kwargs):
    '''
    Returns an array containing timepoints and the data of all saved channels.
    To limit memory-use, not every timepoint is returned (to adjust this
    behaviour see the keyword skip).

    The order of the returned array is as follows:

    index   data
    ---------------------------------------------------
    0       timepoints
    1...    channels


    Parameters
    ----------
    filename : string
        path to the textfile containing the waveforms
    skip : integer (optional)
        number of timepoints to skip after every timepoint. Default: 5
    read_in : string
        which predefined set to use for the read_in of the file. Choose
        either 'new' for the new small oscilloscope, or 'old' for the old
        bigger one. Default: 'new'

    If read_in='old' is used, the file is supposed to contain only one channel,
    e.g. the first 3 columns are empty. These columns are discarded!
    If this behaviour is unwanted, use
        start_col=0
    to get all channels.

    Keyarguments are passed to np.genfromtxt used for loading the
    waveform-data.
    To adjust use:
        delimiter : Symbol indicating fields.
        skip_header :  number of invalid rows in the beginning of the file.
    '''
    # TODO update docstring to incorporate binary files
    filename = os.path.abspath(filename)
    binary_extensions = set([".isf"])
    if os.path.splitext(filename)[-1].lower() in binary_extensions:
        return read_isf_files(filename)[:, ::skip]

    try:
        D = _read_sets[read_in]
    except KeyError:
        raise KeyError('Unknown Read-In-Set "{}"'.format(read_in))
    D.update(kwargs)

    try:
        col = D.pop('start_col')
    except KeyError:
        col = 0

    wave = np.genfromtxt(filename, **D)
    return wave[col:, ::skip]


def find_pulse_maximum(time, signal, radius=1.5e-3):
    """
    Returns the time point of the maximum singal amplitude.

    The maximum amplitude is determined by quadratically fitting the peak
    found in the region determined by radius around t=0.
    """
    sample_rate = np.diff(time[:2])[0]
    zero = np.where(np.isclose(time, 0.0))[0]
    radius = int(radius/sample_rate)
    lower, higher = np.where(
        signal[zero-radius:zero+radius] > 0.01
        )[0][[0, -1]]
    region = slice(zero-radius+lower, zero-radius+higher)
    a, b, c = np.polyfit(time[region], signal[region], deg=2)
    return -b/2/a


def optimal_smoothing_length(times, length=2e-3):
    sample_interval = np.diff(times[:2])[0]
    smoothing = int(length/sample_interval)
    if smoothing % 2 == 0:
        smoothing -= 1
    return smoothing


def get_mean_waveform(filenames, pulse_radius=1.5e-3, sample_int=5e-5,
                      smoothing=True,
                      skip=5, read_in='new', **kwargs):
    """
    Read in a set of waveforms and return the mean signal.

    The waveforms are interpolated to a new set of uniformly spaced timepoints.
    Before averaging the data, the times of the individual waveforms is
    corrected to the maximum of the pulse the aquisition was triggered to.

    Parameters:
    -----------
    filenames : list
        List of filenames of the waveforms
    pulse_radius : float
        Radius around t=0 to search for the maximum of the pulse
    sample_int : float
        Sample interval of the returned timepoints
    smoothing : Bool
        Whether to smooth the waveforms before averaging them
    """
    waves = []
    for f in filenames:
        if not os.path.isfile(f):
            continue
        wave = get_waveform(f, skip=skip, read_in=read_in, **kwargs)
        # correct to maximum of pulse
        t0 = find_pulse_maximum(wave[0], wave[1], radius=pulse_radius)
        wave[0] -= t0
        # smooth waveform
        if smoothing:
            smoothing = optimal_smoothing_length(wave[0])
            for c in range(2, 5):
                wave[c] = smooth(wave[c], window_len=smoothing)
        # save time and optical channels
        waves.append(wave[[0, 2, 3, 4]])
    # create new sample points
    tmin = max([a[0, 0] for a in waves])
    tmax = min([a[0, -1] for a in waves])
    samples = int((tmax - tmin)/sample_int) + 1
    t_sample = np.linspace(tmin, tmax, samples, endpoint=True)
    # create uniformly sampled data
    combined = np.asarray(
        [interp1d(w[0], w[1:], axis=1)(t_sample) for w in waves]
        )
    mean = np.asarray([combined.mean(axis=0), combined.std(axis=0)])
    return t_sample, mean


def calibrate_linear(voltages, positions, full=False):
    """
    Calibrate waveforms by a set of (Voltage, position).

    Calculates slope and intercept of a linear regression through the data
    using scipy.stats.linregress
    The returned parameters assume that 
        
        x = f(U) = a*U + b
        
    Parameters
    ----------
    voltages : array-like
        Values of measured voltage
    positions : array-like
        Corresponding measured positions
    full : bool (optional)
        Whether to return statistical data of the fit.
        
    Returns
    -------
    slope : float
    intercept : float
    r-value : float
        correlation coefficient
    p-value : float
        two-sided p-value
    std : float
        standard error of teh estimate
    """
    sl, inter, r, p, std = linregress(voltages, positions)
    if full:
        return sl, inter, r, p, std
    else:
        return sl, inter
    
    
def scale_linear(waveform, slope, intercept, inplace=False):
    """
    Scale a waveform by applying a linear calibration. 
    
    Parameters
    ----------
    waveform : ndarray
        the waveform to be scaled.
    slope : float
        slope of the linear calibration function.
    intercept : float
        intercept of the linear calibration function.
    inplace : bool
        Whether to return the scaled waveform in a new array (inplace=False),
        this is the default behaviour, or override the unscaled waveform 
        (inplace=True).
        
    Returns
    -------
    scaled : ndarray
        Only if inplace==False. Scaled array of the same shape as waveform.
    """
    result = waveform if inplace==True else waveform.copy()
    result *= slope
    result += intercept
    
    if not inplace: return result
    
    

#def scale_waveform(wave, displacement, urange=500, lrange=None, inplace=True):
#    '''
#    
#    '''
#    #parse slicing
#    if hasattr(urange, '__len__'):
#        if not len(urange)  in [2,3]:
#            raise ValueError()
#        else:
#            upper = slice(*urange)
#    else:
#        upper = slice(0,urange)
#    if lrange == None:
#        lower = slice(-(upper.stop-upper.start), None, upper.step)
#    else:
#        if hasattr(lrange, '__len__'):
#            if not len(lrange) in [2,3]:
#                raise ValueError()
#            else:
#                lower = slice(*lrange)
#        else:
#            lower = slice(lrange)
#
#    result = wave if inplace==True else wave.copy()
#    
#    #scale      
#    y_max = wave[upper].mean()
#    y_min = wave[lower].mean()
#    result -= y_max
#    result *= displacement/(y_max - y_min)
#        
#    if not inplace: return result
#    
#    
#def scale_waveforms(wave, displacement, urange=500, lrange=None, inplace=True):
#    """
#    Scale a set of waveforms
#    """
#    #test displacements
#    try:
#        l_ = len(displacement)
#    except TypeError:
#        displacement = [displacement] * len(wave)
#    else:
#        if l_ != len(wave):
#            displacement = [displacement[0]] * len(wave)
#            print "Length displacement does not match number of channels. Using first displacement for all present channels!"
#            
#    result = wave if inplace==True else wave.copy()
#    
#    for c in range(result.shape[0]):
#        scale_waveform(result[c], displacement[c], urange, lrange, inplace=True)
#        
#    if not inplace: return result
            
                
def smooth(x, window_len=11, window='hanning'):
    '''
    Smooth a waveform. Taken from http://wiki.scipy.org/Cookbook/SignalSmooth
    on 15-04-2014.
    '''   
    if x.ndim != 1:
        raise ValueError, 'smooth only accepts 1-dimensional data!'
        
    if x.size < window_len:
        raise ValueError, 'Input data needs to be bigger than window size!'
        
    if window_len < 3:
        return x
        
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError,\
        "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
  
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]

    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)-1]


def _plot_calculation(t, data, thresholds, opt_func, p, stop):
    "Plot the key-values of the calculation of the starting point."
    fig = plt.figure()
    plt.plot(t, data, color='grey')
    dist = len(t)
    plt.plot(t[dist/4:stop], opt_func[dist/4:stop], color='red')
    for val in thresholds:
        plt.axhline(val, ls='--', color='black')
    plt.axhline(p[-1], ls='-', color='blue')    #y0
    plt.axvline(p[1], ls='-', color='blue')     #t0
    

def calculate_start(t, wave, plower=0.3, pupper=0.055, height_fac=2.,
                    smooth_data=False, full_output=False, verify_plot=False, 
                    verbose=False):
    """
    Estimates the starting point in time of the falling motion of a waveform.
    
    Implementation of Björn Schemmann's algorithm.
    
    Parameters
    ----------
    t : 1darray
        time-points (x-data). usually waveform[0].
    wave : 1darray
        wave-data to be fitted.
    plower : float
        percentage of the falling height to be used for the estimation.
    pupper : float
        percentage of the top part of the falling height not to be used for 
        the estimation.
    smooth_data : bool
        Whether to smooth the y-data when calculating thresholds.
    full_output : bool
        True to return all optional outputs
    verify_plot : bool
        Creates a plot to verify the algorithm.
        
    Returns
    -------
    t0 : float
        The estimated starting point of the falling motion
    a0 : float (Optional)
        The starting acceleration of the falling motion.
    y0 : float (Optional)
        The starting height of the falling motion as used in the 
        estimation.
    """
    
    
    fitfunc = lambda p, x, y0: 0.5*p[0]*((x-p[1])**2.) + y0
    errfunc = lambda p, x, y, y0: fitfunc(p, x, y0) - y
    
    sdata = smooth(wave, 11) if smooth_data else wave
        
    #initialise fit
    ma, mi = sdata.max(), sdata.min()
    upper = ma - np.abs(pupper*ma) #ma * (1-pupper)
    lower = ma - plower*(ma-mi)
    start, stop = np.where(sdata < upper)[0][0], np.where(sdata > lower)[0][-1]
    x = t[start:stop]
    data = wave[start:stop]
    
    #estimate initial guess
    seed = np.polyfit(x, data, 2)
    p0 = (seed[0], -seed[1]/(2*seed[0]))    #a_0, t_0 (apex)
    #estimate starting height
    y0_slice = slice(start-height_fac*(stop-start),
                     start-height_fac*0.2*(stop-start)
                     )
    y0 = np.mean(wave[y0_slice])

    if verbose:
        s = 'min, max = {:.3}, {:.3}'.format(ma, mi)
        s += '\nthresholds: {:.3}, {:.3}'.format(upper, lower)
        s += '\nfitarea: {}({:.3}), {}({:.3})'.format(start, t[start], stop, t[stop])
        s += '\nseed: ({}, {}, {})'.format(*seed)
        print s

    #fit
    popt, _ = optimize.leastsq(errfunc, p0[:], args=(x,data,y0))
    
    if verify_plot:
        dist = stop-start
        drange = slice(start-2.5*dist, stop+dist)
        _plot_calculation(t[drange], wave[drange], 
                          thresholds=(upper, lower), 
                          p=(popt[0], popt[1], y0), 
                          stop=stop,
                          opt_func=fitfunc(popt, t[drange], y0),
                          )
    
    if full_output:
        return popt[1], popt[0], y0
    else:
        return popt[1]
    
        
    
#def calculate_scaling(wave, yrange, const_len=500):
#    '''
#    Calculates scale-factors of a waveform.
#    
#    Calculates:
#        ymax : maximum value of <wave> (mean of the first <const_len> datapoints)
#        ymin : minimum value of <wave> (mean of the last <const_len> datapoints)
#        scalefactor = yrange/(ymin-ymax)
#    
#    Returns (scalefactor, max-value)
#    '''
#    ymax = wave[:const_len].mean()
#    ymin = wave[-const_len:].mean()
#    y = ymax - ymin
#    
#    if y == 0:
#        print 'WARNING: ymax equals ymin, no scaling is applied.'
#        return (yrange, ymax)
#    else:
#        return (yrange/y, ymax)        
#    p[1]*((x-p[0])**2.) + p[2]
#    
#def apply_scaling(wave, scalefactor, offset=0):
#    '''
#    Scales w waveform IN PLACE by adjusting the yrange by <scalefactor> and 
#    substracting <offset>.
#    
#    Calculates (wave - offset) * scalefactor
#    
#    Parameters
#    ----------
#    wave : array
#        the wave to be scaled. The wave is changed in-place, so the original 
#        data is lost.
#    scalefactor : float
#        factor to scale the data with
#    offset : float (optional)
#        offset to substract. Default: 0
#    '''
#    wave = (wave - offset) * scalefactor
    
    

        
        

                         
                         
    
#def get_optic(filename, yrange, skip=5, read_in='new', **kwargs):
#    '''
#    Returns an array containing timepoints and channels 1 to 4 of a saved 
#    waveform as well as smoothed waveforms for channels 2 to 4.
#    The file is supposed to contain the antenna-data as channel 1 and optical 
#    measurements in the following channels.
#    
#    To limit memory-use, not every timepoint is returned (to adjust this 
#    behaviour see the keyword skip).
#    
#    The order of the returned array is as follows:
#    
#    index   data
#    ---------------------------------------------------
#    0       timepoints
#    1       channel 1 (antenna)
#    2-...   channels 2-4 (optical measurement) scaled raw data 
#            (usefull for error estimations)
#    ..-...  channels 2-4 smoothed and normalized
#    
#    
#    Parameters
#    ----------
#    filename : string
#        path to the textfile containing the waveforms
#    yrange : float
#        maximum displacement to normalize the waveform to.
#    skip : integer (optional)
#        number of timepoints to skip after every timepoint. Default: 5
#    read_in : string
#        which predefined set to use for the read_in of the file. Choose
#        either 'new' for the new small oscilloscope, or 'old' for the old 
#        bigger one. Default: 'new'
#    
#    
#    If read_in='old' is used, the file is supposed to contain only one channel,
#    e.g. the first 3 columns are empty. These columns are discarded!
#    If this behaviour is unwanted, use
#        start_col=0
#    to get all channels.
#    
#    Keyarguments are passed to np.genfromtxt used for loading the waveform-data.
#    To adjust use:
#        delimiter : Symbol indicating fields.
#        skip_header :  number of invalid rows in the beginning of the file.
#    '''
#    wave = get_waveform(filename, skip=1, read_in=read_in, *kwargs)
#    
#    wform = np.empty((2*wave.shape[0]-2, wave.shape[1]) )
#    wform[0:2] = wave[0:2].copy()
#    for i in range(2,wform.shape[0]):
#        opt = smooth(wave[i], window_len=201)
#        scale, offset = calculate_scaling(opt, yrange)
#        wform[i] = (wave[i].copy() - offset) * scale
#        wform[i+3] = ((opt - offset) * scale)
#        
#    return wform[:,::skip]
    
    

                             
