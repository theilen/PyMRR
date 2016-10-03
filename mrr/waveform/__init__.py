#version-history
#
#changes in 1.2.1:
#- switched to waveform v1.2.1
#
#chagnes in 1.2:
#- switched to import *
#- switched to waveform v 1.2
#- added __all__
#
#changes in 1.1:
#- switched to waveform version 1.1

__version__ = '1.2.1'
# $Source$

from .waveform import *
from . import calibration

__all__ = ['get_waveform',      #read-in waveforms
           'calibrate_linear',
           'scale_linear',
           'smooth',            #smooth wavefroms
           'calculate_start'    #estimate start of falling motion
           ]