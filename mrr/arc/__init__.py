#

__version__ = '1.2.1'
# $Source$

#version history
#=================
#
#version 1.2.1
#- reimplemented to support mrrcore v1.2
#
#version 1.2
#added parameters G and delta to calc_displacement
#
#version 1.1
#- switched to relative imports


from .arcmr import *

__all__ = ['normalize', 'calc_displacement']



