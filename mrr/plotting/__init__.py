# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:53:47 2014

@author: Sebastian Theilenberg
"""

__version__ = '1.4'
# $Source$


#version history
#=================
#version 1.4
#- changed plot() to work with MRRArray v1.2
#- added keyarg 'labels' to plot()
#
#version 1.3.1
#- bugfix in display_plain() (ValueError when crop=False)
#- added return_index to plot()
#- switched to relative imports
#
#version 1.3
#- added crop_image()
#- added valid_mask()
#- changed display and display_plain() to use crop_image
#
#version 1.2
#- slightly modified display()
#   - added calc_ratio() and crop_limits()
#- added display_plain()
#
#version 1.1
#- modified display():
#   - fixed cropping
#   - added keyword hold to keep figure alive
#   - added support for keyargs passed through to matplotlib.imshow
#
#version 1.0
#---------------
#moved plot() from writy.py to plot.py.
#- moved display from write/write.py to plotting/mrrplot.py
#- rewrote plot()
#- minor changes to display():
#   - added colorbar
#   - turned on minor tics
#   - added zoom based on masked


from .mrrplot import *

__all__ = ['plot', 'display', 'display_plain', 'overview']
