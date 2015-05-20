
__version__ = '1.4.3'
# $Source$

#Version history
#==================
#
#changes in 1.0
#----------------
#- added mrr.read.read_mask to be used for mask-read-in
#
#- switch to mrr.write version 1.1
#   - rewrote plot() and moved from mrr.write to mrr.plotting
#   - changed display() and moved from mrr.write to mrr.plotting
#
#- added mrr.plotting
#
#
#changes in 1.2
#----------------
#- added mrr.waveform version 1.0
#
#- switched to mrr.plot version 1.1
#   - modified display():
#       - fixed cropping
#       - added keyword hold to keep figure alive
#       - added support for keyargs passed through to matplotlib.imshow
#
#- switched to mrr.write version 1.2
#   - added adjust_window()
#       - added window_center(), window_zero() and set_mask()
#
#- switched to mrr.core version 1.1
#   - added min(), max(), zeros(), zeros_like(), load() and save()
#   - minor changes in MRRArray
#
#- added mrr.timeline version 0.8
#
#
#changes in 1.2.1
#----------------
#- switched to mrr.plot version 1.2
#- switched to mrr.write version 1.3
#- added arc version 0.9
#   - added arc to import
#
#
#changes in 1.2.2
#----------------
#- switched to mrr.plot version 1.3
#- switched to mrr.read version 1.1
#
#changes in 1.3.0
#----------------
#- switched to relative imports
#- switched to mrr.arc version 1.1
#- switched to mrr.core version 1.1.1
#- switched to mrr.plotting version 1.3.1
#   - bugfix in display_plain() (ValueError when crop=False)
#   - added return_index to plot()
#- switched to mrr.read version 1.1.1
#- switched to mrr.write version 1.3.1
#- removed unwrap
#- added unwrapping version 1.2
#
#changes in 1.3.1
#----------------
#- switched to unwrapping version 1.3
#- added direct import of unwrap_array
#
#changes in 1.3.2
#----------------
#- switched to read version 1.2
#- switched to unwrapping version 1.3
#
#changes in 1.4
#--------------
#- switched to unwrapping version 1.4
#- switched to core version 1.2
#- switched to plotting version 1.4
#- switched to waveform version 1.2
#
#changes 1.4.1
#-------------
#- removed bugs from arc.arcmr
#- various bugfixes
#- switched to timeline v0.9

import os.path

_ppath = os.path.dirname(os.path.realpath(__file__))
_testpath = os.path.join(_ppath, './test/')


from .mrrcore import *
from .read import *
from .write import *
from .plotting import *
from .unwrapping import *
from .timeline import normalize_image_set, broaden_mask
from . import waveform, unwrapping, bvalue, arc, coordinates, timeline


def return_test_path():
    return os.path.normpath(_testpath)
