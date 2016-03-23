__version__ = '1.2.1'
# $Source$


# Version history
# ==================
# chagnes in 1.2.1
# - switched to cor v1.2.1
#
# changes in 1.2
# - swtiched to core v1.2
# - switched to arithmetics v1.1
# - removed statistics (merged mean() into core)
#
# changes in 1.1.
# - switched to relative imports
#
# changes in 1.1
# ----------------
# - switched to core.py version 1.1
# - added mrr_min and mrr_max as min and max to __all__
# - added load, save, zeros, zeros_like to __all__
#


from .core import *
from .arithmetics import *
from .boolean import *
# from .statistics import *


min = mrr_min
max = mrr_max
mean = mrr_mean

__all__ = ["MRRArray", "empty", "empty_like", "zeros", "zeros_like",
           "min", "max", "load", "save", "mean", "mean_phasor",
           "cond_print"]  # core

__all__.extend(["add", "div", "divide", "mult", "multiply", "sub", "subtract",
                "power"])  # arithmetics

__all__.extend(["all_equal", "any_equal"])  # boolean
