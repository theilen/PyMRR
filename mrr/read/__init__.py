

__version__ = '1.31'
# $Source$

# version history
# =================
# version 1.3
# -----------
# - switched to read version 1.3
# - added parse_dicom v0.8
# - added siemens_csa v0.8
# - switched to explicit import from read
#
# version 1.2
# - switched to read version 1.2

# version 1.1.1
# -------------
# - switched to relative imports
# - added support for different unwrap-algorithms to read_dicom and
#    read-dicom_set to support unwrapping version 1.2
#
# version 1.1
# -----------
# - read_dicom_set: added IOError if nameparser returns empty list
# - read_mask: added a test for invertet masks
# - nameparser: added ValueError if filename does not match pattern
#       (len(items)<2)

from .read import read_dicom, read_dicom_set, read_bitmap, read_mask
from .parse_dicom import parse_parameters, variable_ptft

__all__ = ['read_dicom', 'read_dicom_set', 'read_bitmap', 'read_mask',
           "parse_parameters", "variable_ptft"
           ]
