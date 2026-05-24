# vim: set ft=python tw=88 nu ai et ts=4 sw=4:
# ==============================================================================
#  Trace function
#  Python equivalent of trace.R (lecturer's helper)
# ==============================================================================

import numpy as np


def trace(data):
    """Return the trace (sum of diagonal elements) of a matrix."""
    return float(np.trace(data))
