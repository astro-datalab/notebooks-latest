# Licensed under a 3-clause BSD style license - see LICENSE
"""
MOCPy is a Python library allowing easy creation and manipulation of MOCs (Multi-Order Coverage maps).

MOC is an `IVOA standard <http://ivoa.net/documents/MOC/>`__ enabling description of arbitrary sky regions.
Based on the HEALPix sky tessellation, it maps regions on the sky into hierarchically grouped predefined cells.
An experimental support for TMOC (temporal MOC) has been added since version 0.4.0.
It allows creation, parsing and comparison of TMOCs.
"""

from .fmoc import FrequencyMOC
from .moc import MOC, WCS
from .stmoc import STMOC
from .tmoc import TimeMOC
from .version import __version__

__all__ = [
    "MOC",
    "TimeMOC",
    "FrequencyMOC",
    "STMOC",
    "WCS",
    "__version__",
]
