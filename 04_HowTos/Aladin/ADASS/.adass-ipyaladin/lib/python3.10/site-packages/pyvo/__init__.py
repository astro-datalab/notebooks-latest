# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
PyVO is a package providing access to remote data and services of the
Virtual observatory (VO) using Python.

The pyvo module currently provides these main capabilities:

* find archives that provide particular data of a particular type and/or
  relates to a particular topic

  *  regsearch()

* search an archive for datasets of a particular type

  *  imagesearch(), spectrumsearch()

* do simple searches on catalogs or databases

  *  conesearch(), linesearch(), tablesearch()

Submodules provide additional functions and classes for greater control over
access to these services.

This module also exposes the exception classes raised by the above functions,
of which DALAccessError is the root parent exception.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

from . import registry
from .dal import ssa, sia, sla, scs, tap
from . import auth
from .registry import search as regsearch
from .dal import (
    imagesearch, spectrumsearch, conesearch, linesearch, tablesearch,
    DALAccessError, DALProtocolError, DALFormatError, DALServiceError,
    DALQueryError)
