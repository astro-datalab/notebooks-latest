# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .sia import search as imagesearch
from .sia2 import search as imagesearch2
from .ssa import search as spectrumsearch
from .sla import search as linesearch
from .scs import search as conesearch
from .tap import search as tablesearch

from .query import DALService, DALQuery, DALResults, Record

from .sia import SIAService, SIAQuery, SIAResults, SIARecord
from .sia2 import SIA2Service, SIA2Query, SIA2Results, ObsCoreRecord
from .ssa import SSAService, SSAQuery, SSAResults, SSARecord
from .sla import SLAService, SLAQuery, SLAResults, SLARecord
from .scs import SCSService, SCSQuery, SCSResults, SCSRecord
from .tap import TAPService, TAPQuery, TAPResults, AsyncTAPJob


from .exceptions import (
    DALAccessError, DALProtocolError, DALFormatError, DALServiceError,
    DALQueryError, DALOverflowWarning)

__all__ = [
    "imagesearch", "spectrumsearch", "linesearch", "conesearch", "tablesearch",
    "DALService", "imagesearch2",
    "SIAService", "SIA2Service", "SSAService", "SLAService", "SCSService", "TAPService",
    "DALQuery", "SIAQuery", "SIA2Query", "SSAQuery", "SLAQuery", "SCSQuery", "TAPQuery",
    "DALResults",
    "SIAResults", "SIA2Results", "SSAResults", "SLAResults", "SCSResults", "TAPResults",
    "Record", "ObsCoreRecord",
    "SIARecord", "SSARecord", "SLARecord", "SCSRecord",
    "AsyncTAPJob",
    "DALAccessError", "DALProtocolError", "DALFormatError", "DALServiceError",
    "DALQueryError", "DALOverflowWarning"]
