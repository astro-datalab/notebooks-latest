# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
a package for interacting with registries.

The regtap module supports access to the IVOA Registries
"""

from .regtap import (search, ivoid2service,
    get_RegTAP_query,
    choose_RegTAP_service,
    RegistryResults, RegistryResource)

from .rtcons import (Constraint, SubqueriedConstraint,
                     Freetext, Author, Servicetype, Waveband, Datamodel, Ivoid,
                     UCD, Spatial, Spectral, Temporal, RegTAPFeatureMissing)

__all__ = ["search", "get_RegTAP_query", "Constraint", "SubqueriedConstraint",
           "Freetext", "Author",
           "Servicetype", "Waveband", "Datamodel", "Ivoid", "UCD",
           "Spatial", "Spectral", "Temporal",
           "choose_RegTAP_service", "RegTAPFeatureMissing",
           "RegistryResults", "RegistryResource",]
