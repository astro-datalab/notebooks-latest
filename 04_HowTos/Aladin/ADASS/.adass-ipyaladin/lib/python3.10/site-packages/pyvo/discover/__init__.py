# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Various functions dealing with global data disovery.
"""

import warnings
from pyvo.utils.prototype import PrototypeWarning

# if you remove this warning, also remove the ignorere in test_imagediscovery.
warnings.warn("pyvo.discover's API is still under design in pyVO 1.6 and"
    " may change without prior notice.  Feedback to the authors is most"
    " welcome.", PrototypeWarning)

from .image import images_globally, ImageDiscoverer, ImageFound

__all__ = ['images_globally', "ImageDiscoverer", "ImageFound"]
