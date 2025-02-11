# Licensed under a 3-clause BSD style license - see LICENSE.rst
from astropy.utils.collections import HomogeneousList

from ...utils.xml.elements import xmlelement, Element
from .exceptions import W32, W33, W34, W35

__all__ = ["Availability"]


######################################################################
# FACTORY FUNCTIONS
def _convert_boolean(value, default=None):
    return {
        'false': False,
        '0': False,
        'true': True,
        '1': True
    }.get(value, default)


######################################################################
# ELEMENT CLASSES
class Availability(Element):
    def __init__(self, config=None, pos=None, _name='availability', **kwargs):
        super().__init__(config, pos, _name, **kwargs)

        self._available = None
        self._upsince = None
        self._downat = None
        self._backat = None
        self._notes = HomogeneousList(str)

    @xmlelement(plain=True, multiple_exc=W32)
    def available(self):
        return self._available

    @available.setter
    def available(self, available):
        self._available = _convert_boolean(available)

    @xmlelement(name='upSince', plain=True, multiple_exc=W33)
    def upsince(self):
        return self._upsince

    @upsince.setter
    def upsince(self, upsince):
        self._upsince = upsince

    @xmlelement(name='downAt', plain=True, multiple_exc=W34)
    def downat(self):
        return self._downat

    @downat.setter
    def downat(self, downat):
        self._downat = downat

    @xmlelement(name='backAt', plain=True, multiple_exc=W35)
    def backat(self):
        return self._backat

    @backat.setter
    def backat(self, backat):
        self._backat = backat

    @xmlelement(name='note', plain=True)
    def notes(self):
        return self._notes
