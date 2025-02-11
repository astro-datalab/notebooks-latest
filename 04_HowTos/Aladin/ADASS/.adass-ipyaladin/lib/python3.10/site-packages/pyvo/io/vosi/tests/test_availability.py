#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.io.vosi
"""
import pyvo.io.vosi as vosi

from astropy.utils.data import get_pkg_data_filename


class TestAvailability:
    def test_availability(self):
        availability = vosi.parse_availability(get_pkg_data_filename(
            "data/availability.xml"))

        assert availability.available
        assert availability.upsince == "2000-00-00T00:00:00Z"
        assert availability.downat == "2666-00-00T00:00:00Z"
        assert availability.backat == "2666-23-23T13:37:00Z"
        assert "foo" in availability.notes
        assert "bar" in availability.notes
