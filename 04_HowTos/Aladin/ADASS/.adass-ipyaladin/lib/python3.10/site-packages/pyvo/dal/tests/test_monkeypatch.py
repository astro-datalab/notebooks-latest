#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.datalink
"""
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename


def test_monkeypatch():
    Table.read(get_pkg_data_filename("data/monkeypatch.xml"))
    import pyvo  # noqa: F401
    Table.read(get_pkg_data_filename("data/monkeypatch.xml"))
