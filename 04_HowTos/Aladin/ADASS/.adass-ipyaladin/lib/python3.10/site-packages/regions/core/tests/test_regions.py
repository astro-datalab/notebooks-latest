# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the regions module.
"""

import pytest
from astropy.table import Table

from regions.core import PixCoord, Regions
from regions.shapes import CirclePixelRegion


def test_regions_inputs():
    regs = []
    for radius in range(1, 5):
        center = PixCoord(14, 21)
        regs.append(CirclePixelRegion(center, radius=radius))

    regions = Regions(regs)
    assert len(regions) == 4


def test_regions_no_input():
    regions = Regions()
    assert len(regions) == 0

    regions = Regions([])
    assert len(regions) == 0


def test_regions_invalid_input():
    match = "'int' object is not iterable"
    with pytest.raises(TypeError, match=match):
        Regions(1)

    match = 'Input regions must be a list of Region objects'
    with pytest.raises(TypeError, match=match):
        Regions([1])
    with pytest.raises(TypeError, match=match):
        Regions([1, 2, 3])


def test_regions_append():
    regions = Regions()

    for radius in range(1, 5):
        center = PixCoord(14, 21)
        regions.append(CirclePixelRegion(center, radius=radius))
    assert len(regions) == 4

    match = 'Input region must be a Region object'
    with pytest.raises(TypeError, match=match):
        regions.append(1)


def test_regions_extend():
    regions = Regions()

    regs = []
    for radius in range(1, 5):
        center = PixCoord(14, 21)
        regs.append(CirclePixelRegion(center, radius=radius))

    regions.extend(regs)
    assert len(regions) == 4

    regions.extend(regions)
    assert len(regions) == 8

    match = 'Input regions must be a list of Region objects'
    with pytest.raises(TypeError, match=match):
        regions.extend([1])


def test_regions_methods():
    regions = Regions()

    for radius in range(1, 5):
        center = PixCoord(14, 21)
        regions.append(CirclePixelRegion(center, radius=radius))

    reg_slc = regions[0:2]
    assert len(reg_slc) == 2

    reg = CirclePixelRegion(PixCoord(0, 0), radius=1)
    regions.insert(0, reg)
    assert len(regions) == 5
    assert regions[0] == reg

    regions2 = regions.copy()
    regions.reverse()
    assert regions.regions == regions2.regions[::-1]

    outreg = regions.pop(-1)
    assert outreg == reg


def test_regions_get_formats():
    reg = CirclePixelRegion(PixCoord(0, 0), radius=1)
    regions = Regions([reg])
    tbl = regions.get_formats()
    assert isinstance(tbl, Table)
    assert len(tbl) == 3


def test_regions_repr():
    reg = CirclePixelRegion(PixCoord(0, 0), radius=1)
    regions = Regions([reg])
    regions_str = '[<CirclePixelRegion(center=PixCoord(x=0, y=0), radius=1)>]'
    regions_repr = f'<Regions({regions_str})>'
    assert str(regions) == regions_str
    assert repr(regions) == regions_repr
