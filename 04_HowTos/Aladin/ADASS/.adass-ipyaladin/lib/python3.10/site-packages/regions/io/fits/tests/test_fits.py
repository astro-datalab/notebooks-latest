# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the fits subpackage.
"""

import warnings

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.utils.data import get_pkg_data_filenames
from astropy.utils.exceptions import AstropyUserWarning
from numpy.testing import assert_equal

from regions.core import PixCoord, RegionMeta, Regions
from regions.io.fits.core import FITSParserError
from regions.shapes import (CirclePixelRegion, CircleSkyRegion,
                            LinePixelRegion, RectangleAnnulusPixelRegion,
                            TextPixelRegion)
from regions.tests.helpers import assert_region_allclose


def test_roundtrip(tmpdir):
    filenames = get_pkg_data_filenames('data', pattern='*.fits')

    # AstropyUserWarning will be emitted only for some of the files
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyUserWarning)
        # Check that all test files are readable
        for filename in filenames:
            regions = Regions.read(filename, format='fits')
            assert len(regions) > 0

        tempfile = tmpdir.join('tmp.fits').strpath
        regions.write(tempfile, format='fits', overwrite=True)
        regions2 = Regions.read(tempfile, format='fits')
        for reg1, reg2 in zip(regions, regions2, strict=True):
            assert_region_allclose(reg1, reg2)


def test_only_pixel_regions():
    reg_sky = CircleSkyRegion(SkyCoord(1, 2, unit='deg'), 5 * u.deg)
    reg = Regions([reg_sky])

    match = 'Sky regions cannot be serialized'
    with pytest.warns(AstropyUserWarning, match=match):
        result = reg.serialize(format='fits')
        assert len(result) == 0

    reg_pix = CirclePixelRegion(PixCoord(10, 10), 5)
    reg = Regions([reg_sky, reg_pix])
    match = 'Sky regions cannot be serialized'
    with pytest.warns(AstropyUserWarning, match=match):
        result = reg.serialize(format='fits')
        assert len(result) == 1


def test_valid_columns():
    t = QTable([[1, 2, 3]], names=('a'))
    with pytest.raises(FITSParserError) as excinfo:
        Regions.parse(t, format='fits')
        estr = "'a' is not a valid column name"
        assert estr in str(excinfo.value)


def test_invalid_regions():
    center = PixCoord(42, 43)
    region1 = RectangleAnnulusPixelRegion(center, 4.2, 5.2, 7.2, 8.2,
                                          angle=15 * u.deg)
    region2 = LinePixelRegion(start=center, end=PixCoord(x=52, y=53))
    region3 = TextPixelRegion(center, text='Example Text')
    regions = (region1, region2, region3)

    match = 'cannot be serialized using the FITS format'
    for reg in regions:
        with pytest.warns(AstropyUserWarning, match=match):
            _ = reg.serialize(format='fits')


def test_valid_row():
    x = [1] * u.pix
    y = [2] * u.pix
    shapes = ['CIRCLE']
    tbl = QTable([x, y, shapes], names=('X', 'Y', 'SHAPE'))
    match = 'Table columns are missing'
    with pytest.warns(AstropyUserWarning, match=match):
        Regions.parse(tbl, format='fits')

    # test invalid shape
    shapes = ['INVALID']
    tbl2 = QTable([x, y, shapes], names=('X', 'Y', 'SHAPE'))
    with pytest.raises(FITSParserError) as excinfo:
        Regions.parse(tbl2, format='fits')
        estr = "'invalid' is not a valid FITS region shape"
        assert estr in str(excinfo.value)

    shapes = ['PIE']
    rotang = [[20, 30]] * u.deg
    tbl3 = QTable([x, y, shapes, rotang], names=('X', 'Y', 'SHAPE', 'ROTANG'))
    match = "'pie' is not supported"
    with pytest.warns(AstropyUserWarning, match=match):
        Regions.parse(tbl3, format='fits')


def test_components():
    center = PixCoord(10, 10)
    regions = Regions((CirclePixelRegion(center, 3),
                       CirclePixelRegion(center, 5),
                       CirclePixelRegion(center, 7),
                       CirclePixelRegion(center, 9)))
    tbl1 = regions.serialize(format='fits')
    assert 'COMPONENT' not in tbl1.colnames

    components = [None, None, None, None]
    for region, component in zip(regions, components, strict=True):
        region.meta = RegionMeta({'component': component})
    tbl2 = regions.serialize(format='fits')
    assert 'COMPONENT' not in tbl2.colnames

    components = np.arange(4)
    for region, component in zip(regions, components, strict=True):
        region.meta = RegionMeta({'component': component})
    tbl3 = regions.serialize(format='fits')
    assert 'COMPONENT' in tbl3.colnames
    assert_equal(tbl3['COMPONENT'], components)

    components = [1, 2, None, 4]
    for region, component in zip(regions, components, strict=True):
        region.meta = RegionMeta({'component': component})
    tbl4 = regions.serialize(format='fits')
    assert 'COMPONENT' in tbl4.colnames
    assert_equal(tbl4['COMPONENT'], [1, 2, 5, 4])
