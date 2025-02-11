# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the compound module.
"""

import operator

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from numpy.testing import assert_allclose

from regions.core import CompoundPixelRegion, PixCoord, RegionBoundingBox
from regions.shapes import CirclePixelRegion, CircleSkyRegion
from regions.tests.helpers import make_simple_wcs


class TestCompoundPixel:
    """
    Test compound pixel regions.
    """

    # Two circles that overlap in one column
    c1 = CirclePixelRegion(PixCoord(5, 5), 4)
    c2 = CirclePixelRegion(PixCoord(11, 5), 4)

    def test_copy(self):
        reg = (self.c1 | self.c2).copy()
        assert isinstance(reg, CompoundPixelRegion)
        assert reg.operator is operator.or_

        assert isinstance(reg.region1, CirclePixelRegion)
        assert reg.region1.center.xy == (5, 5)
        assert reg.region1.radius == 4
        assert reg.region1.meta == {}
        assert reg.region1.visual == {}

        assert isinstance(reg.region2, CirclePixelRegion)
        assert reg.region2.center.xy == (11, 5)
        assert reg.region2.radius == 4
        assert reg.region2.meta == {}
        assert reg.region2.visual == {}

    def test_rotate(self):
        reg = (self.c1 | self.c2).rotate(PixCoord(0, 0), 360 * u.deg)
        assert isinstance(reg, CompoundPixelRegion)
        assert reg.operator is operator.or_

        assert isinstance(reg.region1, CirclePixelRegion)
        assert_allclose(reg.region1.center.xy, (5, 5))
        assert reg.region1.radius == 4
        assert reg.region1.meta == {}
        assert reg.region1.visual == {}

        assert isinstance(reg.region2, CirclePixelRegion)
        assert_allclose(reg.region2.center.xy, (11, 5))
        assert reg.region2.radius == 4
        assert reg.region2.meta == {}
        assert reg.region2.visual == {}

    def test_union(self):
        union = self.c1 | self.c2
        assert isinstance(union, CompoundPixelRegion)
        mask = union.to_mask()
        assert_allclose(mask.data[:, 7], [0, 0, 1, 1, 1, 1, 1, 0, 0])
        assert_allclose(mask.data[:, 6], [0, 1, 1, 1, 1, 1, 1, 1, 0])

    def test_intersection(self):
        intersection = self.c1 & self.c2
        mask = intersection.to_mask()
        assert_allclose(mask.data[:, 7], [0, 0, 1, 1, 1, 1, 1, 0, 0])
        assert_allclose(mask.data[:, 6], [0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_symdiff(self):
        symdiff = self.c1 ^ self.c2
        mask = symdiff.to_mask()
        assert_allclose(mask.data[:, 7], [0, 0, 0, 0, 0, 0, 0, 0, 0])
        assert_allclose(mask.data[:, 6], [0, 1, 1, 1, 1, 1, 1, 1, 0])

    def test_to_mask(self):
        # fixes #168
        pixcoord3 = PixCoord(1, 1)
        c3 = CirclePixelRegion(pixcoord3, 4)
        union2 = self.c1 | c3
        mask1 = union2.to_mask()

        pixcoord4 = PixCoord(9, 9)
        c4 = CirclePixelRegion(pixcoord4, 4)
        union3 = self.c1 | c4
        mask2 = union3.to_mask()

        # mask1 and mask2 should be equal
        assert_allclose(mask1.data, mask2.data)

        ref_data = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
                             [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                             [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
        assert_allclose(mask1.data, ref_data)

    def test_bounding_box(self):
        bbox = (self.c1 | self.c2).bounding_box
        assert bbox == RegionBoundingBox(1, 16, 1, 10)

    def test_invalid(self):
        with pytest.raises(ValueError):
            self.c1 | 1.0
        with pytest.raises(ValueError):
            self.c1 | PixCoord(1.0, 2.0)


def test_compound_sky():
    skycoord1 = SkyCoord(0 * u.deg, 0 * u.deg, frame='galactic')
    c1 = CircleSkyRegion(skycoord1, 1 * u.deg)

    skycoord2 = SkyCoord(1 * u.deg, 1 * u.deg, frame='galactic')
    c2 = CircleSkyRegion(skycoord2, 0.5 * u.deg)

    test_coord1 = SkyCoord(1.2 * u.deg, 1.2 * u.deg, frame='galactic')
    test_coord2 = SkyCoord(0.5 * u.deg, 0.5 * u.deg, frame='galactic')
    test_coord3 = SkyCoord(0.7 * u.deg, 0.7 * u.deg, frame='galactic')
    test_coord4 = SkyCoord(2 * u.deg, 5 * u.deg, frame='galactic')

    wcs = make_simple_wcs(skycoord1, 0.1 * u.deg, 20)

    assert c2.contains(test_coord1, wcs) and not c1.contains(test_coord1, wcs)
    assert not c2.contains(test_coord2, wcs) and c1.contains(test_coord2, wcs)
    assert c1.contains(test_coord3, wcs) and c2.contains(test_coord3, wcs)
    assert (not c2.contains(test_coord4, wcs)
            and not c1.contains(test_coord4, wcs))

    coords = SkyCoord([test_coord1, test_coord2, test_coord3, test_coord4],
                      frame='galactic')

    union = c1 | c2
    assert (union.contains(coords, wcs) == [True, True, True, False]).all()

    intersection = c1 & c2
    assert ((intersection.contains(coords, wcs)
             == [False, False, True, False]).all())

    diff = c1 ^ c2
    assert (diff.contains(coords, wcs) == [True, True, False, False]).all()

    c3 = CircleSkyRegion(test_coord4, 0.1 * u.deg)
    union = c1 | c2 | c3
    assert (union.contains(coords, wcs) == [True, True, True, True]).all()

    intersection = c1 & c2 & c3
    assert ((intersection.contains(coords, wcs)
             == [False, False, False, False]).all())

    diff = c1 ^ c2 ^ c3
    assert (diff.contains(coords, wcs) == [True, True, False, True]).all()
    assert 'Compound' in str(union)
