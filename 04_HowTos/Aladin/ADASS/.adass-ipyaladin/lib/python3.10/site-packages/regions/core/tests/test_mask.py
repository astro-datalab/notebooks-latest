# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the mask module.
"""

import astropy.units as u
import numpy as np
import pytest
from astropy.utils import minversion
from numpy.testing import assert_allclose, assert_almost_equal

from regions.core.bounding_box import RegionBoundingBox
from regions.core.mask import RegionMask
from regions.core.pixcoord import PixCoord
from regions.shapes import CircleAnnulusPixelRegion, CirclePixelRegion

NUMPY_LT_2_0 = not minversion(np, '2.0.dev')
COPY_IF_NEEDED = False if NUMPY_LT_2_0 else None
POSITIONS = [(-20, -20), (-20, 20), (20, -20), (60, 60)]


def test_mask_input_shapes():
    with pytest.raises(ValueError):
        mask_data = np.ones((10, 10))
        bbox = RegionBoundingBox(5, 10, 5, 10)
        RegionMask(mask_data, bbox)


def test_mask_array():
    mask_data = np.ones((10, 10))
    bbox = RegionBoundingBox(5, 15, 5, 15)
    mask = RegionMask(mask_data, bbox)
    data = np.array(mask)
    assert_allclose(data, mask.data)


def test_mask_copy():
    bbox = RegionBoundingBox(5, 15, 5, 15)

    mask = RegionMask(np.ones((10, 10)), bbox)
    mask_copy = np.array(mask, copy=True)
    mask_copy[0, 0] = 100.0
    assert mask.data[0, 0] == 1.0

    mask = RegionMask(np.ones((10, 10)), bbox)
    mask_copy = np.array(mask, copy=False)
    mask_copy[0, 0] = 100.0
    assert mask.data[0, 0] == 100.0

    # no copy; copy=None returns a copy only if __array__ returns a copy
    # copy=None was introduced in NumPy 2.0
    mask = RegionMask(np.ones((10, 10)), bbox)
    mask_copy = np.array(mask, copy=COPY_IF_NEEDED)
    mask_copy[0, 0] = 100.0
    assert mask.data[0, 0] == 100.0

    # no copy
    mask = RegionMask(np.ones((10, 10)), bbox)
    mask_copy = np.asarray(mask)
    mask_copy[0, 0] = 100.0
    assert mask.data[0, 0] == 100.0

    # needs to copy because of the dtype change
    mask = RegionMask(np.ones((10, 10)), bbox)
    mask_copy = np.asarray(mask, dtype=int)
    mask_copy[0, 0] = 100.0
    assert mask.data[0, 0] == 1.0


def test_mask_get_overlap_slices():
    aper = CirclePixelRegion(PixCoord(5, 5), radius=10.)
    mask = aper.to_mask()
    slc = ((slice(0, 16, None), slice(0, 16, None)),
           (slice(5, 21, None), slice(5, 21, None)))
    assert mask.get_overlap_slices((25, 25)) == slc


def test_mask_cutout_shape():
    mask_data = np.ones((10, 10))
    bbox = RegionBoundingBox(5, 15, 5, 15)
    mask = RegionMask(mask_data, bbox)

    with pytest.raises(ValueError):
        mask.cutout(np.arange(10))

    with pytest.raises(ValueError):
        mask.to_image((10,))


def test_mask_cutout_copy():
    data = np.ones((50, 50))
    aper = CirclePixelRegion(PixCoord(25, 25), radius=10.)
    mask = aper.to_mask()
    cutout = mask.cutout(data, copy=True)
    data[25, 25] = 100.
    assert cutout[10, 10] == 1.

    # test quantity data
    data2 = np.ones((50, 50)) * u.adu
    cutout2 = mask.cutout(data2, copy=True)
    assert cutout2.unit == data2.unit
    data2[25, 25] = 100. * u.adu
    assert cutout2[10, 10].value == 1.


@pytest.mark.parametrize('position', POSITIONS)
def test_mask_cutout_no_overlap(position):
    data = np.ones((50, 50))
    aper = CirclePixelRegion(PixCoord(position[0], position[1]), radius=10.)
    mask = aper.to_mask()

    cutout = mask.cutout(data)
    assert cutout is None

    weighted_data = mask.multiply(data)
    assert weighted_data is None

    image = mask.to_image(data.shape)
    assert image is None


@pytest.mark.parametrize('position', POSITIONS)
def test_mask_cutout_partial_overlap(position):
    data = np.ones((50, 50))
    aper = CirclePixelRegion(PixCoord(position[0], position[1]), radius=30.)
    mask = aper.to_mask()

    cutout = mask.cutout(data)
    assert cutout.shape == mask.shape

    weighted_data = mask.multiply(data)
    assert weighted_data.shape == mask.shape

    image = mask.to_image(data.shape)
    assert image.shape == data.shape


def test_mask_multiply():
    radius = 10.
    data = np.ones((50, 50))
    region = CirclePixelRegion(PixCoord(25, 25), radius=radius)
    mask = region.to_mask(mode='exact')
    data_weighted = mask.multiply(data)
    assert_almost_equal(np.sum(data_weighted), np.pi * radius**2)

    # test that multiply() returns a copy
    data[25, 25] = 100.
    assert data_weighted[10, 10] == 1.


def test_mask_multiply_quantity():
    radius = 10.
    data = np.ones((50, 50)) * u.adu
    region = CirclePixelRegion(PixCoord(25, 25), radius=radius)
    mask = region.to_mask(mode='exact')
    data_weighted = mask.multiply(data)
    assert data_weighted.unit == u.adu
    assert_almost_equal(np.sum(data_weighted.value), np.pi * radius**2)

    # test that multiply() returns a copy
    data[25, 25] = 100. * u.adu
    assert data_weighted[10, 10].value == 1.


@pytest.mark.parametrize('value', (np.nan, np.inf))
def test_mask_nonfinite_fill_value(value):
    region = CircleAnnulusPixelRegion(PixCoord(0, 0), 10, 20)
    data = np.ones((101, 101)).astype(int)
    cutout = region.to_mask().cutout(data, fill_value=value)
    assert ~np.isfinite(cutout[0, 0])


def test_mask_multiply_fill_value():
    region = CircleAnnulusPixelRegion(PixCoord(0, 0), 10, 20)
    data = np.ones((101, 101)).astype(int)
    cutout = region.to_mask().multiply(data, fill_value=np.nan)
    xypos = ((20, 20), (5, 5), (5, 35), (35, 5), (35, 35))
    for x, y in xypos:
        assert np.isnan(cutout[y, x])


def test_mask_nonfinite_in_bbox():
    """
    Regression test that non-finite data values outside of the mask but
    within the bounding box are set to zero.
    """
    data = np.ones((101, 101))
    data[33, 33] = np.nan
    data[67, 67] = np.inf
    data[33, 67] = -np.inf
    data[22, 22] = np.nan
    data[22, 23] = np.inf

    radius = 20.
    reg1 = CirclePixelRegion(PixCoord(50, 50), radius)
    reg2 = CirclePixelRegion(PixCoord(5, 5), radius)

    wdata1 = reg1.to_mask(mode='exact').multiply(data)
    assert_allclose(np.sum(wdata1), np.pi * radius**2)

    wdata2 = reg2.to_mask(mode='exact').multiply(data)
    assert_allclose(np.sum(wdata2), 561.6040111923013)


@pytest.mark.parametrize('x, y, exp_shape',
                         [(0, 0, 245), (50, 50, 940), (100, 100, 245)])
def test_mask_get_values(x, y, exp_shape):
    aper = CircleAnnulusPixelRegion(PixCoord(x, y), inner_radius=10,
                                    outer_radius=20)
    data = np.ones((101, 101))
    values = aper.to_mask(mode='center').get_values(data)
    assert values.shape == (exp_shape,)
    assert_allclose(np.sum(values), exp_shape)


def test_mask_get_values_no_overlap():
    aper = CirclePixelRegion(PixCoord(-100, -100), radius=3)
    data = np.ones((51, 51))
    values = aper.to_mask().get_values(data)
    assert values.shape == (0,)


def test_mask_get_values_mask():
    aper = CirclePixelRegion(PixCoord(24.5, 24.5), radius=10.)
    data = np.ones((51, 51))
    mask = aper.to_mask(mode='exact')
    with pytest.raises(ValueError):
        mask.get_values(data, mask=np.ones(3))

    arr = mask.get_values(data, mask=None)
    assert_allclose(np.sum(arr), 100. * np.pi)

    data_mask = np.zeros(data.shape, dtype=bool)
    data_mask[25:] = True
    arr2 = mask.get_values(data, mask=data_mask)
    assert_allclose(np.sum(arr2), 100. * np.pi / 2.)
