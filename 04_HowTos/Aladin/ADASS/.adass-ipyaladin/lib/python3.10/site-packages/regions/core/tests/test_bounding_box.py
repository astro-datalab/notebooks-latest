# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the bounding_box module.
"""

import pytest
from numpy.testing import assert_allclose

from regions._utils.optional_deps import HAS_MATPLOTLIB
from regions.core.bounding_box import RegionBoundingBox
from regions.core.pixcoord import PixCoord
from regions.shapes import RectanglePixelRegion


def test_bounding_box_init():
    bbox = RegionBoundingBox(1, 10, 2, 20)
    assert bbox.ixmin == 1
    assert bbox.ixmax == 10
    assert bbox.iymin == 2
    assert bbox.iymax == 20


def test_bounding_box_init_minmax():
    with pytest.raises(ValueError):
        RegionBoundingBox(100, 1, 1, 100)
    with pytest.raises(ValueError):
        RegionBoundingBox(1, 100, 100, 1)


def test_bounding_box_inputs():
    with pytest.raises(TypeError):
        RegionBoundingBox([1], [10], [2], [9])
    with pytest.raises(TypeError):
        RegionBoundingBox([1, 2], 10, 2, 9)
    with pytest.raises(TypeError):
        RegionBoundingBox(1.0, 10.0, 2.0, 9.0)
    with pytest.raises(TypeError):
        RegionBoundingBox(1.3, 10, 2, 9)
    with pytest.raises(TypeError):
        RegionBoundingBox(1, 10.3, 2, 9)
    with pytest.raises(TypeError):
        RegionBoundingBox(1, 10, 2.3, 9)
    with pytest.raises(TypeError):
        RegionBoundingBox(1, 10, 2, 9.3)


def test_bounding_box_from_float():
    # This is the example from the method docstring
    bbox = RegionBoundingBox.from_float(xmin=1.0, xmax=10.0, ymin=2.0,
                                        ymax=20.0)
    assert bbox == RegionBoundingBox(ixmin=1, ixmax=11, iymin=2, iymax=21)

    bbox = RegionBoundingBox.from_float(xmin=1.4, xmax=10.4, ymin=1.6,
                                        ymax=10.6)
    assert bbox == RegionBoundingBox(ixmin=1, ixmax=11, iymin=2, iymax=12)


def test_bounding_box_eq():
    bbox = RegionBoundingBox(1, 10, 2, 20)
    assert bbox == RegionBoundingBox(1, 10, 2, 20)
    assert bbox != RegionBoundingBox(9, 10, 2, 20)
    assert bbox != RegionBoundingBox(1, 99, 2, 20)
    assert bbox != RegionBoundingBox(1, 10, 9, 20)
    assert bbox != RegionBoundingBox(1, 10, 2, 99)

    with pytest.raises(TypeError):
        assert bbox == (1, 10, 2, 20)


def test_bounding_box_repr():
    bbox = RegionBoundingBox(1, 10, 2, 20)
    assert repr(bbox) == ('RegionBoundingBox(ixmin=1, ixmax=10, iymin=2, '
                          'iymax=20)')


def test_bounding_box_shape():
    bbox = RegionBoundingBox(1, 10, 2, 20)
    assert bbox.shape == (18, 9)


def test_bounding_box_center():
    bbox = RegionBoundingBox(1, 10, 2, 20)
    assert bbox.center == (10.5, 5)


def test_bounding_box_get_overlap_slices():
    bbox = RegionBoundingBox(1, 10, 2, 20)
    slc = ((slice(2, 20, None), slice(1, 10, None)),
           (slice(0, 18, None), slice(0, 9, None)))
    assert bbox.get_overlap_slices((50, 50)) == slc

    bbox = RegionBoundingBox(-10, -1, 2, 20)
    assert bbox.get_overlap_slices((50, 50)) == (None, None)

    bbox = RegionBoundingBox(-10, 10, -10, 20)
    slc = ((slice(0, 20, None), slice(0, 10, None)),
           (slice(10, 30, None), slice(10, 20, None)))
    assert bbox.get_overlap_slices((50, 50)) == slc


def test_bounding_box_extent():
    bbox = RegionBoundingBox(1, 10, 2, 20)
    assert_allclose(bbox.extent, (0.5, 9.5, 1.5, 19.5))


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
def test_bounding_box_as_artist():
    bbox = RegionBoundingBox(1, 10, 2, 20)
    patch = bbox.as_artist()

    assert_allclose(patch.get_xy(), (0.5, 1.5))
    assert_allclose(patch.get_width(), 9)
    assert_allclose(patch.get_height(), 18)


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
def test_bounding_box_plot():
    from matplotlib.patches import Patch
    bbox = RegionBoundingBox(1, 10, 2, 20)
    patch = bbox.plot()
    assert isinstance(patch, Patch)


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
def test_bounding_box_to_region():
    bbox = RegionBoundingBox(1, 10, 2, 20)
    region = RectanglePixelRegion(PixCoord(5.0, 10.5), width=9., height=18.)
    bbox_region = bbox.to_region()
    assert_allclose(bbox_region.center.xy, region.center.xy)
    assert bbox_region.width == region.width
    assert bbox_region.height == region.height
    assert bbox_region.angle == region.angle


def test_bounding_box_union():
    bbox1 = RegionBoundingBox(1, 10, 2, 20)
    bbox2 = RegionBoundingBox(5, 21, 7, 32)
    bbox_union_expected = RegionBoundingBox(1, 21, 2, 32)
    bbox_union1 = bbox1 | bbox2
    bbox_union2 = bbox1.union(bbox2)

    assert bbox_union1 == bbox_union_expected
    assert bbox_union1 == bbox_union2

    with pytest.raises(TypeError):
        bbox1.union((5, 21, 7, 32))


def test_bounding_box_intersect():
    bbox1 = RegionBoundingBox(1, 10, 2, 20)
    bbox2 = RegionBoundingBox(5, 21, 7, 32)
    bbox_intersect_expected = RegionBoundingBox(5, 10, 7, 20)
    bbox_intersect1 = bbox1 & bbox2
    bbox_intersect2 = bbox1.intersection(bbox2)

    assert bbox_intersect1 == bbox_intersect_expected
    assert bbox_intersect1 == bbox_intersect2

    with pytest.raises(TypeError):
        bbox1.intersection((5, 21, 7, 32))

    assert bbox1.intersection(RegionBoundingBox(30, 40, 50, 60)) is None
