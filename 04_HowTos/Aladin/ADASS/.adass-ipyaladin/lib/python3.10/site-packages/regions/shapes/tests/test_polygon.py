# Licensed under a 3-clause BSD style license - see LICENSE.rst

from copy import copy

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose, assert_equal

from regions._utils.examples import make_example_dataset
from regions._utils.optional_deps import HAS_MATPLOTLIB
from regions.core import PixCoord, RegionBoundingBox, RegionMeta, RegionVisual
from regions.shapes.polygon import (PolygonPixelRegion, PolygonSkyRegion,
                                    RegularPolygonPixelRegion)
from regions.shapes.tests.test_common import (BaseTestPixelRegion,
                                              BaseTestSkyRegion)
from regions.tests.helpers import make_simple_wcs


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    config = dict(crpix=(18, 9), cdelt=(-10, 10), shape=(18, 36))
    dataset = make_example_dataset(config=config)
    return dataset.wcs


class TestPolygonPixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = PolygonPixelRegion(PixCoord([1, 3, 1], [1, 1, 4]), meta=meta,
                             visual=visual)
    sample_box = [0, 4, 0, 5]
    inside = [(2, 2)]
    outside = [(3, 2), (3, 3)]
    expected_area = 3
    expected_repr = ('<PolygonPixelRegion(vertices=PixCoord(x=[1 3 1], '
                     'y=[1 1 4]))>')
    expected_str = ('Region: PolygonPixelRegion\nvertices: '
                    'PixCoord(x=[1 3 1], y=[1 1 4])')

    # We will be using this polygon for basic tests:
    #
    # (3,0) *
    #       |\
    #       | \
    #       |  \
    # (0,0) *---* (2, 0)
    #

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.vertices.x, [1, 3, 1])
        assert_allclose(reg.vertices.y, [1, 1, 4])
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.vertices.x, self.reg.vertices.x)
        assert_allclose(reg_new.vertices.y, self.reg.vertices.y)
        assert reg_new.meta == self.reg.meta
        assert reg_new.visual == self.reg.visual

        # test that converted region meta and visual are copies and not views
        reg_new.meta['text'] = 'new'
        reg_new.visual['color'] = 'green'
        assert reg_new.meta['text'] != self.reg.meta['text']
        assert reg_new.visual['color'] != self.reg.visual['color']

    def test_bounding_box(self):
        bbox = self.reg.bounding_box
        assert bbox == RegionBoundingBox(ixmin=1, ixmax=4, iymin=1, iymax=5)

    def test_to_mask(self):
        # The true area of this polygon is 3
        # We only do very low-precision asserts below,
        # because results can be unstable with points
        # on the edge of the polygon.
        # Basically we check that mask.data is filled
        # with something useful at all.

        # Bounding box and output shape is independent of subpixels,
        # so we only assert on it once here, not in the other cases below
        mask = self.reg.to_mask(mode='center', subpixels=1)
        assert 2 <= np.sum(mask.data) <= 6
        assert mask.bbox == RegionBoundingBox(ixmin=1, ixmax=4, iymin=1,
                                              iymax=5)
        assert mask.data.shape == (4, 3)

        # Test more cases for to_mask
        # This example is with the default: subpixels=5
        mask = self.reg.to_mask(mode='subpixels')
        assert 2 <= np.sum(mask.data) <= 6

        mask = self.reg.to_mask(mode='subpixels', subpixels=8)
        assert 2 <= np.sum(mask.data) <= 6

        mask = self.reg.to_mask(mode='subpixels', subpixels=9)
        assert 2 <= np.sum(mask.data) <= 6

        mask = self.reg.to_mask(mode='subpixels', subpixels=10)
        assert 2 <= np.sum(mask.data) <= 6

        with pytest.raises(NotImplementedError):
            self.reg.to_mask(mode='exact')

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
    def test_as_artist(self):
        patch = self.reg.as_artist()
        expected = [[1, 1], [3, 1], [1, 4], [1, 1]]
        assert_allclose(patch.xy, expected)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(3, 1), -90 * u.deg)
        assert_allclose(reg.vertices.x, [3, 3, 6])
        assert_allclose(reg.vertices.y, [3, 1, 3])

    def test_origin(self):
        verts = PixCoord([1, 3, 1], [1, 1, 4])
        reg1 = PolygonPixelRegion(verts)

        origin = PixCoord(1, 1)
        relverts = verts - origin
        reg2 = PolygonPixelRegion(relverts, origin=origin)
        assert_equal(reg1.vertices, reg2.vertices)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.vertices = PixCoord([1, 3, 1], [1, 1, 6])
        assert reg != self.reg


class TestPolygonSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = PolygonSkyRegion(SkyCoord([3, 4, 3] * u.deg, [3, 4, 4] * u.deg),
                           meta=meta, visual=visual)

    expected_repr = ('<PolygonSkyRegion(vertices=<SkyCoord (ICRS): (ra, '
                     'dec) in deg\n    [(3., 3.), (4., 4.), (3., 4.)]>)>')
    expected_str = ('Region: PolygonSkyRegion\nvertices: <SkyCoord (ICRS):'
                    ' (ra, dec) in deg\n    [(3., 3.), (4., 4.), (3., 4.)]>')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.vertices.ra.deg, [3, 4, 3])
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_transformation(self, wcs):

        pixpoly = self.reg.to_pixel(wcs)

        assert_allclose(pixpoly.vertices.x, [11.187992, 10.976332, 11.024032],
                        atol=1e-5)
        assert_allclose(pixpoly.vertices.y, [1.999486, 2.039001, 2.077076],
                        atol=1e-5)

        poly = pixpoly.to_sky(wcs)

        # TODO: we should probably assert something about frame
        # attributes, or generally some better way to check if two
        # SkyCoord are the same? For now, we use the following line to
        # transform back to ICRS (`poly` is in Galactic, same as WCS)
        poly = PolygonSkyRegion(
            vertices=poly.vertices.transform_to(self.reg.vertices))
        assert_quantity_allclose(poly.vertices.data.lon,
                                 self.reg.vertices.data.lon, atol=1e-3 * u.deg)
        assert_quantity_allclose(poly.vertices.data.lat,
                                 self.reg.vertices.data.lat, atol=1e-3 * u.deg)

    def test_contains(self, wcs):
        position = SkyCoord([1, 3.25] * u.deg, [2, 3.75] * u.deg)
        # 1,2 is outside, 3.25,3.75 should be inside the triangle...
        assert all(self.reg.contains(position, wcs)
                   == np.array([False, True], dtype='bool'))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.vertices = SkyCoord([3, 4, 3], [3, 4, 6], unit='deg')
        assert reg != self.reg


class TestRegionPolygonPixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RegularPolygonPixelRegion(PixCoord(50, 50), 8, 20, angle=25 * u.deg,
                                    meta=meta, visual=visual)
    inside = [(40, 40)]
    outside = [(20, 20), (80, 90)]
    expected_area = 1131.37085
    expected_repr = ('<RegularPolygonPixelRegion(center=PixCoord(x=50, y=50), '
                     'nvertices=8, radius=20, angle=25.0 deg)>')

    expected_str = ('Region: RegularPolygonPixelRegion\n'
                    'center: PixCoord(x=50, y=50)\n'
                    'nvertices: 8\n'
                    'radius: 20\nangle: 25.0 deg')

    def test_copy(self):
        reg = self.reg.copy()
        x_expected = copy(self.reg.vertices.x)
        y_expected = copy(self.reg.vertices.y)
        assert_allclose(reg.vertices.x, x_expected)
        assert_allclose(reg.vertices.y, y_expected)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_bounding_box(self):
        bbox = self.reg.bounding_box
        assert bbox == RegionBoundingBox(ixmin=31, ixmax=70, iymin=31,
                                         iymax=70)

    def test_to_mask(self):
        # The true area of this polygon is 3
        # We only do very low-precision asserts below,
        # because results can be unstable with points
        # on the edge of the polygon.
        # Basically we check that mask.data is filled
        # with something useful at all.

        # Bounding box and output shape is independent of subpixels,
        # so we only assert on it once here, not in the other cases below
        mask = self.reg.to_mask(mode='center', subpixels=1)
        assert 1130 <= np.sum(mask.data) <= 1135
        assert mask.bbox == RegionBoundingBox(ixmin=31, ixmax=70, iymin=31,
                                              iymax=70)
        assert mask.data.shape == (39, 39)

        # Test more cases for to_mask
        # This example is with the default: subpixels=5
        mask = self.reg.to_mask(mode='subpixels')
        assert 1130 <= np.sum(mask.data) <= 1135

        mask = self.reg.to_mask(mode='subpixels', subpixels=8)
        assert 1130 <= np.sum(mask.data) <= 1135

        mask = self.reg.to_mask(mode='subpixels', subpixels=9)
        assert 1130 <= np.sum(mask.data) <= 1135

        mask = self.reg.to_mask(mode='subpixels', subpixels=10)
        assert 1130 <= np.sum(mask.data) <= 1135

        with pytest.raises(NotImplementedError):
            self.reg.to_mask(mode='exact')

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
    def test_as_artist(self):
        patch = self.reg.as_artist()
        expected = [[41.54763477, 68.12615574], [31.20614758, 56.84040287],
                    [31.87384426, 41.54763477], [43.15959713, 31.20614758],
                    [58.45236523, 31.87384426], [68.79385242, 43.15959713],
                    [68.12615574, 58.45236523], [56.84040287, 68.79385242],
                    [41.54763477, 68.12615574]]
        assert_allclose(patch.xy, expected)

    def test_rotate(self):
        reg = self.reg.rotate(self.reg.center, 20 * u.deg)
        assert reg.angle.value == 45.0
        x_vert = [35.85786438, 30., 35.85786438, 50., 64.14213562, 70.,
                  64.14213562, 50.]
        y_vert = [64.14213562, 50., 35.85786438, 30., 35.85786438, 50.,
                  64.14213562, 70.]
        assert_allclose(reg.vertices.x, x_vert)
        assert_allclose(reg.vertices.y, y_vert)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.radius = 25.
        assert reg != self.reg

    def test_to_polygon(self):
        polyreg = self.reg.to_polygon()
        assert polyreg.vertices == self.reg.vertices
        assert polyreg.meta == self.meta
        assert polyreg.visual == self.visual
        assert polyreg.origin == PixCoord(0, 0)
