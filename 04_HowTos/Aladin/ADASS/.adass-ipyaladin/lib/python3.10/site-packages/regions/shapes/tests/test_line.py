# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from regions._utils.optional_deps import HAS_MATPLOTLIB
from regions.core import PixCoord, RegionMeta, RegionVisual
from regions.shapes.line import LinePixelRegion, LineSkyRegion
from regions.shapes.tests.test_common import (BaseTestPixelRegion,
                                              BaseTestSkyRegion)
from regions.tests.helpers import make_simple_wcs


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestLinePixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = LinePixelRegion(PixCoord(3, 4), PixCoord(4, 4), meta=meta,
                          visual=visual)
    sample_box = [-2, 8, -1, 9]
    inside = []
    outside = [(3.1, 4.2), (5, 4)]
    expected_area = 0
    expected_repr = ('<LinePixelRegion(start=PixCoord(x=3, y=4), '
                     'end=PixCoord(x=4, y=4))>')
    expected_str = ('Region: LinePixelRegion\nstart: PixCoord(x=3, y=4)\n'
                    'end: PixCoord(x=4, y=4)')

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.start.xy == (3, 4)
        assert reg.end.xy == (4, 4)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.start.x, self.reg.start.x)
        assert_allclose(reg_new.start.y, self.reg.start.y)
        assert_allclose(reg_new.end.x, self.reg.end.x)
        assert_allclose(reg_new.end.y, self.reg.end.y)
        assert reg_new.meta == self.reg.meta
        assert reg_new.visual == self.reg.visual

        # test that converted region meta and visual are copies and not views
        reg_new.meta['text'] = 'new'
        reg_new.visual['color'] = 'green'
        assert reg_new.meta['text'] != self.reg.meta['text']
        assert reg_new.visual['color'] != self.reg.visual['color']

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
    def test_as_artist(self):
        patch = self.reg.as_artist()
        assert 'Arrow' in str(patch)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.start.xy, (1, 4))
        assert_allclose(reg.end.xy, (1, 5))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.start = PixCoord(1, 2)
        assert reg != self.reg


class TestLineSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    start = SkyCoord(3 * u.deg, 4 * u.deg, frame='galactic')
    end = SkyCoord(3 * u.deg, 5 * u.deg, frame='galactic')
    reg = LineSkyRegion(start, end, meta=meta, visual=visual)

    expected_repr = ('<LineSkyRegion(start=<SkyCoord (Galactic): (l, b) '
                     'in deg\n    (3., 4.)>, end=<SkyCoord (Galactic): '
                     '(l, b) in deg\n    (3., 5.)>)>')
    expected_str = ('Region: LineSkyRegion\nstart: <SkyCoord (Galactic): '
                    '(l, b) in deg\n    (3., 4.)>\nend: <SkyCoord '
                    '(Galactic): (l, b) in deg\n    (3., 5.)>')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.start.b.deg, 4)
        assert_allclose(reg.end.b.deg, 5)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_transformation(self, wcs):
        pixline = self.reg.to_pixel(wcs)

        assert_allclose(pixline.start.x, -50.5)
        assert_allclose(pixline.start.y, 299.5)
        assert_allclose(pixline.end.x, -50.5)
        assert_allclose(pixline.end.y, 349.5)

        skyline = pixline.to_sky(wcs)

        assert_quantity_allclose(skyline.start.data.lon,
                                 self.reg.start.data.lon)
        assert_quantity_allclose(skyline.start.data.lat,
                                 self.reg.start.data.lat)
        assert_quantity_allclose(skyline.end.data.lon, self.reg.end.data.lon)
        assert_quantity_allclose(skyline.end.data.lat, self.reg.end.data.lat)

    def test_contains(self, wcs):
        position = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        # lines do not contain things
        assert all(self.reg.contains(position, wcs)
                   == np.array([False, False], dtype='bool'))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.start = SkyCoord(1 * u.deg, 2 * u.deg, frame='galactic')
        assert reg != self.reg
        reg.start = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
        assert reg != self.reg
