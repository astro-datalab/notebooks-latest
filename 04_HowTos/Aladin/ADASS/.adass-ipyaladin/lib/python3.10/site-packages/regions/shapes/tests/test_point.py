# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from regions._utils.optional_deps import HAS_MATPLOTLIB
from regions.core import PixCoord, RegionMeta, RegionVisual
from regions.shapes.point import PointPixelRegion, PointSkyRegion
from regions.shapes.tests.test_common import (BaseTestPixelRegion,
                                              BaseTestSkyRegion)
from regions.tests.helpers import make_simple_wcs


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


class TestPointPixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = PointPixelRegion(PixCoord(3, 4), meta=meta, visual=visual)
    sample_box = [-2, 8, -1, 9]
    inside = []
    outside = [(3.1, 4.2), (5, 4)]
    expected_area = 0
    expected_repr = '<PointPixelRegion(center=PixCoord(x=3, y=4))>'
    expected_str = 'Region: PointPixelRegion\ncenter: PixCoord(x=3, y=4)'

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert reg_new.meta == self.reg.meta
        assert reg_new.visual == self.reg.visual

        # test that converted region meta and visual are copies and not views
        reg_new.meta['text'] = 'new'
        reg_new.visual['color'] = 'green'
        assert reg_new.meta['text'] != self.reg.meta['text']
        assert reg_new.visual['color'] != self.reg.visual['color']

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
    def test_as_artist(self):
        artist = self.reg.as_artist()

        assert artist.get_data() == ([3], [4])

        artist = self.reg.as_artist(origin=(1, 1))

        assert artist.get_data() == ([2], [3])

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.center = PixCoord(1, 2)
        assert reg != self.reg


class TestPointSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = PointSkyRegion(SkyCoord(3, 4, unit='deg'), meta=meta, visual=visual)

    expected_repr = ('<PointSkyRegion(center=<SkyCoord (ICRS): (ra, dec) '
                     'in deg\n    (3., 4.)>)>')
    expected_str = ('Region: PointSkyRegion\ncenter: <SkyCoord (ICRS): '
                    '(ra, dec) in deg\n    (3., 4.)>')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_contains(self, wcs):
        position = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)
        # points do not contain things
        assert all(self.reg.contains(position, wcs)
                   == np.array([False, False], dtype='bool'))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.center = SkyCoord(1, 2, unit='deg')
        assert reg != self.reg
