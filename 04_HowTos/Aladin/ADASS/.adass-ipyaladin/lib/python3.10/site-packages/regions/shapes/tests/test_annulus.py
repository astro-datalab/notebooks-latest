# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from numpy.testing import assert_allclose

from regions.core import PixCoord, RegionMeta, RegionVisual
from regions.shapes.annulus import (CircleAnnulusPixelRegion,
                                    CircleAnnulusSkyRegion,
                                    EllipseAnnulusPixelRegion,
                                    EllipseAnnulusSkyRegion,
                                    RectangleAnnulusPixelRegion,
                                    RectangleAnnulusSkyRegion)
from regions.shapes.tests.test_common import (BaseTestPixelRegion,
                                              BaseTestSkyRegion)
from regions.tests.helpers import make_simple_wcs


class TestCircleAnnulusPixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = CircleAnnulusPixelRegion(PixCoord(3, 4), inner_radius=2,
                                   outer_radius=3, meta=meta, visual=visual)
    sample_box = [0, 6, 1, 7]
    inside = [(3, 2)]
    outside = [(3, 0)]
    expected_area = 5 * np.pi
    expected_repr = ('<CircleAnnulusPixelRegion(center=PixCoord(x=3, y=4), '
                     'inner_radius=2, outer_radius=3)>')
    expected_str = ('Region: CircleAnnulusPixelRegion\ncenter: '
                    'PixCoord(x=3, y=4)\ninner_radius: 2\nouter_radius: 3')

    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_init(self):
        assert_quantity_allclose(self.reg.center.x, 3)
        assert_quantity_allclose(self.reg.center.y, 4)
        assert_quantity_allclose(self.reg.inner_radius, 2)
        assert_quantity_allclose(self.reg.outer_radius, 3)

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.inner_radius == 2
        assert reg.outer_radius == 3
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.inner_radius, self.reg.inner_radius)
        assert_allclose(reg_new.outer_radius, self.reg.outer_radius)
        assert reg_new.meta == self.reg.meta
        assert reg_new.visual == self.reg.visual

        # test that converted region meta and visual are copies and not views
        reg_new.meta['text'] = 'new'
        reg_new.visual['color'] = 'green'
        assert reg_new.meta['text'] != self.reg.meta['text']
        assert reg_new.visual['color'] != self.reg.visual['color']

    def test_transformation(self):
        skyannulus = self.reg.to_sky(wcs=self.wcs)
        assert isinstance(skyannulus, CircleAnnulusSkyRegion)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.inner_radius = 3
        assert reg != self.reg


class TestCircleAnnulusSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = CircleAnnulusSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg),
                                 20 * u.arcsec, 30 * u.arcsec, meta=meta,
                                 visual=visual)
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    expected_repr = ('<CircleAnnulusSkyRegion(center=<SkyCoord (ICRS): '
                     '(ra, dec) in deg\n    (3., 4.)>, inner_radius=20.0 '
                     'arcsec, outer_radius=30.0 arcsec)>')
    expected_str = ('Region: CircleAnnulusSkyRegion\ncenter: <SkyCoord '
                    '(ICRS): (ra, dec) in deg\n    (3., 4.)>\ninner_radius: '
                    '20.0 arcsec\nouter_radius: 30.0 arcsec')

    def test_init(self):
        assert_quantity_allclose(self.reg.center.ra, self.skycoord.ra)
        assert_quantity_allclose(self.reg.inner_radius, 20 * u.arcsec)
        assert_quantity_allclose(self.reg.outer_radius, 30 * u.arcsec)

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.inner_radius.to_value('arcsec'), 20)
        assert_allclose(reg.outer_radius.to_value('arcsec'), 30)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_contains(self):
        assert not self.reg.contains(self.skycoord, self.wcs)
        test_coord = SkyCoord(3 * u.deg, 10 * u.deg, frame='icrs')
        assert not self.reg.contains(test_coord, self.wcs)

    def test_transformation(self):
        pixannulus = self.reg.to_pixel(wcs=self.wcs)
        assert isinstance(pixannulus, CircleAnnulusPixelRegion)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.inner_radius = 10 * u.arcsec
        assert reg != self.reg


class TestEllipseAnnulusPixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = EllipseAnnulusPixelRegion(PixCoord(3, 4), 2, 5, 5, 8, meta=meta,
                                    visual=visual)
    sample_box = [1, 6, 0, 9]
    inside = [(3, 7)]
    outside = [(3, 4)]
    expected_area = 7.5 * np.pi
    expected_repr = ('<EllipseAnnulusPixelRegion(center=PixCoord(x=3, y=4), '
                     'inner_width=2, outer_width=5, inner_height=5, '
                     'outer_height=8, angle=0.0 deg)>')
    expected_str = ('Region: EllipseAnnulusPixelRegion\ncenter: '
                    'PixCoord(x=3, y=4)\ninner_width: 2\nouter_width: 5\n'
                    'inner_height: 5\nouter_height: 8\nangle: 0.0 deg')

    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_init(self):
        assert_quantity_allclose(self.reg.center.x, 3)
        assert_quantity_allclose(self.reg.center.y, 4)
        assert_quantity_allclose(self.reg.inner_width, 2)
        assert_quantity_allclose(self.reg.inner_height, 5)
        assert_quantity_allclose(self.reg.outer_width, 5)
        assert_quantity_allclose(self.reg.outer_height, 8)

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.inner_width == 2
        assert reg.inner_height == 5
        assert reg.outer_width == 5
        assert reg.outer_height == 8
        assert_allclose(reg.angle.to_value('deg'), 0)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.inner_width, self.reg.inner_width)
        assert_allclose(reg_new.outer_width, self.reg.outer_width)
        assert_allclose(reg_new.inner_height, self.reg.inner_height)
        assert_allclose(reg_new.outer_height, self.reg.outer_height)
        assert_quantity_allclose(reg_new.angle, self.reg.angle)
        assert reg_new.meta == self.reg.meta
        assert reg_new.visual == self.reg.visual

        # test that converted region meta and visual are copies and not views
        reg_new.meta['text'] = 'new'
        reg_new.visual['color'] = 'green'
        assert reg_new.meta['text'] != self.reg.meta['text']
        assert reg_new.visual['color'] != self.reg.visual['color']

    def test_transformation(self):
        skyannulus = self.reg.to_sky(wcs=self.wcs)
        assert isinstance(skyannulus, EllipseAnnulusSkyRegion)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value('deg'), 90)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.outer_height = 10
        assert reg != self.reg


class TestEllipseAnnulusSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    reg = EllipseAnnulusSkyRegion(skycoord, 20 * u.arcsec, 50 * u.arcsec,
                                  50 * u.arcsec, 80 * u.arcsec, meta=meta,
                                  visual=visual)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    expected_repr = ('<EllipseAnnulusSkyRegion(center=<SkyCoord (ICRS): '
                     '(ra, dec) in deg\n    (3., 4.)>, inner_width=20.0 '
                     'arcsec, outer_width=50.0 arcsec, inner_height=50.0 '
                     'arcsec, outer_height=80.0 arcsec, angle=0.0 deg)>')
    expected_str = ('Region: EllipseAnnulusSkyRegion\ncenter: <SkyCoord '
                    '(ICRS): (ra, dec) in deg\n    (3., 4.)>\ninner_width: '
                    '20.0 arcsec\nouter_width: 50.0 arcsec\n'
                    'inner_height: 50.0 arcsec\nouter_height: '
                    '80.0 arcsec\nangle: 0.0 deg')

    def test_init(self):
        assert_quantity_allclose(self.reg.center.ra, self.skycoord.ra)
        assert_quantity_allclose(self.reg.inner_width, 20 * u.arcsec)
        assert_quantity_allclose(self.reg.inner_height, 50 * u.arcsec)

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.inner_width.to_value('arcsec'), 20)
        assert_allclose(reg.inner_height.to_value('arcsec'), 50)
        assert_allclose(reg.outer_width.to_value('arcsec'), 50)
        assert_allclose(reg.outer_height.to_value('arcsec'), 80)
        assert_allclose(reg.angle.to_value('deg'), 0)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_contains(self):
        assert not self.reg.contains(self.skycoord, self.wcs)
        test_coord = SkyCoord(3 * u.deg, 10 * u.deg, frame='icrs')
        assert not self.reg.contains(test_coord, self.wcs)

    def test_transformation(self):
        pixannulus = self.reg.to_pixel(wcs=self.wcs)
        assert isinstance(pixannulus, EllipseAnnulusPixelRegion)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.outer_height = 85 * u.arcsec
        assert reg != self.reg


class TestRectangleAnnulusPixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RectangleAnnulusPixelRegion(PixCoord(3, 4), 2, 5, 5, 8, meta=meta,
                                      visual=visual)
    sample_box = [1, 6, 0, 9]
    inside = [(3, 7)]
    outside = [(3, 4)]
    expected_area = 30
    expected_repr = ('<RectangleAnnulusPixelRegion(center=PixCoord(x=3, y=4), '
                     'inner_width=2, outer_width=5, inner_height=5, '
                     'outer_height=8, angle=0.0 deg)>')
    expected_str = ('Region: RectangleAnnulusPixelRegion\ncenter: '
                    'PixCoord(x=3, y=4)\ninner_width: 2\nouter_width: 5\n'
                    'inner_height: 5\nouter_height: 8\nangle: 0.0 deg')

    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    def test_init(self):
        assert_quantity_allclose(self.reg.center.x, 3)
        assert_quantity_allclose(self.reg.center.y, 4)
        assert_quantity_allclose(self.reg.inner_width, 2)
        assert_quantity_allclose(self.reg.inner_height, 5)
        assert_quantity_allclose(self.reg.outer_width, 5)
        assert_quantity_allclose(self.reg.outer_height, 8)

    def test_copy(self):
        reg = self.reg.copy()
        assert_quantity_allclose(reg.center.x, 3)
        assert_quantity_allclose(reg.center.y, 4)
        assert_quantity_allclose(reg.inner_width, 2)
        assert_quantity_allclose(reg.inner_height, 5)
        assert_quantity_allclose(reg.outer_width, 5)
        assert_quantity_allclose(reg.outer_height, 8)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.inner_width, self.reg.inner_width)
        assert_allclose(reg_new.outer_width, self.reg.outer_width)
        assert_allclose(reg_new.inner_height, self.reg.inner_height)
        assert_allclose(reg_new.outer_height, self.reg.outer_height)
        assert_quantity_allclose(reg_new.angle, self.reg.angle)
        assert reg_new.meta == self.reg.meta
        assert reg_new.visual == self.reg.visual

        # test that converted region meta and visual are copies and not views
        reg_new.meta['text'] = 'new'
        reg_new.visual['color'] = 'green'
        assert reg_new.meta['text'] != self.reg.meta['text']
        assert reg_new.visual['color'] != self.reg.visual['color']

    def test_transformation(self):
        skyannulus = self.reg.to_sky(wcs=self.wcs)
        assert isinstance(skyannulus, RectangleAnnulusSkyRegion)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value('deg'), 90)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.outer_height = 10
        assert reg != self.reg


class TestRectangleAnnulusSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    skycoord = SkyCoord(3 * u.deg, 4 * u.deg, frame='icrs')
    reg = RectangleAnnulusSkyRegion(skycoord, 20 * u.arcsec, 50 * u.arcsec,
                                    50 * u.arcsec, 80 * u.arcsec, meta=meta,
                                    visual=visual)
    wcs = make_simple_wcs(skycoord, 5 * u.arcsec, 20)

    expected_repr = ('<RectangleAnnulusSkyRegion(center=<SkyCoord (ICRS): '
                     '(ra, dec) in deg\n    (3., 4.)>, inner_width=20.0 '
                     'arcsec, outer_width=50.0 arcsec, inner_height=50.0 '
                     'arcsec, outer_height=80.0 arcsec, angle=0.0 deg)>')
    expected_str = ('Region: RectangleAnnulusSkyRegion\ncenter: <SkyCoord '
                    '(ICRS): (ra, dec) in deg\n    (3., 4.)>\ninner_width: '
                    '20.0 arcsec\nouter_width: 50.0 arcsec\n'
                    'inner_height: 50.0 arcsec\nouter_height: 80.0 arcsec\n'
                    'angle: 0.0 deg')

    def test_init(self):
        assert_quantity_allclose(self.reg.center.ra, self.skycoord.ra)
        assert_quantity_allclose(self.reg.inner_width, 20 * u.arcsec)
        assert_quantity_allclose(self.reg.inner_height, 50 * u.arcsec)

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.inner_width.to_value('arcsec'), 20)
        assert_allclose(reg.inner_height.to_value('arcsec'), 50)
        assert_allclose(reg.outer_width.to_value('arcsec'), 50)
        assert_allclose(reg.outer_height.to_value('arcsec'), 80)
        assert_allclose(reg.angle.to_value('deg'), 0)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_contains(self):
        assert not self.reg.contains(self.skycoord, self.wcs)
        test_coord = SkyCoord(3 * u.deg, 10 * u.deg, frame='icrs')
        assert not self.reg.contains(test_coord, self.wcs)

    def test_transformation(self):
        pixannulus = self.reg.to_pixel(wcs=self.wcs)
        assert isinstance(pixannulus, RectangleAnnulusPixelRegion)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.outer_height = 85 * u.arcsec
        assert reg != self.reg
