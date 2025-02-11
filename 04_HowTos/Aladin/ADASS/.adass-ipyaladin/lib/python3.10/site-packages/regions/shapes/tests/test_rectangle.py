# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from numpy.testing import assert_allclose, assert_equal

from regions._utils.optional_deps import HAS_MATPLOTLIB
from regions.core import PixCoord, RegionMeta, RegionVisual
from regions.shapes.rectangle import RectanglePixelRegion, RectangleSkyRegion
from regions.shapes.tests.test_common import (BaseTestPixelRegion,
                                              BaseTestSkyRegion)
from regions.tests.helpers import make_simple_wcs


@pytest.fixture(scope='session', name='wcs')
def wcs_fixture():
    filename = get_pkg_data_filename('data/example_header.fits')
    header = fits.getheader(filename)
    return WCS(header)


def test_corners():
    xc, yc = 2, 2
    angle = 30 * u.deg
    width = 2
    height = 1
    reg = RectanglePixelRegion(PixCoord(xc, yc), width=width, height=height,
                               angle=angle)

    y1 = yc + np.cos(angle) * (height / 2) + np.sin(angle) * (width / 2)
    x1 = xc + np.cos(angle) * (width / 2) - np.sin(angle) * (height / 2)

    assert (x1, y1) in reg.corners

    reg = RectanglePixelRegion(PixCoord(xc, yc), width=width, height=height,
                               angle=90 * u.deg)
    # simple case: rotate by 90
    np.testing.assert_allclose([(2.5, 1.), (2.5, 3.), (1.5, 3.), (1.5, 1.)],
                               reg.corners)

    reg = RectanglePixelRegion(center=PixCoord(xc, yc), width=width,
                               height=height, angle=0 * u.deg)
    # simpler case: rotate by 0
    np.testing.assert_array_equal([(1, 1.5), (3, 1.5), (3, 2.5), (1, 2.5)],
                                  reg.corners)

    poly = reg.to_polygon()
    assert len(poly.vertices) == 4


class TestRectanglePixelRegion(BaseTestPixelRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RectanglePixelRegion(center=PixCoord(3, 4), width=4, height=3,
                               angle=5 * u.deg, meta=meta, visual=visual)
    sample_box = [-2, 8, -1, 9]
    inside = [(4.5, 4)]
    outside = [(5, 2.5)]
    expected_area = 12
    expected_repr = ('<RectanglePixelRegion(center=PixCoord(x=3, y=4), '
                     'width=4, height=3, angle=5.0 deg)>')
    expected_str = ('Region: RectanglePixelRegion\ncenter: PixCoord(x=3, '
                    'y=4)\nwidth: 4\nheight: 3\nangle: 5.0 deg')

    def test_copy(self):
        reg = self.reg.copy()
        assert reg.center.xy == (3, 4)
        assert reg.width == 4
        assert reg.height == 3
        assert_allclose(reg.angle.to_value('deg'), 5)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_pix_sky_roundtrip(self):
        wcs = make_simple_wcs(SkyCoord(2 * u.deg, 3 * u.deg), 0.1 * u.deg, 20)
        reg_new = self.reg.to_sky(wcs).to_pixel(wcs)
        assert_allclose(reg_new.center.x, self.reg.center.x)
        assert_allclose(reg_new.center.y, self.reg.center.y)
        assert_allclose(reg_new.width, self.reg.width)
        assert_allclose(reg_new.height, self.reg.height)
        assert_quantity_allclose(reg_new.angle, self.reg.angle)
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
        # Note: `reg.center` is the center, `patch.xy` is the lower-left
        # corner
        assert_allclose(patch.xy, (1.138344, 2.331396), atol=1e-3)
        assert_allclose(patch.get_width(), 4)
        assert_allclose(patch.get_height(), 3)
        assert_allclose(patch.angle, 5)

    def test_rotate(self):
        reg = self.reg.rotate(PixCoord(2, 3), 90 * u.deg)
        assert_allclose(reg.center.xy, (1, 4))
        assert_allclose(reg.angle.to_value('deg'), 95)

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.angle = 35 * u.deg
        assert reg != self.reg

    # temporarily disable sync=True test due to random failures
    # @pytest.mark.parametrize('sync', (False, True))
    @pytest.mark.parametrize('sync', (False,))
    def test_as_mpl_selector(self, sync):
        plt = pytest.importorskip('matplotlib.pyplot')
        from matplotlib.testing.widgets import do_event

        rng = np.random.default_rng(0)
        data = rng.random((16, 16))
        mask = np.zeros_like(data)

        ax = plt.subplot(1, 1, 1)
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(mode='subpixels', subpixels=10).to_image(data.shape)

        # For now this will only work with unrotated rectangles. Once
        # this works with rotated rectangles, the following exception
        # check can be removed as well as the ``angle=0 * u.deg`` in the
        # call to copy() below.
        with pytest.raises(NotImplementedError,
                           match=('Cannot create matplotlib selector for rotated rectangle.')):
            self.reg.as_mpl_selector(ax)

        region = self.reg.copy(angle=0 * u.deg)

        selector = region.as_mpl_selector(ax, callback=update_mask, sync=sync)

        do_event(selector, 'press', xdata=7.3, ydata=4.4, button=1)
        do_event(selector, 'onmove', xdata=9.3, ydata=5.4, button=1)
        do_event(selector, 'release', xdata=9.3, ydata=5.4, button=1)

        ax.figure.canvas.draw()

        if sync:
            assert_allclose(region.center.x, 8.3)
            assert_allclose(region.center.y, 4.9)
            assert_allclose(region.width, 2)
            assert_allclose(region.height, 1)
            assert_quantity_allclose(region.angle, 0 * u.deg)

            assert_equal(mask, region.to_mask(mode='subpixels', subpixels=10).to_image(data.shape))

        else:
            assert_allclose(region.center.x, 3)
            assert_allclose(region.center.y, 4)
            assert_allclose(region.width, 4)
            assert_allclose(region.height, 3)
            assert_quantity_allclose(region.angle, 0 * u.deg)

            assert_equal(mask, 0)

        with pytest.raises(AttributeError, match=('Cannot attach more than one selector to a reg')):
            region.as_mpl_selector(ax)

    @pytest.mark.parametrize('anywhere', (False, True))
    def test_mpl_selector_drag(self, anywhere):
        """
        Test dragging of entire region from central handle and anywhere.
        """
        plt = pytest.importorskip('matplotlib.pyplot')
        from matplotlib.testing.widgets import (
            do_event)  # click_and_drag  # MPL_VERSION >= 36

        rng = np.random.default_rng(0)
        data = rng.random((16, 16))
        mask = np.zeros_like(data)

        ax = plt.subplot(1, 1, 1)
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(mode='subpixels', subpixels=10).to_image(data.shape)

        region = self.reg.copy(angle=0 * u.deg)

        selector = region.as_mpl_selector(ax, callback=update_mask,
                                          drag_from_anywhere=anywhere)
        assert selector.drag_from_anywhere is anywhere
        assert region._mpl_selector.drag_from_anywhere is anywhere

        # click_and_drag(selector, start=(3, 4), end=(3.5, 4.5))
        do_event(selector, 'press', xdata=3, ydata=4, button=1)
        do_event(selector, 'onmove', xdata=3.5, ydata=4.5, button=1)
        do_event(selector, 'release', xdata=3.5, ydata=4.5, button=1)

        ax.figure.canvas.draw()

        assert_allclose(region.center.x, 3.5)
        assert_allclose(region.center.y, 4.5)
        assert_allclose(region.width, 4)
        assert_allclose(region.height, 3)

        do_event(selector, 'press', xdata=3.25, ydata=4.25, button=1)
        do_event(selector, 'onmove', xdata=4.25, ydata=5.25, button=1)
        do_event(selector, 'release', xdata=4.25, ydata=5.25, button=1)

        ax.figure.canvas.draw()

        # For drag_from_anywhere=False this will have created a new 1x1 rectangle.
        if anywhere:
            assert_allclose(region.center.x, 4.5)
            assert_allclose(region.center.y, 5.5)
            assert_allclose(region.width, 4)
            assert_allclose(region.height, 3)
        else:
            assert_allclose(region.center.x, 4.5)
            assert_allclose(region.center.y, 5.5)

        assert_equal(mask, region.to_mask(mode='subpixels', subpixels=10).to_image(data.shape))

    @pytest.mark.parametrize('userargs',
                             ({'useblit': True},
                              {'grab_range': 20, 'minspanx': 5, 'minspany': 4},
                              {'props': {'facecolor': 'blue', 'linewidth': 2}},
                              {'twit': 'gumby'}))
    def test_mpl_selector_kwargs(self, userargs):
        """
        Test that additional kwargs are passed to selector.
        """
        plt = pytest.importorskip('matplotlib.pyplot')

        rng = np.random.default_rng(0)
        data = rng.random((16, 16))
        mask = np.zeros_like(data)

        ax = plt.subplot(1, 1, 1)
        ax.imshow(data)

        def update_mask(reg):
            mask[:] = reg.to_mask(mode='subpixels', subpixels=10).to_image(data.shape)

        region = self.reg.copy(angle=0 * u.deg)

        if 'twit' in userargs:
            with pytest.raises(TypeError, match=(r'__init__.. got an unexpected keyword argument')):
                selector = region.as_mpl_selector(ax, callback=update_mask, **userargs)
        else:
            selector = region.as_mpl_selector(ax, callback=update_mask, **userargs)
            assert region._mpl_selector.artists[0].get_edgecolor() == (0, 0, 1, 1)

            if 'props' in userargs:
                assert region._mpl_selector.artists[0].get_facecolor() == (0, 0, 1, 1)
                assert region._mpl_selector.artists[0].get_linewidth() == 2
            else:
                assert region._mpl_selector.artists[0].get_facecolor() == (0, 0, 0, 0)
                assert region._mpl_selector.artists[0].get_linewidth() == 1

                for key, val in userargs.items():
                    assert getattr(region._mpl_selector, key) == val
                    assert getattr(selector, key) == val


def test_rectangular_pixel_region_bbox():
    # odd sizes
    width = 7
    height = 3
    a = RectanglePixelRegion(PixCoord(50, 50), width=width, height=height,
                             angle=0. * u.deg)
    assert a.bounding_box.shape == (height, width)

    a = RectanglePixelRegion(PixCoord(50.5, 50.5), width=width, height=height,
                             angle=0. * u.deg)
    assert a.bounding_box.shape == (height + 1, width + 1)

    a = RectanglePixelRegion(PixCoord(50, 50), width=width, height=height,
                             angle=90. * u.deg)
    assert a.bounding_box.shape == (width, height)

    # even sizes
    width = 8
    height = 4
    a = RectanglePixelRegion(PixCoord(50, 50), width=width, height=height,
                             angle=0. * u.deg)
    assert a.bounding_box.shape == (height + 1, width + 1)

    a = RectanglePixelRegion(PixCoord(50.5, 50.5), width=width, height=height,
                             angle=0. * u.deg)
    assert a.bounding_box.shape == (height, width)

    a = RectanglePixelRegion(PixCoord(50.5, 50.5), width=width, height=height,
                             angle=90. * u.deg)
    assert a.bounding_box.shape == (width, height)


class TestRectangleSkyRegion(BaseTestSkyRegion):
    meta = RegionMeta({'text': 'test'})
    visual = RegionVisual({'color': 'blue'})
    reg = RectangleSkyRegion(center=SkyCoord(3, 4, unit='deg'),
                             width=4 * u.deg, height=3 * u.deg,
                             angle=5 * u.deg, meta=meta, visual=visual)

    expected_repr = ('<RectangleSkyRegion(center=<SkyCoord (ICRS): (ra, dec) '
                     'in deg\n    (3., 4.)>, width=4.0 deg, height=3.0 deg, '
                     'angle=5.0 deg)>')
    expected_str = ('Region: RectangleSkyRegion\ncenter: <SkyCoord '
                    '(ICRS): (ra, dec) in deg\n    (3., 4.)>\nwidth: '
                    '4.0 deg\nheight: 3.0 deg\nangle: 5.0 deg')

    def test_copy(self):
        reg = self.reg.copy()
        assert_allclose(reg.center.ra.deg, 3)
        assert_allclose(reg.width.to_value('deg'), 4)
        assert_allclose(reg.height.to_value('deg'), 3)
        assert_allclose(reg.angle.to_value('deg'), 5)
        assert reg.meta == self.meta
        assert reg.visual == self.visual

    def test_contains(self, wcs):
        position = SkyCoord([1, 3] * u.deg, [2, 4] * u.deg)
        # 1,2 is outside, 3,4 is the center and is inside
        assert all(self.reg.contains(position, wcs)
                   == np.array([False, True], dtype='bool'))

    def test_eq(self):
        reg = self.reg.copy()
        assert reg == self.reg
        reg.angle = 10 * u.deg
        assert reg != self.reg
