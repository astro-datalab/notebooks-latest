# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The tests in this file simply check what functionality is currently
implemented and doesn't check anything about correctness.
"""
import itertools
from collections import OrderedDict

import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from regions.core.core import PixelRegion, Region, SkyRegion
from regions.core.mask import RegionMask
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord
from regions.shapes.annulus import (CircleAnnulusPixelRegion,
                                    CircleAnnulusSkyRegion,
                                    EllipseAnnulusPixelRegion,
                                    EllipseAnnulusSkyRegion,
                                    RectangleAnnulusPixelRegion,
                                    RectangleAnnulusSkyRegion)
from regions.shapes.circle import CirclePixelRegion, CircleSkyRegion
from regions.shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from regions.shapes.point import PointPixelRegion, PointSkyRegion
from regions.shapes.polygon import PolygonPixelRegion, PolygonSkyRegion
from regions.shapes.rectangle import RectanglePixelRegion, RectangleSkyRegion

PIXEL_REGIONS = [
    CirclePixelRegion(PixCoord(3, 4), radius=5),
    CircleAnnulusPixelRegion(PixCoord(3, 4), 5, 7),
    EllipsePixelRegion(PixCoord(3, 4), width=7, height=5, angle=3 * u.deg),
    EllipseAnnulusPixelRegion(PixCoord(3, 4), 7, 9, 8, 10, angle=3 * u.deg),
    PolygonPixelRegion(PixCoord([1, 4, 3], [2, 4, 4])),
    RectanglePixelRegion(PixCoord(6, 5), width=3, height=5),
    RectangleAnnulusPixelRegion(PixCoord(6, 5), 3, 6, 5, 7),
    PointPixelRegion(PixCoord(1, 2))]

SKY_REGIONS = [
    CircleSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), radius=5 * u.deg),
    CircleAnnulusSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 5 * u.deg,
                           7 * u.deg),
    EllipseSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), width=7 * u.deg,
                     height=5 * u.deg, angle=3 * u.deg),
    EllipseAnnulusSkyRegion(SkyCoord(3 * u.deg, 4 * u.deg), 7 * u.deg,
                            9 * u.deg, 5 * u.deg, 9 * u.deg, angle=3 * u.deg),
    PolygonSkyRegion(SkyCoord([1, 4, 3] * u.deg, [2, 4, 4] * u.deg)),
    RectangleSkyRegion(SkyCoord(6 * u.deg, 5 * u.deg), width=3 * u.deg,
                       height=5 * u.deg),
    RectangleAnnulusSkyRegion(SkyCoord(6 * u.deg, 5 * u.deg), 3 * u.deg,
                              5 * u.deg, 5 * u.deg, 7 * u.deg),
    PointSkyRegion(SkyCoord(6 * u.deg, 5 * u.deg))]

MASK_MODES = ['center', 'exact', 'subpixels']
COMMON_WCS = WCS(naxis=2)
COMMON_WCS.wcs.ctype = 'RA---TAN', 'DEC--TAN'


def ids_func(arg):
    if isinstance(arg, Region):
        return arg.__class__.__name__
    else:
        return str(arg)


@pytest.mark.parametrize('region', [PIXEL_REGIONS[i] for i in (0, 2, 5)],
                         ids=ids_func)
def test_pix_in(region):
    assert PixCoord(5, 4) in region


@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_pix_area(region):
    # TODO: remove the pytest.skip once polygon area is implemented
    if isinstance(region, PolygonPixelRegion):
        pytest.skip()

    area = region.area
    assert not isinstance(area, u.Quantity)


@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_pix_to_sky(region):
    try:
        sky_region = region.to_sky(COMMON_WCS)
        assert isinstance(sky_region, SkyRegion)
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.parametrize(('region', 'mode'),
                         itertools.product(PIXEL_REGIONS, MASK_MODES),
                         ids=ids_func)
def test_pix_to_mask(region, mode):
    try:
        mask = region.to_mask(mode=mode)
        assert isinstance(mask, RegionMask)
    except NotImplementedError:
        pytest.xfail()


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_in(region):
    region.contains(SkyCoord(1 * u.deg, 1 * u.deg, frame='icrs'), COMMON_WCS)


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_in_array(region):
    region.contains(SkyCoord([1, 2, 3] * u.deg, [3, 2, 1] * u.deg,
                             frame='icrs'), COMMON_WCS)


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_sky_to_pix(region):
    pix_region = region.to_pixel(wcs=COMMON_WCS)
    assert isinstance(pix_region, PixelRegion)


@pytest.mark.parametrize('region', PIXEL_REGIONS, ids=ids_func)
def test_attribute_validation_pixel_regions(region):
    invalid_values = dict(center=[PixCoord([1, 2], [2, 3]), 1,
                                  SkyCoord(1 * u.deg, 2 * u.deg),
                                  (10, 10), (10 * u.deg, 10 * u.deg)],
                          radius=[u.Quantity('1deg'), [1], PixCoord(1, 2),
                                  3 * u.km, 0.0, -10.],
                          angle=[u.Quantity([1 * u.deg, 2 * u.deg]), 2,
                                 PixCoord(1, 2), 3 * u.km],
                          vertices=[u.Quantity('1'), 2, PixCoord(1, 2),
                                    PixCoord([[1, 2]], [[2, 3]]), 3 * u.km,
                                    (10, 10), (10 * u.deg, 10 * u.deg)])
    invalid_values['width'] = invalid_values['radius']
    invalid_values['height'] = invalid_values['radius']
    invalid_values['inner_height'] = invalid_values['radius']
    invalid_values['inner_width'] = invalid_values['radius']
    invalid_values['outer_height'] = invalid_values['radius']
    invalid_values['outer_width'] = invalid_values['radius']
    invalid_values['inner_radius'] = invalid_values['radius']
    invalid_values['outer_radius'] = invalid_values['radius']
    invalid_values['start'] = invalid_values['center']
    invalid_values['end'] = invalid_values['radius']

    for attr in invalid_values:
        if hasattr(region, attr):
            for val in invalid_values.get(attr, None):
                with pytest.raises(ValueError) as excinfo:
                    setattr(region, attr, val)
                assert f'{attr!r} must' in str(excinfo.value)


@pytest.mark.parametrize('region', SKY_REGIONS, ids=ids_func)
def test_attribute_validation_sky_regions(region):
    invalid_values = dict(center=[PixCoord([1, 2], [2, 3]), 1,
                                  SkyCoord([1 * u.deg], [2 * u.deg]),
                                  (10, 10), (10 * u.deg, 10 * u.deg)],
                          radius=[u.Quantity([1 * u.deg, 5 * u.deg]),
                                  [1], SkyCoord(1 * u.deg, 2 * u.deg),
                                  1, 3 * u.km, 0.0 * u.deg, -10. * u.deg],
                          angle=[u.Quantity([1 * u.deg, 2 * u.deg]), 2,
                                 SkyCoord(1 * u.deg, 2 * u.deg), 3. * u.km],
                          vertices=[u.Quantity('1deg'), 2,
                                    SkyCoord(1 * u.deg, 2 * u.deg),
                                    SkyCoord([[1 * u.deg, 2 * u.deg]],
                                             [[2 * u.deg, 3 * u.deg]]),
                                    3 * u.km, (10, 10),
                                    (10 * u.deg, 10 * u.deg)])

    invalid_values['width'] = invalid_values['radius']
    invalid_values['height'] = invalid_values['radius']
    invalid_values['inner_height'] = invalid_values['radius']
    invalid_values['inner_width'] = invalid_values['radius']
    invalid_values['outer_height'] = invalid_values['radius']
    invalid_values['outer_width'] = invalid_values['radius']
    invalid_values['inner_radius'] = invalid_values['radius']
    invalid_values['outer_radius'] = invalid_values['radius']
    invalid_values['start'] = invalid_values['center']
    invalid_values['end'] = invalid_values['radius']

    for attr in invalid_values:
        if hasattr(region, attr):
            for val in invalid_values.get(attr, None):
                with pytest.raises(ValueError) as excinfo:
                    setattr(region, attr, val)
                assert f'{attr!r} must' in str(excinfo.value)


@pytest.mark.parametrize('region', PIXEL_REGIONS + SKY_REGIONS, ids=ids_func)
def test_metadata(region):
    region.meta = {'text': 'hello'}
    region.visual = {'color': 'blue'}
    assert isinstance(region.meta, RegionMeta)
    assert isinstance(region.visual, RegionVisual)

    # dict subclasses are allowed
    region.meta = OrderedDict({'text': 'hello'})
    region.visual = OrderedDict({'color': 'blue'})
    assert isinstance(region.meta, RegionMeta)
    assert isinstance(region.visual, RegionVisual)

    with pytest.raises(ValueError):
        region.meta = 1
    with pytest.raises(ValueError):
        region.visual = 'blue'
