# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file sets up detailed tests for computing masks with reference
images.
"""
import itertools

import astropy.units as u
import pytest

from regions.core import PixCoord, RegionMask
from regions.shapes.circle import CirclePixelRegion
from regions.shapes.ellipse import EllipsePixelRegion
from regions.shapes.polygon import PolygonPixelRegion
from regions.shapes.rectangle import RectanglePixelRegion

REGIONS = [CirclePixelRegion(PixCoord(3.981987, 4.131378), radius=3.3411),
           EllipsePixelRegion(PixCoord(5.981987, 4.131378), width=10.4466,
                              height=6.6822, angle=12 * u.deg),
           RectanglePixelRegion(PixCoord(5.981987, 4.131378), width=7.2233,
                                height=4.3411, angle=12 * u.deg),
           PolygonPixelRegion(PixCoord([-2.334, 3.631, 1.122, -1.341],
                                       [-3.121, -2.118, 2.987, 1.759]))]

MODES = [{'mode': 'center'},
         {'mode': 'exact'},
         {'mode': 'subpixels', 'subpixels': 1},
         {'mode': 'subpixels', 'subpixels': 5},
         {'mode': 'subpixels', 'subpixels': 18}]


def label(value):
    if isinstance(value, CirclePixelRegion):
        return 'circ'
    elif isinstance(value, EllipsePixelRegion):
        return 'elli'
    elif isinstance(value, RectanglePixelRegion):
        return 'rect'
    elif isinstance(value, PolygonPixelRegion):
        return 'poly'
    else:
        return '-'.join(f'{key}_{value}'
                        for key, value in sorted(value.items()))


@pytest.mark.array_compare(file_format='text', write_kwargs={'fmt': '%12.8e'})
@pytest.mark.parametrize(('region', 'mode'),
                         itertools.product(REGIONS, MODES), ids=label)
def test_to_mask(region, mode):
    try:
        mask = region.to_mask(**mode)
    except NotImplementedError:
        pytest.xfail()
    assert isinstance(mask, RegionMask)
    return mask.data.astype(float)
