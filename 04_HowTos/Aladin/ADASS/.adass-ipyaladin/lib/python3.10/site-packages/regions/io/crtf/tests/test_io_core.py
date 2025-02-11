# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest
from astropy.tests.helper import assert_quantity_allclose
from astropy.units.quantity import Quantity

from regions.io.crtf.io_core import _to_shape_list
from regions.io.crtf.read import _CRTFParser


def test_shape_crtf():
    reg_str = ('circle[[18h12m24s, -23d11m00s], 2.3arcsec], linewidth=2, '
               'coord=J2000, symsize=2')

    parser = _CRTFParser(reg_str)
    shape1 = parser.shapes[0]
    region = parser.shapes.to_regions()
    shape2 = _to_shape_list(region, 'fk5')[0]

    assert shape1.coordsys == shape2.coordsys
    assert shape1.region_type == shape2.region_type
    assert shape1.include == shape2.include
    assert_quantity_allclose(shape1.coord, shape2.coord)
    # assert dict(shape1.meta) == dict(shape2.meta)
    assert shape1.composite == shape2.composite

    # Pixel origin is assumed to be 0 origin. Needs to be checked though.
    reg_str = 'circle[[1pix, 2pix], 3pix], linewidth=2, coord=image, symsize=2'

    parser = _CRTFParser(reg_str)
    shape = parser.shapes[0]
    assert_quantity_allclose(shape.coord,
                             [Quantity('1'), Quantity('2'), Quantity('3')])


def test_valid_shape():
    reg_str = ('circle[[18h12m24s, -23d11m00s], 2.3arcsec], linewidth=2, '
               'coord=J2000, symsize=2')

    shape = _CRTFParser(reg_str).shapes[0]

    with pytest.raises(ValueError) as excinfo:
        shape.region_type = 'box'

    assert '"box" is not a valid region type' in str(excinfo.value)

    shape = _CRTFParser(reg_str).shapes[0]
    with pytest.raises(ValueError) as excinfo:
        shape.coordsys = 'hello'

    estr = '"hello" is not a valid coordinate reference frame'
    assert estr in str(excinfo.value)


def test_valid_ellipse_crtf():
    reg_str1 = 'ellipse[[1pix, 2pix], [3pix, 4pix], 5deg]'
    reg_str2 = 'ellipse[[1deg, 2deg], [3deg, 4deg], 5deg], coord=fk5'

    shape1 = _CRTFParser(reg_str1, 'warn').shapes[0]
    shape2 = _CRTFParser(reg_str2, 'warn').shapes[0]

    assert_quantity_allclose(shape1.coord[:2], [1, 2])
    assert_quantity_allclose(shape2.coord[:2],
                             [Quantity('1deg'), Quantity('2deg')])

    assert_quantity_allclose(shape1.coord[2:-1], [8, 6])
    assert_quantity_allclose(shape2.coord[2:-1],
                             [Quantity('8deg'), Quantity('6deg')])

    assert_quantity_allclose(shape1.coord[-1], Quantity('5deg'))
    assert_quantity_allclose(shape2.coord[-1], Quantity('5deg'))
