# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the pixcoord module.
"""

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from regions._utils.examples import make_example_dataset
from regions.core.pixcoord import PixCoord


@pytest.fixture(scope='session', name='wcs')
def fixture_wcs():
    config = dict(crpix=(18, 9), cdelt=(-10, 10), shape=(18, 36))
    dataset = make_example_dataset(config=config)
    return dataset.wcs


def test_pixcoord_copy_scalar():
    pc1 = PixCoord(x=1, y=2)
    pc2 = pc1.copy()
    pc1.x = 99
    pc1.y = 99
    assert pc2.xy == (1, 2)


def test_pixcoord_copy_array():
    pc1 = PixCoord(x=[1, 2], y=[3, 4])
    pc2 = pc1.copy()
    pc1.x[0] = 99
    pc1.y[0] = 99
    assert_equal(pc2.x, [1, 2])
    assert_equal(pc2.y, [3, 4])


def test_pixcoord_basic_dimension():
    with pytest.raises(ValueError):
        PixCoord(np.array([1, 2]), [3, 4, 5, 6])


def test_pixcoord_basics_scalar():
    pc0 = PixCoord(x=1, y=2)
    pc1 = PixCoord(x=np.array(1), y=2)
    pc2 = PixCoord(x=np.array(1), y=np.array(2))

    assert pc0 == pc1
    assert pc1 == pc2

    assert pc0.x == 1
    assert pc0.y == 2

    assert pc0.isscalar
    assert pc1.isscalar
    assert pc2.isscalar

    assert str(pc0) == 'PixCoord(x=1, y=2)'
    assert repr(pc0) == 'PixCoord(x=1, y=2)'

    with pytest.raises(TypeError):
        _ = len(pc0)

    with pytest.raises(IndexError):
        _ = pc0[0]


def test_pixcoord_basics_array_1d():
    pc0 = PixCoord(x=[1, 2, 3], y=[11, 22, 33])
    pc1 = PixCoord(x=1, y=[1, 2, 3])

    assert_equal(pc0.x, [1, 2, 3])
    assert_equal(pc0.y, [11, 22, 33])

    assert_equal(pc1.x, [1, 1, 1])
    assert_equal(pc1.y, [1, 2, 3])

    assert not pc0.isscalar
    assert not pc1.isscalar

    assert str(pc0) == 'PixCoord(x=[1 2 3], y=[11 22 33])'
    assert str(pc1) == 'PixCoord(x=[1 1 1], y=[1 2 3])'

    assert repr(pc0) == 'PixCoord(x=[1 2 3], y=[11 22 33])'
    assert repr(pc1) == 'PixCoord(x=[1 1 1], y=[1 2 3])'

    assert len(pc0) == 3
    assert len(pc1) == 3

    pc2 = list(iter(pc0))[-1]
    assert pc2.x == 3
    assert pc2.y == 33

    pc3 = pc0[-1]
    assert pc3.isscalar
    assert pc3.x == 3
    assert pc3.y == 33

    pc4 = pc0[1:]
    assert len(pc4) == 2
    assert_equal(pc4.x, [2, 3])
    assert_equal(pc4.y, [22, 33])


def test_pixcoord_basics_array_2d():
    pc0 = PixCoord([[1, 2, 3], [4, 5, 6]], [[11, 12, 13], [14, 15, 16]])
    assert_equal(pc0.x, [[1, 2, 3], [4, 5, 6]])
    assert_equal(pc0.y, [[11, 12, 13], [14, 15, 16]])
    assert pc0.isscalar is False

    expected = ('PixCoord(x=[[1 2 3]\n [4 5 6]], '
                'y=[[11 12 13]\n [14 15 16]])')
    assert str(pc0) == expected
    assert repr(pc0) == expected
    assert len(pc0) == 2

    assert_equal(pc0[0].x, pc0.x[0])
    assert_equal(pc0[0].y, pc0.y[0])

    pc1 = list(iter(pc0))
    assert_equal(pc1[0].x, pc0.x[0])
    assert_equal(pc1[0].y, pc0.y[0])


def test_pixcoord_to_sky_scalar(wcs):
    pc1 = PixCoord(x=18, y=9)

    sc = pc1.to_sky(wcs=wcs)
    assert sc.name == 'galactic'
    assert_allclose(sc.data.lon.deg, 349.88093995435634)
    assert_allclose(sc.data.lat.deg, 10.003028030623508)

    pc2 = PixCoord.from_sky(skycoord=sc, wcs=wcs)
    assert isinstance(pc2.x, float)
    assert pc2.isscalar

    assert pc1 == pc2


def test_pixcoord_to_sky_array_1d(wcs):
    pc1 = PixCoord(x=[17, 18], y=[8, 9])

    sc = pc1.to_sky(wcs=wcs)
    assert sc.name == 'galactic'
    assert_allclose(sc.data.lon.deg, [0, 349.88094])
    assert_allclose(sc.data.lat.deg, [0, 10.003028])

    pc2 = PixCoord.from_sky(skycoord=sc, wcs=wcs)
    assert isinstance(pc2.x, np.ndarray)
    assert pc2.x.shape == (2,)
    assert not pc2.isscalar

    assert pc1 == pc2


def test_pixcoord_to_sky_array_2d(wcs):
    pc0 = PixCoord(x=[[17, 18]], y=[[8, 9]])
    pc1 = PixCoord(x=[[17, 17, 17], [18, 18, 18]], y=[[8, 8, 8], [9, 9, 9]])

    sc0 = pc0.to_sky(wcs=wcs)
    assert sc0.name == 'galactic'

    sc1 = pc1.to_sky(wcs=wcs)
    assert sc1.name == 'galactic'
    assert_allclose(sc1.data.lon.deg, [[0, 0, 0],
                                       [349.88094, 349.88094, 349.88094]])
    assert_allclose(sc1.data.lat.deg, [[0, 0, 0],
                                       [10.003028, 10.003028, 10.003028]])
    assert_allclose(sc0.data.lon.deg, [[0, 349.88094]])
    assert_allclose(sc0.data.lat.deg, [[0, 10.003028]])

    pc2 = PixCoord.from_sky(skycoord=sc0, wcs=wcs)
    assert isinstance(pc2.x, np.ndarray)
    assert pc2.x.shape == (1, 2)
    assert not pc2.isscalar

    assert pc0 == pc2


def test_pixcoord_separation_scalar():
    pc1 = PixCoord(x=1, y=2)
    pc2 = PixCoord(x=4, y=6)
    sep = pc1.separation(pc2)
    assert_allclose(sep, 5)


def test_pixcoord_separation_array_1d():
    pc1 = PixCoord(x=[1, 1], y=[2, 2])
    pc2 = PixCoord(x=[4, 4], y=[6, 6])
    sep = pc1.separation(pc2)
    assert_allclose(sep, [5, 5])


def test_pixcoord_separation_array_2d():
    pc1 = PixCoord(x=[[1, 1]], y=[[2, 2]])
    pc2 = PixCoord(x=[[4, 4]], y=[[6, 6]])
    sep = pc1.separation(pc2)
    assert isinstance(sep, np.ndarray)
    assert sep.shape == (1, 2)
    assert_allclose(sep, [[5, 5]])


def test_equality():
    arr = np.array([1, 2])
    pc1 = PixCoord(arr[0], arr[1])
    pc2 = PixCoord(arr[0] + 0.0000001, arr[1])

    assert pc1 != arr
    assert pc1 == PixCoord(arr[0], arr[1])
    assert pc1 == pc2

    pc3 = PixCoord([[1, 2, 3], [4, 5, 6]], [[11, 12, 13], [14, 15, 16]])
    pc4 = PixCoord([[1, 2, 3], [4, 5, 6]],
                   [[11.0000002, 12, 13], [14, 15, 16]])

    assert pc3 == pc4
    assert pc3 == PixCoord([[1, 2, 3], [4, 5, 6]],
                           [[11, 12, 13], [14, 15, 16]])


def test_pixcoord_xy():
    arr = np.array([1, 2])
    pc = PixCoord(arr[0], arr[1])
    assert pc.xy[0] == pc.x
    assert pc.xy[1] == pc.y


def test_pixcoord_rotate_scalar():
    point = PixCoord(3, 4)
    center = PixCoord(2, 3)
    point = point.rotate(center, 90 * u.deg)
    assert_allclose(point.xy, (1, 4))


def test_pixcoord_rotate_array():
    point = PixCoord([3, 3], [4, 4])
    center = PixCoord(2, 3)
    point = point.rotate(center, 90 * u.deg)
    assert_allclose(point.xy, ([1, 1], [4, 4]))


def test_pixcoord_addition():
    point1 = PixCoord([3, 3], [4, 4])
    point2 = PixCoord([5, 7], [6, 1])
    center = PixCoord(2, 3)

    p1 = point1 + point2
    p2 = point2 + point1
    assert p1 == p2
    assert_equal(p1.x, (8, 10))
    assert_equal(p1.y, (10, 5))

    p3 = point1 + center
    p4 = center + point1
    assert p3 == p4
    assert_equal(p3.x, (5, 5))
    assert_equal(p3.y, (7, 7))

    p5 = point1 + point1 + point1
    assert_equal(p5.x, point1.x * 3)
    assert_equal(p5.y, point1.y * 3)

    with pytest.raises(TypeError):
        point1 + 10

    with pytest.raises(ValueError):
        point1 + PixCoord([1, 1, 1], [2, 2, 2])


def test_pixcoord_subtraction():
    point1 = PixCoord([3, 3], [4, 4])
    point2 = PixCoord([5, 7], [6, 1])
    center = PixCoord(2, 3)

    p1 = point2 - point1
    assert_equal(p1.x, (2, 4))
    assert_equal(p1.y, (2, -3))

    p2 = point1 - center
    assert_equal(p2.x, (1, 1))
    assert_equal(p2.y, (1, 1))

    with pytest.raises(TypeError):
        point1 - 10

    with pytest.raises(ValueError):
        point1 - PixCoord([1, 1, 1], [2, 2, 2])
