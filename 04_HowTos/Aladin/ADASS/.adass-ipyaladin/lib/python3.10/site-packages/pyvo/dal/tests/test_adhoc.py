#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.adhoc
"""
import datetime
from astropy import units as u
from astropy.time import Time
import pytest

from pyvo.dal.adhoc import AxisParamMixin, SodaQuery


def test_pos():
    class TestClass(dict, AxisParamMixin):
        pass
    test_obj = TestClass()
    test_obj.pos.add((1, 2, 3) * u.deg)
    assert len(test_obj._pos) == 1
    assert test_obj['POS'] == ['CIRCLE 1.0 2.0 3.0']

    test_obj.pos.add((1, 2, 3, 4))
    assert len(test_obj._pos) == 2
    assert test_obj['POS'] == ['CIRCLE 1.0 2.0 3.0', 'RANGE 1.0 2.0 3.0 4.0']

    # duplicates are ignored
    test_obj.pos.add((1, 2, 3))
    assert len(test_obj._pos) == 2
    assert test_obj['POS'] == ['CIRCLE 1.0 2.0 3.0', 'RANGE 1.0 2.0 3.0 4.0']

    # polygon
    test_obj.pos.add((1, 2, 3, 4, 5, 6))
    assert len(test_obj._pos) == 3
    assert test_obj['POS'] == ['CIRCLE 1.0 2.0 3.0', 'RANGE 1.0 2.0 3.0 4.0',
                               'POLYGON 1.0 2.0 3.0 4.0 5.0 6.0']

    # deletes
    test_obj.pos.remove((1, 2, 3, 4))
    assert len(test_obj._pos) == 2
    assert test_obj['POS'] == ['CIRCLE 1.0 2.0 3.0',
                               'POLYGON 1.0 2.0 3.0 4.0 5.0 6.0']

    # test borders
    test_obj.pos.discard((1, 2, 3) * u.deg)
    test_obj.pos.discard((1, 2, 3, 4, 5, 6))
    assert (len(test_obj._pos) == 0)
    test_obj.pos.add((0, 90, 90))
    assert len(test_obj._pos) == 1
    assert test_obj['POS'] == ['CIRCLE 0.0 90.0 90.0']
    test_obj.pos.pop()
    test_obj.pos.add((360, -90, 1))
    assert len(test_obj._pos) == 1
    assert test_obj['POS'] == ['CIRCLE 360.0 -90.0 1.0']
    test_obj.pos.pop()
    test_obj.pos.add((0, 360, -90, 90))
    assert len(test_obj._pos) == 1
    assert test_obj['POS'] == ['RANGE 0.0 360.0 -90.0 90.0']
    test_obj.pos.pop()
    test_obj.pos.add((0, 0, 180, 90, 270, -90))
    assert len(test_obj._pos) == 1
    assert test_obj['POS'] == ['POLYGON 0.0 0.0 180.0 90.0 270.0 -90.0']

    # errors
    test_obj.pos.pop()
    with pytest.raises(ValueError):
        test_obj.pos.add(('A', 2, 3))
    with pytest.raises(ValueError):
        test_obj.pos.add((-2, 7, 3))
    with pytest.raises(ValueError):
        test_obj.pos.add((3, 99, 3))
    with pytest.raises(ValueError):
        test_obj.pos.add((2, 7, 91))
    with pytest.raises(ValueError):
        test_obj.pos.add((-1, 7, 3, 4))
    with pytest.raises(ValueError):
        test_obj.pos.add((2, 1, 3, 4))
    with pytest.raises(ValueError):
        test_obj.pos.add((1, 2, 4, 3))
    with pytest.raises(ValueError):
        test_obj.pos.add((-2, 7, 5, 9, 10, 10))
    with pytest.raises(ValueError):
        test_obj.pos.add((2, 99, 5, 9, 10, 10))
    with pytest.raises(ValueError):
        test_obj.pos.add((1, 2, 3, 4, 5, 6, 7))


def test_band():
    class TestClass(dict, AxisParamMixin):
        pass
    test_obj = TestClass()
    assert not hasattr(test_obj, '_band')
    test_obj.band.add(33)
    assert 33 in test_obj.band
    assert test_obj['BAND'] == ['33.0 33.0']
    test_obj.band.add((50 * u.meter, 500))
    assert 33 in test_obj.band
    assert (50 * u.meter, 500) in test_obj.band
    assert test_obj['BAND'] == ['33.0 33.0', '50.0 500.0']
    test_obj.band.discard(33)
    assert (50 * u.meter, 500) in test_obj.band
    assert test_obj['BAND'] == ['50.0 500.0']
    test_obj.band.pop()
    assert not test_obj.band
    assert not test_obj['BAND']
    test_obj.band.add((float('-inf'), 33))
    assert (float('-inf'), 33) in test_obj.band
    assert test_obj.band.dal == ['-inf 33.0']
    test_obj.band.clear()
    test_obj.band.add((33, float('inf')))
    assert (33, float('inf')) in test_obj.band
    assert test_obj.band.dal == ['33.0 inf']
    test_obj.clear()

    # error cases
    with pytest.raises(ValueError):
        test_obj.band.add(())
    with pytest.raises(ValueError):
        test_obj.band.add((1, 2, 3))
    with pytest.raises(ValueError):
        test_obj.band.add(('INVALID', 6))
    with pytest.raises(ValueError):
        test_obj.band.add((3, 1))


def test_time():
    class TestClass(dict, AxisParamMixin):
        pass
    test_obj = TestClass()
    assert not hasattr(test_obj, '_time')
    now = Time(datetime.datetime.now(tz=datetime.timezone.utc))
    test_obj.time.add(now)
    assert now in test_obj.time
    assert test_obj['TIME'] == ['{now} {now}'.format(now=now.mjd)]
    min_time = '2010-01-01T00:00:00.000Z'
    max_time = '2010-01-01T01:00:00.000Z'
    test_obj.time.add((min_time, max_time))
    assert now in test_obj.time
    assert (min_time, max_time) in test_obj.time
    assert test_obj['TIME'] == ['{now} {now}'.format(now=now.mjd),
                                '{min} {max}'.format(min=Time(min_time).mjd,
                                                     max=Time(max_time).mjd)]
    test_obj.time.discard(now)
    assert (min_time, max_time) in test_obj.time
    assert test_obj['TIME'] == ['{min} {max}'.format(min=Time(min_time).mjd,
                                                     max=Time(max_time).mjd)]
    test_obj.time.pop()
    assert not test_obj.time
    assert not test_obj['TIME']

    # error cases
    with pytest.raises(ValueError):
        test_obj.time.add([])
    with pytest.raises(ValueError):
        test_obj.time.add([now, min_time, max_time])
    with pytest.raises(ValueError):
        test_obj.time.add(['INVALID'])
    with pytest.raises(ValueError):
        test_obj.time.add([max_time, min_time])


def test_pol():
    class TestClass(dict, AxisParamMixin):
        pass
    test_obj = TestClass()
    assert not hasattr(test_obj, '_pol')
    test_obj.pol.add('YY')
    assert 'YY' in test_obj.pol
    assert test_obj['POL'] == ['YY']
    test_obj.pol.add('POLI')
    assert 'YY' in test_obj.pol
    assert 'POLI' in test_obj.pol
    assert test_obj['POL'] == ['YY', 'POLI']
    # test duplicate
    test_obj.pol.add('POLI')
    assert 'YY' in test_obj.pol
    assert 'POLI' in test_obj.pol
    assert test_obj['POL'] == ['YY', 'POLI']
    test_obj.pol.remove('YY')
    assert 'POLI' in test_obj.pol
    assert test_obj['POL'] == ['POLI']
    test_obj.pol.pop()
    assert not test_obj._pol
    assert not test_obj['POL']

    # error cases
    with pytest.raises(ValueError):
        test_obj.pol.add(None)
    with pytest.raises(ValueError):
        test_obj.pol.add(['INVALID'])


def test_soda_query():
    test_obj = SodaQuery(baseurl='some/url')
    test_obj.circle = (2, 3, 5)
    assert test_obj._circle == (2, 3, 5)
    assert test_obj['CIRCLE'] == '2.0 3.0 5.0'
    assert test_obj._circle
    assert not hasattr(test_obj, '_polygon')
    assert not hasattr(test_obj, '_range')

    test_obj.range = (8, 9, 3, 4) * u.deg
    assert test_obj['POS'] == 'RANGE 8.0 9.0 3.0 4.0'
    assert test_obj._range is not None
    assert not hasattr(test_obj, '_polygon')
    assert not hasattr(test_obj, '_circle')

    test_obj.polygon = (1, 2, 3, 4, 5, 6)
    assert test_obj['POLYGON'] == '1.0 2.0 3.0 4.0 5.0 6.0'
    assert test_obj._polygon
    assert not hasattr(test_obj, '_range')
    assert not hasattr(test_obj, '_circle')

    del test_obj.polygon
    assert not hasattr(test_obj, '_polygon')
    assert not hasattr(test_obj, '_circle')
    assert not hasattr(test_obj, '_range')

    # error cases
    with pytest.raises(ValueError):
        test_obj.circle = ('A', 1, 2)
    with pytest.raises(ValueError):
        test_obj.circle = (1, 1, 2, 2)
    with pytest.raises(ValueError):
        test_obj.circle = (-1, 1, 2)
    with pytest.raises(ValueError):
        test_obj.circle = (1, 99, 2)
    with pytest.raises(ValueError):
        test_obj.circle = (1, 1, 91)
    with pytest.raises(ValueError):
        test_obj.range = (1, 2, 3)
    with pytest.raises(ValueError):
        test_obj.range = (2, 1, 3, 4)
    with pytest.raises(ValueError):
        test_obj.range = (1, 2, 4, 3)
    with pytest.raises(ValueError):
        test_obj.range = (-1, 2, 3, 4)
    with pytest.raises(ValueError):
        test_obj.range = (2, 1000, 3, 4)
    with pytest.raises(ValueError):
        test_obj.range = (1, 1, -91, 4)
    with pytest.raises(ValueError):
        test_obj.range = (1, 1, 3, 92)
    with pytest.raises(ValueError):
        test_obj.polygon = (1, 2, 3, 4)
    with pytest.raises(ValueError):
        test_obj.polygon = (2, 1, 3, 4, 5, 6, 7)
