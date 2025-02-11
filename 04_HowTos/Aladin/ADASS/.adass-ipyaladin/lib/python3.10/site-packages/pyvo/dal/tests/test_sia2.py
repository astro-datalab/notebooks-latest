#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.sia
"""
from functools import partial
import re
import requests_mock

import pytest

from pyvo.dal.sia2 import search, SIA2Service, SIA2Query, SIAService, SIAQuery

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.utils.data import get_pkg_data_contents
from astropy.utils.exceptions import AstropyDeprecationWarning

get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')

sia_re = re.compile('https://example.com/sia/v2query*')
capabilities_url = 'https://example.com/sia/capabilities'


@pytest.fixture()
def sia(mocker):
    with mocker.register_uri(
        'GET', sia_re, content=get_pkg_data_contents('data/sia2/dataset.xml')
    ) as matcher:
        yield matcher


@pytest.fixture()
def capabilities(mocker):
    with mocker.register_uri(
        'GET', capabilities_url,
            content=get_pkg_data_contents('data/sia2/capabilities.xml')
    ) as matcher:
        yield matcher


def _test_result(record):
    assert record.obs_collection == 'TEST'
    assert record.obs_id == 'TEST-DATASET'
    assert record.instrument_name == 'TEST-INSTR'
    assert record.facility_name == 'TEST-1.6m'


@pytest.mark.usefixtures('sia')
@pytest.mark.usefixtures('capabilities')
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_search():
    results = search('https://example.com/sia',
                     pos=(33.3 * u.deg, 4.2 * u.deg, 0.0166 * u.deg))
    result = results[0]

    _test_result(result)


class TestSIA2Service:

    def test_capabilities(self):
        # this tests the SIA2 capabilities with various combinations:

        with requests_mock.Mocker() as cm:
            cm.get('https://example.com/sia/capabilities',
                   content=get_pkg_data_contents('data/sia2/capabilities.xml'))
            cm.get('https://example.com/sia-basicauth/capabilities',
                   content=get_pkg_data_contents(
                       'data/sia2/capabilities-basicauth.xml'))
            cm.get('https://example.com/sia-newformat/capabilities',
                   content=get_pkg_data_contents(
                       'data/sia2/capabilities-newformat.xml'))
            cm.get('https://example.com/sia-priv/capabilities',
                   content=get_pkg_data_contents(
                       'data/sia2/capabilities-priv.xml')),
            cm.get('https://example.com/sia/myquery/capabilities',
                   content=get_pkg_data_contents('data/sia2/capabilities.xml'))

            # multiple interfaces with single security method each and
            # anonymous access.
            service = SIA2Service('https://example.com/sia')
            assert service.query_ep == 'https://example.com/sia/v2query'

            # one interface with multiple security methods
            service = SIA2Service('https://example.com/sia-newformat')
            assert service.query_ep == 'https://example.com/sia/v2query'

            # multiple interfaces with single security method each (no anon)
            service = SIA2Service('https://example.com/sia-priv')
            assert service.query_ep == 'https://example.com/sia/v2query'

            # any access point will be valid even when it contains query params
            service = SIA2Service('https://example.com/sia/myquery?param=1')
            assert service.query_ep == 'https://example.com/sia/v2query'

            # capabilities checking is bypassed all together with the
            # check_baseurl=False flag
            service = SIA2Service('https://example.com/sia/myquery?param=1&',
                                  check_baseurl=False)
            assert service.query_ep == 'https://example.com/sia/myquery?param=1'

    POSITIONS = [(2, 4, 0.0166 * u.deg),
                 (12, 12.5, 34, 36),
                 (12.0 * u.deg, 34.0 * u.deg,
                  14.0 * u.deg, 35.0 * u.deg,
                  14.0 * u.deg, 36.0 * u.deg,
                  12.0 * u.deg, 35.0 * u.deg),
                 (SkyCoord(2, 4, unit='deg'), 0.166 * u.deg)]

    @pytest.mark.usefixtures('sia')
    @pytest.mark.usefixtures('capabilities')
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W42")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W49")
    @pytest.mark.parametrize("position", POSITIONS)
    def test_search_scalar(self, position):
        service = SIA2Service('https://example.com/sia')

        results = service.search(pos=position)
        result = results[0]
        _test_result(result)

    @pytest.mark.usefixtures('sia')
    @pytest.mark.usefixtures('capabilities')
    def test_search_vector(self, pos=POSITIONS):
        service = SIA2Service('https://example.com/sia')
        results = service.search(pos=pos)
        result = results[0]
        _test_result(result)

    @pytest.mark.usefixtures('sia')
    @pytest.mark.usefixtures('capabilities')
    def test_search_deprecation(self, pos=POSITIONS):
        # test the deprecation
        with pytest.warns(AstropyDeprecationWarning):
            deprecated_service = SIAService('https://example.com/sia')
            deprecated_results = deprecated_service.search(pos=pos)
            result = deprecated_results[0]
            _test_result(result)


class TestSIA2Query():

    def test_query(self):
        query = SIA2Query('someurl')
        query.field_of_view.add((10, 20))
        assert query['FOV'] == ['10.0 20.0']
        query.field_of_view.add((1 * u.rad, 60))
        assert query['FOV'] == ['10.0 20.0', '57.29577951308232 60.0']

        query.spatial_resolution.add((1 * u.arcsec, 2))
        assert query['SPATRES'] == ['1.0 2.0']

        query.spectral_resolving_power.add((3, 5))
        assert query['SPECRP'] == ['3 5']

        query.exptime.add((25, 50))
        assert query['EXPTIME'] == ['25.0 50.0']

        query.timeres.add((1, 3))
        assert query['TIMERES'] == ['1.0 3.0']

        query.publisher_did.add('ID1')
        query.publisher_did.add('ID2')
        assert query['ID'] == ['ID1', 'ID2']

        query.facility.add('TEL1')
        assert query['FACILITY'] == ['TEL1']

        query.collection.add('ABC')
        query.collection.add('EFG')
        assert query['COLLECTION'] == ['ABC', 'EFG']

        query.instrument.add('INST1')
        assert query['INSTRUMENT'] == ['INST1']

        query.data_type.add('TYPEA')
        assert query['DPTYPE'] == ['TYPEA']

        query.calib_level.add(0)
        query.calib_level.add(1)
        assert query['CALIB'] == ['0', '1']

        query.target_name.add('TARGET1')
        assert query['TARGET'] == ['TARGET1']

        query.res_format.add('pdf')
        assert query['FORMAT'] == ['pdf']

        query.maxrec = 1000
        assert query['MAXREC'] == '1000'

        query = SIA2Query('someurl', custom_param=23)
        assert query['custom_param'] == ['23']

        query['custom_param'].append('-Inf 0')
        assert query['custom_param'] == ['23', '-Inf 0']

        query = SIA2Query('someurl', custom_param=[('-Inf', 0), (2, '+Inf')])
        assert query['custom_param'] == ['-Inf 0', '2 +Inf']

        with pytest.warns(AstropyDeprecationWarning):
            deprecated_query = SIAQuery('someurl')
            deprecated_query.field_of_view.add((10, 20))
            assert deprecated_query['FOV'] == ['10.0 20.0']


def test_variable_deprecation():
    # Test this while we are in the deprecation period, as the variable is durectly
    # used at least by astroquery.alma
    with pytest.warns(AstropyDeprecationWarning):
        from pyvo.dal.sia2 import SIA_PARAMETERS_DESC
        assert SIA_PARAMETERS_DESC
