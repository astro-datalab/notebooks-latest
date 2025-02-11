#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.query
"""
from functools import partial

from contextlib import ExitStack

from os import listdir

import pytest

import numpy as np

import platform

from pyvo.dal.query import DALService, DALQuery, DALResults, Record
from pyvo.dal.exceptions import DALServiceError, DALQueryError, DALFormatError, DALOverflowWarning
from pyvo.version import version

from astropy.table import Table, QTable
from astropy.io.votable.tree import VOTableFile

try:
    # Workaround astropy deprecation, remove try/except once >=6.0 is required
    from astropy.io.votable.tree import TableElement
except ImportError:
    from astropy.io.votable.tree import Table as TableElement

from astropy.io.fits import HDUList
from astropy.utils.data import get_pkg_data_contents

get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')


@pytest.fixture()
def register_mocks(mocker):
    with ExitStack() as stack:
        matchers = [
            stack.enter_context(mocker.register_uri(
                'GET', '//example.com/query/basic',
                content=get_pkg_data_contents('data/query/basic.xml')
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/missingtable',
                content=get_pkg_data_contents('data/query/missingtable.xml')
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/missingresource',
                content=get_pkg_data_contents('data/query/missingresource.xml')
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/missingcolumns',
                content=get_pkg_data_contents('data/query/missingcolumns.xml')
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/firstresource',
                content=get_pkg_data_contents('data/query/firstresource.xml')
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/rootinfo',
                content=get_pkg_data_contents('data/query/rootinfo.xml')
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/tableinfo',
                content=get_pkg_data_contents('data/query/tableinfo.xml')
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/dataset',
                content=get_pkg_data_contents('data/query/dataset.xml')
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/querydata/image.fits',
                content=get_pkg_data_contents('data/querydata/image.fits')
            )),
            # mocker.register_uri(
            #     'GET', 'http://example.com/querydata/votable.xml',
            #     content=get_pkg_data_contents('data/querydata/votable.xml')
            # ),
            # mocker.register_uri(
            #     'GET', 'http://example.com/querydata/votable-datalink.xml',
            #     content=get_pkg_data_contents('data/querydata/votable-datalink.xml')
            # ),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/nonexistant',
                text='Not Found', status_code=404
            )),
            stack.enter_context(mocker.register_uri(
                'GET', '//example.com/query/errornous',
                text='Error', status_code=500
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/errorstatus',
                content=get_pkg_data_contents('data/query/errorstatus.xml')
            )),
            stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/query/overflowstatus',
                content=get_pkg_data_contents('data/query/overflowstatus.xml')
            )),
        ]

        def verbosetest_callback(request, context):
            assert 'VERBOSE' in request.qs and '1' in request.qs['VERBOSE']
            return get_pkg_data_contents('data/query/basic.xml')

        matchers.append(stack.enter_context(mocker.register_uri(
            'GET', 'http://example.com/query/verbosetest',
            content=verbosetest_callback
        )))

        def useragent_callback(request, context):
            assert 'User-Agent' in request.headers
            assert request.headers['User-Agent'] == 'pyVO/{} Python/{} ({})'.format(
                version, platform.python_version(), platform.system())
            return get_pkg_data_contents('data/query/basic.xml')

        matchers.append(stack.enter_context(mocker.register_uri(
            'GET', 'http://example.com/query/useragent',
            content=useragent_callback
        )))

        yield matchers


def _test_results(results):
    """Regression test result columns for correctnes"""
    assert len(results) == 3

    assert results['1', 0] == 23
    assert results['1', 1] == 42
    assert results['1', 2] == 1337

    truth = 'Illuminatus'
    assert results['2', 0] == truth
    truth = "Don't panic, and always carry a towel"
    assert results['2', 1] == truth
    truth = 'Elite'
    assert results['2', 2] == truth


def _test_records(records):
    """ Regression test dal records for correctness"""
    assert len(records) == 3

    assert all([isinstance(record, Record) for record in records])

    assert records[0]['1'] == 23

    truth = 'Illuminatus'
    assert records[0]['2'] == truth

    assert records[1]['1'] == 42
    truth = "Don't panic, and always carry a towel"
    assert records[1]['2'] == truth

    assert records[2]['1'] == 1337
    truth = 'Elite'
    assert records[2]['2'] == truth


@pytest.fixture
def url():
    return "http://example.com/query/basic"


@pytest.fixture
def description():
    return "An example service."


@pytest.fixture
def basic_service(url, description):
    return DALService(url, capability_description=description)


@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
@pytest.mark.usefixtures('register_mocks')
class TestDALService:

    def test_init(self, url, description, basic_service):
        """Test if baseurl and description are passed correctly"""
        assert basic_service.baseurl == url
        assert basic_service.capability_description == "An example service."

    def test__repr__(self, basic_service):
        assert str(basic_service) == (f"DALService(baseurl : '{basic_service.baseurl}',"
                                      f" description : '{basic_service.capability_description}')")

    def test_search(self):
        """
        Test (in conjunction with mocker) that parameters arrive serverside,
        while also ensuring data consistency
        """
        service = DALService('http://example.com/query/verbosetest')
        dalresults = service.search(VERBOSE=1)

        _test_results(dalresults)
        _test_records(dalresults)

    def test_useragent(self):
        service = DALService('http://example.com/query/useragent')
        service.search()

    def test_http_exception_404(self):
        service = DALService('http://example.com/query/nonexistant')

        try:
            service.search()
        except DALServiceError as exc:
            assert exc.code == 404
        else:
            assert False

    def test_http_exception_500(self):
        service = DALService('http://example.com/query/errornous')

        try:
            service.search()
        except DALServiceError as exc:
            assert exc.code == 500
        else:
            assert False

    def test_query_exception(self):
        service = DALService('http://example.com/query/errorstatus')

        with pytest.raises(DALQueryError):
            service.search()

    def test_query_warning(self):
        service = DALService('http://example.com/query/overflowstatus')

        with pytest.warns(DALOverflowWarning):
            service.search()

    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W53")
    def test_format_exception(self):
        with pytest.raises(DALFormatError):
            service = DALService('http://example.com/query/missingtable')
            service.search()

        with pytest.raises(DALFormatError):
            service = DALService('http://example.com/query/missingresource')
            service.search()

        with pytest.raises(DALFormatError):
            service = DALService('http://example.com/query/missingcolumns')
            service.search()


@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
@pytest.mark.usefixtures('register_mocks')
class TestDALQuery:
    def test_url(self):
        queries = (
            DALQuery('http://example.com/query/basic'),
            DALQuery(b'http://example.com/query/basic'),
        )

        assert all(
            q.queryurl == 'http://example.com/query/basic' for q in queries
        )

    def test_params(self):
        query = DALQuery(
            'http://example.com/query/basic', verbose=1, foo='BAR')

        assert query['VERBOSE'] == 1
        assert query['FOO'] == 'BAR'

    def test_execute(self):
        query = DALQuery('http://example.com/query/basic')
        dalresults = query.execute()

        assert dalresults.queryurl == 'http://example.com/query/basic'

        _test_results(dalresults)
        _test_records(dalresults)

    def test_execute_raw(self):
        query = DALQuery('http://example.com/query/basic')
        raw = query.execute_raw()

        assert raw.startswith(b'<?xml')
        assert raw.strip().endswith(b'</VOTABLE>')


@pytest.mark.filterwarnings('ignore::astropy.io.votable.exceptions.W03')
@pytest.mark.filterwarnings('ignore::astropy.io.votable.exceptions.W06')
@pytest.mark.usefixtures('register_mocks')
class TestDALResults:
    def test_init(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/basic')

        assert dalresults.queryurl == 'http://example.com/query/basic'
        assert isinstance(dalresults.votable, VOTableFile)
        assert isinstance(dalresults.resultstable, TableElement)

        assert dalresults.fieldnames == ('1', '2')
        assert (
            dalresults.fielddescs[0].name, dalresults.fielddescs[1].name
        ) == ('1', '2')

        assert dalresults.status == ('OK', 'OK')

    def test_from_result_url(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/basic')
        assert dalresults.status == ('OK', 'OK')

    def test_init_errorstatus(self):
        with pytest.raises(DALQueryError):
            DALResults.from_result_url('http://example.com/query/errorstatus')

    def test_init_overflowstatus(self):
        with pytest.warns(DALOverflowWarning):
            DALResults.from_result_url('http://example.com/query/overflowstatus')

    def test_init_missingtable(self):
        with pytest.raises(DALFormatError):
            DALResults.from_result_url('http://example.com/query/missingtable')

    @pytest.mark.filterwarnings('ignore::astropy.io.votable.exceptions.W53')
    def test_init_missingresource(self):
        with pytest.raises(DALFormatError):
            DALResults.from_result_url(
                'http://example.com/query/missingresource')

    def test_init_missingcolumns(self):
        with pytest.raises(DALFormatError):
            DALResults.from_result_url(
                'http://example.com/query/missingcolumns')

    def test_init_firstresource(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/firstresource')
        assert dalresults.status == ('OK', 'OK')

    def test_init_tableinfo(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/tableinfo')
        assert dalresults.status == ('OK', 'OK')

    def test_init_rootinfo(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/rootinfo')
        assert dalresults.status == ('OK', 'OK')

    def test_repr(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/basic')

        assert repr(dalresults)[0:26] == "<DALResultsTable length=3>"

    def test_iter(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/basic')

        records = list(iter(dalresults))

        _test_results(dalresults)
        _test_records(records)

    def test_dataconsistency(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/basic')

        assert isinstance(dalresults['1'], np.ndarray)
        assert isinstance(dalresults['2'], np.ndarray)

        _test_results(dalresults)
        _test_records(dalresults)

    def test_table_conversion(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/basic')

        assert isinstance(dalresults.to_table(), Table)
        assert isinstance(dalresults.to_qtable(), QTable)
        assert len(dalresults) == len(dalresults.to_table())
        assert len(dalresults) == len(dalresults.to_qtable())

    def test_id_over_name(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/basic')

        assert isinstance(dalresults['_1'], np.ndarray)
        assert isinstance(dalresults['_2'], np.ndarray)

        table = dalresults.to_table()
        with pytest.raises(KeyError):
            assert table['_1']
        with pytest.raises(KeyError):
            assert table['_2']

    def test_nosuchcolumn(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/basic')

        with pytest.raises(KeyError):
            dalresults['nosuchcolumn']

        with pytest.raises(KeyError):
            dalresults.getdesc('nosuchcolumn')

    def test_columnaliases(self):
        dalresults = DALResults.from_result_url(
            'http://example.com/query/basic')

        assert dalresults.fieldname_with_ucd('foo') == '1'
        assert dalresults.fieldname_with_ucd('bar') == '1'

        assert dalresults.fieldname_with_utype('foobar') == '2'

        assert dalresults.fieldname_with_ucd('baz') is None
        assert dalresults.fieldname_with_utype('foobaz') is None


@pytest.mark.filterwarnings('ignore::astropy.io.votable.exceptions.W03')
@pytest.mark.filterwarnings('ignore::astropy.io.votable.exceptions.W06')
@pytest.mark.usefixtures('register_mocks')
class TestRecord:
    def test_itemaccess(self):
        record = DALResults.from_result_url(
            'http://example.com/query/basic')[0]

        assert record['1'] == 23
        truth = 'Illuminatus'
        assert record['2'] == truth

        assert record['_1'] == 23
        assert record['_2'] == truth

    def test_nosuchcolumn(self):
        record = DALResults.from_result_url(
            'http://example.com/query/basic')[0]

        with pytest.raises(KeyError):
            record['nosuchcolumn']

    def test_iter(self):
        record = DALResults.from_result_url(
            'http://example.com/query/basic')[0]

        record = list(iter(record))

        assert record[0] == '1'
        assert record[1] == '2'

    def test_len(self):
        record = DALResults.from_result_url(
            'http://example.com/query/basic')[0]

        assert len(record) == 2

    def test_repr(self):
        record = DALResults.from_result_url(
            'http://example.com/query/basic')[0]
        truth = 'Illuminatus'
        assert repr(record) == repr((23, truth))

    def test_get(self):
        record = DALResults.from_result_url(
            'http://example.com/query/basic')[0]

        assert record.get('2', decode=True) == 'Illuminatus'

    def test_columnaliases(self):
        record = DALResults.from_result_url(
            'http://example.com/query/basic')[0]

        assert record.getbyucd('foo') == 23
        assert record.getbyucd('bar') == 23

        truth = 'Illuminatus'
        assert record.getbyutype('foobar') == truth

        record.getbyucd('baz') is None
        record.getbyutype('foobaz') is None

    def test_datasets(self):
        records = DALResults.from_result_url(
            'http://example.com/query/dataset')

        record = records[0]
        assert record.getdataurl() == 'http://example.com/querydata/image.fits'
        dataset = record.getdataset()
        HDUList.fromstring(dataset.read())

    def test_nodataset(self):
        record = DALResults.from_result_url(
            'http://example.com/query/basic')[0]

        assert record.getdataurl() is None

        with pytest.raises(KeyError):
            record.getdataset().read()

    def test_cachedataset(self, tmpdir):
        tmpdir = str(tmpdir)

        record = DALResults.from_result_url(
            'http://example.com/query/dataset')[0]

        record.cachedataset(dir=tmpdir)

        assert "dataset.dat" in listdir(tmpdir)


class TestUpload:
    pass
