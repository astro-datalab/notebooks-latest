# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.tap
"""
from functools import partial
from contextlib import ExitStack
import datetime
import re
from io import BytesIO
from urllib.parse import parse_qsl
import tempfile

import pytest
import requests_mock

from pyvo.dal.tap import escape, search, AsyncTAPJob, TAPService
from pyvo.dal import DALQueryError, DALServiceError

from pyvo.io.uws import JobFile
from pyvo.io.uws.tree import Parameter, Result, ErrorSummary, Message
from pyvo.io.vosi.exceptions import VOSIError
from pyvo.utils import prototype

from astropy.time import Time, TimeDelta

from astropy.utils.data import get_pkg_data_contents

get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')

job_re_path_full = re.compile('^http://example.com/tap/async/([0-9]+)')
job_re_path = re.compile('^/tap/async/([0-9]+)')

job_re_phase_full = re.compile('^http://example.com/tap/async/([0-9]+)/phase')

job_re_parameters_full = re.compile(
    '^http://example.com/tap/async/([0-9]+)/parameters')

job_re_result_full = re.compile(
    '^http://example.com/tap/async/([0-9]+)/results/result')


def _test_image_results(results):
    assert len(results) == 10


@pytest.fixture()
def sync_fixture(mocker):
    def callback(request, context):
        return get_pkg_data_contents('data/tap/obscore-image.xml')

    with mocker.register_uri(
        'POST', 'http://example.com/tap/sync', content=callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def create_fixture(mocker):
    def match_request(request):
        data = request.text.read()
        if b'VOSITable' in data:
            assert request.headers['Content-Type'] == 'text/xml', 'Wrong file format'
        elif b'VOTable' in data:
            assert request.headers['Content-Type'] == \
                'application/x-votable+xml', 'Wrong file format'
        else:
            assert False, 'BUG'
        return True

    with mocker.register_uri(
        'PUT', 'https://example.com/tap/tables/abc',
        additional_matcher=match_request, status_code=201
    ) as matcher:
        yield matcher


@pytest.fixture()
def delete_fixture(mocker):
    with mocker.register_uri(
        'DELETE', 'https://example.com/tap/tables/abc', status_code=200,
    ) as matcher:
        yield matcher


@pytest.fixture()
def load_fixture(mocker):
    def match_request(request):
        data = request.text.read()
        if b',' in data:
            assert request.headers['Content-Type'] == 'text/csv', 'Wrong file format'
        elif b'\t' in data:
            assert request.headers['Content-Type'] == \
                'text/tab-separated-values', 'Wrong file format'
        elif b'FITSTable' in data:
            assert request.headers['Content-Type'] == \
                'application/fits', 'Wrong file format'
        else:
            assert False, 'BUG'
        return True

    with mocker.register_uri(
        'POST', 'https://example.com/tap/load/abc',
        additional_matcher=match_request, status_code=200,
    ) as matcher:
        yield matcher


def get_index_job(phase):
    return """<?xml version="1.0" encoding="UTF-8"?>
    <uws:job xmlns:uws="http://www.ivoa.net/xml/UWS/v1.0"
    xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1">
        <uws:jobId>v3njuz4k1ebpdb5q</uws:jobId>
        <uws:runId />
        <uws:ownerId>user</uws:ownerId>
        <uws:phase>{}</uws:phase>
        <uws:quote>2021-10-29T17:34:19.638</uws:quote>
        <uws:creationTime>2021-10-28T17:34:19.638</uws:creationTime>
        <uws:startTime xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true" />
        <uws:endTime xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:nil="true" />
        <uws:executionDuration>14400</uws:executionDuration>
        <uws:destruction>2021-11-04T17:34:19.638</uws:destruction>
        <uws:parameters>
            <uws:parameter id="index">article</uws:parameter>
            <uws:parameter id="table">cadcauthtest1.pyvoTestTable</uws:parameter>
            <uws:parameter id="unique">true</uws:parameter>
        </uws:parameters>
        <uws:results />
    </uws:job>""".format(phase).encode('utf-8')


class MockAsyncTAPServer:
    def __init__(self):
        self._jobs = dict()

    def validator(self, request):
        pass

    def use(self, mocker):
        with ExitStack() as stack:
            matchers = {
                'create': stack.enter_context(mocker.register_uri(
                    'POST', 'http://example.com/tap/async',
                    content=self.create
                )),
                'job': stack.enter_context(mocker.register_uri(
                    requests_mock.ANY, job_re_path_full, content=self.job
                )),
                'phase': stack.enter_context(mocker.register_uri(
                    requests_mock.ANY, job_re_phase_full, content=self.phase
                )),
                'parameters': stack.enter_context(mocker.register_uri(
                    requests_mock.ANY, job_re_parameters_full,
                    content=self.parameters
                )),
                'result': stack.enter_context(mocker.register_uri(
                    'GET', job_re_result_full, content=self.result
                )),
                'get_job': stack.enter_context(mocker.register_uri(
                    'GET', 'http://example.com/tap/async/111',
                    content=self.get_job
                )),
                'get_job_list': stack.enter_context(mocker.register_uri(
                    'GET', 'http://example.com/tap/async',
                    content=self.get_job_list
                ))
            }
            yield matchers

    def create(self, request, context):
        if request.method == 'GET':
            return self.get_job_list(request, context)
        self.validator(request)
        newid = max(list(self._jobs.keys()) or [0]) + 1
        data = dict(parse_qsl(request.body))

        job = JobFile()
        job.version = "1.1"
        job.jobid = newid
        if 'test_erroneus_submit.non_existent' in request.text:
            job.phase = 'ERROR'
            job._errorsummary = ErrorSummary()
            job.errorsummary.message = Message()
            job.errorsummary.message.content =\
                'test_erroneus_submit.non_existent not found'
        else:
            job.phase = 'PENDING'
        job.quote = Time.now() + TimeDelta(1, format='sec')
        job.creationtime = Time.now()
        job.executionduration = TimeDelta(3600, format='sec')
        job.destruction = Time.now() + TimeDelta(3600, format='sec')

        for key, value in data.items():
            param = Parameter(id=key)
            param.content = value
            job.parameters.append(param)

        context.status_code = 303
        context.reason = 'See other'
        context.headers['Location'] = (
            'http://example.com/tap/async/{}'.format(newid))

        self._jobs[newid] = job

    def job(self, request, context):
        self.validator(request)
        jobid = int(job_re_path.match(request.path).group(1))

        if request.method == 'GET':
            job = self._jobs[jobid]
            io = BytesIO()
            job.to_xml(io)
            return io.getvalue()
        elif request.method == 'POST':
            data = dict(parse_qsl(request.body))
            action = data.get('ACTION')

            if action == 'DELETE':
                del self._jobs[jobid]

    def phase(self, request, context):
        self.validator(request)
        jobid = int(job_re_path.match(request.path).group(1))

        if request.method == 'GET':
            phase = self._jobs[jobid].phase
            return phase
        elif request.method == 'POST':
            newphase = request.body.split('=')[-1]
            job = self._jobs[jobid]
            result = get_pkg_data_contents('data/tap/obscore-image.xml')

            if newphase == 'RUN':
                newphase = 'COMPLETED'
                result = Result(**{
                    'id': 'result',
                    'size': len(result),
                    'mime-type': 'application/x-votable+xml',
                    'xlink:href': (
                        'http://example.com/tap/async/{}/results/result'
                    ).format(jobid)
                })

                try:
                    job.results[0] = result
                except (IndexError, TypeError):
                    job.results.append(result)

            job.phase = newphase

    def parameters(self, request, context):
        self.validator(request)
        jobid = int(job_re_path.match(request.path).group(1))
        job = self._jobs[jobid]

        if request.method == 'GET':
            pass
        elif request.method == 'POST':
            data = dict(parse_qsl(request.body))

            if 'QUERY' in data:
                assert data['QUERY'] == 'SELECT TOP 42 * FROM ivoa.obsCore'
                for param in job.parameters:
                    if param.id_.lower() == 'query':
                        param.content = data['QUERY']
            if 'UPLOAD' in data:
                for param in job.parameters:
                    if param.id_.lower() == 'upload':
                        uploads1 = {data[0]: data[1] for data in [
                            data.split(',') for data
                            in data['UPLOAD'].split(';')
                        ]}

                        uploads2 = {data[0]: data[1] for data in [
                            data.split(',') for data
                            in param.content.split(';')
                        ]}

                        uploads1.update(uploads2)

                        param.content = ';'.join([
                            '{}={}'.format(key, value) for key, value
                            in uploads1.items()
                        ])

    def result(self, request, context):
        self.validator(request)
        return get_pkg_data_contents('data/tap/obscore-image.xml')

    def get_job(self, request, context):
        self.validator(request)
        jobid = int(job_re_path.match(request.path).group(1))

        job = JobFile()
        job.jobid = jobid
        job.phase = 'EXECUTING'
        job.ownerid = '222'
        job.creationtime = Time.now()
        io = BytesIO()
        job.to_xml(io)
        return io.getvalue()

    def _get_jobref_rep(self, jobid, phase, runid, ownerid, creation_time):
        doc = ('    <uws:jobref id="{}">\n'
               '        <uws:phase>{}</uws:phase>\n'
               '        <uws:runId>{}</uws:runId>\n'
               '        <uws:ownerId>{}</uws:ownerId>\n'
               '        <uws:creationTime>{}</uws:creationTime>\n'
               '    </uws:jobref>\n')
        return doc.format(jobid, phase, runid, ownerid, creation_time)

    def get_job_list(self, request, context):
        self.validator(request)
        fields = parse_qsl(request.query)
        phases = []
        last = None
        after = None
        for arg, val in fields:
            if arg == 'PHASE':
                phases.append(val)
            elif arg == 'LAST':
                last = int(val)
            elif arg == 'AFTER':
                after = val

        doc = '<?xml version="1.0" encoding="UTF-8"?>\n' +\
              '<uws:jobs xmlns:uws="http://www.ivoa.net/xml/UWS/v1.0" ' +\
              'xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1">\n'

        if phases:
            doc += self._get_jobref_rep('abc1', 'EXECUTING', 'def1', '21',
                                        '2018-12-20T00:23:15.79')
            doc += self._get_jobref_rep('abc2', 'EXECUTING', 'def2', '21',
                                        '2018-12-20T00:23:15.79')

        if after:
            doc += self._get_jobref_rep('abc3', 'EXECUTING', 'def3', '21',
                                        '2018-12-20T00:23:15.79')

        if last:
            doc += self._get_jobref_rep('abc4', 'EXECUTING', 'def4', '21',
                                        '2018-12-20T00:23:15.79')
            doc += self._get_jobref_rep('abc5', 'EXECUTING', 'def5', '21',
                                        '2018-12-20T00:23:15.79')
            doc += self._get_jobref_rep('abc6', 'EXECUTING', 'def6', '21',
                                        '2018-12-20T00:23:15.79')
        doc += '</uws:jobs>'

        return doc.encode('UTF-8')


@pytest.fixture()
def async_fixture(mocker):
    mock_server = MockAsyncTAPServer()
    yield from mock_server.use(mocker)


@pytest.fixture()
def tables(mocker):
    def callback_tables(request, context):
        return get_pkg_data_contents('data/tap/tables.xml')

    def callback_table1(request, context):
        return get_pkg_data_contents('data/tap/lazy-table1.xml')

    def callback_table2(request, context):
        return get_pkg_data_contents('data/tap/lazy-table2.xml')

    with ExitStack() as stack:
        matchers = {
            'tables': stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/tap/tables', content=callback_tables
            )),
            'table1': stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/tap/tables/test.table1',
                content=callback_table1
            )),
            'table2': stack.enter_context(mocker.register_uri(
                'GET', 'http://example.com/tap/tables/test.table2',
                content=callback_table2
            )),
        }

        yield matchers


@pytest.fixture()
def examples(mocker):
    def callback_examplesXHTML(request, context):
        uri = f"://{request.netloc}{request.path}"
        if uri == '://example.com/tap/examples':
            return get_pkg_data_contents('data/tap/examples.htm')
        elif uri == '://example.org/obscore-examples.html':
            return get_pkg_data_contents('data/tap/obscore-examples.html')
        else:
            assert False, f"Unexpected examples URI: {uri}"

    with mocker.register_uri(
        'GET', requests_mock.ANY, content=callback_examplesXHTML
    ) as matcher:
        yield matcher


@pytest.fixture()
def capabilities(mocker):
    def callback(request, context):
        return get_pkg_data_contents('data/tap/capabilities.xml')

    with mocker.register_uri(
        'GET', 'http://example.com/tap/capabilities', content=callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def tapservice(capabilities):
    """
    preferably use this fixture when you need a generic TAP service;
    it saves a bit of parsing overhead.

    (but of course make sure you don't modify it).
    """
    return TAPService('http://example.com/tap')


def test_escape():
    query = 'SELECT * FROM ivoa.obscore WHERE dataproduct_type = {}'
    query = query.format(escape("'image'"))
    assert query == (
        "SELECT * FROM ivoa.obscore WHERE dataproduct_type = ''image''")


@pytest.mark.usefixtures('sync_fixture')
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_search():
    results = search('http://example.com/tap', "SELECT * FROM ivoa.obscore")

    _test_image_results(results)


class TestTAPService:
    def test_init(self):
        url = 'http://example.com/tap'
        description = "An example TAP service."
        service = TAPService(url, capability_description=description)

        assert service.baseurl == url
        assert service.capability_description == description
        assert str(service) == f"""TAPService(baseurl : '{url}', description : '{description}')"""

    def _test_tables(self, table1, table2):
        assert table1.description == 'Lazy Test Table 1'
        assert table1.title == 'Test table 1'

        assert table2.description == 'Lazy Test Table 2'
        assert table2.title == 'Test table 2'

    @pytest.mark.usefixtures('tables')
    def test_tables(self):
        service = TAPService('http://example.com/tap')
        vositables = service.tables

        assert list(vositables.keys()) == ['test.table1', 'test.table2']

        assert "test.table1" in vositables
        assert "any.random.stuff" not in vositables

        table1, table2 = list(vositables)
        self._test_tables(table1, table2)

    def _test_examples(self, parsed_examples):
        assert len(parsed_examples) == 6
        assert "SELECT * FROM rosmaster" in parsed_examples[0]['QUERY']
        # the last query is from the continuation
        assert parsed_examples[-1]['QUERY'].startswith(
            "\nSELECT access_url, t_exptime, t_min FROM ivoa.obscore")

    @pytest.mark.usefixtures('examples')
    def test_examples(self):
        service = TAPService('http://example.com/tap')
        service_examples = service.examples

        self._test_examples(service_examples)

    @pytest.mark.usefixtures('capabilities')
    def test_maxrec(self):
        service = TAPService('http://example.com/tap')

        assert service.maxrec == 20000

    @pytest.mark.usefixtures('capabilities')
    def test_hardlimit(self):
        service = TAPService('http://example.com/tap')

        assert service.hardlimit == 10000000

    @pytest.mark.usefixtures('capabilities')
    def test_upload_methods(self):
        service = TAPService('http://example.com/tap')
        upload_methods = service.upload_methods

        assert upload_methods[0].ivo_id == (
            'ivo://ivoa.net/std/TAPRegExt#upload-https')
        assert upload_methods[1].ivo_id == (
            'ivo://ivoa.net/std/TAPRegExt#upload-ftp')
        assert upload_methods[2].ivo_id == (
            'ivo://ivoa.net/std/TAPRegExt#upload-inline')
        assert upload_methods[3].ivo_id == (
            'ivo://ivoa.net/std/TAPRegExt#upload-http')

    @pytest.mark.usefixtures('sync_fixture')
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
    def test_run_sync(self):
        service = TAPService('http://example.com/tap')
        results = service.run_sync("SELECT * FROM ivoa.obscore")
        _test_image_results(results)

    @pytest.mark.usefixtures('sync_fixture')
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
    def test_search(self):
        service = TAPService('http://example.com/tap')
        results = service.search("SELECT * FROM ivoa.obscore")
        _test_image_results(results)

    @pytest.mark.usefixtures('async_fixture')
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
    def test_run_async(self):
        service = TAPService('http://example.com/tap')
        results = service.run_async("SELECT * FROM ivoa.obscore")
        _test_image_results(results)

    @pytest.mark.usefixtures('async_fixture')
    def test_submit_job(self):
        service = TAPService('http://example.com/tap')
        job = service.submit_job("SELECT * FROM ivoa.obscore")

        assert job.url == 'http://example.com/tap/async/' + job.job_id
        assert job.phase == 'PENDING'
        assert job.execution_duration == TimeDelta(3600, format='sec')
        assert isinstance(job.destruction, Time)
        assert isinstance(job.quote, Time)
        assert job.query == "SELECT * FROM ivoa.obscore"

        job.run()
        job.wait()
        job.delete()

    @pytest.mark.usefixtures('async_fixture')
    def test_erroneus_submit_job(self):
        service = TAPService('http://example.com/tap')
        job = service.submit_job(
            "SELECT * FROM test_erroneus_submit.non_existent")
        with pytest.raises(DALQueryError) as e:
            job.raise_if_error()
        assert 'test_erroneus_submit.non_existent not found' in str(e)

    @pytest.mark.usefixtures('async_fixture')
    def test_submit_job_case(self):
        """Test using mixed case in the QUERY parameter to a job.

        DALI requires that query parameter names be case-insensitive, and
        some TAP servers reflect the input case into the job record, so the
        TAP client has to be prepared for any case for the QUERY parameter
        name.
        """
        service = TAPService('http://example.com/tap')

        # This has to be tested manually, bypassing the normal client layer,
        # in order to force a mixed-case parameter name.
        response = service._session.post(
            "http://example.com/tap/async",
            data={
                "REQUEST": "doQuery",
                "LANG": "ADQL",
                "quERy": "SELECT * FROM ivoa.obscore",
            }
        )
        response.raw.read = partial(response.raw.read, decode_content=True)
        job = AsyncTAPJob(response.url, session=service._session)

        assert job.url == 'http://example.com/tap/async/' + job.job_id
        assert job.query == "SELECT * FROM ivoa.obscore"

    @pytest.mark.usefixtures('async_fixture')
    def test_modify_job(self):
        service = TAPService('http://example.com/tap')
        job = service.submit_job(
            "SELECT * FROM ivoa.obscore", uploads={
                'one': 'http://example.com/uploads/one'
            })
        job.query = "SELECT TOP 42 * FROM ivoa.obsCore"
        job.upload(two='http://example.com/uploads/two')

        for parameter in job._job.parameters:
            if parameter.id_ == 'query':
                assert parameter.content == 'SELECT TOP 42 * FROM ivoa.obsCore'
                break
            elif parameter.id_ == 'upload':
                assert (
                    'one=http://example.com/uploads/one' in parameter.content)
                assert (
                    'two=http://example.com/uploads/two' in parameter.content)

    @pytest.mark.usefixtures('async_fixture')
    def test_get_job(self):
        service = TAPService('http://example.com/tap')
        job = service.get_job('111')
        assert job.jobid == '111'
        assert job.phase == 'EXECUTING'
        assert job.ownerid == '222'

    @pytest.mark.usefixtures('async_fixture')
    @pytest.mark.remote_data
    def test_get_job_list(self):
        service = TAPService('http://example.com/tap')
        # server returns:
        #       - 3 jobs for last atribute
        #       - 2 jobs for phase attribute
        #       - 1 job for after attribute
        # Tests consists in counting the cumulative number of jobs as per
        # above rules
        after = datetime.datetime.now(tz=datetime.timezone.utc)
        assert len(service.get_job_list()) == 0
        assert len(service.get_job_list(last=3)) == 3
        assert len(service.get_job_list(after='2018-04-25T17:46:01Z')) == 1
        assert len(service.get_job_list(phases=['EXECUTING'])) == 2
        assert len(service.get_job_list(after=after,
                                        phases=['EXECUTING'])) == 3
        assert len(service.get_job_list(after='2018-04-25T17:46:01.123Z',
                                        last=3)) == 4
        assert len(service.get_job_list(phases=['EXECUTING'], last=3)) == 5
        assert len(service.get_job_list(phases=['EXECUTING'], last=3,
                                        after=datetime.datetime.now(tz=datetime.timezone.utc))) == 6

    @pytest.mark.usefixtures('create_fixture')
    def test_create_table(self):
        prototype.activate_features('cadc-tb-upload')
        try:
            buffer = BytesIO(b'table definition in VOSITable format')
            service = TAPService('https://example.com/tap')
            service.create_table(name='abc', definition=buffer)
            tmpfile = tempfile.NamedTemporaryFile('w+b', delete=False)
            tmpfile.write(b'table definition in VOTable format here')
            tmpfile.close()
            with open(tmpfile.name, 'rb') as f:
                service.create_table('abc', definition=f, format='VOTable')

            with pytest.raises(ValueError):
                service.create_table('abc', definition=buffer, format='Unknown')
            with pytest.raises(ValueError):
                service.create_table('abc', definition=None, format='VOSITable')
        finally:
            prototype.deactivate_features('cadc-tb-upload')

    @pytest.mark.usefixtures('delete_fixture')
    def test_remove_table(self):
        prototype.activate_features('cadc-tb-upload')
        try:
            service = TAPService('https://example.com/tap')
            service.remove_table(name='abc')
        finally:
            prototype.deactivate_features('cadc-tb-upload')

    @pytest.mark.usefixtures('load_fixture')
    def test_load_table(self):
        # csv content in buffer
        prototype.activate_features('cadc-tb-upload')
        try:
            service = TAPService('https://example.com/tap')
            table_content = BytesIO(b'article,count\nart1,1\nart2,2\nart3,3')
            service.load_table(name='abc', source=table_content, format='csv')

            # tsv content in file
            tmpfile = tempfile.NamedTemporaryFile('w+b', delete=False)
            tmpfile.write(b'article\tcount\nart1\t1\nart2\t2\nart3\t3')
            tmpfile.close()
            with open(tmpfile.name, 'rb') as f:
                service.load_table('abc', source=f, format='tsv')

            # FITSTable content in file
            tmpfile = tempfile.NamedTemporaryFile('w+b', delete=False)
            tmpfile.write(b'FITSTable content here')
            tmpfile.close()
            with open(tmpfile.name, 'rb') as f:
                service.load_table('abc', source=f, format='FITSTable')

            with pytest.raises(ValueError):
                service.load_table('abc', source=table_content, format='Unknown')

            with pytest.raises(ValueError):
                service.load_table('abc', source=None, format='tsv')
        finally:
            prototype.deactivate_features('cadc-tb-upload')

    def test_create_index(self):
        prototype.activate_features('cadc-tb-upload')
        try:
            service = TAPService('https://example.com/tap')

            def match_request_text(request):
                # check details of index are present
                return 'table=abc&index=col1&unique=true' in request.text

            with requests_mock.Mocker() as rm:
                # mock initial post to table-update and the subsequent calls to
                # get, run and check status of the job
                rm.post('https://example.com/tap/table-update',
                        additional_matcher=match_request_text,
                        status_code=303,
                        headers={'Location': 'https://example.com/tap/uws'})
                rm.get('https://example.com/tap/uws',
                       [{'content': get_index_job("PENDING")},
                        {'content': get_index_job("COMPLETED")}])
                rm.post('https://example.com/tap/uws/phase', status_code=200)
                # finally the call
                service.create_index(table_name='abc', column_name='col1',
                                     unique=True)
            # test wrong return status code
                with requests_mock.Mocker() as rm:
                    # mock initial post to table-update and the subsequent calls to
                    # get, run and check status of the job
                    rm.post('https://example.com/tap/table-update',
                            additional_matcher=match_request_text,
                            status_code=200,  # NOT EXPECTED!
                            headers={'Location': 'https://example.com/tap/uws'})
                    with pytest.raises(RuntimeError):
                        service.create_index(table_name='abc', column_name='col1',
                                             unique=True)
        finally:
            prototype.deactivate_features('cadc-tb-upload')


@pytest.mark.usefixtures("tapservice")
class TestTAPCapabilities:
    def test_no_tap_cap(self):
        svc = TAPService('http://example.com/tap')
        svc.capabilities = []
        with pytest.raises(DALServiceError) as excinfo:
            svc.get_tap_capability()
        assert str(excinfo.value) == ("Invalid TAP service:"
            " Does not expose a tr:TableAccess capability")

    def test_no_adql(self):
        svc = TAPService('http://example.com/tap')
        svc.get_tap_capability()._languages = []
        with pytest.raises(VOSIError) as excinfo:
            svc.get_tap_capability().get_adql()
        assert str(excinfo.value) == ("Invalid TAP service:"
            " Does not declare an ADQL language")

    def test_get_adql(self, tapservice):
        assert tapservice.get_tap_capability().get_adql().name == "ADQL"

    def test_missing_featurelist(self, tapservice):
        assert (
            tapservice.get_tap_capability().get_adql().get_feature_list("fump")
            == [])

    def test_get_featurelist(self, tapservice):
        features = tapservice.get_tap_capability().get_adql().get_feature_list(
            "ivo://ivoa.net/std/TAPRegExt#features-adqlgeo")
        assert set(f.form for f in features) == {
            'CENTROID', 'CONTAINS', 'COORD1', 'POLYGON',
            'INTERSECTS', 'COORD2', 'BOX', 'AREA', 'DISTANCE',
            'REGION', 'CIRCLE', 'POINT'}

    def test_get_missing_feature(self, tapservice):
        assert tapservice.get_tap_capability().get_adql().get_feature(
            "ivo://ivoa.net/std/TAPRegExt#features-adqlgeo",
            "Garage") is None

    def test_get_feature(self, tapservice):
        feature = tapservice.get_tap_capability().get_adql().get_feature(
            "ivo://ivoa.net/std/TAPRegExt#features-adqlgeo",
            "AREA")
        assert feature.form == "AREA"
        assert feature.description is None

    def test_missing_udf(self, tapservice):
        assert tapservice.get_tap_capability().get_adql().get_udf("duff function") is None

    def test_get_udf(self, tapservice):
        func = tapservice.get_tap_capability().get_adql().get_udf("IVO_hasword")  # case insensitive!
        assert func.form == "ivo_hasword(haystack TEXT, needle TEXT) -> INTEGER"


def test_get_endpoint_candidates():
    # Directly instantiate the TAPService with a known base URL
    svc = TAPService("http://astroweb.projects.phys.ucl.ac.uk:8000/tap")

    # Check if the correct endpoint candidates are generated
    expected_urls = [
        "http://astroweb.projects.phys.ucl.ac.uk:8000/tap/capabilities",
        "http://astroweb.projects.phys.ucl.ac.uk:8000/capabilities"
    ]
    assert svc._get_endpoint_candidates("capabilities") == expected_urls
