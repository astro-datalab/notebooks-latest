# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module for accessing remote source and observation catalogs
"""
from functools import partial
from datetime import datetime
from time import sleep

import requests
from urllib.parse import urlparse, urljoin

from astropy.io.votable import parse as votableparse

from .query import (
    DALResults, DALQuery, DALService, Record, UploadList,
    DALServiceError, DALQueryError)
from .vosi import AvailabilityMixin, CapabilityMixin, VOSITables
from .adhoc import DatalinkResultsMixin, DatalinkRecordMixin, SodaRecordMixin

from ..io import vosi, uws
from ..io.vosi import tapregext as tr

from ..utils.formatting import para_format_desc
from ..utils.http import use_session
from ..utils.prototype import prototype_feature
import xml.etree.ElementTree
import io

__all__ = [
    "search", "escape", "TAPService", "TAPQuery", "AsyncTAPJob", "TAPResults"]

IVOA_DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%S.%fZ"

# file formats supported by table upload and their corresponding MIME types
TABLE_UPLOAD_FORMAT = {'tsv': 'text/tab-separated-values',
                       'csv': 'text/csv',
                       'FITSTable': 'application/fits'}
# file formats supported by table create and their corresponding MIME types
TABLE_DEF_FORMAT = {'VOSITable': 'text/xml',
                    'VOTable': 'application/x-votable+xml'}


def _from_ivoa_format(datetime_str):
    """
    parses an ivoa date in ISO 8601 format: YYYY-MM-DDTHH:MM:SS.[mmm]Z

    :param datetime_str:
    :return: corresponding datetime object
    """
    # TODO Replace with datetime.fromisoformat(date_string) in Python3.7+
    try:
        # with fraction of seconds first
        return datetime.strptime(datetime_str, IVOA_DATETIME_FORMAT)
    except ValueError:
        # and without
        return datetime.strptime(datetime_str, "%Y-%m-%dT%H:%M:%SZ")


def escape(term):
    """
    escapes a term for use in ADQL
    """
    return str(term).replace("'", "''")


def search(url, query, *, language="ADQL", maxrec=None, uploads=None, **keywords):
    """
    submit a Table Access query that returns rows matching the criteria given.

    Parameters
    ----------
    url : str
        the base URL of the query service.
    query : str, dict
        The query string / parameters
    language : str
        specifies the query language, default ADQL.
        useful for services which allow to use the backend query language.
    maxrec : int
        the maximum records to return. defaults to the service default
    uploads : dict
        a mapping from table names to file like objects containing a votable

    Returns
    -------
    TAPResults
        a container holding a table of matching catalog records

    Raises
    ------
    DALServiceError
        for errors connecting to or
        communicating with the service.
    DALQueryError
        if the service responds with
        an error, including a query syntax error.
    """
    service = TAPService(url)
    return service.search(query, language=language, maxrec=maxrec, uploads=uploads, **keywords)


class TAPService(DALService, AvailabilityMixin, CapabilityMixin):
    """
    a representation of a Table Access Protocol service
    """

    _tables = None
    _examples = None

    def __init__(self, baseurl, *, capability_description=None, session=None):
        """
        instantiate a Table Access Protocol service

        Parameters
        ----------
        baseurl :  str
           the base URL that should be used for forming queries to the service.
        session : object
           optional session to use for network requests
        """
        super().__init__(baseurl, session=session, capability_description=capability_description)

        # Check if the session has an update_from_capabilities attribute.
        # This means that the session is aware of IVOA capabilities,
        # and can use this information in processing network requests.
        # One such use case for this is auth.
        if hasattr(self._session, 'update_from_capabilities'):
            self._session.update_from_capabilities(self.capabilities)

    def get_tap_capability(self):
        """
        returns the (first) TAP capability of this service.

        Returns
        -------
        A `~pyvo.io.vosi.tapregext.TableAccess` instance.
        """
        for capa in self.capabilities:
            if isinstance(capa, tr.TableAccess):
                return capa
        raise DALServiceError("Invalid TAP service: Does not"
            " expose a tr:TableAccess capability")

    @property
    def tables(self):
        """
        returns tables as a dict-like object
        """
        if self._tables is None:
            tables_url = '{}/tables'.format(self.baseurl)

            response = self._session.get(tables_url, stream=True)

            try:
                response.raise_for_status()
            except requests.RequestException as ex:
                raise DALServiceError.from_except(ex, tables_url)

            # requests doesn't decode the content by default
            response.raw.read = partial(response.raw.read, decode_content=True)

            self._tables = VOSITables(
                vosi.parse_tables(response.raw.read), tables_url)
        return self._tables

    def _parse_examples(self, examples_uri, *, depth=0):
        """returns the TAP queries from a DALI examples URI.
        """
        if depth > 5:
            raise DALServiceError("Suspecting endless recursion when"
                " parsing TAP examples")

        response = self._session.get(examples_uri, stream=True)
        if response.status_code == 404:
            return []

        try:
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, examples_uri)

        try:
            root = xml.etree.ElementTree.parse(
                io.BytesIO(response.content)).getroot()
            exampleElements = root.findall('.//*[@property="query"]')
        except Exception as ex:
            raise DALServiceError.from_except(ex, examples_uri)

        examples = [TAPQuery(self.baseurl, example.text)
            for example in exampleElements]

        for continuation in root.findall('.//*[@property="continuation"]'):
            examples.extend(
                self._parse_examples(continuation.get("href"), depth=depth + 1))

        return examples

    @property
    def examples(self):
        """
        returns examples as a list of TAPQuery objects
        """
        if self._examples is None:
            examples_url = '{}/examples'.format(self.baseurl)

            self._examples = self._parse_examples(examples_url)
        return self._examples

    @property
    def maxrec(self):
        """
        the default output limit.

        Raises
        ------
        DALServiceError
            if the property is not exposed by the service
        """
        try:
            return self.get_tap_capability().outputlimit.default.content
        except AttributeError:
            pass
        raise DALServiceError("Default limit not exposed by the service")

    @property
    def hardlimit(self):
        """
        the hard output limit.

        Raises
        ------
        DALServiceError
            if the property is not exposed by the service
        """
        try:
            return self.get_tap_capability().outputlimit.hard.content
        except AttributeError:
            pass
        raise DALServiceError("Hard limit not exposed by the service")

    @property
    def upload_methods(self):
        """
        a list of upload methods in form of
        :py:class:`~pyvo.io.vosi.tapregext.UploadMethod` objects
        """
        return self.get_tap_capability().uploadmethods

    def run_sync(
            self, query, *, language="ADQL", maxrec=None, uploads=None,
            **keywords):
        """
        runs sync query and returns its result

        Parameters
        ----------
        query : str
            The query
        language : str
            specifies the query language, default ADQL.
            useful for services which allow to use the backend query language.
        maxrec : int
            the maximum records to return. defaults to the service default
        uploads : dict
            a mapping from table names to objects containing a votable

        Returns
        -------
        TAPResults
            the query result

        See Also
        --------
        TAPResults
        """
        return self.create_query(
            query, language=language, maxrec=maxrec, uploads=uploads,
            **keywords).execute()

    # alias for service discovery
    search = run_sync

    def run_async(
            self, query, *, language="ADQL", maxrec=None, uploads=None,
            **keywords):
        """
        runs async query and returns its result

        Parameters
        ----------
        query : str, dict
            the query string / parameters
        language : str
            specifies the query language, default ADQL.
            useful for services which allow to use the backend query language.
        maxrec : int
            the maximum records to return. defaults to the service default
        uploads : dict
            a mapping from table names to objects containing a votable

        Returns
        -------
        TAPResult
            the query instance

        Raises
        ------
        DALServiceError
           for errors connecting to or communicating with the service
        DALQueryError
           for errors either in the input query syntax or
           other user errors detected by the service
        DALFormatError
           for errors parsing the VOTable response

        See Also
        --------
        AsyncTAPJob
        """
        job = AsyncTAPJob.create(
            self.baseurl, query, language=language, maxrec=maxrec, uploads=uploads,
            session=self._session, **keywords)
        job = job.run().wait()
        job.raise_if_error()
        result = job.fetch_result()
        job.delete()

        return result

    def submit_job(
            self, query, *, language="ADQL", maxrec=None, uploads=None,
            **keywords):
        """
        submit a async query without starting it and returns a AsyncTAPJob
        object

        Parameters
        ----------
        query : str
            the query string / parameters
        language : str
            specifies the query language, default ADQL.
            useful for services which allow to use the backend query language.
        maxrec : int
            the maximum records to return. defaults to the service default
        uploads : dict
            a mapping from table names to objects containing a votable

        Returns
        -------
        AsyncTAPJob
            the query instance

        See Also
        --------
        AsyncTAPJob
        """
        return AsyncTAPJob.create(
            self.baseurl, query, language=language, maxrec=maxrec, uploads=uploads,
            session=self._session, **keywords)

    def create_query(
            self, query=None, *, mode="sync", language="ADQL", maxrec=None,
            uploads=None, **keywords):
        """
        create a query object that constraints can be added to and then
        executed.  The input arguments will initialize the query with the
        given values.

        Parameters
        ----------
        baseurl : str
            the base URL for the TAP service
        query : str
            the query string / parameters
        mode : str
            the query mode (sync | async). default "sync"
        language : str
            specifies the query language, default ADQL.
            useful for services which allow to use the backend query language.
        maxrec : int
            specifies the maximum records to return.
            defaults to the service default.
        uploads : dict
            a mapping from table names to objects containing a votable.
        """
        return TAPQuery(
            self.baseurl, query, mode=mode, language=language, maxrec=maxrec,
            uploads=uploads, session=self._session, **keywords)

    def get_job(self, job_id):
        """
        Returns the job corresponding to an ID. Note that the caller must be
        able to see the job in the current security context.

        Parameters
        ----------
        job_id : str
            ID of the job to view

        Returns
        -------
        `~pyvo.io.uws.tree.JobSummary` corresponding to the job ID
        """
        response = self._session.get(
            self.baseurl + '/async/' + job_id,
            stream=True)
        response.raw.read = partial(response.raw.read,
                                    decode_content=True)
        return uws.parse_job(response.raw.read)

    def get_job_list(self, *, phases=None, after=None, last=None,
                     short_description=True):
        """
        lists jobs that the caller can see in the current security context.
        The list can be filtered on the server side by the phases of the jobs,
        creation date time or
        Note that by default jobs in 'ARCHIVED' phase are not returned.

        Parameters
        ----------
        phases: list of str
            Union of job phases to filter the results by.
        after: datetime
            Return only jobs created after this datetime
        last: int
            Return only the most recent number of jobs
        short_description: flag - True or False
            If True, the jobs in the list will contain only the information
            corresponding to the TAP ShortJobDescription object (job ID, phase,
            run ID, owner ID and creation ID) whereas if False, a separate GET
            call to each job is performed for the complete job description.

        Returns
        -------
        list of `~pyvo.io.uws.tree.JobSummary`
        """

        params = {'PHASE': phases, 'LAST': last}

        if after:
            if isinstance(after, str):
                after = _from_ivoa_format(after)
            params['AFTER'] = after.strftime(IVOA_DATETIME_FORMAT)

        response = self._session.get('{}/async'.format(self.baseurl),
                                     params=params,
                                     stream=True)
        response.raw.read = partial(response.raw.read, decode_content=True)

        jobs = uws.parse_job_list(response.raw.read)
        if not short_description:
            dj = []
            for job in jobs:
                dj.append(self.get_job(job.jobid))
            return dj
        else:
            return list(jobs)

    def describe(self, width=None):
        """
        Print a summary description of this service.

        This includes the interface capabilities, and the content description
        if it doesn't contains multiple data collections (in other words, it is
        not a TAP service).
        """
        if len(self.tables) == 1:
            description = next(self.tables.values()).description

            if width:
                description = para_format_desc(description, width)

            print(description)
            print()

        capabilities = filter(
            lambda x: not str(x.standardid).startswith(
                'ivo://ivoa.net/std/VOSI'),
            self.capabilities
        )

        for cap in capabilities:
            cap.describe()
            print()

    @prototype_feature('cadc-tb-upload')
    def create_table(self, name, definition, *, format='VOSITable'):
        """
        Creates a table in the catalog service.

        Parameters
        ----------
        name: str
            Name of the table in the TAP service
        definition: stream (object with a read method)
            Definition of the table
        format: str
            Format of the table definition (VOSITable or VOTable).
        """
        if not name or not definition:
            raise ValueError(
                'table name and definition required in create: {}/{}'.
                format(name, definition))
        if format not in TABLE_DEF_FORMAT.keys():
            raise ValueError(
                'Table definition file format {} not supported ({})'.
                format(format, ' '.join(TABLE_DEF_FORMAT.keys())))

        headers = {'Content-Type': TABLE_DEF_FORMAT[format]}
        response = self._session.put('{}/tables/{}'.format(self.baseurl, name),
                                     headers=headers,
                                     data=definition)
        response.raise_for_status()

    @prototype_feature('cadc-tb-upload')
    def remove_table(self, name):
        """
        Remove a table from the catalog service (Equivalent to drop command
        in DB).

        Parameters
        ----------
        name: str
            Name of the table in the TAP service
        """
        if not name:
            raise ValueError(
                'table name required in : {}'.
                format(name))

        response = self._session.delete(
            '{}/tables/{}'.format(self.baseurl, name))
        response.raise_for_status()

    @prototype_feature('cadc-tb-upload')
    def load_table(self, name, source, *, format='tsv'):
        """
        Loads content to a table

        Parameters
        ----------
        name: str
            Name of the table
        source: stream with a read method
            Stream containing the data to be loaded
        format: str
            Format of the data source: tab-separated values(tsv),
            comma-separated values (csv) or FITS table (FITSTable)
        """
        if not name or not source:
            raise ValueError(
                'table name and source required in upload: {}/{}'.
                format(name, source))
        if format not in TABLE_UPLOAD_FORMAT.keys():
            raise ValueError(
                'Table content file format {} not supported ({})'.
                format(format, ' '.join(TABLE_UPLOAD_FORMAT.keys())))

        headers = {'Content-Type': TABLE_UPLOAD_FORMAT[format]}
        response = self._session.post(
            '{}/load/{}'.format(self.baseurl, name),
            headers=headers,
            data=source)
        response.raise_for_status()

    @prototype_feature('cadc-tb-upload')
    def create_index(self, table_name, column_name, *, unique=False):
        """
        Creates a table index in the catalog service.

        Parameters
        ----------
        table_name: str
            Name of the table
        column_name: str
            Name of the column in the table
        unique: bool
            True for unique index, False otherwise
        """
        if not table_name or not column_name:
            raise ValueError(
                'table and column names are required in index: {}/{}'.
                format(table_name, column_name))

        result = self._session.post('{}/table-update'.format(self.baseurl),
                                    data={'table': table_name,
                                          'index': column_name,
                                          'unique': 'true' if unique
                                          else 'false'},
                                    allow_redirects=False)

        if result.status_code == 303:
            job_url = result.headers['Location']
            if not job_url:
                raise RuntimeError(
                    'table update job location missing in response')
            # run the job
            job = AsyncTAPJob(job_url, session=self._session)
            job = job.run().wait()
            job.raise_if_error()
            # TODO job.delete()
        else:
            raise RuntimeError(
                'BUG: table update expected status 303 received {}'.
                format(result.status_code))


class AsyncTAPJob:
    """
    This class represents a UWS TAP Job.
    """

    _job = {}

    @classmethod
    def create(
            cls, baseurl, query, *, language="ADQL", maxrec=None, uploads=None,
            session=None, **keywords):
        """
        creates a async tap job on the server under ``baseurl``

        Parameters
        ----------
        baseurl : str
            the TAP baseurl
        query : str
            the query string
        language : str
            specifies the query language, default ADQL.
            useful for services which allow to use the backend query language.
        maxrec : int
            the maximum records to return. defaults to the service default
        uploads : dict
            a mapping from table names to objects containing a votable
        session : object
           optional session to use for network requests
        """
        tapquery = TAPQuery(
            baseurl, query, mode="async", language=language, maxrec=maxrec,
            uploads=uploads, session=session, **keywords)
        response = tapquery.submit()
        job = cls(response.url, session=session)
        return job

    def __init__(self, url, *, session=None):
        """
        initialize the job object with the given url and fetch remote values

        Parameters
        ----------
        url : str
            the job url
        """
        self._url = url
        self._session = use_session(session)
        self._update()

    def __enter__(self):
        """
        Enters the context
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exits the context. The job is silently deleted.
        """
        try:
            self.delete()
        except Exception:
            pass

    def _update(self, wait_for_statechange=False, timeout=10.):
        """
        updates local job infos with remote values
        """
        try:
            if wait_for_statechange:
                response = self._session.get(
                    self.url, stream=True, timeout=timeout, params={
                        "WAIT": "-1"
                    }
                )
            else:
                response = self._session.get(self.url, stream=True, timeout=timeout)
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, self.url)

        # requests doesn't decode the content by default
        response.raw.read = partial(response.raw.read, decode_content=True)

        self._job = uws.parse_job(response.raw.read)

    @property
    def job(self):
        """
        all up-to-date uws job infos as dictionary
        """
        # keep it up to date
        self._update()
        return self._job

    @property
    def url(self):
        """
        the job url
        """
        return self._url

    @property
    def job_id(self):
        """
        the job id
        """
        return self._job.jobid

    @property
    def phase(self):
        """
        the current query phase
        """
        self._update()
        return self._job.phase

    @property
    def execution_duration(self):
        """
        maximum execution duration as `~astropy.time.TimeDelta`.
        """
        self._update()
        return self._job.executionduration

    @execution_duration.setter
    def execution_duration(self, value):
        try:
            response = self._session.post(
                "{}/executionduration".format(self.url),
                data={"EXECUTIONDURATION": str(value)})
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, self.url)

        self._update()

    @property
    def destruction(self):
        """
        datetime after which the job results are deleted automatically.
        read-write
        """
        self._update()
        return self._job.destruction

    @destruction.setter
    def destruction(self, value):
        """
        datetime after which the job results are deleted automatically.
        read-write

        Parameters
        ----------
        value : datetime
            datetime after which the job results are deleted automatically
        """
        if isinstance(value, str):
            value = _from_ivoa_format(value)

        try:
            response = self._session.post(
                "{}/destruction".format(self.url),
                data={"DESTRUCTION": value.strftime(IVOA_DATETIME_FORMAT)})
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, self.url)

        self._update()

    @property
    def quote(self):
        """
        estimated runtime
        """
        self._update()
        return self._job.quote

    @property
    def owner(self):
        """
        job owner (if applicable)
        """
        self._update()
        return self._job.ownerid

    @property
    def query(self):
        """
        the job query
        """
        self._update()
        for parameter in self._job.parameters:
            if parameter.id_.lower() == 'query':
                return parameter.content
        return ''

    @query.setter
    def query(self, query):
        try:
            response = self._session.post(
                '{}/parameters'.format(self.url),
                data={"QUERY": query})
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, self.url)

        self._update()

    def upload(self, **kwargs):
        """
        upload a table to the job. the job must not been started.
        """
        uploads = UploadList.fromdict(kwargs)
        files = {
            upload.name: upload.fileobj()
            for upload in uploads
            if upload.is_inline
        }

        try:
            response = self._session.post(
                '{}/parameters'.format(self.url),
                data={'UPLOAD': uploads.param()},
                files=files
            )
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, self.url)

        self._update()

    @property
    def results(self):
        """
        The job results if exists
        """
        return self._job.results

    @property
    def result(self):
        """
        The job result if exists
        """
        try:
            for r in self._job.results:
                if r.id_ == 'result':
                    return r

            return self._job.results[0]
        except IndexError:
            return None

    @property
    def result_uris(self):
        """
        a list of the last result uri's
        """
        return [result.href for result in self._job.results]

    @property
    def result_uri(self):
        """
        the uri of the result
        """
        try:
            uri = self.result.href
            if not urlparse(uri).netloc:
                uri = urljoin(self.url, uri)
            return uri
        except IndexError:
            return None

    @property
    def uws_version(self):
        """
        the version of the UWS serving this async job

        Asynchronous TAP jobs are managed using a standard called Universal
        Worker Service (UWS). For instance, starting version 1.1, you can
        have long polls, which save on monitoring requests.  Normal users
        generally will not have to look at this.

        """
        self._update()
        return self._job.version

    def run(self):
        """
        starts the job / change phase to RUN
        """
        try:
            response = self._session.post(
                '{}/phase'.format(self.url), data={"PHASE": "RUN"})
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, self.url)

        return self

    def abort(self):
        """
        aborts the job / change phase to ABORT
        """
        try:
            response = self._session.post(
                '{}/phase'.format(self.url), data={"PHASE": "ABORT"})
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, self.url)

        return self

    def wait(self, *, phases=None, timeout=600.):
        """
        waits for the job to reach the given phases.

        Parameters
        ----------
        phases : list
            phases to wait for

        Raises
        ------
        DALServiceError
            if the job is in a state that won't lead to an result
        """
        if not phases:
            phases = {"COMPLETED", "ABORTED", "ERROR"}

        interval = 1.0
        increment = 1.2

        active_phases = {
            "QUEUED", "EXECUTING", "RUN", "COMPLETED", "ERROR", "UNKNOWN"}

        while True:
            self._update(wait_for_statechange=True, timeout=timeout)
            # use the cached value
            cur_phase = self._job.phase

            if cur_phase not in active_phases:
                raise DALServiceError(
                    "Cannot wait for job completion. Job is not active!")

            if cur_phase in phases:
                break

            # fallback for uws 1.0 or unsupported WAIT parameter
            sleep(interval)
            interval = min(120, interval * increment)

        return self

    def delete(self):
        """
        deletes the job. this object will become invalid.
        """
        try:
            response = self._session.delete(
                self.url,
                allow_redirects=False)
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, self.url)

        self._url = None

    def raise_if_error(self):
        """
        raise a exception if theres an error

        Raises
        ------
        DALQueryError
            if theres an error
        """
        if self.phase in {"ERROR", "ABORTED"}:
            msg = ""
            if self._job and self._job.errorsummary:
                msg = self._job.errorsummary.message.content
            msg = msg or "<No useful error from server>"
            raise DALQueryError("Query Error: " + msg, self.url)

    def fetch_result(self):
        """
        returns the result votable if query is finished
        """
        try:
            response = self._session.get(self.result_uri, stream=True)
            response.raise_for_status()
        except requests.RequestException as ex:
            self._update()
            # we propably got a 404 because query error. raise with error msg
            self.raise_if_error()
            raise DALServiceError.from_except(ex, self.url)

        response.raw.read = partial(
            response.raw.read, decode_content=True)
        return TAPResults(votableparse(response.raw.read), url=self.result_uri, session=self._session)


class TAPQuery(DALQuery):
    """
    a class for preparing a query to a TAP service.  Query constraints
    are added via service type-specific methods.  The various execute()
    functions will submit the query and return the results.

    The base URL for the query, which controls where the query will be sent
    when one of the execute functions is called, is typically set at
    construction time; however, it can be updated later via the
    :py:attr:`~pyvo.dal.DALQuery.baseurl` to send a configured
    query to another service.

    In addition to the search constraint attributes described below, search
    parameters can be set generically by name via dict semantics.

    The typical function for submitting the query is ``execute()``; however,
    alternate execute functions provide the response in different forms,
    allowing the caller to take greater control of the result processing.
    """

    def __init__(
            self, baseurl, query, *, mode="sync", language="ADQL", maxrec=None,
            uploads=None, session=None, **keywords):
        """
        initialize the query object with the given parameters

        Parameters
        ----------
        baseurl : str
            the TAP baseurl
        query : str
            the query string
        mode : str
            the query mode (sync | async). default "sync"
        language : str
            the query language. defaults to ADQL
        maxrec : int
            the amount of records to fetch
        uploads : dict
            Files to upload. Uses table name as key and table content as value.
        session : object
           optional session to use for network requests
        """
        baseurl = baseurl.rstrip("?")

        super().__init__(baseurl, session=session, **keywords)

        self._mode = mode if mode in ("sync", "async") else "sync"
        self._uploads = UploadList.fromdict(uploads or {})

        self["REQUEST"] = "doQuery"
        self["LANG"] = language

        if maxrec:
            self["MAXREC"] = maxrec

        self["QUERY"] = query

        if self._uploads:
            self["UPLOAD"] = self._uploads.param()

    @property
    def queryurl(self):
        """
        the URL to which to submit queries

        In TAP, that varies depending on whether we run sync or async
        queries.

        """
        return '{baseurl}/{mode}'.format(baseurl=self.baseurl, mode=self._mode)

    def execute_stream(self, *, post=False):
        """
        submit the query and return the raw VOTable XML as a file stream

        Raises
        ------
        DALServiceError
           for errors connecting to or communicating with the service
        DALQueryError
           for errors in the input query syntax
        """

        # theres nothing to execute in non-sync queries
        if self._mode != "sync":
            raise DALServiceError(
                "Cannot execute a non-synchronous query. Use submit instead")

        return super().execute_stream(post=post)

    def execute(self):
        """
        submit the query and return the results as a TAPResults instance

        Raises
        ------
        DALServiceError
           for errors connecting to or communicating with the service
        DALQueryError
           for errors either in the input query syntax or
           other user errors detected by the service
        DALFormatError
           for errors parsing the VOTable response
        """
        return TAPResults(self.execute_votable(), url=self.queryurl, session=self._session)

    def submit(self, *, post=False):
        """
        Does the request part of the TAP query.
        This function is separated from response parsing because async queries
        return no votable but behave like sync queries in terms of request.
        It returns the requests response.
        """
        url = self.queryurl

        files = {
            upload.name: upload.fileobj()
            for upload in self._uploads
            if upload.is_inline
        }

        response = self._session.post(
            url, data=self, stream=True, files=files)
        # requests doesn't decode the content by default
        response.raw.read = partial(response.raw.read, decode_content=True)
        return response


class TAPResults(DatalinkResultsMixin, DALResults):
    """
    The list of matching images resulting from an image (SIA) query.
    Each record contains a set of metadata that describes an available
    image matching the query constraints.  The number of records in
    the results is available by passing it to the Python built-in ``len()`` function.

    This class supports iterable semantics; thus,
    individual records (in the form of
    :py:class:`~pyvo.dal.Record` instances) are typically
    accessed by iterating over an ``TAPResults`` instance.

    Alternatively, records can be accessed randomly via
    :py:meth:`getrecord` or through a Python Database API (v2)
    Cursor (via :py:meth:`~pyvo.dal.DALResults.cursor`).
    Column-based data access is possible via the
    :py:meth:`~pyvo.dal.query.DALResults.getcolumn` method.

    ``TAPResults`` is essentially a wrapper around an Astropy
    :py:mod:`~astropy.io.votable`
    :py:class:`~astropy.io.votable.tree.TableElement` instance where the
    columns contain the various metadata describing the images.
    One can access that VOTable directly via the
    :py:attr:`pyvo.dal.DALResults.votable` attribute.  Thus,
    when one retrieves a whole column via
    :py:meth:`~pyvo.dal.DALResults.getcolumn`, the result is
    a Numpy array.  Alternatively, one can manipulate the results
    as an Astropy :py:class:`astropy.table.Table` via the
    following conversion:

    ``table = results.to_table``

    ``SIAResults`` supports the array item operator ``[...]`` in a
    read-only context.  When the argument is numerical, the result
    is an
    :py:class:`~pyvo.dal.Record` instance, representing the
    record at the position given by the numerical index.  If the
    argument is a string, it is interpreted as the name of a column,
    and the data from the column matching that name is returned as
    a Numpy array.
    """

    @property
    def infos(self):
        """
        return the info element as dictionary
        """
        return getattr(self, "_infos", {})

    @property
    def query_status(self):
        """
        return the query status
        """
        return getattr(self, "_infos", {}).get("QUERY_STATUS", None)

    def getrecord(self, index):
        """
        return a representation of a tap result record that follows
        dictionary semantics. The keys of the dictionary are those returned by
        this instance's fieldnames attribute. The returned record has
        additional image-specific properties

        Parameters
        ----------
        index : int
           the integer index of the desired record where 0 returns the first
           record

        Returns
        -------
        REc
           a dictionary-like wrapper containing the result record metadata.

        Raises
        ------
        IndexError
           if index is negative or equal or larger than the number of rows in
           the result table.

        See Also
        --------
        Record
        """
        return TAPRecord(self, index, session=self._session)


class TAPRecord(SodaRecordMixin, DatalinkRecordMixin, Record):
    pass
