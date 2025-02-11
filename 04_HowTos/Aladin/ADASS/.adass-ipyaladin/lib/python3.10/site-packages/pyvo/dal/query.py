# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module for walking through the query response of VO data access layer
(DAL) queries and general VOTable-based datasets.

Most data queries in the VO return a table as a result, usually
formatted as a VOTable.  Each row of the table describes a single
physical or virtual dataset which can be retrieved.  For uniformity,
datasets are described via standard metadata defined by a data model
specific to the type of data being queried.  The fields of the data
model are identified most generally by their VOClient alias as defined
in this interface, or at a lower level by the Utype or UCD of the
specific standard and version of the standard being queried.  While
the data model differs depending upon the type of data being queried,
the form of the query response is the same for all classes of data,
allowing a common query response interface to be used.

An exception to this occurs when querying an astronomical catalog or
other externally defined table.  In this case there is no VO defined
standard data model.  Usually the field names are used to uniquely
identify table columns.
"""
__all__ = ["DALService", "DALQuery", "DALResults", "Record"]

import os
import shutil
import re
import requests
from collections.abc import Mapping

import collections

from warnings import warn

from astropy.table import Table, QTable
from astropy.io.votable import parse as votableparse
from astropy.io.votable.ucd import parse_ucd
from astropy.utils.exceptions import AstropyDeprecationWarning

from .mimetype import mime_object_maker
from .exceptions import (DALFormatError, DALServiceError, DALQueryError,
                         DALOverflowWarning)

from .. import samp

from ..utils.decorators import stream_decode_content
from ..utils.http import use_session


class DALService:
    """
    an abstract base class representing a DAL service located a particular
    endpoint.
    """

    def __init__(self, baseurl, *, session=None, capability_description=None,):
        """
        instantiate the service connecting it to a base URL

        Parameters
        ----------
        baseurl :  str
           the base URL that should be used for forming queries to the service.
        session : object
           optional session to use for network requests
        description : str, optional
           the description of the service.
        """
        self._baseurl = baseurl
        self._capability_description = capability_description
        self._session = use_session(session)

    @property
    def baseurl(self):
        """
        the base URL identifying the location of the service and where
        queries are submitted (read-only)
        """
        return self._baseurl

    @property
    def capability_description(self):
        """
        The service description.
        """
        return self._capability_description

    def __repr__(self) -> str:
        return (f"{type(self).__name__}(baseurl : '{self.baseurl}',"
                f" description : '{self.capability_description}')")

    def search(self, **keywords):
        """
        send a search query to this service.

        This implementation has no knowledge of the type of service being
        queried.  The query parameters are given as arbitrary keywords which
        will be assumed to be understood by the service (i.e. there is no
        argument checking).  The response is a generic DALResults object.

        Raises
        ------
        DALServiceError
           for errors connecting to or communicating with the service
        DALQueryError
           for errors either in the input query syntax or other user errors
           detected by the service
        DALFormatError
           for errors parsing the VOTable response
        """
        q = self.create_query(**keywords)
        return q.execute()

    def create_query(self, **keywords):
        """
        create a query object that constraints can be added to and then
        executed.

        Returns
        -------
        DALQuery
           a generic query object
        """
        q = DALQuery(self.baseurl, session=self._session, **keywords)
        return q

    def describe(self):
        """
        describe the general information about the DAL service
        """
        print('DAL Service at {}'.format(self.baseurl))


class DALQuery(dict):
    """
    a class for preparing a query to a particular service.  Query constraints
    are added via its service type-specific methods. The various execute()
    functions will submit the query and return the results.

    The base URL for the query can be changed via the baseurl property.

    A session can also optionally be passed in that will be used for
    network transactions made by this object to remote services.
    """

    _ex = None

    def __init__(self, baseurl, *, session=None, **keywords):
        """
        initialize the query object with a baseurl
        """
        if isinstance(baseurl, bytes):
            baseurl = baseurl.decode("utf-8")

        self._baseurl = baseurl.rstrip("?")
        self._session = use_session(session)

        self.update({key.upper(): value for key, value in keywords.items()})

    @property
    def baseurl(self):
        """
        the base URL that this query will be sent to when one of the
        execute functions is called.
        """
        return self._baseurl

    def execute(self):
        """
        submit the query and return the results as a Results subclass instance

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
        return DALResults(self.execute_votable(), url=self.queryurl, session=self._session)

    def execute_raw(self):
        """
        submit the query and return the raw response as a string.

        No exceptions are raised here because non-2xx responses might still
        contain payload. They can be raised later by calling ``raise_if_error``
        """
        f = self.execute_stream()
        out = None
        try:
            out = f.read()
        finally:
            f.close()
        return out

    @stream_decode_content
    def execute_stream(self, *, post=False):
        """
        Submit the query and return the raw response as a file stream.

        No exceptions are raised here because non-2xx responses might still
        contain payload. They can be raised later by calling ``raise_if_error``
        """
        response = self.submit(post=post)

        try:
            response.raise_for_status()
        except requests.RequestException as ex:
            # save for later use
            self._ex = ex
        finally:
            return response.raw

    def submit(self, *, post=False):
        """
        does the actual request
        """
        url = self.queryurl
        params = {k: v for k, v in self.items()}

        if post:
            response = self._session.post(url, data=params, stream=True,
                                          allow_redirects=True)
        else:
            response = self._session.get(url, params=params, stream=True,
                                         allow_redirects=True)
        return response

    def execute_votable(self, *, post=False):
        """
        Submit the query and return the results as an AstroPy votable instance.
        As this is the level where qualified error messages are available,
        they are raised here instead of in the underlying execute_stream.

        Returns
        -------
        astropy.io.votable.tree.VOTableFile
           an Astropy votable object

        Raises
        ------
        DALServiceError
           for errors connecting to or communicating with the service
        DALFormatError
           for errors parsing the VOTable response

        See Also
        --------
        astropy.io.votable
        DALServiceError
        DALFormatError
        DALQueryError
        """
        try:
            return votableparse(self.execute_stream(post=post).read)
        except Exception as e:
            self.raise_if_error()
            raise DALFormatError(e, self.queryurl)

    def raise_if_error(self):
        """
        Raise if there was an error on http level.
        """
        if self._ex:
            e = self._ex
            raise DALServiceError.from_except(e, self.queryurl)

    @property
    def queryurl(self):
        """
        The URL that encodes the current query. This is the
        URL that the execute functions will use if called next.
        """
        return self.baseurl


class DALResults:
    """
    Results from a DAL query.  It provides random access to records in
    the response.  Alternatively, it can provide results via a Cursor
    (compliant with the Python Database API) or an iterable.
    """
    @classmethod
    @stream_decode_content
    def _from_result_url(cls, result_url, session):
        return session.get(result_url, stream=True).raw

    @classmethod
    def from_result_url(cls, result_url, *, session=None):
        """
        Create a result object from a url.

        Uses the optional session to make the request.
        """
        session = use_session(session)
        return cls(
            votableparse(cls._from_result_url(result_url, session).read),
            url=result_url,
            session=session)

    def __init__(self, votable, *, url=None, session=None):
        """
        initialize the cursor.  This constructor is not typically called
        by directly applications; rather an instance is obtained from calling
        a DALQuery's execute().

        Parameters
        ----------
        votable : str
           the service response parsed into an
           astropy.io.votable.tree.VOTableFile instance.
        url : str
           the URL that produced the response
        session : object
           optional session to use for network requests

        Raises
        ------
        DALFormatError
           if the response VOTable does not contain a response table

        See Also
        --------
        pyvo.dal.DALFormatError
        """
        self._votable = votable

        self._url = url
        self._session = use_session(session)

        self._status = self._findstatus(votable)
        if self._status[0].lower() not in ("ok", "overflow"):
            raise DALQueryError(self._status[1], self._status[0], url)

        if self._status[0].lower() == "overflow":
            warn("Partial result set. Potential causes MAXREC, async storage space, etc.",
                 category=DALOverflowWarning)

        self._resultstable = self._findresultstable(votable)
        if not self._resultstable:
            raise DALFormatError(
                reason="VOTable response missing results table", url=url)

        self._fldnames = tuple(
            field.name for field in self._resultstable.fields)

        if not self._fldnames:
            raise DALFormatError(
                reason="response table missing column descriptions.", url=url)

        self._infos = self._findinfos(votable)

    def _findresultstable(self, votable):
        # this can be overridden to specialize for a particular DAL protocol
        res = self._findresultsresource(votable)
        if not res or len(res.tables) < 1:
            return None
        return res.tables[0]

    def _findresultsresource(self, votable):
        # this can be overridden to specialize for a particular DAL protocol
        if len(votable.resources) < 1:
            return None
        for res in votable.resources:
            if res.type.lower() == "results":
                return res
        return votable.resources[0]

    def _findstatus(self, votable):
        # this can be overridden to specialize for a particular DAL protocol

        # look first in the result resource
        res = self._findresultsresource(votable)
        if res:
            # should be a RESOURCE/INFO
            info = self._findstatusinfo(res.infos)
            if info:
                return (info.value, info.content)

            # if not there, check inside first table
            if len(res.tables) > 0:
                info = self._findstatusinfo(res.tables[0].infos)
                if info:
                    return (info.value, info.content)

        # otherwise, look just below the root element
        info = self._findstatusinfo(votable.infos)
        if info:
            return (info.value, info.content)

        # assume it's okay
        return ("OK", "QUERY_STATUS not specified")

    def _findstatusinfo(self, infos):
        # this can be overridden to specialize for a particular DAL protocol
        status = None
        # return the last status to catch potential overflow or error sent
        # after the results.
        for info in infos:
            if info.name.lower() == 'query_status':
                status = info
        return status

    def _findinfos(self, votable):
        # this can be overridden to specialize for a particular DAL protocol
        infos = {}
        res = self._findresultsresource(votable)
        for info in res.infos:
            infos[info.name] = info.value
        for info in votable.infos:
            infos[info.name] = info.value
        return infos

    def __repr__(self):
        return f"<DALResults{repr(self.to_table())[1:]}"

    @property
    def queryurl(self):
        """
        the URL query that produced these results.  None is returned if unknown
        """
        return self._url

    @property
    def votable(self):
        """
        The complete votable XML Document `astropy.io.votable.tree.VOTableFile`
        """
        return self._votable

    @property
    def resultstable(self):
        """
        The votable XML element `astropy.io.votable.tree.TableElement`
        """
        return self._resultstable

    def to_table(self):
        """
        Returns a astropy Table object.

        Returns
        -------
        `astropy.table.Table`
        """
        return self.resultstable.to_table(use_names_over_ids=True)

    def to_qtable(self):
        """
        Returns a astropy QTable object containing quantities instead of simple
        values

        Returns
        -------
        `astropy.table.QTable`
        """
        return QTable(self.resultstable.to_table(use_names_over_ids=True))

    @property
    def table(self):
        warn(AstropyDeprecationWarning(
            'Using the table property is deprecated. '
            'Please use the to_table() instead.'
        ))
        return self.to_table()

    def __len__(self):
        """
        return the record count
        """
        return len(self.resultstable.array)

    def __getitem__(self, indx):
        """
        if indx is a string, r[indx] will return the field with the name of
        indx; if indx is an integer, r[indx] will return the indx-th record.
        """
        if isinstance(indx, int):
            return self.getrecord(indx)
        elif isinstance(indx, tuple):
            return self.getvalue(*indx)
        else:
            return self.getcolumn(indx)

    @property
    def fieldnames(self):
        """
        return the names of the columns.  These are the names that are used
        to access values from the dictionaries returned by getrecord().  They
        correspond to the column name.
        """
        return self._fldnames

    @property
    def fielddescs(self):
        """
        return the full metadata the columns as a list of Field instances,
        a simple object with attributes corresponding the the VOTable FIELD
        attributes, namely: name, id, type, ucd, utype, arraysize, description
        """
        return self.resultstable.fields

    @property
    def status(self):
        """
        The query status as a 2-element tuple e.g. ('OK', 'Everythings fine')
        """
        return self._status

    def fieldname_with_ucd(self, ucd):
        """
        return the field name that has a given UCD value or None if the UCD
        is not found.
        """
        search_ucds = set(parse_ucd(ucd, has_colon=True))

        for field in (field for field in self.fielddescs if field.ucd):
            field_ucds = set(parse_ucd(field.ucd, has_colon=True))

            if search_ucds & field_ucds:
                return field.name

        return None

    def fieldname_with_utype(self, utype):
        """
        return the field name that has a given UType value or None if the UType
        is not found.
        """
        utype = utype.lower()
        try:
            iterchain = (
                self.getdesc(fieldname) for fieldname in self.fieldnames)
            iterchain = (field for field in iterchain if (field.utype or "").lower() == utype)
            return next(iterchain).name
        except StopIteration:
            return None

    def getcolumn(self, name):
        """
        return a numpy array containing the values for the column with the
        given name
        """
        try:
            if name not in self.fieldnames:
                name = self.resultstable.get_field_by_id(name).name

            return self.resultstable.array[name]
        except KeyError:
            raise KeyError("No such column: {}".format(name))

    def getrecord(self, index):
        """
        return a representation of a result record that follows dictionary
        semantics.  The keys of the dictionary are those returned by this
        instance's fieldnames attribute.The returned record may have additional
        accessor methods for getting at stardard DAL response metadata
        (e.g. ra, dec).

        Parameters
        ----------
        index : int
           the integer index of the desired record where 0 returns the first
           record

        Returns
        -------
        Record
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
        return Record(self, index, session=self._session)

    def getvalue(self, name, index):
        """
        return the value of a record attribute--a value from a column and row.

        Parameters
        ----------
        name : str
           the name of the attribute (column)
        index : int
           the zero-based index of the record

        Raises
        ------
        IndexError
           if index is negative or equal or larger than the
           number of rows in the result table.
        KeyError
           if name is not a recognized column name
        """
        return self.getrecord(index)[name]

    def getdesc(self, name):
        """
        return the field description for the record attribute (column) with
        the given name

        Parameters
        ----------
        name : str
           the name of the attribute (column)

        Returns
        -------
        object
           with attributes (name, id, datatype, unit, ucd, utype, arraysize)
           which describe the column

        """
        if name not in self._fldnames:
            raise KeyError(name)
        return self.resultstable.get_field_by_id_or_name(name)

    def __iter__(self):
        """
        return a python iterable for stepping through the records in this
        result
        """
        pos = 0

        while True:
            try:
                out = self.getrecord(pos)
            except IndexError:
                break

            yield out
            pos += 1

    def broadcast_samp(self, *, client_name=None):
        """
        Broadcast the table to ``client_name`` via SAMP
        """
        with samp.connection() as conn:
            samp.send_table_to(
                conn, self.to_table(),
                client_name=client_name, name=self.queryurl)

    def cursor(self):
        """
        return a cursor that is compliant with the Python Database API's
        ``Cursor`` interface.  See PEP 249 for details.
        """
        from .dbapi2 import Cursor
        return Cursor(self)


class Record(Mapping):
    """
    one record from a DAL query result.  The column values are accessible
    as dictionary items.  It also provides special added functions for
    accessing the dataset the record corresponds to.  Subclasses may provide
    additional functions for access to service type-specific data.
    """

    def __init__(self, results, index, *, session=None):
        self._results = results
        self._index = index
        self._session = use_session(session)
        self._mapping = collections.OrderedDict(
            zip(
                results.fieldnames,
                results.resultstable.array.data[index]
            )
        )

    def __getitem__(self, key):
        try:
            if key not in self._mapping:
                key = self._results.resultstable.get_field_by_id(key).name

            return self._mapping[key]
        except KeyError:
            raise KeyError("No such column: {}".format(key))

    def __iter__(self):
        return iter(self._mapping)

    def __len__(self):
        return len(self._mapping)

    def __repr__(self):
        return repr(tuple(self.values()))

    def get(self, key, default=None, decode=False):
        """
        This method mimics the dict get method and adds a decode parameter
        to allow decoding of binary strings.
        """
        out = self._mapping.get(key, default)

        if decode and isinstance(out, bytes):
            out = out.decode('ascii')

        return out

    def getbyucd(self, ucd, default=None, decode=False):
        """
        return the column with the given ucd.
        """
        return self.get(
            self._results.fieldname_with_ucd(ucd), default, decode)

    def getbyutype(self, utype, default=None, decode=False):
        """
        return the column with the given utype.

        Raises
        ------
        KeyError
            if theres no column with the given utype.
        """
        return self.get(
            self._results.fieldname_with_utype(utype), default, decode)

    def getdataformat(self):
        """
        return the mimetype of the dataset described by this record.
        """
        return self.getbyucd('meta.code.mime', decode=True)

    def getdataurl(self):
        """
        return the URL contained in the access URL column which can be used
        to retrieve the dataset described by this record.  None is returned
        if no such column exists.
        """
        for fieldname in self._results.fieldnames:
            field = self._results.getdesc(fieldname)
            if (field.utype and "access.reference" in field.utype.lower()) or (
                    field.ucd and "meta.dataset" in field.ucd
                    and "meta.ref.url" in field.ucd
            ):
                out = self[fieldname]
                if isinstance(out, bytes):
                    out = out.decode('utf-8')
                return out
        return None

    def getdataobj(self):
        """
        return the appropriate data object suitable for the data content behind
        this record.
        """
        return mime_object_maker(self.getdataurl(), self.getdataformat())

    @stream_decode_content
    def getdataset(self, timeout=None):
        """
        Get the dataset described by this record from the server.

        Parameters
        ----------
        timeout : float
           the time in seconds to allow for a successful
           connection with server before failing with an
           IOError (specifically, socket.timeout) exception

        Returns
        -------
            A file-like object which may be read to retrieve the referenced
            dataset.

        Raises
        ------
        KeyError
           if no datast access URL is included in the record
        URLError
           if the dataset access URL is invalid (note: subclass of IOError)
        HTTPError
           if an HTTP error occurs while accessing the dataset
           (note: subclass of IOError)
        socket.timeout
           if the timeout is exceeded before a connection is established.
           (note: subclass of IOError)
        IOError
           if some other error occurs while establishing the data stream.
        """
        url = self.getdataurl()
        if not url:
            raise KeyError("no dataset access URL recognized in record")

        if timeout:
            response = self._session.get(url, stream=True, timeout=timeout)
        else:
            response = self._session.get(url, stream=True)
        try:
            response.raise_for_status()
        except requests.RequestException as ex:
            raise DALServiceError.from_except(ex, url)

        return response.raw

    def cachedataset(self, *, filename=None, dir=".", timeout=None, bufsize=None):
        """
        retrieve the dataset described by this record and write it out to
        a file with the given name.  If the file already exists, it will be
        over-written.

        Parameters
        ----------
        filename : str
           the name of the file to write dataset to.  If the
           value represents a relative path, it will be taken
           to be relative to the value of the ``dir``
           parameter.  If None, a default name is attempted
           based on the record title and format.
        dir : str
           the directory to write the file into.  This value
           will be ignored if filename is an absolute path.
        timeout : int
           the time in seconds to allow for a successful
           connection with server before failing with an
           IOError (specifically, socket.timeout) exception
        bufsize : int
           a buffer size in bytes for copying the data to disk
           (default: 0.5 MB)

        Raises
        ------
        KeyError
            if no datast access URL is included in the record
        URLError
           if the dataset access URL is invalid
        HTTPError
           if an HTTP error occurs while accessing the dataset
        socket.timeout
           if the timeout is exceeded before a connection is established.
           (note: subclass of IOError)
        IOError
            if an error occurs while writing out the dataset
        """
        if not bufsize:
            bufsize = 524288

        if not filename:
            filename = self.make_dataset_filename(dir=dir)

        inp = self.getdataset(timeout)
        try:
            with open(filename, 'wb') as out:
                shutil.copyfileobj(inp, out)
        finally:
            inp.close()

    _dsname_no = 0  # used by make_dataset_filename

    def make_dataset_filename(self, *, dir=".", base=None, ext=None):
        """
        create a viable pathname in a given directory for saving the dataset
        available via getdataset().  The pathname that is returned is
        guaranteed not to already exist (under single-threaded conditions).

        This implementation will first try combining the base name with the
        file extension (with a dot).  If this file already exists in the
        directory, a name that appends an integer suffix ("-#") to the base
        before joining with the extension will be tried.  The integer will
        be incremented until a non-existent filename is created.

        Parameters
        ----------
        dir : str
           the directory to save the dataset under.  This must already exist.
        base : str
           a basename to use to as the base of the filename.  If None, the
           result of ``suggest_dataset_basename()`` will be used.
        ext : str
           the filename extension to use.  If None, the result of
           ``suggest_extension()`` will be used.
        """
        if not dir:
            raise ValueError(
                "make_dataset_filename(): no dir parameter provided")
        if not os.path.exists(dir):
            os.mkdir(dir)
        if not os.path.isdir(dir):
            raise ValueError("{}: not a directory".format(dir))

        if not base:
            base = self.suggest_dataset_basename()
        if not ext:
            ext = self.suggest_extension(default="dat")

        base = base.replace("/", "_"
            ).replace("\\", "_")

        # be efficient when writing a bunch of files into the same directory
        # in succession
        n = self._dsname_no

        def mkpath(i):
            return os.path.join(dir, "{}-{}.{}".format(base, i, ext))

        if n > 0:
            # find the last file written of the form, base-n.ext
            while n > 0 and not os.path.exists(mkpath(n)):
                n -= 1
        if n > 0:
            n += 1
        if n == 0:
            # never wrote a file of form, base-n.ext; try base.ext
            path = os.path.join(dir, "{}.{}".format(base, ext))
            if not os.path.exists(path):
                return path
            n += 1
        # find next available name
        while os.path.exists(mkpath(n)):
            n += 1
        self._dsname_no = n
        return mkpath(n)

    def suggest_dataset_basename(self):
        """
        return a default base filename that the dataset available via
        ``getdataset()`` can be saved as.  This function is
        specialized for a particular service type this record originates from
        so that it can be used by ``cachedataset()`` via
        ``make_dataset_filename()``.
        """
        # abstract; specialized for the different service types
        return "dataset"

    def suggest_extension(self, *, default=None):
        """
        returns a recommended filename extension for the dataset described
        by this record.  Typically, this would look at the column describing
        the format and choose an extension accordingly.  This function is
        specialized for a particular service type this record originates from
        so that it can be used by ``cachedataset()`` via
        ``make_dataset_filename()``.
        """
        # abstract; specialized for the different service types
        return default


class Iter:
    def __init__(self, res):
        self.resultset = res
        self.pos = 0

    def __iter__(self):
        return self

    def __next__(self):
        try:
            out = self.resultset.getrecord(self.pos)
            self.pos += 1
            return out
        except IndexError:
            raise StopIteration()

    next = __next__


class Upload:
    """
    This class represents a DALI Upload as described in
    http://www.ivoa.net/documents/DALI/20161101/PR-DALI-1.1-20161101.html#tth_sEc3.4.5
    """

    def __init__(self, name, content):
        """
        Initialise the Upload object with the given parameters

        Parameters
        ----------
        name : str
            Tablename for use in queries
        content : object
            If its a file-like object, a string pointing to a local file,
            a `DALResults` object or a astropy table, `is_inline` will be true
            and it will expose a file-like object under `fileobj`

            Otherwise it exposes a URI under `uri`
        """
        try:
            self._is_file = os.path.isfile(content)
        except Exception:
            self._is_file = False
        self._is_fileobj = hasattr(content, "read")
        self._is_table = isinstance(content, Table)
        self._is_resultset = isinstance(content, DALResults)

        self._inline = any((
            self._is_file,
            self._is_fileobj,
            self._is_table,
            self._is_resultset,
        ))

        self._name = name
        self._content = content

    @property
    def is_inline(self):
        """
        True if the upload can be inlined
        """
        return self._inline

    @property
    def name(self):
        return self._name

    def fileobj(self):
        """
        A file-like object for a local resource

        Raises
        ------
        ValueError
            if theres no valid local resource
        """

        if not self.is_inline:
            raise ValueError(
                "Upload {name} doesn't refer to a local resource".format(
                    name=self.name))

        # astropy table
        if isinstance(self._content, Table):
            from io import BytesIO
            fileobj = BytesIO()

            self._content.write(output=fileobj, format="votable")
            fileobj.seek(0)

            return fileobj
        elif isinstance(self._content, DALResults):
            from io import BytesIO
            fileobj = BytesIO()

            table = self._content.to_table()
            table.write(output=fileobj, format="votable")
            fileobj.seek(0)

            return fileobj

        fileobj = self._content
        try:
            fileobj = open(self._content)
        finally:
            return fileobj

    def uri(self):
        """
        The URI pointing to the result
        """

        # TODO: use a async job base class instead of hasattr for inspection
        if hasattr(self._content, "result_uri"):
            self._content.raise_if_error()
            uri = self._content.result_uri
        else:
            uri = str(self._content)

        return uri

    def query_part(self):
        """
        The query part for use in DALI requests
        """

        if self.is_inline:
            value = "{name},param:{name}"
        else:
            value = "{name},{uri}"

        return value.format(name=self.name, uri=self.uri())


class UploadList(list):
    """
    This class extends the native python list with utility functions for
    upload handling
    """

    @classmethod
    def fromdict(cls, dct):
        """
        Constructs a upload list from a dictionary with table_name: content
        """
        return cls(Upload(key, value) for key, value in dct.items())

    def param(self):
        """
        Returns a string suitable for use in UPLOAD parameters
        """
        return ";".join(upload.query_part() for upload in self)


_image_mt_re = re.compile(r'^image/(\w+)')
_text_mt_re = re.compile(r'^text/(\w+)')
_votable_mt_re = re.compile(r'^(\w+)/(x-)?votable(\+\w+)?')
