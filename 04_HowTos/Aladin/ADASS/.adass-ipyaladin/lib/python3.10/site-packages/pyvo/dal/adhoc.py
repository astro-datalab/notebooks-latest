# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Datalink classes and mixins
"""
import numpy as np
import warnings
import copy
import requests
from collections import OrderedDict

from .query import DALResults, DALQuery, DALService, Record
from .exceptions import DALServiceError
from .vosi import AvailabilityMixin, CapabilityMixin
from .params import find_param_by_keyword, get_converter

from astropy.io.votable.tree import Param
from astropy import units as u
from astropy.units import Quantity, Unit
from astropy.units import spectral as spectral_equivalencies

from astropy.io.votable.tree import Resource, Group
from astropy.utils.collections import HomogeneousList

from ..utils.decorators import stream_decode_content
from ..utils import vocabularies
from .params import PosQueryParam, IntervalQueryParam, TimeQueryParam, EnumQueryParam
from ..dam.obscore import POLARIZATION_STATES

# calls to DataLink from the results pages are batched for performance
# reasons. This is the size of a batch
DATALINK_BATCH_CALL_SIZE = 50

SODA_SYNC_IVOID = 'ivo://ivoa.net/std/SODA#sync-1'
DATALINK_IVOID = 'ivo://ivoa.net/std/datalink'

# MIME types
DATALINK_MIME_TYPE = 'application/x-votable+xml;content=datalink'


# monkeypatch astropy with group support in RESOURCE
def _monkeypath_astropy_resource_groups():
    old_group_unknown_tag = Group._add_unknown_tag

    def new_group_unknown_tag(self, iterator, tag, data, config, pos):
        if tag == "PARAM":
            return self._add_param(self, iterator, tag, data, config, pos)
        else:
            old_group_unknown_tag(self, iterator, tag, data, config, pos)

    old_init = Resource.__init__

    def new_init(self, *args, **kwargs):
        old_init(self, *args, **kwargs)
        self._groups = HomogeneousList(Group)
    Resource.__init__ = new_init

    def resource_groups(self):
        return self._groups
    Resource.groups = property(resource_groups)

    def resource_add_group(self, iterator, tag, data, config, pos):
        group = Group(self, config=config, pos=pos, **data)
        self.groups.append(group)
        group.parse(iterator, config)
    Resource._add_group = resource_add_group

    old_resource_unknown_tag = Resource._add_unknown_tag

    def new_resource_unknown_tag(self, iterator, tag, data, config, pos):
        if tag == "GROUP":
            return self._add_group(iterator, tag, data, config, pos)
        else:
            old_resource_unknown_tag(self, iterator, tag, data, config, pos)
    Resource._add_unknown_tag = new_resource_unknown_tag


if not hasattr(Resource, 'groups'):
    _monkeypath_astropy_resource_groups()

__all__ = [
    "AdhocServiceResultsMixin", "DatalinkResultsMixin", "DatalinkRecordMixin",
    "DatalinkService", "DatalinkQuery", "DatalinkResults", "DatalinkRecord",
    "SodaRecordMixin", "SodaQuery"]


def _get_input_params_from_resource(resource):
    # get the group with name inputParams
    group_input_params = next(
        group for group in resource.groups if group.name == "inputParams")
    # get only Param elements from the group
    return {
        _.name: _ for _ in group_input_params.entries if isinstance(_, Param)}


def _get_params_from_resource(resource):
    return {_.name: _ for _ in resource.params}


def _get_accessurl_from_params(params):
    if "accessURL" not in params:
        raise DALServiceError("Datalink has no accessURL")
    return params["accessURL"].value


class AdhocServiceResultsMixin:
    """
    Mixin for adhoc:service functionality for results classes.
    """

    def __init__(self, votable, *, url=None, session=None):
        super().__init__(votable, url=url, session=session)
        self._adhocservices = list(
            resource for resource in votable.resources
            if resource.type == "meta" and resource.utype == "adhoc:service"
        )

    def iter_adhocservices(self):
        yield from self._adhocservices

    def get_adhocservice_by_ivoid(self, ivo_id):
        """
        Return the adhoc service starting with the given ivo_id.

        Parameters
        ----------
        ivoid : str
           the ivoid of the service we want to have.

        Returns
        -------
        Resource
            The resource element describing the service.
        """
        if isinstance(ivo_id, bytes):
            ivo_id = ivo_id.decode('utf-8')
        for adhocservice in self.iter_adhocservices():
            if any(
                    all((
                        param.name == "standardID",
                        param.value.lower().startswith(ivo_id.lower())
                    )) for param in adhocservice.params
            ):
                return adhocservice
        raise DALServiceError(
            "No Adhoc Service with ivo-id {}!".format(ivo_id))

    def get_adhocservice_by_id(self, id_):
        """
        Return the adhoc service starting with the given service_def id.

        Parameters
        ----------
        id_ : str
           the service_def id of the service we want to have.

        Returns
        -------
        Resource
            The resource element describing the service.
        """
        for adhocservice in self.iter_adhocservices():
            if adhocservice.ID == id_:
                return adhocservice
        raise DALServiceError(
            "No Adhoc Service with service_def id {}!".format(id_))


class DatalinkResultsMixin(AdhocServiceResultsMixin):
    """
    Mixin for datalink functionallity for results classes.
    """
    def _iter_datalinks_from_dlblock(self, datalink_service):
        """yields datalinks from the current rows using a datalink
        service RESOURCE.
        """
        remaining_ids = []  # remaining IDs to processed
        current_batch = None  # retrieved but not returned yet
        current_ids = []  # retrieved but not returned
        processed_ids = []  # retrived and returned IDs
        batch_size = None  # size of the batch

        for row in self:
            if not current_ids:
                if batch_size is None:
                    # first call.
                    self.query = DatalinkQuery.from_resource(
                        [_ for _ in self],
                        self._datalink,
                        session=self._session,
                        original_row=row)
                    remaining_ids = self.query['ID']
                if not remaining_ids:
                    # we are done
                    return
                if batch_size:
                    # subsequent calls are limitted to batch size
                    self.query['ID'] = remaining_ids[:batch_size]
                current_batch = self.query.execute(post=True)
                current_ids = list(OrderedDict.fromkeys(
                    [_ for _ in current_batch.to_table()['ID']]))
                if not current_ids:
                    raise DALServiceError(
                        'Could not retrieve datalinks for: {}'.format(
                            ', '.join([_ for _ in remaining_ids])))
                batch_size = len(current_ids)
            id1 = current_ids.pop(0)
            processed_ids.append(id1)
            remaining_ids.remove(id1)
            yield current_batch.clone_byid(
                id1,
                original_row=row)

    @staticmethod
    def _guess_access_format(row):
        """returns a guess for the format of what __guess_access_url will
        return.

        This tries a few heuristics based on how obscore or SIA records might
        be marked up.  If will return None if row does not look as if
        it contained an access format.  Note that the heuristics are
        tried in sequence; for now, we do not define the sequence of
        heuristics.
        """
        if hasattr(row, "access_format"):
            return row.access_format

        if "access_format" in row:
            return row["access_format"]

        access_format = row.getbyutype("obscore:access.format"
            ) or row.getbyutype("ssa:Access.Format")
        if access_format:
            return access_format

        access_format = row.getbyucd("meta.code.mime"
            ) or row.getbyucd("VOX:Image_Format")
        if access_format:
            return access_format

    @staticmethod
    def _guess_access_url(row):
        """returns a guess for a URI to a data product in row.

        This tries a few heuristics based on how obscore or SIA records might
        be marked up.  If will return None if row does not look as if
        it contained a product access url.  Note that the heuristics are
        tried in sequence; for now, we do not define the sequence of
        heuristics.
        """
        if hasattr(row, "access_url"):
            return row.access_url

        if "access_url" in row:
            return row["access_url"]

        access_url = row.getbyutype("obscore:access.reference"
            ) or row.getbyutype("ssa:Access.Reference")
        if access_url:
            return access_url

        access_url = row.getbyucd("meta.ref.url"
            ) or row.getbyucd("VOX:Image_AccessReference")
        if access_url:
            return access_url

    @staticmethod
    def _guess_datalink(row, **kwargs):
        # TODO: we should be more careful about whitespace, case
        # and perhaps more parameters in the following comparison
        if row._results._guess_access_format(row) == DATALINK_MIME_TYPE:
            access_url = row._results._guess_access_url(row)
            if access_url is not None:
                return DatalinkResults.from_result_url(access_url, **kwargs)

    def _iter_datalinks_from_product_rows(self):
        """yield datalinks from self's rows if they describe datalink-valued
        products.
        """
        for row in self:
            # TODO: we should be more careful about whitespace, case
            # and perhaps more parameters in the following comparison
            if self._guess_access_format(row) == DATALINK_MIME_TYPE:
                access_url = self._guess_access_url(row)
                if access_url is not None:
                    yield DatalinkResults.from_result_url(
                        access_url,
                        original_row=row)

    def iter_datalinks(self):
        """
        Iterates over all datalinks in a DALResult.
        """
        # To reduce the number of calls to the Datalink service, multiple
        # IDs are sent in batches. The appropriate batch size is not available
        # as each service has its own criteria for determining the maximum
        # number of IDs processed per request. To overcome this limitation,
        # this method initially sends all the available IDs and if the
        # response is partial, uses the size of the response as the batch
        # size.
        if not hasattr(self, '_datalink'):
            try:
                self._datalink = self.get_adhocservice_by_ivoid(DATALINK_IVOID)
            except DALServiceError:
                self._datalink = None

        if self._datalink is None:
            yield from self._iter_datalinks_from_product_rows()

        else:
            yield from self._iter_datalinks_from_dlblock(self._datalink)


class DatalinkRecordMixin:
    """
    Mixin for record classes, providing functionallity for datalink.

    - ``getdataset()`` considers datalink.
    """

    def getdatalink(self):
        """
        Retrieve the datalink information for this record.

        Returns
        -------
        DatalinkResults
            The datalink results for this record.

        Raises
        ------
        DALServiceError
            If no datalink information is found for this record.
        """
        try:
            datalink = self._results.get_adhocservice_by_ivoid(DATALINK_IVOID)

            query = DatalinkQuery.from_resource(self, datalink, session=self._session)
            return query.execute()
        except DALServiceError as error:
            datalink = self._results._guess_datalink(self, session=self._session)
            if datalink is not None:
                return datalink
            else:
                # re-raise the original error if nothing works
                raise DALServiceError("No datalink found for record.") from error

    @stream_decode_content
    def getdataset(self, timeout=None):
        try:
            url = next(self.getdatalink().bysemantics('#this')).access_url
            response = self._session.get(url, stream=True, timeout=timeout)
            try:
                response.raise_for_status()
            except requests.RequestException as ex:
                raise DALServiceError.from_except(ex, url)
            return response.raw
        except (DALServiceError, ValueError, StopIteration):
            # this should go to Record.getdataset()
            return super().getdataset(timeout=timeout)


class DatalinkService(DALService, AvailabilityMixin, CapabilityMixin):
    """
    a representation of a Datalink service
    """

    def __init__(self, baseurl, session=None):
        """
        instantiate a Datalink service

        Parameters
        ----------
        baseurl :  str
           the base URL that should be used for forming queries to the service.
        session : object
           optional session to use for network requests
        """
        super().__init__(baseurl, session=session)

        # Check if the session has an update_from_capabilities attribute.
        # This means that the session is aware of IVOA capabilities,
        # and can use this information in processing network requests.
        # One such usecase for this is auth.
        if hasattr(self._session, 'update_from_capabilities'):
            self._session.update_from_capabilities(self.capabilities)

    def run_sync(self, id, *, responseformat=None, **keywords):
        """
        runs sync query and returns its result

        Parameters
        ----------
        id : str
            the dataset identifier
        responseformat : str
            the output format

        Returns
        -------
        DatalinkResults
            the query result

        See Also
        --------
        DatalinkResults
        """
        return self.create_query(id, responseformat=responseformat, **keywords).execute()

    # alias for service discovery
    search = run_sync

    def create_query(self, id, *, responseformat=None, **keywords):
        """
        create a query object that constraints can be added to and then
        executed.  The input arguments will initialize the query with the
        given values.

        Parameters
        ----------
        baseurl : str
            the base URL for the Datalink service
        id : str
            the dataset identifier
        responseformat : str
            the output format
        """
        return DatalinkQuery(
            self.baseurl, id=id, responseformat=responseformat, **keywords)


class DatalinkQuery(DALQuery):
    """
    A class for preparing a query to a Datalink service.  Query constraints
    are added via its service type-specific methods.  The various execute()
    functions will submit the query and return the results.

    The base URL for the query, which controls where the query will be sent
    when one of the execute functions is called, is typically set at
    construction time; however, it can be updated later via the
    :py:attr:`~pyvo.dal.DALQuery.baseurl` to send a configured
    query to another service.

    A session can also optionally be passed in that will be used for
    network transactions made by this object to remote services.

    In addition to the search constraint attributes described below, search
    parameters can be set generically by name via dict semantics.
    The typical function for submitting the query is ``execute()``; however,
    alternate execute functions provide the response in different forms,
    allowing the caller to take greater control of the result processing.
    """
    @classmethod
    def from_resource(cls, rows, resource, *, session=None, **kwargs):
        """
        Creates a instance from a number of records and a Datalink Resource.

        XML Hierarchy:

        .. code:: xml

            <FIELD id="ID">
            <FIELD id="srcGroup">

            <GROUP name="inputParams">
                <PARAM name="ID" datatype="char" arraysize="*" value=""
                    ref="primaryID"/>
                <PARAM name="CALIB" datatype="char" arraysize="*"
                    value="FLUX"/>
                <PARAM name="GROUP" datatype="char" arraysize="*" value=""
                    ref="srcGroup"/>
            </GROUP>
        """
        original_row = kwargs.pop("original_row", None)

        input_params = _get_input_params_from_resource(resource)
        # get params outside of any group
        dl_params = _get_params_from_resource(resource)
        accessurl = _get_accessurl_from_params(dl_params)

        query_params = dict()
        for name, input_param in input_params.items():
            if input_param.ref:
                if isinstance(rows, list):
                    query_params[name] = []
                    for r in rows:
                        query_params[name].append(r[input_param.ref])
                else:
                    # scalars are also accepted for backwards compatibility
                    query_params[name] = rows[input_param.ref]
            elif np.isscalar(input_param.value) and input_param.value:
                query_params[name] = input_param.value
            elif (
                    not np.isscalar(input_param.value)
                    and input_param.value.all()
                    and len(input_param.value)
            ):
                query_params[name] = " ".join(
                    str(_) for _ in input_param.value)

        for name, query_param in kwargs.items():
            try:
                input_param = find_param_by_keyword(name, input_params)
                if input_param and query_param is None:
                    del query_params[input_param.name]
                converter = get_converter(input_param)
                query_params[input_param.name] = converter.serialize(
                    query_param)
            except KeyError:
                query_params[name] = query_param

        return cls(
            accessurl,
            session=session,
            original_row=original_row,
            **query_params)

    def __init__(
            self, baseurl, *, id=None, responseformat=None, session=None, **keywords):
        """
        initialize the query object with the given parameters

        Parameters
        ----------
        baseurl : str
            the Datalink baseurl
        id : str
            the dataset identifier
        responseformat : str
            the output format
        session : object
            optional session to use for network requests
        """
        self.original_row = keywords.pop("original_row", None)

        super().__init__(baseurl, session=session, **keywords)

        if id is not None:
            self["ID"] = id
        if responseformat is not None:
            self["RESPONSEFORMAT"] = responseformat

    def execute(self, post=False):
        """
        submit the query and return the results as a DatalinkResults instance

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
        return DatalinkResults(
            self.execute_votable(post=post),
            url=self.queryurl,
            original_row=self.original_row,
            session=self._session)


class DatalinkResults(DatalinkResultsMixin, DALResults):
    """
    The list of matching records resulting from an datalink query.
    Each record contains a set of metadata that describes an available
    record matching the query constraints.  The number of records in
    the results is available by passing it to the Python built-in ``len()``
    function.

    This class supports iterable semantics; thus,
    individual records (in the form of
    :py:class:`~pyvo.dal.Record` instances) are typically
    accessed by iterating over an ``DatalinkResults`` instance.

    Alternatively, records can be accessed randomly via
    :py:meth:`getrecord` or through a Python Database API (v2)
    Cursor (via :py:meth:`~pyvo.dal.DALResults.cursor`).
    Column-based data access is possible via the
    :py:meth:`~pyvo.dal.DALResults.getcolumn` method.

    ``DatalinkResults`` is essentially a wrapper around an Astropy
    :py:mod:`~astropy.io.votable`
    :py:class:`~astropy.io.votable.tree.TableElement` instance where the
    columns contain the various metadata describing the images.
    One can access that VOTable directly via the
    :py:attr:`~pyvo.dal.DALResults.votable` attribute.  Thus,
    when one retrieves a whole column via
    :py:meth:`~pyvo.dal.DALResults.getcolumn`, the result is
    a Numpy array.  Alternatively, one can manipulate the results
    as an Astropy :py:class:`~astropy.table.Table` via the
    following conversion:

    ``table = results.to_table()``

    ``DatalinkResults`` supports the array item operator ``[...]`` in a
    read-only context.  When the argument is numerical, the result
    is an
    :py:class:`~pyvo.dal.Record` instance, representing the
    record at the position given by the numerical index.  If the
    argument is a string, it is interpreted as the name of a column,
    and the data from the column matching that name is returned as
    a Numpy array.
    """

    def __init__(self, *args, **kwargs):
        self.original_row = kwargs.pop("original_row", None)
        super().__init__(*args, **kwargs)

    def getrecord(self, index):
        """
        return a representation of a datalink result record that follows
        dictionary semantics. The keys of the dictionary are those returned by
        this instance's fieldnames attribute. The returned record has the
        additional function ``getdataset``

        Parameters
        ----------
        index : int
           the integer index of the desired record where 0 returns the first
           record

        Returns
        -------
        Rec
           a dictionary-like wrapper containing the result record metadata.

        Raises
        ------
        IndexError
           if index is negative or equal or larger than the number of rows in
           the result table.

        See Also
        --------
        pyvo.dal.Record
        """
        return DatalinkRecord(self, index, session=self._session)

    def bysemantics(self, semantics, *, include_narrower=True):
        """
        return the rows with the dataset identified by the given semantics

        Parameters
        ----------
        semantics: str or list
            One or more term(s) from the datalink vocabulary
            (http://www.ivoa.net/rdf/datalink/core).  datalink/core
            terms may be passed in with or without a leading hash.
            Free URIs may also be passed an and will be compared literally,
            i.e., without any URI normalisation.
        include_narrower: boolean
            If true, the result will include matches for any term
            that is narrower than the term passed in.

        Returns
        -------
        Sequence of DatalinkRecord
            a sequence of dictionary-like wrappers containing the result record
        """
        # If the URL juggling gets any more complicated here, we ought
        # to bite the bullet and only deal with full URLs.  Sigh.
        if isinstance(semantics, str):
            semantics = [semantics]

        core_terms, other_terms = [], []
        for term in semantics:
            if ":" in term:
                # it's a full URI, see if it's ours
                if term.startswith("http://www.ivoa.net/rdf/datalink/core#"):
                    core_terms.append(term.split("#", 1)[-1])
                else:
                    other_terms.append(term)
            else:
                # it's a local term
                core_terms.append(term.lstrip("#"))

        if include_narrower:
            additional_terms = []
            voc = vocabularies.get_vocabulary("datalink/core")
            for term in core_terms:
                if term in voc["terms"]:
                    additional_terms.extend(voc["terms"][term]["narrower"])
            core_terms = core_terms + additional_terms

        semantics = set("#" + term for term in core_terms) | set(other_terms)
        for record in self:
            if record.semantics in semantics:
                yield record

    def clone_byid(self, id, *, original_row=None):
        """
        return a clone of the object with results and corresponding
        resources matching a given id

        Returns
        -------
        Sequence of DatalinkRecord
            a sequence of dictionary-like wrappers containing the result record
        """

        copy_tb = copy.deepcopy(self.votable)
        votable = copy_tb.get_first_table()
        # find index of ID column
        id_index = None
        for index, field in enumerate(votable.fields):
            if field.name == 'ID':
                id_index = index
        rows = [x for x in votable.array if x[id_index] == id]
        votable.create_arrays(len(rows))
        for index, row in enumerate(rows):
            votable.array[index] = row
        # now remove unreferenced services from resources
        referenced_serviced = [x for x in votable.array['service_def'] if x]
        # remove customized that are not referenced by the current results
        for x in copy_tb.resources:
            if x.ID and x.ID not in referenced_serviced:
                copy_tb.resources.remove(x)
        return DatalinkResults(copy_tb, original_row=original_row)

    def getdataset(self, *, timeout=None):
        """
        return the first row with the dataset identified by semantics #this

        Returns
        -------
        DatalinkRecord
            a dictionary-like wrapper containing the result record.
        """
        try:
            return next(self.bysemantics("#this")).getdataset(timeout=timeout)
        except StopIteration:
            raise ValueError("No row with semantics #this found!")

    def iter_procs(self):
        """
        iterate over all rows with a processing service
        """
        for row in self:
            if row.service_def:
                yield row

    def get_first_proc(self):
        """
        returns the first datalink row with a processing service.
        """
        for proc in self.iter_procs():
            return proc
        raise IndexError("No processing service found in datalink result")

    @classmethod
    def from_result_url(cls, result_url, *, session=None, original_row=None):
        res = super().from_result_url(result_url, session=session)
        res.original_row = original_row
        return res


class SodaRecordMixin:
    """
    Mixin for soda functionality for record classes.
    If used, it's result class must have
    `pyvo.dal.adhoc.AdhocServiceResultsMixin` mixed in.
    """

    def _get_soda_resource(self):
        try:
            return self._results.get_adhocservice_by_ivoid(SODA_SYNC_IVOID)
        except DALServiceError:
            pass

        dataformat = self.getdataformat()
        if dataformat is None:
            raise DALServiceError(
                "No SODA Resouces available and no dataformat in record. "
                "Maybe you forgot to include it into the TAP Query?")

        if "content=datalink" in dataformat:
            dataurl = self.getdataurl()
            if not dataurl:
                raise DALServiceError(
                    "No dataurl in record, but dataformat contains "
                    "'content=datalink'. Maybe you forgot to include it into "
                    "the TAP Query?")

            try:
                datalink_result = DatalinkResults.from_result_url(dataurl)
                return datalink_result.get_adhocservice_by_ivoid(
                    SODA_SYNC_IVOID)
            except DALServiceError:
                pass

        return None

    def processed(
            self, *, circle=None, range=None, polygon=None, band=None, **kwargs):
        """
        Returns processed dataset.

        Parameters
        ----------
        circle : `astropy.units.Quantity`
            latitude, longitude and radius
        range : `astropy.units.Quantity`
            two longitude + two latitude values describing a rectangle
        polygon : `astropy.units.Quantity`
            multiple (at least three) pairs of longitude and latitude points
        band : `astropy.units.Quantity`
            two bandwidth or frequency values
        """
        soda_resource = self._get_soda_resource()

        if soda_resource:
            soda_query = SodaQuery.from_resource(
                self, soda_resource, circle=circle, range=range,
                polygon=polygon, band=band, **kwargs)

            soda_stream = soda_query.execute_stream()
            soda_query.raise_if_error()
            return soda_stream
        else:
            return self.getdataset()


class DatalinkRecord(DatalinkRecordMixin, SodaRecordMixin, Record):
    """
    a dictionary-like container for data in a record from the results of an
    datalink query,

    The commonly accessed metadata which are stadardized by the datalink
    standard are available as attributes.  If the metadatum accessible
    via an attribute is not available, the value of that attribute
    will be None.  All metadata, including non-standard metadata, are
    acessible via the ``get(`` *key* ``)`` function (or the [*key*]
    operator) where *key* is table column name.
    """

    @property
    def id(self):
        """
        Input identifier
        """
        return self.get("ID", decode=True)

    @property
    def access_url(self):
        """
        Link to data or processing service
        """
        row_url = self.get("access_url", decode=True)

        if not row_url:
            proc_resource = self._results.get_adhocservice_by_id(
                self.service_def)
            dl_params = _get_params_from_resource(proc_resource)
            row_url = _get_accessurl_from_params(dl_params)

        return row_url

    @property
    def service_def(self):
        """
        reference to the service descriptor resource
        """
        return self.get("service_def", decode=True)

    @property
    def error_message(self):
        """
        Error if an access_url cannot be created
        """
        return self.get("error_message", decode=True)

    @property
    def description(self):
        """
        Human-readable text describing this link
        """
        return self.get("description", decode=True)

    @property
    def semantics(self):
        """
         Term from a controlled vocabulary describing the link
        """
        return self.get("semantics", decode=True)

    @property
    def content_type(self):
        """
        Mime-type of the content the link returns
        """
        return self.get("content_type", decode=True)

    @property
    def content_length(self):
        """
        Size of the download the link returns
        """
        return int(self["content_length"])

    def getdataurl(self):
        """
        return the URL contained in the access URL column which can be used
        to retrieve the dataset described by this record. Raises
        :py:class:`~pyvo.dal.DALServiceError` if theres an error.
        """
        if self.error_message:
            raise DALServiceError(self.error_message)

        return self.access_url

    def process(self, **kwargs):
        """
        calls the processing service and returns it's result as a file-like
        object
        """
        proc_resource = self._results.get_adhocservice_by_id(self.service_def)
        proc_query = DatalinkQuery.from_resource(self, proc_resource, **kwargs)
        proc_stream = proc_query.execute_stream()
        return proc_stream

    @property
    def params(self):
        """
        the access parameters of the service behind this datalink row.
        """
        proc_resource = self._results.get_adhocservice_by_id(self.service_def)
        return proc_resource.params

    @property
    def input_params(self):
        """
        a list of input parameters for the service behind this datalink row.
        """
        proc_resource = self._results.get_adhocservice_by_id(self.service_def)
        return list(_get_input_params_from_resource(proc_resource).values())


class AxisParamMixin():
    """
    Stores the axis parameters (pos, band, time and pol) used in SODA
    or SIA2 queries
    """
    @property
    def pos(self):
        if not hasattr(self, '_pos'):
            self._pos = PosQueryParam()
            self['POS'] = self._pos.dal
        return self._pos

    @property
    def band(self):
        if not hasattr(self, '_band'):
            self._band = IntervalQueryParam(unit=u.meter,
                                            equivalencies=u.spectral())
            self['BAND'] = self.band.dal
        return self._band

    @property
    def time(self):
        if not hasattr(self, '_time'):
            self._time = TimeQueryParam()
            self['TIME'] = self.time.dal
        return self._time

    @property
    def pol(self):
        if not hasattr(self, '_pol'):
            self._pol = EnumQueryParam(POLARIZATION_STATES)
            self['POL'] = self.pol.dal
        return self._pol


class SodaQuery(DatalinkQuery, AxisParamMixin):
    """
    a class for preparing a query to a SODA Service.
    """

    def __init__(
            self, baseurl, *, circle=None, range=None, polygon=None, band=None,
            **kwargs):
        super().__init__(baseurl, **kwargs)

        if circle is not None:
            self.circle = circle

        if range is not None:
            self.range = range

        if polygon is not None:
            self.polygon = polygon

        if band is not None:
            self.band = band

    @property
    def circle(self):
        """
        The CIRCLE parameter defines a spatial region using the circle xtype
        defined in DALI.
        """
        return getattr(self, '_circle', None)

    @circle.setter
    def circle(self, circle):
        if len(circle) != 3:
            raise ValueError(
                "Range must be a sequence with exactly three values")
        self['CIRCLE'] = PosQueryParam().get_dal_format(circle).\
            replace('CIRCLE ', '')
        setattr(self, '_circle', circle)
        del self.range
        del self.polygon

    @circle.deleter
    def circle(self):
        if hasattr(self, '_circle'):
            delattr(self, '_circle')
        if 'CIRCLE' in self:
            del self['CIRCLE']

    @property
    def range(self):
        """
        A rectangular range.
        """
        warnings.warn(
            "Use pos attribute instead",
            DeprecationWarning
        )
        return getattr(self, '_circle', None)

    @range.setter
    def range(self, range):
        if len(range) != 4:
            raise ValueError(
                "Range must be a sequence with exactly four values")
        self['POS'] = PosQueryParam().get_dal_format(range)
        setattr(self, '_range', range)
        del self.circle
        del self.polygon

    @range.deleter
    def range(self):

        if hasattr(self, '_range'):
            delattr(self, '_range')
        if 'POS' in self and self['POS'].startswith('RANGE'):
            del self['POS']

    @property
    def polygon(self):
        """
        The POLYGON parameter defines a spatial region using the polygon xtype
        defined in DALI.
        """
        return getattr(self, '_polygon', None)

    @polygon.setter
    def polygon(self, polygon):
        if len(polygon) < 6:
            raise ValueError(
                "Polygon must be a sequence of at least six values")
        self['POLYGON'] = PosQueryParam().get_dal_format(polygon).\
            replace('POLYGON ', '')
        setattr(self, '_polygon', polygon)
        del self.circle
        del self.range

    @polygon.deleter
    def polygon(self):
        if hasattr(self, '_polygon'):
            delattr(self, '_polygon')
        if 'POLYGON' in self:
            del self['POLYGON']

    @property
    def band(self):
        """
        The BAND parameter defines the wavelength interval(s) to be extracted
        from the data using a floating point interval
        """
        return getattr(self, "_band", None)

    @band.setter
    def band(self, band):
        setattr(self, "_band", band)

        if not isinstance(band, Quantity):
            valerr = ValueError(
                'Band must be specified with exactly two values, ',
                'expressing a frequency or wavelength range'
            )

            try:
                # assume meters
                band = band * Unit("meter")
                if len(band) != 2:
                    raise valerr
            except (ValueError, TypeError):
                raise valerr

        # transform to meters
        band = band.to(Unit("m"), equivalencies=spectral_equivalencies())
        # frequency is counter-proportional to wavelength, so we just sort
        # it to have the right order again
        band.sort()

        self["BAND"] = "{start} {end}".format(
            start=band.value[0], end=band.value[1])

    @band.deleter
    def band(self):
        if hasattr(self, '_band'):
            delattr(self, '_band')
        if 'BAND' in self:
            del self['BAND']
