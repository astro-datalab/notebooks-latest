# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module for searching for spectral line metadata in a remote database.

A Simple Line Access (SLA) service allows a client to search for
metadata describing atomic and molecular transitions that can result
in spectral line emission and absorption.  The service responds to a
search query with a table in which each row represents a transition
that matches the query constraints.  The columns provide the metadata
describing the transition.  This module provides an interface for
accessing an SLA service.  It is implemented as a specialization of
the DAL Query interface.

The ``search()`` function support the simplest and most common types
of queries, returning an SLAResults instance as its results which
represents the matching imagess from the archive.  The SLAResults
supports access to and iterations over the individual records; these
are provided as SLARecord instances, which give easy access to key
metadata in the response, such as the transition title.

The SLAService class can represent a specific service available at a URL
endpoint.
"""
from pyvo.io.vosi.vodataservice import TableParam

from astropy.units import Quantity, Unit
from astropy.units import spectral as spectral_equivalencies
from astropy.io.votable.tree import Field
from astropy.table import Table

from .query import DALResults, DALQuery, DALService, Record

__all__ = ["search", "SLAService", "SLAQuery", "SLAResults", "SLARecord"]


def search(baseurl, wavelength, **keywords):
    """
    submit a simple SLA query that requests spectral lines within a
    wavelength range

    Parameters
    ----------
    baseurl : str
       the base URL for the SLA service
    wavelength : `~astropy.units.Quantity` class or sequence of two floats
        the bandwidth range the observations belong to.
        assuming meters if unit is not specified.
    **keywords :
        additional parameters can be given via arbitrary
        case insensitive keyword arguments. Where there is overlap
        with the parameters set by the other arguments to
        this function, these keywords will override.

    Returns
    -------
    SLAResults
        a container holding a table of matching spectral lines

    Raises
    ------
    DALServiceError
       for errors connecting to or communicating with the service
    DALQueryError
       if the service responds with an error, including a query syntax error.
    """
    service = SLAService(baseurl)
    return service.search(wavelength, **keywords)


class SLAService(DALService):
    """
    a representation of an spectral line catalog (SLA) service
    """
    def __init__(self, baseurl, *, capability_description=None, session=None):
        """
        instantiate an SLA service

        Parameters
        ----------
        baseurl : str
           the base URL for submitting search queries to the service.
        session : object
           optional session to use for network requests
        """
        super().__init__(baseurl, capability_description=capability_description, session=session)

    def _get_metadata(self):
        """
        download the metadata resource
        """
        if not hasattr(self, "_metadata"):
            query = self.create_query(request='getCapabilities')
            metadata = query.execute_votable()

            setattr(self, "_metadata", metadata)

    @property
    def description(self):
        """
        the service description.

        If this is not provided during instantiation, this method will download
        a sample from the service and read the description in the sample's
        metadata instead.
        """
        if self._description is not None:
            return self._description

        self._get_metadata()
        try:
            return getattr(self, "_metadata", None).description
        except AttributeError:
            return None

    @property
    def columns(self):
        """
        the available columns on this service
        """
        self._get_metadata()
        fields = filter(
            lambda field_or_param: isinstance(field_or_param, Field),
            self._metadata.iter_fields_and_params()
        )

        try:
            return [
                TableParam.from_field(field) for field in fields]
        except AttributeError:
            return []

    def search(self, wavelength, **keywords):
        """
        submit a simple SLA query to this service with the given constraints.

        This method is provided for a simple but typical SLA queries.  For
        more complex queries, one should create an SLAQuery object via
        create_query()

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity` class or sequence of two floats
            the bandwidth range the observations belong to.
            assuming meters if unit is not specified.
        **keywords :
            additional parameters can be given via arbitrary
            case insensitive keyword arguments. Where there is overlap
            with the parameters set by the other arguments to
            this function, these keywords will override.

        Returns
        -------
        SLAResults
            a container holding a table of matching spectral lines

        Raises
        ------
        DALServiceError
           for errors connecting to or communicating with the service
        DALQueryError
           if the service responds with
           an error, including a query syntax error.

        See Also
        --------
        SLAResults
        pyvo.dal.DALServiceError
        pyvo.dal.DALQueryError
        """
        return self.create_query(wavelength, **keywords).execute()

    def create_query(self, wavelength=None, *, request="queryData", **keywords):
        """
        create a query object that constraints can be added to and then
        executed.  The input arguments will initialize the query with the
        given values.

        Parameters
        ----------
        wavelength : `~astropy.units.Quantity` class or sequence of two floats
            the bandwidth range the observations belong to.
            assuming meters if unit is not specified.
        **keywords :
            additional parameters can be given via arbitrary
            case insensitive keyword arguments. Where there is overlap
            with the parameters set by the other arguments to
            this function, these keywords will override.

        Returns
        -------
        SLAQuery
           the query instance

        See Also
        --------
        SLAQuery
        """
        return SLAQuery(baseurl=self.baseurl, wavelength=wavelength, request=request,
                        session=self._session, **keywords)

    def describe(self):
        print(self.description)
        print()

        rows = [(
            col.name,
            col.description,
            col.unit,
            col.ucd,
            col.utype,
            col.datatype.arraysize,
            col.datatype.content,
        ) for col in self.columns]

        names = (
            'name',
            'description',
            'unit',
            'ucd',
            'utype',
            'arraysize',
            'datatype',
        )

        table = Table(rows=rows, names=names)
        table.pprint(
            max_lines=-1, max_width=-1, show_unit=False, show_dtype=False)


class SLAQuery(DALQuery):
    """
    a class for preparing a query to an SLA service.  Query constraints
    are added via its service type-specific methods.  The various execute()
    functions will submit the query and return the results.

    The base URL for the query, which controls where the query will be sent
    when one of the execute functions is called, is typically set at
    construction time; however, it can be updated later via the
    :py:attr:`~pyvo.dal.DALQuery.baseurl` to send a configured
    query to another service.

    In addition to the search constraint attributes described below, search
    parameters can be set generically by name via the dict semantics.

    The typical function for submitting the query is ``execute()``; however,
    alternate execute functions provide the response in different forms,
    allowing the caller to take greater control of the result processing.
    """

    def __init__(
            self, baseurl, wavelength=None, *, request="queryData", session=None):
        """
        initialize the query object with a baseurl and the given parameters

        Parameters
        ----------
        baseurl : str
            the base URL for the SLA service
        wavelength : `~astropy.units.Quantity` class or sequence of two floats
            the bandwidth range the observations belong to.
            assuming meters if unit is not specified.
        session : object
           optional session to use for network requests
        """
        super().__init__(baseurl, session=session)

        if wavelength is not None:
            self.wavelength = wavelength

        self.request = request

    @property
    def wavelength(self):
        """
        the frequency/wavelength range the observations belong to.
        """
        return getattr(self, "_wavelength", None)

    @wavelength.setter
    def wavelength(self, wavelength):
        setattr(self, "_wavelength", wavelength)

        if not isinstance(wavelength, Quantity):
            valerr = ValueError(
                'Wavelength range must be a sequence with exactly two values',
                'expressing a frequency or wavelength range')

            try:
                # assume meters
                wavelength = wavelength * Unit("meter")
            except ValueError:
                raise valerr

            try:
                if len(wavelength) != 2:
                    raise valerr
            except TypeError:
                raise valerr

        # transform to meters
        wavelength = wavelength.to(
            Unit("m"), equivalencies=spectral_equivalencies())
        # frequency is counter-proportional to wavelength, so we just sort it
        # to have the right order again
        wavelength.sort()

        self["WAVELENGTH"] = "{start}/{end}".format(
            start=wavelength.value[0], end=wavelength.value[1])

    @wavelength.deleter
    def wavelength(self):
        delattr(self, "_wavelength")
        del self["WAVELENGTH"]

    @property
    def request(self):
        """
        the type of service operation which is being performed
        """
        return getattr(self, "_request", None)

    @request.setter
    def request(self, val):
        setattr(self, "_request", val)
        self["REQUEST"] = val

    @request.deleter
    def request(self):
        delattr(self, "_request")
        del self["REQUEST"]

    def execute(self):
        """
        submit the query and return the results as a SLAResults instance

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
        return SLAResults(self.execute_votable(), url=self.queryurl, session=self._session)


class SLAResults(DALResults):
    """
    The list of matching spectral lines resulting from a spectal line
    catalog (SLA) query.
    Each record contains a set of metadata that describes a source or
    observation within the requested circular region (i.e. a "cone").  The
    number of records in the results is available by passing it to the Python built-in
    ``len()`` function.

    This class supports iterable semantics; thus,
    individual records (in the form of
    :py:class:`~pyvo.dal.sla.SLARecord` instances) are typically
    accessed by iterating over an ``SLAResults`` instance.

    Alternatively, records can be accessed randomly via
    :py:meth:`getrecord` or through a Python Database API (v2)
    Cursor (via :py:meth:`~pyvo.dal.DALResults.cursor`).
    Column-based data access is possible via the
    :py:meth:`~pyvo.dal.DALResults.getcolumn` method.

    ``SLAResults`` is essentially a wrapper around an Astropy
    :py:mod:`~astropy.io.votable`
    :py:class:`~astropy.io.votable.tree.TableElement` instance where the
    columns contain the various metadata describing the images.
    One can access that VOTable directly via the
    :py:attr:`~pyvo.dal.DALResults.votable` attribute.  Thus,
    when one retrieves a whole column via
    :py:meth:`~pyvo.dal.DALResults.getcolumn`, the result is
    a Numpy array.  Alternatively, one can manipulate the results
    as an Astropy :py:class:`~astropy.table.table.Table` via the
    following conversion:

    ``table = results.votable.to_table()``

    ``SLAResults`` supports the array item operator ``[...]`` in a
    read-only context.  When the argument is numerical, the result
    is an
    :py:class:`~pyvo.dal.sla.SLARecord` instance, representing the
    record at the position given by the numerical index.  If the
    argument is a string, it is interpreted as the name of a column,
    and the data from the column matching that name is returned as
    a Numpy array.
    """

    def getrecord(self, index):
        """
        return a representation of a sla result record that follows
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
        SLARecord
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
        return SLARecord(self, index, session=self._session)


class SLARecord(Record):
    """
    a dictionary-like container for data in a record from the results of an
    spectral line (SLA) query, describing a spectral line transition.

    The commonly accessed metadata which are stadardized by the SLA
    protocol are available as attributes.  All metadata, particularly
    non-standard metadata, are acessible via the ``get(`` *key* ``)``
    function (or the [*key*] operator) where *key* is table column name.
    """

    @property
    def title(self):
        """
        a title/small description of the line transition
        """
        return self.getbyutype("ssldm:Line.title", decode=True)

    @property
    def wavelength(self):
        """
        the vacuum wavelength of the line in meters.
        """
        return self.getbyutype("ssldm:Line.wavelength.value") * Unit("m")

    @property
    def species_name(self):
        """
        the name of the chemical species that produces the transition.
        """
        return self.getbyutype("ssldm:Line.species.name")

    @property
    def status(self):
        """
        the name of the chemical species that produces the transition.
        """
        return self.getbyutype("ssldm:Line.identificationStatus")

    @property
    def initial_level(self):
        """
        a description of the initial (higher energy) quantum level
        """
        return self.getbyutype("ssldm:Line.initialLevel.name", decode=True)

    @property
    def final_level(self):
        """
        a description of the final (higher energy) quantum level
        """
        return self.getbyutype("ssldm:Line.finalLevel.name")
