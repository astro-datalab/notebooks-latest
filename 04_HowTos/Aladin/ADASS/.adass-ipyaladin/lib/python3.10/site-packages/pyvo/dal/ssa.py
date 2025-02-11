# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module for searching for spectra in a remote archive.

A Simple Spectral Access (SSA) service allows a client to search for
spectra in an archive whose field of view overlaps with a given cone
on the sky.  The service responds to a search query with a table in
which each row represents an image that is available for download.
The columns provide metadata describing each image and one column in
particular provides the image's download URL (also called the *access
reference*, or *acref*).  Some SSA services can create spectra
on-the-fly from underlying data (e.g. image cubes); in this case, the
query result is a table of images whose aperture matches the
requested cone and which will be created when accessed via the
download URL.

This module provides an interface for accessing an SSA service.  It is
implemented as a specialization of the DAL Query interface.

The ``search()`` function support the simplest and most common types
of queries, returning an SSAResults instance as its results which
represents the matching imagess from the archive.  The SSAResults
supports access to and iterations over the individual records; these
are provided as SSARecord instances, which give easy access to key
metadata in the response, such as the position of the spectrum's
aperture, the spectrum format, its frequency range, and its download
URL.

The SSAService class can represent a specific service available at a URL
endpoint.
"""
import re

from pyvo.io.vosi.vodataservice import TableParam

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.units import Quantity, Unit
from astropy.units import spectral as spectral_equivalencies
from astropy.io.votable.tree import Field
from astropy.table import Table
import numpy as np

from .query import DALResults, DALQuery, DALService, Record
from .mimetype import mime2extension
from .adhoc import DatalinkResultsMixin, DatalinkRecordMixin, SodaRecordMixin

from .. import samp

__all__ = ["search", "SSAService", "SSAQuery", "SSAResults", "SSARecord"]


def search(
        baseurl, pos=None, *, diameter=None, band=None, time=None, format=None,
        **keywords):
    """
    submit a simple SSA query that requests spectra overlapping a given region

    Parameters
    ----------
    baseurl : str
        the base URL for the SSA service
    pos : `~astropy.coordinates.SkyCoord` class or sequence of two floats
        the position of the center of the circular search region.
        assuming icrs decimal degrees if unit is not specified.
    diameter : `~astropy.units.Quantity` class or scalar float
        the diameter of the circular region around pos in which to search.
        assuming icrs decimal degrees if unit is not specified.
    band : `~astropy.units.Quantity` class or sequence of two floats
        the bandwidth range the observations belong to.
        assuming meters if unit is not specified.
    time : `~astropy.time.Time` class or sequence of two strings
        the datetime range the observations were made in.
        assuming iso 8601 if format is not specified.
    format : str
        the image format(s) of interest.  "all" indicates
        all available formats; "graphic" indicates
        graphical images (e.g. jpeg, png, gif; not FITS);
        "metadata" indicates that no images should be
        returned--only an empty table with complete metadata.
    **keywords :
        additional case insensitive parameters can be given via arbitrary
        case insensitive keyword arguments. Where there is overlap
        with the parameters set by the other arguments to
        this function, these keywords will override.

    Returns
    -------
    SSAResults
        a container holding a table of matching spectrum records

    Raises
    ------
    DALServiceError
       for errors connecting to or communicating with the service
    DALQueryError
       if the service responds with an error, including a query syntax error.

    See Also
    --------
    SSAResults
    pyvo.dal.DALServiceError
    pyvo.dal.DALQueryError
    """
    return SSAService(baseurl).search(
        pos=pos, diameter=diameter, band=band, time=time, format=format, **keywords)


class SSAService(DALService):
    """
    a representation of an SSA service
    """
    def __init__(self, baseurl, *, capability_description=None, session=None):
        """
        instantiate an SSA service

        Parameters
        ----------
        baseurl : str
           the base URL for submitting search queries to the service.
        """
        super().__init__(baseurl, session=session, capability_description=capability_description)

    def _get_metadata(self):
        """
        the metadata resource element
        """
        if not hasattr(self, "_metadata"):
            query = self.create_query(format='metadata')
            metadata = query.execute_votable()

            setattr(self, "_metadata", metadata)

    @property
    def description(self):
        """
        the service description.
        """
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

    def search(
            self, pos=None, *, diameter=None, band=None, time=None, format=None,
            **keywords):
        """
        submit a SSA query to this service with the given constraints.

        Parameters
        ----------
        pos : `~astropy.coordinates.SkyCoord` class or sequence of two floats
            the position of the center of the circular search region.
            assuming icrs decimal degrees if unit is not specified.
        diameter : `~astropy.units.Quantity` class or scalar float
            the diameter of the circular region around pos in which to search.
            assuming icrs decimal degrees if unit is not specified.
        band : `~astropy.units.Quantity` class or sequence of two floats
            the bandwidth range the observations belong to.
            assuming meters if unit is not specified.
        time : `~astropy.time.Time` class or sequence of two strings
            the datetime range the observations were made in.
            assuming iso 8601 if format is not specified.
        format : str
           the image format(s) of interest.  "all" indicates
           all available formats; "graphic" indicates
           graphical images (e.g. jpeg, png, gif; not FITS);
           "metadata" indicates that no images should be
           returned--only an empty table with complete metadata.
        **keywords :
           additional case insensitive parameters can be given via arbitrary
           case insensitive keyword arguments. Where there is overlap
           with the parameters set by the other arguments to
           this function, these keywords will override.

        Returns
        -------
        SSAResults
           a container holding a table of matching catalog records

        Raises
        ------
        DALServiceError
           for errors connecting to or communicating with the service
        DALQueryError
           if the service responds with an error, including query syntax errors

        See Also
        --------
        SSAResults
        pyvo.dal.DALServiceError
        pyvo.dal.DALQueryError
        """
        return self.create_query(
            pos=pos, diameter=diameter, band=band, time=time, format=format, **keywords).execute()

    def create_query(
            self, pos=None, *, diameter=None, band=None, time=None, format=None,
            request="queryData", **keywords):
        """
        create a query object that constraints can be added to and then
        executed.  The input arguments will initialize the query with the
        given values.

        Parameters
        ----------
        pos : `~astropy.coordinates.SkyCoord` class or sequence of two floats
            the position of the center of the circular search region.
            assuming icrs decimal degrees if unit is not specified.
        diameter : `~astropy.units.Quantity` class or scalar float
            the diameter of the circular region around pos in which to search.
            assuming icrs decimal degrees if unit is not specified.
        band : `~astropy.units.Quantity` class or sequence of two floats
            the bandwidth range the observations belong to.
            assuming meters if unit is not specified.
        time : `~astropy.time.Time` class or sequence of two strings
            the datetime range the observations were made in.
            assuming iso 8601 if format is not specified.
        format : str
           the image format(s) of interest.  "all" indicates
           all available formats; "graphic" indicates
           graphical images (e.g. jpeg, png, gif; not FITS);
           "metadata" indicates that no images should be
           returned--only an empty table with complete metadata.
        **keywords :
           additional case insensitive parameters can be given via arbitrary
           case insensitive keyword arguments. Where there is overlap
           with the parameters set by the other arguments to
           this function, these keywords will override.

        Returns
        -------
        SSAQuery
            the query instance

        See Also
        --------
        SSAQuery
        """
        return SSAQuery(
            self.baseurl, pos=pos, diameter=diameter, band=band, time=time, format=format, request=request,
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


class SSAQuery(DALQuery):
    """
    a class for preparing a query to an SSA service.  Query constraints
    are added via its service type-specific properties and methods.  Once
    all the constraints are set, one of the various execute() functions
    can be called to submit the query and return the results.

    The base URL for the query, which controls where the query will be sent
    when one of the execute functions is called, is typically set at
    construction time; however, it can be updated later via the
    :py:attr:`~pyvo.dal.DALQuery.baseurl` to send a configured
    query to another service.

    The typical function for submitting the query is ``execute()``; however,
    alternate execute functions provide the response in different forms,
    allowing the caller to take greater control of the result processing.

    """

    def __init__(
            self, baseurl, pos=None, *, diameter=None, band=None, time=None,
            format=None, request="queryData", session=None, **keywords):
        """
        initialize the query object with a baseurl and the given parameters

        Parameters
        ----------
        baseurl : str
            the base URL for the SSA service
        pos : `~astropy.coordinates.SkyCoord` class or sequence of two floats
            the position of the center of the circular search region.
            assuming icrs decimal degrees if unit is not specified.
        diameter : `~astropy.units.Quantity` class or scalar float
            the diameter of the circular region around pos in which to search.
            assuming icrs decimal degrees if unit is not specified.
        band : `~astropy.units.Quantity` class or sequence of two floats
            the bandwidth range the observations belong to.
            assuming meters if unit is not specified.
        time : `~astropy.time.Time` class or sequence of two strings
            the datetime range the observations were made in.
            assuming iso 8601 if format is not specified.
        format : str
           the image format(s) of interest.  "all" indicates
           all available formats; "graphic" indicates
           graphical images (e.g. jpeg, png, gif; not FITS);
           "metadata" indicates that no images should be
           returned--only an empty table with complete metadata.
        session : object
           optional session to use for network requests
        **keywords :
           additional case insensitive parameters can be given via arbitrary
           case insensitive keyword arguments. Where there is overlap
           with the parameters set by the other arguments to
           this function, these keywords will override.
        """
        super().__init__(baseurl, session=session)

        if pos is not None:
            self.pos = pos

        if diameter is not None:
            self.diameter = diameter

        if band is not None:
            self.band = band

        if time is not None:
            self.time = time

        if format is not None:
            self.format = format

        self.request = request

        self.update({key.upper(): value for key, value in keywords.items()})

    @property
    def pos(self):
        """
        the position of the center of the circular search region as a
        `~astropy.coordinates.SkyCoord` instance.
        """
        return getattr(self, "_pos", None)

    @pos.setter
    def pos(self, pos):
        setattr(self, "_pos", pos)

        if not isinstance(pos, SkyCoord):
            try:
                ra, dec = pos
            except (TypeError, ValueError):
                raise ValueError(
                    'Pos must be a sequence with exactly two values, '
                    'expressing ra and dec in icrs degrees'
                )

            # assume degrees
            pos = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")

        self["POS"] = "{ra},{dec}".format(
            ra=pos.icrs.ra.deg, dec=pos.icrs.dec.deg)

    @pos.deleter
    def pos(self):
        delattr(self, "_pos")
        del self["POS"]

    @property
    def diameter(self):
        """
        the diameter of the circular region around pos as a
        `~astropy.units.Quantity` instance.
        """
        return getattr(self, "_diameter", None)

    @diameter.setter
    def diameter(self, diameter):
        setattr(self, "_diameter", diameter)

        if not isinstance(diameter, Quantity):
            valerr = ValueError(
                'Diameter must be exactly one value, expressing degrees')

            try:
                # assume degrees
                diameter = diameter * Unit("deg")
            except ValueError:
                raise valerr

            try:
                if len(diameter):
                    raise valerr
            except TypeError:
                pass  # len 1

        self["SIZE"] = diameter.to(Unit("deg")).value

    @diameter.deleter
    def diameter(self):
        delattr(self, "_diameter")
        del self["SIZE"]

    @property
    def band(self):
        """
        the bandwidth range the observations belong to.
        """
        return getattr(self, "_band", None)

    @band.setter
    def band(self, band):
        setattr(self, "_band", band)

        if not isinstance(band, Quantity):
            valerr = ValueError(
                'Band must be a sequence with exactly two values',
                'expressing a frequency or wavelength range')

            try:
                # assume meters
                band = band * Unit("meter")
            except ValueError:
                raise valerr

            try:
                if len(band) != 2:
                    raise valerr
            except TypeError:
                raise valerr

        # transform to meters
        band = band.to(Unit("m"), equivalencies=spectral_equivalencies())
        # frequency is counter-proportional to wavelength, so we just sort
        # it to have the right order again
        band.sort()

        self["BAND"] = "{start}/{end}".format(
            start=band.value[0], end=band.value[1])

    @band.deleter
    def band(self):
        delattr(self, "_band")
        del self["BAND"]

    @property
    def time(self):
        """
        the datetime range the observations were made in.
        """
        return getattr(self, "_time", None)

    @time.setter
    def time(self, time):
        setattr(self, "_time", time)

        if not isinstance(time, Time):
            valerr = ValueError(
                'Time must be a sequence with exactly two values, '
                'expressing a datetime in ISO 8601'
            )

            try:
                # assume iso8601
                time = Time(time, format="isot")
            except ValueError:
                raise valerr

            try:
                if len(time) != 2:
                    raise valerr
            except TypeError:
                raise valerr

        # It seems astropy either has seconds and microseconds (the date_hms
        # subformat) or no seconds at all (the date_hm subformat).  SSAP
        # probably doesn't allow microseconds.  Rather than fix this
        # via a new astropy subformat, let's get by with local string
        # operations.
        literals = time.to_value('isot')
        self["TIME"] = "{start}/{end}".format(
            start=literals[0].split(".")[0],
            end=literals[1].split(".")[0])

    @time.deleter
    def time(self):
        delattr(self, "_time")
        del self["TIME"]

    @property
    def format(self):
        """
        the image format(s) of interest.  "all" indicates
        all available formats; "graphic" indicates
        graphical images (e.g. jpeg, png, gif; not FITS);
        "metadata" indicates that no images should be
        returned--only an empty table with complete metadata.
        """
        return getattr(self, "_format", None)

    @format.setter
    def format(self, val):
        setattr(self, "_format", val)

        if isinstance(val, (str, bytes)):
            val = [val]

        self["FORMAT"] = ",".join(val)

    @format.deleter
    def format(self):
        delattr(self, "_format")
        del self["FORMAT"]

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
        submit the query and return the results as a SSAResults instance

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
        return SSAResults(self.execute_votable(), url=self.queryurl, session=self._session)


class SSAResults(DatalinkResultsMixin, DALResults):
    """
    The list of matching images resulting from a spectrum (SSA) query.
    Each record contains a set of metadata that describes an available
    spectrum matching the query constraints.  The number of records in
    the results is by passing it to the Python built-in ``len()`` function.

    This class supports iterable semantics; thus,
    individual records (in the form of
    :py:class:`~pyvo.dal.ssa.SSARecord` instances) are typically
    accessed by iterating over an ``SSAResults`` instance.

    Alternatively, records can be accessed randomly via
    :py:meth:`getrecord` or through a Python Database API (v2)
    Cursor (via :py:meth:`~pyvo.dal.DALResults.cursor`).
    Column-based data access is possible via the
    :py:meth:`~pyvo.dal.DALResults.getcolumn` method.

    ``SSAResults`` is essentially a wrapper around an Astropy
    :py:mod:`~astropy.io.votable`
    :py:class:`~astropy.io.votable.tree.TableElement` instance where the
    columns contain the various metadata describing the spectra.
    One can access that VOTable directly via the
    :py:attr:`~pyvo.dal.DALResults.votable` attribute.  Thus,
    when one retrieves a whole column via
    :py:meth:`~pyvo.dal.DALResults.getcolumn`, the result is
    a Numpy array.  Alternatively, one can manipulate the results
    as an Astropy :py:class:`~astropy.table.table.Table` via the
    following conversion:

    ``table = results.votable.to_table()``

    ``SSAResults`` supports the array item operator ``[...]`` in a
    read-only context.  When the argument is numerical, the result
    is an
    :py:class:`~pyvo.dal.ssa.SSARecord` instance, representing the
    record at the position given by the numerical index.  If the
    argument is a string, it is interpreted as the name of a column,
    and the data from the column matching that name is returned as
    a Numpy array.
    """

    def getrecord(self, index):
        """
        return a representation of a sia result record that follows
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
        SIARecord
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
        return SSARecord(self, index, session=self._session)


class SSARecord(SodaRecordMixin, DatalinkRecordMixin, Record):
    """
    a dictionary-like container for data in a record from the results of an
    SSA query, describing an available spectrum.

    The commonly accessed metadata which are stadardized by the SSA
    protocol are available as attributes.  If the metadatum accessible
    via an attribute is not available, the value of that attribute
    will be None.  All metadata, including non-standard metadata, are
    acessible via the ``get(`` *key* ``)`` function (or the [*key*]
    operator) where *key* is table column name.
    """

    @property
    def ra(self):
        """
        return the right ascension of the center of the spectrum
        """
        return self.getbyutype("ssa:Target.Pos")[0]

    @property
    def dec(self):
        """
        return the declination of the center of the spectrum
        """
        return self.getbyutype("ssa:Target.Pos")[1]

    @property
    def title(self):
        """
        return the title of the spectrum
        """
        return self.getbyutype("ssa:DataID.Title", decode=True)

    @property
    def format(self):
        """
        return the file format that this the spectrum is stored in
        """
        return self.getbyutype("ssa:Access.Format", decode=True)

    @property
    def dateobs(self):
        """
        return the modified Julien date (MJD) of the mid-point of the
        observational data that went into the spectrum
        """
        dateobs = self.getbyutype("ssa:DataID.Date", decode=True)
        try:
            if not dateobs or np.isnan(dateobs):
                return None
        except TypeError:
            # np.isnan can only check floats. If can't check for nan, pass it along
            pass
        return Time(dateobs, format="iso")

    @property
    def instr(self):
        """
        return the name of the instrument (or instruments) that produced the
        data that went into this spectrum.
        """
        return self.getbyutype("ssa:DataID.Instrument", decode=True)

    @property
    def acref(self):
        """
        return the URL that can be used to retrieve the spectrum.
        """
        return self.getbyutype("ssa:Access.Reference", decode=True)

    @property
    def filesize(self):
        """
        The (estimated) size of the image in bytes
        """
        return self.getbyutype("ssa:Access.Size")

    def getdataurl(self):
        """
        return the URL contained in the access URL column which can be used
        to retrieve the dataset described by this record.  None is returned
        if no such column exists.
        """
        dataurl = super().getdataurl()
        if dataurl is None:
            return self.acref
        else:
            return dataurl

    def suggest_dataset_basename(self):
        """
        return a default base filename that the dataset available via
        ``getdataset()`` can be saved as.  This function is
        specialized for a particular service type this record originates from
        so that it can be used by ``cachedataset()`` via
        ``make_dataset_filename()``.
        """
        out = self.title
        if isinstance(out, bytes):
            out = out.decode('utf-8')

        if not out:
            out = "spectrum"
        else:
            out = re.sub(r'\s+', '_', out.strip())
        return out

    def suggest_extension(self, *, default=None):
        """
        returns a recommended filename extension for the dataset described
        by this record.  Typically, this would look at the column describing
        the format and choose an extension accordingly.
        """
        return mime2extension(self.format, default)

    def broadcast_samp(self, *, client_name=None):
        """
        Broadcast the spectrum to ``client_name`` via SAMP
        """
        with samp.connection() as conn:
            samp.send_spectrum_to(
                conn, self.getdataurl(), client_name,
                name=self.suggest_dataset_basename())
