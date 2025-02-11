# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module for searching for images in a remote archive.

A Simple Image Access (SIA) service allows a client to search for
images in an archive whose field of view overlaps with a given
rectangular region on the sky.  The service responds to a search query
with a table in which each row represents an image that is available
for download.  The columns provide metadata describing each image and
one column in particular provides the image's download URL (also
called the *access reference*, or *acref*).  Some SIA services act as
a cut-out service; in this case, the query result is a table of images
whose field of view matches the requested region and which will be
created when accessed via the download URL.

This module provides an interface for accessing an SIA service.  It is
implemented as a specialization of the DAL Query interface.

The ``search()`` function support the simplest and most common types
of queries, returning an SIAResults instance as its results which
represents the matching images from the archive.  The SIAResults
supports access to and iterations over the individual records; these
are provided as SIARecord instances, which give easy access to key
metadata in the response, such as the position of the image's center,
the image format, the size and shape of the image, and its download
URL.

The ``SIAService`` class can represent a specific service available at a URL
endpoint.
"""
import re

from pyvo.io.vosi.vodataservice import TableParam

from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.units import Quantity, Unit
import numpy as np

from .query import DALResults, DALQuery, DALService, Record
from .mimetype import mime2extension
from .adhoc import DatalinkResultsMixin, DatalinkRecordMixin, SodaRecordMixin

from .. import samp

__all__ = ["search", "SIAService", "SIAQuery", "SIAResults", "SIARecord"]


def search(
        url, pos, size=1.0, *, format=None, intersect=None, verbosity=2,
        **keywords):
    """
    submit a simple SIA query that requests images overlapping a given region

    Parameters
    ----------
    url : str
        the base URL for the SIA service
    pos : `~astropy.coordinates.SkyCoord` class or sequence of two floats
        the position of the center of the rectangular search region.
        assuming icrs decimal degrees if unit is not specified.
    size : `~astropy.units.Quantity` class or up to 2 floats.
        the full rectangular size of the search region along the
        RA and Dec directions.
        converted if it's a iterable containing scalars,
        assuming decimal degrees.
    format : str
        the image format(s) of interest.  "all" (server-side default)
        indicates all available formats; "graphic" indicates
        graphical images (e.g. jpeg, png, gif; not FITS);
        "metadata" indicates that no images should be
        returned--only an empty table with complete metadata;
        "image/\\*" indicates a particular image format where
        * can have values like "fits", "jpeg", "png", etc.
    intersect : str
        a token indicating how the returned images should
        intersect with the search region; recognized values include:

        ========= ======================================================
        COVERS    select images that completely cover the search region
        ENCLOSED  select images that are complete enclosed by the region
        OVERLAPS  select any image that overlaps the search region (server-side default)
        CENTER    select images whose center is within the search region
        ========= ======================================================

    verbosity : int
        an integer value that indicates the volume of columns
        to return in the result table.  0 means the minimum
        set of columsn, 3 means as many columns as are  available.
        (client-side default == 2)
    **keywords :
        additional parameters can be given via arbitrary
        case insensitive keyword arguments. Where there is overlap
        with the parameters set by the other arguments to
        this function, these keywords will override.

    Returns
    -------
    SIAResults
        a container holding a table of matching image records

    Raises
    ------
    DALServiceError
        for errors connecting to or communicating with the service
    DALQueryError
        if the service responds with an error,
        including a query syntax error.

    See Also
    --------
    SIAResults
    pyvo.dal.DALServiceError
    pyvo.dal.DALQueryError
    """
    service = SIAService(url)
    return service.search(pos=pos, size=size, format=format, intersect=intersect,
                          verbosity=verbosity, **keywords)


class SIAService(DALService):
    """
    a representation of an SIA service
    """

    def __init__(self, baseurl, *, capability_description=None, session=None):
        """
        instantiate an SIA service

        Parameters
        ----------
        baseurl : str
           the base URL for submitting search queries to the service.
        session : object
           optional session to use for network requests
        """
        super().__init__(baseurl, capability_description=capability_description, session=session)
        self._description = capability_description

    def _get_metadata(self):
        """
        the metadata resource element
        """
        if not hasattr(self, "_metadata"):
            query = self.create_query(format='metadata')
            metadata = query.execute_votable()

            setattr(self, "_metadata", metadata)
            try:
                setattr(self, "_metadata_resource", metadata.resources[0])
            except IndexError:
                setattr(self, "_metadata_resource", None)

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
    def params(self):
        """
        the service parameters.
        """
        self._get_metadata()

        try:
            return getattr(self, "_metadata_resource", None).params
        except AttributeError:
            return None

    @property
    def columns(self):
        """
        the available columns on this service
        """
        self._get_metadata()
        fields = getattr(self, '_metadata', None).get_first_table().fields

        try:
            return [
                TableParam.from_field(field) for field in fields]
        except AttributeError:
            return []

    def search(
            self, pos, size=1.0, *, format=None, intersect=None,
            verbosity=2, **keywords):
        """
        submit a SIA query to this service with the given parameters.

        Parameters
        ----------
        pos : `~astropy.coordinates.SkyCoord` class or sequence of two floats
            the position of the center of the rectangular search region.
            assuming icrs decimal degrees if unit is not specified.
        size : `~astropy.units.Quantity` class or up to 2 floats.
            the full rectangular size of the search region along the
            RA and Dec directions.
            converted if it's a iterable containing scalars,
            assuming decimal degrees.
        size : `~astropy.units.Quantity` class or scalar float
            the size of the rectangular region around pos.
            assuming icrs decimal degrees if unit is not specified.
        format : str
            the image format(s) of interest.  "all" (server-side default)
            indicates all available formats; "graphic" indicates
            graphical images (e.g. jpeg, png, gif; not FITS);
            "metadata" indicates that no images should be
            returned--only an empty table with complete metadata;
            "image/\\*" indicates a particular image format where
            * can have values like "fits", "jpeg", "png", etc.
        intersect : str
            a token indicating how the returned images should
            intersect with the search region; recognized values include:

            ========= ======================================================
            COVERS    select images that completely cover the search region
            ENCLOSED  select images that are complete enclosed by the region
            OVERLAPS  select any image that overlaps the search region (server-side default)
            CENTER    select images whose center is within the search region
            ========= ======================================================

        verbosity : int
            an integer value that indicates the volume of columns
            to return in the result table.  0 means the minimum
            set of columns, 3 means as many columns as are  available.
            (client-side default == 2)
        **keywords :
            additional parameters can be given via arbitrary
            case insensitive keyword arguments. Where there is overlap
            with the parameters set by the other arguments to
            this function, these keywords will override.

        Returns
        -------
        SIAResults
           a container holding a table of matching catalog records

        Raises
        ------
        DALServiceError
            for errors connecting to or communicating with the service
        DALQueryError
            if the service responds with an error,
            including a query syntax error.

        See Also
        --------
        SIAResults
        pyvo.dal.DALServiceError
        pyvo.dal.DALQueryError
        """
        return self.create_query(
            pos=pos, size=size, format=format, intersect=intersect, verbosity=verbosity, **keywords).execute()

    def create_query(
            self, pos=None, size=None, *, format=None, intersect=None,
            verbosity=None, **keywords):
        """
        create a query object that constraints can be added to and then
        executed.  The input arguments will initialize the query with the
        given values.

        Parameters
        ----------
        pos : `~astropy.coordinates.SkyCoord` class or sequence of two floats
            the position of the center of the rectangular search region.
            assuming icrs decimal degrees if unit is not specified.
        size : `~astropy.units.Quantity` class or up to 2 floats.
            the full rectangular size of the search region along the
            RA and Dec directions.
            converted if it's a iterable containing scalars,
            assuming decimal degrees.
        size : `~astropy.units.Quantity` class or scalar float
            the size of the rectangular region around pos.
            assuming icrs decimal degrees if unit is not specified.
        format : str
            the image format(s) of interest.  "all" (server-side default)
            indicates all available formats; "graphic" indicates
            graphical images (e.g. jpeg, png, gif; not FITS);
            "metadata" indicates that no images should be
            returned--only an empty table with complete metadata;
            "image/\\*" indicates a particular image format where
            * can have values like "fits", "jpeg", "png", etc.
        intersect : str
            a token indicating how the returned images should
            intersect with the search region; recognized values include:

            ========= ======================================================
            COVERS    select images that completely cover the search region
            ENCLOSED  select images that are complete enclosed by the region
            OVERLAPS  select any image that overlaps the search region (server-side default)
            CENTER    select images whose center is within the search region
            ========= ======================================================

        verbosity : int
            an integer value that indicates the volume of columns
            to return in the result table.  0 means the minimum
            set of columsn, 3 means as many columns as are  available.
        **keywords :
            additional parameters can be given via arbitrary
            case insensitive keyword arguments. Where there is overlap
            with the parameters set by the other arguments to
            this function, these keywords will override.

        Returns
        -------
        SIAQuery
            the query instance

        See Also
        --------
        SIAQuery
        """
        return SIAQuery(
            self.baseurl, pos=pos, size=size, format=format, intersect=intersect,
            verbosity=verbosity, session=self._session, **keywords)

    def describe(self):
        print(self.description)
        print()


class SIAQuery(DALQuery):
    """
    a class for preparing a query to an SIA service.  Query constraints
    are added via its service type-specific methods.  The various execute()
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
            self, baseurl, pos=None, size=None, *, format=None, intersect=None,
            verbosity=None, session=None, **keywords):
        """
        initialize the query object with a baseurl and the given parameters

        Parameters
        ----------
        baseurl : str
            the base URL for the SIA service
        pos : `~astropy.coordinates.SkyCoord` class or sequence of two floats
            the position of the center of the rectangular search region.
            assuming icrs decimal degrees if unit is not specified.
        size : `~astropy.units.Quantity` class or up to 2 floats.
            the full rectangular size of the search region along the
            RA and Dec directions.
            converted if it's a iterable containing scalars,
            assuming decimal degrees.
        size : `~astropy.units.Quantity` class or scalar float
            the size of the rectangular region around pos.
            assuming icrs decimal degrees if unit is not specified.
        format : str
            the image format(s) of interest.  "all" (server-side default)
            indicates all available formats; "graphic" indicates
            graphical images (e.g. jpeg, png, gif; not FITS);
            "metadata" indicates that no images should be
            returned--only an empty table with complete metadata;
            "image/\\*" indicates a particular image format where
            * can have values like "fits", "jpeg", "png", etc.
        intersect : str
            a token indicating how the returned images should
            intersect with the search region; recognized values include:

            ========= ======================================================
            COVERS    select images that completely cover the search region
            ENCLOSED  select images that are complete enclosed by the region
            OVERLAPS  select any image that overlaps the search region (server-side default)
            CENTER    select images whose center is within the search region
            ========= ======================================================

        verbosity : int
            an integer value that indicates the volume of columns
            to return in the result table.  0 means the minimum
            set of columsn, 3 means as many columns as are  available.
        session : object
           optional session to use for network requests
        **keywords :
            additional parameters can be given via arbitrary
            case insensitive keyword arguments. Where there is overlap
            with the parameters set by the other arguments to
            this function, these keywords will override.
        """
        super().__init__(baseurl, session=session, **keywords)

        if pos is not None:
            self.pos = pos

        if size is not None:
            self.size = size

        if format is not None:
            self.format = format

        if intersect is not None:
            self.intersect = intersect

        if verbosity is not None:
            self.verbosity = verbosity

    @property
    def pos(self):
        """
        the position of the center of the rectangular search region as a
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
    def size(self):
        """
        the size of the rectangular region around pos as a
        `~astropy.units.Quantity` instance.
        """
        return getattr(self, "_size", None)

    @size.setter
    def size(self, size):
        setattr(self, "_size", size)

        if not isinstance(size, Quantity):
            valerr = ValueError(
                'Size must be either a single value or a sequence with two'
                'values, expressing degrees'
            )

            try:
                # assume degrees
                size = size * Unit("deg")
            except ValueError:
                raise valerr

            try:
                if len(size) > 2:
                    raise valerr
            except TypeError:
                pass  # len 1

        try:
            self["SIZE"] = ",".join(
                str(deg) for deg in size.to(Unit("deg")).value)
        except TypeError:
            self["SIZE"] = str(size.to(Unit("deg")).value)

    @size.deleter
    def size(self):
        delattr(self, "_size")
        del self["SIZE"]

    @property
    def format(self):
        """
        the image format(s) of interest.  "all" (default)
        indicates all available formats; "graphic" indicates
        graphical images (e.g. jpeg, png, gif; not FITS);
        "metadata" indicates that no images should be
        returned--only an empty table with complete metadata;
        "image/\\*" indicates a particular image format where
        * can have values like "fits", "jpeg", "png", etc.
        """
        return getattr(self, "_format", None)

    @format.setter
    def format(self, format_):
        setattr(self, "_format", format_)

        if not isinstance(format_, list):
            format_ = format_.split(',')
        normalized_formats = []
        for user_input in format_:
            if user_input.upper() in ['ALL', 'METADATA', 'GRAPHIC', 'GRAPHIC-ALL']:
                normalized_formats.append(user_input.upper())
            elif user_input.split('-')[0].upper() == 'GRAPHIC':
                normalized_formats.append(user_input.split('-')[0].upper()+"-"+user_input.split('-')[1])
            else:
                normalized_formats.append(user_input)

        self["FORMAT"] = ",".join(normalized_formats)

    @format.deleter
    def format(self):
        delattr(self, "_format")
        del self["FORMAT"]

    @property
    def intersect(self):
        """
        a token indicating how the returned images should
        intersect with the search region; recognized values include:

        ========= ======================================================
        COVERS    select images that completely cover the search region
        ENCLOSED  select images that are complete enclosed by the region
        OVERLAPS  select any image that overlaps with the search region
        CENTER    select images whose center is within the search region
        ========= ======================================================
        """
        return getattr(self, "_intersect", None)

    @intersect.setter
    def intersect(self, intersect):
        setattr(self, "_intersect", intersect)
        self["INTERSECT"] = intersect.upper()

    @intersect.deleter
    def intersect(self):
        delattr(self, "_intersect")
        del self["INTERSECT"]

    @property
    def verbosity(self):
        """
        an integer value that indicates the volume of columns
        to return in the result table.  0 means the minimum
        set of columsn, 3 means as many columns as are  available.
        """
        return getattr(self, "_verbosity", None)

    @verbosity.setter
    def verbosity(self, verbosity):
        setattr(self, "_verbosity", verbosity)
        self["VERB"] = verbosity

    @verbosity.deleter
    def verbosity(self):
        delattr(self, "_verbosity")
        del self["VERB"]

    def execute(self):
        """
        submit the query and return the results as a SIAResults instance

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
        return SIAResults(self.execute_votable(), url=self.queryurl, session=self._session)


class SIAResults(DatalinkResultsMixin, DALResults):
    """
    The list of matching images resulting from an image (SIA) query.
    Each record contains a set of metadata that describes an available
    image matching the query constraints.  The number of records in
    the results is available by passing it to the Python built-in ``len()`` function.

    This class supports iterable semantics; thus,
    individual records (in the form of
    :py:class:`~pyvo.dal.sia.SIARecord` instances) are typically
    accessed by iterating over an ``SIAResults`` instance.

    Alternatively, records can be accessed randomly via
    :py:meth:`getrecord` or through a Python Database API (v2)
    Cursor (via :py:meth:`~pyvo.dal.DALResults.cursor`).
    Column-based data access is possible via the
    :py:meth:`~pyvo.dal.DALResults.getcolumn` method.

    ``SIAResults`` is essentially a wrapper around an Astropy
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

    ``SIAResults`` supports the array item operator ``[...]`` in a
    read-only context.  When the argument is numerical, the result
    is an
    :py:class:`~pyvo.dal.sia.SIARecord` instance, representing the
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
        return SIARecord(self, index, session=self._session)


class SIARecord(SodaRecordMixin, DatalinkRecordMixin, Record):
    """
    a dictionary-like container for data in a record from the results of an
    image (SIA) search, describing an available image.

    The commonly accessed metadata which are stadardized by the SIA
    protocol are available as attributes.  If the metadatum accessible
    via an attribute is not available, the value of that attribute
    will be None.  All metadata, including non-standard metadata, are
    acessible via the ``get(`` *key* ``)`` function (or the [*key*]
    operator) where *key* is table column name.
    """

    def getdataformat(self):
        """
        return the mimetype of the dataset described by this record.
        """
        return self.getbyucd("VOX:Image_Format", decode=True)

    @property
    def pos(self):
        """
        the position of the object or observation described by this record.
        """
        return SkyCoord(
            ra=self.getbyucd("POS_EQ_RA_MAIN"),
            dec=self.getbyucd("POS_EQ_DEC_MAIN"),
            unit="deg", frame="icrs")

    # Image Metadata
    @property
    def title(self):
        """
        the title of the image
        """
        return self.getbyucd("VOX:Image_Title", decode=True)

    @property
    def instr(self):
        """
        the name of the instrument (or instruments) that produced the data that
        went into this image.
        """
        return self.getbyucd("INST_ID", decode=True)

    @property
    def dateobs(self):
        """
        the modified Julian date (MJD) of the mid-point of the
        observational data that went into the image,
        as an astropy.time.Time instance
        """
        dateobs = self.getbyucd("VOX:Image_MJDateObs")
        try:
            if not dateobs or np.isnan(dateobs):
                return None
        except TypeError:
            # np.isnan can only check floats. If can't check for nan, pass it along
            pass
        return Time(dateobs, format="mjd")

    @property
    def naxes(self):
        """
        the number of axes in this image.
        """
        return self.getbyucd("VOX:Image_Naxes")

    @property
    def naxis(self):
        """
        the lengths of the sides along each axis, in pix,
        as a astropy Quantity pix
        """
        return self.getbyucd("VOX:Image_Naxis") * Unit("pix")

    @property
    def scale(self):
        """
        the scale of the pixels in each image axis, in degrees/pixel,
        as a astropy Quantity deg / pix
        """
        return self.getbyucd("VOX:Image_Scale") * (Unit("deg") / Unit("pix"))

    @property
    def format(self):
        """
        the format of the image
        """
        return self.getbyucd("VOX:Image_Format", decode=True)

    # Coordinate System Metadata
    @property
    def coord_frame(self):
        """
        the coordinate system reference frame, one of the following:
        "ICRS", "FK5", "FK4", "ECL", "GAL", and "SGAL".
        """
        return self.getbyucd("VOX:STC_CoordRefFrame", decode=True)

    @property
    def coord_equinox(self):
        """
        the equinox of the used coordinate system
        """
        return self.getbyucd("VOX:STC_CoordEquinox")

    @property
    def coord_projection(self):
        """
        the celestial projection (TAN / ARC / SIN / etc.)
        """
        return self.getbyucd("VOX:WCS_CoordProjection", decode=True)

    @property
    def coord_refpixel(self):
        """
        the image pixel coordinates of the WCS reference pixel
        """
        return self.getbyucd("VOX:WCS_CoordRefPixel")

    @property
    def coord_refvalue(self):
        """
        the world coordinates of the WCS reference pixel.
        """
        return self.getbyucd("VOX:WCS_CoordRefValue")

    @property
    def cdmatrix(self):
        """
        the WCS CD matrix defining the scale and rotation (among other things)
        of the image. ordered as CD[i,j] = [0,0], [0,1], [1,0], [1,1].
        """
        return self.getbyucd("VOX:WCS_CDMatrix").reshape((2, 2))

    # Spectral Bandpass Metadata
    @property
    def bandpass_id(self):
        """
        the bandpass by name (e.g., "V", "SDSS_U", "K", "K-Band", etc.)
        """
        return self.getbyucd("VOX:BandPass_ID", decode=True)

    @property
    def bandpass_unit(self):
        """
        the astropy unit used to represent spectral values.
        """
        sia_unit = self.getbyucd("VOX:BandPass_Unit", decode=True)

        # the standard says this should be `"meters", "hertz", and "keV"',
        # which in practice everyone has ignored.  Still, let's
        # translate the standard terms.  Somewhat more dangerously,
        # we replace a missing unit with metres; in theory, we should
        # reject anything but the three terms above, but then 99%
        # of SIA services would break.  And without the assumption of
        # a metre default, 20% of 2024 SIA1 services would have unusable
        # bandpasses.

        astropy_unit = {
            None: "m",
            "": "m",
            "meters": "m",
            "hertz": "Hz"}.get(sia_unit, sia_unit)

        return Unit(astropy_unit)

    @property
    def bandpass_refvalue(self):
        """
        the characteristic (reference) wavelength, frequency or energy
        for the bandpass model, as an astropy Quantity of bandpass_unit
        """
        value = self.getbyucd("VOX:BandPass_RefValue")
        if value is not None and self.bandpass_unit:
            return Quantity(value, self.bandpass_unit)

    @property
    def bandpass_hilimit(self):
        """
        the upper limit of the bandpass as astropy Quantity in bandpass_unit
        """
        value = self.getbyucd("VOX:BandPass_HiLimit")
        if value is not None and self.bandpass_unit:
            return Quantity(value, self.bandpass_unit)

    @property
    def bandpass_lolimit(self):
        """
        the lower limit of the bandpass as astropy Quantity in bandpass_unit
        """
        value = self.getbyucd("VOX:BandPass_LoLimit")
        if value is not None and self.bandpass_unit:
            return Quantity(value, self.bandpass_unit)

    # Processig Metadata
    @property
    def pixflags(self):
        """
        the type of processing done by the image service to produce an output
        image pixel

        a string of one or more of the following values:

        * C -- The image pixels were copied from a source image without change,
               as when an atlas image or cutout is returned.
        * F -- The image pixels were computed by resampling an existing image,
               e.g., to rescale or reproject the data,
               and were filtered by an interpolator.
        * X -- The image pixels were computed by the service directly from a
               primary data set hence were not filtered by an interpolator.
        * Z -- The image pixels contain valid flux (intensity) values, e.g., if
               the pixels were resampled with a flux-preserving interpolator.
        * V -- The image pixels contain some unspecified visualization of the
               data, hence are suitable for display but not for numerical
               analysis.
        """
        return self.getbyucd("VOX:Image_PixFlags", decode=True)

    # Access Metadata
    @property
    def acref(self):
        """
        the URL that can be used to retrieve the image
        """
        return self.getbyucd("VOX:Image_AccessReference", decode=True)

    @property
    def acref_ttl(self):
        """
        the minimum time to live in seconds of the access reference
        """
        return self.getbyucd("VOX:Image_AccessRefTTL")

    @property
    def filesize(self):
        """
        the (estimated) size of the image in bytes
        """
        return self.getbyucd("VOX:Image_FileSize")

    def getdataurl(self):
        """
        return the URL contained in the access URL column which can be used
        to retrieve the dataset described by this record.  None is returned
        if no such column exists.
        """
        if self.acref is not None:
            return self.acref
        else:
            return super().getdataurl()

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
            out = "image"
        else:
            out = re.sub(r'\s+', '_', out.strip())
        return out

    def suggest_extension(self, *, default='dat'):
        """
        returns a recommended filename extension for the dataset described
        by this record.  Typically, this would look at the column describing
        the format and choose an extension accordingly.
        """
        return mime2extension(self.format, default)

    def broadcast_samp(self, *, client_name=None):
        """
        Broadcast the image to ``client_name`` via SAMP
        """
        with samp.connection() as conn:
            samp.send_image_to(
                conn, self.getdataurl(), client_name,
                name=self.suggest_dataset_basename())
