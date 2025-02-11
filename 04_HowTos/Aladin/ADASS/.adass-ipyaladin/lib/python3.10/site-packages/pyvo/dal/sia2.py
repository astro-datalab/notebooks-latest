# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module for searching for images in a remote archive.

A Simple Image Access v2 (SIA2) service allows a client to search for
images based on a number of criteria/parameters. The results are
represented in ``pyvo.dam.ObsCoreMetadata`` format.

The ``SIA2Service`` class can represent a specific service available at a URL
endpoint.
"""

import warnings
import numpy as np

from astropy import units as u
from astropy.time import Time
from astropy.utils.decorators import deprecated
from astropy.utils.exceptions import AstropyDeprecationWarning

from .query import DALResults, DALQuery, DALService, Record
from .adhoc import DatalinkResultsMixin, AxisParamMixin, SodaRecordMixin, DatalinkRecordMixin
from .params import IntervalQueryParam, StrQueryParam, EnumQueryParam
from .vosi import AvailabilityMixin, CapabilityMixin
from ..dam import ObsCoreMetadata, CALIBRATION_LEVELS


__all__ = ["search", "SIA2Service", "SIA2Query", "SIA2Results", "ObsCoreRecord"]

SIA2_STANDARD_ID = 'ivo://ivoa.net/std/SIA#query-2.0'


SIA2_PARAMETERS_DESC = """
pos : single or list of tuples
    angle units (default: deg)
    the positional region(s) to be searched for data. Each region can
    be expressed as a tuple representing a CIRCLE, RANGE or POLYGON as
    follows:
    (ra, dec, radius) - for CIRCLE. (angle units - defaults to)
    (long1, long2, lat1, lat2) - for RANGE (angle units required)
    (ra, dec, ra, dec, ra, dec ... ) ra/dec points for POLYGON all
    in angle units
band : scalar, tuple(interval) or list of tuples
    (spectral units (default: meter)
    the energy interval(s) to be searched for data.
time : single or list of `~astropy.time.Time` or compatible strings
    the time interval(s) to be searched for data.
pol : single or list of str from ``pyvo.dam.obscore.POLARIZATION_STATES``
    the polarization state(s) to be searched for data.
field_of_view : single or list of tuples
    angle units (default arcsec)
    the range(s) of field of view (size) to be searched for data
spatial_resolution : single or list of tuples
    angle units required
    the range(s) of spatial resolution to be searched for data
spectral_resolving_power : single or list of tuples
    the range(s) of spectral resolving power to be searched for data
exptime : single or list of tuples
    time units (default: second)
    the range(s) of exposure times to be searched for data
timeres : single of list of tuples
    time units (default: second)
    the range(s) of temporal resolution to be searched for data
publisher_did : single or list of str
    specifies the unique identifier of dataset(s). It is global because
    it must include information regarding the publisher
    (obs_publisher_did in ObsCore)
collection : single or list of str
    name of the collection that the data belongs to
facility : single or list of str
    specifies the name of the facility (usually telescope) where
    the data was acquired.
instrument : single or list of str
    specifies the name of the instrument with which the data was
    acquired.
data_type : 'image'|'cube'
    specifies the type of the data
calib_level : single or list from enum
    ``pyvo.dam.obscore.CALIBRATION_LEVELS``
    specifies the calibration level of the data. Can be a single value
    or a list of values
target_name : single or list of str
    specifies the name of the target (e.g. the intention of the
    original science program or observation)
res_format : single or list of strings
    specifies response format(s).
max_records : int
    allows the client to limit the number or records in the response
**kwargs : custom query parameters
    single or a list of values (or tuples for intervals) custom query
    parameters that a specific service accepts. The values of the
    parameters need to follow the SIA2 format and represent the
    appropriate quantities (where applicable).
"""


def __getattr__(name):
    if name == 'SIA_PARAMETERS_DESC':
        warnings.warn("The name SIA_PARAMETERS_DESC is deprecated in v1.5 for SIA v2 services, "
                      "use SIA2_PARAMETERS_DESC instead.", AstropyDeprecationWarning)
        return SIA2_PARAMETERS_DESC
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def search(url, pos=None, *, band=None, time=None, pol=None,
           field_of_view=None, spatial_resolution=None,
           spectral_resolving_power=None, exptime=None,
           timeres=None, publisher_did=None, facility=None, collection=None,
           instrument=None, data_type=None, calib_level=None,
           target_name=None, res_format=None, maxrec=None, session=None,
           **kwargs):
    """
    submit a simple SIA query to a SIA2 compatible service

    Parameters
    ----------

    url : str
       url of the SIA service (base or endpoint)
    _SIA2_PARAMETERS

    """
    service = SIA2Service(url, session=session)
    return service.search(pos=pos, band=band, time=time, pol=pol,
                          field_of_view=field_of_view,
                          spatial_resolution=spatial_resolution,
                          spectral_resolving_power=spectral_resolving_power,
                          exptime=exptime, timeres=timeres,
                          publisher_did=publisher_did,
                          facility=facility, collection=collection,
                          instrument=instrument, data_type=data_type,
                          calib_level=calib_level, target_name=target_name,
                          res_format=res_format, maxrec=maxrec, **kwargs)


search.__doc__ = search.__doc__.replace('_SIA2_PARAMETERS',
                                        SIA2_PARAMETERS_DESC)


def _tolist(value):
    # return value as a list - is there something in Python to do that?
    if not value:
        return []
    if isinstance(value, list):
        return value
    return [value]


class SIA2Service(DALService, AvailabilityMixin, CapabilityMixin):
    """
    a representation of an SIA2 service

    Note that within pyVO, SIA2 services are associated with the
    (non-existing) standard id ivo://ivoa.net/std/sia2 rather than
    ivo://ivoa.net/std/sia#query-2.0 as in the standard; users should
    generally not notice that, though.
    """

    def __init__(self, baseurl, *, capability_description=None, session=None, check_baseurl=True):
        """
        instantiate an SIA2 service

        Parameters
        ----------
        url : str
           url - URL of the SIA2service (base or query endpoint)
        session : object
           optional session to use for network requests
        check_baseurl : bool
           True - use the capabilities end point of the service to get the
           query end point, False - baseurl is the query end point
        """

        super().__init__(baseurl, capability_description=capability_description, session=session)

        # Check if the session has an update_from_capabilities attribute.
        # This means that the session is aware of IVOA capabilities,
        # and can use this information in processing network requests.
        # One such usecase for this is auth.
        if hasattr(self._session, 'update_from_capabilities'):
            self._session.update_from_capabilities(self.capabilities)

        self.query_ep = baseurl.strip('&')  # default service end point
        if check_baseurl:
            for cap in self.capabilities:
                # assumes that the access URL is the same regardless of the
                # authentication method except BasicAA which is not supported
                # in pyvo. So pick any access url as long as it's not
                if cap.standardid.lower() == SIA2_STANDARD_ID.lower():
                    for interface in cap.interfaces:
                        if interface.accessurls and \
                                not (len(interface.securitymethods) == 1
                                     and interface.securitymethods[0].standardid
                                     == 'ivo://ivoa.net/sso#BasicAA'):
                            self.query_ep = interface.accessurls[0].content
                            break
                    else:
                        continue
                    break

    def search(self, pos=None, *, band=None, time=None, pol=None,
               field_of_view=None, spatial_resolution=None,
               spectral_resolving_power=None, exptime=None,
               timeres=None, publisher_did=None, facility=None, collection=None,
               instrument=None, data_type=None, calib_level=None,
               target_name=None, res_format=None, maxrec=None, **kwargs):
        """
        Performs a SIA2 search against a SIA2 service

        See Also
        --------
        pyvo.dal.sia2.SIA2Query

        """
        return SIA2Query(self.query_ep, pos=pos, band=band,
                         time=time, pol=pol,
                         field_of_view=field_of_view,
                         spatial_resolution=spatial_resolution,
                         spectral_resolving_power=spectral_resolving_power,
                         exptime=exptime, timeres=timeres,
                         publisher_did=publisher_did,
                         facility=facility, collection=collection,
                         instrument=instrument, data_type=data_type,
                         calib_level=calib_level, target_name=target_name,
                         res_format=res_format, maxrec=maxrec,
                         session=self._session, **kwargs).execute()


class SIA2Query(DALQuery, AxisParamMixin):
    """
    a class very similar to :py:attr:`~pyvo.dal.SIAQuery` class but
    used to interact with SIA2 services.
    """

    def __init__(self, url, pos=None, *, band=None, time=None, pol=None,
                 field_of_view=None, spatial_resolution=None,
                 spectral_resolving_power=None, exptime=None,
                 timeres=None, publisher_did=None,
                 facility=None, collection=None,
                 instrument=None, data_type=None, calib_level=None,
                 target_name=None, res_format=None, maxrec=None,
                 session=None, **kwargs):
        """
        initialize the query object with a url and the given parameters

        Note: The majority of the attributes represent constraints used to
        query the SIA2 service and are represented through lists. Multiple value
        attributes are OR-ed in the query, however the values of different
        attributes are AND-ed. Intervals are represented with tuples and
        open-ended intervals should be expressed with float("-inf") or
        float("inf"). Eg. For all values less than or equal to 600 use
        (float(-inf), 600)

        Additional attribute constraints can be specified (or removed) after
        this object has been created using the *.add and *_del methods.

        Parameters
        ----------
        url : url where to send the query request to
        _SIA2_PARAMETERS
        session : object
           optional session to use for network requests

        Returns
        -------
        SIA2Results
            a container holding a table of matching image records. Records are
            represented in IVOA ObsCore format

        Raises
        ------
        DALServiceError
            for errors connecting to or communicating with the service
        DALQueryError
            if the service responds with an error,
            including a query syntax error.

        See Also
        --------
        SIA2Results
        pyvo.dal.DALServiceError
        pyvo.dal.DALQueryError

        """
        super().__init__(url, session=session)

        for pp in _tolist(pos):
            self.pos.add(pp)

        for bb in _tolist(band):
            self.band.add(bb)

        for tt in _tolist(time):
            self.time.add(tt)

        for pp in _tolist(pol):
            self.pol.add(pp)

        for ff in _tolist(field_of_view):
            self.field_of_view.add(ff)

        for sp in _tolist(spatial_resolution):
            self.spatial_resolution.add(sp)

        for sr in _tolist(spectral_resolving_power):
            self.spectral_resolving_power.add(sr)

        for et in _tolist(exptime):
            self.exptime.add(et)

        for tr in _tolist(timeres):
            self.timeres.add(tr)

        for ii in _tolist(publisher_did):
            self.publisher_did.add(ii)

        for ff in _tolist(facility):
            self.facility.add(ff)

        for col in _tolist(collection):
            self.collection.add(col)

        for inst in _tolist(instrument):
            self.instrument.add(inst)

        for dt in _tolist(data_type):
            self.data_type.add(dt)

        for cal in _tolist(calib_level):
            self.calib_level.add(cal)

        for tt in _tolist(target_name):
            self.target_name.add(tt)

        for rf in _tolist(res_format):
            self.res_format.add(rf)

        for name, value in kwargs.items():
            custom_arg = []
            for kw in _tolist(value):
                if isinstance(kw, tuple):
                    val = '{} {}'.format(kw[0], kw[1])
                else:
                    val = str(kw)
                custom_arg.append(val)
            self[name] = custom_arg

        self.maxrec = maxrec

    __init__.__doc__ = \
        __init__.__doc__.replace('_SIA2_PARAMETERS', SIA2_PARAMETERS_DESC)

    @property
    def field_of_view(self):
        if not hasattr(self, '_fov'):
            self._fov = IntervalQueryParam(u.deg)
            self['FOV'] = self._fov.dal
        return self._fov

    @property
    def spatial_resolution(self):
        if not hasattr(self, '_spatres'):
            self._spatres = IntervalQueryParam(u.arcsec)
            self['SPATRES'] = self._spatres.dal
        return self._spatres

    @property
    def spectral_resolving_power(self):
        if not hasattr(self, '_specrp'):
            self._specrp = IntervalQueryParam()
            self['SPECRP'] = self._specrp.dal
        return self._specrp

    @property
    def exptime(self):
        if not hasattr(self, '_exptime'):
            self._exptime = IntervalQueryParam(u.second)
            self['EXPTIME'] = self._exptime.dal
        return self._exptime

    @property
    def timeres(self):
        if not hasattr(self, '_timeres'):
            self._timeres = IntervalQueryParam(u.second)
            self['TIMERES'] = self._timeres.dal
        return self._timeres

    @property
    def publisher_did(self):
        if not hasattr(self, '_publisher_did'):
            self._publisher_did = StrQueryParam()
            self['ID'] = self._publisher_did.dal
        return self._publisher_did

    @property
    def facility(self):
        if not hasattr(self, '_facility'):
            self._facility = StrQueryParam()
            self['FACILITY'] = self._facility.dal
        return self._facility

    @property
    def collection(self):
        if not hasattr(self, '_collection'):
            self._collection = StrQueryParam()
            self['COLLECTION'] = self._collection.dal
        return self._collection

    @property
    def instrument(self):
        if not hasattr(self, '_instrument'):
            self._instrument = StrQueryParam()
            self['INSTRUMENT'] = self._instrument.dal
        return self._instrument

    @property
    def data_type(self):
        if not hasattr(self, '_data_type'):
            self._data_type = StrQueryParam()
            self['DPTYPE'] = self._data_type.dal
        return self._data_type

    @property
    def calib_level(self):
        if not hasattr(self, '_cal'):
            self._cal = EnumQueryParam(CALIBRATION_LEVELS)
            self['CALIB'] = self._cal.dal
        return self._cal

    @property
    def target_name(self):
        if not hasattr(self, '_target'):
            self._target_name = StrQueryParam()
            self['TARGET'] = self._target_name.dal
        return self._target_name

    @property
    def res_format(self):
        if not hasattr(self, '_res_format'):
            self._res_format = StrQueryParam()
            self['FORMAT'] = self._res_format.dal
        return self._res_format

    @property
    def maxrec(self):
        return self._maxrec

    @maxrec.setter
    def maxrec(self, val):
        if not isinstance(val, int):
            if not val:
                return
            raise ValueError(f'maxrec {val} must be non-negative int')
        self._maxrec = val
        self['MAXREC'] = str(val)

    def execute(self):
        """
        submit the query and return the results as a SIA2Results instance

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
        return SIA2Results(self.execute_votable(), url=self.queryurl, session=self._session)


class SIA2Results(DatalinkResultsMixin, DALResults):
    """
    The list of matching images resulting from an image (SIA2) query.
    Each record contains a set of metadata that describes an available
    image matching the query constraints.  The number of records in
    the results is available by passing it to the Python built-in ``len()`` function.

    This class supports iterable semantics; thus,
    individual records (in the form of
    :py:class:`~pyvo.dal.sia2.ObsCoreRecord` instances) are typically
    accessed by iterating over an ``SIA2Results`` instance.

    Alternatively, records can be accessed randomly via
    :py:meth:`getrecord` or through a Python Database API (v2)
    Cursor (via :py:meth:`~pyvo.dal.DALResults.cursor`).
    Column-based data access is possible via the
    :py:meth:`~pyvo.dal.DALResults.getcolumn` method.

    ``SIA2Results`` is essentially a wrapper around an Astropy
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

    ``SIA2Results`` supports the array item operator ``[...]`` in a
    read-only context.  When the argument is numerical, the result
    is an
    :py:class:`~pyvo.dal.sia2.ObsCoreRecord` instance, representing the
    record at the position given by the numerical index.  If the
    argument is a string, it is interpreted as the name of a column,
    and the data from the column matching that name is returned as
    a Numpy array.
    """

    def getrecord(self, index):
        """
        return a representation of a SIA2 result record that follows
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
        ObsCoreMetadataRecord
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
        return ObsCoreRecord(self, index, session=self._session)


class ObsCoreRecord(SodaRecordMixin, DatalinkRecordMixin, Record,
                    ObsCoreMetadata):
    """
    a dictionary-like container for data in a record from the results of an
    image (SIA2) search, describing an available image in ObsCore format.

    The commonly accessed metadata which are stadardized by the SIA2
    protocol are available as attributes.  If the metadatum accessible
    via an attribute is not available, the value of that attribute
    will be None.  All metadata, including non-standard metadata, are also
    acessible via the ``get(`` *key* ``)`` function (or the [*key*]
    operator) where *key* is table column name.
    """

    #          OBSERVATION INFO
    @property
    def dataproduct_type(self):
        """
        Data product (file content) primary type. This is coded as a string
        that conveys a general idea of the content and organization of a
        dataset.
        """
        return self.get('dataproduct_type', decode=True)

    @property
    def dataproduct_subtype(self):
        """
        Data product more specific type
        """
        return self.get('dataproduct_subtype', decode=True, default=None)

    @property
    def calib_level(self):
        """
        Calibration level of the observation: in {0, 1, 2, 3, 4}
        """
        return self.get('calib_level')

    #          TARGET INFO
    @property
    def target_name(self):
        """
        The target_name attribute contains the name of the target of the
        observation, if any. This is typically the name of an astronomical
        object, but could be the name of a survey field.
        The target name is most useful for output, to identify the target of
        an observation to the user. In queries it is generally better to refer
        to astronomical objects by position, using a name resolver to convert
        the target name into a coordinate (when possible).
        """
        return self.get('target_name', decode=True)

    @property
    def target_class(self):
        """
        This field indicates the type of object that was pointed for this
        observation. It is a string with possible values defined in a special
        vocabulary set to be defined: list of object classes (or types) used
        by the SIMBAD database, NED or defined in another IVOA vocabulary.
        """
        return self.get('target_class', decode=True, default=None)

    #          DATA DESCRIPTION
    @property
    def obs_id(self):
        """
        Collection specific internal ID given by the ObsTAP service
        """
        return self.get('obs_id', decode=True)

    @property
    def obs_title(self):
        """
        Brief description of dataset in free format
        """
        return self.get('obs_title', decode=True, default=None)

    @property
    def obs_collection(self):
        """
        The name of the collection (DataID.Collection) identifies the data
        collection to which the data product belongs. A data collection can be
        any collection of datasets which are alike in some fashion. Typical
        data collections might be all the data from a particular telescope,
        instrument, or survey. The value is either the registered shortname
        for the data collection, the full registered IVOA identifier for the
        collection, or a data provider defined short name for the collection.
        Examples: HST/WFPC2, VLT/FORS2, CHANDRA/ACIS-S, etc.
        """
        return self.get('obs_collection', decode=True)

    @property
    def obs_create_date(self):
        """
        Date when the dataset was created
        """
        cd = self.get('obs_create_date', default=None)
        return cd if not cd else Time(cd)

    @property
    def obs_creator_name(self):
        """
        The name of the institution or entity which created the dataset.
        """
        return self.get('obs_creator_name', decode=True, default=None)

    @property
    def obs_creator_did(self):
        """
        IVOA dataset identifier given by its creator.
        """
        return self.get('obs_creator_did', decode=True, default=None)

    #         CURATION INFORMATION
    @property
    def obs_release_date(self):
        """
        Observation release date
        """
        rt = self.get('obs_release_date', default=None, decode=True)
        return rt if not rt else Time(rt)

    @property
    def obs_publisher_did(self):
        """
        ID for the Dataset assigned by the publisher. Note that data from
        a source (creator_did) can be published by multiple publishers
        and have assigned multiple publisher data IDs.
        """
        return self.get('obs_publisher_did', decode=True)

    @property
    def publisher_id(self):
        """
        IVOA-ID for the Publisher. It will also be globally unique since each
        publisher has a unique registered publisher ID
        """
        return self.get('publisher_id', decode=True, default=None)

    @property
    def bib_reference(self):
        """
        URL or bibcode for documentation. This is a forward link to major
        publications which reference the dataset.
        """
        return self.get('bib_reference', decode=True, default=None)

    @property
    def data_rights(self):
        """
        This parameter allows mentioning the availability of a dataset.
        Possible values are: public, secure, or proprietary.
        """
        return self.get('data_rights', decode=True, default=None)

    #           ACCESS INFORMATION
    @property
    def access_url(self):
        """
        The access_url column contains a URL that can be used to download the
        data product (as a file of some sort). Access URLs are not guaranteed
        to remain valid and unchanged indefinitely. To access a specific data
        product after a period of time (e.g., days or weeks) a query should be
        performed to obtain a fresh access URL.
        """
        return self.get('access_url', decode=True)

    @property
    def access_format(self):
        """
        Content format of the dataset. The value of access_format should be a
        MIME type, either a standard MIME type, an extended MIME type from
        the above table, or a new custom MIME type defined by the data
        provider.
        """
        return self.get('access_format', decode=True)

    @property
    def access_estsize(self):
        """
        The approximate size (in kilobytes) of the file available via the
        access_url. This is used only to gain some idea of the size of a data
        product before downloading it, hence only an approximate value is
        required. Provision of dataset size estimates is important whenever it
        is possible that datasets can be very large.
        """
        return self.get('access_estsize') * 1000 * u.byte

    #           SPATIAL CHARACTERISATION
    @property
    def s_ra(self):
        """
        ICRS Right Ascension of the center of the observation
        """
        return self.get('s_ra') * u.deg

    @property
    def s_dec(self):
        """
        CRS Declination of the center of the observation
        """
        return self.get('s_dec') * u.deg

    @property
    def s_fov(self):
        """
        Approximate size of the covered region as the diameter of a containing
        circle. For most data products the value given should be large enough
        to include the entire area of the observation; coverage within the
        bounded region need not be complete, for example if the specified
        radius encompasses a rotated rectangular region. For observations
        which do not have a well-defined boundary, e.g. radio or
        high energy observations, a characteristic value should be given.
        The radius attribute provides a simple way to characterize and use
        (e.g. for discovery computations) the approximate spatial coverage of a
        data product. The spatial coverage of a data product can be more
        precisely specified using the region attribute.
        """
        return self.get('s_fov') * u.deg

    @property
    def s_region(self):
        """
        Sky region covered by the data product (expressed in ICRS frame).
        It can be used to precisely specify the covered spatial region of a
        data product.
        It is often an exact, or almost exact, representation of the
        illumination region of a given observation defined in a standard way
        by the concept of Support in the Characterisation data model.
        """
        return self.get('s_region', decode=True)

    @property
    def s_resolution(self):
        """
        Spatial resolution of data specifies a reference value chosen by the
        data provider for the estimated spatial resolution of the data product
        in arcseconds. This refers to the smallest spatial feature in the
        observed signal that can be resolved.
        In cases where the spatial resolution varies across the field the best
        spatial resolution (smallest resolvable spatial feature) should be
        specified. In cases where the spatial frequency sampling of an
        observation is complex (e.g., interferometry) a typical value for
        spatial resolution estimate should be given; additional
        characterisation may be necessary to fully specify the spatial
        characteristics of the data.
        """
        return self.get('s_resolution') * u.arcsec

    @property
    def s_xel1(self):
        """
        Number of elements along the first spatial axis
        """
        return self.get('s_xel1')

    @property
    def s_xel2(self):
        """
        Number of elements along the second spatial axis
        """
        return self.get('s_xel2')

    @property
    def s_ucd(self):
        """
        UCD for the nature of the spatial axis (pos or u,v data)
        """
        return self.get('s_ucd', decode=True, default=None)

    @property
    def s_unit(self):
        """
        Unit used for spatial axis
        """
        return self.get('s_unit', decode=True, default=None)

    @property
    def s_resolution_min(self):
        """
        Resolution min value on spatial axis (FHWM of PSF)
        """
        rmin = self.get('s_resolution_min', default=None)
        return rmin if not rmin else rmin * u.arcsec

    @property
    def s_resolution_max(self):
        """
        Resolution max value on spatial axis (FHWM of PSF)
        """
        rmax = self.get('s_resolution_max', default=None)
        return rmax if not rmax else rmax * u.arcsec

    @property
    def s_calib_status(self):
        """
        A string to encode the calibration status along the spatial axis
        (astrometry). Possible values could be {uncalibrated, raw, calibrated}
        """
        return self.get('s_calib_status', decode=True, default=None)

    @property
    def s_stat_error(self):
        """
        This parameter gives an estimate of the astrometric statistical error
        after the astrometric calibration phase.
        """
        return self.get('s_stat_error', decode=True, default=None)

    @property
    def s_pixel_scale(self):
        """
        This corresponds to the sampling precision of the data along the
        spatial axis. It is stored as a real number corresponding to the
        spatial sampling period, i.e., the distance in world coordinates
        system units between two pixel centers. It may contain two values if
        the pixels are rectangular.
        """
        return self.get('s_pixel_scale', decode=True, default=None)

    #           TIME CHARACTERISATION
    @property
    def t_xel(self):
        """
        Number of elements along the time axis
        """
        return self.get('t_xel')

    @property
    def t_ref_pos(self):
        """
        Time Axis Reference Position as defined in STC REC, Section 4.4.1.1.1
        """
        return self.get('t_ref_pos', decode=True, default=None)

    @property
    def t_min(self):
        """
        The start time of the observation specified in MJD as an
        `~astropy.time.Time` instance. In case of data products result of the
        combination of multiple frames, min time must be the minimum of the
        start times. ``None`` is used for NaN response values.
        """
        t_min = self.get('t_min')
        if np.isfinite(t_min):
            return Time(t_min, format='mjd')
        else:
            return None

    @property
    def t_max(self):
        """
        The stop time of the observation specified in MJD as an
        `~astropy.time.Time` instance. In case of data products result of the
        combination of multiple frames, t_max must be the maximum of the
        stop times. ``None`` is used for NaN response values.
        """
        t_max = self.get('t_max')
        if np.isfinite(t_max):
            return Time(t_max, format='mjd')
        else:
            return None

    @property
    def t_exptime(self):
        """
        Total exposure time. For simple exposures, this is just the time_bounds
         size expressed in seconds. For data where the detector is not active
         at all times (e.g. data products made by combining exposures taken at
         different times), the t_exptime will be smaller than the time_bounds
         interval. For data where the xptime is not constant over the entire
         data product, the median exposure time per pixel is a good way to
         characterize the typical value. In some cases, exptime is generally
         used as an indicator of the relative sensitivity (depth) within a
         single data collection (e.g. obs_collection); data providers should
         supply a suitable relative value when it is not feasible to define or
         compute the true exposure time.

        In case of targeted observations, on the contrary the exposure time is
        often adjusted to achieve similar signal to noise ratio for different
        targets.
        """
        return self.get('t_exptime') * u.second

    @property
    def t_resolution(self):
        """
        Estimated or average value of the temporal resolution.
        """
        return self.get('t_resolution') * u.second

    @property
    def t_calib_status(self):
        """
        Type of time coordinate calibration. Possible values are principally
        {uncalibrated, calibrated, raw, relative}. This may be extended for
        specific time domain collections.
        """
        return self.get('t_calib_status', decode=True, default=None)

    @property
    def t_stat_error(self):
        """
        Time coord statistical error on the time measurements in seconds
        """
        ter = self.get('t_stat_error', default=None)
        return ter if not ter else ter * u.second

    #           SPECTRAL CHARACTERISATION
    @property
    def em_xel(self):
        """
        Number of elements along the spectral axis
        """
        return self.get('em_xel')

    @property
    def em_ucd(self):
        """
        Nature of the spectral axis
        """
        return self.get('em_ucd', decode=True, default=None)

    @property
    def em_unit(self):
        """
        Units along the spectral axis
        """
        return self.get('em_unit', decode=True, default=None)

    @property
    def em_calib_status(self):
        """
        This attribute of the spectral axis indicates the status of the data
        in terms of spectral calibration. Possible values are defined in the
        Characterisation Data Model and belong to {uncalibrated , calibrated,
        relative, absolute}.
        """
        return self.get('em_calib_status', decode=True, default=None)

    @property
    def em_min(self):
        """
        Minimum of the spectral interval covered by the observation
        """
        return self.get('em_min') * u.meter

    @property
    def em_max(self):
        """
        Maximum of the spectral interval covered by the observation
        """
        return self.get('em_max') * u.meter

    @property
    def em_res_power(self):
        """
        Average estimation for the spectral resolution power stored as a
        double value, with no unit.
        """
        return self.get("em_res_power")

    @property
    def em_res_power_min(self):
        """
        Resolving power min value on spectral axis
        """
        return self.get('em_res_power_min', None)

    @property
    def em_res_power_max(self):
        """
        Resolving power max value on spectral axis
        """
        return self.get('em_res_power_max', None)

    @property
    def em_resolution(self):
        """
        A mean estimate of the resolution, e.g. Full Width at Half Maximum
        (FWHM) of the Line Spread Function (or LSF). This can be used for
        narrow range spectra whereas in the majority of cases, the resolution
        power is preferable due to the LSF variation along the spectral axis.
        """
        if 'em_resolution' in self.keys():
            return self.get('em_resolution') * u.meter
        return None

    @property
    def em_stat_error(self):
        """
        Spectral coord statistical error (accuracy along the spectral axis)
        """
        if 'em_stat_error' in self.keys():
            return self.get('em_stat_error') * u.meter
        return None

    #           OBSERVABLE AXIS
    @property
    def o_ucd(self):
        """
        Nature of the observable axis within the data product
        """
        return self.get('o_ucd', decode=True)

    @property
    def o_unit(self):
        """
        Units along the observable axis
        """
        return self.get('o_unit', decode=True, default=None)

    @property
    def o_calib_status(self):
        """
        Type of calibration applied on the Flux observed (or other observable
        quantity).
        """
        return self.get('o_calib_status', decode=True, default=None)

    @property
    def o_stat_error(self):
        """
        Statistical error on the Observable axis.
        Note: the return value has the units defined in unit
        """
        return self.get('o_stat_error', decode=True, default=None)

    #           POLARIZATION CHARACTERISATION
    @property
    def pol_xel(self):
        """
        Number of different polarization states present in the data. The
        default value is 0, indicating that polarization was not explicitly
        observed. Corresponding values are stored in the ``pol`` property
        """
        return self.get('pol_xel')

    @property
    def pol_states(self):
        """
        List of polarization states present in the data file. Possible values
        are: {I Q U V RR LL RL LR XX YY XY YX POLI POLA}. Values in the
        set are separated by the '/' character. A leading / character must
        start the list and a trailing / character must end it. It should be
        ordered following the above list, compatible with the FITS list table
        for polarization definition.
        """
        return self.get('pol_states', decode=True, default=None)

    #           PROVENANCE
    @property
    def instrument_name(self):
        """
        The name of the instrument used for the acquisition of the data
        """
        return self.get('instrument_name', decode=True)

    @property
    def facility_name(self):
        """
        Name of the facility or observatory used to collect the data
        """
        return self.get('facility_name', decode=True, default=None)

    @property
    def proposal_id(self):
        """
        Identifier of proposal to which observation belongs
        """
        return self.get('proposal_id', default=None, decode=True)


@deprecated("1.5", alternative="SIA2Service")
class SIAService(SIA2Service):
    pass


@deprecated("1.5", alternative="SIA2Query")
class SIAQuery(SIA2Query):
    pass


@deprecated("1.5", alternative="SIA2Results")
class SIAResults(SIA2Results):
    pass
