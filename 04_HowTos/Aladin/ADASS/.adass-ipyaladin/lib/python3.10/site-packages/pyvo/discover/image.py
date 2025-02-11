# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Global image discovery.

This module is intended to cover the much-requested use case "give me all
images with properties x, y, and z, where these properties for now is the
location in space, time, and spectrum; sensitivity clearly would be great,
but that needs standard work.

The code here looks for data in SIA1, SIA2, and Obscore services for now.
"""

# TODOs:
# * Evaluate and log Overflow conditions in the three protocols
# * do more (s_region!) on the basis of the WCS parts in SIA1
# * In obscore, we probably should make the dataproduct type constraint
#   configurable (this should really also work for cubes)
# * Obscore query generation *might* want some extra logic for inclusive,
#   (as in OR whatever IS NULL) -- but is that a good idea?
# * It would be nice if we preserved datalink availability and perhaps
#   even let people do automatic cutouts to the RoI.

import functools
import threading

import requests

from astropy import units as u
from astropy import table
from astropy import time

from ..dam import obscore
from .. import dal
from .. import registry
from ..registry import regtap


# imports for type hints
from typing import Callable, Generator, List, Optional, Set, Tuple
from astropy.units import quantity

__all__ = ["ImageFound", "ImageDiscoverer", "images_globally"]


# We should probably have a general way to set query timeouts in pyVO.
# For now, we don't, but for global discovery there's really no way
# around them.  So, let me hack something here.

class SessionWithTimeout(requests.Session):
    def __init__(self, *args, default_timeout=None, **kwargs):
        self.default_timeout = default_timeout
        super().__init__(*args, **kwargs)

    def request(self, *args, **kwargs):
        if "timeout" not in kwargs:
            kwargs["timeout"] = self.default_timeout
        return super().request(*args, **kwargs)


class Queriable:
    """A facade for a queriable service.

    We keep these rather than work on
    `pyvo.registry.RegistryResource`-s directly because the latter
    actually live in VOTables which are very hard to manipulate.

    They are constructed with a resource record.
    """
    def __init__(self, res_rec):
        self.res_rec = res_rec
        self.ivoid = res_rec.ivoid
        self.title = self.res_rec.res_title

    def __str__(self):
        return f"<{self.ivoid}>"

    def __repr__(self):
        return str(self)


@functools.lru_cache(maxsize=None)
def obscore_column_names():
    """returns the names of obscore columns.

    For lack of alternatives, I pull them out of ObsCoreMetadata for now.
    """
    return [name for name in dir(obscore.ObsCoreMetadata())
        if not name.startswith("_")]


class ImageFound(obscore.ObsCoreMetadata):
    """Obscore metadata for a found record.

    The mandatory obscore fields are kept as attributes.

    We are pulling these out from the various VOTables that we
    retrieve because we need to do some manipulation of them
    that is simpler to do if they are proper Python objects.

    This also keeps track of the service the record came from,
    which is available from the origin_service attribute.
    """
    attr_names = set(obscore_column_names())

    def __init__(self, origin_service, obscore_record):
        self.origin_service = origin_service
        for k, v in obscore_record.items():
            if k not in self.attr_names:
                raise TypeError(
                    f"ImageFound does not accept {k} keys")
            setattr(self, k, v)

    @classmethod
    def from_obscore_recs(cls, origin_service, obscore_result):
        if not obscore_result:
            return

        ap_table = obscore_result.to_table()
        our_keys = [n for n in ap_table.colnames if n in cls.attr_names]
        for row in obscore_result:
            yield cls(
                origin_service,
                dict(zip(our_keys, (row[n] for n in our_keys))))

    @classmethod
    def from_sia1_recs(cls, origin_service, sia1_result, filter_func):
        for rec in sia1_result:
            if not filter_func(rec):
                continue
            mapped = {
                "dataproduct_type": "image" if rec.naxes == 2 else "cube",
                "access_url": rec.acref,
                "em_max": None if rec.bandpass_hilimit is None
                    else rec.bandpass_hilimit.to(u.m).value,
                "em_min": None if rec.bandpass_lolimit is None
                    else rec.bandpass_lolimit.to(u.m).value,
                # Sigh.  Try to guess exposure time?
                "t_min": None if rec.dateobs is None
                    else rec.dateobs.mjd,
                "t_max": None if rec.dateobs is None
                    else rec.dateobs.mjd,
                "access_estsize": None if rec.filesize is None
                    else rec.filesize/1024,
                "access_format": rec.format,
                "instrument_name": rec.instr,
                "s_xel1": rec.naxis[0].to(u.pix).value,
                "s_xel2": rec.naxis[1].to(u.pix).value,
                "s_ra": rec.pos.icrs.ra.to(u.deg).value,
                "s_dec": rec.pos.icrs.dec.to(u.deg).value,
                "obs_title": rec.title,
            }
            yield cls(origin_service, mapped)


def _clean_for(records: List[Queriable], ivoids_to_remove: Set[str]):
    """returns the Queriables in records the ivoids of which are
    not in ivoids_to_remove.
    """
    return [r for r in records if r.ivoid not in ivoids_to_remove]


class ImageDiscoverer:
    """A management class for VO global image discovery.

    This encapsulates all of constraints, service lists, results, and
    diagnostics.  This probably should not be considered API but
    rather as an implementation detail of discover_images.

    The normal usage is do call discover_services(), which will locate
    all VO services that may have relevant data.  Alternatively, call
    set_services(registry_results) with some result of a registry.search()
    call.  ImageDiscoverer will then pick capabilities it can use out
    of the resource records.  Records without usable capabilities are
    silently ignored.

    Then call query_services to execute the discovery query on these
    services.

    See images_globally for a discussion of its constructor parameters.
    """
    # Constraint defaults
    # a float in metres
    spectrum = None
    # MJD floats
    time_min = time_max = None
    # a center as a 2-tuple in ICRS degrees
    center = None
    # a radius as a float in degrees
    radius = None

    def __init__(self, *,
            space=None, spectrum=None, time=None,
            inclusive=False,
            watcher=None,
            timeout=20):
        self.session = SessionWithTimeout(default_timeout=timeout)

        if space:
            self.center = (space[0], space[1])
            self.radius = space[2]

        if spectrum is not None:
            self.spectrum = spectrum.to(u.m, equivalencies=u.spectral()).value

        if time is not None:
            if isinstance(time, tuple):
                self.time_min, self.time_max = time[0].mjd, time[1].mjd
            else:
                self.time_min = self.time_max = time.mjd

        self.inclusive = inclusive
        self.already_queried, self.failed_services = 0, 0
        self.results: List[obscore.ObsCoreMetadata] = []
        self.watcher = watcher
        self.log_messages: List[str] = []
        self.known_access_urls: Set[str] = set()

        self._service_list_lock = threading.Lock()
        with self._service_list_lock:
            # only write to these while holding the lock
            self.sia1_recs, self.sia2_recs, self.obscore_recs = [], [], []

    def _info(self, message: str) -> None:
        """sends message to our watcher (if there is any)
        """
        if self.watcher is not None:
            self.watcher(self, message)

    def _log(self, message: str) -> None:
        """logs message.

        This will also do whatever _info does.
        """
        self.log_messages.append(message)
        self._info(message)

    def _purge_redundant_services(self):
        """removes services querying data already covered by more capable
        services from our current services lists.
        """
        def ids(recs):
            return set(r.ivoid for r in recs)

        self.sia1_recs = _clean_for(self.sia1_recs,
                ids(self.sia2_recs) | ids(self.obscore_recs))
        self.sia2_recs = _clean_for(self.sia2_recs, ids(self.obscore_recs))

        # In addition, now throw out all services that have an
        # IsServedBy relationship to another service we will also
        # query.  That's particularly valuable if there are large
        # obscore services covering data from many SIA1 services.
        ids_present = table.Table([
            table.Column(
                name="id",
                data=list(
                    sorted(ids(self.sia1_recs)
                        | ids(self.sia2_recs)
                        | ids(self.obscore_recs))),
                description="ivoids of candiate services",
                meta={"ucd": "meta.ref.ivoid"}),])
        if len(ids_present) == 0:
            return

        services_for = regtap.get_RegTAP_service().run_sync(
            """SELECT ivoid, related_id
                FROM rr.relationship
                  JOIN tap_upload.ids AS leftids ON (ivoid=leftids.id)
                  JOIN tap_upload.ids AS rightids ON (related_id=rightids.id)
                WHERE relationship_type='isservedby'
            """, uploads={'ids': ids_present})
        for rec in services_for:
            self._log(f"Skipping {rec['ivoid']} because"
                f" it is served by {rec['related_id']}")

        collections_to_remove = set(r["ivoid"] for r in services_for)
        self.sia1_recs = _clean_for(self.sia1_recs, collections_to_remove)
        self.sia2_recs = _clean_for(self.sia2_recs, collections_to_remove)
        self.obscore_recs = _clean_for(self.obscore_recs, collections_to_remove)

    def _discover_obscore_services(self, *constraints):
        # For obscore, we currently have a defunct discovery pattern
        # ("obscore" in the Datamodel constraint).  There is obscore-new,
        # which fixes the problem, but until that's adopted by all the
        # obscore services, we have to try both and the pick the
        # more suitable version.
        # Once we move obscore-new to obscore, remove this function
        # and put
        # self.obscore_recs = [Queriable(r) for r in registry.search(
        #   registry.Datamodel("obscore"), *constraints)]
        # back into discover_services.
        obscore_services = registry.search(
            registry.Datamodel("obscore_new"), *constraints)
        tap_services_with_obscore = registry.search(
            registry.Datamodel("obscore"), *constraints)

        new_style_access_urls = set()
        for rec in obscore_services:
            new_style_access_urls |= set(
                i.access_url for i in rec.list_interfaces("tap"))

        for tap_rec in tap_services_with_obscore:
            access_urls = set(
                i.baseurl for i in rec.list_services("tap"))
            if new_style_access_urls.isdisjoint(access_urls):
                obscore_services.append(obscore_services)

        return [Queriable(r) for r in obscore_services]

    def discover_services(self):
        """fills the X_recs attributes with resources declaring coverage
        for our constraints.

        X at this point is sia1, sia2, and obscore.

        It tries to filter out as many duplicates (i.e., services operating on
        the same data collections) as it can.  The order of preference is
        Obscore, SIA2, SIA.
        """
        constraints = []
        if self.center is not None:
            constraints.append(
                registry.Spatial(list(self.center)+[self.radius],
                    inclusive=self.inclusive))
        if self.spectrum is not None:
            constraints.append(
                registry.Spectral(self.spectrum*u.m,
                    inclusive=self.inclusive))
        if self.time_min is not None or self.time_max is not None:
            constraints.append(
                registry.Temporal(
                    (self.time_min, self.time_max),
                    inclusive=self.inclusive))

        with self._service_list_lock:
            self.sia1_recs = [Queriable(r) for r in registry.search(
                registry.Servicetype("sia"), *constraints)]
            self._info("Found {} SIA1 service(s)".format(len(self.sia1_recs)))

            self.sia2_recs = [Queriable(r) for r in registry.search(
                registry.Servicetype("sia2"), *constraints)]
            self._info("Found {} SIA2 service(s)".format(len(self.sia2_recs)))

            self.obscore_recs = self._discover_obscore_services(*constraints)
            self._info("Found {} Obscore service(s)".format(
                len(self.obscore_recs)))

            self._purge_redundant_services()

    def set_services(self,
            registry_results: registry.RegistryResults,
            purge_redundant: bool = True) -> None:
        """as an alternative to discover_services, this sets the services
        to be queried to the result of a custom registry query.

        This will pick the "most capabable" interface from each record
        and ignore records without image discovery capabilities.

        If you set purge_redundant to false, this will not attempt
        to eliminate services that are alternative interfaces to services
        that are already called.  There are very few situations in which
        you would want to do that, mostly related to debugging.
        """
        with self._service_list_lock:
            for rsc in registry_results:
                if "tap" in rsc.access_modes():
                    # TODO: we ought to ensure there's an obscore
                    # table on this; but by the time this sees users,
                    # I hope to have fixed obscore discovery.
                    self.obscore_recs.append(Queriable(rsc))
                elif "sia2" in rsc.access_modes():
                    self.sia2_recs.append(Queriable(rsc))
                elif "sia" in rsc.access_modes():
                    self.sia1_recs.append(Queriable(rsc))
                # else ignore this record

            if purge_redundant:
                self._purge_redundant_services()

    def reset_services(self):
        """clears the queues of services to query.

        This is the only supported way to interrupt operations once
        query_services has been called.

        This will not interrupt the query currently running.  There is
        currently no way to do that.
        """
        with self._service_list_lock:
            self.sia1_recs, self.sia2_recs, self.obscore_recs = [], [], []
            self._log("Cancelling queries with {} service(s) queried"
                .format(self.already_queried))

    def _add_records(self,
            recgen: Generator[ImageFound, None, None]) -> int:
        """adds records from regen to the global results.

        This will skip datasets the access urls of which we have already seen
        and will return the number of datasets actually added.
        """
        n_added = 0

        for obscore_record in recgen:
            if obscore_record.access_url in self.known_access_urls:
                continue
            self.known_access_urls.add(obscore_record.access_url)
            self.results.append(obscore_record)
            n_added += 1

        return n_added

    def _query_one_sia1(self, rec: Queriable):
        """runs our query against a SIA1 capability of rec.

        Since SIA1 cannot do spectral and temporal constraints itself,
        we do them client-side provided sufficient metadata is present.
        If metadata is missing, we keep or discard accoding to
        self.inclusive.
        """
        def non_spatial_filter(sia1_rec):
            if self.spectrum and not self.inclusive and (
                    sia1_rec.bandpass_hilimit is not None
                    and sia1_rec.bandpass_lolimit is not None):
                if not (sia1_rec.bandpass_lolimit
                        <= self.spectrum*u.m
                        <= sia1_rec.bandpass_hilimit):
                    return False

            # Regrettably, the exposure time is not part of SIA1 standard
            # metadata.  We fudge things a bit; this should not
            # increase false positives very badly.
            if (self.time_min is not None or self.time_max is not None) \
                    and not self.inclusive and sia1_rec.dateobs:
                if not (self.time_min-0.1
                        < sia1_rec.dateobs.mjd
                        < self.time_max+0.1):
                    return False
            return True

        self._info("Querying SIA1 {}...".format(rec.title))
        svc = rec.res_rec.get_service("sia", session=self.session, lax=True)
        n_found = self._add_records(
            ImageFound.from_sia1_recs(
                rec.ivoid,
                svc.search(
                    pos=self.center, size=self.radius, intersect='overlaps'),
                non_spatial_filter))
        self._log(f"SIA1 {rec.title} {n_found} records")

    def _query_sia1(self):
        """runs the SIA1 part of our discovery.

        This will be a no-op without a space constraint due to
        limitations of SIA1.
        """
        if self.sia1_recs and self.center is None:
            self._log("SIA1 service(s) skipped due to missing space"
                " constraint")
            return

        # we don't do a for loop here because we want to react to changes
        # in self.sia1_recs
        while self.sia1_recs:
            rec = self.sia1_recs.pop()
            try:
                self._query_one_sia1(rec)
            except Exception as msg:
                self._log(f"SIA1 {rec.ivoid} skipped: {msg}")
                self.failed_services += 1
            self.already_queried += 1

    def _query_one_sia2(self, rec: Queriable):
        """runs our query against a SIA2 capability of rec.
        """
        self._info("Querying SIA2 {}...".format(rec.title))

        svc = rec.res_rec.get_service("sia2", session=self.session, lax=True)
        constraints = {}
        if self.center is not None:
            constraints["pos"] = self.center+(self.radius,)
        if self.spectrum is not None:
            constraints["band"] = self.spectrum
        if self.time_min is not None or self.time_max is not None:
            constraints["time"] = (
                time.Time(self.time_min, format="mjd"),
                time.Time(self.time_max, format="mjd"))

        n_found = self._add_records(
            ImageFound.from_obscore_recs(
                rec.ivoid,
                svc.search(**constraints)))
        self._log(f"SIA2 {rec.title}: {n_found} records")

    def _query_sia2(self):
        """runs the SIA2 part of our discovery.
        """
        while self.sia2_recs:
            rec = self.sia2_recs.pop()
            try:
                self._query_one_sia2(rec)
            except Exception as msg:
                self._log(f"SIA2 {rec.ivoid} skipped: {msg}")
                self.failed_services += 1
            self.already_queried += 1

    def _query_one_obscore(self, rec: Queriable, where_clause: str):
        """runs our query against a Obscore capability of rec.
        """
        self._info("Querying Obscore {}...".format(rec.title))
        svc = rec.res_rec.get_service("tap", session=self.session, lax=True)

        n_found = self._add_records(
            ImageFound.from_obscore_recs(
                rec.ivoid,
                svc.run_sync("select * from ivoa.obscore "+where_clause)))
        self._log(f"Obscore {rec.title}: {n_found} records")

    def _query_obscore(self):
        """runs the Obscore part of our discovery.
        """
        where_parts = ["dataproduct_type='image'"]
        if self.center is not None:
            where_parts.append(
                "(1=contains(point('ICRS', s_ra, s_dec),"
                " circle('ICRS', {}, {}, {}))".format(
                    self.center[0], self.center[1], self.radius)
                + " or 1=intersects(circle({}, {}, {}), s_region))".format(
                    self.center[0], self.center[1], self.radius))
        if self.spectrum is not None:
            where_parts.append(
                f"(em_min<={self.spectrum} AND em_max>={self.spectrum})")

        if self.time_min is not None or self.time_max is not None:
            where_parts.append(
                "({h1}>={l2} AND {h2}>={l1}"
                " AND {l1}<={h1} AND {l2}<={h2})".format(
                    l1="t_min", h1="t_max",
                    l2=self.time_min, h2=self.time_max))

        where_clause = "WHERE "+(" AND ".join(where_parts))

        while self.obscore_recs:
            rec = self.obscore_recs.pop()
            try:
                self._query_one_obscore(rec, where_clause)
            except Exception as msg:
                self._log("Obscore {} skipped: {}".format(
                    rec.res_rec['ivoid'], msg))
                self.failed_services += 1
            self.already_queried += 1

    def get_query_stats(self):
        """returns a tuple of n(total to query), n(already queried)
        """
        total_to_query = len(self.obscore_recs)\
            + len(self.sia1_recs)\
            + len(self.sia2_recs)
        return total_to_query, self.already_queried, self.failed_services

    def query_services(self):
        """queries the discovered image services according to our
        constraints.

        This creates and fills the results and the log_messages attributes.
        """
        if (not self.sia1_recs
                and not self.sia2_recs
                and not self.obscore_recs):
            raise dal.DALQueryError("No services to query.  Unless"
                " you overrode service selection, you will have to"
                " loosen your constraints.")
        self._query_obscore()
        self._query_sia2()
        self._query_sia1()


def images_globally(
        *,
        space: Optional[Tuple[float, float, float]] = None,
        spectrum: Optional[quantity.Quantity] = None,
        time: Optional[time.Time] = None,
        inclusive: bool = False,
        watcher: Optional[Callable[['ImageDiscoverer', str], None]] = None,
        timeout: float = 20,
        services: Optional[registry.RegistryResults] = None)\
        -> Tuple[List[obscore.ObsCoreMetadata], List[str]]:
    """returns a collection of ObsCoreMetadata-s matching certain constraints
    and a list of log lines.

    Parameters
    ----------

    space :
        An optional tuple of ra, dec, and the search radius, all in degrees;
        images returned must intersect a spherical circle described in this
        way.
    spectrum :
        An astropy quantity convertible to a (vacuum) wavelength; images
        must cover this point in the (electromagnetic) spectrum.
    time :
        An astropy time that must be in the observation time of the image
        (if it declares a time).
    inclusive :
        Set to True to incluse services that do not declare their
        STC coverage.  As of 2024, it may be advisable to do that as many
        relevant archives do not do that.
    watcher :
        A callable that will be called with the ImageDiscoverer instance and
        a string perhaps suitable for displaying to a human.
    services :
        An optional `~pyvo.registry.RegistryResults` instance to
        override automatic services detection.

    When an image has insufficient metadata to evaluate a constraint, it
    is excluded; this mimics the behaviour of SQL engines that consider
    comparisons with NULL-s false.

    Returns
    -------
    discovered_images, log_lines
        The images found are returned in a list of `pyvo.discover.ImageFound`
        instances.
        log_lines is a list of strings reporting which services were
        queried with which result (and possibly more).  So far, this
        is not considered machine-readable.
    """
    discoverer = ImageDiscoverer(
        space=space, spectrum=spectrum, time=time,
        inclusive=inclusive,
        watcher=watcher,
        timeout=timeout)

    if services is None:
        discoverer.discover_services()
    else:
        discoverer.set_services(services)

    discoverer.query_services()
    return discoverer.results, discoverer.log_messages
