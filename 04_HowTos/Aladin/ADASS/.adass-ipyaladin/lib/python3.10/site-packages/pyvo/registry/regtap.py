# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
a module for basic VO Registry interactions.

A VO registry is a database of VO resources--data collections and
services--that are available for VO applications.  Typically, it is
aware of the resources from all over the world.  A registry can find
relevant data collections and services through search
queries--typically, subject-based.  The registry responds with a list
of records describing matching resources.  With a record in hand, the
application can use the information in the record to access the
resource directly.  Most often, the resource is a data service that
can be queried for individual datasets of interest.

This module provides basic, low-level access to the RegTAP Registries using
standardized TAP-based services.
"""

import functools
import itertools
import os
import textwrap
import warnings

from astropy import table
from astropy.utils.decorators import deprecated
from astropy.utils.exceptions import AstropyDeprecationWarning

import numpy

from . import rtcons
from ..dal import scs, sia, sia2, ssa, sla, tap, query as dalq
from ..io.vosi import vodataservice
from ..utils.formatting import para_format_desc


__all__ = ["search", "get_RegTAP_query", "Interface",
           "RegistryResource", "RegistryResults", "ivoid2service"]

REGISTRY_BASEURL = os.environ.get("IVOA_REGISTRY", "http://reg.g-vo.org/tap"
                                  ).rstrip("/")


# ADQL only has string_agg, where we need string arrays.  We fake arrays
# by joining elements with a token separator that we think shouldn't
# turn up in the things joined.  Of course, people could create
# resources that break us; let's assume there's nothing be gained
# from that ever.
TOKEN_SEP = ":::py VO sep:::"


def shorten_stdid(s):
    """removes leading ivo://ivoa.net/std/ from s if present.

    We're using this to make the display and naming of standard ivoids
    less ugly in several places.

    Nones remain Nones.
    """
    if s and s.startswith("ivo://ivoa.net/std/"):
        return s[19:]
    return s


def expand_stdid(s):
    """returns s if it already looks like a URI, and it prepends
    ivo://ivoa.net/std otherwise.

    This is the (approximate) reverse of shorten_stdid.
    """
    if s is None or "://" in s:
        return s
    return "ivo://ivoa.net/std/" + s


def regularize_SIA2_id(standard_id):
    """returns standard_id with SIA2 standard ids modified to what they should
    have been.

    Regrettably, SIA2 uses the same standardID as SIA1; sure, they use
    different fragments, but that doesn't really help with the logic
    we have here.

    To make up for that, we replace them with the sia2 ids they should
    have had on input.  This function assumes lowercased ids as they
    come from RegTAP services.
    """
    if standard_id.startswith("ivo://ivoa.net/std/sia#query-2"):
        return "ivo://ivoa.net/std/sia2"
    elif standard_id.startswith("ivo://ivoa.net/std/sia#query-aux-2"):
        # query-aux-2 is mentioned in discovering data collections,
        # which isn't really the place to define this.  But then
        # it's endorsed, and SIA2 doesn't say anything about it.
        return "ivo://ivoa.net/std/sia2#aux"
    else:
        return standard_id


@functools.lru_cache(1)
def get_RegTAP_service():
    """
    a lazily created TAP service offering the RegTAP services.

    Always get the TAP service there using this function to avoid
    re-creating the server and profit from caching of capabilties,
    tables, etc.

    To switch to a different RegTAP service, use
    :py:func:`choose_RegTAP_service`.
    """
    return tap.TAPService(REGISTRY_BASEURL)


def choose_RegTAP_service(access_url):
    """
    changes the RegTAP service used by :py:func:`search`
    to the one at access_url.

    By default, pyVO uses whatever is given in the environment variable
    ``IVOA_REGISTRY``, defaulting to GAVO's TAP service.  In order to
    change the service used on the fly, always use this function in order
    to clear caches that need clearing.

    Parameters
    ----------
    access_url : str
        The TAP access URL of the new RegTAP endpoints.
        To find alternate endpoints, try ``regsearch(datamodel='regtap')``
        and look at ``.get_interface("tap").access_url`` of the results.
    """
    global REGISTRY_BASEURL
    get_RegTAP_service.cache_clear()
    REGISTRY_BASEURL = access_url


def get_RegTAP_query(*constraints: rtcons.Constraint,
                     includeaux=False,
                     service=None,
                     **kwargs):
    """returns SQL for a RegTAP query for constraints and keywords.

    This function's parameters are as for search; this is basically
    a wrapper for rtcons.build_regtap_query maintaining the legacy
    keyword-based interface.
    """
    # we don't document the service parameter -- it's probably not useful
    # to users and is just the conscequence of having retrofitted service
    # sensing into the API.
    if service is None:
        service = get_RegTAP_service()

    constraints = list(constraints) + rtcons.keywords_to_constraints(kwargs)

    # maintain legacy includeaux by locating any Servicetype constraints
    # and replacing them with ones that includes auxiliaries.
    if includeaux:
        for index, constraint in enumerate(constraints):
            if isinstance(constraint, rtcons.Servicetype):
                constraints[index] = constraint.include_auxiliary_services()

    return rtcons.build_regtap_query(constraints, service)


def search(*constraints: rtcons.Constraint,
        includeaux: bool = False,
        maxrec: int = None,
        **kwargs):

    """
    execute a simple query to the RegTAP registry.

    The function accepts query constraints either as Constraint objects
    passed in as positional arguments or as their associated keywords.
    For what constraints are available, see :ref:`registry-basic-interface`.


    The values of keyword arguments may be tuples or lists when the associated
    Constraint objects take multiple arguments.

    All constraints, whether passed in directly or via keywords, are
    evaluated as a conjunction (i.e., in an AND clause).

    Parameters
    ----------
    *constraints : `~pyvo.registry.Constraint` instances
        The constraints (keywords to match, positions to cover, ...)
        that the returned records need to satisfy.
        The accepted constraints are:

            - keywords: one or more freetext words, mached in the title,
              description or subject of the resource.
            - servicetype: constrain to one of tap, ssa, sia, conesearch
              (or full ivoids for other service types). This is the
              constraint you want to use for service discovery.
            - ucd: constrain by one or more UCD patterns; resources match
              when they serve columns having a matching UCD
              (e.g., phot.mag;em.ir.% for "any infrared magnitude").
            - waveband: one or more terms from the vocabulary at
              http://www.ivoa.net/rdf/messenger giving the rough spectral
              location of the resource.
            - author: an author ("creator"). This is a single SQL pattern,
              and given the sloppy practices in the VO for how to write
              author names, you should probably generously use wildcards.
            - datamodel: one of obscore, epntap, or regtap: only return TAP
              services having tables of this kind.
            - ivoid: exactly match a single IVOA identifier (that is,
              in effect, the primary key in the VO).
            - spatial: match resources covering a certain geometry (point,
              circle, polygon, or MOC). RegTAP 1.2 Extension.
            - spectral: match resources covering a certain part of the spectrum
              (usually, but not limited to, the electromagnetic spectrum).
              RegTAP 1.2 Extension
            - temporal: match resources covering a some point or interval in
              time. RegTAP 1.2 Extension

        Multiple constraints are combined conjunctively ("AND").


    includeaux : bool
        Flag for whether to include auxiliary capabilities in results.
        This may result in duplicate capabilities being returned,
        especially if the servicetype is not specified.

    maxrec : int
        Overrides the RegTAP server's default limit on the number of rows to
        return.  You may need to use this if you want to retrieve more
        than a few thousand matches.  The server may also have a hard limit
        that ``maxrec`` cannot override.  Note that truncated search results
        are not reproducible.

    **kwargs : strings, mostly
        shorthands for ``constraints``; see the documentation of
        a specific constraint for what keyword it uses and what literal
        it expects.

    Returns
    -------
    ~pyvo.registry.RegistryResults`
       a container holding a table of matching resource (e.g. services)

    See Also
    --------
    RegistryResults
    """
    service = get_RegTAP_service()
    query = RegistryQuery(
        service.baseurl,
        get_RegTAP_query(*constraints,
            includeaux=includeaux,
            service=service,
            **kwargs),
        maxrec=maxrec)
    return query.execute()


class RegistryQuery(tap.TAPQuery):
    def execute(self):
        """
        submit the query and return the results as a RegistryResults instance

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
        return RegistryResults(self.execute_votable(), url=self.queryurl)


class RegistryResults(dalq.DALResults):
    """
    an iterable set of results from a registry query. Each record is
    returned as RegistryResults

    You can iterate over these, or access them by (numeric) index; note,
    however, that these indexes will not be stable across different
    executions and thus should only be used in interactive sessions.
    Alternatively, you can use short names as indexes; there *might*
    be clashes for these, as they are not unique VO-wide.  Where this
    matters, you need to use full ivoids as index.

    """

    def getrecord(self, index):
        """
        return all the attributes of a resource record with the given index
        as SimpleResource instance (a dictionary-like object).

        Parameters
        ----------
        index  : int
            the zero-based index of the record
        """
        return RegistryResource(self, index)

    def get_summary(self):
        """
        returns a brief overview of the matched results as an astropy table.

        This is mainly intended for interactive use, where people would
        like to inspect the matches in, perhaps, notebooks.
        """
        return table.Table([
            list(range(len(self))),
            [r.short_name for r in self],
            [r.res_title for r in self],
            [r.res_description for r in self],
            [", ".join(sorted(r.access_modes())) for r in self]],
            names=("index", "short_name", "title", "description", "interfaces"),
            descriptions=(
                "Index to access the resource within self",
                "Short name",
                "Resource title",
                "Resource description",
                "Access modes offered"))

    @functools.lru_cache(maxsize=None)
    def _get_ivo_index(self):
        return dict((r.ivoid, index)
                    for index, r in enumerate(self))

    @functools.lru_cache(maxsize=None)
    def _get_short_name_index(self):
        return dict((r.short_name, index)
                    for index, r in enumerate(self))

    def __getitem__(self, item):
        """
        returns a record by numeric index, short names, or ivoid.

        This will raise an IndexError or a KeyError when item does
        not match a record returned.
        """
        if isinstance(item, int):
            return self.getrecord(item)

        elif isinstance(item, str):
            if item.startswith("ivo://"):
                return self.getrecord(self._get_ivo_index()[item])
            else:
                return self.getrecord(self._get_short_name_index()[item])

        else:
            raise IndexError(f"No resource matching {item}")


class _BrowserService:
    """A pseudo-service class just opening a web browser for browser-based
    services.
    """

    def __init__(self, access_url, capability_description=None):
        self.baseurl = access_url
        self.capability_description = capability_description

    def __repr__(self):
        return f"BrowserService(baseurl : '{self.baseurl}', description : '{self.capability_description}')"

    def search(self):
        import webbrowser
        webbrowser.open(self.baseurl, 2)


class Interface:
    """
    a service interface.

    These consist of an access URL, a standard id for the capability
    (typically the ivoid of an IVOA standard, or None for free services),
    an interface type (something like vs:paramhttp or vr:webbrowser)
    and an indication if the interface is the "standard" interface
    of the capability.

    Such interfaces can be turned into services using the ``to_service``
    method if pyvo knows how to talk to the interface.

    Note that the constructor arguments are assumed to be normalised
    as in regtap (e.g., lowercased for the standardIDs).
    """
    service_for_standardid = {
        "ivo://ivoa.net/std/conesearch": scs.SCSService,
        "ivo://ivoa.net/std/sia": sia.SIAService,
        "ivo://ivoa.net/std/sia2": sia2.SIA2Service,
        "ivo://ivoa.net/std/ssa": ssa.SSAService,
        "ivo://ivoa.net/std/sla": sla.SLAService,
        "ivo://ivoa.net/std/tap": tap.TAPService}

    def __init__(self, access_url, *, standard_id=None,
                 intf_type=None, intf_role=None, capability_description=None):
        self.access_url = access_url
        self.standard_id = standard_id
        self.type = intf_type
        self.role = intf_role
        self.capability_description = capability_description
        self.is_standard = self.role == "std"

        if self.standard_id is not None:
            self.is_vosi = self.standard_id.startswith("ivo://ivoa.net/std/vosi")
        else:
            self.is_vosi = False

        # a service is user visible if it has a corresponding service class
        if self.standard_id is not None and self.standard_id != "":
            service_type = self.standard_id.split("#")[0]  # remove possible suffixes/standard keys
            self.is_user_visible = service_type in self.service_for_standardid
        # or if it is a webpage
        else:
            self.is_user_visible = self.type == "vr:webbrowser"

    def __repr__(self):
        if not self.standard_id or self.standard_id is None:
            return (f"Interface(type={self.type!r}, "
                    f"description={self.capability_description!r}, "
                    f"url={self.access_url!r})")
        return (f"Interface(type={self.standard_id.rsplit('/')[-1]!r}, "
                f"description={self.capability_description!r}, "
                f"url={self.access_url!r})")

    def to_service(self, *, session=None):
        if self.type == "vr:webbrowser":
            return _BrowserService(self.access_url, self.capability_description)

        if self.standard_id is None or not self.is_standard:
            raise ValueError("This is not a standard interface.  PyVO"
                             " cannot speak to it.")

        service_class = self.service_for_standardid.get(
            self.standard_id.split("#")[0])
        if service_class is None:
            raise ValueError("PyVO has no support for interfaces with"
                             f" standard id {self.standard_id}.")

        if service_class == sia2.SIA2Service:
            return service_class(self.access_url,
                                 capability_description=self.capability_description,
                                 check_baseurl=False,
                                 session=session)
        else:
            return service_class(
                self.access_url,
                capability_description=self.capability_description,
                session=session)

    def supports(self, standard_id):
        """returns true if we believe the interface should be able to talk
        standard_id.

        At this point, we naively check if the interfaces's standard_id
        has standard_id as a prefix.  At this point, we cut off standard_id
        fragments for this purpose.  This works for all current DAL
        standards but would, for instance, not work for VOSI.  Hence,
        this may need further logic if we wanted to extend our service
        generation to VOSI or, perhaps, VOSpace.

        Parameters
        ----------

        standard_id : str
            The ivoid of a standard.
        """
        if not self.standard_id:
            return False
        standard_id = regularize_SIA2_id(standard_id)
        return self.standard_id.split("#")[0] == standard_id.split("#")[0]


class RegistryResource(dalq.Record):
    """
    a dictionary for the resource metadata returned in one record of a
    registry query.

    A SimpleResource acts as a dictionary, so in general, all attributes can
    be accessed by name via the [] operator, and the attribute names can
    by returned via the keys() function.  For convenience, it also stores
    key values as properties; these include:
    """

    _service = None

    # the following attribute is used by datasearch._build_regtap_query
    # to figure build the select clause; it is maintained here
    # because this class knows what it expects to get.
    #
    # Each item is either a plain string for a column name, or
    # a 2-tuple for an as clause; all plain strings are used
    # used in the group by, and so it is assumed they are
    # 1:1 to ivoid.
    expected_columns = [
        "ivoid",
        "res_type",
        "short_name",
        "res_title",
        "content_level",
        "res_description",
        "reference_url",
        "creator_seq",
        "created",
        "updated",
        "rights",
        "content_type",
        "source_format",
        "source_value",
        "region_of_regard",
        "waveband",
        (f"\n  ivo_string_agg(COALESCE(access_url, ''), '{TOKEN_SEP}')",
            "access_urls"),
        (f"\n  ivo_string_agg(COALESCE(standard_id, ''), '{TOKEN_SEP}')",
            "standard_ids"),
        (f"\n  ivo_string_agg(COALESCE(intf_type, ''), '{TOKEN_SEP}')",
            "intf_types"),
        (f"\n  ivo_string_agg(COALESCE(intf_role, ''), '{TOKEN_SEP}')",
            "intf_roles"),
        (f"\n  ivo_string_agg(COALESCE(cap_description, ''), '{TOKEN_SEP}')",
            "cap_descriptions")]

    def __init__(self, results, index, *, session=None):
        dalq.Record.__init__(self, results, index, session=session)

        self._mapping["access_urls"
                      ] = self._parse_pseudo_array(self._mapping["access_urls"])
        self._mapping["standard_ids"] = [
            regularize_SIA2_id(id) for id in
            self._parse_pseudo_array(self._mapping["standard_ids"])]
        self._mapping["intf_types"
                      ] = self._parse_pseudo_array(self._mapping["intf_types"])
        self._mapping["intf_roles"
                      ] = self._parse_pseudo_array(self._mapping["intf_roles"])
        self._mapping["cap_descriptions"
                      ] = self._parse_pseudo_array(self._mapping["cap_descriptions"])

        self.interfaces = [Interface(props[0], standard_id=props[1], intf_type=props[2],
                                     intf_role=props[3], capability_description=props[4])
                           for props in itertools.zip_longest(
            self["access_urls"],
            self["standard_ids"],
            self["intf_types"],
            self["intf_roles"],
            self["cap_descriptions"])]

    @staticmethod
    def _parse_pseudo_array(literal):
        """
        parses RegTAP pseudo-arrays into lists.

        Parameters
        ----------
        literal : str
            the result of an ivo_string_agg call with TOKEN_SEP

        Returns
        -------
        A list of strings corresponding to the orginal, database-side
        aggregate.
        """
        if not literal:
            # As VOTable, we don't distinguish between None and ""
            return []
        return literal.split(TOKEN_SEP)

    @property
    def ivoid(self):
        """
        the IVOA identifier for the resource.
        """
        return self.get("ivoid", decode=True)

    @property
    def res_type(self):
        """
        the resource types that characterize this resource.
        """
        return self.get("res_type", decode=True)

    @property
    def short_name(self):
        """
        the short name for the resource
        """
        return self.get("short_name", decode=True)

    @property
    def res_title(self):
        """
        the title of the resource
        """
        return self.get("res_title", default=None, decode=True)

    @property
    def content_levels(self):
        """
        a list of content level labels that describe the intended audience
        for this resource.
        """
        return self.get("content_level", default="", decode=True).split("#")

    @property
    def res_description(self):
        """
        the textual description of the resource.

        """
        return self.get("res_description", decode=True)

    @property
    def reference_url(self):
        """
        URL pointing to a human-readable document describing this resource.
        """
        return self.get("reference_url", decode=True)

    @property
    def creators(self):
        """
        The creator(s) of the resource
        in the order given by the resource record author
        """
        return self.get("creator_seq", default="", decode=True).split(";")

    @property
    def created(self):
        """Date of creation of the resource."""
        return self.get("created", decode=True)

    @property
    def updated(self):
        """Date of last modification of the resource."""
        return self.get("updated", decode=True)

    @property
    def rights(self):
        """A statement of usage conditions for the content of the resource.

        This information is often incomplete in the registry, you
        might get more information at the ``reference_url``.
        """
        return self.get("rights", decode=True)

    @property
    def content_types(self):
        """
        list of natures or genres of the content of the resource.
        """
        return self.get("content_type", decode=True).split("#")

    @property
    def source_format(self):
        """
        The format of source_value.
        """
        return self.get("source_format", decode=True)

    @property
    def source_value(self):
        """
        The bibliographic source for this resource (typically a bibcode
        or a DOI).
        """
        return self.get("source_value", decode=True)

    @property
    def region_of_regard(self):
        """
        numeric value representing the angle, given in decimal degrees,
        by which a positional query against this resource should be "blurred"
        in order to get an appropriate match.
        """
        # we get NULLs as NaNs here
        val = self["region_of_regard"]
        if numpy.isnan(val):
            return None
        return val

    @property
    def waveband(self):
        """
        a list of names of the wavebands that the resource provides data for
        """
        return self.get("waveband", default="", decode=True).split("#")

    @property
    def access_url(self):
        """
        the URL that can be used to access the service resource.
        """
        # some services declare some data models using multiple
        # identifiers; in this case, we'll get the same access URL
        # multiple times in here.  Be cool about that situation:
        access_urls = list(sorted(set(self["access_urls"])))

        if len(access_urls) == 0:
            raise dalq.DALQueryError(
                f"The resource {self.ivoid} has no queriable interfaces.")

        elif len(access_urls) > 1:
            warnings.warn(AstropyDeprecationWarning(
                f"The resource {self.ivoid} has multiple capabilities. "
                " You should explicitly pick one using get_service. "
                " Returning some access_url now, but this behaviour "
                " may change in the future."))
        return access_urls[0]

    @property
    def standard_id(self):
        """
        the IVOA standard identifier
        """
        standard_ids = list(set(self["standard_ids"]))
        if len(standard_ids) == 1:
            return standard_ids[0]
        else:
            raise dalq.DALQueryError(
                "This resource supports several standards ({})."
                "  Use get_service or restrict your query using Servicetype."
                .format(", ".join(sorted(self.access_modes()))))

    def access_modes(self):
        """
        returns a set of interface identifiers available on
        this resource.

        For standard interfaces, get_service will return a service
        suitable for querying if you pass in an identifier from this
        list as the service_type.

        This will ignore VOSI (infrastructure) services.
        """
        return set(shorten_stdid(intf.standard_id) or "web"
                   for intf in self.interfaces
                   if (intf.standard_id or intf.type == "vr:webbrowser")
                   and not intf.is_vosi)

    def get_interface(self, *,
                      service_type: str,
                      lax: bool = False,
                      std_only: bool = False,
                      keyword: str = None):
        """returns a regtap.Interface class for service_type.

        The meaning of the parameters is as for get_service.  This
        method does not return services, though, so you can use it to
        obtain access URLs and such for interfaces that pyVO does
        not (directly) support.

        Parameters
        ----------

        service_type : str
            If you leave out ``service_type``, this will return a service
            for "the" standard interface of the resource.  If a resource
            has multiple standard capabilities (e.g., both TAP and SSAP
            endpoints), this will raise a DALQueryError.

            Otherwise, a service of the given service type will be returned.
            Pass in an ivoid of a standard or one of the shorthands from
            rtcons.SERVICE_TYPE_MAP, or "web" for a web page (the "service"
            for this will be an object opening a web browser when you call
            its query method).

        lax : bool
            If there are multiple capabilities for service_type, the
            function choose the first matching capability by default
            Pass lax=False to instead raise a DALQueryError.

        std_only : bool
            Only return interfaces declared as "std".  This is what you
            want when you want to construct pyVO service objects later.
            This parameter is ignored for the "web" service type.


        Returns
        -------

        `~pyvo.registry.regtap.Interface`
        """
        if service_type == "web":
            # this works very much differently in the Registry
            # than do the proper services
            candidates = [intf for intf in self.interfaces
                          if intf.type == "vr:webbrowser"]

        else:
            service_type = expand_stdid(
                rtcons.SERVICE_TYPE_MAP.get(
                    service_type, service_type))

            candidates = [intf for intf in self.interfaces
                          if ((not std_only) or intf.is_standard)
                          and not intf.is_vosi
                          and ((not service_type) or intf.supports(service_type))]
        if not candidates:
            raise ValueError("No matching interface.")

        if keyword is not None:
            candidates = [candidate for candidate in candidates
                          if candidate.capability_description is not None
                          and keyword in candidate.capability_description]

        if len(candidates) > 1 and not lax:
            raise ValueError("Multiple matching interfaces found. If you know that they are all equivalent, "
                             "use 'lax=True'.\n"
                             "Otherwise, the search can be constrained by specifying a 'service_type' "
                             "or a 'keyword' to be found in the service description."
                             " You might also want to see all the available services"
                             " with `pyvo.registry.regtap.RegistryResource.list_interfaces()`.")

        return candidates[0]

    def get_service(self, service_type: str = None, *,
                    lax: bool = False,
                    keyword: str = None,
                    session: object = None):
        """
        return an appropriate DALService subclass for this resource that
        can be used to search the resource using service_type.

        Raise a ValueError if the service_type is not offered on
        the resource (or no standard service is offered).  With
        lax=False, also raise a ValueError if multiple interfaces
        exist for the given service_type.

        VOSI (infrastructure) services are always ignored here.

        A magic service_type "web" can be passed in to get non-standard,
        browser-based interfaces.  The service in this case is an
        object that opens a web browser if its query() method is called.

        Parameters
        ----------
        service_type : str
            If you leave out ``service_type``, this will return a service
            for "the" standard interface of the resource.  If a resource
            has multiple standard capabilities (e.g., both TAP and SSAP
            endpoints), this will raise a DALQueryError.

            Otherwise, a service of the given service type will be returned.
            Pass in an ivoid of a standard or one of the shorthands from
            rtcons.SERVICE_TYPE_MAP, or "web" for a web page (the "service"
            for this will be an object opening a web browser when you call
            its query method).

        lax : bool
            If there are multiple capabilities for service_type, the
            function choose the first matching capability by default
            Pass lax=False to instead raise a DALQueryError.

        keyword : str
            A keyword that should be in the service description.
            Some resources have multiple capabilities ("services") of the same
            type (see list_interfaces for a discussion).  get_service will
            usually raise a ValueError in that case. Passing a
            ``keyword`` that uniquely identifies the capability to
            query will make get_service predictably return the desired
            service.  Use list_interfaces to find such a unique description
            fragment.

        session : object
            optional requests session to use to communicate with the service
            constructed.

        Returns
        -------
        `pyvo.dal.DALService`
            For standard service types, a specific DAL service instance
            (e.g., a `pyvo.dal.tap.TAPService` when requesting ``tap``
            services) is returned.  For ``web`` services, what is returned is
            an opaque service object that has a ``search()`` method simply
            opening a web browser on the access URL.

        See Also
        --------
        list_interfaces : return a list with all the available services.
        """
        return self.get_interface(service_type=service_type, lax=lax, std_only=True,
                                  keyword=keyword).to_service(session=session)

    @property
    def service(self):
        """
        return a service for this resource.

        This will in general only work if the registry query has
        constrained the service type; otherwise, many resources will
        have multiple capabilities.  Use get_service instead in
        such cases.

        See Also
        --------
        get_service : select a service of a specific service type.
        list_interfaces : return a list with all the available interfaces to services.
        """
        if self._service is not None:
            return self._service
        self._service = self.get_service(service_type=None, lax=True)
        return self._service

    def list_interfaces(self, service_type: str = None):
        """List the interfaces to services available for this registry record.

        The interface objects in the list can then be used to instantiate the service with
        the  `~pyvo.registry.regtap.Interface.to_service` method.

        Parameters
        ----------
        service_type : str, optional
            Restricts the list to interfaces corresponding to services of this service type.
            The list of possible values can be printed with
            ``pyvo.registry.rtcons.SERVICE_TYPE_MAP.keys()``.

        Returns
        -------
        list
            A list of `~pyvo.registry.regtap.Interface` objects.

        See Also
        --------
        get_service : when you already know that there is only one service of a specific service type.

        """
        interfaces = [interface for interface in self.interfaces
                      if interface.is_user_visible]

        if service_type is not None:
            service_type = expand_stdid(rtcons.SERVICE_TYPE_MAP.get(service_type, service_type))
            interfaces = [interface for interface in interfaces
                          if interface.is_standard and interface.supports(service_type)]

        return sorted(interfaces, key=lambda interface: interface.access_url)

    def search(self, *args, **keys):
        """
        assuming this resource refers to a searchable service, execute a
        search against the resource.  This is equivalent to:

        .. code:: python

           self.to_service().search(*args, **keys)

        The arguments provided should be appropriate for the service that
        the DAL service type would expect.  See the documentation for the
        appropriate service type:

        ============  ===========================================
        Service type  Use the argument syntax for
        ============  ===========================================
        catalog       :py:meth:`pyvo.dal.scs.SCSService.search`
        sia           :py:meth:`pyvo.dal.sia.SIAService.search`
        sia2          :py:meth:`pyvo.dal.sia2.SIA2Service.search`
        ssa           :py:meth:`pyvo.dal.ssa.SSAService.search`
        line          :py:meth:`pyvo.dal.sla.SLAService.search`
        database      *not yet supported*
        ============  ===========================================

        Raises
        ------
        DALServiceError
           if the resource does not describe a searchable service.
        """
        try:
            return self.service.search(*args, **keys)
        except ValueError:
            # I blindly assume the ValueError comes out of get_interface.
            # But then that's likely enough.
            raise dalq.DALServiceError(
                f"Resource {self.ivoid} is not a searchable service")

    def describe(self, *, verbose=False, width=78, file=None):
        """
        Print a summary description of this resource.

        Parameters
        ----------
        verbose : bool
            If false (default), only user-oriented information is
            printed.
            If true, additional information -- reference url,
            reference to the related article, and alternative identifier
            (often a DOI) -- will be printed if available.
        width : int
            Format the description with given character-width.
        file : writable file-like object
            If provided, write information to this output stream.
            Otherwise, it is written to standard out.
        """
        print(para_format_desc(self.res_title), file=file)
        print("Short Name: " + self.short_name, file=file)
        print("IVOA Identifier: " + self.ivoid, file=file)
        print("Access modes: " + ", ".join(sorted(self.access_modes())),
              file=file)
        for interface in self.list_interfaces():
            if interface.is_user_visible:
                if interface.type == "vr:webbrowser":
                    print(f"- webpage: {interface.access_url}", file=file)
                else:
                    service_item = (f"- {interface.standard_id.rsplit('/')[-1]}: "
                                   f"{interface.access_url}")
                    if interface.capability_description:
                        service_item += f", description: {interface.capability_description}"
                    print(textwrap.fill(service_item, width, subsequent_indent=" "), file=file)

        if self.res_description:
            print(file=file)
            print(para_format_desc(self.res_description), file=file)
            print(file=file)

        if self.waveband:
            val = (str(v) for v in self.waveband)
            print(
                para_format_desc("Waveband Coverage: " + ", ".join(val)),
                file=file)

        if verbose:
            if self.source_value:
                print(f"\nSource: {self.source_value}", file=file)
            if self.creators:
                # if any creator has a name longer than 70 characters, we
                # truncate it.
                creators = [f"{creator[:70]}..." if len(creator) > 70
                            else creator for creator in self.creators]
                nmax_authors = 5
                if len(creators) <= nmax_authors:
                    print(f"Authors: {', '.join(creators)}", file=file)
                else:
                    print(f"Authors: {', '.join(creators[:nmax_authors])} et al.\n"
                    "See creators attribute for the complete list of authors.", file=file)

            alt_identifiers = self.get_alt_identifiers()
            if alt_identifiers:
                print(
                    "Alternative identifier(s): {}".format(
                        ", ".join(alt_identifiers)),
                    file=file)

            if self.reference_url:
                print("More info: " + self.reference_url, file=file)

    def get_contact(self):
        """
        return contact information for this resource in a string.

        Use this to report bugs or unexpected downtime.
        """
        res = get_RegTAP_service().run_sync("""
            SELECT role_name, email, telephone
            FROM rr.res_role
            WHERE
                base_role='contact'
                AND ivoid={}""".format(
            rtcons.make_sql_literal(self.ivoid)))

        contacts = []
        for row in res:
            contact = row["role_name"]
            if row["telephone"]:
                contact += f" ({row['telephone']})"
            if row["email"]:
                contact += f" <{row['email']}>"
            contacts.append(contact)

        return "\n".join(contacts)

    def get_alt_identifiers(self):
        """return a sequence of non-ivoid identifiers for the resource.

        This is typically used to provide a DOI for the resource.
        """
        res = get_RegTAP_service().run_sync("""
            SELECT alt_identifier
            FROM rr.alt_identifier
            WHERE ivoid={}""".format(rtcons.make_sql_literal(self.ivoid)))
        return [r["alt_identifier"] for r in res]

    def _build_vosi_column(self, column_row):
        """
        return a io.vosi.vodataservice.Column element for a
        query result from get_tables.
        """
        res = vodataservice.TableParam()
        for att_name in ["name", "ucd", "unit", "utype"]:
            setattr(res, att_name, column_row[att_name])
        res.description = column_row["column_description"]

# TODO: be more careful with the type; this isn't necessarily a
# VOTable type (regrettably)
        res.datatype = vodataservice.VOTableType(
            arraysize=column_row["arraysize"],
            extendedType=column_row["extended_type"])
        res.datatype.content = column_row["datatype"]

        return res

    def _build_vosi_table(self, table_row, columns):
        """
        return a io.vosi.vodataservice.VODataServiceTable element for a
        query result from get_tables.
        """
        res = vodataservice.VODataServiceTable()
        res.name = table_row["table_name"]
        res.title = table_row["table_title"]
        res.description = table_row["table_description"]
        res.utype = table_row["table_utype"]
        res._columns = [
            self._build_vosi_column(row)
            for row in columns]

        res.origin = self

        return res

    def get_tables(self, *, table_limit=20):
        """
        return the structure of the tables underlying the service.

        This returns a dict with table names as keys and vodataservice.VODataServiceTable
        objects as values (pretty much what tables returns for a TAP
        service).  The table instances will have an ``origin`` attribute
        pointing back to the registry record.

        Note that not only TAP services can (and do) define table
        structures.  The meaning of non-TAP tables is not always
        as clear.

        Also note that resources do not need to define tables at all.
        You will receive an empty dictionary if they don't.
        """
        svc = get_RegTAP_service()

        tables = svc.run_sync(
            """SELECT table_name, table_description, table_index, table_title,
                table_utype
            FROM rr.res_table
            WHERE ivoid={}""".format(
                rtcons.make_sql_literal(self.ivoid)))
        if len(tables) > table_limit:
            raise dalq.DALQueryError(f"Resource {self.ivoid} reports"
                                     f" {len(tables)} tables.  Pass a higher table_limit"
                                     " to see them all.")

        res = {}
        for table_row in tables:
            columns = svc.run_sync(
                """
                SELECT name, ucd, unit, utype, datatype, arraysize,
                    extended_type, column_description
                FROM rr.table_column
                WHERE ivoid={}
                    AND table_index={}""".format(
                    rtcons.make_sql_literal(self.ivoid),
                    rtcons.make_sql_literal(table_row["table_index"])))
            res[table_row["table_name"]] = self._build_vosi_table(
                table_row, columns)

        return res


@deprecated("1.5", "ivoid2service does not work in the presence of"
    " multiple capabilities.  Use"
    " registry.search(ivoid=...)[0].get_service('capname') instead.")
def ivoid2service(ivoid, servicetype=None):
    """
    return service(s) for a given IVOID.

    The servicetype option specifies the kind of service requested
    (conesearch, sia, sia2, ssa, slap, or tap).  By default, if none is
    given, a list of all matching services is returned.

    """
    constraints = [rtcons.Ivoid(ivoid)]
    if servicetype is not None:
        constraints.append(rtcons.Servicetype(servicetype))
    resources = search(*constraints)
    if len(resources) == 0:
        if servicetype:
            raise dalq.DALQueryError(f"No resource {ivoid} with"
                                     f" {servicetype} capability.")
        else:
            raise dalq.DALQueryError(f"No resource {ivoid}")

    # We're grouping by ivoid in search, so if there's a result
    # there is only one.
    resource = resources[0]

    return resource.get_service(servicetype, lax=True)
