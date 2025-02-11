# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Constraints for doing registry searches.

The Constraint class encapsulates a query fragment in a RegTAP query, e.g., a
keyword, a sky location, an author name, a class of services.  They are used
either directly as arguments to registry.search, or by passing keyword
arguments into registry.search.  The mapping from keyword arguments to
constraint classes happens through the _keyword attribute in Constraint-derived
classes.
"""

import datetime
import warnings

from astropy import units as u
from astropy import constants
from astropy.coordinates import SkyCoord
from astropy.utils.exceptions import AstropyDeprecationWarning
import numpy

from ..dal import query as dalq
from ..utils import vocabularies
from .import regtap


# Classes from this module are exposed at the higher level namespace, not listing them here
__all__ = ["build_regtap_query"]


# a mapping of service type shorthands to the ivoids of the
# corresponding standards.  This is mostly to keep legacy APIs.
# In the future, preferably rely on shorten_stdid and expand_stdid
# from regtap.
SERVICE_TYPE_MAP = dict((k, "ivo://ivoa.net/std/" + v)
                        for k, v in [
    ("image", "sia"),
    ("sia", "sia"),
    ("sia1", "sia"),
    # SIA2 is irregular
    # funky scheme used by SIA2 without breaking everything else
    ("spectrum", "ssa"),
    ("ssap", "ssa"),
    ("ssa", "ssa"),
    ("scs", "conesearch"),
    ("conesearch", "conesearch"),
    ("line", "slap"),
    ("slap", "slap"),
    ("table", "tap"),
    ("tap", "tap"),
])


class RegTAPFeatureMissing(dalq.DALQueryError):
    """
    Raised when the current RegTAP server does not support a feature
    needed for a constraint.

    This could be that it is missing some ADQL feature indispensible
    to write that constraint, or because it is missing a table or column.

    To recover, choose another RegTAP service. Search constraining
    ``datamodel="regtap"``, and then use `pyvo.registry.choose_RegTAP_service`
    with a TAP access URL discovered in this way.
    """


class _AsIs(str):
    """a sentinel class make `make_sql_literal` not escape a string.
    """


def make_sql_literal(value):
    """makes a SQL literal from a python value.

    This is not suitable as a device to ward against SQL injections;
    in what we produce, callers could produce arbitrary SQL anyway.
    The point of this function is to minimize surprises when building
    constraints.

    Parameters
    ----------

    value : object
        Conceptually, the function should produces SQL literals
        for anything that might reasonably add up in a registry
        query.  In reality, a ValueError will be raised for anything
        we do not know about.

    Returns
    -------

    str
        A SQL literal.
    """
    if isinstance(value, _AsIs):
        return value

    if isinstance(value, str):
        return "'{}'".format(value.replace("'", "''"))

    elif isinstance(value, bytes):
        return "'{}'".format(value.decode("ascii").replace("'", "''"))

    elif isinstance(value, (int, numpy.integer)):
        return "{:d}".format(value)

    elif isinstance(value, (float, numpy.floating)):
        return repr(value)

    elif isinstance(value, datetime.datetime):
        return "'{}'".format(value.isoformat())

    else:
        raise ValueError("Cannot format {} as a SQL literal"
                         .format(repr(value)))


def format_function_call(func_name, args):
    """make an ADQL literal for a function call with arguments.

    Parameters
    ----------
    func_name : str
        the name of the function to call.

    args : sequence of anything
        python values for the arguments for the function.

    Returns
    -------
    str
        ADQL ready for inclusion into a query.
    """
    return "{}({})".format(
        func_name,
        ", ".join(make_sql_literal(a) for a in args))


class Constraint:
    """an abstract base class for data discovery contraints.

    These, essentially, are configurable RegTAP query fragments,
    consisting of a where clause, parameters for filling that,
    and possibly additional tables.

    Users construct concrete constraints with whatever they would like
    to constrain things with.

    To implement a new constraint, in the constructor set ``_condition`` to a
    string with {}-type replacement fields (assume all parameters are strings),
    and define ``fillers`` to be a dictionary with values for the _condition
    template.  Don't worry about SQL-serialising the values, Constraint takes
    care of that.  If you need your Constraint to be "lazy"
    (cf. Servicetype), it's ok to overrride get_search_condition without
    an upcall to Constraint.

    If your constraints need extra tables, give them in a list
    in _extra_tables.

    For the legacy x_search with keywords, define a _keyword
    attribute containing the name of the parameter that should
    generate such a constraint.  When pickung up such keywords,
    sequence values will in general be unpacked and turned into
    sequences of constraints.  Constraints that want to the all
    arguments in the constructor can set takes_sequence to True.
    """
    # TODO: _extra_tables is only used in the legacy leg of
    # the fulltext constraint any more, and there it's wrong, too.
    # Let's do away with extra_tables and tell people to use
    # SubqueriedConstraint whenever they think they need it.
    _extra_tables = []
    _condition = None
    _fillers = None
    _keyword = None

    takes_sequence = False

    def get_search_condition(self, service):
        """
        Formats this constraint to an ADQL fragment.

        This takes the service the constraint is being executed on
        as an argument because constraints may be written differently
        depending on the service's features or refuse to run altogether.

        Parameters
        ----------
        service : `~pyvo.dal.TAPService`
            The RegTAP service the query is supposed to be run on
            (that is relevant because we adapt to the features available
            on given services).

        Returns
        -------
        str
            A string ready for inclusion into a WHERE clause.
        """
        if self._condition is None:
            raise NotImplementedError("{} is an abstract Constraint"
                                      .format(self.__class__.__name__))

        return self._condition.format(**self._get_sql_literals())

    def _get_sql_literals(self):
        """
        returns self._fillers as a dictionary of properly SQL-escaped
        literals.
        """
        if self._fillers:
            return {k: make_sql_literal(v) for k, v in self._fillers.items()}
        return {}


class SubqueriedConstraint(Constraint):
    """
    An (abstract) constraint for when the constraint is over a table
    other than rr.resource.

    We need to be careful with these, because they will in general have
    1:n relationships to rr.resource, and these will lead to
    duplicated interfaces if we just do the NATURAL JOIN we normally
    do.

    Instead, we have to resort to ivoid in (subquery) conditions.  In
    particular, extra_tables will always be empty for these.

    To configure those, give the table to query in _subquery_table
    and _condition, _fillers, and _keyword as usual.  The rest is taken
    care of by get_search_condition.
    """
    _subquery_table = None

    def get_search_condition(self, service):
        if self._condition is None or self._subquery_table is None:
            raise NotImplementedError("{} is an abstract Constraint"
                                      .format(self.__class__.__name__))

        return ("ivoid IN (SELECT DISTINCT ivoid FROM {subquery_table}"
            " WHERE {condition})".format(
                subquery_table=self._subquery_table,
                condition=self._condition.format(**self._get_sql_literals())))


class Freetext(Constraint):
    """
    A constraint using plain text to match against title, description,
    subjects, and person names.
    """
    _keyword = "keywords"

    def __init__(self, *words: str):
        """

        Parameters
        ----------
        *words : str
            One or more string arguments.
            It is recommended to pass multiple words in multiple strings
            arguments.  You can pass in phrases (i.e., multiple words
            separated by space), but behaviour might then vary quite
            significantly between different registries.

        Examples
        --------
        >>> from pyvo import registry
        >>> registry.Freetext("Gamma", "Ray", "Burst") # doctest: +IGNORE_OUTPUT
        """
        self.words = words

    def get_search_condition(self, service):
        # cross-table ORs kill the query planner.  We therefore
        # write the constraint as an IN condition on a UNION
        # of subqueries if we can (i.e., the service has UNION);
        # It may look as if this has to be really slow, but in fact it's almost
        # always a lot faster than direct ORs.
        if service.get_tap_capability().get_adql().get_feature(
                "ivo://ivoa.net/std/TAPRegExt#features-adql-sets", "UNION"):
            return self._get_union_condition(service)
        else:
            self._extra_tables = ["rr.res_subject"]
            return self._get_or_condition(service)

    def _get_union_condition(self, service):
        base_queries = [
            "SELECT DISTINCT ivoid FROM rr.resource WHERE"
            " 1=ivo_hasword(res_description, {{{parname}}})",
            "SELECT DISTINCT ivoid FROM rr.resource WHERE"
            " 1=ivo_hasword(res_title, {{{parname}}})",
            "SELECT DISTINCT ivoid FROM rr.res_subject WHERE"
            " rr.res_subject.res_subject ILIKE {{{parpatname}}}"]
        self._fillers, subqueries = {}, []

        for index, word in enumerate(self.words):
            parname = "fulltext{}".format(index)
            parpatname = "fulltextpar{}".format(index)
            self._fillers[parname] = word
            self._fillers[parpatname] = '%' + word + '%'
            args = locals()
            subqueries.append(" UNION ALL ".join(
                q.format(**args) for q in base_queries))

        self._condition = " AND ".join(
            f"ivoid IN ({part})" for part in subqueries)

        return super().get_search_condition(service)

    def _get_or_condition(self, service):
        base_queries = [
            " 1=ivo_hasword(res_description, {{{parname}}})",
            " 1=ivo_hasword(res_title, {{{parname}}})",
            " rr.res_subject.res_subject ILIKE {{{parpatname}}}"]
        self._fillers, conditions = {}, []

        for index, word in enumerate(self.words):
            parname = "fulltext{}".format(index)
            parpatname = "fulltextpar{}".format(index)
            self._fillers[parname] = word
            self._fillers[parpatname] = '%' + word + '%'
            args = locals()
            conditions.append(" OR ".join(
                q.format(**args) for q in base_queries))

        self._condition = " AND ".join(f"({part})"
            for part in conditions)

        return super().get_search_condition(service)


class Author(SubqueriedConstraint):
    """
    A constraint for creators (“authors”) of a resource; you can use SQL
    patterns here.

    The match is case-sensitive.
    """
    _keyword = "author"
    _subquery_table = "rr.res_role"

    def __init__(self, name: str):
        """

        Parameters
        ----------
        name : str
            Note that regrettably there are no guarantees as to how authors
            are written in the VO.  This means that you will generally have
            to write things like ``%Hubble%`` (% being “zero or more
            characters” in SQL) here.
        """
        self._condition = "role_name LIKE {auth} AND base_role='creator'"
        self._fillers = {"auth": name}


class Servicetype(Constraint):
    """
    A constraint for for the availability of a certain kind of service
    on the result.

    The constraint normally is a custom keyword, one of:


    * ``sia``, ``sia1`` (SIAP version 1 services; prefer ``sia1`` for symmetry,
      although ``sia`` will be kept as the official IVOA short name for SIA1)
    * ``sia2`` (SIAP version 2 services)
    * ``ssa``, ``ssap`` (synonymous for SSAP services)
    * ``scs``, ``conesearch`` (synonymous for cone search services, prefer
      ``scs``)
    * ``line`` (for SLAP services)
    * ``tap``, ``table`` (synonymous for TAP services, prefer ``tap``)
    * ``image`` (a deprecated alias for sia)
    * ``spectrum`` (a deprecated alias for ssap)

    You can also pass in the standards' ivoid (which
    generally looks like
    ``ivo://ivoa.net/std/<standardname>`` (except for SIA2) and have
    to be URIs with a scheme part in any case); note, however, that for
    standards pyVO does not know about it will not build service instances
    for you.

    Multiple service types can be passed in; a match in that case
    is for records having any of the service types passed in.

    Contrary to what the names of the service types may suggest,
    “image” and "spectrum" do *not* select all resources serving images
    or spectra.  For instance, at the moment you would have to search for
    images served in three different ways: SIAv1, SIAv2, and ObsTAP. Against
    that, “image” only selects SIAv1, and "spectrum" similarly omits ObsTAP.  A
    global discovery module is under active development to provide global
    dataset discovery. When it is ready, these two misleading options may be
    removed.

    The match is literal (i.e., no patterns are allowed); this means
    that you will not receive records that only have auxiliary
    services, which is what you want when enumerating all services
    of a certain type in the VO.  In data discovery, you
    can use ``Servicetype(...).include_auxiliary_services()`` or
    use registry.search's ``includeaux`` parameter; but, really, there
    is little point using this constraint in data discovery in the first
    place.
    """
    _keyword = "servicetype"

    def __init__(self, *stds):
        """

        Parameters
        ----------
        *stds : str
            one or more standards identifiers.  The constraint will
            match records that have any of them.
        """
        self.stdids = set()
        self.extra_fragments = []

        for std in stds:
            if std in ('image', 'spectrum'):
                warnings.warn(AstropyDeprecationWarning(
                    f"The '{std}' servicetype is deprecated.  See"
                    " https://pyvo.readthedocs.io/en/latest/api/pyvo.registry"
                    ".Servicetype.html#pyvo.registry.Servicetype"
                    " for more information."))

            if std in SERVICE_TYPE_MAP:
                self.stdids.add(SERVICE_TYPE_MAP[std])
            elif "://" in std:
                self.stdids.add(std)
            elif std == 'sia2':
                self.extra_fragments.append(
                    "standard_id like 'ivo://ivoa.net/std/sia#query-2.%'")
            else:
                raise dalq.DALQueryError("Service type {} is neither a full"
                                         " standard URI nor one of the bespoke identifiers"
                                         " {}, sia2".format(std, ", ".join(SERVICE_TYPE_MAP)))

    def clone(self):
        """returns a copy of this servicetype constraint.
        """
        new_constraint = Servicetype()
        new_constraint.stdids = set(self.stdids)
        new_constraint.extra_fragments = self.extra_fragments[:]
        return new_constraint

    def get_search_condition(self, service):
        # we sort the stdids to make it easy for tests (and it's
        # virtually free for the small sets we have here).
        fragments = []
        std_ids = ", ".join(make_sql_literal(s)
            for s in sorted(self.stdids))
        if std_ids:
            fragments.append(f"standard_id IN ({std_ids})")

        return " OR ".join(fragments + self.extra_fragments)

    def include_auxiliary_services(self):
        """returns a Servicetype constraint that has self's
        service types but includes the associated auxiliary services.

        This is a convenience to maintain registry.search's signature.
        """
        expanded = self.clone()
        expanded.stdids |= set(
            std + '#aux' for std in expanded.stdids)
        if "standard_id like 'ivo://ivoa.net/std/sia#query-2.%'" in expanded.extra_fragments:
            expanded.extra_fragments.append(
                "standard_id like 'ivo://ivoa.net/std/sia#query-aux-2.%'")
        return expanded


class Waveband(Constraint):
    """
    A constraint on messenger particles.

    This builds a constraint against rr.resource.waveband, i.e.,
    a verbal indication of the messenger particle, coming
    from the IVOA vocabulary http://www.ivoa.net/rdf/messenger.

    The :py:class:`pyvo.registry.Spectral` constraint enables selections by particle energy,
    but few resources actually give the necessary metadata (in 2021).

    Multiple wavebands can be given (and are effectively combined with OR).
    """
    _keyword = "waveband"
    _legal_terms = None

    def __init__(self, *bands):
        """

        Parameters
        ----------
        *bands : strings
            One or more of the terms given in http://www.ivoa.net/rdf/messenger.
            The constraint matches when a resource declares at least
            one of the messengers listed.
        """
        if self.__class__._legal_terms is None:
            self.__class__._legal_terms = {w.lower() for w in
                                           vocabularies.get_vocabulary("messenger")["terms"]}

        bands = [band.lower() for band in bands]
        for band in bands:
            if band not in self._legal_terms:
                raise dalq.DALQueryError(
                    f"Waveband {band} is not in the IVOA messenger"
                    " vocabulary http://www.ivoa.net/rdf/messenger.")

        self.bands = list(bands)
        self._condition = " OR ".join(
            "1 = ivo_hashlist_has(rr.resource.waveband, {})".format(
                make_sql_literal(band))
            for band in self.bands)


class Datamodel(SubqueriedConstraint):
    """
    A constraint on the adherence to a data model.

    This constraint only lets resources pass that declare support for
    one of several well-known data models; the SQL produced depends
    on the data model identifier.

    Records returned with these constraints will have a TAP (or
    auxiliary TAP) capability.

    Known data models at this point include:

    * obscore -- generic observational data
    * epntap -- solar system data
    * regtap -- the VO registry.
    * obscore-new -- the table-based new-style obscore discovery.  Don't
      use in code built to last: this will become normal obscore
      when we have migrated the VO.

    DM names are matched case-insensitively here mainly for
    historical reasons.
    """
    _keyword = "datamodel"

    # if you add to this list, you have to define a method
    # _make_<dmname>_constraint.
    _known_dms = {"obscore", "epntap", "regtap", "obscore_new"}

    def __init__(self, dmname):
        """

        Parameters
        ----------
        dmname : string
            A well-known name; currently one of obscore, epntap, and regtap.
        """
        dmname = dmname.lower()
        if dmname not in self._known_dms:
            raise dalq.DALQueryError("Unknown data model id {}.  Known are: {}."
                                     .format(dmname, ", ".join(sorted(self._known_dms))))
        self._subquery_table, self._condition = getattr(
            self, f"_make_{dmname}_constraint")()

    def _make_obscore_constraint(self):
        # There was a bit of chaos with the DM ids for Obscore.
        # Be lenient here
        obscore_pat = 'ivo://ivoa.net/std/obscore%'
        return "rr.res_detail", (
            "detail_xpath = '/capability/dataModel/@ivo-id'"
            f" AND 1 = ivo_nocasematch(detail_value, '{obscore_pat}')")

    def _make_obscore_new_constraint(self):
        return "rr.res_table NATURAL JOIN rr.resource", (
            "table_utype LIKE 'ivo://ivoa.net/std/obscore#table-1.%'"
            # Only use catalogresource-typed records to keep out
            # full TAP services that may have the table in their
            # tablesets.
            " AND res_type = 'vs:catalogresource'")

    def _make_epntap_constraint(self):
        # we include legacy, pre-IVOA utypes for matches; lowercase
        # any new identifiers (utypes case-fold).
        return "rr.res_table", " OR ".join(
            f"table_utype LIKE '{pat}'" for pat in
            ['ivo://vopdc.obspm/std/epncore#schema-2.%',
             'ivo://ivoa.net/std/epntap#table-2.%'])

    def _make_regtap_constraint(self):
        regtap_pat = 'ivo://ivoa.net/std/RegTAP#1.%'
        return "rr.res_detail", (
            "detail_xpath = '/capability/dataModel/@ivo-id'"
            f" AND 1 = ivo_nocasematch(detail_value, '{regtap_pat}')")


class Ivoid(Constraint):
    """
    A constraint selecting a single resource by its IVOA identifier.
    """
    _keyword = "ivoid"

    def __init__(self, ivoid, *more_ivoids):
        """

        Parameters
        ----------

        ivoid : string
            One or more string arguments.
            The IVOA identifier of the resource to match.  As RegTAP
            requires lowercasing ivoids on ingestion, the constraint
            lowercases the ivoid passed in, too.

        more_ivoids : strings
            You can pass in multiple ivoids to match.  As usual,
            they are combined by an or.
        """
        self.ivoids = [id.lower() for id in (ivoid,) + more_ivoids]

    def get_search_condition(self, service):
        return " OR ".join(
            f"ivoid={make_sql_literal(id)}" for id in self.ivoids)


class UCD(SubqueriedConstraint):
    """
    A constraint selecting resources having tables with columns having
    UCDs matching a SQL pattern (% as wildcard).
    """
    _keyword = "ucd"
    _subquery_table = "rr.table_column"

    def __init__(self, *patterns):
        """

        Parameters
        ----------

        patterns : strings
            SQL patterns (i.e., ``%`` is 0 or more characters) for
            UCDs.  The constraint will match when a resource has
            at least one column matching one of the patterns.
        """
        self._condition = " OR ".join(
            f"ucd LIKE {{ucd{i}}}" for i in range(len(patterns)))
        self._fillers = dict((f"ucd{index}", pattern)
                             for index, pattern in enumerate(patterns))


class Spatial(SubqueriedConstraint):
    """
    A RegTAP constraint selecting resources covering a geometry in
    space.

    This is a RegTAP 1.2 extension not yet available on all Registries
    (in 2022).  Also note that not all data providers give spatial coverage
    for their resources.

    To find resources having data for RA/Dec 347.38/8.6772::

        >>> from pyvo import registry
        >>> resources = registry.Spatial((347.38, 8.6772))

    To find resources claiming to have data for a spherical circle 2 degrees
    around that point::

        >>> resources = registry.Spatial((347.38, 8.6772, 2))

    To find resources claiming to have data for a polygon described by
    the vertices (23, -40), (26, -39), (25, -43) in ICRS RA/Dec::

        >>> resources = registry.Spatial([23, -40, 26, -39, 25, -43])

    To find resources claiming to cover a MOC_, pass an ASCII MOC::

        >>> resources = registry.Spatial("0/1-3 3/")

    .. _MOC: https://www.ivoa.net/documents/MOC/

    To find resources which coverage is enclosed in a region,

        >>> enclosed = registry.Spatial("0/0-11", intersect="enclosed")

    To find resources which coverage intersects a region,

        >>> overlaps = registry.Spatial("0/0-11", intersect="overlaps")

    When you already have an astropy SkyCoord::

        >>> from astropy.coordinates import SkyCoord
        >>> resources = registry.Spatial(SkyCoord("23d +3d"))

    SkyCoords also work as circle centers (plain floats for the radius
    are interpreted in degrees)::

        >>> resources = registry.Spatial((SkyCoord("23d +3d"), 3))

    Or you can provide the radius angle as an Astropy Quantity:

        >>> resources = registry.Spatial((SkyCoord("23d +3d"), 1*u.rad))
    """
    _keyword = "spatial"
    _subquery_table = "rr.stc_spatial"

    takes_sequence = True

    def __init__(self, geom_spec,
            *, order=6, intersect="covers", inclusive=False):
        """

        Parameters
        ----------
        geom_spec : object
            For now, this is DALI-style: a 2-sequence is interpreted
            as a DALI point, a 3-sequence as a DALI circle, a 2n sequence
            as a DALI polygon.  Additionally, strings are interpreted
            as ASCII MOCs, SkyCoords as points, and a pair of a
            SkyCoord and a float or Quantity as a circle.  Other types (proper
            geometries or MOCPy objects) might be supported in the
            future.
        order : int, optional
            Non-MOC geometries are converted to MOCs before comparing
            them to the resource coverage.  By default, this constraint
            uses order 6, which corresponds to about a degree of resolution
            and is what RegTAP recommends as a sane default for the
            order actually used for the coverages in the database.
        intersect : str, optional
            Allows to specify the connection between the resource coverage
            and the *geom_spec*. The possible values are 'covers' for services
            that completely cover the *geom_spec* region, 'enclosed' for services
            completely enclosed in the region and 'overlaps' for services which
            coverage intersect the region.
        inclusive : bool, optional
            Normally, this constraint will remove all resources that do
            not declare their spatial coverage.  Pass inclusive=True to
            retain these.
        """
        self.inclusive = inclusive

        def tomoc(s):
            return _AsIs("MOC({}, {})".format(order, s))

        if isinstance(geom_spec, str):
            geom = _AsIs("MOC({})".format(
                make_sql_literal(geom_spec)))

        elif isinstance(geom_spec, SkyCoord):
            geom = tomoc(format_function_call("POINT",
                                              (geom_spec.ra.value, geom_spec.dec.value)))

        elif len(geom_spec) == 2:
            if isinstance(geom_spec[0], SkyCoord):
                # If radius given is astropy quantity, then convert to degrees
                if isinstance(geom_spec[1], u.Quantity):
                    if geom_spec[1].unit.physical_type != 'angle':
                        raise ValueError("Radius quantity is not of type angle.")
                    radius = geom_spec[1].to(u.deg).value
                else:
                    radius = geom_spec[1]
                geom = tomoc(format_function_call("CIRCLE",
                                                  [geom_spec[0].ra.value, geom_spec[0].dec.value,
                                                   radius]))
            else:
                geom = tomoc(format_function_call("POINT", geom_spec))

        elif len(geom_spec) == 3:
            geom = tomoc(format_function_call("CIRCLE", geom_spec))

        elif len(geom_spec) % 2 == 0:
            geom = tomoc(format_function_call("POLYGON", geom_spec))

        else:
            raise ValueError("This constraint needs DALI-style geometries.")

        if intersect == "covers":
            self._condition = f"1 = CONTAINS({geom}, coverage)"
        elif intersect == "enclosed":
            self._condition = f"1 = CONTAINS(coverage, {geom})"
        elif intersect == "overlaps":
            self._condition = f"1 = INTERSECTS(coverage, {geom})"
        else:
            raise ValueError("'intersect' should be one of 'covers', 'enclosed', or 'overlaps' "
                             f"but its current value is '{intersect}'.")

    def get_search_condition(self, service):
        # we *could* make this a bit less demanding on the server
        # if we MOC-ified the geometries locally -- but then we'd
        # have to depend on pymoc, and that's too high a price for
        # something as esoteric as a server that understands
        # MOC-based geometries but does not have a MOC function.
        if not service.get_tap_capability().get_adql().get_feature(
                "ivo://org.gavo.dc/std/exts#extra-adql-keywords", "MOC"):
            raise RegTAPFeatureMissing(
                "Current RegTAP service does not support MOC.")

        # We should compare case-insensitively here, but then we don't
        # with delimited identifiers -- in the end, that would have to
        # be handled in dal.vosi.VOSITables.
        if "rr.stc_spatial" not in service.tables:
            raise RegTAPFeatureMissing(
                "stc_spatial missing on current RegTAP service")

        cond = super().get_search_condition(service)
        if self.inclusive:
            return cond[:-1]+" OR coverage IS NULL)"
        else:
            return cond


class Spectral(SubqueriedConstraint):
    """
    A RegTAP constraint on the spectral coverage of resources.

    This is a RegTAP 1.2 extension not yet available on all Registries
    (in 2022).  Worse, not too many resources bother declaring this
    at this point.  For robustness, it might be preferable to use
    the `Waveband` constraint for the time being..

    This constraint accepts quantities, i.e., values with units, and will
    convert them to RegTAP's representation (which is Joule of particle energy)
    if it can. This ought to work for wavelengths, frequencies, and energies.
    Plain numbers are interpreted as particle energies in Joule.

    RegTAP uses the observer frame at the solar system barycenter, but
    it is probably wise to use constraints suitably relaxed such that
    frame and reference position (within reason) do not matter.

    To find resources covering the messenger particle energy 5 eV::

        >>> from astropy import units as u
        >>> from pyvo import registry
        >>> resources =  registry.Spectral(5*u.eV)

    To find resources overlapping the band between 5000 and 6000 Ångström::

        >>> resources =  registry.Spectral((5000*u.Angstrom, 6000*u.Angstrom))

    To find resources having data in the FM band::

        >>> resources =  registry.Spectral((88*u.MHz, 102*u.MHz))
    """
    _keyword = "spectral"
    _subquery_table = "rr.stc_spectral"

    takes_sequence = True

    def __init__(self, spec, *, inclusive=False):
        """

        Parameters
        ----------
        spec : astropy.Quantity or a 2-tuple of astropy.Quantity-s
            A spectral point or interval to cover.  This must be a wavelength,
            a frequency, or an energy, or a pair of such quantities,
            in which case the argument is interpreted as an interval.
            All resources *overlapping* the interval are returned.
            Plain floats are interpreted as messenger energy in Joule.

        inclusive : bool, optional
            Normally, this constraint will remove all resources that do
            not declare their spectral coverage.  Pass inclusive=True to
            retain these.
        """
        self.inclusive = inclusive

        if isinstance(spec, tuple):
            self._fillers = {
                "spec_lo": self._to_joule(spec[0]),
                "spec_hi": self._to_joule(spec[1])}
            self._condition = ("1 = ivo_interval_overlaps("
                               "spectral_start, spectral_end, {spec_lo}, {spec_hi})")

        else:
            self._fillers = {
                "spec": self._to_joule(spec)}
            self._condition = "{spec} BETWEEN spectral_start AND spectral_end"

    def _to_joule(self, quant):
        """returns a spectral quantity as a float in joule.

        A plain float is returned as-is.
        """
        if isinstance(quant, (float, int)):
            return quant

        try:
            # is it an energy?
            return quant.to(u.Joule).value
        except u.UnitConversionError:
            pass  # try next

        try:
            # is it a wavelength?
            return (constants.h * constants.c / quant.to(u.m)).value
        except u.UnitConversionError:
            pass  # try next

        try:
            # is it a frequency?
            return (constants.h * quant.to(u.Hz)).value
        except u.UnitConversionError:
            pass  # fall through to give up

        raise ValueError(f"Cannot make a spectral quantity out of {quant}")

    def get_search_condition(self, service):
        if "rr.stc_spectral" not in service.tables:
            raise RegTAPFeatureMissing(
                "stc_spectral missing on current RegTAP service")
        cond = super().get_search_condition(service)
        if self.inclusive:
            return (f"({cond}) OR NOT EXISTS("
                "SELECT 1 FROM rr.stc_spectral AS inner_s WHERE"
                    " inner_s.ivoid=rr.resource.ivoid)")
        else:
            return cond


class Temporal(SubqueriedConstraint):
    """
    A RegTAP constraint on the temporal coverage of resources.

    This is a RegTAP 1.2 extension not yet available on all Registries
    (in 2022).  Worse, not too many resources bother declaring this
    at this point.  Until this changes, you will probably have a lot of false
    negatives (i.e., resources that should match but do not because they
    are not declaring their time coverage) if you use this constraint.

    This constraint accepts astropy Time instances or pairs of Times
    when specifying intervals.  Plain numbers will be interpreted as
    MJD.  RegTAP uses TDB times at the solar system barycenter, and it is
    probably wise to relax constraints such that such details do not matter.
    This constraint does not attempt any conversions of time scales or
    reference positions.

    To find resources claiming to have data for Jan 10, 2022::

        >>> from pyvo import registry
        >>> from astropy.time import Time
        >>> resources = registry.Temporal(Time('2022-01-10'))

    To find resources claiming to have data for some time between
    MJD 54130 and 54200::

        >>> resources = registry.Temporal((54130, 54200))
    """
    _keyword = "temporal"
    _subquery_table = "rr.stc_temporal"

    takes_sequence = True

    def __init__(self, times, *, inclusive=False):
        """

        Parameters
        ----------
        spec : astropy.Time or a 2-tuple of astropy.Time-s
            A point in time or time interval to cover.  Plain numbers
            are interpreted as MJD.  All resources *overlapping* the
            interval are returned.

        inclusive : bool, optional
            Normally, this constraint will remove all resources that do
            not declare their temproal coverage.  Pass inclusive=True to
            retain these.
        """
        self.inclusive = inclusive

        if isinstance(times, tuple):
            self._fillers = {
                "time_lo": self._to_mjd(times[0]),
                "time_hi": self._to_mjd(times[1])}
            self._condition = ("1 = ivo_interval_overlaps("
                               "time_start, time_end, {time_lo}, {time_hi})")

        else:
            self._fillers = {
                "time": self._to_mjd(times)}
            self._condition = "{time} BETWEEN time_start AND time_end"

    def _to_mjd(self, quant):
        """returns a time specification in MJD.

        Times not corresponding to a single point in time are rejected.

        A plain float is returned as-is.
        """
        if isinstance(quant, (float, int)):
            return quant

        val = quant.to_value('mjd')
        if not isinstance(val, numpy.number):
            raise ValueError("RegTAP time constraints must be made from"
                             " single time instants.")
        return val

    def get_search_condition(self, service):
        if "rr.stc_temporal" not in service.tables:
            raise RegTAPFeatureMissing(
                "stc_temporal missing on current RegTAP service")

        cond = super().get_search_condition(service)
        if self.inclusive:
            return (f"({cond}) OR NOT EXISTS("
                "SELECT 1 FROM rr.stc_temporal AS inner_t WHERE"
                    " inner_t.ivoid=rr.resource.ivoid)")
        else:
            return cond


# NOTE: If you add new Contraint-s, don't forget to add them in
# registry.__init__, in docs/registry/index.rst and in the docstring
# of regtap.query.


def build_regtap_query(constraints, service):
    """returns a RegTAP query ready for submission from a list of
    Constraint instances.

    Parameters
    ----------
    constraints : sequence of `~pyvo.registry.Constraint`-s
        A sequence of constraints for a RegTAP query.  All of them
        will become part of a conjunction (i.e., all of them have
        to be satisfied for a record to match).

    service : `~pyvo.dal.TAPService`
        The RegTAP service the query is supposed to be run on
        (that is relevant because we adapt to the features available
        on given services).

    Returns
    -------
    str
        An ADQL literal ready for submission to a RegTAP service.
    """
    if not constraints:
        raise dalq.DALQueryError(
            "No search parameters passed to registry search")

    serialized, extra_tables = [], set()
    for constraint in constraints:
        if isinstance(constraint, str):
            constraint = Freetext(constraint)

        serialized.append(
            "(" + constraint.get_search_condition(service) + ")")
        extra_tables |= set(constraint._extra_tables)

    joined_tables = ["rr.resource", "rr.capability", "rr.interface"
                     ] + list(sorted(extra_tables))

    # see comment in regtap.RegistryResource for the following
    # oddity
    select_clause, plain_columns = [], []
    for col_desc in regtap.RegistryResource.expected_columns:
        if isinstance(col_desc, str):
            select_clause.append(col_desc)
            plain_columns.append(col_desc)
        else:
            select_clause.append("{} AS {}".format(*col_desc))

    fragments = ["SELECT",
                 ", ".join(select_clause),
                 "FROM",
                 "\nNATURAL LEFT OUTER JOIN ".join(joined_tables),
                 "WHERE",
                 "\n  AND ".join(serialized),
                 "GROUP BY",
                 ", ".join(plain_columns)]

    return "\n".join(fragments)


def keywords_to_constraints(keywords):
    """returns constraints expressed as keywords as Constraint instances.

    Parameters
    ----------
    keywords : dict
        regsearch arguments as a kwargs-style dictionary.

    Returns
    -------
    sequence of `Constraint`-s

    Raises
    ------
    DALQueryError
        if an unknown keyword is encountered.
    """
    constraints = []
    for keyword, value in keywords.items():
        if keyword not in _KEYWORD_TO_CONSTRAINT:
            raise TypeError(f"{keyword} is not a valid registry"
                            " constraint keyword.  Use one of {}.".format(
                                ", ".join(sorted(_KEYWORD_TO_CONSTRAINT))))

        constraint_class = _KEYWORD_TO_CONSTRAINT[keyword]
        if (isinstance(value, (tuple, list))
                and not constraint_class.takes_sequence):
            constraints.append(constraint_class(*value))
        else:
            constraints.append(constraint_class(value))
    return constraints


def _make_constraint_map():
    """returns a map of _keyword to constraint classes.

    This is used in module initialisation.
    """
    keyword_to_constraint = {}
    for att_name, obj in globals().items():
        if (isinstance(obj, type)
                and issubclass(obj, Constraint)
                and obj._keyword):
            keyword_to_constraint[obj._keyword] = obj
    return keyword_to_constraint


_KEYWORD_TO_CONSTRAINT = _make_constraint_map()
