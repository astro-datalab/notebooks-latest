#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.registry.rtcons, i.e. RegTAP constraints and query building.
"""

import datetime

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils.exceptions import AstropyDeprecationWarning

import numpy
import pytest

from pyvo import registry
from pyvo.registry import rtcons
from pyvo.dal import query as dalq

from .commonfixtures import messenger_vocabulary, FAKE_GAVO, FAKE_PLAIN  # noqa: F401


def _build_regtap_query_with_fake(
        *args,
        service=FAKE_GAVO,
        **kwargs):
    return rtcons.build_regtap_query(
        *args, service=service, **kwargs)


def _make_subquery(table, condition):
    """returns how condition would show up in something produced by
    SubqueriedConstraint.
    """
    return f"ivoid IN (SELECT DISTINCT ivoid FROM {table} WHERE {condition})"


class TestAbstractConstraint:
    def test_no_search_condition(self):
        with pytest.raises(NotImplementedError):
            rtcons.Constraint().get_search_condition(FAKE_GAVO)


class TestSQLLiterals:
    @pytest.fixture(scope="class", autouse=True)
    def literals(self):
        class _WithFillers(rtcons.Constraint):
            _fillers = {
                "aString": "some harmless stuff",
                "nastyString": "that's not nasty",
                "bytes": b"keep this ascii for now",
                "anInt": 210,
                "aFloat": 5e7,
                "numpyStuff": numpy.float64(23.7),
                "timestamp": datetime.datetime(2021, 6, 30, 9, 1), }

        return _WithFillers()._get_sql_literals()

    def test_strings(self, literals):
        assert literals["aString"] == "'some harmless stuff'"
        assert literals["nastyString"] == "'that''s not nasty'"

    def test_bytes(self, literals):
        assert literals["bytes"] == "'keep this ascii for now'"

    def test_int(self, literals):
        assert literals["anInt"] == "210"

    def test_float(self, literals):
        assert literals["aFloat"] == "50000000.0"

    def test_numpy(self, literals):
        assert float(literals["numpyStuff"]) - 23.7 < 1e-10

    def test_timestamp(self, literals):
        assert literals["timestamp"] == "'2021-06-30T09:01:00'"

    def test_odd_type_rejected(self):
        with pytest.raises(ValueError) as excinfo:
            rtcons.make_sql_literal({})
        assert str(excinfo.value) == "Cannot format {} as a SQL literal"


class TestFreetextConstraint:
    def test_basic(self):
        assert rtcons.Freetext("star").get_search_condition(FAKE_GAVO) == (
            "ivoid IN (SELECT DISTINCT ivoid FROM rr.resource"
            " WHERE 1=ivo_hasword(res_description, 'star') "
            "UNION ALL SELECT DISTINCT ivoid FROM rr.resource"
            " WHERE 1=ivo_hasword(res_title, 'star') "
            "UNION ALL SELECT DISTINCT ivoid FROM rr.res_subject"
            " WHERE rr.res_subject.res_subject ILIKE '%star%')")

    def test_interesting_literal(self):
        assert rtcons.Freetext("α Cen's planets").get_search_condition(FAKE_GAVO) == (
            "ivoid IN (SELECT DISTINCT ivoid FROM rr.resource"
            " WHERE 1=ivo_hasword(res_description, 'α Cen''s planets') "
            "UNION ALL SELECT DISTINCT ivoid FROM rr.resource"
            " WHERE 1=ivo_hasword(res_title, 'α Cen''s planets') "
            "UNION ALL SELECT DISTINCT ivoid FROM rr.res_subject"
            " WHERE rr.res_subject.res_subject"
            " ILIKE '%α Cen''s planets%')")

    def test_multipleLiterals(self):
        assert rtcons.Freetext("term1", "term2").get_search_condition(FAKE_GAVO) == (
            "ivoid IN (SELECT DISTINCT ivoid FROM rr.resource"
            " WHERE 1=ivo_hasword(res_description, 'term1') "
            "UNION ALL SELECT DISTINCT ivoid FROM rr.resource"
            " WHERE 1=ivo_hasword(res_title, 'term1') "
            "UNION ALL SELECT DISTINCT ivoid FROM rr.res_subject"
            " WHERE rr.res_subject.res_subject ILIKE '%term1%')"
            " AND "
            "ivoid IN (SELECT DISTINCT ivoid FROM rr.resource"
            " WHERE 1=ivo_hasword(res_description, 'term2') "
            "UNION ALL SELECT DISTINCT ivoid FROM rr.resource"
            " WHERE 1=ivo_hasword(res_title, 'term2') "
            "UNION ALL SELECT DISTINCT ivoid FROM rr.res_subject"
            " WHERE rr.res_subject.res_subject ILIKE '%term2%')")

    def test_adaption_to_service(self):
        assert rtcons.Freetext("term1", "term2").get_search_condition(FAKE_PLAIN) == (
            "( 1=ivo_hasword(res_description, 'term1') OR  1=ivo_hasword(res_title, 'term1')"
            " OR  rr.res_subject.res_subject ILIKE '%term1%')"
            " AND ( 1=ivo_hasword(res_description, 'term2') OR  1=ivo_hasword(res_title, 'term2')"
            " OR  rr.res_subject.res_subject ILIKE '%term2%')")


class TestAuthorConstraint:
    def test_basic(self):
        assert (rtcons.Author("%Hubble%").get_search_condition(FAKE_GAVO)
                == _make_subquery("rr.res_role", "role_name LIKE '%Hubble%' AND base_role='creator'"))


class TestServicetypeConstraint:
    def test_standardmap(self):
        assert (rtcons.Servicetype("scs").get_search_condition(FAKE_GAVO)
                == "standard_id IN ('ivo://ivoa.net/std/conesearch')")

    def test_fulluri(self):
        assert (rtcons.Servicetype("http://extstandards/invention"
                                   ).get_search_condition(FAKE_GAVO)
                == "standard_id IN ('http://extstandards/invention')")

    def test_multi(self):
        assert (rtcons.Servicetype("http://extstandards/invention", "sia"
                                   ).get_search_condition(FAKE_GAVO)
                == "standard_id IN ('http://extstandards/invention',"
                " 'ivo://ivoa.net/std/sia')")

    def test_includeaux(self):
        assert (rtcons.Servicetype("http://extstandards/invention", "sia1"
                                   ).include_auxiliary_services().get_search_condition(FAKE_GAVO)
                == "standard_id IN ('http://extstandards/invention',"
                " 'http://extstandards/invention#aux',"
                " 'ivo://ivoa.net/std/sia',"
                " 'ivo://ivoa.net/std/sia#aux')")

    def test_junk_rejected(self):
        with pytest.raises(dalq.DALQueryError) as excinfo:
            rtcons.Servicetype("junk")
        assert str(excinfo.value) == ("Service type junk is neither"
                                      " a full standard URI nor one of the bespoke identifiers"
                                      " image, sia, sia1, spectrum, ssap, ssa, scs, conesearch, line, slap,"
                                      " table, tap, sia2")

    def test_legacy_term(self):
        assert (rtcons.Servicetype("conesearch").get_search_condition(FAKE_GAVO)
                == "standard_id IN ('ivo://ivoa.net/std/conesearch')")

    def test_sia2(self):
        assert (
            rtcons.Servicetype("conesearch", "sia2").get_search_condition(FAKE_GAVO)
            == ("standard_id IN ('ivo://ivoa.net/std/conesearch')"
                " OR standard_id like 'ivo://ivoa.net/std/sia#query-2.%'"))

    def test_sia2_aux(self):
        constraint = rtcons.Servicetype("conesearch", "sia2").include_auxiliary_services()
        assert (constraint.get_search_condition(FAKE_GAVO)
                == ("standard_id IN ('ivo://ivoa.net/std/conesearch', 'ivo://ivoa.net/std/conesearch#aux')"
                    " OR standard_id like 'ivo://ivoa.net/std/sia#query-2.%'"
                    " OR standard_id like 'ivo://ivoa.net/std/sia#query-aux-2.%'"))

    def test_image_deprecated(self):
        with pytest.warns(AstropyDeprecationWarning):
            assert (rtcons.Servicetype("image").get_search_condition(FAKE_GAVO)
                == "standard_id IN ('ivo://ivoa.net/std/sia')")

    def test_spectrum_deprecated(self):
        with pytest.warns(AstropyDeprecationWarning):
            assert (rtcons.Servicetype("spectrum").get_search_condition(FAKE_GAVO)
                == "standard_id IN ('ivo://ivoa.net/std/ssa')")


@pytest.mark.usefixtures('messenger_vocabulary')
class TestWavebandConstraint:
    def test_basic(self):
        assert (rtcons.Waveband("Infrared", "EUV").get_search_condition(FAKE_GAVO)
                == "1 = ivo_hashlist_has(rr.resource.waveband, 'infrared')"
                " OR 1 = ivo_hashlist_has(rr.resource.waveband, 'euv')")

    def test_junk_rejected(self):
        with pytest.raises(dalq.DALQueryError) as excinfo:
            rtcons.Waveband("junk")
        assert str(excinfo.value) == (
            "Waveband junk is not in the IVOA messenger vocabulary http://www.ivoa.net/rdf/messenger.")

    def test_normalisation(self):
        assert (rtcons.Waveband("oPtIcAl").get_search_condition(FAKE_GAVO)
                == "1 = ivo_hashlist_has(rr.resource.waveband, 'optical')")


class TestDatamodelConstraint:
    def test_junk_rejected(self):
        with pytest.raises(dalq.DALQueryError) as excinfo:
            rtcons.Datamodel("junk")
        assert str(excinfo.value) == (
            "Unknown data model id junk.  Known are: epntap, obscore, obscore_new, regtap.")

    def test_obscore_new(self):
        cons = rtcons.Datamodel("obscore_new")
        assert (cons.get_search_condition(FAKE_GAVO)
            == "ivoid IN (SELECT DISTINCT ivoid FROM rr.res_table"
               " NATURAL JOIN rr.resource"
               " WHERE table_utype LIKE 'ivo://ivoa.net/std/obscore#table-1.%'"
               " AND res_type = 'vs:catalogresource')")
        assert (cons._extra_tables == [])

    def test_obscore(self):
        cons = rtcons.Datamodel("ObsCore")
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.res_detail",
                    "detail_xpath = '/capability/dataModel/@ivo-id'"
                    " AND 1 = ivo_nocasematch(detail_value,"
                    " 'ivo://ivoa.net/std/obscore%')"))
        assert cons._extra_tables == []

    def test_epntap(self):
        cons = rtcons.Datamodel("epntap")
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.res_table",
                    "table_utype LIKE 'ivo://vopdc.obspm/std/epncore#schema-2.%'"
                    " OR table_utype LIKE 'ivo://ivoa.net/std/epntap#table-2.%'"))
        assert cons._extra_tables == []

    def test_regtap(self):
        cons = rtcons.Datamodel("regtap")
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.res_detail",
                    "detail_xpath = '/capability/dataModel/@ivo-id'"
                    " AND 1 = ivo_nocasematch(detail_value,"
                    " 'ivo://ivoa.net/std/RegTAP#1.%')"))
        assert cons._extra_tables == []


class TestIvoidConstraint:
    def test_basic(self):
        cons = rtcons.Ivoid("ivo://example/some_path")
        assert (cons.get_search_condition(FAKE_GAVO)
                == "ivoid='ivo://example/some_path'")

    def test_multiple(self):
        cons = rtcons.Ivoid(
            "ivo://org.gavo.dc/tap",
            "ivo://org.gavo.dc/__system__/siap2/sitewide")
        assert (cons.get_search_condition(FAKE_GAVO)
                == ("ivoid='ivo://org.gavo.dc/tap'"
                    " OR ivoid='ivo://org.gavo.dc/__system__/siap2/sitewide'"))


class TestUCDConstraint:
    def test_basic(self):
        cons = rtcons.UCD("phot.mag;em.opt.%", "phot.mag;em.ir.%")
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.table_column",
                    "ucd LIKE 'phot.mag;em.opt.%' OR ucd LIKE 'phot.mag;em.ir.%'"))


class TestSpatialConstraint:
    def test_point(self):
        cons = registry.Spatial([23, -40])
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spatial", "1 = CONTAINS(MOC(6, POINT(23, -40)), coverage)")
        assert cons._extra_tables == []

    def test_circle_and_order(self):
        cons = registry.Spatial([23, -40, 0.25], order=7)
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spatial", "1 = CONTAINS(MOC(7, CIRCLE(23, -40, 0.25)), coverage)")

    def test_polygon(self):
        cons = registry.Spatial([23, -40, 26, -39, 25, -43])
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spatial", "1 = CONTAINS(MOC(6, POLYGON(23, -40, 26, -39, 25, -43)), coverage)")

    def test_moc(self):
        cons = registry.Spatial("0/1-3 3/")
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spatial", "1 = CONTAINS(MOC('0/1-3 3/'), coverage)")

    def test_moc_and_inclusive(self):
        cons = registry.Spatial("0/1-3 3/", inclusive=True)
        assert cons.get_search_condition(FAKE_GAVO) == (
            "ivoid IN (SELECT DISTINCT ivoid FROM rr.stc_spatial WHERE 1 = "
            "CONTAINS(MOC('0/1-3 3/'), coverage) OR coverage IS NULL)")

    def test_SkyCoord(self):
        cons = registry.Spatial(SkyCoord(3 * u.deg, -30 * u.deg))
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spatial", "1 = CONTAINS(MOC(6, POINT(3.0, -30.0)), coverage)")

    def test_SkyCoord_Circle(self):
        cons = registry.Spatial((SkyCoord(3 * u.deg, -30 * u.deg), 3))
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spatial", "1 = CONTAINS(MOC(6, CIRCLE(3.0, -30.0, 3)), coverage)")

    def test_SkyCoord_Circle_RadiusQuantity(self):
        for radius in [3*u.deg, 180*u.Unit('arcmin'), 10800*u.Unit('arcsec')]:
            cons = registry.Spatial((SkyCoord(3 * u.deg, -30 * u.deg), radius))
            assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
                "rr.stc_spatial", "1 = CONTAINS(MOC(6, CIRCLE(3.0, -30.0, 3.0)), coverage)")

        with pytest.raises(ValueError, match="is not of type angle."):
            cons = registry.Spatial((SkyCoord(3 * u.deg, -30 * u.deg), (1 * u.m)))

    def test_enclosed(self):
        cons = registry.Spatial("0/1-3", intersect="enclosed")
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spatial", "1 = CONTAINS(coverage, MOC('0/1-3'))")

    def test_overlaps(self):
        cons = registry.Spatial("0/1-3", intersect="overlaps")
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spatial", "1 = INTERSECTS(coverage, MOC('0/1-3'))")

    def test_not_an_intersect_mode(self):
        with pytest.raises(ValueError, match="'intersect' should be one of 'covers', 'enclosed',"
                           " or 'overlaps' but its current value is 'wrong'."):
            registry.Spatial("0/1-3", intersect="wrong")

    def test_no_MOC(self):
        cons = registry.Spatial((SkyCoord(3 * u.deg, -30 * u.deg), 3))
        with pytest.raises(rtcons.RegTAPFeatureMissing) as excinfo:
            cons.get_search_condition(FAKE_PLAIN)
        assert (str(excinfo.value)
            == "Current RegTAP service does not support MOC.")

    def test_no_spatial_table(self):
        cons = registry.Spatial((SkyCoord(3 * u.deg, -30 * u.deg), 3))
        previous = FAKE_GAVO.tables.pop("rr.stc_spatial")
        try:
            with pytest.raises(rtcons.RegTAPFeatureMissing) as excinfo:
                cons.get_search_condition(FAKE_GAVO)
            assert (str(excinfo.value)
                == "stc_spatial missing on current RegTAP service")
        finally:
            FAKE_GAVO.tables["rr.spatial"] = previous


class TestSpectralConstraint:
    # These tests might need some float literal fuzziness.  I'm just
    # too lazy at this point to see if pytest has something on board
    # that would be useful there.
    def test_energy_float(self):
        cons = registry.Spectral(1e-19)
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spectral", "1e-19 BETWEEN spectral_start AND spectral_end")

    def test_energy_eV(self):
        cons = registry.Spectral(5 * u.eV)
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.stc_spectral",
                    "8.01088317e-19 BETWEEN spectral_start AND spectral_end"))

    def test_energy_interval(self):
        cons = registry.Spectral((1e-10 * u.erg, 2e-10 * u.erg))
        assert cons.get_search_condition(FAKE_GAVO) == _make_subquery(
            "rr.stc_spectral",
            "1 = ivo_interval_overlaps(spectral_start, spectral_end, 1e-17, 2e-17)")

    def test_wavelength(self):
        cons = registry.Spectral(5000 * u.Angstrom)
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.stc_spectral",
                    "3.9728917142978567e-19 BETWEEN spectral_start AND spectral_end"))

    def test_wavelength_interval(self):
        cons = registry.Spectral((20 * u.cm, 22 * u.cm))
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.stc_spectral",
                    "1 = ivo_interval_overlaps(spectral_start, spectral_end,"
                    " 9.932229285744642e-25, 9.029299350676949e-25)"))

    def test_frequency(self):
        cons = registry.Spectral(2 * u.GHz)
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.stc_spectral",
                    "1.32521403e-24 BETWEEN spectral_start AND spectral_end"))

    def test_frequency_and_inclusive(self):
        cons = registry.Spectral(2 * u.GHz, inclusive=True)
        assert (cons.get_search_condition(FAKE_GAVO)
                == "(ivoid IN (SELECT DISTINCT ivoid FROM rr.stc_spectral"
                " WHERE 1.32521403e-24 BETWEEN spectral_start AND spectral_end))"
                " OR NOT EXISTS(SELECT 1 FROM rr.stc_spectral AS inner_s"
                " WHERE inner_s.ivoid=rr.resource.ivoid)")

    def test_frequency_interval(self):
        cons = registry.Spectral((88 * u.MHz, 102 * u.MHz))
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.stc_spectral",
                    "1 = ivo_interval_overlaps("
                    "spectral_start, spectral_end, 5.830941732e-26, 6.758591553e-26)"))

    def test_no_spectral(self):
        cons = registry.Spectral((88 * u.MHz, 102 * u.MHz))
        with pytest.raises(rtcons.RegTAPFeatureMissing) as excinfo:
            cons.get_search_condition(FAKE_PLAIN)
        assert (str(excinfo.value)
            == "stc_spectral missing on current RegTAP service")


class TestTemporalConstraint:
    def test_plain_float(self):
        cons = registry.Temporal((54130, 54200))
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.stc_temporal",
                    "1 = ivo_interval_overlaps(time_start, time_end, 54130, 54200)"))

    def test_plain_float_and_inclusive(self):
        cons = registry.Temporal((54130, 54200), inclusive=True)
        assert (cons.get_search_condition(FAKE_GAVO)
                == "(ivoid IN (SELECT DISTINCT ivoid FROM rr.stc_temporal"
                " WHERE 1 = ivo_interval_overlaps(time_start, time_end, 54130, 54200)))"
                " OR NOT EXISTS(SELECT 1 FROM rr.stc_temporal AS inner_t"
                " WHERE inner_t.ivoid=rr.resource.ivoid)")

    def test_single_time(self):
        cons = registry.Temporal(Time('2022-01-10'))
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery("rr.stc_temporal", "59589.0 BETWEEN time_start AND time_end"))

    def test_time_interval(self):
        cons = registry.Temporal((Time(2459000, format='jd'),
                                  Time(59002, format='mjd')))
        assert (cons.get_search_condition(FAKE_GAVO)
                == _make_subquery(
                    "rr.stc_temporal",
                    "1 = ivo_interval_overlaps(time_start, time_end, 58999.5, 59002.0)"))

    def test_multi_times_rejected(self):
        with pytest.raises(ValueError) as excinfo:
            registry.Temporal(Time(['1999-01-01', '2010-01-01']))
        assert (str(excinfo.value) == "RegTAP time constraints must"
                " be made from single time instants.")

    def test_no_temporal(self):
        cons = registry.Temporal((Time(2459000, format='jd'),
                                  Time(59002, format='mjd')))
        with pytest.raises(rtcons.RegTAPFeatureMissing) as excinfo:
            cons.get_search_condition(FAKE_PLAIN)
        assert (str(excinfo.value)
            == "stc_temporal missing on current RegTAP service")


class TestWhereClauseBuilding:
    @staticmethod
    def where_clause_for(*args, **kwargs):
        cons = list(args) + rtcons.keywords_to_constraints(kwargs)
        return _build_regtap_query_with_fake(cons).split("\nWHERE\n", 1)[1].split("\nGROUP BY\n")[0]

    @pytest.mark.usefixtures('messenger_vocabulary')
    def test_from_constraints(self):
        assert self.where_clause_for(
            rtcons.Waveband("EUV"),
            rtcons.Author("%Hubble%"))\
            == ("(1 = ivo_hashlist_has(rr.resource.waveband, 'euv'))\n"
              "  AND ({})".format(
                _make_subquery(
                    "rr.res_role",
                    "role_name LIKE '%Hubble%' AND base_role='creator'")))

    @pytest.mark.usefixtures('messenger_vocabulary')
    def test_from_keywords(self):
        assert self.where_clause_for(
            waveband="EUV",
            author="%Hubble%")\
            == ("(1 = ivo_hashlist_has(rr.resource.waveband, 'euv'))\n  AND ({})".format(
                _make_subquery("rr.res_role", "role_name LIKE '%Hubble%' AND base_role='creator'")))

    @pytest.mark.usefixtures('messenger_vocabulary')
    def test_mixed(self):
        assert self.where_clause_for(
            rtcons.Waveband("EUV"),
            author="%Hubble%")\
            == ("(1 = ivo_hashlist_has(rr.resource.waveband, 'euv'))\n"
              "  AND ({})".format(
                _make_subquery("rr.res_role", "role_name LIKE '%Hubble%' AND base_role='creator'")))

    def test_bad_keyword(self):
        with pytest.raises(TypeError) as excinfo:
            _build_regtap_query_with_fake(
                *rtcons.keywords_to_constraints({"foo": "bar"}))
        # the following assertion will fail when new constraints are
        # defined (or old ones vanish).  I'd say that's a convenient
        # way to track changes; so, let's update the assertion as we
        # go.
        assert str(excinfo.value) == ("foo is not a valid registry"
                                      " constraint keyword.  Use one of"
                                      " author, datamodel, ivoid, keywords, servicetype,"
                                      " spatial, spectral, temporal, ucd, waveband.")

    def test_with_legacy_keyword(self):
        assert self.where_clause_for(
            "plain", "string"
        ) == (
            '(ivoid IN (SELECT DISTINCT ivoid FROM rr.resource WHERE '
            "1=ivo_hasword(res_description, 'plain') UNION ALL SELECT DISTINCT ivoid FROM rr.resource "
            "WHERE 1=ivo_hasword(res_title, 'plain') UNION ALL SELECT DISTINCT ivoid FROM "
            "rr.res_subject WHERE rr.res_subject.res_subject ILIKE '%plain%'))\n"
            '  AND (ivoid IN (SELECT DISTINCT ivoid FROM rr.resource WHERE '
            "1=ivo_hasword(res_description, 'string') UNION ALL SELECT DISTINCT ivoid FROM rr.resource "
            "WHERE 1=ivo_hasword(res_title, 'string') UNION ALL SELECT DISTINCT ivoid FROM "
            "rr.res_subject WHERE rr.res_subject.res_subject ILIKE '%string%'))")


class TestSelectClause:
    def test_expected_columns(self):
        # This will break as regtap.RegistryResource.expected_columns
        # is changed.  Just update the assertion then.
        assert _build_regtap_query_with_fake(
            rtcons.keywords_to_constraints({"author": "%Hubble%"})
        ).split("\nFROM\nrr.resource\n")[0] == (
            "SELECT\n"
            "ivoid, "
            "res_type, "
            "short_name, "
            "res_title, "
            "content_level, "
            "res_description, "
            "reference_url, "
            "creator_seq, "
            "created, "
            "updated, "
            "rights, "
            "content_type, "
            "source_format, "
            "source_value, "
            "region_of_regard, "
            "waveband, "
            "\n  ivo_string_agg(COALESCE(access_url, ''), ':::py VO sep:::') AS access_urls, "
            "\n  ivo_string_agg(COALESCE(standard_id, ''), ':::py VO sep:::') AS standard_ids, "
            "\n  ivo_string_agg(COALESCE(intf_type, ''), ':::py VO sep:::') AS intf_types, "
            "\n  ivo_string_agg(COALESCE(intf_role, ''), ':::py VO sep:::') AS intf_roles, "
            "\n  ivo_string_agg(COALESCE(cap_description, ''), ':::py VO sep:::') AS cap_descriptions")

    def test_group_by_columns(self):
        # Again, this will break as regtap.RegistryResource.expected_columns
        # is changed.  Just update the assertion then.
        assert (_build_regtap_query_with_fake([rtcons.Author("%Hubble%")]).split("\nGROUP BY\n")[-1]
                == ("ivoid, "
                    "res_type, "
                    "short_name, "
                    "res_title, "
                    "content_level, "
                    "res_description, "
                    "reference_url, "
                    "creator_seq, "
                    "created, "
                    "updated, "
                    "rights, "
                    "content_type, "
                    "source_format, "
                    "source_value, "
                    "region_of_regard, "
                    "waveband"))

    def test_joined_tables(self):
        expected_tables = [
            # from author constraint
            "rr.res_role",
            # default tables
            "rr.resource",
            "rr.capability",
            "rr.interface",
        ]
        assert all(table in _build_regtap_query_with_fake([rtcons.Author("%Hubble%")])
                   for table in expected_tables)


@pytest.mark.remote_data
def test_all_constraints():
    text = rtcons.Freetext("star")
    author = rtcons.Author(r"%ESA%")
    servicetype = rtcons.Servicetype("tap")
    waveband = rtcons.Waveband("optical")
    datamodel = rtcons.Datamodel("obscore")
    ivoid = rtcons.Ivoid(r"ivoid")
    ucd = rtcons.UCD(r"pos.eq.ra")
    moc = rtcons.Spatial("0/0-11", intersect="overlaps")
    spectral = rtcons.Spectral((5000 * u.Angstrom, 6000 * u.Angstrom))
    time = rtcons.Temporal((50000, 60000))
    result = registry.search(
        text, author, servicetype, waveband, datamodel,
        ivoid, ucd, moc, spectral, time
    )
    assert result.fieldnames == (
        'ivoid', 'res_type', 'short_name',
        'res_title', 'content_level', 'res_description',
        'reference_url', 'creator_seq', 'created', 'updated',
        'rights', 'content_type', 'source_format', 'source_value',
        'region_of_regard', 'waveband', 'access_urls', 'standard_ids',
        'intf_types', 'intf_roles', 'cap_descriptions')
