#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.io.vosi
"""
import io
from operator import eq as equals

import pytest

import pyvo.io.vosi as vosi
import pyvo.io.vosi.vodataservice as vs
import pyvo.io.vosi.tapregext as tr
from pyvo.io.vosi.exceptions import W06

from astropy.utils.data import get_pkg_data_filename


@pytest.fixture(name='parsed_caps')
def _parsed_caps():
    return vosi.parse_capabilities(get_pkg_data_filename(
        "data/capabilities.xml"))


@pytest.mark.usefixtures("parsed_caps")
class TestCapabilities:
    def test_availability(self, parsed_caps):
        assert equals(
            parsed_caps[0].standardid,
            "ivo://ivoa.net/std/VOSI#availability")
        assert isinstance(parsed_caps[0].interfaces[0], vs.ParamHTTP)
        assert parsed_caps[0].interfaces[0].accessurls[0].use == "full"
        assert equals(
            parsed_caps[0].interfaces[0].accessurls[0].content,
            "http://example.org/tap/availability")

    def test_capendpoint(self, parsed_caps):
        assert equals(
            parsed_caps[1].standardid,
            "ivo://ivoa.net/std/VOSI#capabilities")
        assert isinstance(parsed_caps[1].interfaces[0], vs.ParamHTTP)
        assert parsed_caps[1].interfaces[0].accessurls[0].use == "full"
        assert equals(
            parsed_caps[1].interfaces[0].accessurls[0].content,
            "http://example.org/tap/capabilities")

    def test_tablesendpoint(self, parsed_caps):
        assert parsed_caps[2].standardid == "ivo://ivoa.net/std/VOSI#tables"
        assert isinstance(parsed_caps[2].interfaces[0], vs.ParamHTTP)
        assert parsed_caps[2].interfaces[0].accessurls[0].use == "full"
        assert equals(
            parsed_caps[2].interfaces[0].accessurls[0].content,
            "http://example.org/tap/tables")

    def test_type_parsed(self, parsed_caps):
        assert isinstance(parsed_caps[3], tr.TableAccess)

    def test_stdid_parsed(self, parsed_caps):
        assert parsed_caps[3].standardid == "ivo://ivoa.net/std/TAP"

    def test_dm_parsed(self, parsed_caps):
        assert equals(
            parsed_caps[3].datamodels[0].ivo_id,
            "ivo://ivoa.net/std/ObsCore#table-1.1")
        assert parsed_caps[3].datamodels[0].content == "Obscore-1.1"

        assert equals(
            parsed_caps[3].datamodels[1].ivo_id,
            "ivo://ivoa.net/std/RegTAP#1.0")
        assert parsed_caps[3].datamodels[1].content == "Registry 1.0"

    def test_language_parsed(self, parsed_caps):
        assert parsed_caps[3].languages[0].name == "ADQL"
        assert equals(
            parsed_caps[3].languages[0].versions[0].ivo_id,
            "ivo://ivoa.net/std/ADQL#v2.0")
        assert parsed_caps[3].languages[0].versions[0].content == "2.0"
        assert parsed_caps[3].languages[0].description == "ADQL 2.0"

    def test_udfs(self, parsed_caps):
        assert equals(
            parsed_caps[3].languages[0].languagefeaturelists[0].type,
            "ivo://ivoa.net/std/TAPRegExt#features-udf")
        assert equals(
            parsed_caps[3].languages[0].languagefeaturelists[0][0].form,
            "form 1")
        assert equals(
            parsed_caps[3].languages[0].languagefeaturelists[0][
                0].description,
            "description 1")
        assert equals(
            parsed_caps[3].languages[0].languagefeaturelists[0][1].form,
            "form 2")
        assert equals(
            parsed_caps[3].languages[0].languagefeaturelists[0].features[
                1].description,
            "description 2")

    def test_adqlgeos(self, parsed_caps):
        assert equals(
            parsed_caps[3].languages[0].languagefeaturelists[1].type,
            "ivo://ivoa.net/std/TAPRegExt#features-adqlgeo")
        assert equals(
            parsed_caps[3].languages[0].languagefeaturelists[1].features[
                0].form,
            "BOX")
        assert equals(
            parsed_caps[3].languages[0].languagefeaturelists[1].features[
                1].form,
            "POINT")

    def test_outputformats(self, parsed_caps):
        assert equals(
            parsed_caps[3].outputformats[0].ivo_id,
            "ivo://ivoa.net/std/TAPRegExt#output-votable-binary")
        assert parsed_caps[3].outputformats[0].mime == "text/xml"

        assert parsed_caps[3].outputformats[1].ivo_id is None
        assert parsed_caps[3].outputformats[1].mime == "text/html"

    def test_uploadmethods(self, parsed_caps):
        assert equals(
            parsed_caps[3].uploadmethods[0].ivo_id,
            "ivo://ivoa.net/std/TAPRegExt#upload-https")
        assert equals(
            parsed_caps[3].uploadmethods[1].ivo_id,
            "ivo://ivoa.net/std/TAPRegExt#upload-inline")

    def test_temporal_limits(self, parsed_caps):
        assert parsed_caps[3].retentionperiod.default == 172800
        assert parsed_caps[3].executionduration.default == 3600

    def test_spatial_limits(self, parsed_caps):
        assert parsed_caps[3].outputlimit.default.unit == "row"
        assert parsed_caps[3].outputlimit.default.content == 2000

        assert parsed_caps[3].outputlimit.hard.unit == "row"
        assert parsed_caps[3].outputlimit.hard.content == 10000000

        assert parsed_caps[3].uploadlimit.hard.unit == "byte"
        assert parsed_caps[3].uploadlimit.hard.content == 100000000

    def test_multiple_capa_descriptions(self):
        with pytest.warns(W06):
            vosi.parse_capabilities(get_pkg_data_filename(
                'data/capabilities/multiple_capa_descriptions.xml'))

        with pytest.raises(W06):
            vosi.parse_capabilities(get_pkg_data_filename(
                'data/capabilities/multiple_capa_descriptions.xml'),
                pedantic=True)


@pytest.mark.usefixtures("parsed_caps")
class TestInterface:
    def test_interface_parsed(self, parsed_caps):
        assert parsed_caps[3].interfaces[0].accessurls[0].use == "base"
        assert equals(
            parsed_caps[3].interfaces[0].accessurls[0].content,
            "http://example.org/tap")

    def test_mirrors_parsed(self, parsed_caps):
        assert len(parsed_caps[3].interfaces[0].mirrorurls) == 2

    def test_mirrors_have_titles(self, parsed_caps):
        assert ([m.title for m in parsed_caps[3].interfaces[0].mirrorurls]
                == ["https version", "Paris mirror"])

    def test_mirrors_have_urls(self, parsed_caps):
        assert ([m.content for m in parsed_caps[3].interfaces[0].mirrorurls]
                == ['https://example.org/tap', 'https://paris.example.org/tap'])

    def test_testquerystring_parsed(self, parsed_caps):
        assert (parsed_caps[3].interfaces[0].testquerystring.content
                == 'QUERY=SELECT%20*%20FROM%20tap_schema.tables&LANG=ADQL')


@pytest.fixture(name='cap_with_free_prefix')
def _cap_with_free_prefix(recwarn):
    caps = vosi.parse_capabilities(io.BytesIO(b"""
<ns2:capabilities xmlns:ns2="http://www.ivoa.net/xml/VOSICapabilities/v1.0">
  <capability
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xmlns:ns4="http://www.ivoa.net/xml/TAPRegExt/v1.0"
    xsi:type="ns4:TableAccess"
    standardID="ivo://ivoa.net/std/TAP">
    <interface xmlns:ns5="http://www.ivoa.net/xml/VODataService/v1.1"
        xsi:type="ns5:ParamHTTP" version="1.1" role="std">
      <accessURL use="base">https://archive.eso.org/tap_obs</accessURL>
    </interface>
    <dataModel ivo-id="ivo://ivoa.net/std/ObsCore#core-1.1">ObsCore-1.1</dataModel>
    <outputFormat/>
    <language>
      <name>ADQL</name>
      <version ivo-id="ivo://ivoa.net/std/ADQL#v2.0">2.0</version>
      <description>ADQL-2.0</description>
    </language>
  </capability>
</ns2:capabilities>"""))
    return recwarn, caps


# this is a test for when people ignore the canonical prefixes for
# registry documents.
class TestFreePrefixes:
    def test_parses_without_warning(self, cap_with_free_prefix):
        warnings, _ = cap_with_free_prefix
        assert len(warnings) == 0

    def test_parses_as_tapregext(self, cap_with_free_prefix):
        _, cap = cap_with_free_prefix
        assert isinstance(cap[0], tr.TableAccess)
