#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.registry.regtap
"""

import io
import re
from functools import partial
from urllib.parse import parse_qsl

import pytest

from astropy import time

from pyvo.registry import regtap
from pyvo.registry import rtcons
from pyvo.registry.regtap import REGISTRY_BASEURL
from pyvo.registry import search as regsearch
from pyvo.dal import DALOverflowWarning
from pyvo.dal import query as dalq
from pyvo.dal import tap, sia2

from astropy.utils.data import get_pkg_data_contents

from .commonfixtures import messenger_vocabulary, FAKE_GAVO  # noqa: F401


get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')


@pytest.fixture(name='capabilities')
def _capabilities(mocker, scope="session"):
    def callback(request, context):
        return get_pkg_data_contents('data/capabilities.xml')

    with mocker.register_uri(
        'GET', REGISTRY_BASEURL + '/capabilities', content=callback
    ) as matcher:
        yield matcher


# to update this, run
# import requests
# from pyvo.registry import regtap
#
# with open("data/regtap.xml", "wb") as f:
# 	f.write(requests.get(regtap.REGISTRY_BASEURL+"/sync", {
# 		"LANG": "ADQL",
# 		"QUERY": regtap.get_RegTAP_query(keywords="pulsar", ucd=["pos.distance"])}).content)


@pytest.fixture(name='regtap_pulsar_distance_response')
def _regtap_pulsar_distance_response(mocker, capabilities, scope="session"):
    with mocker.register_uri(
        'POST', REGISTRY_BASEURL + '/sync',
            content=get_pkg_data_contents('data/regtap.xml')) as matcher:
        yield matcher


@pytest.fixture()
def keywords_fixture(mocker, capabilities, scope="session"):
    def keywordstest_callback(request, context):
        data = dict(parse_qsl(request.body))
        query = data['QUERY']

        assert "res_subject ILIKE '%vizier%'" in query
        assert "ivo_hasword(res_description, 'vizier')" in query
        assert "1=ivo_hasword(res_title, 'vizier')" in query

        assert ".res_subject ILIKE '%pulsar%'" in query
        assert "1=ivo_hasword(res_description, 'pulsar')" in query
        assert "1=ivo_hasword(res_title, 'pulsar')" in query

        return get_pkg_data_contents('data/regtap.xml')

    with mocker.register_uri(
        'POST', REGISTRY_BASEURL + '/sync',
        content=keywordstest_callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def single_keyword_fixture(mocker, capabilities, scope="session"):
    def keywordstest_callback(request, context):
        data = dict(parse_qsl(request.body))
        query = data['QUERY']

        assert (" UNION ALL SELECT DISTINCT ivoid FROM rr.res_subject WHERE "
                "rr.res_subject.res_subject ILIKE '%single%'") in query
        assert "1=ivo_hasword(res_description, 'single') " in query
        assert "1=ivo_hasword(res_title, 'single')" in query

        return get_pkg_data_contents('data/regtap.xml')

    with mocker.register_uri(
        'POST', REGISTRY_BASEURL + '/sync',
        content=keywordstest_callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def servicetype_fixture(mocker, capabilities, scope="session"):
    def servicetypetest_callback(request, context):
        data = dict(parse_qsl(request.body))
        query = data['QUERY']

        assert "'ivo://ivoa.net/std/conesearch'" not in query
        assert "'ivo://ivoa.net/std/sia'" not in query
        assert "'ivo://ivoa.net/std/ssa'" not in query
        assert "'ivo://ivoa.net/std/slap'" not in query
        assert "'ivo://ivoa.net/std/tap'" in query

        return get_pkg_data_contents('data/regtap.xml')

    with mocker.register_uri(
        'POST', REGISTRY_BASEURL + '/sync',
        content=servicetypetest_callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def waveband_fixture(mocker, capabilities, scope="session"):
    def wavebandtest_callback(request, content):
        data = dict(parse_qsl(request.body))
        query = data['QUERY']

        assert "1 = ivo_hashlist_has(rr.resource.waveband, 'optical'" in query

        return get_pkg_data_contents('data/regtap.xml')

    with mocker.register_uri(
        'POST', REGISTRY_BASEURL + '/sync',
        content=wavebandtest_callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def datamodel_fixture(mocker, capabilities, scope="session"):
    def datamodeltest_callback(request, content):
        data = dict(parse_qsl(request.body))
        query = data['QUERY']

        assert (
            "detail_xpath = '/capability/dataModel/@ivo-id'" in query)

        assert (
            "ivo_nocasematch(detail_value, 'ivo://ivoa.net/std/obscore%'))"
            in query)

        return get_pkg_data_contents('data/regtap.xml')

    with mocker.register_uri(
        'POST', REGISTRY_BASEURL + '/sync',
        content=datamodeltest_callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def aux_fixture(mocker, capabilities, scope="session"):
    def auxtest_callback(request, context):
        data = dict(parse_qsl(request.body))
        query = data['QUERY']

        assert "ivo://ivoa.net/std/tap#aux" in query

        return get_pkg_data_contents('data/regtap.xml')

    with mocker.register_uri(
        'POST', REGISTRY_BASEURL + '/sync',
        content=auxtest_callback,
    ) as matcher:
        yield matcher


@pytest.fixture(name='multi_interface_fixture')
def _multi_interface_fixture(mocker, capabilities, scope="session"):
    # to update this, run
    # import requests
    # from pyvo.registry import regtap
    #
    # with open("data/multi-interface.xml", "wb") as f:
    # 	f.write(requests.get(regtap.REGISTRY_BASEURL+"/sync", {
    # 		"LANG": "ADQL",
    # 		"QUERY": regtap.get_RegTAP_query(
    # 			ivoid="ivo://org.gavo.dc/flashheros/q/ssa")}).content)
    with mocker.register_uri(
        'POST', REGISTRY_BASEURL + '/sync',
        content=get_pkg_data_contents('data/multi-interface.xml')
    ) as matcher:
        yield matcher


@pytest.fixture(name='flash_service')
def _flash_service(multi_interface_fixture, scope="session"):
    return regtap.search(
        ivoid="ivo://org.gavo.dc/flashheros/q/ssa")[0]


class TestInterfaceClass:
    def test_basic(self):
        intf = regtap.Interface("http://example.org")
        assert intf.access_url == "http://example.org"
        assert intf.standard_id is None
        assert intf.type is None
        assert intf.role is None
        assert intf.is_standard is False
        assert not intf.is_vosi
        assert intf.capability_description is None

    def test_repr(self):
        description = "An example description."
        intf = regtap.Interface("http://example.org", standard_id="ivo://gavo/std/a",
                                intf_type="vs:paramhttp", intf_role="std",
                                capability_description=description)
        assert (repr(intf) == "Interface(type='a', description='An example description.',"
                " url='http://example.org')")
        intf = regtap.Interface("http://example.org", capability_description=description)
        assert repr(intf) == ("Interface(type=None, description='An example description.',"
                " url='http://example.org')")

    def test_unknown_standard(self):
        intf = regtap.Interface("http://example.org", standard_id="ivo://gavo/std/a",
                                intf_type="vs:paramhttp", intf_role="std")
        assert intf.is_standard
        with pytest.raises(ValueError) as excinfo:
            intf.to_service()

        assert str(excinfo.value) == (
            "PyVO has no support for interfaces with standard"
            " id ivo://gavo/std/a.")

    def test_known_standard(self):
        intf = regtap.Interface("http://example.org", standard_id="ivo://ivoa.net/std/tap#aux",
                                intf_type="vs:paramhttp", intf_role="std")
        assert isinstance(intf.to_service(), tap.TAPService)
        assert not intf.is_vosi

    def test_sia2_standard(self):
        intf = regtap.Interface("http://example.org", standard_id="ivo://ivoa.net/std/sia2",
                                intf_type="vs:paramhttp", intf_role="std")
        assert isinstance(intf.to_service(), sia2.SIA2Service)
        assert not intf.is_vosi

    def test_secondary_interface(self):
        intf = regtap.Interface("http://example.org", standard_id="ivo://ivoa.net/std/tap#aux",
                                intf_type="vs:webbrowser", intf_role="web")

        with pytest.raises(ValueError) as excinfo:
            intf.to_service()

        assert str(excinfo.value) == (
            "This is not a standard interface.  PyVO cannot speak to it.")

    def test_VOSI(self):
        intf = regtap.Interface("http://example.org", standard_id="ivo://ivoa.net/std/vosi#capabilities",
                                intf_type="vs:ParamHTTP", intf_role="std")
        assert intf.is_vosi


# The following tests have their assertions in the fixtures.
# It would certainly not hurt to refactor this so they are
# in the tests (we could also just rely on the rtcons tests
# that exercise about the same thing).

@pytest.mark.usefixtures('keywords_fixture', 'capabilities')
def test_keywords():
    regsearch(keywords=['vizier', 'pulsar'])


@pytest.mark.usefixtures('single_keyword_fixture', 'capabilities')
def test_single_keyword():
    regsearch(keywords=['single'])
    regsearch(keywords='single')


@pytest.mark.usefixtures('servicetype_fixture', 'capabilities')
def test_servicetype():
    regsearch(servicetype='table')


@pytest.mark.usefixtures(
    'waveband_fixture',
    'capabilities',
    'messenger_vocabulary')
def test_waveband():
    regsearch(waveband='optical')


@pytest.mark.usefixtures('datamodel_fixture', 'capabilities')
def test_datamodel():
    regsearch(datamodel='ObsCore')


@pytest.mark.usefixtures('aux_fixture', 'capabilities')
def test_servicetype_aux():
    regsearch(servicetype='table', includeaux=True)


@pytest.mark.usefixtures('aux_fixture', 'capabilities')
def test_bad_servicetype_aux():
    with pytest.raises(dalq.DALQueryError):
        regsearch(servicetype='bad_servicetype', includeaux=True)


class _NS:
    """a namespace exposing its keyword arguments as attributes.

    We need this here to let us conveniently construct _FakeResults.
    """
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


class _FakeResults:
    """
    a minimal standin for dal.query.Results to be used with dal.query.Record.

    It is constructed with a dictionary that should eventually be
    used as the mapping in the Record.
    """
    def __init__(self, valdict):
        self.fieldnames = list(valdict.keys())
        self.resultstable = _NS(array=_NS(data=[list(valdict.values())]))


def get_regtap_results(**kwargs):
    """
    return a RegTAP result as expected by RegistryResult with all values
    empty, completed with what's in kwargs.
    """
    res = {}
    for key in regtap.RegistryResource.expected_columns:
        if isinstance(key, str):
            res[key] = None
        else:
            res[key[-1]] = None

    res.update(kwargs)
    return _FakeResults(res)


def test_spatial():
    assert (rtcons.keywords_to_constraints({
            "spatial": (23, -40)})[0].get_search_condition(FAKE_GAVO)
            == "ivoid IN (SELECT DISTINCT ivoid FROM rr.stc_spatial"
            " WHERE 1 = CONTAINS(MOC(6, POINT(23, -40)), coverage))")


def test_spectral():
    assert (rtcons.keywords_to_constraints({
            "spectral": (1e-17, 2e-17)})[0].get_search_condition(FAKE_GAVO)
            == "ivoid IN (SELECT DISTINCT ivoid FROM rr.stc_spectral WHERE"
            " 1 = ivo_interval_overlaps(spectral_start, spectral_end, 1e-17, 2e-17))")


def test_to_table(multi_interface_fixture, capabilities):
    t = regtap.search(
        ivoid="ivo://org.gavo.dc/flashheros/q/ssa").get_summary()
    assert (set(t.columns.keys())
            == {'index', 'short_name', 'title', 'description', 'interfaces'})
    assert t["index"][0] == 0
    assert t["title"][0] == 'Flash/Heros SSAP'
    assert (t["description"][0][:40]
            == 'Spectra from the Flash and Heros Echelle')
    assert (t["interfaces"][0]
            == 'datalink#links-1.1, soda#sync-1.0, ssa, tap#aux, web')


@pytest.fixture(name='rt_pulsar_distance')
def _rt_pulsar_distance(regtap_pulsar_distance_response, capabilities):
    return regsearch(keywords="pulsar", ucd=["pos.distance"])


def test_record_fields(rt_pulsar_distance):
    rec = rt_pulsar_distance["VII/156"]
    assert rec.ivoid == "ivo://cds.vizier/vii/156"
    assert rec.res_type == "vs:catalogservice"
    assert rec.short_name == "VII/156"
    assert rec.res_title == "Catalog of 558 Pulsars"
    assert rec.content_levels == ['research']
    assert rec.res_description[:20] == "The catalogue is an up-to-date"[:20]
    assert rec.reference_url == "https://cdsarc.cds.unistra.fr/viz-bin/cat/VII/156"
    assert rec.creators == ['Taylor J.H.', ' Manchester R.N.', ' Lyne A.G.']
    assert rec.content_types == ['catalog']
    assert rec.source_format == "bibcode"
    assert rec.source_value == "1993ApJS...88..529T"
    assert rec.region_of_regard is None
    assert rec.waveband == ['radio']
    assert rec.created == "1999-03-02T10:21:53"
    # updated might break when regenerating data/regtap.xml,
    # replace by the new date
    assert rec.updated == "2021-10-21T00:00:00"
    assert rec.rights == "https://cds.unistra.fr/vizier-org/licences_vizier.html"
    # access URL, standard_id and friends exercised in TestInterfaceSelection


class TestResultIndexing:
    def test_get_with_index(self, rt_pulsar_distance):
        # this is expecte to break when the fixture is updated
        assert (rt_pulsar_distance[0].res_title
                == 'Pulsar Timing for Fermi Gamma-ray Space Telescope')

    def test_get_with_short_name(self, rt_pulsar_distance):
        assert (rt_pulsar_distance["ATNF"].res_title
                == 'ATNF Pulsar Catalog')

    def test_get_with_ivoid(self, rt_pulsar_distance):
        assert (rt_pulsar_distance["ivo://nasa.heasarc/atnfpulsar"
                                   ].res_title == 'ATNF Pulsar Catalog')

    def test_out_of_range(self, rt_pulsar_distance):
        with pytest.raises(IndexError) as excinfo:
            rt_pulsar_distance[40320]
        assert (str(excinfo.value)
                == f"index 40320 is out of bounds for axis 0 with size {len(rt_pulsar_distance)}")

    def test_bad_key(self, rt_pulsar_distance):
        with pytest.raises(KeyError) as excinfo:
            rt_pulsar_distance["hunkatunka"]
        assert (str(excinfo.value) == "'hunkatunka'")

    def test_not_indexable(self, rt_pulsar_distance):
        with pytest.raises(IndexError) as excinfo:
            rt_pulsar_distance[None]
        assert (str(excinfo.value)
                == "No resource matching None")


@pytest.mark.usefixtures('multi_interface_fixture', 'capabilities',
                         'flash_service')
class TestInterfaceSelection:
    """
    tests for the selection and generation of services within
    RegistryResource.
    """

    def test_exactly_one_result(self):
        results = regtap.search(
            ivoid="ivo://org.gavo.dc/flashheros/q/ssa")
        assert len(results) == 1

    def test_access_modes(self, flash_service):
        assert set(flash_service.access_modes()) == {
            'datalink#links-1.1', 'soda#sync-1.0', 'ssa', 'tap#aux',
            'web'}

    def test_standard_id_multi(self, flash_service):
        with pytest.raises(dalq.DALQueryError) as excinfo:
            _ = flash_service.standard_id

        assert str(excinfo.value) == ("This resource supports several"
                                      " standards (datalink#links-1.1, soda#sync-1.0, ssa,"
                                      " tap#aux, web).  Use get_service or restrict your query"
                                      " using Servicetype.")

    def test_get_web_interface(self, flash_service):
        svc = flash_service.get_service(service_type="web")
        assert isinstance(svc,
                          regtap._BrowserService)
        assert (svc.baseurl
                == "http://dc.zah.uni-heidelberg.de/flashheros/q/web/form")
        assert str(svc) == ("BrowserService(baseurl : "
                            "'http://dc.zah.uni-heidelberg.de/flashheros/q/web/form',"
                            " description : '')")

        import webbrowser
        orig_open = webbrowser.open
        try:
            open_args = []
            webbrowser.open = lambda *args: open_args.append(args)
            svc.search()
            assert open_args == [
                ("http://dc.zah.uni-heidelberg.de/flashheros/q/web/form", 2)]
        finally:
            webbrowser.open = orig_open

    def test_get_aux_interface(self, flash_service):
        svc = flash_service.get_service(service_type="tap#aux")
        assert (svc._baseurl
                == "http://dc.zah.uni-heidelberg.de/tap")

    def test_get_aux_as_main(self, flash_service):
        assert (flash_service.get_service(service_type="tap")._baseurl
                == "http://dc.zah.uni-heidelberg.de/tap")

    def test_get__main_from_aux(self, flash_service):
        assert (flash_service.get_service(service_type="tap")._baseurl
                == "http://dc.zah.uni-heidelberg.de/tap")

    def test_get_by_alias(self, flash_service):
        assert (flash_service.get_service(service_type="spectrum")._baseurl
                == "http://dc.zah.uni-heidelberg.de/fhssa?")

    def test_get_unsupported_standard(self, flash_service):
        with pytest.raises(ValueError) as excinfo:
            flash_service.get_service(service_type="soda#sync-1.0")

        assert str(excinfo.value) == (
            "PyVO has no support for interfaces with standard id"
            " ivo://ivoa.net/std/soda#sync-1.0.")

    def test_get_nonexisting_standard(self, flash_service):
        with pytest.raises(ValueError) as excinfo:
            flash_service.get_service(service_type="http://nonsense#fancy")

        assert str(excinfo.value) == (
            "No matching interface.")

    def test_unconstrained(self, flash_service):
        with pytest.raises(ValueError, match="Multiple matching interfaces found.*"):
            flash_service.get_service(lax=False)

    def test_interface_without_role(self):
        # There's an ugly corner case in our array simulation for
        # capabilities and interfaces: if there's a single untyped
        # interface, the returned type (or role) will be an empty
        # string, and the split() will return an empty list.
        # This swallowed the interface in pyVO 1.3.
        rec = get_regtap_results(
            access_urls="http://example.org/tap",
            standard_ids="ivo://ivoa.net/std/TAP",
            intf_types="vr:webbrowser",
            intf_roles="")

        resource = regtap.RegistryResource(rec, 0)
        assert len(resource.interfaces) == 1
        assert resource.interfaces[0].standard_id == 'ivo://ivoa.net/std/TAP'

        # get_service still won't work because it needs a paramhttp
        # interface (and a role="std").
        with pytest.raises(ValueError) as excinfo:
            resource.get_service(service_type='tap')
        assert (str(excinfo.value) == "No matching interface.")

    def test_get_service_by_keyword(self):
        rec = _makeRegistryRecord(
            access_urls=["http://example1.org/tap", "http://example2.org/tap",
                         "http://example3.org/tap"],
            standard_ids=["ivo://ivoa.net/std/tap"] * 3,
            intf_types=["vs:paramhttp"] * 3,
            intf_roles=["std"] * 3,
            cap_descriptions=["Example 1 TAP", "Example 2 TAP"])  # no description for the third one
        service = rec.get_service(service_type="tap", keyword="Example 2")
        assert service.capability_description == "Example 2 TAP"

    def test_sia2_query(self):
        rec = _makeRegistryRecord(
            access_urls=["http://sia2.example.com", "http://sia.example.com"],
            standard_ids=[
                "ivo://ivoa.net/std/sia#query-2.0",
                "ivo://ivoa.net/std/sia"],
            intf_roles=["std"] * 2,
            intf_types=["vs:paramhttp"] * 2,
            cap_descriptions=["A mock SIA2 Service."] * 2)
        assert rec.access_modes() == {"sia", "sia2"}
        assert rec.get_interface(service_type="sia2").access_url == 'http://sia2.example.com'
        assert rec.get_interface(service_type="sia").access_url == 'http://sia.example.com'

    def test_sia2_aux(self):
        rec = _makeRegistryRecord(
            access_urls=["http://sia2.example.com", "http://sia.example.com"],
            standard_ids=[
                "ivo://ivoa.net/std/sia#query-aux-2.0",
                "ivo://ivoa.net/std/sia#aux"],
            intf_roles=["std"] * 2,
            intf_types=["vs:paramhttp"] * 2,
            cap_descriptions=["A mock service."] * 2)
        assert rec.access_modes() == {"sia#aux", "sia2#aux"}
        assert rec.get_interface(service_type="sia2").access_url == 'http://sia2.example.com'
        assert rec.get_interface(service_type="sia").access_url == 'http://sia.example.com'

    def test_non_standard_interface(self):
        intf = regtap.Interface("http://url", standard_id="", intf_type="", intf_role="")
        assert intf.supports("ivo://ivoa.net/std/sia") is False

    def test_supports_none(self):
        intf = regtap.Interface("http://url", standard_id="", intf_type="", intf_role="")
        assert intf.supports(None) is False

    def test_non_searchable_service(self):
        rec = _makeRegistryRecord()
        with pytest.raises(dalq.DALServiceError) as excinfo:
            rec.search()

        assert str(excinfo.value) == (
            "Resource ivo://pyvo/test_regtap.py is not a searchable service")

    def test_list_interfaces(self):
        rec = _makeRegistryRecord(
            access_urls=["http://sia.example.com", "http://sia.example.com",
                         "http://tap.example.com", "http://vosi.example.com",
                         "http://website.com"
                         ],
            standard_ids=["ivo://ivoa.net/std/sia#aux",
                          "ivo://ivoa.net/std/sia#aux",
                          "ivo://ivoa.net/std/tap",
                          "ivo://ivoa.net/std/vosi", ""],
            intf_roles=["std"] * 4 + ["non standard"],
            intf_types=["vs:paramhttp"] * 4 + ["vr:webbrowser"],
            cap_descriptions=["A mock service."] * 5)
        # webpages should be there, but not vosi
        assert len(rec.list_interfaces()) == 4
        # only sia ones
        assert len(rec.list_interfaces("sia")) == 2


class _FakeResult:
    """A fake class just sufficient for giving dal.query.Record enough
    to pull in the dict this is constructed.

    As an extra service, list values are stringified with
    regtap.TOKEN_SEP -- this is how they ought to come in from
    RegTAP services.
    """

    def __init__(self, d):
        self.fieldnames = list(d.keys())
        vals = [regtap.TOKEN_SEP.join(v) if isinstance(v, list) else v
                for v in d.values()]

        class _:
            class array:
                data = [vals]

        self.resultstable = _


def _makeRegistryRecord(**overrides):
    """returns a minimal RegistryResource instance, overriding
    some built-in defaults with the dict overrides.
    """
    defaults = {
        "access_urls": "",
        "standard_ids": "",
        "intf_types": "",
        "intf_roles": "",
        "cap_descriptions": "",
        "ivoid": "ivo://pyvo/test_regtap.py",
    }
    defaults.update(overrides)
    return regtap.RegistryResource(_FakeResult(defaults), 0)


class TestInterfaceRejection:
    """tests for various artificial corner cases where interface selection
    should fail (or just not fail).
    """

    def test_nonunique(self):
        rsc = _makeRegistryRecord(
            access_urls=["http://a", "http://b"],
            standard_ids=["ivo://ivoa.net/std/tap"] * 2,
            intf_types=["vs:paramhttp"] * 2,
            intf_roles=["std"] * 2)
        with pytest.raises(ValueError, match="Multiple matching interfaces found.*"):
            rsc.get_service(service_type="tap", lax=False)

    def test_nonunique_lax(self):
        rsc = _makeRegistryRecord(
            access_urls=["http://a", "http://b"],
            standard_ids=["ivo://ivoa.net/std/tap"] * 2,
            intf_types=["vs:paramhttp"] * 2,
            intf_roles=["std"] * 2,
            cap_descriptions=["Test TAP service."] * 2)

        assert (rsc.get_service(service_type="tap", lax=True)._baseurl
                == "http://a")

    def test_nonstd_ignored(self):
        rsc = _makeRegistryRecord(
            access_urls=["http://a", "http://b"],
            standard_ids=["ivo://ivoa.net/std/tap"] * 2,
            intf_types=["vs:paramhttp"] * 2,
            intf_roles=["std", ""])
        assert (rsc.get_service(service_type="tap", lax=False)._baseurl
                == "http://a")

    def test_select_single_matching_service(self):
        rsc = _makeRegistryRecord(
            access_urls=["http://a", "http://b"],
            standard_ids=["", "ivo://ivoa.net/std/tap"],
            intf_types=["vs:webbrowser", "vs:paramhttp"],
            intf_roles=["", "std"])

        assert (rsc.service._baseurl == "http://b")
        # this makes sure caching the service obtained doesn't break
        # things
        assert (rsc.service._baseurl == "http://b")

    def test_capless(self):
        rsc = _makeRegistryRecord()

        with pytest.raises(ValueError) as excinfo:
            rsc.service._baseurl

        assert str(excinfo.value) == (
            "No matching interface.")


@pytest.mark.remote_data
def test_maxrec():
    with pytest.warns(DALOverflowWarning) as w:
        _ = regsearch(servicetype="tap", maxrec=1)
    assert len(w) == 1
    assert str(w[0].message).startswith("Partial result set.")


@pytest.mark.remote_data
def test_get_contact():
    rsc = _makeRegistryRecord(
        ivoid="ivo://org.gavo.dc/flashheros/q/ssa")
    assert (rsc.get_contact()
            == "GAVO Data Centre Team (+49 6221 54 1837)"
            " <gavo@ari.uni-heidelberg.de>")


@pytest.mark.remote_data
def test_get_alt_identifier():
    rsc = _makeRegistryRecord(ivoid="ivo://cds.vizier/i/337")
    assert set(rsc.get_alt_identifiers()) == {
        'doi:10.26093/cds/vizier.1337',
        'bibcode:doi:10.5270/esa-ogmeula',
        'bibcode:2016yCat.1337....0G'}


@pytest.mark.remote_data
class TestDatamodelQueries:
    # right now, the data model queries are all rather sui generis, and
    # rather fickly on top.  Let's make sure they actually return something
    # against the live registry.  Admittedly, this is about as much
    # a test of the VO infrastructure as of pyvo.

    def test_obscore(self):
        assert len(regsearch(rtcons.Datamodel('obscore'))) > 0

    def test_epntap(self):
        assert len(regsearch(rtcons.Datamodel('epntap'))) > 0

    def test_regtap(self):
        assert len(regsearch(rtcons.Datamodel('regtap'))) > 0


@pytest.mark.usefixtures('multi_interface_fixture', 'capabilities',
                         'flash_service')
class TestExtraResourceMethods:
    """
    tests for methods of RegistryResource containing some non-trivial
    logic (except service selection, which is in TestInterfaceSelection,
    and get_tables, which is in TestGetTables).
    """

    def test_unique_standard_id(self):
        rsc = _makeRegistryRecord(
            access_urls=["http://a"],
            standard_ids=["ivo://ivoa.net/std/tap"],
            intf_types=["vs:paramhttp"],
            intf_roles=["std"])
        assert rsc.standard_id == "ivo://ivoa.net/std/tap"

    def test_describe_multi(self, flash_service):
        out = io.StringIO()
        flash_service.describe(verbose=True, file=out)
        output = out.getvalue()

        assert "Flash/Heros SSAP" in output
        assert ("Access modes: datalink#links-1.1, soda#sync-1.0,"
                " ssa, tap#aux, web" in output)
        assert "- webpage: http://dc.zah.uni-heidelberg.de/flashheros/q/web/form" in output
        assert "- tap#aux: http://dc.zah.uni-heidelberg.de/tap" in output
        # datalink, soda and vosi are not printed here
        assert "- datalink" not in output
        assert "Source: 1996A&A...312..539S" in output
        assert "Authors: Wolf" in output
        assert "Alternative identifier(s): doi:10.21938/" in output
        assert "More info: http://dc.zah" in output

    def test_describe_long_authors_list(self):
        """Check that long list of authors use et al.."""
        rsc = _makeRegistryRecord(
            access_urls=[],
            standard_ids=["ivo://pyvo/test"],
            short_name=["name"],
            intf_types=[],
            intf_roles=[],
            creator_seq=["a;" * 6],
            res_title=["title"]
        )
        out = io.StringIO()
        rsc.describe(verbose=True, file=out)
        output = out.getvalue()
        # output should cut at 5 authors
        assert "Authors: a, a, a, a, a et al." in output

    def test_describe_long_author_name(self):
        """Check that long author names are truncated."""
        rsc = _makeRegistryRecord(
            access_urls=[],
            standard_ids=["ivo://pyvo/test"],
            short_name=["name"],
            intf_types=[],
            intf_roles=[],
            creator_seq=["a" * 71],
            res_title=["title"]
        )
        out = io.StringIO()
        rsc.describe(verbose=True, file=out)
        output = out.getvalue()
        # should cut the long author name at 70 characters
        assert f"Authors: {'a' * 70}..." in output

    def test_no_access_url(self):
        rsc = _makeRegistryRecord(
            access_urls=[],
            standard_ids=[],
            intf_types=[],
            intf_roles=[])
        with pytest.raises(dalq.DALQueryError) as excinfo:
            rsc.access_url

        assert str(excinfo.value) == ("The resource"
                                      " ivo://pyvo/test_regtap.py has no queriable interfaces.")

    def test_unique_access_url(self):
        rsc = _makeRegistryRecord(
            access_urls=["http://a"],
            standard_ids=["ivo://ivoa.net/std/tap"],
            intf_types=["vs:paramhttp"],
            intf_roles=[""])
        assert rsc.access_url == "http://a"

    def test_ambiguous_access_url_warns(self, recwarn):
        rsc = _makeRegistryRecord(
            access_urls=["http://a", "http://b"],
            standard_ids=["ivo://ivoa.net/std/tap"] * 2,
            intf_types=["vs:paramhttp"] * 2,
            intf_roles=["std"] * 2)
        assert rsc.access_url == "http://a"
        assert ('The resource ivo://pyvo/test_regtap.py has multipl' in
                [str(w.message)[:50] for w in recwarn])


# TODO: While I suppose the contact test should keep requiring network,
# I think we should can the network responses involved in the following;
# the stuff might change upstream any time and then break our unit tests.
@pytest.fixture(name='obscore_tables')
def _obscore_tables(scope="session"):
    rsc = _makeRegistryRecord(
        ivoid="ivo://org.gavo.dc/__system__/obscore/obscore")
    return rsc.get_tables()


@pytest.mark.usefixtures("obscore_tables")
class TestGetTables:
    @pytest.mark.remote_data
    def test_get_tables_limit_enforced(self):
        rsc = _makeRegistryRecord(
            ivoid="ivo://org.gavo.dc/tap")
        with pytest.raises(dalq.DALQueryError) as excinfo:
            rsc.get_tables()

        assert re.match(r"Resource ivo://org.gavo.dc/tap reports \d+ tables."
                        "  Pass a higher table_limit to see them all.", str(excinfo.value))

    @pytest.mark.remote_data
    def test_get_tables_names(self, obscore_tables):
        assert (list(sorted(obscore_tables.keys()))
                == ["ivoa.obscore"])

    @pytest.mark.remote_data
    def test_get_tables_table_instance(self, obscore_tables):
        assert (obscore_tables["ivoa.obscore"].name
                == "ivoa.obscore")
        assert (obscore_tables["ivoa.obscore"].description[:42]
                == "The IVOA-defined obscore table, containing")
        assert (obscore_tables["ivoa.obscore"].title
                == "GAVO Data Center Obscore Table")

        assert (obscore_tables["ivoa.obscore"].origin.ivoid
                == "ivo://org.gavo.dc/__system__/obscore/obscore")

    @pytest.mark.remote_data
    def test_get_tables_column_meta(self, obscore_tables):
        def getcol(name):
            for col in obscore_tables["ivoa.obscore"].columns:
                if name == col.name:
                    return col
            raise KeyError(name)

        assert getcol("access_url").datatype.content == "char"
        assert getcol("access_url").datatype.arraysize == "*"

        assert getcol("s_region").datatype._extendedtype == "adql:region"

        assert getcol("access_format").ucd == 'meta.code.mime'

        assert getcol("em_min").unit == "m"

        assert (getcol("t_max").utype
                == "obscore:char.timeaxis.coverage.bounds.limits.stoptime")

        assert (getcol("t_exptime").description
                == "Total exposure time")

    @pytest.mark.remote_data
    def test_get_tables_utype(self, obscore_tables):
        assert (obscore_tables["ivoa.obscore"].utype
                == "ivo://ivoa.net/std/obscore#table-1.1")


@pytest.mark.remote_data
def test_sia2_service_operation():
    svcs = regsearch(
        servicetype='sia2',
        ivoid='ivo://cadc.nrc.ca/sia')
    assert len(svcs) == 1

    res = svcs[0].search(pos=(30, 40, 0.1),
        time=(time.Time(58794.9, format="mjd"),
            time.Time(58795, format="mjd")))
    assert len(res) > 10
    assert "s_dec" in res.to_table().columns


@pytest.mark.remote_data
def test_endpoint_switching():
    alt_svc = "http://vao.stsci.edu/RegTAP/TapService.aspx"
    previous_url = regtap.REGISTRY_BASEURL
    try:
        regtap.choose_RegTAP_service(alt_svc)
        assert regtap.get_RegTAP_service()._baseurl == alt_svc

        res = regtap.search(keywords="wirr")
        assert len(res) > 0
    finally:
        regtap.choose_RegTAP_service(previous_url)
