#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.datalink
"""
from functools import partial
import re

import pytest

import pyvo as vo
from pyvo.dal.adhoc import DatalinkResults, DALServiceError
from pyvo.dal.sia2 import SIA2Results
from pyvo.dal.tap import TAPResults
from pyvo.utils import testing, vocabularies
from pyvo.dal.sia import search

from astropy.utils.data import get_pkg_data_contents, get_pkg_data_filename

get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')


@pytest.fixture()
def ssa_datalink(mocker):
    def callback(request, context):
        return get_pkg_data_contents('data/datalink/datalink-ssa.xml')

    with mocker.register_uri(
        'GET', 'http://example.com/ssa_datalink', content=callback
    ) as matcher:
        yield matcher


sia_re = re.compile("http://example.com/sia.*")


@pytest.fixture()
def register_mocks(mocker):
    with mocker.register_uri(
        "GET",
        "http://example.com/querydata/image.fits",
        content=get_pkg_data_contents("data/querydata/image.fits"),
    ) as matcher:
        yield matcher


@pytest.fixture()
def sia(mocker):
    with mocker.register_uri(
        "GET", sia_re, content=get_pkg_data_contents("data/sia/dataset.xml")
    ) as matcher:
        yield matcher


@pytest.fixture()
def datalink(mocker):
    def callback(request, context):
        return get_pkg_data_contents('data/datalink/datalink.xml')

    with mocker.register_uri(
        'POST', 'http://example.com/datalink', content=callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def datalink_product(mocker):
    def callback(request, context):
        return get_pkg_data_contents('data/datalink/datalink.xml')

    with mocker.register_uri(
        'GET', 'http://example.com/datalink.xml', content=callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def obscore_datalink(mocker):
    def callback(request, context):
        return get_pkg_data_contents('data/datalink/datalink-obscore.xml')

    with mocker.register_uri(
        'GET', 'http://example.com/obscore', content=callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def res_datalink(mocker):

    first_batch = True

    def callback(request, context):
        nonlocal first_batch
        if first_batch:
            first_batch = False
            return get_pkg_data_contents('data/datalink/cutout1.xml')
        else:
            return get_pkg_data_contents('data/datalink/cutout2.xml')

    with mocker.register_uri(
        'POST', 'https://example.com/obscore-datalink', content=callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def proc(mocker):
    def callback(request, context):
        return get_pkg_data_contents('data/datalink/proc.xml')

    with mocker.register_uri(
        'GET', 'http://example.com/proc', content=callback
    ) as matcher:
        yield matcher


@pytest.fixture()
def datalink_vocabulary(mocker):
    # astropy download_file (which get_vocabluary uses) does not use
    # requests, so we can't mock this as we can mock the others.  We
    # replace the entire function for a while
    dl_voc_uri = 'http://www.ivoa.net/rdf/datalink/core'

    def fake_download_file(src_url, *args, **kwargs):
        assert src_url == dl_voc_uri
        return get_pkg_data_filename('data/datalink/datalink.desise')

    real_download_file = vocabularies.download_file
    try:
        vocabularies.download_file = fake_download_file
        yield
    finally:
        vocabularies.download_file = real_download_file


@pytest.mark.usefixtures('ssa_datalink', 'datalink')
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.E02")
def test_datalink():
    results = vo.spectrumsearch(
        'http://example.com/ssa_datalink', (30, 30))

    datalinks = next(results.iter_datalinks())

    assert datalinks.original_row["accsize"] == 100800

    row = datalinks[0]
    assert row.semantics == "#progenitor"

    row = datalinks[1]
    assert row.semantics == "#proc"

    row = datalinks[2]
    assert row.semantics == "#this"

    row = datalinks[3]
    assert row.semantics == "#preview"


@pytest.mark.usefixtures('obscore_datalink', 'res_datalink')
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.E02")
def test_datalink_batch():
    results = vo.dal.imagesearch(
        'http://example.com/obscore', (30, 30))

    dls = list(results.iter_datalinks())
    assert len(dls) == 3
    assert dls[0].original_row["obs_collection"] == "MACHO"


@pytest.mark.usefixtures('proc', 'datalink_vocabulary')
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.E02")
class TestSemanticsRetrieval:
    def test_access_with_string(self):
        datalinks = DatalinkResults.from_result_url('http://example.com/proc')

        assert datalinks.original_row is None
        res = [r["access_url"] for r in datalinks.bysemantics("#this")]
        assert len(res) == 1
        assert res[0].endswith("eq010000ms/20100927.comb_avg.0001.fits.fz")

    def test_access_with_list(self):
        datalinks = DatalinkResults.from_result_url('http://example.com/proc')
        res = [r["access_url"]
               for r in datalinks.bysemantics(["#this", "#preview-image"])]
        assert len(res) == 2
        assert res[0].endswith("eq010000ms/20100927.comb_avg.0001.fits.fz")
        assert res[1].endswith("20100927.comb_avg.0001.fits.fz?preview=True")

    def test_access_with_expansion(self):
        datalinks = DatalinkResults.from_result_url('http://example.com/proc')
        res = [r["access_url"]
               for r in datalinks.bysemantics(["#this", "#preview"])]
        assert len(res) == 3
        assert res[0].endswith("eq010000ms/20100927.comb_avg.0001.fits.fz")
        assert res[1].endswith("20100927.comb_avg.0001.fits.fz?preview=True")
        assert res[2].endswith("http://dc.zah.uni-heidelberg.de/wider.dat")

    def test_access_without_expansion(self):
        datalinks = DatalinkResults.from_result_url('http://example.com/proc')
        res = [r["access_url"] for r in datalinks.bysemantics(
            ["#this", "#preview"], include_narrower=False)]
        assert len(res) == 2
        assert res[0].endswith("eq010000ms/20100927.comb_avg.0001.fits.fz")
        assert res[1].endswith("http://dc.zah.uni-heidelberg.de/wider.dat")

    def test_with_full_url(self):
        datalinks = DatalinkResults.from_result_url('http://example.com/proc')
        res = [r["access_url"]
               for r in datalinks.bysemantics("urn:example:rdf/dlext#oracle")]
        assert len(res) == 1
        assert res[0].endswith("when-will-it-be-back")

    def test_all_mixed(self):
        datalinks = DatalinkResults.from_result_url('http://example.com/proc')
        res = [r["access_url"]
               for r in datalinks.bysemantics([
                   "urn:example:rdf/dlext#oracle",
                   'http://www.ivoa.net/rdf/datalink/core#preview',
                   '#this',
                   'non-existing-term'])]
        assert len(res) == 4
        assert res[0].endswith("eq010000ms/20100927.comb_avg.0001.fits.fz")
        assert res[1].endswith("comb_avg.0001.fits.fz?preview=True")
        assert res[2].endswith("http://dc.zah.uni-heidelberg.de/wider.dat")
        assert res[3].endswith("when-will-it-be-back")


@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.E02")
@pytest.mark.usefixtures('datalink_product', 'datalink_vocabulary')
class TestIterDatalinksProducts:
    """Tests for producing datalinks from tables containing links to
    datalink documents.
    """
    def test_no_access_format(self):
        res = testing.create_dalresults([
            {"name": "access_url", "datatype": "char", "arraysize": "*",
                "utype": "obscore:access.reference"}],
            [("http://foo.bar/baz.jpeg",)],
            resultsClass=TAPResults)
        assert list(res.iter_datalinks()) == []

    def test_obscore_utype(self):
        res = testing.create_dalresults([
            {"name": "data_product", "datatype": "char", "arraysize": "*",
                "utype": "obscore:access.reference"},
            {"name": "content_type", "datatype": "char", "arraysize": "*",
                "utype": "obscore:access.format"},],
            [("http://example.com/datalink.xml",
                "application/x-votable+xml;content=datalink")],
            resultsClass=TAPResults)
        links = list(res.iter_datalinks())
        assert len(links) == 1
        assert (next(links[0].bysemantics("#this"))["access_url"]
            == "http://dc.zah.uni-heidelberg.de/getproduct/flashheros/data/ca90/f0011.mt")

    def test_sia2_record(self):
        res = testing.create_dalresults([
            {"name": "access_url", "datatype": "char", "arraysize": "*",
                "utype": "obscore:access.reference"},
            {"name": "access_format", "datatype": "char", "arraysize": "*",
                "utype": "obscore:access.format"},],
            [("http://example.com/datalink.xml",
                "application/x-votable+xml;content=datalink")],
            resultsClass=SIA2Results)
        links = list(res.iter_datalinks())
        assert len(links) == 1
        assert (next(links[0].bysemantics("#this"))["access_url"]
            == "http://dc.zah.uni-heidelberg.de/getproduct/flashheros/data/ca90/f0011.mt")

    def test_sia1_record(self):
        res = testing.create_dalresults([
            {"name": "product", "datatype": "char", "arraysize": "*",
                "ucd": "VOX:Image_AccessReference"},
            {"name": "mime", "datatype": "char", "arraysize": "*",
                "ucd": "VOX:Image_Format"},],
            [("http://example.com/datalink.xml",
                "application/x-votable+xml;content=datalink")],
            resultsClass=TAPResults)
        links = list(res.iter_datalinks())
        assert len(links) == 1
        assert (next(links[0].bysemantics("#this"))["access_url"]
            == "http://dc.zah.uni-heidelberg.de/getproduct/flashheros/data/ca90/f0011.mt")

    def test_ssap_record(self):
        res = testing.create_dalresults([
            {"name": "product", "datatype": "char", "arraysize": "*",
                "utype": "ssa:access.reference"},
            {"name": "mime", "datatype": "char", "arraysize": "*",
                "utype": "ssa:access.format"},],
            [("http://example.com/datalink.xml",
                "application/x-votable+xml;content=datalink")],
            resultsClass=TAPResults)
        links = list(res.iter_datalinks())
        assert len(links) == 1
        assert (next(links[0].bysemantics("#this"))["access_url"]
            == "http://dc.zah.uni-heidelberg.de/getproduct/flashheros/data/ca90/f0011.mt")

    def test_generic_record(self):
        # The meta.code.mime and meta.ref.url UCDs are perhaps too
        # generic.  To ensure a somewhat predictable behaviour,
        # we at least make sure we pick the first of possibly multiple
        # pairs (not that this would preclude arbitrary amounts of
        # chaos).
        res = testing.create_dalresults([
            {"name": "access_url", "datatype": "char", "arraysize": "*",
                "ucd": "meta.ref.url"},
            {"name": "access_format", "datatype": "char", "arraysize": "*",
                "utype": "meta.code.mime"},
            {"name": "alt_access_url", "datatype": "char", "arraysize": "*",
                "ucd": "meta.ref.url"},
            {"name": "alt_access_format", "datatype": "char", "arraysize": "*",
                "utype": "meta.code.mime"},],
            [("http://example.com/datalink.xml",
                "application/x-votable+xml;content=datalink",
                "http://example.com/bad-pick.xml",
                "application/x-votable+xml;content=datalink",)],
            resultsClass=TAPResults)
        links = list(res.iter_datalinks())
        assert len(links) == 1
        assert (next(links[0].bysemantics("#this"))["access_url"]
            == "http://dc.zah.uni-heidelberg.de/getproduct/flashheros/data/ca90/f0011.mt")


@pytest.mark.usefixtures("sia", "register_mocks")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_no_datalink():
    # for issue #328 getdatalink() exits messily when there isn't a datalink

    results = search("http://example.com/sia", pos=(288, 15), format="all")
    result = results[0]
    with pytest.raises(DALServiceError, match="No datalink found for record."):
        result.getdatalink()
