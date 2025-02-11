#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.sia
"""
from functools import partial
import os
import re

import pytest

from pyvo.dal.sia import search, SIAService, SIAQuery

from astropy.io.fits import HDUList
from astropy.coordinates import SkyCoord
from astropy.utils.data import get_pkg_data_contents

get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')

sia_re = re.compile('http://example.com/sia.*')


@pytest.fixture()
def register_mocks(mocker):
    with mocker.register_uri(
        'GET', 'http://example.com/querydata/image.fits',
        content=get_pkg_data_contents('data/querydata/image.fits')
    ) as matcher:
        yield matcher


@pytest.fixture()
def sia(mocker):
    with mocker.register_uri(
        'GET', sia_re, content=get_pkg_data_contents('data/sia/dataset.xml')
    ) as matcher:
        yield matcher


def _test_result(result):
    assert result.getdataurl() == 'http://example.com/querydata/image.fits'
    assert isinstance(result.getdataobj(), HDUList)
    assert result.filesize == 153280


@pytest.mark.usefixtures('sia')
@pytest.mark.usefixtures('register_mocks')
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
@pytest.mark.parametrize("position", ((288, 15), SkyCoord(288, 15, unit="deg")))
@pytest.mark.parametrize("format", ("IMAGE/JPEG", "all"))
def test_search(position, format):
    results = search('http://example.com/sia', pos=position, format=format)
    result = results[0]

    _test_result(result)


class TestSIAService:
    @pytest.mark.usefixtures('sia')
    @pytest.mark.usefixtures('register_mocks')
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W42")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W49")
    def test_search(self):
        description = "A SIA service."
        url = 'http://example.com/sia'
        service = SIAService(url, capability_description=description)
        assert service.baseurl == url
        assert service.capability_description == description

        results = service.search(pos=(288, 15))
        result = results[0]

        _test_result(result)

        assert results[1].dateobs is None

    @pytest.mark.usefixtures('sia')
    @pytest.mark.usefixtures('register_mocks')
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W42")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W49")
    def test_formatter(self):
        service = SIAQuery('http://example.com/sia')
        service.format = "image"
        assert service["FORMAT"] == "image"
        service.format = "all"
        assert service["FORMAT"] == "ALL"
        service.format = "Graphic-png"
        assert service["FORMAT"] == "GRAPHIC-png"
        service.format = "Unsupported"
        assert service["FORMAT"] == "Unsupported"


@pytest.mark.usefixtures('sia')
class TestNameMaking:
    def test_slash_replace(self):
        res = SIAService('http://example.com/sia').search(pos=(288, 15))
        res["imageTitle"][0] = "ogsa/dai output"
        assert res[0].make_dataset_filename() == os.path.join(".", "ogsa_dai_output.fits")

    def test_default_for_broken_media_type(self):
        res = SIAService('http://example.com/sia').search(pos=(288, 15))
        res["mime"][0] = "application/x-youcannotknowthis"
        assert res[0].make_dataset_filename() == os.path.join(".", "Test_Observation.dat")

    def test_default_media_type_adaption(self):
        res = SIAService('http://example.com/sia').search(pos=(288, 15))
        res["mime"][0] = "image/png"
        assert res[0].make_dataset_filename() == os.path.join(".", "Test_Observation.png")
