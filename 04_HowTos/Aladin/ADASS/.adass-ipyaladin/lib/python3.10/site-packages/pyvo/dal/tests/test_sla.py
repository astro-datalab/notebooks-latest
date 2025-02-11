#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.sla
"""
from functools import partial
import re

import pytest

from pyvo.dal.sla import search, SLAService

from astropy.utils.data import get_pkg_data_contents

get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')

sla_re = re.compile('http://example.com/sla.*')


@pytest.fixture()
def sla(mocker):
    with mocker.register_uri(
        'GET', sla_re, content=get_pkg_data_contents('data/sla/dataset.xml')
    ) as matcher:
        yield matcher


@pytest.mark.usefixtures('sla')
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W42")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W49")
def test_search():
    results = search('http://example.com/sla', wavelength=(7.6e-6, 1.e-5))

    assert len(results) == 21


class TestSLAService:

    def test_service_init(self):
        url = "http://test"
        description = "A SLA service."
        service = SLAService(url, capability_description=description)
        assert service.baseurl == url
        assert service.capability_description == description
        assert str(service) == f"SLAService(baseurl : '{url}', description : '{description}')"

    @pytest.mark.usefixtures('sla')
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W42")
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W49")
    def test_search(self):
        service = SLAService('http://example.com/sla')

        results = service.search(wavelength=(7.6e-6, 1.e-5))

        assert len(results) == 21
