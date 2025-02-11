#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.scs
"""
from functools import partial
import re

import pytest

from pyvo.dal.scs import search, SCSService

from astropy.utils.data import get_pkg_data_contents

get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')

scs_re = re.compile('http://example.com/scs.*')


@pytest.fixture()
def scs(mocker):
    def callback(request, context):
        return get_pkg_data_contents('data/scs/result.xml')

    with mocker.register_uri(
        'GET', scs_re, content=callback
    ) as matcher:
        yield matcher


@pytest.mark.usefixtures('scs')
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_search():
    results = search('http://example.com/scs', pos=(78, 2), radius=0.5)

    assert len(results) == 1273


class TestSCSService:
    def test_init(self):
        service = SCSService('http://example.com/scs', capability_description="SCS")

        assert service.baseurl == 'http://example.com/scs'
        assert service.capability_description == "SCS"

    @pytest.mark.usefixtures('scs')
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
    def test_search(self):
        service = SCSService('http://example.com/scs')

        results = service.search(pos=(78, 2), radius=0.5)

        assert len(results) == 1273
