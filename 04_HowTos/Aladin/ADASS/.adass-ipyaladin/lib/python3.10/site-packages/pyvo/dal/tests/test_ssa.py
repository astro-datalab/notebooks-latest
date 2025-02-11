#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.ssa
"""
from functools import partial
import re

import pytest

from pyvo.dal.ssa import search, SSAService

from astropy.utils.data import get_pkg_data_contents

get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')

ssa_re = re.compile('http://example.com/ssa.*')


@pytest.fixture()
def ssa(mocker):
    def callback(request, context):
        return get_pkg_data_contents('data/ssa/result.xml')

    with mocker.register_uri(
        'GET', ssa_re, content=callback
    ) as matcher:
        yield matcher


@pytest.mark.usefixtures('ssa')
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W42")
def test_search():
    results = search('http://example.com/ssa', pos=(0.0, 0.0), diameter=1.0)
    assert len(results) == 36


class TestSSAService:
    @pytest.mark.usefixtures('ssa')
    @pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W42")
    def test_search(self):
        service = SSAService('http://example.com/ssa')

        assert str(service) == "SSAService(baseurl : 'http://example.com/ssa', description : 'None')"

        results = service.search(pos=(0.0, 0.0), diameter=1.0)

        assert len(results) == 36
        assert results[35].dateobs is None
