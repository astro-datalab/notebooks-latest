#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.mimetype
"""

from functools import partial
import pytest
import requests_mock

from astropy.utils.data import get_pkg_data_contents

from pyvo.dal.mimetype import mime_object_maker

mime_url = 'https://someurl.com/'

get_pkg_data_contents = partial(
    get_pkg_data_contents, package=__package__, encoding='binary')

try:
    from PIL import Image  # noqa: F401
    HAS_PILLOW = True
except ImportError:
    HAS_PILLOW = False


@pytest.fixture()
def mime(mocker):
    def callback(request, context):
        if 'mime-text' in request.url:
            return b'Text content'
        elif 'image' in request.url:
            return get_pkg_data_contents('data/mimetype/ivoa_logo.jpg')
        elif 'fits' in request.url:
            return get_pkg_data_contents('data/mimetype/test.fits')

    with mocker.register_uri(
        'GET', requests_mock.ANY, content=callback
    ) as matcher:
        yield matcher


@pytest.mark.usefixtures('mime')
@pytest.mark.skipif('not HAS_PILLOW')
def test_mime_object_maker():

    assert 'Text content' == mime_object_maker(mime_url + 'mime-text',
                                               'text/csv')

    img = mime_object_maker(mime_url + 'image', 'image/jpeg')
    assert img
    assert 'JPEG' == img.format

    fits = mime_object_maker(mime_url + 'fits', 'application/fits')
    assert 2 == len(fits)

    # error cases
    with pytest.raises(ValueError):
        mime_object_maker(None, "not/a/mime/type")
    with pytest.raises(ValueError):
        mime_object_maker(None, None)
