# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.utils.http
"""

import platform

from pyvo.utils.http import create_session
from pyvo.version import version


def test_create_session():
    test_session = create_session()
    assert (test_session.headers['User-Agent']
            == f'pyVO/{version} Python/{platform.python_version()} ({platform.system()})')
