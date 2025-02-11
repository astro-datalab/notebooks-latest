# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
HTTP utils
"""
import platform
import requests
from ..version import version

DEFAULT_USER_AGENT = f'pyVO/{version} Python/{platform.python_version()} ({platform.system()})'


def use_session(session):
    """
    Return the session passed in, or create a default
    session to use for this network request.
    """
    if session:
        return session
    else:
        return create_session()


def create_session():
    """
    Create a new empty requests session with a pyvo
    user agent.
    """
    session = requests.Session()
    session.headers['User-Agent'] = DEFAULT_USER_AGENT
    return session
