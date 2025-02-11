#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.auth
"""
import base64
from requests.cookies import RequestsCookieJar
from string import Template

import pytest

from astropy.utils.data import get_pkg_data_contents

import pyvo.dal
from pyvo.auth.authsession import AuthSession

from pyvo.dal.tests.test_tap import MockAsyncTAPServer


@pytest.fixture(name='security_methods')
def _security_methods():
    return [None]


@pytest.fixture()
def auth_capabilities(mocker, security_methods):
    anon_element = '<securityMethod />'
    sm_element = Template('<securityMethod standardID="$m" />')
    method_elements = []

    for m in security_methods:
        if m is None:
            method_elements.append(anon_element)
        else:
            method_elements.append(sm_element.substitute(m=m))

    def callback(request, context):
        t = Template(get_pkg_data_contents('data/tap/capabilities.xml'))
        capabilities = t.substitute(security_methods=method_elements)
        return capabilities.encode('utf-8')

    with mocker.register_uri(
        'GET', 'http://example.com/tap/capabilities', content=callback
    ) as matcher:
        yield matcher


class MockAnonAuthTAPServer(MockAsyncTAPServer):
    def validator(self, request):
        assert request.cert is None
        assert 'Authorization' not in request.headers
        assert 'Cookie' not in request.headers


@pytest.fixture()
def anon_auth_service(mocker):
    yield from MockAnonAuthTAPServer().use(mocker)


class MockCookieAuthTAPServer(MockAsyncTAPServer):
    def validator(self, request):
        assert request.cert is None
        assert 'Authorization' not in request.headers
        assert request.headers.get('Cookie', None) == 'TEST_COOKIE=BADCOOKIE'


@pytest.fixture()
def cookie_auth_service(mocker):
    yield from MockCookieAuthTAPServer().use(mocker)


class MockCertificateAuthTAPServer(MockAsyncTAPServer):
    def validator(self, request):
        assert request.cert == 'client-certificate.pem'
        assert 'Authorization' not in request.headers
        assert 'Cookie' not in request.headers


@pytest.fixture()
def certificate_auth_service(mocker):
    yield from MockCertificateAuthTAPServer().use(mocker)


class MockBasicAuthTAPServer(MockAsyncTAPServer):
    def validator(self, request):
        pw = 'testuser:hunter2'.encode('ascii')
        basic_encoded = 'Basic ' + base64.b64encode(pw).decode('ascii')

        assert request.cert is None
        assert request.headers.get('Authorization', None) == basic_encoded
        assert 'Cookie' not in request.headers


@pytest.fixture()
def basic_auth_service(mocker):
    yield from MockBasicAuthTAPServer().use(mocker)


@pytest.mark.usefixtures('cookie_auth_service', 'auth_capabilities')
@pytest.mark.parametrize('security_methods', [[None, 'ivo://ivoa.net/sso#cookie']])
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_cookies_auth():
    session = AuthSession()
    session.credentials.set_cookie('TEST_COOKIE', 'BADCOOKIE')
    service = pyvo.dal.TAPService('http://example.com/tap', session=session)
    service.run_async("SELECT * FROM ivoa.obscore")


@pytest.mark.usefixtures('cookie_auth_service', 'auth_capabilities')
@pytest.mark.parametrize('security_methods', [[None, 'ivo://ivoa.net/sso#cookie']])
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_cookie_jar_auth():
    session = AuthSession()
    jar = RequestsCookieJar()
    jar.set('TEST_COOKIE', 'BADCOOKIE')
    session.credentials.set_cookie_jar(jar)
    service = pyvo.dal.TAPService('http://example.com/tap', session=session)
    service.run_async("SELECT * FROM ivoa.obscore")


@pytest.mark.usefixtures('certificate_auth_service', 'auth_capabilities')
@pytest.mark.parametrize('security_methods', [[None, 'ivo://ivoa.net/sso#tls-with-certificate']])
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_certificate_auth():
    session = AuthSession()
    session.credentials.set_client_certificate('client-certificate.pem')
    service = pyvo.dal.TAPService('http://example.com/tap', session=session)
    service.run_async("SELECT * FROM ivoa.obscore")


@pytest.mark.usefixtures('basic_auth_service', 'auth_capabilities')
@pytest.mark.parametrize('security_methods', [[None, 'ivo://ivoa.net/sso#BasicAA']])
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_basic_auth():
    session = AuthSession()
    session.credentials.set_password('testuser', 'hunter2')
    service = pyvo.dal.TAPService('http://example.com/tap', session=session)
    service.run_async("SELECT * FROM ivoa.obscore")


@pytest.mark.usefixtures('basic_auth_service', 'auth_capabilities')
@pytest.mark.parametrize('security_methods', [[None, 'ivo://ivoa.net/sso#BasicAA']])
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_negotiation():
    session = AuthSession()
    session.credentials.set_password('testuser', 'hunter2')
    session.credentials.set_client_certificate('client-certificate.pem')
    session.credentials.set_cookie('TEST_COOKIE', 'BADCOOKIE')
    service = pyvo.dal.TAPService('http://example.com/tap', session=session)
    service.run_async("SELECT * FROM ivoa.obscore")


@pytest.mark.usefixtures('anon_auth_service', 'auth_capabilities')
@pytest.mark.parametrize('security_methods', [[None, 'ivo://ivoa.net/sso#FancyAuth']])
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W27")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W48")
@pytest.mark.filterwarnings("ignore::astropy.io.votable.exceptions.W06")
def test_no_common_auth_negotiation():
    session = AuthSession()
    session.credentials.set_password('testuser', 'hunter2')
    session.credentials.set_client_certificate('client-certificate.pem')
    session.credentials.set_cookie('TEST_COOKIE', 'BADCOOKIE')
    service = pyvo.dal.TAPService('http://example.com/tap', session=session)
    service.run_async("SELECT * FROM ivoa.obscore")
