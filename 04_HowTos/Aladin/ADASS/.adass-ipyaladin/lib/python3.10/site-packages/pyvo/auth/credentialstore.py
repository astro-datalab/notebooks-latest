import logging

from pyvo.utils.http import create_session

from . import securitymethods

__all__ = ["CredentialStore"]


class CredentialStore(object):
    """
    The credential store takes user credentials, and uses them
    to create appropriate requests sessions for dispatching
    requests using those credentials.

    Different types of credentials can be passed in, such as
    cookies, a jar of cookies, certificates, and basic auth.

    A session can also be associated with a security method
    URI by calling the set function.

    Before a request is to be dispatched, the AuthSession
    calls the get method to retrieve the appropriate
    requests.Session for making that HTTP request.
    """

    def __init__(self):
        self.credentials = {}
        self.set(securitymethods.ANONYMOUS, create_session())

    def negotiate_method(self, allowed_methods):
        """
        Compare the credentials provided by the user against the
        security methods passed in, and determine which method is
        to be used for making this request.

        Parameters
        ----------
        allowed_methods : list(str)
            list of allowed security methods to return

        Raises
        ------
        Raises an exception if a common method could not be negotiated.
        """
        available_methods = set(self.credentials.keys())
        methods = available_methods.intersection(allowed_methods)
        logging.debug('Available methods: %s', methods)

        # If we have no common auth mechanism, then fail.
        if not methods:
            msg = 'Negotiation failed.  Server supports %s, client supports %s' % \
                (allowed_methods, available_methods)
            raise Exception(msg)

        # If there are more than 1 method to pick, don't pick
        # anonymous over an actual method.
        if len(methods) > 1:
            methods.discard(securitymethods.ANONYMOUS)

        # Pick a random method.
        return methods.pop()

    def set(self, method_uri, session):
        """
        Associate a security method URI with a requests.Session like object.

        Parameters
        ----------
        method_uri : str
            URI representing the security method
        session : object
            the requests.Session like object that will dispatch requests
            for the authentication method provided by method_uri
        """
        self.credentials[method_uri] = session

    def get(self, method_uri):
        """
        Retrieve the requests.Session like object associated with a security
        method URI.

        Parameters
        ----------
        method_uri : str
            URI representing the security method
        """
        return self.credentials[method_uri]

    def set_cookie(self, cookie_name, cookie_value, domain='', path='/'):
        """
        Add a cookie to use as authentication.

        More than one call to set_cookie will add multiple cookies into
        the same cookie jar used for the request.

        Parameters
        ----------
        cookie_name : str
            name of the cookie
        cookie_value : str
            value of the cookie
        domain : str
            restrict usage of this cookie to this domain
        path : str
            restrict usage of this cookie to this path
        """
        cookie_session = self.credentials.setdefault(securitymethods.COOKIE, create_session())
        cookie_session.cookies.set(cookie_name, cookie_value, domain=domain, path=path)

    def set_cookie_jar(self, cookie_jar):
        """
        Set the cookie jar to use for authentication.

        Any previous cookies or cookie jars set will be removed.

        Parameters
        ----------
        cookie_jar : obj
            the cookie jar to use.
        """
        cookie_session = self.credentials.setdefault(securitymethods.COOKIE, create_session())
        cookie_session.cookies = cookie_jar

    def set_client_certificate(self, certificate_path):
        """
        Add a client certificate to use for authentication.

        Parameters
        ----------
        certificate_path : str
            path to the file of the client certificate
        """
        cert_session = create_session()
        cert_session.cert = certificate_path
        self.set(securitymethods.CLIENT_CERTIFICATE, cert_session)

    def set_password(self, username, password):
        """
        Add a username / password for basic authentication.

        Parameters
        ----------
        username : str
            username to use
        password : str
            password to use
        """
        basic_session = create_session()
        basic_session.auth = (username, password)
        self.set(securitymethods.BASIC, basic_session)

    def __repr__(self):
        return 'Support for %s' % self.credentials.keys()
