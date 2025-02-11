import logging

from .authurls import AuthURLs
from .credentialstore import CredentialStore

__all__ = ["AuthSession"]


class AuthSession:
    """
    A requests-like session that pyvo can use to dispatch its
    network calls with authentication.

    The user adds their credentials to the credentials object,
    such as adding a cookie, certificate, or password.

    The network requests made by pyvo pass through here, and
    the URL of the request is matched against the capabilities
    of the service.  Based on what credentials have been
    provided and the capabilities of the service, appropriate
    credentials are added to the request before it is sent.
    """

    def __init__(self):
        super(AuthSession, self).__init__()
        self.credentials = CredentialStore()
        self._auth_urls = AuthURLs()

    def add_security_method_for_url(self, url, security_method, exact=False):
        """
        Add a security method for a url.
        This is additive with update_from_capabilities.  This
        can be useful to set additional security methods that
        aren't set in the capabilities for whatever reason.

        Parameters
        ----------
        url : str
            URL to set a security method for
        security_method : str
            URI of the security method to set
        exact : bool
            If True, match only this URL.  If false, match all URLs that
            match this as a base URL.
        """
        self._auth_urls.add_security_method_for_url(url, security_method, exact=exact)

    def update_from_capabilities(self, capabilities):
        """
        Update the URL to security method mapping using the
        capabilities provided.

        Parameters
        ----------
        capabilities : object
            List of `~pyvo.io.vosi.voresource.Capability`
        """
        self._auth_urls.update_from_capabilities(capabilities)

    def get(self, url, **kwargs):
        """
        Wrapper to make a HTTP GET request with authentication.
        """
        return self._request('GET', url, **kwargs)

    def post(self, url, **kwargs):
        """
        Wrapper to make a HTTP POST request with authentication.
        """
        return self._request('POST', url, **kwargs)

    def put(self, url, **kwargs):
        """
        Wrapper to make a HTTP PUT request with authentication.
        """
        return self._request('PUT', url, **kwargs)

    def delete(self, url, **kwargs):
        """
        Wrapper to make a HTTP DELETE request with authentication.
        """
        return self._request('DELETE', url, **kwargs)

    def _request(self, http_method, url, **kwargs):
        """
        Make an HTTP request with authentication.

        This function looks at the url of the request, determines
        what credentials it should attach to the request to
        authenticate, and then dispatches the request to the
        underlying requests library using the session that
        has been configured with the credentials.

        Parameters
        ----------
        http_method : str
            the HTTP verb of the request.
        url : str
            the URL to request
        """
        auth_methods = self._auth_urls.allowed_auth_methods(url)
        logging.debug('Possible auth methods: %s', auth_methods)

        negotiated_method = self.credentials.negotiate_method(auth_methods)
        logging.debug('Using auth method: %s', negotiated_method)

        session = self.credentials.get(negotiated_method)
        return session.request(http_method, url, **kwargs)

    def __repr__(self):
        return '\n'.join([repr(self.credentials), repr(self._auth_urls)])
