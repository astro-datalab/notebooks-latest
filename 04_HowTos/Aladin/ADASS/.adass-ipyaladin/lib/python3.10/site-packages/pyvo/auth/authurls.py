import collections
import logging

from . import securitymethods

__all__ = ["AuthURLs"]


class AuthURLs():
    """
    AuthURLs helps determine which security method should be used
    with a given URL.  It learns the security methods through the
    VOSI capabilities, which are passed in via update_from_capabilities.
    """

    def __init__(self):
        self.full_urls = collections.defaultdict(set)
        self.base_urls = collections.defaultdict(set)

    def update_from_capabilities(self, capabilities):
        """
        Update the URL to security method mapping using the
        capabilities provided.

        Parameters
        ----------
        capabilities : object
            List of `~pyvo.io.vosi.voresource.Capability`
        """
        for c in capabilities:
            for i in c.interfaces:
                for u in i.accessurls:
                    url = u.content
                    exact = u.use == 'full'

                    if not i.securitymethods:
                        self.add_security_method_for_url(url, securitymethods.ANONYMOUS, exact)

                    for sm in i.securitymethods:
                        method = sm.standardid or securitymethods.ANONYMOUS
                        self.add_security_method_for_url(url, method, exact)

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
        if exact:
            self.full_urls[url].add(security_method)
        else:
            self.base_urls[url].add(security_method)

    def allowed_auth_methods(self, url):
        """
        Return the authentication methods allowed for a particular URL.
        The methods are returned as URIs that represent security methods.

        Parameters
        ----------
        url : str
            the URL to determine authentication methods
        """
        logging.debug('Determining auth method for %s', url)

        if url in self.full_urls:
            methods = self.full_urls[url]
            logging.debug('Matching full url %s, methods %s', url, methods)
            return methods

        for base_url, methods in self._iterate_base_urls():
            if url.startswith(base_url):
                logging.debug('Matching base url %s, methods %s', base_url, methods)
                return methods

        logging.debug('No match, using anonymous auth')
        return {securitymethods.ANONYMOUS}

    def _iterate_base_urls(self):
        """
        A generator to sort the base URLs in the correct way
        to determine the most specific base_url.  This is done
        by returning them longest to shortest.
        """
        def sort_by_len(x):
            return len(x[0])

        # Sort the base path matching URLs, so that
        # the longest URLs (the most specific ones, if
        # there is a tie) are used to determine the
        # auth method.
        for url, method in sorted(self.base_urls.items(),
                                  key=sort_by_len,
                                  reverse=True):
            yield url, method

    def __repr__(self):
        urls = []

        for url, methods in self.full_urls.items():
            urls.append('Full match:' + url + ':' + str(methods))
        for url, methods in self._iterate_base_urls():
            urls.append('Base match:' + url + ':' + str(methods))

        return '\n'.join(urls)
