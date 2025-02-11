"""
DAL Exceptions.
"""

__all__ = [
    "DALAccessError", "DALProtocolError", "DALFormatError", "DALServiceError",
    "DALQueryError", "DALOverflowWarning"]

import re

import requests

from astropy.utils.exceptions import AstropyUserWarning


class DALAccessError(Exception):
    """
    a base class for failures while accessing a DAL service
    """
    _defreason = "Unknown service access error"

    def __init__(self, reason=None, url=None):
        """
        initialize the exception with an error message

        Parameters
        ----------
        reason : str
           a message describing the cause of the error
        url : str
           the query URL that produced the error
        """
        if not reason:
            reason = self._defreason
        super().__init__(reason)
        self._reason = reason
        self._url = url

    @classmethod
    def _typeName(cls, exc):
        try:
            return exc.__qualname__
        except AttributeError:
            return re.sub(
                r"'>$", '',
                re.sub(r"<(type|class) '(.*\.)?", '', str(type(exc)))
            )

    def __str__(self):
        return self._reason

    def __repr__(self):
        return "{}: {}".format(self._typeName(self), self._reason)

    @property
    def reason(self):
        """
        a string description of what went wrong
        """
        return self._reason

    @property
    def url(self):
        """
        the URL that produced the error.  If None, the URL is unknown or unset
        """
        return self._url


class DALProtocolError(DALAccessError):
    """
    a base exception indicating that a DAL service responded with an error.
    This can be either an HTTP protocol error or a response format error;
    both of these are handled by separate subclasses. This base class captures
    an underlying exception clause.
    """
    _defreason = "Unknown DAL Protocol Error"

    def __init__(self, reason=None, cause=None, url=None):
        """
        initialize with a string message and an optional HTTP response code

        Parameters
        ----------
        reason : str
           a message describing the cause of the error
        code : int
           the HTTP error code (as an integer)
        cause : str
           an exception issued as the underlying cause.  A value
           of None indicates that no underlying exception was
           caught.
        url : str
           the query URL that produced the error
        """
        super().__init__(reason, url)
        self._cause = cause

    @property
    def cause(self):
        """
        a string description of what went wrong
        """
        return self._cause


class DALFormatError(DALProtocolError):
    """
    an exception indicating that a DAL response contains fatal format errors.
    This would include XML or VOTable format errors.
    """
    _defreason = "Unknown VOTable Format Error"

    def __init__(self, cause=None, url=None, reason=None):
        """
        create the exception

        Parameters
        ----------
        cause : str
           an exception issued as the underlying cause.  A value
           of None indicates that no underlying exception was caught.
        url
           the query URL that produced the error
        reason
           a message describing the cause of the error
        """
        if cause and not reason:
            reason = "{}: {}".format(
                DALAccessError._typeName(cause), str(cause))

        super().__init__(reason, cause, url)


class DALServiceError(DALProtocolError):
    """
    an exception indicating a failure communicating with a DAL
    service.  Most typically, this is used to report DAL queries that result
    in an HTTP error.
    """
    _defreason = "Unknown service error"

    def __init__(self, reason=None, code=None, cause=None, url=None):
        """
        initialize with a string message and an optional HTTP response code

        Parameters
        ----------
        reason : str
           a message describing the cause of the error
        code : int
           the HTTP error code (as an integer)
        cause : str
           an exception issued as the underlying cause.  A value
           of None indicates that no underlying exception was
           caught.
        url : str
           the query URL that produced the error
        """
        super().__init__(reason, cause, url)
        self._code = code

    @property
    def code(self):
        """
        the HTTP error code that resulted from the DAL service query,
        indicating the error.  If None, the service did not produce an HTTP
        response.
        """
        return self._code

    @classmethod
    def from_except(cls, exc, url=None):
        """
        create and return DALServiceError exception appropriate
        for the given exception that represents the underlying cause.
        """
        if isinstance(exc, requests.exceptions.RequestException):
            try:
                response = exc.response
            except AttributeError:
                response = None
            code = 0
            message = str(exc)

            # if there is a response, refine the error message
            if response is not None:
                code = response.status_code
                content_type = response.headers.get('content-type', None)
                if content_type and 'text/plain' in content_type:
                    message = '{} for {}'.format(response.text, url)

            # TODO votable handling

            return DALServiceError(message, code, exc, url)
        elif isinstance(exc, Exception):
            return DALServiceError("{}: {}".format(cls._typeName(exc), str(exc)),
                                   cause=exc, url=url)
        else:
            raise TypeError("from_except: expected Exception")


class DALQueryError(DALAccessError):
    """
    an exception indicating an error by a working DAL service while processing
    a query.  Generally, this would be an error that the service successfully
    detected and consequently was able to respond with a legal error response--
    namely, a VOTable document with an INFO element contains the description
    of the error.  Possible errors will include bad usage by the client, such
    as query-syntax errors.
    """
    _defreason = "Unknown DAL Query Error"

    def __init__(self, reason=None, label=None, url=None):
        """
        Parameters
        ----------
        reason : str
           a message describing the cause of the error.  This should
           be set to the content of the INFO error element.
        label : str
           the identifying name of the error.  This should be the
           value of the INFO element's value attribute within the
           VOTable response that describes the error.
        url : str
           the query URL that produced the error
        """
        super().__init__(reason, url)
        self._label = label

    @property
    def label(self):
        """
        the identifing name for the error given in the DAL query response.
        DAL queries that produce an error which is detectable on the server
        will respond with a VOTable containing an INFO element that contains
        the description of the error.  This property contains the value of
        the INFO's value attribute.
        """
        return self._label


class PyvoUserWarning(AstropyUserWarning):
    pass


class DALOverflowWarning(UserWarning):
    pass
