# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains a contains the high-level functions to read the various
VOSI Endpoints.
"""

from astropy.utils.xml.writer import XMLWriter
from astropy.io.votable.util import convert_to_writable_filelike

from ...utils.xml.elements import xmlattribute, parse_for_object
from .tree import JobSummary, Jobs

__all__ = ["parse_job", "parse_job_list", "JobFile"]


def parse_job_list(
    source, pedantic=None, filename=None, _debug_python_based_parser=False
):
    """
    Parses a job xml file (or file-like object), and returns a
    `~pyvo.io.uws.tree.Jobs` object.

    Parameters
    ----------
    source : str or readable file-like object
        Path or file object containing a tableset xml file.
    pedantic : bool, optional
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        the standard Python mechanisms.  See the `warnings`
        module in the Python standard library for more information.
        Defaults to False.
    filename : str, optional
        A filename, URL or other identifier to use in error messages.
        If *filename* is None and *source* is a string (i.e. a path),
        then *source* will be used as a filename for error messages.
        Therefore, *filename* is only required when source is a
        file-like object.

    Returns
    -------
    `~pyvo.io.uws.tree.Jobs` object

    See also
    --------
    pyvo.io.vosi.exceptions : The exceptions this function may raise.
    """
    return parse_for_object(source, Jobs, pedantic, filename,
                            _debug_python_based_parser).joblist


def parse_job(
    source, pedantic=None, filename=None, _debug_python_based_parser=False
):
    """
    Parses a job xml file (or file-like object), and returns a
    `~pyvo.io.uws.endpoint.JobFile` object.

    Parameters
    ----------
    source : str or readable file-like object
        Path or file object containing a tableset xml file.
    pedantic : bool, optional
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        the standard Python mechanisms.  See the `warnings`
        module in the Python standard library for more information.
        Defaults to False.
    filename : str, optional
        A filename, URL or other identifier to use in error messages.
        If *filename* is None and *source* is a string (i.e. a path),
        then *source* will be used as a filename for error messages.
        Therefore, *filename* is only required when source is a
        file-like object.

    Returns
    -------
    `~pyvo.io.uws.endpoint.JobFile` object

    See also
    --------
    pyvo.io.vosi.exceptions : The exceptions this function may raise.
    """
    return parse_for_object(source, JobFile, pedantic, filename,
                            _debug_python_based_parser)


class JobFile(JobSummary):
    """
    availability element: represents an entire file.
    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    def __init__(self, config=None, pos=None, **kwargs):
        super().__init__(config=config, pos=pos, **kwargs)
        self._version = None

    @xmlattribute
    def version(self):
        return self._version

    @version.setter
    def version(self, version):
        self._version = version

    def parse(self, iterator, config):
        for start, tag, data, pos in iterator:
            if start and tag == 'xml':
                pass
            elif start and tag == 'job':
                # version was not required in v1.0, so default to that.
                self._version = data.get('version', '1.0')
                break

        return super().parse(iterator, config)

    def to_xml(self, fd):
        with convert_to_writable_filelike(fd) as _fd:
            w = XMLWriter(_fd)

            xml_header = (
                '<?xml version="1.0" encoding="utf-8"?>\n'
                '<!-- Produced with pyvo.io.uws -->\n'
            )

            w.write(xml_header)

            super().to_xml(w)
