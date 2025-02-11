# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module for parsing and working with mimetypes
"""

import mimetypes
from email.message import Message

from astropy.io.fits import HDUList

from ..utils.http import use_session


mimetypes.add_type('application/fits', 'fits')
mimetypes.add_type('application/x-fits', 'fits')
mimetypes.add_type('image/fits', 'fits')
mimetypes.add_type('text/plain', 'txt')


def mime2extension(mimetype, default=None):
    """
    return a recommended file extension for a file with a given MIME-type.

    This function provides some generic mappings that can be leveraged in
    implementations of ``suggest_extension()`` in ``Record`` subclasses.

      >>> mime2extension('application/fits')
      'fits'
      >>> mime2extension('image/jpeg')
      'jpg'
      >>> mime2extension('application/x-zed', 'dat')
      'dat'

    Parameters
    ----------
    mimetype : str
       the file MIME-type byte-string to convert
    default : str
       the default extension to return if one could not be
       recommended based on ``mimetype``.  By convention,
       this should not include a preceding '.'

    Returns
    -------
    str
       the recommended extension without a preceding '.', or the
       value of ``default`` if no recommendation could be made.
    """
    if not mimetype:
        return default

    if isinstance(mimetype, bytes):
        mimetype = mimetype.decode('utf-8')

    ext = mimetypes.guess_extension(mimetype, strict=False)
    if ext is None:
        return default

    return ext.lstrip(".")


def mime_object_maker(url, mimetype, *, session=None):
    """
    return a data object suitable for the mimetype given.
    this will either return a astropy fits object or a pyvo DALResults object,
    a PIL object for conventional images or string for text content.

    Parameters
    ----------
    url : str
        the object download url
    mimetype : str
        the content mimetype
    session : object
        optional session to use for network requests

    Raises
    ------
    ValueError if the mimetype is missing or cannot be parsed correctly
    """
    if not mimetype:
        raise ValueError('mimetype required')
    session = use_session(session)
    msg = Message()
    msg['content-type'] = mimetype
    pp = msg.get_params()
    full_type = pp[0][0]
    params = pp[1:]
    mtype = [x.strip() for x in full_type.split('/')] if '/' in full_type else None
    if not mtype or len(mtype) > 2:
        raise ValueError("Can't parse mimetype \"{}\"".format(full_type))

    if mtype[0] == 'text':
        return session.get(url).text

    if mtype[1] == 'fits' or mtype[1] == 'x-fits':
        response = session.get(url)
        return HDUList.fromstring(response.content)

    if mtype[0] == 'image':
        from PIL import Image
        from io import BytesIO
        response = session.get(url)
        bio = BytesIO(response.content)
        return Image.open(bio)

    if mtype[1] == 'x-votable' or mtype[1] == 'x-votable+xml':
        # As soon as there are some kind of recursive data structures,
        # things start to get messy
        for param in params:
            if (param[0].lower() == 'content') and (param[1].lower() == 'datalink'):
                from .adhoc import DatalinkResults
                return DatalinkResults.from_result_url(url)
        from .query import DALResults
        return DALResults.from_result_url(url)
