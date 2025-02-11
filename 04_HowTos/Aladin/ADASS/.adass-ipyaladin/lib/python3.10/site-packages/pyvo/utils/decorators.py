# Licensed under a 3-clause BSD style license - see LICENSE.rst

from functools import partial, wraps


def stream_decode_content(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        raw = func(*args, **kwargs)
        raw.read = partial(raw.read, decode_content=True)
        return raw

    return wrapper


def response_decode_content(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        response = func(*args, **kwargs)
        response.raw.read = partial(response.raw.read, decode_content=True)
        return response

    return wrapper
