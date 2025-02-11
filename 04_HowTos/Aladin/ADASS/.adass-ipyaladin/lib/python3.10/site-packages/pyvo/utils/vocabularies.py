"""
A shallow interface to IVOA vocabularies.

See http://ivoa.net/documents/Vocabularies/ (>= version 2) for the
larger background.  In this module, we essentially wrap the retrieval
and caching of the desise files.
"""

import functools
import json
import os
import time
import warnings

from astropy.utils.data import download_file, clear_download_cache

from pyvo.dal.exceptions import PyvoUserWarning


IVOA_VOCABULARY_ROOT = "http://www.ivoa.net/rdf/"


class VocabularyError(Exception):
    """A generic error that occurred when interacting with the IVOA
    vocabulary repository.
    """


@functools.lru_cache()
def get_vocabulary(voc_name, force_update=False):
    """returns an IVOA vocabulary in its "desise" form.

    See Vocabularies in the VO 2 to see what is inside of this.

    This will use a cache to avoid repeated updates, but it
    will attempt to re-download if the cached copy is older than 6 months.
    """
    src_url = IVOA_VOCABULARY_ROOT + voc_name
    if force_update:
        clear_download_cache(src_url)

    try:
        src_name = download_file(
            src_url,
            cache=True,
            show_progress=False,
            http_headers={"accept": "application/x-desise+json"})
    except Exception as msg:
        raise VocabularyError("No such vocabulary: {} ({})".format(
            voc_name, msg))

    if time.time() - os.path.getmtime(src_name) > 3600 * 60 * 150:
        # attempt a re-retrieval, but ignore failure
        try:
            src_name = download_file(
                IVOA_VOCABULARY_ROOT + voc_name,
                cache="update", show_progress=False,
                http_headers={"accept": "application/x-desise+json"})
        except Exception as msg:
            warnings.warn("Updating cache for the vocabulary"
                          f" {voc_name} failed: {msg}",
                          category=PyvoUserWarning)

    with open(src_name, "r", encoding="utf-8") as f:
        return json.load(f)


def get_label(voc, term, default=None):
    """returns the label of term if it's in the desise vocabulary voc,
    term capitalised otherwise.
    """
    if term in voc["terms"]:
        return voc["terms"][term]["label"]
    else:
        return default


# vi:et:sw=4:sta
