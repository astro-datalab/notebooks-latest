# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.utils.vocabularies

It's hard to write meaningful tests for those that don't require network
connectivity because essentially it's all just wrapping downloads.  Hence,
I'm just giving in rather than bother with a mock server.
"""

import os
import pathlib
import time

import pytest

from astropy.utils import data

from pyvo.dal.exceptions import PyvoUserWarning
from pyvo.utils import vocabularies


@pytest.mark.remote_data
class TestVocabularies:

    def test_basic_getting(self):
        # clear the lru cache in case someone else has already used
        # datalink/core.
        vocabularies.get_vocabulary.cache_clear()
        voc = vocabularies.get_vocabulary("datalink/core")
        assert "progenitor" in voc["terms"]
        assert data.is_url_in_cache("http://www.ivoa.net/rdf/datalink/core")

    def test_label_getting(self):
        voc = vocabularies.get_vocabulary("datalink/core")
        assert (vocabularies.get_label(voc, "coderived")
                == "Coderived Data")

    def test_label_getting_default(self):
        voc = vocabularies.get_vocabulary("datalink/core")
        assert vocabularies.get_label(voc, "oov", "Missing") == "Missing"

    def test_refreshing(self):
        voc = vocabularies.get_vocabulary("datalink/core", force_update=True)

        # first make sure that things didn't break
        assert "progenitor" in voc["terms"]

        # now guess that a download has actually happened; we don't want
        # to reflect cache name generation here, so we just check if there's
        # a recent download in the cache directory
        dldir = data._get_download_cache_loc()
        with os.scandir(dldir) as entries:
            last_change = 0
            for entry in entries:
                last_change = max(last_change, entry.stat().st_mtime)
        assert time.time() - last_change < 2

    def test_non_existing_voc(self):
        with pytest.raises(vocabularies.VocabularyError):
            vocabularies.get_vocabulary("not_an_ivoa_vocabulary")

    def test_failed_update(self):
        # Create a fake vocabulary and make it so old the machine
        # will want to refresh it.
        fake_voc = "http://www.ivoa.net/rdf/astropy-test-failure"

        cache_dir = pathlib.Path(data._get_download_cache_loc()
                                 ) / data._url_to_dirname(fake_voc)
        cache_dir.mkdir(exist_ok=True)

        cache_name = cache_dir / "contents"
        with open(cache_name, "w") as f:
            f.write("{}")
        with open(cache_dir / "url", "w") as f:
            f.write(fake_voc)
        os.utime(cache_name, (1000000000, 1000000000))

        with pytest.warns(PyvoUserWarning) as msgs:
            vocabularies.get_vocabulary("astropy-test-failure")
        # this sometimes catches a warning about an unclosed socket that,
        # I think, originates somewhere else; let me work around it for
        # the moment.
        for msg in msgs:
            if str(msg.message) == ("Updating cache for the vocabulary"
                                    " astropy-test-failure failed: HTTP Error 404: Not Found"):
                break
        else:
            raise AssertionError("No warning about failed cache update")
