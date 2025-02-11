#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.io.vosi
"""

import pyvo.io.uws as uws

from astropy.utils.data import get_pkg_data_filename


class TestJob:
    def test_job(self):
        job = uws.parse_job(get_pkg_data_filename(
            "data/job.xml"))

        assert job.jobid == '1337'
        assert job.version == '1.1'

        job = uws.parse_job(get_pkg_data_filename(
            "data/job-implicit-v1.0.xml"))

        assert job.version == '1.0'

    def test_error_job(self):
        job = uws.parse_job(get_pkg_data_filename(
            "data/job-error.xml"))

        assert job.jobid == '1337'
        assert job.version == '1.1'

        assert not job.errorsummary.has_detail
        assert job.errorsummary.type_ == 'fatal'
        assert job.errorsummary.message.content == 'We have problem'
