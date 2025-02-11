# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from numpy.testing import assert_allclose

from regions._utils.examples import make_example_dataset


class TestExampleSimulatedDataset:
    def setup_method(self):
        self.dataset = make_example_dataset(data='simulated')

    def test_source_table(self):
        source_table = self.dataset.source_table
        assert len(source_table) == 3

    def test_event_table(self):
        # source_table = self.dataset.source_table
        event_table = self.dataset.event_table
        assert len(event_table) == 300

    def test_wcs(self):
        wcs = self.dataset.wcs
        assert_allclose(wcs.wcs.crval, (0, 0))
        assert_allclose(wcs.wcs.crpix, (180, 90))
        assert_allclose(wcs.wcs.cdelt, (-1, 1))
        assert wcs.wcs.ctype[0] == 'GLON-AIT'
        assert wcs.wcs.ctype[1] == 'GLAT-AIT'

    def test_image(self):
        image = self.dataset.image
        assert image.data.shape == (180, 360)
        assert np.nansum(image.data) == 300

    def test_hdu_list(self):
        hdu_list = self.dataset.hdu_list

        # Check that all data is present
        assert hdu_list[1] == hdu_list['sources']
        assert hdu_list[2] == hdu_list['events']
        assert hdu_list[3] == hdu_list['image']
