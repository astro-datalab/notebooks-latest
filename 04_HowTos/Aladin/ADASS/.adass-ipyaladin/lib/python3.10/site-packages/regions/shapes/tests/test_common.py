# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Base class for all shape tests.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from regions._utils.optional_deps import HAS_MATPLOTLIB
from regions.core.pixcoord import PixCoord


class BaseTestRegion:

    def test_repr(self):
        assert repr(self.reg) == self.expected_repr

    def test_str(self):
        assert str(self.reg) == self.expected_str

    def test_del(self):
        # test that Region shape attributes (descriptors) cannot be
        # deleted)
        with pytest.raises(AttributeError):
            delattr(self.reg, self.reg._params[0])


class BaseTestPixelRegion(BaseTestRegion):

    def test_area(self):
        assert_allclose(self.reg.area, self.expected_area)

    def test_mask_area(self):
        try:
            mask = self.reg.to_mask(mode='exact')
            assert_allclose(np.sum(mask.data), self.expected_area)
        except NotImplementedError:
            try:
                mask = self.reg.to_mask(mode='subpixels', subpixels=30)
                assert_allclose(np.sum(mask.data), self.expected_area,
                                rtol=0.005)
            except NotImplementedError:
                pytest.skip()

    def test_contains_scalar(self):
        if len(self.inside) > 0:
            pixcoord = PixCoord(*self.inside[0])
            assert self.reg.contains(pixcoord)
            assert pixcoord in self.reg

        if len(self.outside) > 0:
            pixcoord = PixCoord(*self.outside[0])
            assert not self.reg.contains(pixcoord)
            assert pixcoord not in self.reg

    def test_contains_array_1d(self):
        pixcoord = PixCoord(*zip(*(self.inside + self.outside), strict=True))
        actual = self.reg.contains(pixcoord)
        assert_equal(actual[:len(self.inside)], True)
        assert_equal(actual[len(self.inside):], False)

        with pytest.raises(ValueError) as excinfo:
            assert pixcoord in self.reg
        assert 'coord must be scalar' in str(excinfo.value)

    def test_contains_array_2d(self):
        x, y = zip(*(self.inside + self.outside), strict=True)
        pixcoord = PixCoord([x] * 3, [y] * 3)

        actual = self.reg.contains(pixcoord)
        assert actual.shape == (3, len(x))
        assert_equal(actual[:, :len(self.inside)], True)
        assert_equal(actual[:, len(self.inside):], False)

    @pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
    def test_plot(self):
        from matplotlib.artist import Artist
        artist = self.reg.plot()
        assert isinstance(artist, Artist)


class BaseTestSkyRegion(BaseTestRegion):
    # TODO: here we should add inside/outside tests as above
    pass
