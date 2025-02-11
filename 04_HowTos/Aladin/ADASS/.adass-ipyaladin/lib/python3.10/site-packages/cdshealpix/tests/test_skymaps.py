from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pytest

from ..skymap import Skymap
from .. import cdshealpix

path_to_test_skymap = Path(__file__).parent.resolve() / "resources" / "skymap.fits"


def test_read():
    values = Skymap.from_fits(path_to_test_skymap).values
    assert values.dtype == np.int32
    assert len(values) == 49152


def test_read_write_read_conservation():
    skymap = Skymap.from_fits(path_to_test_skymap)
    with NamedTemporaryFile() as fp:
        skymap.to_fits(fp.name)
        skymap2 = Skymap.from_fits(fp.name)
        assert all(skymap.values == skymap2.values)


def test_plot():
    skymap = Skymap.from_array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    img = cdshealpix.pixels_skymap(skymap.values, 256, True)
    assert img.shape == (256, 512, 4)


def test_depth():
    n = 12
    skymap = Skymap.from_array(np.zeros(12 * 4**n))
    assert skymap.depth == n


def test_from_array():
    with pytest.raises(ValueError, match="The length of values should be*"):
        Skymap.from_array(np.zeros(3))
    with pytest.raises(ValueError, match="The accepted types are*"):
        Skymap.from_array(["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"])
    with pytest.raises(ValueError, match="Skymap values should be one-dimensional*"):
        Skymap.from_array([[1, 2, 3], [1, 2, 3]])
