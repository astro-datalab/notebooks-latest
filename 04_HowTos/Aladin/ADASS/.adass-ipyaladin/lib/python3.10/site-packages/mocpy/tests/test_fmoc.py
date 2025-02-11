import numpy as np
import pytest
from astropy import units as u

from ..fmoc import FrequencyMOC


def test_new_empty():
    fmoc = FrequencyMOC.new_empty(FrequencyMOC.MAX_ORDER)
    assert fmoc.empty()


def test_max_order():
    fmoc = FrequencyMOC.new_empty(32)
    assert fmoc.max_order == 32


def test_n_cells():
    assert FrequencyMOC.n_cells(0) == 2
    with pytest.raises(
        ValueError,
        match=f"The depth should be comprised between 0 and {FrequencyMOC.MAX_ORDER}*",
    ):
        FrequencyMOC.n_cells(-1)
        FrequencyMOC.n_cells(FrequencyMOC.MAX_ORDER + 1)
    assert FrequencyMOC.n_cells(5) == 2 * FrequencyMOC.n_cells(4)


def test_to_depth59_ranges():
    fmoc = FrequencyMOC.new_empty(FrequencyMOC.MAX_ORDER).complement()
    # 2^60 = 1152921504606846976
    assert (
        fmoc.to_depth59_ranges
        == np.array([[0, FrequencyMOC.MAX_INDEX_EXCLUSIVE]], dtype=np.uint64)
    ).all()


def test_to_hz_ranges():
    fmoc = FrequencyMOC.new_empty(FrequencyMOC.MAX_ORDER).complement()
    assert np.allclose(
        fmoc.to_hz_ranges(),
        np.array(
            [[FrequencyMOC.FREQ_MIN_HZ, FrequencyMOC.FREQ_MAX_HZ]],
            dtype=np.float64,
        ),
    )


def test_degrade_to_order():
    fmoc = FrequencyMOC.new_empty(FrequencyMOC.MAX_ORDER).complement()
    fmoc = fmoc.degrade_to_order(0)
    # No need for more than this since the same degrade code is used for smoc and tmoc
    assert fmoc.max_order == 0
    assert (
        fmoc.to_depth59_ranges
        == np.array([[0, FrequencyMOC.MAX_INDEX_EXCLUSIVE]], dtype=np.uint64)
    ).all()


def test_from_depth59_ranges():
    ranges = np.array([[0, FrequencyMOC.MAX_INDEX_EXCLUSIVE]], dtype=np.uint64)
    fmoc = FrequencyMOC.from_depth59_ranges(FrequencyMOC.MAX_ORDER, ranges)
    assert fmoc == FrequencyMOC.new_empty(FrequencyMOC.MAX_ORDER).complement()


def test_from_frequencies():
    freq = np.array([FrequencyMOC.FREQ_MIN_HZ, FrequencyMOC.FREQ_MAX_HZ]) * u.Hz
    fmoc = FrequencyMOC.from_frequencies(0, freq)
    assert fmoc == FrequencyMOC.new_empty(0).complement()


def test_from_frequency_ranges():
    fmin = np.array([FrequencyMOC.FREQ_MIN_HZ]) * u.Hz
    fmax = np.array([FrequencyMOC.FREQ_MAX_HZ]) * u.Hz
    fmoc = FrequencyMOC.from_frequency_ranges(58, fmin, fmax)
    assert fmoc == FrequencyMOC.new_empty(58).complement()


def test_min_freq():
    fmoc = FrequencyMOC.new_empty(FrequencyMOC.MAX_ORDER).complement()
    assert np.allclose(fmoc.min_freq, FrequencyMOC.FREQ_MIN_HZ * u.Hz)


def test_max_freq():
    fmoc = FrequencyMOC.new_empty(FrequencyMOC.MAX_ORDER).complement()
    assert np.allclose(fmoc.max_freq, (FrequencyMOC.FREQ_MAX_HZ) * u.Hz)


def test_contains():
    fmoc = FrequencyMOC.new_empty(12).complement()
    freq = (
        np.array([0.5 * (FrequencyMOC.FREQ_MIN_HZ + FrequencyMOC.FREQ_MAX_HZ)]) * u.Hz
    )
    assert np.all(fmoc.contains(freq))


def test_order_to_relative_precision():
    for order in range(FrequencyMOC.MAX_ORDER):
        v = FrequencyMOC.order_to_relative_precision(order)
        p = 0.5 * (v[0] + v[1])
        o = FrequencyMOC.relative_precision_to_order(p)
        assert order == o
