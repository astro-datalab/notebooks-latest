import json
from pathlib import Path

import pytest
from astropy.io import fits
from astropy.time import Time

from ..abstract_moc import AbstractMOC
from ..fmoc import FrequencyMOC
from ..moc import MOC
from ..stmoc import STMOC
from ..tmoc import TimeMOC

# -----------------------------------------------------------fixtures


@pytest.fixture
def all_moc_types():
    """Create a fixture with a MOC of each type."""
    moc = MOC.from_str("0/0-11")
    tmoc = TimeMOC.from_str("0/0-1")
    fmoc = FrequencyMOC.from_str("0/0-1")
    stmoc = STMOC.from_spatial_coverages(Time("2023-11-13"), Time("2023-11-14"), moc)
    return [moc, tmoc, fmoc, stmoc]


@pytest.fixture
def moc():
    return MOC.from_str("0/0-11")


@pytest.fixture
def path(tmp_path):
    return tmp_path / "path"


# ------------------------------------------------------instantiation


def test_failing_instantiation():
    with pytest.raises(
        TypeError,
        match="Can't instantiate abstract class AbstractMOC*",
    ):
        AbstractMOC()


# ---------------------------------------------------------------save


def test_failing_save(moc, path):
    """All the pitfalls in the save method."""
    moc.save(path)
    assert path.is_file()
    with pytest.raises(OSError, match=r"File *"):
        moc.save(path)
    with pytest.raises(ValueError, match=r"The ``fold`` argument*"):
        moc.save(path, format="fits", overwrite=True, fold=80)
    with pytest.raises(ValueError, match=r"The ``fits_keyword`` argument*"):
        moc.save(path, format="json", fits_keywords={"TEST", "keyword"}, overwrite=True)
    with pytest.raises(ValueError, match=r"'format' should be *"):
        moc.save(path, overwrite=True, format="test")


def test_passing_save(all_moc_types, path):
    for moc in all_moc_types:
        # ----
        # fits
        # ----
        moc.save(path, format="fits", fits_keywords={"TEST": "written"})
        assert moc == moc.load(path, format="fits")
        # if it's a spatial moc it's the all sky
        if isinstance(moc, MOC):
            assert moc.sky_fraction == 1
        # test that astroquery can also deserialize
        with fits.open(path) as fits_moc:
            # and that the additional keyword is there
            assert "TEST    = 'written '" in str(fits_moc[1].header)
        # delete the moc we just created
        path.unlink()
        # ----
        # json
        # ----
        # without fold
        moc.save(path, format="json")
        assert moc == moc.load(path, format="json")
        # test standard module can also read this
        with Path.open(path) as f:
            moc_json = json.load(f)
        # test that it is a valid dictionary
        indices = moc_json[0]["s"]["0"] if isinstance(moc, STMOC) else moc_json["0"]
        assert {0, 1}.issubset(set(indices))
        # delete the moc we just created
        path.unlink()
        # -----
        # ascii
        # -----
        fold = 20  # will only be triggered by the STMOC
        moc.save(path, format="ascii", fold=fold)
        assert moc == moc.load(path, format="ascii")
        with Path.open(path) as f:
            moc_string = f.read()
        moc_string_no_whitespaces = moc_string.replace(" ", "")
        assert len(max(moc_string_no_whitespaces.split("\n"), key=len)) <= fold
        assert moc_string.replace("\n", "") == moc.to_string(
            "ascii",
            fold=fold,
        ).replace("\n", "")
        # delete the moc we just created
        path.unlink()


# --------------------------------------------------------------write


def test_write(moc, path, mocker):
    # test that the deprecation warning is raised
    with pytest.warns(
        DeprecationWarning,
        match=r"This method is deprecated. Use MOC.save*",
    ):
        mocked_save = mocker.patch("mocpy.abstract_moc.AbstractMOC.save")
        moc.write(path)
        # test that save was called by write
        mocked_save.assert_called_once_with(
            path,
            "fits",
            False,  # noqa: FBT003
            False,  # noqa: FBT003
            fits_keywords=None,
        )


# --------------------------------------------------------operations


def test_operations(all_moc_types):
    # nothing too fancy here, we just check that the operations can be called
    for moc in all_moc_types:
        # intersection
        assert moc.intersection(moc) == moc
        assert moc.intersection(moc, moc) == moc
        # union
        assert moc.union(moc) == moc
        assert moc.union(moc, moc) == moc
        # symmetric difference
        if isinstance(moc, STMOC):
            with pytest.raises(
                NotImplementedError,
                match=(
                    "Symmetric difference is not implemented yet for Space-Time MOCs"
                ),
            ):
                moc.symmetric_difference(moc)
        else:
            assert moc.symmetric_difference(moc).empty()
            assert moc.symmetric_difference(moc, moc) == moc
        # difference
        assert moc.difference(moc).empty()
        assert moc.difference(moc, moc).empty()
