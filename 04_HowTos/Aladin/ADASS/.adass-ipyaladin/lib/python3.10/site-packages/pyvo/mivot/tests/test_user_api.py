"""
Check the API exposed to the user works properly
This unit test can also be used as a code example
Created on 26 Feb. 2024
@author: michel
"""

import os
import pytest
from urllib.request import urlretrieve
import astropy.units as u
from astropy.coordinates import SkyCoord
from pyvo.dal.scs import SCSService
from pyvo.mivot.version_checker import check_astropy_version
from pyvo.mivot import MivotViewer
from astropy.io.votable import parse
from pyvo.mivot.utils.dict_utils import DictUtils

ref_ra = [
    0.04827189,
    0.16283175,
    0.29222255,
    0.42674592,
    359.5190115,
    359.94372764,
]
ref_dec = [
    -0.36042119,
    0.22293899,
    -0.07592034,
    -0.21749947,
    -0.1281483,
    -0.28005255,
]
ref_pmdec = [
    -11.67,
    -3.09,
    -73.28,
    -114.08,
    -19.05,
    -25.43
]
ref_pmra = [
    61.75,
    39.02,
    54.94,
    20.73,
    -45.19,
    -5.14
]


@pytest.fixture
def data_path():
    return os.path.dirname(os.path.realpath(__file__))


@pytest.fixture
def data_sample_url():
    return "https://raw.githubusercontent.com/ivoa/dm-usecases/main/pyvo-ci-sample/"


@pytest.fixture
def vizier_url():
    return "https://cdsarc.cds.unistra.fr/beta/viz-bin/mivotconesearch/I/239/hip_main"


@pytest.fixture
def delt_coo():
    """acceptable delta for coordinate value comparisons"""
    return 0.0000001


@pytest.fixture
def path_to_votable(data_path, data_sample_url):

    votable_name = "vizier_for_user_api.xml"
    votable_path = os.path.join(data_path, "data", votable_name)
    urlretrieve(data_sample_url + votable_name, votable_path)

    yield votable_path
    os.remove(votable_path)


@pytest.fixture
def path_to_full_mapped_votable(data_path, data_sample_url):

    votable_name = "gaia_epoch_propagation_full.xml"
    votable_path = os.path.join(data_path, "data", votable_name)
    urlretrieve(data_sample_url + votable_name, votable_path)

    yield votable_path
    os.remove(votable_path)


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_mivot_viewer_next(path_to_votable):
    """
    Check that ten MivotViewer iterating over the data rows provides the expected values
    """
    mivot_viewer = MivotViewer(path_to_votable)
    mivot_instance = mivot_viewer.dm_instance
    assert mivot_instance.dmtype == "mango:EpochPosition"
    assert mivot_instance.coordSys.spaceRefFrame.value == "ICRS"
    ra = []
    dec = []
    pmra = []
    pmdec = []
    while mivot_viewer.next():
        ra.append(mivot_instance.longitude.value)
        dec.append(mivot_instance.latitude.value)
        pmra.append(mivot_instance.pmLongitude.value)
        pmdec.append(mivot_instance.pmLatitude.value)
    assert ra == ref_ra
    assert dec == ref_dec
    assert pmra == ref_pmra
    assert pmdec == ref_pmdec


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_mivot_tablerow_next(path_to_votable):
    """
    Check that the MIVOT interpreter can be applied to a classical table readout.
    - The MivotViewer is initialized on the first data row (behind the stage)
      to be able to build a MivoInstance that will provide a model view on the data
    - The data table is read in a classical way
    - The MivotInstance is updated with data row providing so a model view on it
    """
    votable = parse(path_to_votable)
    table = votable.resources[0].tables[0]
    mivot_viewer = MivotViewer(votable)

    mivot_instance = mivot_viewer.dm_instance
    assert mivot_instance.dmtype == "mango:EpochPosition"
    assert mivot_instance.coordSys.spaceRefFrame.value == "ICRS"
    ra = []
    dec = []
    pmra = []
    pmdec = []
    for rec in table.array:
        mivot_instance.update(rec)
        ra.append(mivot_instance.longitude.value)
        dec.append(mivot_instance.latitude.value)
        pmra.append(mivot_instance.pmLongitude.value)
        pmdec.append(mivot_instance.pmLatitude.value)
    assert ra == ref_ra
    assert dec == ref_dec
    assert pmra == ref_pmra
    assert pmdec == ref_pmdec


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_external_iterator(path_to_votable, delt_coo):
    """Checks that the values returned by MIVOT are
    the same as those read by Astropy
    - The MivotViewer is initialized on the first data row (behind the stage)
      to be able to build a MivoInstance that will provide a model view on the data
    - The data table is read in a classical way
    - The MivotInstance is updated with data row providing so a model view on it
    - The attribute values are then checked against the table data
    """
    # parse the VOTable outside of the viewer
    votable = parse(path_to_votable)
    table = votable.resources[0].tables[0]
    # init the viewer
    mivot_viewer = MivotViewer(votable)
    mivot_instance = mivot_viewer.dm_instance
    for rec in table.array:
        mivot_instance.update(rec)
        assert rec["RAICRS"] == mivot_instance.longitude.value
        assert rec["DEICRS"] == mivot_instance.latitude.value
        assert rec["pmRA"] == mivot_instance.pmLongitude.value
        assert rec["pmDE"] == mivot_instance.pmLatitude.value


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_with_withstatement(path_to_votable):
    """check that the values read by a MivotViewer used with a with statement
    are those expected
    """

    read_ra = []
    read_dec = []
    read_pmra = []
    read_pmdec = []
    with MivotViewer(path_to_votable) as mivot_viewer:
        mivot_object = mivot_viewer.dm_instance
        while mivot_viewer.next():
            read_ra.append(mivot_object.longitude.value)
            read_dec.append(mivot_object.latitude.value)
            read_pmra.append(mivot_object.pmLongitude.value)
            read_pmdec.append(mivot_object.pmLatitude.value)
    assert read_ra == ref_ra
    assert read_dec == ref_dec
    assert read_pmra == ref_pmra
    assert read_pmdec == ref_pmdec


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_with_dict(path_to_votable):
    """check that the MIVOT object user dictionary match the read data"""

    with MivotViewer(path_to_votable) as mivot_viewer:
        mivot_object = mivot_viewer.dm_instance
        # let's focus on the last data row
        while mivot_viewer.next():
            pass
    DictUtils.print_pretty_json(mivot_object.to_hk_dict())

    # check the slim (user friendly) dictionary
    assert mivot_object.to_dict() == {
        "dmtype": "mango:EpochPosition",
        "longitude": {"value": 359.94372764, "unit": "deg"},
        "latitude": {"value": -0.28005255, "unit": "deg"},
        "pmLongitude": {"value": -5.14, "unit": "mas/yr"},
        "pmLatitude": {"value": -25.43, "unit": "mas/yr"},
        "epoch": {"value": 1991.25, "unit": "year"},
        "coordSys": {
            "dmtype": "coords:SpaceSys",
            "dmid": "SpaceFrame_ICRS",
            "dmrole": "coords:Coordinate.coordSys",
            "spaceRefFrame": {"value": "ICRS"},
        },
    }
    # check the whole dictionary
    assert mivot_object.to_hk_dict() == {
        "dmtype": "mango:EpochPosition",
        "longitude": {
            "dmtype": "ivoa:RealQuantity",
            "value": 359.94372764,
            "unit": "deg",
            "ref": "RAICRS",
        },
        "latitude": {
            "dmtype": "ivoa:RealQuantity",
            "value": -0.28005255,
            "unit": "deg",
            "ref": "DEICRS",
        },
        "pmLongitude": {
            "dmtype": "ivoa:RealQuantity",
            "value": -5.14,
            "unit": "mas/yr",
            "ref": "pmRA",
        },
        "pmLatitude": {
            "dmtype": "ivoa:RealQuantity",
            "value": -25.43,
            "unit": "mas/yr",
            "ref": "pmDE",
        },
        "epoch": {
            "dmtype": "ivoa:RealQuantity",
            "value": 1991.25,
            "unit": "year",
            "ref": None,
        },
        "coordSys": {
            "dmtype": "coords:SpaceSys",
            "dmid": "SpaceFrame_ICRS",
            "dmrole": "coords:Coordinate.coordSys",
            "spaceRefFrame": {
                "dmtype": "coords:SpaceFrame",
                "value": "ICRS",
                "unit": None,
                "ref": None,
            },
        },
    }


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_with_full_dict(path_to_full_mapped_votable):
    """check that the MIVOT object user dictionary match the read data
    - The data sample used here maps the whole EpochPosition model except TimeSys
    """

    with MivotViewer(path_to_full_mapped_votable) as mivot_viewer:
        mivot_object = mivot_viewer.dm_instance
        # let's focus on the second data row
        while mivot_viewer.next():
            # check the slim (user friendly) dictionary
            assert mivot_object.to_dict() == {
                "dmtype": "mango:EpochPosition",
                "longitude": {"value": 307.79115807079, "unit": "deg"},
                "latitude": {"value": 20.43108005561, "unit": "deg"},
                "parallax": {"value": 0.4319, "unit": "mas"},
                "radialVelocity": {"value": None, "unit": "km/s"},
                "pmLongitude": {"value": -2.557, "unit": "mas/yr"},
                "pmLatitude": {"value": -5.482, "unit": "mas/yr"},
                "epoch": {"value": "2016.5"},
                "pmCosDeltApplied": {"value": True},
                "errors": {
                    "dmrole": "mango:EpochPosition.errors",
                    "dmtype": "mango:EpochPositionErrors",
                    "parallax": {
                        "dmrole": "mango:EpochPositionErrors.parallax",
                        "dmtype": "mango:ErrorTypes.PropertyError1D",
                        "sigma": {"value": 0.06909999996423721, "unit": "mas"},
                    },
                    "radialVelocity": {
                        "dmrole": "mango:EpochPositionErrors.radialVelocity",
                        "dmtype": "mango:ErrorTypes.PropertyError1D",
                        "sigma": {"value": None, "unit": "km/s"},
                    },
                    "position": {
                        "dmrole": "mango:EpochPositionErrors.position",
                        "dmtype": "mango:ErrorTypes.ErrorMatrix",
                        "sigma1": {"value": 0.0511, "unit": "mas"},
                        "sigma2": {"value": 0.0477, "unit": "mas"},
                    },
                    "properMotion": {
                        "dmrole": "mango:EpochPositionErrors.properMotion",
                        "dmtype": "mango:ErrorTypes.ErrorMatrix",
                        "sigma1": {"value": 0.06400000303983688, "unit": "mas/yr"},
                        "sigma2": {"value": 0.06700000166893005, "unit": "mas/yr"},
                    },
                },
                "correlations": {
                    "dmrole": "mango:EpochPosition.correlations",
                    "dmtype": "mango:EpochPositionCorrelations",
                    "positionPm": {
                        "dmrole": "mango:EpochPositionCorrelations.positionPm",
                        "dmtype": "mango:Correlation22",
                        "isCovariance": {"value": True},
                        "a2b1": {"value": -0.0085},
                        "a2b2": {"value": -0.2983},
                        "a1b1": {"value": -0.4109},
                        "a1b2": {"value": -0.0072},
                    },
                    "parallaxPm": {
                        "dmrole": "mango:EpochPositionCorrelations.parallaxPm",
                        "dmtype": "mango:Correlation12",
                        "isCovariance": {"value": True},
                        "a1b1": {"value": -0.2603},
                        "a1b2": {"value": -0.0251},
                    },
                    "positionParallax": {
                        "dmrole": "mango:EpochPositionCorrelations.positionParallax",
                        "dmtype": "mango:Correlation21",
                        "isCovariance": {"value": True},
                        "a2b1": {"value": 0.0069},
                        "a1b1": {"value": 0.1337},
                    },
                    "positionPosition": {
                        "dmrole": "mango:EpochPositionCorrelations.positionPosition",
                        "dmtype": "mango:Correlation22",
                        "isCovariance": {"value": True},
                        "a2b1": {"value": 0.1212},
                        "a1b2": {"value": 0.1212},
                    },
                    "properMotionPm": {
                        "dmrole": "mango:EpochPositionCorrelations.properMotionPm",
                        "dmtype": "mango:Correlation22",
                        "isCovariance": {"value": True},
                        "a2b1": {"value": 0.2688},
                        "a1b2": {"value": 0.2688},
                    },
                },
                "coordSys": {
                    "dmid": "_spacesys_icrs",
                    "dmrole": "mango:EpochPosition.coordSys",
                    "dmtype": "coords:SpaceSys",
                    "frame": {
                        "dmrole": "coords:PhysicalCoordSys.frame",
                        "dmtype": "coords:SpaceFrame",
                        "spaceRefFrame": {"value": "ICRS"},
                    },
                },
            }
            break


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_cone_search(vizier_url):
    """
    Test that data returned by Vizier CS (1 row) and read through MIVOt are those expected
    """
    scs_srv = SCSService(vizier_url)
    m_viewer = MivotViewer(
        scs_srv.search(
            pos=SkyCoord(ra=52.26708 * u.degree, dec=59.94027 * u.degree, frame="icrs"),
            radius=0.05,
        )
    )
    mivot_instance = m_viewer.dm_instance
    assert mivot_instance.dmtype == "mango:EpochPosition"
    assert mivot_instance.coordSys.spaceRefFrame.value == "ICRS"
    ra = []
    dec = []
    pmra = []
    pmdec = []
    while m_viewer.next():
        ra.append(mivot_instance.longitude.value)
        dec.append(mivot_instance.latitude.value)
        pmra.append(mivot_instance.pmLongitude.value)
        pmdec.append(mivot_instance.pmLatitude.value)
    assert ra == [52.26722684]
    assert dec == [59.94033461]
    assert pmra == [-0.82]
    assert pmdec == [-1.85]
