'''
The first service in operation the annotates query responses in the fly is Vizier
https://cds/viz-bin/mivotconesearch/VizierParams
Data are mapped o the EPochPropagtion model as it is implemented in the current code.
This test case is based on 2 VOTables:

- The Vizier native (vizier_cs_withname.xml) where all ATTRIBUTE@ref are
  based on FIELD@name even when a field has an ID.
- The patched vizier (vizier_cs_withid.xml) where all ATTRIBUTE@ref are
  based on FIELD@ID or FIELD@name if it exists.
The test checks that:
- The position fields can be retrieved through the mapping.
- Both cases give the same results
A third test checks the case where a referenced column does not exist.

Created on 26 janv. 2024
@author: michel
'''
import os
import pytest
from urllib.request import urlretrieve
from pyvo.mivot.version_checker import check_astropy_version
from pyvo.mivot import MivotViewer
from pyvo.mivot.utils.exceptions import MivotError


@pytest.fixture
def data_path():
    return os.path.dirname(os.path.realpath(__file__))


@pytest.fixture
def data_sample_url():
    return "https://raw.githubusercontent.com/ivoa/dm-usecases/main/pyvo-ci-sample/"


@pytest.fixture
def delt_coo():
    """ acceptable delta for coordinate value comparisons
    """
    return 0.0000005


@pytest.fixture
def path_to_withname(data_path, data_sample_url):

    votable_name = "vizier_cs_withname.xml"
    votable_path = os.path.join(data_path, "data", votable_name)
    urlretrieve(data_sample_url + votable_name,
                votable_path)

    yield votable_path
    os.remove(votable_path)


@pytest.fixture
def path_to_withid(data_path, data_sample_url):

    votable_name = "vizier_cs_withid.xml"
    votable_path = os.path.join(data_path, "data", votable_name)
    urlretrieve(data_sample_url + votable_name,
                votable_path)

    yield votable_path
    os.remove(votable_path)


@pytest.fixture
def path_to_badref(data_path, data_sample_url):

    votable_name = "vizier_cs_badref.xml"
    votable_path = os.path.join(data_path, "data", votable_name)
    urlretrieve(data_sample_url + votable_name,
                votable_path)

    yield votable_path
    os.remove(votable_path)


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_with_name(path_to_withname, delt_coo):
    """ Test that the epoch propagation works with all FIELDs referenced by name or by ID
    """
    m_viewer = MivotViewer(votable_path=path_to_withname)
    m_viewer.next()
    mivot_object = m_viewer.dm_instance

    assert abs(mivot_object.longitude.value - 52.2340018) < delt_coo
    assert abs(mivot_object.latitude.value - 59.8937333) < delt_coo
    assert abs(mivot_object.pmLongitude.value - 1.5) < delt_coo
    assert abs(mivot_object.pmLatitude.value - -12.30000019) < delt_coo
    assert str(mivot_object.epoch.value) == '2013.418'
    assert str(mivot_object.coordSys.spaceRefFrame.value) == 'ICRS'

    m_viewer.next()

    assert abs(mivot_object.longitude.value - 32.2340018) < delt_coo
    assert abs(mivot_object.latitude.value - 49.8937333) < delt_coo
    assert abs(mivot_object.pmLongitude.value - 1.5) < delt_coo
    assert abs(mivot_object.pmLatitude.value - -12.30000019) < delt_coo
    assert str(mivot_object.epoch.value) == '2013.418'
    assert str(mivot_object.coordSys.spaceRefFrame.value) == 'ICRS'


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_with_id(path_to_withid, delt_coo):
    """ Test that the epoch propagation works with all FIELDs referenced by name or by ID
    """
    # Test with all FILELDs referenced by names
    m_viewer = MivotViewer(votable_path=path_to_withid)
    m_viewer.next()
    mivot_instance = m_viewer.dm_instance
    assert abs(mivot_instance.longitude.value - 52.2340018) < delt_coo
    assert abs(mivot_instance.latitude.value - 59.8937333) < delt_coo
    assert abs(mivot_instance.pmLongitude.value - 1.5) < delt_coo
    assert abs(mivot_instance.pmLatitude.value - -12.30000019) < delt_coo
    assert str(mivot_instance.epoch.value) == '2013.418'


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_bad_ref(path_to_badref, delt_coo):
    """ Test that the epoch propagation works with all FIELDs referenced by name or by ID
    """
    # Test with all FILELDs referenced by names
    with (pytest.raises(MivotError, match="Attribute mango:EpochPosition.epoch can not be set.*")):
        MivotViewer(votable_path=path_to_badref)
