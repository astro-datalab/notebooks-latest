# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test for mivot.viewer.mivot_viewer.py
"""
import os
import pytest
import re
from astropy.utils.data import get_pkg_data_filename
from pyvo.mivot.utils.vocabulary import Constant
from pyvo.mivot.utils.dict_utils import DictUtils
from pyvo.mivot.utils.exceptions import MappingError
from pyvo.mivot.version_checker import check_astropy_version
from pyvo.mivot import MivotViewer
from astropy import version as astropy_version


@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_get_first_instance_dmtype(path_to_first_instance):
    """
    Test the function get_first_instance_dmtype() which is
    used to find the first INSTANCE/COLLECTION in TEMPLATES.
    """
    m_viewer = MivotViewer(votable_path=path_to_first_instance)
    assert m_viewer.get_first_instance_dmtype("one_instance") == "one_instance"
    assert m_viewer.get_first_instance_dmtype("coll_and_instances") == "first"
    assert m_viewer.get_first_instance_dmtype("one_collection") == Constant.ROOT_COLLECTION
    assert m_viewer.get_first_instance_dmtype("only_collection") == Constant.ROOT_COLLECTION
    with pytest.raises(Exception, match="Can't find the first INSTANCE/COLLECTION in TEMPLATES"):
        m_viewer.get_first_instance_dmtype("empty")


@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_table_ref(m_viewer):
    """
    Test if the mivot_viewer can find each table_ref and connect to the right table_ref.
    Test if the mivot_viewer can find each models.
    """
    assert m_viewer._mapped_tables == ['_PKTable', 'Results']
    with pytest.raises(Exception,
                       match=re.escape(r"The table first_table doesn't match with any mapped_table "
                                       r"(['_PKTable', 'Results']) encountered in TEMPLATES")):
        m_viewer._connect_table("wrong_tableref")
    assert m_viewer.connected_table_ref == Constant.FIRST_TABLE
    assert (m_viewer.get_models()
            == {'mango': 'file:/Users/sao/Documents/IVOA/GitHub/ivoa-dm-examples/tmp/Mango-v1.0.vo-dml.xml',
                'cube': 'https://volute.g-vo.org/svn/trunk/projects/dm/Cube/vo-dml/Cube-1.0.vo-dml.xml',
                'ds': 'https://volute.g-vo.org/svn/trunk/projects/dm/'
                      'DatasetMetadata/vo-dml/DatasetMetadata-1.0.vo-dml.xml',
                'meas': 'https://www.ivoa.net/xml/Meas/20200908/Meas-v1.0.vo-dml.xml',
                'coords': 'https://www.ivoa.net/xml/STC/20200908/Coords-v1.0.vo-dml.xml',
                'ivoa': 'https://www.ivoa.net/xml/VODML/IVOA-v1.vo-dml.xml'})


@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_global_getters(m_viewer):
    """
    Test each getter for TEMPLATES of the model_viewer.
    """
    assert m_viewer.get_table_ids() == ['_PKTable', 'Results']
    assert m_viewer.get_globals_models() == DictUtils.read_dict_from_file(
        get_pkg_data_filename("data/reference/globals_models.json"))
    assert m_viewer.get_templates_models() == DictUtils.read_dict_from_file(
        get_pkg_data_filename("data/reference/templates_models.json"))
    m_viewer._connect_table('_PKTable')
    row = m_viewer.next_table_row()
    assert row[0] == '5813181197970338560'
    assert row[1] == 'G'
    row = m_viewer.next_table_row()
    assert row[0] == '5813181197970338560'
    assert row[1] == 'BP'
    m_viewer.rewind()
    row = m_viewer.next_table_row()
    assert row[0] == '5813181197970338560'
    assert row[1] == 'G'


@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_no_mivot(path_no_mivot):
    """
    Test each getter for GLOBALS of the model_viewer specific .
    """
    m_viewer = MivotViewer(path_no_mivot)
    assert m_viewer.get_table_ids() is None
    assert m_viewer.get_globals_models() is None

    assert m_viewer.get_templates_models() is None
    with pytest.raises(MappingError):
        m_viewer._connect_table('_PKTable')
    with pytest.raises(MappingError):
        m_viewer._connect_table()

    assert m_viewer.next_table_row() is None


def test_check_version(path_to_viewer):
    if not check_astropy_version():
        with pytest.raises(Exception,
                           match=f"Astropy version {astropy_version.version} "
                                 f"is below the required version 6.0 for the use of MIVOT."):
            MivotViewer(votable_path=path_to_viewer)
    if astropy_version.version is None:
        assert not check_astropy_version()
    elif astropy_version.version < '6.0':
        assert not check_astropy_version()
    else:
        assert check_astropy_version() is True


@pytest.fixture
def m_viewer():
    if not check_astropy_version():
        pytest.skip("MIVOT test skipped because of the astropy version.")

    votable_name = "test.mivot_viewer.xml"
    votable_path = get_pkg_data_filename(os.path.join("data", votable_name))
    return MivotViewer(votable_path=votable_path)


@pytest.fixture
def path_to_viewer():
    if not check_astropy_version():
        pytest.skip("MIVOT test skipped because of the astropy version.")

    votable_name = "test.mivot_viewer.xml"
    return get_pkg_data_filename(os.path.join("data", votable_name))


@pytest.fixture
def path_to_first_instance():
    votable_name = "test.mivot_viewer.first_instance.xml"
    return get_pkg_data_filename(os.path.join("data", votable_name))


@pytest.fixture
def path_no_mivot():
    votable_name = "test.mivot_viewer.no_mivot.xml"
    return get_pkg_data_filename(os.path.join("data", votable_name))
