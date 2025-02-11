# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test for mivot.viewer.model_viewer_level2.py
"""
import pytest
try:
    from defusedxml.ElementTree import Element as element
except ImportError:
    from xml.etree.ElementTree import Element as element
from astropy.utils.data import get_pkg_data_filename
from pyvo.mivot.version_checker import check_astropy_version
from pyvo.mivot import MivotViewer
from pyvo.mivot.utils.exceptions import MivotError


@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_xml_viewer(m_viewer):

    m_viewer.next()
    xml_viewer = m_viewer.xml_viewer
    with pytest.raises(MivotError,
                       match="Cannot find dmrole wrong_role in any instances of the VOTable"):
        xml_viewer.get_instance_by_role("wrong_role")

    with pytest.raises(MivotError,
                       match="Cannot find dmrole wrong_role in any instances of the VOTable"):
        xml_viewer.get_instance_by_role("wrong_role", all_instances=True)

    with pytest.raises(MivotError,
                       match="Cannot find dmtype wrong_dmtype in any instances of the VOTable"):
        xml_viewer.get_instance_by_type("wrong_dmtype")

    with pytest.raises(MivotError,
                       match="Cannot find dmtype wrong_dmtype in any instances of the VOTable"):
        xml_viewer.get_instance_by_type("wrong_dmtype", all_instances=True)

    with pytest.raises(MivotError,
                       match="Cannot find dmrole wrong_role in any collections of the VOTable"):
        xml_viewer.get_collection_by_role("wrong_role")

    with pytest.raises(MivotError,
                       match="Cannot find dmrole wrong_role in any collections of the VOTable"):
        xml_viewer.get_collection_by_role("wrong_role", all_instances=True)

    instances_list_role = xml_viewer.get_instance_by_role("cube:MeasurementAxis.measure")
    assert isinstance(instances_list_role, element)

    instances_list_role = xml_viewer.get_instance_by_role("cube:MeasurementAxis.measure", all_instances=True)
    assert len(instances_list_role) == 3

    instances_list_type = xml_viewer.get_instance_by_type("cube:Observable")
    assert isinstance(instances_list_type, element)

    instances_list_type = xml_viewer.get_instance_by_type("cube:Observable", all_instances=True)
    assert len(instances_list_type) == 3

    collections_list_role = xml_viewer.get_collection_by_role("cube:NDPoint.observable")
    assert isinstance(collections_list_role, element)

    collections_list_role = xml_viewer.get_collection_by_role("cube:NDPoint.observable", all_instances=True)
    assert len(collections_list_role) == 1


@pytest.fixture
def m_viewer():
    return MivotViewer(get_pkg_data_filename("data/test.mivot_viewer.xml"),
                       tableref="Results")
