# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test for mivot.viewer.model_viewer_level3.py and mivot.viewer.mivot_time.py
"""
import os
import pytest
from urllib.request import urlretrieve
from pyvo.mivot.version_checker import check_astropy_version
from pyvo.mivot import MivotViewer
from pyvo.mivot.utils.mivot_utils import MivotUtils


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_model_viewer3(votable_test, simple_votable):
    """
    Recursively compare an XML element with an element of MIVOT
    class with the function recursive_xml_check.
    This test run on 2 votables : votable_test and simple_votable.
    """
    m_viewer_simple_votable = MivotViewer(votable_path=simple_votable)
    MivotInstance = m_viewer_simple_votable.dm_instance
    xml_simple_votable = m_viewer_simple_votable.xml_view
    assert xml_simple_votable.tag == 'TEMPLATES'
    recusive_xml_check(xml_simple_votable, MivotInstance)
    m_viewer_votable_test = MivotViewer(votable_path=votable_test)
    m_viewer_votable_test.next()
    mivot_instance = m_viewer_votable_test.dm_instance
    xml_votable_test = m_viewer_votable_test.xml_view
    assert xml_simple_votable.tag == 'TEMPLATES'
    recusive_xml_check(xml_votable_test, mivot_instance)


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def recusive_xml_check(xml_simple_votable, MivotInstance):
    if xml_simple_votable.tag == 'TEMPLATES':
        recusive_xml_check(xml_simple_votable[0], MivotInstance)
    else:
        for child in xml_simple_votable:
            if child.tag == 'INSTANCE':
                for key, value in child.attrib.items():
                    if key == 'dmrole':
                        if value == '':
                            if child.tag == 'ATTRIBUTE':
                                recusive_xml_check(child,
                                                   getattr(MivotInstance,
                                                           MivotInstance._remove_model_name(
                                                               child.get('dmrole'))))
                            elif child.tag == 'INSTANCE':
                                recusive_xml_check(child, getattr(MivotInstance,
                                                                  MivotInstance._remove_model_name
                                                                  (child.get('dmrole'))))
                        else:
                            if child.tag == 'ATTRIBUTE':
                                recusive_xml_check(child, getattr(MivotInstance,
                                                                  MivotInstance._remove_model_name(
                                                                      child.get('dmrole'))))
                            elif child.tag == 'INSTANCE':
                                recusive_xml_check(child, getattr(MivotInstance,
                                                                  MivotInstance._remove_model_name(
                                                                      child.get('dmrole'))))
                            elif child.tag == 'COLLECTION':
                                recusive_xml_check(child, getattr(MivotInstance,
                                                                  MivotInstance._remove_model_name(
                                                                      child.get('dmrole'))))
            elif child.tag == 'COLLECTION':
                for key, value in child.attrib.items():
                    assert len(getattr(MivotInstance,
                                       MivotInstance._remove_model_name(child.get('dmrole')))) == len(child)
                    i = 0
                    for child2 in child:
                        recusive_xml_check(child2, getattr(MivotInstance, MivotInstance._remove_model_name
                        (child.get('dmrole')))[i])
                        i += 1
            elif child.tag == 'ATTRIBUTE':
                MivotInstance_attribute = getattr(MivotInstance,
                                                  MivotInstance._remove_model_name(child.get('dmrole')))
                for key, value in child.attrib.items():
                    if key == 'dmtype':
                        assert MivotInstance_attribute.dmtype in value
                    elif key == 'value':
                        if (MivotInstance_attribute.value is not None
                                and not isinstance(MivotInstance_attribute.value, bool)):
                            if isinstance(MivotInstance_attribute.value, float):
                                pytest.approx(float(value), MivotInstance_attribute.value, 0.0001)
                            else:
                                assert value == MivotInstance_attribute.value
            else:
                assert False


@pytest.mark.remote_data
@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_dict_model_viewer3(votable_test, simple_votable):
    """
    To test the generation of the MIVOT class, the function builds a ModelViewerLevel3
    with his MIVOT class and his previous dictionary from XML.
    Then, it calls the function recursive_check which recursively compares an element of MIVOT class
    with the dictionary on which it was built.
    MIVOT class is itself a dictionary with only essential information of the ModelViewerLevel3._dict.
    This test run on 2 votables : votable_test and simple_votable.
    """
    m_viewer_votable_test = MivotViewer(votable_path=votable_test)
    m_viewer_votable_test.next()
    mivot_instance = m_viewer_votable_test.dm_instance
    _dict = MivotUtils.xml_to_dict(m_viewer_votable_test.xml_viewer.view)
    recursive_check(mivot_instance, **_dict)
    mivot_instance = m_viewer_votable_test.dm_instance
    _dict = MivotUtils.xml_to_dict(m_viewer_votable_test.xml_view)
    recursive_check(mivot_instance, **_dict)


def recursive_check(MivotInstance, **kwargs):
    for key, value in kwargs.items():
        # the root instance ha no role: this makes an empty value in the unpacked dict
        if key == '':
            continue
        if isinstance(value, list):
            nbr_item = 0
            for item in value:
                if isinstance(item, dict):
                    assert 'dmtype' in item.keys()
                    recursive_check(getattr(MivotInstance,
                                            MivotInstance._remove_model_name(key))[nbr_item],
                                    **item
                                    )
                    nbr_item += 1
        elif isinstance(value, dict) and 'value' not in value:
            # for INSTANCE of INSTANCEs dmrole needs model_name
            assert MivotInstance._remove_model_name(key, True) in vars(MivotInstance).keys()
            recursive_check(getattr(MivotInstance, MivotInstance._remove_model_name(key, True)), **value)
        else:
            if isinstance(value, dict) and MivotInstance._is_leaf(**value):
                assert value.keys().__contains__('dmtype' and 'value' and 'unit' and 'ref')
                lower_dmtype = value['dmtype'].lower()
                if "real" in lower_dmtype or "double" in lower_dmtype or "float" in lower_dmtype:
                    assert isinstance(value['value'], float)
                elif "bool" in lower_dmtype:
                    assert isinstance(value['value'], bool)
                elif value['dmtype'] is None:
                    assert (value['value'] in
                            ('notset', 'noset', 'null', 'none', 'NotSet', 'NoSet', 'Null', 'None'))
                else:
                    if value['value'] is not None:
                        assert isinstance(value['value'], str)
                recursive_check(getattr(MivotInstance, MivotInstance._remove_model_name(key)), **value)
            else:
                assert key == 'dmtype' or 'value'


@pytest.fixture
def votable_test(data_path, data_sample_url):
    votable_name = "vizier_csc2_gal.annot.xml"
    votable_path = os.path.join(data_path, "data", votable_name)
    urlretrieve(data_sample_url + votable_name,
                votable_path)
    yield votable_path
    os.remove(votable_path)


@pytest.fixture
def simple_votable(data_path, data_sample_url):
    votable_name = "simple-annotation-votable.xml"
    votable_path = os.path.join(data_path, "data", votable_name)
    urlretrieve(data_sample_url + votable_name,
                votable_path)
    yield votable_path
    os.remove(votable_path)


@pytest.fixture
def data_path():
    return os.path.dirname(os.path.realpath(__file__))


@pytest.fixture
def data_sample_url():
    return "https://raw.githubusercontent.com/ivoa/dm-usecases/main/pyvo-ci-sample/"
