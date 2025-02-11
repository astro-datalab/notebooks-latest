# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test for mivot.features.static_reference_resolver.py
"""
import pytest
from astropy.utils.data import get_pkg_data_filename
from pyvo.mivot.seekers.annotation_seeker import AnnotationSeeker
from pyvo.mivot.features.static_reference_resolver import StaticReferenceResolver
from pyvo.mivot.version_checker import check_astropy_version
from pyvo.mivot import MivotViewer
from . import XMLOutputChecker


@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_static_reference_resolve(a_seeker, instance):
    StaticReferenceResolver.resolve(a_seeker, None, instance)
    XMLOutputChecker.assertXmltreeEqualsFile(
        instance.getroot(),
        get_pkg_data_filename("data/reference/static_reference_resolved.xml")
    )


@pytest.fixture
def instance():
    return XMLOutputChecker.xmltree_from_file(
        get_pkg_data_filename("data/static_reference.xml"))


@pytest.fixture
def a_seeker():
    m_viewer = MivotViewer(
        get_pkg_data_filename("data/test.mivot_viewer.xml"),
        tableref="Results")
    return AnnotationSeeker(m_viewer._mapping_block)
