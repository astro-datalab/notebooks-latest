# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test for mivot.seekers.resource_seeker.py
"""
import pytest
from astropy.io.votable import parse
from astropy.utils.data import get_pkg_data_filename
from pyvo.mivot.seekers.resource_seeker import ResourceSeeker
from pyvo.mivot.version_checker import check_astropy_version


@pytest.mark.skipif(not check_astropy_version(), reason="need astropy 6+")
def test_id_table(rseeker):
    """
    Checks the IDs of tables found by the RessourceSeeker,
    checks the IDs of the field of the table concerned.
    """
    assert rseeker.get_table_ids() == ['_PKTable', 'Results']
    assert (rseeker.get_id_index_mapping('_PKTable')
            == {'pksrcid': {'ID': '_pksrcid', 'indx': 0},
                'pkband': {'ID': '_pkband', 'indx': 1}
                }
            )
    assert (rseeker.get_id_index_mapping('Results')
            == {'source_id': {'ID': '_srcid', 'indx': 0},
                'transit_id': {'ID': 'transit_id', 'indx': 1},
                'band': {'ID': '_band', 'indx': 2},
                'time': {'ID': '_obstime', 'indx': 3},
                'mag': {'ID': '_mag', 'indx': 4},
                'flux': {'ID': '_flux', 'indx': 5},
                'flux_error': {'ID': '_fluxerr', 'indx': 6},
                'flux_over_error': {'ID': 'flux_over_error', 'indx': 7},
                'rejected_by_photometry': {'ID': 'rejected_by_photometry', 'indx': 8},
                'rejected_by_variability': {'ID': 'rejected_by_variability', 'indx': 9},
                'other_flags': {'ID': 'other_flags', 'indx': 10},
                'solution_id': {'ID': 'solution_id', 'indx': 11}
                }
            )
    table = rseeker.get_table('_PKTable')
    for field in table.fields:
        field.ID = None
    assert (rseeker.get_id_index_mapping('_PKTable')
            == {'pksrcid': {'indx': 0, 'ID': 'pksrcid'},
                'pkband': {'indx': 1, 'ID': 'pkband'}
                }
            )
    for table in rseeker._resource.tables:
        table.ID = None
    assert rseeker.get_table_ids() == ['AnonymousTable', 'AnonymousTable']

    for table in rseeker._resource.tables:
        table.name = "any_name"

    assert rseeker.get_table_ids() == ['any_name', 'any_name']
    assert (rseeker.get_id_index_mapping('any_name')
            == {'pksrcid': {'indx': 0, 'ID': 'pksrcid'},
                'pkband': {'indx': 1, 'ID': 'pkband'}}
            )


@pytest.fixture
def rseeker():

    votable_path = get_pkg_data_filename("data/test.mivot_viewer.xml")

    votable = parse(votable_path)
    for resource in votable.resources:
        return ResourceSeeker(resource)
