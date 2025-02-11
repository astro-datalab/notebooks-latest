#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for pyvo.dal.sia2 against remote services
"""

import pytest

import astropy.units as u
from astropy.utils.exceptions import AstropyDeprecationWarning

from pyvo.dal.sia2 import search, SIA2Service
from pyvo.dal.adhoc import DatalinkResults
from pyvo import regsearch


CADC_SIA_URL = 'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/sia'


@pytest.mark.remote_data
class TestSIACadc():
    # Tests the SIA2 client against the CADC SIA service

    def test_service(self):
        cadc = SIA2Service(baseurl=CADC_SIA_URL)

        with pytest.raises(AstropyDeprecationWarning):
            assert cadc.availability

        with pytest.raises(AstropyDeprecationWarning):
            assert cadc.availability.available

        with pytest.raises(AstropyDeprecationWarning):
            assert cadc.availability.notes

        with pytest.raises(AstropyDeprecationWarning):
            assert cadc.availability.notes[0] == 'service is accepting queries'

        assert cadc.capabilities

    @pytest.mark.xfail(reason="https://github.com/astropy/pyvo/issues/361")
    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_datalink_batch(self):
        # Maximum batch size in CADC SIA is around 25
        # Test whether multiple batches can be retrieved
        results = search(CADC_SIA_URL, pos=(2.8425, 74.4846, 10), maxrec=55)
        ids = []
        for i in results.iter_datalinks():
            assert i.to_table()[0]['ID'] not in ids
            ids.append(i.to_table()[0]['ID'])
        assert len(ids) == 55

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_pos(self):
        results = search(CADC_SIA_URL, pos=(2.8425, 74.4846, 0.001))
        assert len(results) > 10

        # check that results are datalink
        assert isinstance(results[0].getdataobj(), DatalinkResults)

        # limit results to 5 to expedite tests
        results = search(CADC_SIA_URL, pos=(2.8425, 74.4846, 0.001), maxrec=5)
        assert len(results) == 5
        # check attributes of a record
        record = results[0]
        record.dataproduct_type
        record.dataproduct_subtype
        record.calib_level

        #          TARGET INFO
        record.target_name
        record.target_class

        #           DATA DESCRIPTION
        record.obs_id
        record.obs_title
        record.obs_collection
        record.obs_create_date
        record.obs_creator_name
        record.obs_creator_did

        #          CURATION INFORMATION
        record.obs_release_date
        record.obs_publisher_did
        record.publisher_id
        record.bib_reference
        record.data_rights

        #            ACCESS INFORMATION
        record.access_url
        record.access_format
        record.access_estsize

        #            SPATIAL CHARACTERISATION
        record.s_ra
        record.s_dec
        record.s_fov
        record.s_region
        record.s_resolution
        record.s_xel1
        record.s_xel2
        record.s_ucd
        record.s_unit
        record.s_resolution_min
        record.s_resolution_max
        record.s_calib_status
        record.s_stat_error
        record.s_pixel_scale

        #            TIME CHARACTERISATION
        record.t_xel
        record.t_ref_pos
        record.t_min
        record.t_max
        record.t_exptime
        record.t_resolution
        record.t_calib_status
        record.t_stat_error

        #            SPECTRAL CHARACTERISATION
        record.em_xel
        record.em_ucd
        record.em_unit
        record.em_calib_status
        record.em_min
        record.em_max
        record.em_res_power
        record.em_res_power_min
        record.em_res_power_max
        record.em_resolution
        record.em_stat_error

        #            OBSERVABLE AXIS
        record.o_ucd
        record.o_unit
        record.o_calib_status
        record.o_stat_error

        #            POLARIZATION CHARACTERISATION
        record.pol_xel
        record.pol_states

        #            PROVENANCE
        record.instrument_name
        record.facility_name
        record.proposal_id

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_band(self):
        results = search(CADC_SIA_URL, band=(0.0002, 0.0003), maxrec=5)
        # TODO - correctness
        assert len(results) == 5

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_time(self):
        results = search(CADC_SIA_URL,
                         time=('2002-01-01T00:00:00.00',
                               '2002-01-02T00:00:00.00'),
                         maxrec=5)
        assert len(results) == 5

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_pol(self):
        results = search(CADC_SIA_URL, pol=['YY', 'U'], maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert 'YY' in rr.pol_states or 'U' in rr.pol_states

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_fov(self):
        results = search(CADC_SIA_URL, field_of_view=(10, 20), collection='TESS', maxrec=5)
        assert len(results) == 5
        # how to test values

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_spatial_res(self):
        results = search(CADC_SIA_URL, spatial_resolution=(1, 2), collection='MACHO', maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert 1 * u.arcsec <= rr.s_resolution <= 2 * u.arcsec

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_spec_resp(self):
        results = search(CADC_SIA_URL, spectral_resolving_power=(1, 2), collection='NEOSSAT', maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert 1 <= rr.em_res_power <= 2

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_exptime(self):
        results = search(CADC_SIA_URL, exptime=(1, 2), collection='GEMINI', maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert 1 * u.second <= rr.t_exptime <= 2 * u.second

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_timeres(self):
        results = search(CADC_SIA_URL, timeres=(1, 2),
                         maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert 1 * u.second <= rr.t_resolution <= 2 * u.second

    def test_publisher_did(self):
        ids = ['ivo://cadc.nrc.ca/CFHT?447231/447231o',
               'ivo://cadc.nrc.ca/CFHT?447232/447232o']
        results = search(CADC_SIA_URL, publisher_did=ids)
        assert len(results) == 2
        assert results[0].obs_publisher_did in ids
        assert results[1].obs_publisher_did in ids

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_facility(self):
        results = search(CADC_SIA_URL, facility='JCMT', maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert rr.facility_name == 'JCMT'

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_collection(self):
        results = search(CADC_SIA_URL, collection='IRIS', maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert rr.obs_collection == 'IRIS'

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_instrument(self):
        results = search(CADC_SIA_URL, instrument='IRAS', maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert rr.instrument_name == 'IRAS'

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_dataproduct_type(self):
        results = search(CADC_SIA_URL, data_type='image', maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert rr.dataproduct_type == 'image'

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_target_name(self):
        results = search(CADC_SIA_URL, target_name='BarnardStar', collection='NEOSSAT', maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert rr.target_name == 'BarnardStar'

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_res_format(self):
        results = search(
            CADC_SIA_URL,
            res_format='application/x-votable+xml;content=datalink',
            maxrec=5)
        assert len(results) == 5
        for rr in results:
            assert rr.access_format == \
                'application/x-votable+xml;content=datalink'

    @pytest.mark.filterwarnings("ignore::pyvo.dal.exceptions.DALOverflowWarning")
    def test_reg_sia2(self):
        image_services = regsearch(servicetype='sia2')
        irsa_seip = \
            [s for s in image_services if
             'irsa' in s.ivoid and 'seip' in s.ivoid][0]
        result = irsa_seip.search(pos=(31.8425, 77.4846, 0.1), maxrec=1)
        assert len(result) == 1
