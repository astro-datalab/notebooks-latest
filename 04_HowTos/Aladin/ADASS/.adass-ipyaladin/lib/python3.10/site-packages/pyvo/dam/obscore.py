# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A module for representing data in ObsCore format

"""

__all__ = ['ObsCoreMetadata', 'POLARIZATION_STATES', 'CALIBRATION_LEVELS']


# to be moved to ObsCore
POLARIZATION_STATES = ['I', 'Q', 'U', 'V', 'RR', 'LL', 'RL', 'LR',
                       'XX', 'YY', 'XY', 'YX', 'POLI', 'POLA']
CALIBRATION_LEVELS = [0, 1, 2, 3, 4]


class ObsCoreMetadata():
    """
    Representation of an ObsCore observation

    TBD setters to do validation and unit check.
    """
    def __init__(self):

        #          OBSERVATION INFO
        self.dataproduct_type = None
        self.dataproduct_subtype = None
        self.calib_level = None

        #          TARGET INFO
        self.target_name = None
        self.target_class = None

        #           DATA DESCRIPTION
        self.obs_id = None
        self.obs_title = None
        self.obs_collection = None
        self.obs_create_date = None
        self.obs_creator_name = None
        self.obs_creator_did = None

        #          CURATION INFORMATION
        self.obs_release_date = None
        self.obs_publisher_did = None
        self.publisher_id = None
        self.bib_reference = None
        self.data_rights = None

        #            ACCESS INFORMATION
        self.access_url = None
        self.access_format = None
        self.access_estsize = None

        #            SPATIAL CHARACTERISATION
        self.s_ra = None
        self.s_dec = None
        self.s_fov = None
        self.s_region = None
        self.s_resolution = None
        self.s_xel1 = None
        self.s_xel2 = None
        self.s_ucd = None
        self.s_unit = None
        self.s_resolution_min = None
        self.s_resolution_max = None
        self.s_calib_status = None
        self.s_stat_error = None
        self.s_pixel_scale = None

        #            TIME CHARACTERISATION
        self.t_xel = None
        self.t_ref_pos = None
        self.t_min = None
        self.t_max = None
        self.t_exptime = None
        self.t_resolution = None
        self.t_calib_status = None
        self.t_stat_error = None

        #            SPECTRAL CHARACTERISATION
        self.em_xel = None
        self.em_ucd = None
        self.em_unit = None
        self.em_calib_status = None
        self.em_min = None
        self.em_max = None
        self.em_res_power = None
        self.em_res_power_min = None
        self.em_res_power_max = None
        self.em_resolution = None
        self.em_stat_error = None

        #            OBSERVABLE AXIS
        self.o_ucd = None
        self.o_unit = None
        self.o_calib_status = None
        self.o_stat_error = None

        #            POLARIZATION CHARACTERISATION
        self.pol_xel = None
        self.pol_states = None

        #            PROVENANCE
        self.instrument_name = None
        self.facility_name = None
        self.proposal_id = None
