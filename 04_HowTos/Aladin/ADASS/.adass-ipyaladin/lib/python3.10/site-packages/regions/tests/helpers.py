# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from astropy.wcs import WCS
from numpy.testing import assert_allclose

from regions.core import PixCoord, Region


def make_simple_wcs(skycoord, resolution, size):
    crpix = (size + 1) / 2
    cdelt = resolution.to(u.deg).value
    skycoord_icrs = skycoord.transform_to('icrs')
    ra = skycoord_icrs.ra.degree
    dec = skycoord_icrs.dec.degree

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [crpix, crpix]
    wcs.wcs.cdelt = np.array([-cdelt, cdelt])
    wcs.wcs.crval = [ra, dec]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    return wcs


def assert_skycoord_allclose(skycoord1, skycoord2, **kwargs):
    """
    Test that two SkyCoord objects are nearly equal.
    """
    for attr in (skycoord1._extra_frameattr_names
                 | skycoord2._extra_frameattr_names):
        if not skycoord1.frame._frameattr_equiv(getattr(skycoord1, attr),
                                                getattr(skycoord2, attr)):
            raise ValueError('Extra frame attributes are not equivalent')

    assert skycoord1.is_equivalent_frame(skycoord2)

    for attr in (set(skycoord1._sky_coord_frame._data.components)
                 | set(skycoord2._sky_coord_frame._data.components)):
        sdata1 = skycoord1._sky_coord_frame._data
        sdata2 = skycoord2._sky_coord_frame._data
        assert_quantity_allclose(getattr(sdata1, attr), getattr(sdata2, attr),
                                 **kwargs)


def assert_region_allclose(region1, region2, **kwargs):
    """
    Test that two Region objects have parameters which are nearly equal.

    Meta and Visual properties are matched identically.
    """
    if not (isinstance(region1, Region) and isinstance(region2, Region)):
        raise TypeError('Both inputs must be Region instances')

    # check that both have identical parameters
    region1_params = list(region1._params)
    region2_params = list(region2._params)
    if region1_params != region2_params:
        raise ValueError('Inputs do not have the same parameters')

    assert region1.meta == region2.meta
    assert region1.visual == region2.visual

    # now check the parameter values
    # Note that Quantity comparisons allow for different units
    # if they directly convertible (e.g., 1. * u.deg == 60. * u.arcmin)
    for param in region1_params:
        value1 = getattr(region1, param)
        value2 = getattr(region2, param)
        if isinstance(value1, SkyCoord):
            assert_skycoord_allclose(value1, value2, **kwargs)
        elif isinstance(value1, PixCoord):
            assert_allclose(value1.xy, value2.xy)
        elif isinstance(value1, u.Quantity):
            assert_quantity_allclose(value1, value2, **kwargs)
        elif isinstance(value1, str):
            assert value1 == value2
        else:
            assert_allclose(value1, value2)
