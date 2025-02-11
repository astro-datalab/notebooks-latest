# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the crtf subpackage.
"""

import astropy.units as u
import pytest
from astropy.coordinates import Angle, SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filename

from regions.core import Regions
from regions.io.crtf.core import CRTFRegionParserError
from regions.io.crtf.read import _CRTFParser
from regions.shapes.circle import CircleSkyRegion
from regions.shapes.ellipse import EllipseSkyRegion

implemented_region_types = ('ellipse', 'circle', 'rectangle', 'poly', 'point',
                            'text', 'symbol')


def test_global_parser():
    """
    Checks that the global_parser does what's expected.
    """
    global_test_str = ('global coord=B1950_VLA, frame=BARY, corr=[I, Q], '
                       'color=blue')
    global_parser = _CRTFParser(global_test_str)
    expected = {'coord': 'B1950_VLA', 'frame': 'BARY',
                'corr': ['I', 'Q'], 'color': 'blue'}
    assert dict(global_parser.global_meta) == expected


def test_valid_crtf_line():
    """
    Checks whether a the line is valid CRTF format.
    """
    line_str = 'coord=B1950_VLA, frame=BARY, corr=[I, Q], color=blue'

    with pytest.raises(CRTFRegionParserError) as excinfo:
        Regions.parse(line_str, format='crtf')

    assert 'Not a valid CRTF line:' in str(excinfo.value)


def test_valid_region_type():
    """
    Checks whether the region type is valid in CRTF format.
    """
    reg_str = 'hyperbola[[18h12m24s, -23d11m00s], 2.3arcsec]'

    with pytest.raises(CRTFRegionParserError) as excinfo:
        Regions.parse(reg_str, format='crtf')

    assert 'Not a valid CRTF Region type: "hyperbola"' in str(excinfo.value)


def test_valid_global_meta_key():
    """
    Checks whether the global key is valid or not.
    """
    global_test_str = ('global label=B1950_VLA, frame=BARY, corr=[I, Q], '
                       'color=blue')
    with pytest.raises(CRTFRegionParserError) as excinfo:
        Regions.parse(global_test_str, format='crtf')
    assert '"label" is not a valid global meta key' in str(excinfo.value)


def test_valid_meta_key():
    """
    Checks whether the key is valid or not.
    """
    meta_test_str = ('annulus[[17h51m03.2s, -45d17m50s], [0.10deg, 4.12deg]], '
                     'hello="My label here"')
    with pytest.raises(CRTFRegionParserError) as excinfo:
        Regions.parse(meta_test_str, format='crtf')
    assert '"hello" is not a valid meta key' in str(excinfo.value)


def test_valid_region_syntax():
    """
    Checks whether the region has valid parameters.
    """
    reg_str1 = 'circle[[18h12m24s, -23d11m00s], [2.3arcsec,4.5arcsec]'
    with pytest.raises(CRTFRegionParserError) as excinfo:
        Regions.parse(reg_str1, format='crtf')

    estr = ("Not in proper format: ('2.3arcsec', '4.5arcsec') should be "
            'a single length')
    assert estr in str(excinfo.value)

    reg_str2 = ('symbol[[32.1423deg, 12.1412deg], 12deg], linewidth=2, '
                'coord=J2000, symsize=2')
    with pytest.raises(CRTFRegionParserError) as excinfo:
        Regions.parse(reg_str2, format='crtf')
    estr = 'Not in proper format: "12deg" should be a symbol'
    assert estr in str(excinfo.value)

    reg_str3 = 'circle[[18h12m24s, -23d11m00s]'
    with pytest.raises(CRTFRegionParserError) as excinfo:
        Regions.parse(reg_str3, format='crtf')
    estr = ('Does not contain expected number of parameters for the region '
            '"circle"')
    assert estr in str(excinfo.value)

    reg_str4 = 'poly[[1, 2], [4, 5]]'
    with pytest.raises(CRTFRegionParserError) as excinfo:
        Regions.parse(reg_str4, format='crtf')
    assert 'polygon should have >= 3 coordinates' in str(excinfo.value)

    reg_str6 = 'rotbox[[12h01m34.1s, 12d23m33s], [3arcmin,], 12deg]'
    with pytest.raises(CRTFRegionParserError) as excinfo:
        Regions.parse(reg_str6, format='crtf')
    assert "('3arcmin', '') should be a pair of lengths" in str(excinfo.value)


def test_issue_312_regression():
    """
    Make sure there is no trailing comma when writing a CRTF string with
    no metadata.
    """
    reg = EllipseSkyRegion(center=SkyCoord(279.174990 * u.deg,
                                           -21.257123 * u.deg, frame='fk5'),
                           width=0.001571 * u.deg, height=0.001973 * u.deg,
                           angle=111.273322 * u.deg)
    crtfstr = reg.serialize(format='crtf', coordsys='fk5', fmt='.6f',
                            radunit='deg')
    assert crtfstr.strip()[-1] != ','


@pytest.mark.parametrize(('filename', 'outname', 'coordsys', 'fmt'),
                         [('data/CRTFgeneral.crtf',
                           'data/CRTFgeneraloutput.crtf',
                           'fk4', '.3f'),
                          ('data/CRTFgeneraloutput.crtf',
                           'data/CRTFgeneraloutput.crtf',
                           'fk4', '.3f'),
                          ('data/CRTF_labelcolor.crtf',
                           'data/CRTF_labelcolor_output.crtf',
                           'fk5', '.6f')])
def test_file_crtf(filename, outname, coordsys, fmt):
    """
    The "labelcolor" example is a regression test for Issue 405 The
    others are just a general serialization self-consistency check.
    """
    filename = get_pkg_data_filename(filename)
    regs = Regions.read(filename, errors='warn', format='crtf')
    actual_output = regs.serialize(format='crtf', coordsys=coordsys,
                                   fmt=fmt).strip()

    with open(get_pkg_data_filename(outname)) as fh:
        ref_output = fh.read().strip()

    # since metadata is not required to preserve order, we have to do a more
    # complex comparison
    desired_lines = [set(line.split(',')) for line in ref_output.split('\n')]
    actual_lines = [set(line.split(',')) for line in actual_output.split('\n')]
    for split_line in actual_lines:
        assert split_line in desired_lines

    for split_line in desired_lines:
        assert split_line in actual_lines


def test_crtf_header():
    """
    Regression test for issue #239.
    """
    crtf_str = ('#CRTFv0 CASA Region Text Format version 0\n'
                'circle[[42deg, 43deg], 3deg], coord=J2000, color=green')
    reg = Regions.parse(crtf_str, format='crtf')[0]
    assert isinstance(reg, CircleSkyRegion)
    assert reg.center.ra.value == 42.0
    assert reg.center.ra.unit == 'deg'
    assert reg.center.dec.value == 43.0
    assert reg.center.dec.unit == 'deg'
    assert reg.radius.value == 3.0
    assert reg.radius.unit == 'deg'


def test_space_after_regname():
    """
    Regression test for #271: space is allowed.
    """
    reg_str = 'circle [[42deg, 43deg], 3deg], coord=J2000, color=green'
    reg = Regions.parse(reg_str, format='crtf')[0]
    assert isinstance(reg, CircleSkyRegion)


def test_no_comma_after_region():
    reg_str = 'circle [[42deg, 43deg], 3deg] coord=J2000, color=green'
    reg = Regions.parse(reg_str, format='crtf')[0]
    assert isinstance(reg, CircleSkyRegion)
    assert reg.center.ra.value == 42.0
    assert reg.center.ra.unit == 'deg'
    assert reg.center.dec.value == 43.0
    assert reg.center.dec.unit == 'deg'
    assert reg.radius.value == 3.0
    assert reg.radius.unit == 'deg'


def test_casa_file_crtf():
    filename = get_pkg_data_filename('data/CRTF_CARTA.crtf')
    regions = Regions.read(filename, format='crtf')
    assert len(regions) == 2


def test_angle_serialization():
    """
    Regression test for issue #223 to ensure Angle arcsec inputs are
    correctly converted to degrees.
    """
    reg = Regions([CircleSkyRegion(SkyCoord(10, 20, unit='deg'),
                                   Angle(1, 'arcsec'))])
    regstr = reg.serialize(format='crtf')
    expected = ('#CRTFv0\nglobal coord=J2000\ncircle[[10.000009deg, '
                '20.000002deg], 0.000278deg]\n')
    assert regstr == expected


def test_read_sexagesimal_regression407():
    filename = 'data/crtf_carta_sexagesimal.crtf'
    filename = get_pkg_data_filename(filename)
    regs = Regions.read(filename, errors='warn', format='crtf')

    assert_quantity_allclose(regs[0].center.ra, 281.93342096 * u.deg, rtol=1e-9)
