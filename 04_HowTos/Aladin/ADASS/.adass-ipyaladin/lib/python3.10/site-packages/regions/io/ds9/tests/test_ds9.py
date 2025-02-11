# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the ds9 subpackage.
"""

import os
import warnings

import astropy.units as u
import pytest
from astropy.coordinates import Angle, SkyCoord
from astropy.tests.helper import assert_quantity_allclose
from astropy.utils.data import get_pkg_data_filenames
from astropy.utils.exceptions import AstropyUserWarning
from numpy.testing import assert_allclose, assert_equal

from regions._utils.optional_deps import HAS_MATPLOTLIB
from regions.core import PixCoord, Regions, RegionVisual
from regions.shapes import (CirclePixelRegion, CircleSkyRegion,
                            PointPixelRegion, RegularPolygonPixelRegion,
                            TextPixelRegion)
from regions.tests.helpers import assert_region_allclose


def test_roundtrip(tmpdir):
    filenames = get_pkg_data_filenames('data', pattern='*.reg')

    # AstropyUserWarning will be emitted only for some of the files
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyUserWarning)
        # Check that all test files are readable
        for filename in filenames:
            regions = Regions.read(filename, format='ds9')
            assert len(regions) > 0

            tempfile = tmpdir.join('tmp.ds9').strpath
            regions.write(tempfile, format='ds9', overwrite=True, precision=20)
            regions2 = Regions.read(tempfile, format='ds9')
            assert len(regions2) > 0
            for reg1, reg2 in zip(regions, regions2, strict=True):
                assert_region_allclose(reg1, reg2)


def test_serialize():
    """
    Simple test for serialize.
    """
    center = SkyCoord(42, 43, unit='deg')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)
    expected = ('# Region file format: DS9 astropy/regions\nicrs\n'
                'circle(42.0000,43.0000,3.0000)\n')
    actual = region.serialize(format='ds9', precision=4)
    assert actual == expected


def test_serialize_parse_text():
    """
    Test serialization of Text region.
    """
    text = 'Example Text'
    region = TextPixelRegion(PixCoord(42, 43), text=text)
    expected = ('# Region file format: DS9 astropy/regions\nglobal '
                f'text={{{text}}}\nimage\ntext(43.0000,44.0000)\n')
    actual = region.serialize(format='ds9', precision=4)
    assert actual == expected

    reg = Regions.parse(actual, format='ds9')[0]
    assert reg.text == text
    assert 'text' not in reg.meta


def test_parse():
    """
    Simple test for parse.
    """
    region_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
                  'circle(42.0000,43.0000,3.0000)\n')
    region = Regions.parse(region_str, format='ds9')[0]

    assert_allclose(region.center.ra.deg, 42)
    assert_allclose(region.center.dec.deg, 43)
    assert_quantity_allclose(region.radius, 3 * u.deg)


def test_exclude():
    """
    Simple parse test for an excluded region.
    """
    region_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
                  ' -circle(42.0000,43.0000,3.0000)\n')
    region = Regions.parse(region_str, format='ds9')[0]
    assert not region.meta['include']


def test_read_write(tmpdir):
    """
    Simple test for write and read.
    """
    center = SkyCoord(42, 43, unit='deg', frame='fk5')
    radius = Angle(3, 'deg')
    region = CircleSkyRegion(center, radius)
    region.meta['text'] = 'ExampleText'

    filename = os.path.join(str(tmpdir), 'ds9.reg')
    region.write(filename, format='ds9')
    region2 = Regions.read(filename, format='ds9')[0]

    assert_quantity_allclose(region.center.ra, 42 * u.deg)
    assert_quantity_allclose(region.center.dec, 43 * u.deg)
    assert_quantity_allclose(region.radius, 3 * u.deg)
    assert 'text' in region2.meta
    assert region2.meta['text'] == 'ExampleText'


def test_invalid_region_warns():
    ds9_str = ('# Region file format: DS9 astropy/regions\nfk5\n'
               'circle(42.0000,43.0000,3.0000)\ninvalidregion(blah)')

    with pytest.warns(AstropyUserWarning) as warn_results:
        regions = Regions.parse(ds9_str, format='ds9')
    assert len(regions) == 1
    assert len(warn_results) == 1


def test_global_parser():
    """
    Test parsing global metadata.
    """
    region_str = ('global color=green dashlist=8 3 width=1 '
                  'font="helvetica 10 normal roman" select=1 highlite=1 '
                  'dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
                  'icrs;circle(42.0000,43.0000,3.0000)')
    region = Regions.parse(region_str, format='ds9')[0]

    meta_expected = {'select': 1, 'highlite': 1, 'fixed': 0, 'edit': 1,
                     'move': 1, 'delete': 1, 'include': 1, 'source': 1}
    visual_expected = {'linewidth': 1, 'fontname': 'helvetica',
                       'fontsize': 10, 'fontweight': 'normal',
                       'fontstyle': 'normal', 'facecolor': 'green',
                       'edgecolor': 'green', 'default_style': 'ds9'}

    assert region.meta == meta_expected
    assert region.visual == visual_expected


def test_meta_color():
    region_str = ('# Region file format: DS9 astropy/regions\nicrs\n'
                  'circle(42.0000,43.0000,3.0000) # color=green\n'
                  'circle(43.0000,43.0000,3.0000) # color=orange\n')
    regions = Regions.parse(region_str, format='ds9')

    assert regions[0].visual['facecolor'] == 'green'
    assert regions[0].visual['edgecolor'] == 'green'
    assert regions[1].visual['facecolor'] == 'orange'
    assert regions[1].visual['edgecolor'] == 'orange'


def test_meta_color_override_global():
    """
    Color parsing test in the presence of a global color.
    """
    region_str = ('# Region file format: DS9 astropy/regions\n'
                  'global color=blue dashlist=8 3 width=1 '
                  'font="helvetica 10 normal roman" select=1 highlite=1 '
                  'dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
                  'icrs\ncircle(42.0000,43.0000,3.0000) # color=green\n'
                  'circle(42.0000,43.0000,5.0000) # color=orange\n'
                  'circle(42.0000,43.0000,7.0000)')
    regions = Regions.parse(region_str, format='ds9')
    assert regions[0].visual['facecolor'] == 'green'
    assert regions[1].visual['facecolor'] == 'orange'
    assert regions[2].visual['facecolor'] == 'blue'


def test_issue134_regression():
    region_str = 'galactic; circle(+0:14:26.064,+0:00:45.206,30.400")'
    region = Regions.parse(region_str, format='ds9')[0]
    assert_quantity_allclose(region.radius, 30.4 * u.arcsec)


def test_issue65_regression():
    region_str = 'J2000; circle 188.5557102 12.0314056 1" # color=red'
    region = Regions.parse(region_str, format='ds9')[0]
    assert_quantity_allclose(region.center.ra, 188.5557102 * u.deg)
    assert_quantity_allclose(region.center.dec, 12.0314056 * u.deg)
    assert_quantity_allclose(region.radius, 1 * u.arcsec)


def test_frame_uppercase():
    """
    Regression test for issue #236 (PR #237).
    """
    region_str = ('GALACTIC\ncircle(188.5557102,12.0314056,7.1245) # '
                  'text={This message has both a " and : in it} '
                  'textangle=30')
    region = Regions.parse(region_str, format='ds9')[0]
    assert_quantity_allclose(region.center.l, 188.5557102 * u.deg)
    assert_quantity_allclose(region.center.b, 12.0314056 * u.deg)
    assert_quantity_allclose(region.radius, 7.1245 * u.deg)


def test_pixel_angle():
    """
    Check whether angle in PixelRegions is a u.Quantity object.
    """
    region_str = 'image\nbox(1.5,2,2,1,0)'
    region = Regions.parse(region_str, format='ds9')[0]
    assert isinstance(region.angle, u.Quantity)
    assert_quantity_allclose(region.angle, 0. * u.deg)


def test_parse_formats():
    region1_str = 'image\ncircle(1.5, 2, 2)'
    region1 = Regions.parse(region1_str, format='ds9')[0]
    region2_str = 'image\ncircle(1.5i, 2i, 2i)'
    region2 = Regions.parse(region2_str, format='ds9')[0]
    assert region1 == region2

    region3_str = 'icrs\ncircle(1.5, 2, 2)'
    region3 = Regions.parse(region3_str, format='ds9')[0]
    region4_str = 'icrs\ncircle(1.5d, 2d, 2d)'
    region4 = Regions.parse(region4_str, format='ds9')[0]
    assert region3 == region4

    region5_str = 'icrs\ncircle(1.5r, 1r, 2r)'
    region5 = Regions.parse(region5_str, format='ds9')[0]
    assert_quantity_allclose(region5.center.ra, 1.5 * u.radian)
    assert_quantity_allclose(region5.center.dec, 1.0 * u.radian)
    assert_quantity_allclose(region5.radius, 2.0 * u.radian)

    region6_str = 'icrs\ncircle(1:20:30, 2:3:7, 2)'
    region6 = Regions.parse(region6_str, format='ds9')[0]
    assert_quantity_allclose(region6.center.ra, 20.124999999999996 * u.deg)
    assert_quantity_allclose(region6.center.dec, 2.051944444444444 * u.deg)
    assert_quantity_allclose(region6.radius, 2. * u.deg)

    region7_str = 'icrs\ncircle(1h20m30s, 2d3m7s, 2d)'
    region7 = Regions.parse(region7_str, format='ds9')[0]
    assert_quantity_allclose(region7.center.ra, 20.124999999999996 * u.deg)
    assert_quantity_allclose(region7.center.dec, 2.051944444444444 * u.deg)
    assert_quantity_allclose(region7.radius, 2. * u.deg)


def test_angle_serialization():
    """
    Regression test for issue #223 to ensure Angle arcsec inputs are
    correctly converted to degrees.
    """
    region = CircleSkyRegion(SkyCoord(10, 20, unit='deg'), Angle(1, 'arcsec'))
    region_str = region.serialize(format='ds9', precision=6)

    expected = ('# Region file format: DS9 astropy/regions\nicrs\n'
                'circle(10.000000,20.000000,0.000278)\n')
    assert region_str == expected


def test_semicolon():
    """
    Regression test for issue #238 to allow semicolons in the text
    field.
    """
    region_str = ('galactic\n'
                  'circle(0.003, 0.1, 206.696") # text={S17; test} color=red\n'
                  'circle(0.1, -0.5, 360.148") # text={S19} color=red')

    regions = Regions.parse(region_str, format='ds9')
    assert len(regions) == 2

    regstr = ('galactic;'
              'circle(202.4,47.2,10.9) # text={test; test; A1} color=red;'
              'circle(103.8,40.6,3.6) # text="A, \'B\', ; {C}, D"'
              ' color=orange;'
              'circle(103.8,30.6,3.6) # text="A, \'B, ; C}, D" color=orange;'
              'circle(303.8,40.6,3.6) # text=\'A, "B", ; {C}, D\' color=cyan;'
              'circle(303.8,40.6,3.11);'
              'circle(303.8,30.6,3.6) # text=\'A, "B, ; {C, D\' color=cyan;'
              'circle(503.8,40.6,3.6) # text={A, "B", ; \'C\', D}'
              ' color=yellow\n'
              'circle(503.8,30.6,3.6) # text={A, "B, ; \'C, D} color=yellow;'
              '-ellipse(104.84,68.55,7.93",3.96",158.39) # dashlist=8 3'
              ' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1'
              ' source=1 text={Ellipse;Text} background color=#0ff width=1'
              ' font="helvetica 10 normal italic";'
              "ellipse(104.84,68.55,7.93',3.96',158.39) # dashlist=8 3"
              ' text="hello; world;";'
              'circle(1.5, 2, 2) # text={this_is_text};'
              'text(151.1,51.2) # text={text ; hi ; there;} textangle=45'
              ' color=orange font="helvetica 12 normal roman";'
              'text(151.1,51.2) # text  =  {text ; hi ; there;} textangle=45'
              ' color=orange font="helvetica 12 normal roman";'
              'circle(1.5, 2, 2) # text={text="hi;world";text;"tes\'t};'
              "circle(1.5, 2, 2) # text={text='hi;world';text;\"tes't};")

    regions = Regions.parse(regstr, format='ds9')
    assert len(regions) == 15


def test_parser_no_metadata():
    """
    Regression test for issue #259 to ensure regions without metadata
    are parsed correctly.
    """
    region1_str = 'galactic;circle(42,43,3)'
    region2_str = 'galactic;circle 42 43 3'
    region1 = Regions.parse(region1_str, format='ds9')[0]
    region2 = Regions.parse(region2_str, format='ds9')[0]

    assert isinstance(region1, CircleSkyRegion)
    assert isinstance(region2, CircleSkyRegion)
    assert_quantity_allclose(region1.radius, 3.0 * u.deg)
    assert_quantity_allclose(region2.radius, 3.0 * u.deg)


def test_spaces_metadata():
    """
    Regression test for spaces before or after "=" in region metadata.
    """
    regstr = ('galactic;'
              'circle(202.4,47.2,10.9) # text={a; test; test; A1} color=red;'
              'circle(202.4,47.2,10.9) # text ={a; test; test; A1} color= red;'
              'circle(202.4,47.2,10.9) # text= {a; test; test; A1} color =red;'
              'circle(202.4,47.2,10.9) # text = {a; test; test; A1}'
              ' color = red;')
    regions = Regions.parse(regstr, format='ds9')
    assert len(regions) == 4
    for i in range(1, len(regions)):
        assert regions[0].visual == regions[i].visual
        assert regions[0].meta == regions[i].meta

    regstr2 = ('galactic;'
               'text(151.1,51.2) # text={text ; hi ; there;} textangle= 45'
               ' color=orange font="helvetica 12 normal roman";'
               'text(151.1,51.2) # text  ={text ; hi ; there;} textangle =45'
               ' color=orange font="helvetica 12 normal roman";'
               'text(151.1,51.2) # text=  {text ; hi ; there;} textangle  =45'
               ' color   =   orange font="helvetica 12 normal roman";'
               'text(151.1,51.2) # text   =  {text ; hi ; there;} textangle=45'
               ' color   =orange font="helvetica 12 normal roman";')

    regions = Regions.parse(regstr2, format='ds9')
    assert len(regions) == 4
    for i in range(1, len(regions)):
        assert regions[0].visual == regions[i].visual
        assert regions[0].meta == regions[i].meta


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
def test_point_boxcircle():
    import matplotlib.path as mpath

    region_str = ('# Region file format: DS9 astropy/regions\n'
                  'image\n'
                  'point(101.000000,101.000000) # color=red\n'
                  'point(101.000000,201.000000) # point=boxcircle color=blue\n'
                  'point(101.000000,301.000000)')
    regions = Regions.parse(region_str, format='ds9')

    assert isinstance(regions[0].as_artist().get_marker(), mpath.Path)
    assert regions[0].as_artist().get_markeredgecolor() == 'red'
    assert isinstance(regions[1].visual['marker'], mpath.Path)
    assert regions[1].as_artist().get_markeredgecolor() == 'blue'
    assert isinstance(regions[2].as_artist().get_marker(), mpath.Path)
    assert regions[2].as_artist().get_markeredgecolor() == '#00ff00'


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
def test_compound_color():
    region_str = ('# Region file format: DS9 astropy/regions\n'
                  'image\n'
                  'annulus(651.0,301.0,60.0,90.0) # color=red')
    regions = Regions.parse(region_str, format='ds9')
    assert regions[0].as_artist().get_edgecolor() == (1., 0., 0., 1.)


@pytest.mark.skipif(not HAS_MATPLOTLIB, reason='matplotlib is required')
def test_default_mpl_kwargs():
    region_str = ('# Region file format: DS9 astropy/regions\n'
                  'image\n'
                  'circle(101.0,101.0,3.0) # color=red\n'
                  'circle(101.0,301.0,3.0)\n'
                  'point(101.0,101.0) # color=red\n'
                  'point(101.0,301.0)\n'
                  'text(101.0,101.0) # text={Text} color=red\n'
                  'text(101.0,101.0) # text={Text}')
    regions = Regions.parse(region_str, format='ds9')

    # Patch
    assert regions[0].visual['default_style'] == 'ds9'
    assert regions[0].as_artist().get_edgecolor() == (1, 0, 0, 1)
    assert regions[1].as_artist().get_edgecolor() == (0, 1, 0, 1)
    for i in (0, 1):
        regions[i].visual['default_style'] = None
    assert regions[0].as_artist().get_edgecolor() == (1, 0, 0, 1)
    assert regions[1].as_artist().get_edgecolor() == (0, 0, 0, 1)

    # Line2D
    assert regions[2].visual['default_style'] == 'ds9'
    assert regions[2].as_artist().get_markeredgecolor() == 'red'
    assert regions[3].as_artist().get_markeredgecolor() == '#00ff00'
    for i in (2, 3):
        regions[i].visual['default_style'] = None
    assert regions[2].as_artist().get_markeredgecolor() == 'red'
    assert regions[3].as_artist().get_markeredgecolor() == 'C0'

    # Text
    assert regions[4].visual['default_style'] == 'ds9'
    assert regions[4].as_artist().get_color() == 'red'
    assert regions[5].as_artist().get_color() == '#00ff00'
    for i in (4, 5):
        assert regions[i].as_artist().get_va() == 'center'
        assert regions[i].as_artist().get_ha() == 'center'
    for i in (4, 5):
        regions[i].visual['default_style'] = None
        assert regions[i].as_artist().get_va() == 'baseline'
        assert regions[i].as_artist().get_ha() == 'left'
    assert regions[4].as_artist().get_color() == 'red'
    assert regions[5].as_artist().get_color() == 'black'


def test_annulus():
    region1_str = 'image\nannulus(1, 2, 3, 4)'
    regions1 = Regions.parse(region1_str, format='ds9')
    assert len(regions1) == 1
    assert_equal(regions1[0].center.xy, (0, 1))
    assert_equal(regions1[0].inner_radius, 3.0)
    assert_equal(regions1[0].outer_radius, 4.0)

    # multiple annuli
    region2_str = 'image\nannulus(1, 2, 3, 4, 5, 6)'
    regions2 = Regions.parse(region2_str, format='ds9')
    assert len(regions2) == 3
    assert_equal(regions2[0].center.xy, (0, 1))
    assert_equal(regions2[0].inner_radius, 3.0)
    assert_equal(regions2[0].outer_radius, 4.0)
    assert_equal(regions2[1].center.xy, (0, 1))
    assert_equal(regions2[1].inner_radius, 4.0)
    assert_equal(regions2[1].outer_radius, 5.0)
    assert_equal(regions2[2].center.xy, (0, 1))
    assert_equal(regions2[2].inner_radius, 5.0)
    assert_equal(regions2[2].outer_radius, 6.0)


def test_ellipse():
    region1_str = 'image\nellipse (1, 2, 3, 4, 5)'
    region1 = Regions.parse(region1_str, format='ds9')[0]
    assert_equal(region1.center.xy, (0, 1))
    assert_equal(region1.width, 6.0)
    assert_equal(region1.height, 8.0)
    assert_equal(region1.angle, 5.0 * u.deg)

    region2_str = 'icrs\nellipse (1, 2, 3, 4, 5)'
    region2 = Regions.parse(region2_str, format='ds9')[0]
    assert_equal(region2.center.ra, 1. * u.deg)
    assert_equal(region2.center.dec, 2. * u.deg)
    assert_equal(region2.width, 6.0 * u.deg)
    assert_equal(region2.height, 8.0 * u.deg)
    assert_equal(region2.angle, 5.0 * u.deg)

    # elliptical annulus
    region3_str = 'image\nellipse (1, 2, 3, 4, 5, 6, 7)'
    region3 = Regions.parse(region3_str, format='ds9')[0]
    assert_equal(region3.center.xy, (0, 1))
    assert_equal(region3.inner_width, 6.0)
    assert_equal(region3.inner_height, 8.0)
    assert_equal(region3.outer_width, 10.0)
    assert_equal(region3.outer_height, 12.0)
    assert_equal(region3.angle, 7.0 * u.deg)

    # multiple elliptical annuli
    region4_str = 'image\nellipse (1, 2, 3, 4, 5, 6, 7, 8, 9)'
    region4 = Regions.parse(region4_str, format='ds9')
    assert len(region4) == 2
    assert_equal(region4[0].center.xy, (0, 1))
    assert_equal(region4[0].inner_width, 6.0)
    assert_equal(region4[0].inner_height, 8.0)
    assert_equal(region4[0].outer_width, 10.0)
    assert_equal(region4[0].outer_height, 12.0)
    assert_equal(region4[0].angle, 9.0 * u.deg)
    assert_equal(region4[1].center.xy, (0, 1))
    assert_equal(region4[1].inner_width, 10.0)
    assert_equal(region4[1].inner_height, 12.0)
    assert_equal(region4[1].outer_width, 14.0)
    assert_equal(region4[1].outer_height, 16.0)
    assert_equal(region4[1].angle, 9.0 * u.deg)


def test_box():
    region1_str = 'image\nbox (1, 2, 3, 4, 5)'
    region1 = Regions.parse(region1_str, format='ds9')[0]
    assert_equal(region1.center.xy, (0, 1))
    assert_equal(region1.width, 3.0)
    assert_equal(region1.height, 4.0)
    assert_equal(region1.angle, 5.0 * u.deg)

    region2_str = 'icrs\nbox (1, 2, 3, 4, 5)'
    region2 = Regions.parse(region2_str, format='ds9')[0]
    assert_equal(region2.center.ra, 1. * u.deg)
    assert_equal(region2.center.dec, 2. * u.deg)
    assert_equal(region2.width, 3.0 * u.deg)
    assert_equal(region2.height, 4.0 * u.deg)
    assert_equal(region2.angle, 5.0 * u.deg)

    # box (rectangle) annulus
    region3_str = 'image\nbox (1, 2, 3, 4, 5, 6, 7)'
    region3 = Regions.parse(region3_str, format='ds9')[0]
    assert_equal(region3.center.xy, (0, 1))
    assert_equal(region3.inner_width, 3.0)
    assert_equal(region3.inner_height, 4.0)
    assert_equal(region3.outer_width, 5.0)
    assert_equal(region3.outer_height, 6.0)
    assert_equal(region3.angle, 7.0 * u.deg)

    # multiple box (rectangle) annuli
    region4_str = 'image\nbox (1, 2, 3, 4, 5, 6, 7, 8, 9)'
    region4 = Regions.parse(region4_str, format='ds9')
    assert len(region4) == 2
    assert_equal(region4[0].center.xy, (0, 1))
    assert_equal(region4[0].inner_width, 3.0)
    assert_equal(region4[0].inner_height, 4.0)
    assert_equal(region4[0].outer_width, 5.0)
    assert_equal(region4[0].outer_height, 6.0)
    assert_equal(region4[0].angle, 9.0 * u.deg)
    assert_equal(region4[1].center.xy, (0, 1))
    assert_equal(region4[1].inner_width, 5.0)
    assert_equal(region4[1].inner_height, 6.0)
    assert_equal(region4[1].outer_width, 7.0)
    assert_equal(region4[1].outer_height, 8.0)
    assert_equal(region4[1].angle, 9.0 * u.deg)


def test_invalid_metadata():
    # test that invalid ds9 metadata raises warnings
    regstr = ('# Region file format: DS9 version 4.1\n'
              'global color=green dashlist=8 3 width=1 '
              'font="helvetica 10 normal roman" select=1 highlite=1 dash=1 '
              'fixed=0 edit=1 move=1 rotate=1 delete=1 include=1 source=1.1 '
              'background=0 text={hello world} fill=1.8 point=diamond 11.2 '
              'textrotate=1 textangle=34 tag={Tag 1} tag={Tag 2} line=3 1\n'
              'image; circle(100.4,47.2,10.9) # color=red\n'
              'image; circle(202.4,47.2,10.9) # color=blue point=junk 10\n'
              'image; circle(302.4,147.2,10.9) # color=cyan point=invalid\n')

    with pytest.warns(AstropyUserWarning) as record:
        regs = Regions.parse(regstr, format='ds9')
        assert len(record) == 12
        assert 'source=1.1' in record[0].message.args[0]
        assert 'fill=1.8' in record[1].message.args[0]
        assert 'point=diamond 11.2' in record[2].message.args[0]
        assert 'line=3 1' in record[3].message.args[0]
        assert 'point=junk 10' in record[6].message.args[0]
        assert 'point=invalid' in record[10].message.args[0]

    assert len(regs) == 3


def test_unsupported_metadata():
    regstr = ('image; point(335.5,415.6) # point=diamond 11 color=yellow '
              'width=2 text={example} dash=1 tag={Tag 1} tag={Tag 2}')
    with pytest.warns(UserWarning) as record:
        Regions.parse(regstr, format='ds9')
        assert len(record) == 1
        assert 'dashed lines are unsupported' in record[0].message.args[0]


def test_text_metadata():
    # Regression test for issue #233: make sure that text metadata is
    # parsed and stored appropriately.
    region_str = 'image\ncircle(1.5, 2, 2) # text={this_is_text}'
    regions = Regions.parse(region_str, format='ds9')
    assert len(regions) == 1
    assert regions[0].meta['text'] == 'this_is_text'

    rstr = 'image; circle(503.6,490.6,31.1) # color=blue text={A, {B}, C}'
    reg = Regions.parse(rstr, format='ds9')[0]
    assert reg.meta['text'] == 'A, {B'

    rstr = 'image; circle(503.6,490.6,31.1) # color=blue text={A, {{B}}, C}'
    reg = Regions.parse(rstr, format='ds9')[0]
    assert reg.meta['text'] == 'A, {{B'

    rstr = "image; circle(503.6,490.6,31.1) # color=blue text='A, 'B', C'"
    reg = Regions.parse(rstr, format='ds9')[0]
    assert reg.meta['text'] == 'A, '

    rstr = 'image; circle(503.6,490.6,31.1) # color=blue text="A, "B", C"'
    reg = Regions.parse(rstr, format='ds9')[0]
    assert reg.meta['text'] == 'A, '


def test_mixed_coord():
    with pytest.raises(ValueError):
        CirclePixelRegion(PixCoord(10, 20), Angle(1, 'arcsec'))


def test_unsupported_marker():
    """
    Test that warning is issued when serializing a valid matplotlib
    marker, but unsupported by DS9.
    """
    region = PointPixelRegion(PixCoord(2, 2), visual=RegionVisual(marker='Z'))
    with pytest.warns(AstropyUserWarning):
        region.serialize(format='ds9')

    region = PointPixelRegion(PixCoord(2, 2),
                              visual=RegionVisual(marker='$f_{init}$'))
    with pytest.warns(AstropyUserWarning):
        region.serialize(format='ds9')


def test_serialize_regularpolygon():
    region = RegularPolygonPixelRegion(PixCoord(10, 10), 4, 20)
    poly_region = region.to_polygon()
    result1 = region.serialize(format='ds9')
    result2 = poly_region.serialize(format='ds9')
    assert result1 == result2


def test_serialize_empty_list():
    regions = Regions([])
    assert regions.serialize(format='ds9') == ''
