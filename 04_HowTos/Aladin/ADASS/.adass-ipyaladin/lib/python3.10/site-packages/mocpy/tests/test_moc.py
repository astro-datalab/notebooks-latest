import copy
import math
import re

import astropy.units as u
import cdshealpix
import numpy as np
import pytest
import regions
from astropy.coordinates import Angle, Latitude, Longitude, SkyCoord, angular_separation
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.table import QTable

from ..moc import MOC, WCS

rng = np.random.default_rng()


@pytest.fixture
def isets():
    a = MOC.from_depth29_ranges(
        29,
        np.array(
            [
                [49, 73],
                [53, 54],
                [33, 63],
                [65, 80],
                [51, 80],
                [100, 126],
                [38, 68],
                [61, 72],
                [74, 102],
                [27, 43],
            ],
            dtype=np.uint64,
        ),
    )
    b = MOC.from_depth29_ranges(
        29,
        np.array(
            [
                [17, 26],
                [17, 41],
                [12, 31],
                [32, 61],
                [68, 90],
                [77, 105],
                [18, 27],
                [12, 35],
                [9, 37],
                [87, 97],
            ],
            dtype=np.uint64,
        ),
    )
    return {"a": a, "b": b}


def test_interval_set_consistency(isets):
    assert isets["a"] == MOC.from_depth29_ranges(
        29,
        np.array([[27, 126]], dtype=np.uint64),
    )
    assert isets["b"] == MOC.from_depth29_ranges(
        29,
        np.array([[9, 61], [68, 105]], dtype=np.uint64),
    )


def test_uniq_hpx():
    moc = MOC.from_depth29_ranges(29, np.array([[0, 1]], dtype=np.uint64))
    uniq_hpx = np.array([4 * 4**29], dtype=np.uint64)
    assert moc.uniq_hpx == uniq_hpx

    moc = MOC.from_depth29_ranges(29, np.array([[7, 76]], dtype=np.uint64))
    uniq_hpx = np.array(
        [
            1 + 4 * 4**27,
            2 + 4 * 4**27,
            3 + 4 * 4**27,
            2 + 4 * 4**28,
            3 + 4 * 4**28,
            16 + 4 * 4**28,
            17 + 4 * 4**28,
            18 + 4 * 4**28,
            7 + 4 * 4**29,
        ],
        dtype=np.uint64,
    )
    assert (np.sort(moc.uniq_hpx) == uniq_hpx).all()


def test_to_depth29_ranges(isets):
    l = isets["a"].to_depth29_ranges  # noqa: E741
    r = np.array([[27, 126]], dtype=np.uint64)
    assert np.array_equal(l, r)
    l = isets["b"].to_depth29_ranges  # noqa: E741
    r = np.array([[9, 61], [68, 105]], dtype=np.uint64)
    assert np.array_equal(l, r)


def test_n_cells():
    assert MOC.n_cells(0) == 12
    with pytest.raises(
        ValueError,
        match=f"The depth should be comprised between 0 and {MOC.MAX_ORDER}*",
    ):
        MOC.n_cells(-2)
    with pytest.raises(
        ValueError,
        match=f"The depth should be comprised between 0 and {MOC.MAX_ORDER}*",
    ):
        MOC.n_cells(MOC.MAX_ORDER + 1)
    assert MOC.n_cells(6) == 4 * MOC.n_cells(5)


def test_interval_set_union(isets):
    assert isets["a"].union(isets["b"]) == MOC.from_depth29_ranges(
        29,
        np.array([[9, 126]], dtype=np.uint64),
    )
    assert isets["a"].union(MOC.new_empty(29)) == MOC.from_depth29_ranges(
        29,
        np.array([[27, 126]], dtype=np.uint64),
    )
    assert MOC.new_empty(29).union(isets["a"]) == MOC.from_depth29_ranges(
        29,
        np.array([[27, 126]], dtype=np.uint64),
    )


def test_interval_set_intersection(isets):
    assert isets["a"].intersection(isets["b"]) == MOC.from_depth29_ranges(
        29,
        np.array([[27, 61], [68, 105]], dtype=np.uint64),
    )
    assert isets["a"].intersection(MOC.new_empty(29)) == MOC.new_empty(29)
    assert MOC.new_empty(29).intersection(isets["a"]) == MOC.new_empty(29)


def test_interval_set_difference(isets):
    assert isets["a"].difference(isets["b"]) == MOC.from_depth29_ranges(
        29,
        np.array([[61, 68], [105, 126]], dtype=np.uint64),
    )
    assert isets["b"].difference(isets["a"]) == MOC.from_depth29_ranges(
        29,
        np.array([[9, 27]], dtype=np.uint64),
    )
    assert MOC.new_empty(29).difference(isets["a"]) == MOC.new_empty(29)
    assert isets["a"].difference(MOC.new_empty(29)) == isets["a"]


def test_interval_set_min(isets):
    assert isets["a"].min_index == np.uint64(27)
    assert isets["b"].min_index == np.uint64(9)
    assert isets["a"].union(isets["b"]).min_index == np.uint64(9)


def test_interval_set_max(isets):
    assert isets["a"].max_index == np.uint64(126)
    assert isets["b"].max_index == np.uint64(105)
    assert isets["a"].union(isets["b"]).max_index == np.uint64(126)


def test_interval_min_depth():
    big_cells = np.array([[0, 4**29]], dtype=np.uint64)
    itv_result = MOC.from_depth29_ranges(29, big_cells)

    small_cells = np.array(
        [
            [0, 4**28],
            [4**28, 2 * 4**28],
            [2 * 4**28, 3 * 4**28],
            [3 * 4**28, 4**29],
        ],
        dtype=np.uint64,
    )
    itv_small_cells = MOC.from_depth29_ranges(29, small_cells)
    assert itv_result == itv_small_cells


def test_complement():
    assert MOC.from_depth29_ranges(
        max_depth=29,
        ranges=None,
    ).complement() == MOC.from_depth29_ranges(
        max_depth=29,
        ranges=np.array([[0, 12 * 4**29]], dtype=np.uint64),
    )
    assert MOC.new_empty(
        max_depth=29,
    ).complement().complement() == MOC.from_depth29_ranges(max_depth=29, ranges=None)
    assert MOC.from_depth29_ranges(
        29,
        np.array([[1, 2], [6, 8], [5, 6]], dtype=np.uint64),
    ).complement() == MOC.from_depth29_ranges(
        29,
        np.array([[0, 1], [2, 5], [8, 12 * 4**29]], dtype=np.uint64),
    )


# --- TESTING MOC creation ---#


def test_new_empty_serialization():
    # regression test for https://github.com/cds-astro/mocpy/issues/146
    empty = MOC.new_empty(max_depth=0)
    assert empty.serialize("json") == {"0": []}


def get_random_skycoords(size):
    return SkyCoord(
        ra=rng.random(size) * 360,
        dec=rng.random(size) * 180 - 90,
        unit="deg",
    )


@pytest.fixture
def skycoords_gen_f():
    def gen_f(size):
        return SkyCoord(
            ra=rng.random(size) * 360,
            dec=rng.random(size) * 180 - 90,
            unit="deg",
        )

    return gen_f


@pytest.fixture
def lonlat_gen_f():
    def gen_f(size):
        return ((rng.random(size) * 360) * u.deg, (rng.random(size) * 180 - 90) * u.deg)

    return gen_f


def test_moc_from_skycoords(skycoords_gen_f):
    skycoords = skycoords_gen_f(10)
    moc = MOC.from_skycoords(skycoords, max_norder=7)
    assert all(moc.contains_skycoords(skycoords))


def test_moc_from_lonlat(lonlat_gen_f):
    lon, lat = lonlat_gen_f(10)
    moc_lon_lat = MOC.from_lonlat(lon=lon, lat=lat, max_norder=6)
    moc_skycoo = MOC.from_skycoords(SkyCoord(lon, lat), max_norder=6)
    assert moc_lon_lat == moc_skycoo


def test_from_healpix_cells():
    ipix = np.array([40, 87, 65])
    depth = np.array([3, 3, 3])
    moc = MOC.from_healpix_cells(max_depth=3, ipix=ipix, depth=depth)
    # the MOC is correct
    assert moc == MOC.from_str("3/40 87 65")
    # we can also give depth as an integer if the cells have only one order
    assert moc == MOC.from_healpix_cells(ipix=ipix, depth=3, max_depth=3)
    # negative indices are ignored
    with pytest.warns(
        UserWarning,
        match="The list of indices contain negative values*",
    ):
        moc = MOC.from_healpix_cells(ipix=[40, -1, 65], depth=depth, max_depth=3)
    assert moc == MOC.from_str("3/40 65")
    # also allow order zero (regression for issue #157)
    assert MOC.from_healpix_cells(np.array([0]), 0, 0) == MOC.from_str("0/0")


def test_from_polygons():
    list_vertices = [
        SkyCoord([-4, 4, 4, -4], [4, 4, -4, -4], unit="deg"),
        SkyCoord([0, 6, 0, -6], [6, 0, -6, 0], unit="deg"),
    ]
    list_mocs = MOC.from_polygons(list_vertices)
    assert all(
        moc.barycenter().separation(SkyCoord(0, 0, unit="deg")) < 0.1 * u.arcmin
        for moc in list_mocs
    )
    list_mocs_no_skycoord = MOC.from_polygons(
        [[356, 4, 4, 356], [4, 4, -4, -4], [0, 6, 0, 354], [6, 0, -6, 0]],
    )
    assert list_mocs == list_mocs_no_skycoord


def test_from_vizier():
    # deprecated nside should still work (nside=8 means order=3)
    with pytest.warns(
        DeprecationWarning,
        match="'nside' is deprecated in favor of 'max_depth'.*",
    ):
        moc = MOC.from_vizier_table("I/355", nside=8)
    assert moc.max_order == 3
    # gaia is the whole sky at order 3
    assert moc.sky_fraction == 1
    with pytest.warns(
        DeprecationWarning,
        match="'nside' is deprecated in favor of 'max_depth'.*",
    ), pytest.raises(ValueError, match="Bad value for nside.*"):
        moc = MOC.from_vizier_table("I/355", nside=1)
    moc = MOC.from_vizier_table("I/355")
    # default order is 10 for catalogs
    assert moc.max_order == 10
    # non valid table or catalog
    with pytest.raises(ValueError, match="No catalog/table was found for 'abc'*"):
        moc = MOC.from_vizier_table("abc")


def test_from_ivorn():
    with pytest.warns(
        DeprecationWarning,
        match="'nside' is deprecated in favor of 'max_depth'.*",
    ):
        moc = MOC.from_ivorn("ivo://CDS/J/A+AS/133/387/table5", nside=8)
    assert moc.max_order == 3

    with pytest.warns(UserWarning, match="This MOC is empty.*"):
        MOC.from_ivorn("abc")


def test_moc_from_fits():
    fits_path = "resources/P-GALEXGR6-AIS-FUV.fits"
    MOC.load(fits_path, "fits")


def test_moc_from_fits_url():
    url = "http://skies.esac.esa.int/Spitzer/IRAC1_bright_ISM/Moc.fits"
    MOC.from_fits(url)


def test_moc_consistent_with_aladin():
    truth = MOC.load("resources/CDS-I-125A-catalog_MOC.fits", "fits")
    table = parse_single_table("resources/I_125A_catalog.vot").to_table()

    moc = MOC.from_lonlat(
        table["_RAJ2000"].T * u.deg,
        table["_DEJ2000"].T * u.deg,
        max_norder=8,
    )

    assert moc == truth


def test_moc_from_fits_images():
    image_path = "resources/image_with_mask.fits.gz"
    MOC.from_fits_images([image_path], max_norder=15, hdu_index=-1)


def test_from_fits_images_2():
    MOC.from_fits_images(["resources/u_gal.fits"], max_norder=10)


def test_from_fits_image_without_cdelt():
    MOC.from_fits_images(["resources/horsehead.fits"], max_norder=5, hdu_index=-1)


@pytest.fixture
def moc_from_fits_image():
    image_path = "resources/image_with_mask.fits.gz"

    with fits.open(image_path) as hdulist:
        return MOC.from_fits_image(hdu=hdulist[0], max_norder=7, mask=hdulist[0].data)


@pytest.fixture
def moc_from_json():
    return MOC.from_json({"8": [45, 78], "4": [42, 57]})


def test_moc_from_fits_image(moc_from_fits_image):
    assert isinstance(moc_from_fits_image, MOC)


def test_moc_serialize_and_from_json(moc_from_json):
    ipix_d = moc_from_json.serialize(format="json")
    moc2 = MOC.from_json(ipix_d)
    assert moc_from_json == moc2


@pytest.mark.parametrize(
    ("expected", "moc_str"),
    [
        (
            MOC.from_json(
                {
                    "5": [8, 9, 10, 42, 43, 44, 45, 54, 46],
                    "6": [4500],
                    "7": [],
                    "8": [45],
                },
            ),
            "5/8-10 42-46 54\n\r 6/4500 8/45",
        ),
        (MOC.from_json({}), "0/"),
        (MOC.from_json({"29": [101]}), "29/101"),
        (MOC.from_json({"0": [1, 0, 9]}), "0/0-1 9"),
        (MOC.from_json({"0": [2, 9]}), "0/2 9"),
        (MOC.from_json({"0": [2], "8": [8, 9, 10], "11": []}), "0/2\r \n 8/8-10\n 11/"),
    ],
    # (MOC.from_json({"0": [2, 9], "1": [9]}), "0/2 9"),
)
def test_from_str(expected, moc_str):
    assert MOC.from_str(moc_str) == expected


@pytest.mark.parametrize(
    ("expected", "moc_str"),
    [
        (
            MOC.from_json(
                {
                    "5": [8, 9, 10, 42, 43, 44, 45, 54, 46],
                    "6": [4500],
                    "7": [],
                    "8": [45],
                },
            ),
            "5/8-10 42-46 54\n\r 6/4500 8/45",
        ),
        (MOC.from_json({}), "0/"),
        (MOC.from_json({"29": [101]}), "29/101"),
        (MOC.from_json({"0": [1, 0, 9]}), "0/0-1 9"),
        (MOC.from_json({"0": [2, 9]}), "0/2 9"),
        (MOC.from_json({"0": [2], "8": [8, 9, 10], "11": []}), "0/2\r \n 8/8-10\n 11/"),
    ],
    # (MOC.from_json({"0": [2, 9], "1": [9]}), "0/2 9"),
)
def test_from_string(expected, moc_str):
    assert MOC.from_string(moc_str, "ascii") == expected


def test_moc_full_skyfraction():
    moc = MOC.from_json({"0": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]})
    assert moc.sky_fraction == 1.0


def test_moc_skyfraction():
    moc = MOC.from_json({"0": [0, 1, 2, 3, 4, 5]})
    assert moc.sky_fraction == 0.5


def test_sky_fraction_on_empty_coverage():
    moc = MOC.new_empty(max_depth=29)
    assert moc.sky_fraction == 0


# --- TESTING MOC serialization ---#
def test_moc_serialize_to_fits(moc_from_fits_image):
    hdulist = moc_from_fits_image.serialize(format="fits")
    assert isinstance(hdulist, fits.hdu.hdulist.HDUList)


def test_moc_serialize_to_json(moc_from_fits_image):
    moc_json = moc_from_fits_image.serialize(format="json")
    assert isinstance(moc_json, dict)


@pytest.mark.parametrize(
    ("moc", "expected"),
    [
        (
            MOC.from_json(
                {
                    "5": [8, 9, 10, 42, 43, 44, 45, 54, 46],
                    "6": [4500],
                    "7": [],
                    "8": [45],
                },
            ),
            "5/8-10 42-46 54 6/4500 8/45",
        ),
        (MOC.from_json({}), "0/"),
        (MOC.from_json({"29": [101]}), "29/101"),
        (MOC.from_json({"0": [1, 0, 9]}), "0/0-1 9"),
        (MOC.from_json({"0": [2, 9]}), "0/2 9"),
    ],
    #  (MOC.from_json({"0": [2, 9], "1": [9]}), "0/2 9"),
)
def test_serialize_to_str(moc, expected):
    assert moc.serialize(format="str") == expected


# --- TESTING MOC plot functions ---#
def test_mpl_fill():
    fits_path = "resources/P-GALEXGR6-AIS-FUV.fits"
    moc = MOC.load(fits_path, "fits")

    import matplotlib.pyplot as plt

    fig = plt.figure(111, figsize=(10, 10))
    with WCS(
        fig,
        fov=50 * u.deg,
        center=SkyCoord(0, 20, unit="deg", frame="icrs"),
        coordsys="icrs",
        rotation=Angle(0, u.degree),
        projection="AIT",
    ) as wcs:
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc.fill(ax=ax, wcs=wcs, alpha=0.5, color="r")


def test_mpl_border():
    fits_path = "resources/P-GALEXGR6-AIS-FUV.fits"
    moc = MOC.load(fits_path, "fits")

    import matplotlib.pyplot as plt

    fig = plt.figure(111, figsize=(10, 10))
    with WCS(
        fig,
        fov=50 * u.deg,
        center=SkyCoord(0, 20, unit="deg", frame="icrs"),
        coordsys="icrs",
        rotation=Angle(0, u.degree),
        projection="AIT",
    ) as wcs:
        ax = fig.add_subplot(1, 1, 1, projection=wcs)
        moc.border(ax=ax, wcs=wcs, color="g")


# --- TESTING MOC features ---#
@pytest.mark.parametrize("order", [4, 5, 6, 15, 20, 28])
def test_moc_contains(order):
    # defines 20 random healpix cells of the required order
    size = 20
    healpix_arr = rng.integers(0, 12 * 4**order, size, dtype="uint64", endpoint=False)
    # defines a moc containing the 20 points
    moc = MOC.from_json(json_moc={str(order): np.unique(healpix_arr).tolist()})
    # the complementary should not contain them
    moc_complement = moc.complement()
    # coordinates of the 20 random points
    lon, lat = cdshealpix.healpix_to_lonlat(healpix_arr, order)
    # tests
    should_be_inside_arr = moc.contains_lonlat(lon=lon, lat=lat)
    assert should_be_inside_arr.all()
    should_be_outside_arr = moc_complement.contains_lonlat(lon=lon, lat=lat)
    assert not should_be_outside_arr.any()
    # test keep_inside field
    should_be_outside_arr = moc.contains_lonlat(lon=lon, lat=lat, keep_inside=False)
    assert not should_be_outside_arr.any()
    should_be_inside_arr = moc_complement.contains_lonlat(
        lon=lon,
        lat=lat,
        keep_inside=False,
    )
    assert should_be_inside_arr.all()
    # test only floats in arguments, this is a regression test from #108
    moc = MOC.from_string("2/4")
    assert moc.contains_lonlat(164.43 * u.deg, 45.54 * u.deg) == [False]


# test 2d-arrays as lon lat input
def test_moc_contains_2d_parameters():
    """Test that not only 1d arrays are accepted."""
    lon = Angle(np.array([[1, 2, 3], [2, 40, -5]]), unit=u.deg)
    lat = Angle(np.array([[20, 25, 10], [-60, 80, 0]]), unit=u.deg)
    lat2 = Angle(np.array([[20, 25, 10, 22], [-60, 80, 0, 10]]), unit=u.deg)
    moc = MOC.from_polygon(lon=lon, lat=lat, max_depth=12)
    should_be_inside = moc.contains_lonlat(lon=lon, lat=lat)
    complement_moc = MOC.from_polygon(lon=lon, lat=lat, max_depth=12, complement=True)
    assert complement_moc.sky_fraction > moc.sky_fraction
    assert should_be_inside.all()

    # test mismatched
    with pytest.raises(
        ValueError,
        match=re.escape(
            "'lon' and 'lat' should have the same shape but are of shapes (2, 3) and (2, 4)",
        ),
    ):
        moc.contains_lonlat(lon=lon, lat=lat2)


def test_degrade_to_order():
    hst_fits_path = "resources/hst.fits"
    hst_moc = MOC.load(hst_fits_path, "fits")

    max_depth = hst_moc.max_order

    for order in reversed(range(max_depth)):
        hst_moc = hst_moc.degrade_to_order(order)
        assert hst_moc.sky_fraction <= 1.0


def test_from_ring():
    MOC.from_ring(
        lon=0 * u.deg,
        lat=0 * u.deg,
        internal_radius=Angle(5, u.deg),
        external_radius=Angle(10, u.deg),
        max_depth=10,
    )


def test_from_zone():
    moc = MOC.from_zone(SkyCoord([[-50, -50], [50, 50]], unit="deg"), max_depth=5)
    # test the diagonal
    for coordinate in range(-50, 40, 10):  ## (50,50) not included
        assert moc.contains_skycoords(SkyCoord(coordinate, coordinate, unit="deg"))
    # regression for #180
    MOC.from_zone(SkyCoord([(180, 30), (360, 50)], unit="deg"), max_depth=3)


def test_from_box():
    a = Angle("10d")
    b = Angle("2d")
    moc = MOC.from_box(
        lon=Longitude("0d"),
        lat=Latitude("0d"),
        a=a,
        b=b,
        angle=Angle("30deg"),
        max_depth=10,
    )
    area = moc.sky_fraction * 4 * math.pi * u.steradian
    # the moc covers a slightly bigger area than the region defined by the
    # parameters
    assert area.to(u.deg**2).value > 80
    assert area.to(u.deg**2).value < 90
    # test from_boxes
    list_boxes_same = MOC.from_boxes(
        lon=[0, 0] * u.deg,
        lat=[0, 0] * u.deg,
        a=a,
        b=b,
        angle=30 * u.deg,
        max_depth=10,
    )
    assert len(list_boxes_same) == 2
    assert list_boxes_same[0] == moc
    # union strategies
    union_boxes_same = MOC.from_boxes(
        lon=[0, 0] * u.deg,
        lat=[0, 0] * u.deg,
        a=a,
        b=b,
        angle=30 * u.deg,
        max_depth=10,
        union_strategy="large_boxes",
    )
    assert union_boxes_same == MOC.from_box(
        lon=0 * u.deg, lat=0 * u.deg, a=a, b=b, angle=30 * u.deg, max_depth=10
    )

    # not same boxes
    a = [Angle("10d"), Angle("20d")]
    b = [Angle("2d"), Angle("4d")]
    list_boxes_different = MOC.from_boxes(
        lon=[0, 0] * u.deg,
        lat=[0, 0] * u.deg,
        a=a,
        b=b,
        angle=[30, 45] * u.deg,
        max_depth=10,
    )
    assert len(list_boxes_different) == 2
    assert list_boxes_different[0] == moc
    # mixed iterables and scalars raise an error
    with pytest.raises(ValueError, match="'a', 'b' and 'angle' should*"):
        MOC.from_boxes(
            lon=[0, 0] * u.deg,
            lat=[0, 0] * u.deg,
            a=a,
            b=b,
            angle=30 * u.deg,
            max_depth=10,
        )
    # union strategy possible choices
    with pytest.raises(
        ValueError,
        match="'union_strategy' can only be None, 'large_boxes', or 'small_boxes'.",
    ):
        MOC.from_boxes(
            lon=[0, 0] * u.deg,
            lat=[0, 0] * u.deg,
            a=a,
            b=b,
            angle=[30, 30] * u.deg,
            max_depth=10,
            union_strategy="large_cones",  # voluntary confusion between cones and boxes
        )


def test_from_astropy_regions():
    center = SkyCoord(42, 43, unit="deg", frame="fk5")
    # circle
    circle = regions.CircleSkyRegion(center, radius=3 * u.deg)
    moc = MOC.from_astropy_regions(circle, max_depth=10)
    assert round(moc.barycenter().ra.value) == 42
    assert round(moc.barycenter().dec.value) == 43
    assert round(moc.largest_distance_from_coo_to_vertices(center).to("deg").value) == 3
    # ring
    ring = regions.CircleAnnulusSkyRegion(center, 3 * u.deg, 4 * u.deg)
    moc = MOC.from_astropy_regions(ring, max_depth=9)
    assert not moc.contains_skycoords(center)
    # ellipse
    ellipse = regions.EllipseSkyRegion(center, 3 * u.deg, 6 * u.deg, 10 * u.deg)
    moc = MOC.from_astropy_regions(ellipse, max_depth=9)
    assert moc.contains_skycoords(center)
    assert round(moc.largest_distance_from_coo_to_vertices(center).to("deg").value) == 3
    # invert axes in the ellipse
    ellipse_swapped = regions.EllipseSkyRegion(
        center,
        6 * u.deg,
        3 * u.deg,
        (-90 + 10) * u.deg,
    )
    moc_swapped = MOC.from_astropy_regions(ellipse_swapped, max_depth=9)
    assert moc_swapped == moc
    # rectangle
    box = regions.RectangleSkyRegion(center, 12 * u.deg, 6 * u.deg, 10 * u.deg)
    moc = MOC.from_astropy_regions(box, max_depth=8)
    assert all(
        moc.contains_skycoords(SkyCoord([42, 45], [44, 44], unit="deg", frame="icrs")),
    )
    # polygons
    vertices = SkyCoord([1, 2, 2], [1, 1, 2], unit="deg", frame="fk5")
    polygon = regions.PolygonSkyRegion(vertices)
    moc_polygon = MOC.from_astropy_regions(polygon, max_depth=10)
    assert all(moc_polygon.contains_skycoords(vertices))
    # points
    point = SkyCoord("+23h20m48.3s +61d12m06s")
    region_point = regions.PointSkyRegion(point)
    moc_point = MOC.from_astropy_regions(region_point, max_depth=10)
    assert moc_point.max_order == 10
    assert moc_point.contains_skycoords(point)
    # multi regions
    multi_region = regions.Regions([region_point, polygon])
    moc = MOC.from_astropy_regions(multi_region, max_depth=10)
    assert moc == moc_polygon + moc_point


# TODO: IMPROVE THE ALGO
"""
def test_boundaries():
    fits_path = 'resources/P-GALEXGR6-AIS-FUV.fits'
    moc = MOC.load(fits_path, 'fits')
    moc = moc.degrade_to_order(6)
    boundaries_l = moc.get_boundaries()
"""


def test_from_elliptical_cone():
    MOC.from_elliptical_cone(
        lon=0 * u.deg,
        lat=0 * u.deg,
        a=Angle(10, u.deg),
        b=Angle(5, u.deg),
        pa=Angle(0, u.deg),
        max_depth=10,
    )


def test_from_cone():
    with pytest.raises(ValueError, match="'MOC.from_cone' only works with one cone.*"):
        MOC.from_cone([2, 4] * u.deg, [5, 6] * u.deg, radius=2 * u.arcmin, max_depth=2)


def test_from_cones():
    # same radius
    radius = 2 * u.arcmin
    lon = [2, 4] * u.deg
    lat = [5, 6] * u.deg
    cones = MOC.from_cones(lon, lat, radius=radius, max_depth=14)
    for cone, lon_unique, lat_unique in zip(cones, lon, lat):
        barycenter = cone.barycenter()
        assert angular_separation(
            lon_unique,
            lat_unique,
            barycenter.ra,
            barycenter.dec,
        ) < Angle(1 * u.arcsec)
    moc = MOC.from_cones(
        lon,
        lat,
        radius=radius,
        max_depth=14,
        union_strategy="small_cones",
    )
    assert isinstance(moc, MOC)
    moc2 = MOC.from_cones(
        lon,
        lat,
        radius=radius,
        max_depth=14,
        union_strategy="large_cones",
    )
    assert moc == moc2
    # different radii
    radii = [5, 6] * u.arcmin
    cones = MOC.from_cones(lon, lat, radius=radii, max_depth=14)
    for cone, radius in zip(cones, radii):
        # we check their area (Pi simplifies)
        assert (
            cone.sky_fraction
            > (((radius) ** 2).to(u.steradian) / 4 * u.steradian).value
        )
    moc = MOC.from_cones(
        lon,
        lat,
        radius=radii,
        max_depth=14,
        union_strategy="small_cones",
    )
    assert isinstance(moc, MOC)
    # test error
    with pytest.raises(ValueError, match="'union_strategy'*"):
        MOC.from_cones(lon, lat, radius=radii, max_depth=14, union_strategy="big_cones")


@pytest.fixture
def mocs():
    moc1 = {"1": [0]}
    moc1_increased = {"0": [0], "1": [17, 19, 22, 23, 35]}
    moc2 = {"1": [30]}
    moc2_increased = {"0": [7], "1": [8, 9, 25, 43, 41]}

    return {
        "moc1": MOC.from_json(moc1),
        "moc1_increased": MOC.from_json(moc1_increased),
        "moc2": MOC.from_json(moc2),
        "moc2_increased": MOC.from_json(moc2_increased),
    }


def test_add_neighbours(mocs):
    mocs["moc1"].add_neighbours()
    assert mocs["moc1"] == mocs["moc1_increased"]
    mocs["moc2"].add_neighbours()
    assert mocs["moc2"] == mocs["moc2_increased"]


def test_remove_neighbours(mocs):
    mocs["moc1_increased"].remove_neighbours()
    mocs["moc2_increased"].remove_neighbours()
    assert mocs["moc1_increased"] == mocs["moc1"]
    assert mocs["moc2_increased"] == mocs["moc2"]


def test_neighbours(mocs):
    moc1 = copy.deepcopy(mocs["moc1"])
    moc2 = copy.deepcopy(mocs["moc2"])
    moc1.add_neighbours().remove_neighbours()
    moc2.add_neighbours().remove_neighbours()
    assert moc1 == mocs["moc1"]
    assert moc2 == mocs["moc2"]


# --- TESTING MOC operations ---#
@pytest.fixture
def mocs_op():
    moc1 = MOC.from_json({"0": [0, 2, 3, 4, 5]})
    moc2 = MOC.from_json({"0": [0, 1, 7, 4, 3]})
    return {"first": moc1, "second": moc2}


def test_moc_union(mocs_op):
    assert mocs_op["first"].union(mocs_op["second"]) == MOC.from_json(
        {"0": [0, 1, 2, 3, 4, 5, 7]},
    )
    assert mocs_op["first"] + mocs_op["second"] == MOC.from_json(
        {"0": [0, 1, 2, 3, 4, 5, 7]},
    )
    assert mocs_op["first"] | mocs_op["second"] == MOC.from_json(
        {"0": [0, 1, 2, 3, 4, 5, 7]},
    )


def test_sum(mocs_op):
    assert sum([mocs_op["first"], mocs_op["second"]]) == MOC.from_json(
        {"0": [0, 1, 2, 3, 4, 5, 7]},
    )


def test_moc_intersection(mocs_op):
    assert mocs_op["first"].intersection(mocs_op["second"]) == MOC.from_json(
        {"0": [0, 3, 4]},
    )
    assert mocs_op["first"] & mocs_op["second"] == MOC.from_json({"0": [0, 3, 4]})


def test_moc_difference(mocs_op):
    assert mocs_op["first"].difference(mocs_op["second"]) == MOC.from_json(
        {"0": [2, 5]},
    )
    assert mocs_op["first"] - mocs_op["second"] == MOC.from_json({"0": [2, 5]})


def test_moc_complement_consistency():
    moc = MOC.load("resources/P-GALEXGR6-AIS-FUV.fits", "fits")
    assert moc.complement().complement() == moc


def test_from_fits_old():
    MOC.from_fits("resources/V_147_sdss12.moc.fits")


@pytest.mark.parametrize(
    ("input_MOC", "expected"),
    [
        (
            MOC.from_json({"0": [1, 3]}),
            MOC.from_json({"0": [0, 2, 4, 5, 6, 7, 8, 9, 10, 11]}),
        ),
    ],
)
def test_moc_complement(input_MOC, expected):
    assert input_MOC.complement() == expected
    assert ~input_MOC == expected


def test_spatial_res_to_order():
    order = np.arange(14)

    res = MOC.order_to_spatial_resolution(order)
    output = MOC.spatial_resolution_to_order(res)

    assert (order == output).all()


def test_from_valued_healpix_cells():
    uniq = np.array([8, 12])
    values = np.array([0.4, 0.6])
    moc_1 = MOC.from_valued_healpix_cells(
        uniq,
        values,
        cumul_from=0,
        cumul_to=1,
        max_depth=12,
    )
    assert moc_1 == MOC.from_string("0/4 8\n12/")
    moc_0p7 = MOC.from_valued_healpix_cells(
        uniq,
        values,
        cumul_from=0,
        cumul_to=0.7,
        max_depth=12,
    )
    assert moc_0p7 == MOC.from_string("0/8\n12/")

    # different sizes
    uniq = np.array([500])
    values = np.array([])
    with pytest.raises(
        ValueError,
        match="`uniq` and values do not have the same size.",
    ):
        MOC.from_valued_healpix_cells(uniq, values, 12)
    # failed comparison between cumul_from and cumul_to
    values = np.array([1.0])
    with pytest.raises(ValueError, match="`cumul_from` has to be < to `cumul_to`."):
        MOC.from_valued_healpix_cells(uniq, values, 12, cumul_from=0.8, cumul_to=-5.0)
    # negative uniq
    uniq = [-2, 270]
    values = [0.2, 0.2]
    with pytest.warns(
        UserWarning,
        match="The list of indices contain negative values*",
    ):
        assert MOC.from_valued_healpix_cells(
            uniq,
            values,
            max_depth=3,
        ) == MOC.from_string("3/14")


@pytest.mark.parametrize(
    ("cumul_from", "cumul_to"),
    [(-5.0, 1.0), (np.nan, np.inf), (np.nan, np.nan), (np.inf, np.nan), (-10.0, -5.0)],
)
def test_from_valued_healpix_cells_weird_values(cumul_from, cumul_to):
    uniq = np.array([500])
    values = np.array([-1.0])

    MOC.from_valued_healpix_cells(
        uniq,
        values,
        12,
        cumul_from=cumul_from,
        cumul_to=cumul_to,
    )


def test_from_valued_healpix_cells_bayestar():
    from astropy.io import fits

    fits_image_filename = "./resources/bayestar.multiorder.fits"

    with fits.open(fits_image_filename) as hdul:
        data = hdul[1].data

    uniq = data["UNIQ"]
    probdensity = data["PROBDENSITY"]

    import astropy.units as u

    orders = (np.log2(uniq // 4)) // 2
    area = 4 * np.pi / np.array([MOC.n_cells(int(order)) for order in orders]) * u.sr

    prob = probdensity * area

    cumul_to = np.linspace(0.01, 2.0, num=10)

    for b in cumul_to:
        MOC.from_valued_healpix_cells(uniq, prob, 12, cumul_from=0.0, cumul_to=b)


def test_from_valued_healpix_cells_bayestar_and_split():
    fits_mom_filename = "./resources/bayestar.multiorder.fits"
    moc = MOC.from_multiordermap_fits_file(
        fits_mom_filename,
        cumul_from=0.0,
        cumul_to=0.9,
    )
    count = moc.split_count()
    assert count == 2
    mocs = list(moc.split())
    assert len(mocs) == 2
    for moc in mocs:
        assert moc.max_order == 11


def test_probability_in_multiordermap():
    # from path
    moc = MOC.from_str("0/4")
    fits_mom_filename = "./resources/bayestar.multiorder.fits"
    proba = moc.probability_in_multiordermap(fits_mom_filename)
    assert np.isclose(proba, 0.20877154164727782)

    # has no mask
    mom = QTable()
    mom["UNIQ"] = [4 + x for x in range(20)]
    mom["PROBDENSITY"] = [x / 10 for x in range(20)]
    proba = moc.probability_in_multiordermap(mom)
    assert np.isclose(proba, 0.41887902047863906)

    # is not a valid mom or path
    with pytest.raises(
        ValueError,
        match="An argument of type 'str', 'pathlib.Path', or "
        "'astropy.table.Table' is expected. Got '<class 'int'>'",
    ):
        moc.probability_in_multiordermap(1)


def test_probabilities_in_multiordermap():
    moc = MOC.from_str("0/4")
    mocs = [moc, moc]

    mom = QTable()
    mom["UNIQ"] = [4 + x for x in range(20)]
    mom["PROBDENSITY"] = [x / 10 for x in range(20)]
    proba = moc.probability_in_multiordermap(mom)
    list_probas = MOC.probabilities_in_multiordermap(mocs, mom)
    assert proba == list_probas[0]


def test_sum_in_multiordermap():
    all_sky = MOC.from_str("0/0-11")
    mom = QTable()
    range_20 = range(20)
    mom["UNIQ"] = [4 * 4**3 + x for x in range_20]
    mom["TO_SUM"] = range_20
    assert all_sky.sum_in_multiordermap(mom, "TO_SUM") == sum(range_20)


def test_values_and_weights_in_multiordermap():
    all_sky = MOC.from_str("0/0-11")
    mom = QTable()
    range_20 = range(20)
    uniq = np.array([4 * 4**3 + x for x in range_20])
    mom["UNIQ"] = uniq
    mom["values"] = range_20
    values, weights = all_sky.values_and_weights_in_multiordermap(mom, "values")
    assert all(values == mom["values"])
    assert all(np.isclose(weight, 4 * math.pi / (12 * 4**3)) for weight in weights)

    one_cell = MOC.from_str("3/0")
    mom = QTable()
    mom["values"] = [0, 1, 2, 3]
    mom["UNIQ"] = [4 * 4**2 + x for x in [0, 1, 2, 3]]  # corresponds to "1/0"
    values, weights = one_cell.values_and_weights_in_multiordermap(mom, "values")
    assert all(values == np.array([0]))
    assert all(np.isclose(weights, one_cell.sky_fraction * 4 * math.pi))


def test_mask_uniq():
    uniq = [4 * 4**3 + x for x in range(8)]
    moc = MOC.from_str("3/4-20")
    assert all(
        moc.mask_uniq(uniq) == [False, False, False, False, True, True, True, True],
    )

    # fully covered should have less matches
    cone1 = MOC.from_cone(20 * u.deg, 20 * u.deg, radius=2 * u.deg, max_depth=10)
    cone2 = MOC.from_cone(21 * u.deg, 21 * u.deg, radius=2 * u.deg, max_depth=10)
    assert sum(cone1.mask_uniq(cone2.uniq_hpx)) > sum(
        cone1.mask_uniq(cone2.uniq_hpx, fully_covered_only=True),
    )


def test_from_stcs():
    moc1 = MOC.from_stcs("Circle ICRS 147.6 69.9 0.4", 14, 2)
    moc2 = MOC.from_cone(
        lon=147.6 * u.deg,
        lat=69.9 * u.deg,
        radius=Angle(0.4, u.deg),
        max_depth=14,
    )
    assert moc1 == moc2
