# Standard Library
import pathlib
import re

# Astropy tools
# General Astronomy tools
import astropy.units as u
from astropy.coordinates import Angle, Latitude, Longitude, SkyCoord

import numpy as np
import pytest

from .. import from_ring, to_ring
from ..nested.healpix import (
    bilinear_interpolation,
    cone_search,
    box_search,
    zone_search,
    elliptical_cone_search,
    external_neighbours,
    healpix_to_lonlat,
    healpix_to_skycoord,
    healpix_to_xy,
    lonlat_to_healpix,
    lonlat_to_xy,
    neighbours,
    polygon_search,
    skycoord_to_healpix,
    vertices,
    xy_to_lonlat,
)


@pytest.mark.parametrize("size", [1, 10, 100, 1000, 10000, 100000])
def test_lonlat_to_healpix(size):
    depth = np.random.randint(30)
    lon = Longitude(np.random.rand(size) * 360, u.deg)
    lat = Latitude(np.random.rand(size) * 180 - 90, u.deg)

    ipixels, dx, dy = lonlat_to_healpix(
        lon=lon, lat=lat, depth=depth, return_offsets=True
    )
    ipixels, dx, dy = skycoord_to_healpix(
        SkyCoord(lon, lat, frame="icrs"),
        depth=depth,
        return_offsets=True,
    )

    npix = 12 * 4 ** (depth)
    assert ((ipixels >= 0) & (ipixels < npix)).all()
    assert ((dx >= 0) & (dx <= 1)).all()
    assert ((dy >= 0) & (dy <= 1)).all()


@pytest.mark.parametrize("size", [1, 10, 100, 1000, 10000, 100000])
def test_healpix_to_lonlat(size):
    depth = np.random.randint(30)
    ipixels = np.random.randint(12 * 4**depth, size=size, dtype=np.uint64)

    lon, lat = healpix_to_lonlat(ipix=ipixels, depth=depth)
    assert lon.shape == lat.shape


@pytest.mark.parametrize(
    "ipix",
    [
        np.arange(100)[::2],
        np.broadcast_to(np.array([1, 2, 3]), (10000, 3)),
        np.broadcast_to(3, (3, 4, 5)),
    ],
)
def test_healpix_to_lonlat_on_broadcasted_arrays(ipix):
    depth = 12
    lon, _ = healpix_to_lonlat(ipix=ipix, depth=depth)
    assert lon.shape == ipix.shape


def test_healpix_to_lonlat_on_broadcasted_arrays2():
    depth = np.arange(2)
    ipix = np.arange(3)
    lon, lat = healpix_to_lonlat(ipix[:, np.newaxis], depth[np.newaxis, :])
    assert lon.shape == lat.shape and lon.shape == (3, 2)


def test_invalid_depth_exception():
    size = 1000
    ipix = np.zeros(size)
    match = re.escape("Depth must be in the [0, 29] closed range")
    with pytest.raises(ValueError, match=match):
        healpix_to_lonlat(ipix, -2)
    with pytest.raises(ValueError, match=match):
        healpix_to_lonlat(ipix, 30)
    with pytest.raises(ValueError, match=match):
        lonlat_to_healpix(Longitude(0, u.deg), Latitude(0, u.deg), -2)
    with pytest.raises(ValueError, match=match):
        vertices(ipix, -2)
    with pytest.raises(ValueError, match=match):
        neighbours(ipix, -2)
    with pytest.raises(ValueError, match=match):
        cone_search(Longitude(0, u.deg), Latitude(0, u.deg), 15 * u.deg, -2)


def test_lonlat_shape_exception():
    lon = Longitude([2, 5], u.deg)
    lat = Latitude([5], u.deg)
    with pytest.raises(
        ValueError,
        match=re.escape(
            "'lon' and 'lat' should have the same shape but are of shapes (2,) and (1,)"
        ),
    ):
        lonlat_to_healpix(lon, lat, 12)


@pytest.mark.parametrize("depth", [0, 12, 24])
def test_invalid_ipix_exception(depth):
    npix = 12 * 4**depth
    invalid_ipix1 = np.array([-20, 0, 11])
    invalid_ipix2 = np.array([0, npix + 1, 11])
    match = "The input HEALPix array contains values out of *"
    with pytest.raises(ValueError, match=match):
        healpix_to_lonlat(invalid_ipix1, depth)
    with pytest.raises(ValueError, match=match):
        healpix_to_lonlat(invalid_ipix2, depth)
    with pytest.raises(ValueError, match=match):
        vertices(invalid_ipix1, depth)
    with pytest.raises(ValueError, match=match):
        neighbours(invalid_ipix1, depth)


def test_healpix_to_skycoord():
    ipix = np.array([0, 2, 4])
    skycoord = healpix_to_skycoord(ipix=ipix, depth=0)
    assert skycoord.icrs.ra.shape == skycoord.icrs.dec.shape


def test_vertices_lonlat():
    depth = 12
    size = 10000
    ipixels = np.random.randint(12 * 4**depth, size=size)

    lon, lat = vertices(ipix=ipixels, depth=depth)
    assert lon.shape == lat.shape
    assert lon.shape == (size, 4)


def test_neighbours():
    depth = 0
    size = 10000
    ipixels = np.random.randint(12 * 4 ** (depth), size=size)

    neigh = neighbours(ipix=ipixels, depth=depth)
    assert neigh.shape == (size, 9)

    npix = 12 * 4 ** (depth)
    assert ((neigh >= -1) & (neigh < npix)).all()


def test_cone_search():
    lon = Longitude(np.random.rand(1)[0] * 360, u.deg)
    lat = Latitude(np.random.rand(1)[0] * 180 - 90, u.deg)
    radius = (np.random.rand(1)[0] * 45) * u.deg
    max_depth = 5

    ipix, depth, _ = cone_search(
        lon=lon, lat=lat, radius=radius, depth=max_depth, flat=True
    )

    npix = 12 * 4 ** (max_depth)
    assert ((depth >= 0) & (depth <= max_depth)).all()
    assert ((ipix >= 0) & (ipix < npix)).all()

    with pytest.raises(
        ValueError, match="The longitude, latitude and radius must be scalar objects"
    ):
        cone_search(
            Longitude([5, 4], u.deg), Latitude([5, 4], u.deg), 15 * u.deg, depth=12
        )


def test_box_search():
    lon = Longitude(0 * u.deg)
    lat = Latitude(0 * u.deg)
    _, depth, fully_covered = box_search(lon, lat, 1 * u.deg, 1 * u.deg, 0 * u.deg)
    assert all(depth <= 14)
    assert not all(fully_covered)


def test_zone_search():
    lon_min = Longitude(-5 * u.deg)
    lat_min = Latitude(-5 * u.deg)
    lon_max = Longitude(5 * u.deg)
    lat_max = Latitude(5 * u.deg)
    ipix, depth, fully_covered = zone_search(
        lon_min, lat_min, lon_max, lat_max, depth=6
    )
    centers_lon, center_lat = healpix_to_lonlat(
        ipix[fully_covered], depth[fully_covered]
    )
    assert all(
        centers_lon.wrap_at(180 * u.deg).radian >= lon_min.wrap_at(180 * u.deg).radian
    )


@pytest.mark.parametrize("size", [0, 1, 2, 3, 5, 6, 9])
def test_polygon_search(size):
    max_depth = 12
    lon = Longitude(np.random.rand(size) * 360, u.deg)
    lat = Latitude(np.random.rand(size) * 180 - 90, u.deg)

    if size < 3:
        with pytest.raises(
            ValueError,
            match="There must be at least 3 vertices in order to form a polygon",
        ):
            polygon_search(lon=lon, lat=lat, depth=max_depth)
    else:
        ipix, depth, _ = polygon_search(lon=lon, lat=lat, depth=max_depth)

        npix = 12 * 4 ** (max_depth)
        assert ((depth >= 0) & (depth <= max_depth)).all()
        assert ((ipix >= 0) & (ipix < npix)).all()


# From https://github.com/cds-astro/cds-healpix-python/issues/10
def test_polygon_search_issue10():
    coords = SkyCoord(
        [
            (353.8156714, -56.33202193),
            (6.1843286, -56.33202193),
            (5.27558041, -49.49378172),
            (354.72441959, -49.49378172),
        ],
        unit=u.deg,
    )
    polygon_search(Longitude(coords.ra), Latitude(coords.dec), 12)


# From https://github.com/cds-astro/mocpy/issues/57
def test_polygon_search_issue57():
    coords = np.load(pathlib.Path(__file__).parent / "moc_coords.npy")
    polygon_search(Longitude(coords[:, 0], u.deg), Latitude(coords[:, 1], u.deg), 9)


def test_polygon_search_not_enough_vertices_exception():
    # 4 total vertices but only 2 distincts. This should fail
    with pytest.raises(
        ValueError,
        match="There must be at least 3 distinct vertices in order to form a polygon",
    ):
        polygon_search(
            Longitude([1, 1, 2, 1], u.deg), Latitude([1, 1, 3, 1], u.deg), depth=12
        )

    # 4 total vertices and 3 distincts. This should pass
    polygon_search(
        Longitude([1, 1, 2, 1], u.deg), Latitude([1, 1, 3, 2], u.deg), depth=12
    )


# Following an error spotted by the hips2fits service that usually deal with big polygons
def test_polygon_search_big_polygon():
    ipix, _, _ = polygon_search(
        Longitude([268.84102386, 299.40278164, 66.0951825, 96.66953469], u.deg),
        Latitude([-45.73283624, 38.47909742, 38.76894431, -45.43470273], u.deg),
        depth=0,
        flat=True,
    )
    assert (ipix == np.asarray([0, 3, 4, 5, 7, 8, 9, 10, 11])).all()


def test_elliptical_cone_search():
    lon = Longitude(0, u.deg)
    lat = Latitude(0, u.deg)
    a = Angle(50, unit="deg")
    b = Angle(5, unit="deg")
    pa = Angle(0, unit="deg")
    max_depth = 12

    ipix, depth, _ = elliptical_cone_search(lon, lat, a, b, pa, max_depth)

    npix = 12 * 4 ** (max_depth)
    assert ((depth >= 0) & (depth <= max_depth)).all()
    assert ((ipix >= 0) & (ipix < npix)).all()


@pytest.mark.parametrize(
    "depth,ipix,expected_border_cells,expected_corner_cells",
    [
        (
            0,
            0,
            np.array([90, 91, 94, 95, 26, 27, 30, 31, 53, 55, 61, 63, 69, 71, 77, 79]),
            np.array([143, -1, 47, -1]),
        ),
        (
            27,
            0,
            np.array(
                [
                    1633305464859699882,
                    1633305464859699883,
                    1633305464859699886,
                    1633305464859699887,
                    16,
                    18,
                    24,
                    26,
                    32,
                    33,
                    36,
                    37,
                    1248998296657417557,
                    1248998296657417559,
                    1248998296657417565,
                    1248998296657417567,
                ]
            ),
            np.array(
                [2594073385365405695, 1633305464859699898, 48, 1248998296657417589]
            ),
        ),
    ],
)
def test_external_neighbours(depth, ipix, expected_border_cells, expected_corner_cells):
    delta_depth = 2
    ipix_border_cells, ipix_corner_cells = external_neighbours(ipix, depth, delta_depth)
    assert (expected_border_cells == ipix_border_cells).all()
    assert (expected_corner_cells == ipix_corner_cells).all()


@pytest.mark.parametrize(
    "ipix, depth, expected_x, expected_y",
    [
        (
            np.arange(12),
            0,
            np.array(
                [1.0, 3.0, 5.0, 7.0, 0.0, 2.0, 4.0, 6.0, 1.0, 3.0, 5.0, 7.0],
                dtype=np.float64,
            ),
            np.array(
                [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0],
                dtype=np.float64,
            ),
        ),
        ([], 0, np.array([], dtype=np.float64), np.array([], dtype=np.float64)),
        (0, 0, np.array([1.0], dtype=np.float64), np.array([1.0], dtype=np.float64)),
    ],
)
def test_healpix_to_xy(ipix, depth, expected_x, expected_y):
    x, y = healpix_to_xy(ipix, depth)

    assert (x == expected_x).all()
    assert (y == expected_y).all()


def test_healpix_to_xy_expection():
    with pytest.raises(ValueError):
        healpix_to_xy(np.array([]), -5)

    with pytest.raises(ValueError):
        healpix_to_xy(np.array([-5]), 12)


def assert_equal_array(first, second, tol=1e-8):
    assert ((first - second) < tol).all()


@pytest.mark.parametrize(
    "lon, lat, expected_x, expected_y",
    [
        (
            Longitude(np.array([0.0, 0.78539816, 1.57079633]), u.rad),
            Latitude(np.array([-0.72972766, 0.0, 0.72972766]), u.rad),
            np.array([0.0, 1.0, 2.0]),
            np.array([-1.0, 0.0, 1.0]),
        )
    ],
)
def test_lonlat_to_xy(lon, lat, expected_x, expected_y):
    x, y = lonlat_to_xy(lon, lat)
    assert_equal_array(x, expected_x)
    assert_equal_array(y, expected_y)


@pytest.mark.parametrize(
    "x, y, expected_lon, expected_lat",
    [
        (
            np.array([0.0, 1.0, 2.0]),
            np.array([-1.0, 0.0, 1.0]),
            np.array([0.0, 0.78539816, 1.57079633]) * u.rad,
            np.array([-0.72972766, 0.0, 0.72972766]) * u.rad,
        )
    ],
)
def test_xy_to_lonlat(x, y, expected_lon, expected_lat):
    lon, lat = xy_to_lonlat(x, y)
    assert_equal_array(lon.to_value(u.rad), expected_lon.to_value(u.rad))
    assert_equal_array(lat.to_value(u.rad), expected_lat.to_value(u.rad))


# -- Test to_ring & from_ring
@pytest.mark.parametrize(
    "pix, depth, expected_ring_pix",
    [
        (np.arange(12), 0, np.arange(12)),
        (4, 1, 15),
        (np.array([4, 5, 6, 7]), 1, np.array([15, 7, 6, 1])),
        (np.array([28, 29, 30, 31]), 1, np.array([34, 26, 25, 18])),
    ],
)
def test_to_ring(pix, depth, expected_ring_pix):
    ring_pix = to_ring(pix, depth)
    assert (ring_pix == expected_ring_pix).all()


@pytest.mark.parametrize(
    "pix, depth, expected_nested_pix",
    [
        (np.arange(12), 0, np.arange(12)),
        (15, 1, 4),
        (np.array([15, 7, 6, 1]), 1, np.array([4, 5, 6, 7])),
        (np.array([34, 26, 25, 18]), 1, np.array([28, 29, 30, 31])),
    ],
)
def test_from_ring(pix, depth, expected_nested_pix):
    nested_pix = from_ring(pix, depth)
    assert (nested_pix == expected_nested_pix).all()


@pytest.mark.parametrize("size", [1, 10, 100, 1000, 10000, 100000])
def test_from_vs_to_ring(size):
    depth = np.random.randint(30)
    nside = 1 << depth

    ipixels = np.random.randint(12 * (nside**2), size=size, dtype="uint64")
    ring = to_ring(ipix=ipixels, depth=depth)
    ipixels_result = from_ring(ipix=ring, depth=depth)
    assert (ipixels == ipixels_result).all()


@pytest.mark.parametrize("depth", [5, 0, 7, 12, 20, 29])
def test_bilinear_interpolation(depth):
    size = 1000

    lon = Longitude(np.random.rand(size) * 360, u.deg)
    lat = Latitude(np.random.rand(size) * 180 - 90, u.deg)

    ipix, weights = bilinear_interpolation(lon, lat, depth)

    assert ((weights >= 0.0) & (weights <= 1.0)).all()
    assert weights.sum() == ipix.shape[0]
    assert ((ipix >= 0) & (ipix < 12 * (4**depth))).all()


def test_bilinear_interpolation_accepts_nan():
    lon = Longitude([10, np.nan], unit="deg")
    lat = Latitude([5, np.nan], unit="deg")
    depth = 5

    ipix, _ = bilinear_interpolation(lon, lat, depth)

    assert np.all(ipix.mask == [[False, False, False, False], [True, True, True, True]])
