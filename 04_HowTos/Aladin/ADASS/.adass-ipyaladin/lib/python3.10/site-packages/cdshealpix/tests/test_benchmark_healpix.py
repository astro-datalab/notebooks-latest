# Astropy tools
import astropy.units as u
import astropy_healpix
from astropy.coordinates import Latitude, Longitude

import numpy as np
import pytest

from ..nested.healpix import (
    cone_search,
    healpix_to_lonlat,
    healpix_to_skycoord,
    lonlat_to_healpix,
    neighbours,
    vertices,
)


@pytest.mark.benchmark(group="lonlat_to_healpix_broadcast")
def test_lonlat_to_healpix_broadcast(benchmark):
    size = 10000
    depth = np.random.randint(low=0, high=29, size=10)

    lon = Longitude(np.random.rand(size) * 360, u.deg)
    lat = Latitude(np.random.rand(size) * 180 - 90, u.deg)

    benchmark(
        lonlat_to_healpix,
        lon=lon[:, np.newaxis],
        lat=lat[:, np.newaxis],
        depth=depth[np.newaxis, :],
        num_threads=1,
    )


@pytest.mark.benchmark(group="lonlat_to_healpix_broadcast")
def test_lonlat_to_healpix_broadcast_multithread(benchmark):
    size = 10000
    depth = np.random.randint(low=0, high=29, size=10)

    lon = Longitude(np.random.rand(size) * 360, u.deg)
    lat = Latitude(np.random.rand(size) * 180 - 90, u.deg)

    benchmark(
        lonlat_to_healpix,
        lon=lon[:, np.newaxis],
        lat=lat[:, np.newaxis],
        depth=depth[np.newaxis, :],
    )


@pytest.mark.benchmark(group="lonlat_to_healpix_broadcast")
def test_lonlat_to_healpix_astropy_broadcast(benchmark):
    depth = np.random.randint(low=0, high=29, size=10)
    nside = 2**depth

    size = 10000
    lon = np.random.rand(size) * 360 * u.deg
    lat = (np.random.rand(size) * 180 - 90) * u.deg

    benchmark(
        astropy_healpix.core.lonlat_to_healpix,
        lon=lon[:, np.newaxis],
        lat=lat[:, np.newaxis],
        nside=nside[np.newaxis, :],
        order="nested",
    )


@pytest.mark.benchmark(group="lonlat_to_healpix")
def test_lonlat_to_healpix(benchmark):
    size = 10000
    depth = 12
    lon = Longitude(np.random.rand(size) * 360, u.deg)
    lat = Latitude(np.random.rand(size) * 180 - 90, u.deg)

    benchmark(lonlat_to_healpix, lon=lon, lat=lat, depth=depth, num_threads=1)


@pytest.mark.benchmark(group="lonlat_to_healpix")
def test_lonlat_to_healpix_multithread(benchmark):
    size = 10000
    depth = 12
    lon = Longitude(np.random.rand(size) * 360, u.deg)
    lat = Latitude(np.random.rand(size) * 180 - 90, u.deg)

    benchmark(lonlat_to_healpix, lon=lon, lat=lat, depth=depth)


@pytest.mark.benchmark(group="lonlat_to_healpix")
def test_lonlat_to_healpix_astropy(benchmark):
    size = 10000
    depth = 12
    nside = 1 << depth
    lon = np.random.rand(size) * 360 * u.deg
    lat = (np.random.rand(size) * 180 - 90) * u.deg

    benchmark(
        astropy_healpix.core.lonlat_to_healpix,
        lon=lon,
        lat=lat,
        nside=nside,
        order="nested",
    )


@pytest.mark.benchmark(group="healpix_to_lonlat_broadcast")
def test_healpix_to_lonlat_broadcast(benchmark):
    size = 1000
    depth = np.random.randint(low=0, high=29, size=10)
    ipixels = np.random.randint(12 * 4 ** np.min(depth), size=size)

    benchmark(
        healpix_to_lonlat,
        ipix=ipixels[:, np.newaxis],
        depth=depth[np.newaxis, :],
        num_threads=1,
    )


@pytest.mark.benchmark(group="healpix_to_lonlat_broadcast")
def test_healpix_to_lonlat_broadcast_multithread(benchmark):
    size = 1000
    depth = np.random.randint(low=0, high=29, size=10)
    ipixels = np.random.randint(12 * 4 ** np.min(depth), size=size)

    benchmark(
        healpix_to_lonlat, ipix=ipixels[:, np.newaxis], depth=depth[np.newaxis, :]
    )


@pytest.mark.benchmark(group="healpix_to_lonlat_broadcast")
def test_healpix_to_lonlat_astropy_broadcast(benchmark):
    size = 1000
    depth = np.random.randint(low=0, high=29, size=10)
    nside = 2**depth
    ipixels = np.random.randint(12 * 4 ** np.min(depth), size=size)

    benchmark(
        astropy_healpix.core.healpix_to_lonlat,
        healpix_index=ipixels[:, np.newaxis],
        nside=nside[np.newaxis, :],
        order="nested",
    )


@pytest.mark.benchmark(group="healpix_to_lonlat")
def test_healpix_to_lonlat(benchmark):
    size = 10000
    depth = 12
    ipixels = np.random.randint(12 * 4 ** (depth), size=size)

    benchmark(healpix_to_lonlat, ipix=ipixels, depth=depth, num_threads=1)


@pytest.mark.benchmark(group="healpix_to_lonlat")
def test_healpix_to_lonlat_multithread(benchmark):
    size = 10000
    depth = 12
    ipixels = np.random.randint(12 * 4 ** (depth), size=size)

    benchmark(healpix_to_lonlat, ipix=ipixels, depth=depth)


@pytest.mark.benchmark(group="healpix_to_lonlat")
def test_healpix_to_lonlat_astropy(benchmark):
    size = 10000
    depth = 12
    nside = 1 << depth
    ipixels = np.random.randint(12 * 4 ** (depth), size=size)

    benchmark(
        astropy_healpix.core.healpix_to_lonlat,
        healpix_index=ipixels,
        nside=nside,
        order="nested",
    )


def test_healpix_to_skycoord():
    healpix_to_skycoord(ipix=[0, 2, 4], depth=0)


@pytest.mark.benchmark(group="vertices")
def test_healpix_vertices_lonlat(benchmark):
    depth = 12
    size = 100000
    ipixels = np.random.randint(12 * 4 ** (depth), size=size)

    lon, lat = benchmark(vertices, ipix=ipixels, depth=depth, num_threads=1)


@pytest.mark.benchmark(group="vertices")
def test_healpix_vertices_lonlat_multithread(benchmark):
    depth = 12
    size = 100000
    ipixels = np.random.randint(12 * 4 ** (depth), size=size)

    lon, lat = benchmark(vertices, ipix=ipixels, depth=depth)


@pytest.mark.benchmark(group="vertices")
def test_healpix_vertices_lonlat_astropy(benchmark):
    depth = 12
    size = 100000
    ipixels = np.random.randint(12 * 4 ** (depth), size=size)

    lon2, lat2 = benchmark(
        astropy_healpix.core.boundaries_lonlat,
        healpix_index=ipixels,
        nside=(1 << depth),
        step=1,
        order="nested",
    )


@pytest.mark.benchmark(group="neighbours")
def test_healpix_neighbours(benchmark):
    depth = 12
    size = 100000
    ipixels = np.random.randint(12 * 4 ** (depth), size=size)

    benchmark(neighbours, ipix=ipixels, depth=depth)


lon_c = np.random.rand(1)[0] * 360 * u.deg
lat_c = (np.random.rand(1)[0] * 180 - 90) * u.deg
radius_c = np.random.rand(1)[0] * 45 * u.deg
depth_c = 5


@pytest.mark.benchmark(group="cone_search")
def test_cone_search(benchmark):
    benchmark.pedantic(
        cone_search,
        kwargs={
            "lon": Longitude(lon_c),
            "lat": Latitude(lat_c),
            "radius": radius_c,
            "depth": depth_c,
        },
        iterations=10,
        rounds=100,
    )


@pytest.mark.benchmark(group="cone_search")
def test_cone_search_astropy(benchmark):
    nside = 1 << depth_c
    benchmark.pedantic(
        astropy_healpix.core.healpix_cone_search,
        kwargs={
            "lon": lon_c,
            "lat": lat_c,
            "radius": radius_c,
            "nside": nside,
            "order": "nested",
        },
        iterations=10,
        rounds=100,
    )
