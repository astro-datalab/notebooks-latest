"""Manipulation of HEALPix in ring configuration."""

# Astropy tools
import astropy.units as u
from astropy.coordinates import Latitude, Longitude, SkyCoord

import numpy as np

from .. import cdshealpix
from ..utils import _validate_lonlat

# Do not fill by hand :)
# > egrep "^ *def" healpix.py | cut -c 5- | egrep -v '^_'| cut -d '(' -f 1 | sed -r "s/^(.*)$/ '\1'/" | tr '\n' ','
__all__ = [
    "lonlat_to_healpix",
    "skycoord_to_healpix",
    "healpix_to_lonlat",
    "healpix_to_skycoord",
    "healpix_to_xy",
    "vertices",
    "vertices_skycoord",
]

# Raise a ValueError exception if the input
# HEALPix cells array contains invalid values
def _check_ipixels(data, nside):
    npix = 12 * (nside**2)
    if (data >= npix).any() or (data < 0).any():
        valid_ipix = np.stack((np.zeros(npix.shape), npix)).T
        raise ValueError(
            f"The input HEALPix array contains values out of {valid_ipix}."
        )


@_validate_lonlat
def lonlat_to_healpix(lon, lat, nside, return_offsets=False, num_threads=0):
    r"""Get the HEALPix indexes that contains specific sky coordinates.

    The ``nside`` of the returned HEALPix cell indexes must be specified. This
    method is wrapped around the `hash <https://docs.rs/cdshealpix/0.1.5/cdshealpix/nested/struct.Layer.html#method.hash>`__
    method from the `cdshealpix Rust crate <https://crates.io/crates/cdshealpix>`__.

    Parameters
    ----------
    lon : `astropy.coordinates.Longitude`
        The longitudes of the sky coordinates.
    lat : `astropy.coordinates.Latitude`
        The latitudes of the sky coordinates.
    nside : `numpy.ndarray`
        The nside of the returned HEALPix cell indexes.
    return_offsets : bool, optional
        If set to `True`, returns a tuple made of 3 elements, the HEALPix cell
        indexes and the dx, dy arrays telling where the (``lon``, ``lat``) coordinates
        passed are located on the cells. ``dx`` and ``dy`` are :math:`\in [0, 1]`
    num_threads : int, optional
        Specifies the number of threads to use for the computation. Default to 0 means
        it will choose the number of threads based on the RAYON_NUM_THREADS environment variable (if set),
        or the number of logical CPUs (otherwise)

    Returns
    -------
    ipix : `numpy.ndarray`
        A numpy array containing all the HEALPix cell indexes stored as `np.uint64`.

    Raises
    ------
    ValueError
        When the number of longitudes and latitudes given do not match.

    Examples
    --------
    >>> from cdshealpix.ring import lonlat_to_healpix
    >>> from astropy.coordinates import Longitude, Latitude
    >>> import astropy.units as u
    >>> import numpy as np
    >>> lon = Longitude([0, 50, 25], u.deg)
    >>> lat = Latitude([6, -12, 45], u.deg)
    >>> depth = np.array([12, 14])
    >>> nside = 2 ** depth
    >>> ipix = lonlat_to_healpix(lon[:, np.newaxis], lat[:, np.newaxis], nside[np.newaxis, :])
    """
    # Check arrays
    #
    # We could have continued to use `.to_value(u.rad)` instead of `.rad`.
    # Although `to_value` is more generical (method of Quantity),
    # Longitude/Latitude ensure that the values the contain are in the correct ranges.
    lon = np.atleast_1d(lon.rad)
    lat = np.atleast_1d(lat.rad)
    nside = np.atleast_1d(nside)

    if (nside < 1).any() or (nside > (1 << 29)).any():
        raise ValueError("nside must be in the [1, (1 << 29)[ closed range")

    # Broadcasting
    lon, lat, nside = np.broadcast_arrays(lon, lat, nside)

    # Allocation of the array containing the resulting coordinates
    num_ipix = lon.shape
    ipix = np.empty(num_ipix, dtype=np.uint64)
    dx = np.empty(num_ipix, dtype=np.float64)
    dy = np.empty(num_ipix, dtype=np.float64)

    # Call the Rust extension
    nside = nside.astype(np.uint32)
    num_threads = np.uint16(num_threads)
    cdshealpix.lonlat_to_healpix_ring(nside, lon, lat, ipix, dx, dy, num_threads)

    if return_offsets:
        return ipix, dx, dy
    return ipix


def skycoord_to_healpix(skycoord, nside, return_offsets=False, num_threads=0):
    r"""Get the HEALPix indexes that contains specific sky coordinates.

    The ``nside`` of the returned HEALPix cell indexes must be specified.
    This method is wrapped around the
    `hash <https://docs.rs/cdshealpix/0.1.5/cdshealpix/nested/struct.Layer.html#method.hash>`__
    method from the `cdshealpix Rust crate <https://crates.io/crates/cdshealpix>`__.

    Parameters
    ----------
    skycoord : `astropy.coordinates.SkyCoord`
        The sky coordinates.
    nside : `numpy.ndarray`
        The nside of the returned HEALPix cell indexes.
    return_offsets : bool, optional
        If set to `True`, returns a tuple made of 3 elements, the HEALPix cell
        indexes and the dx, dy arrays telling where the (``lon``, ``lat``) coordinates
        passed are located in the cells. ``dx`` and ``dy`` are :math:`\in [0, 1]`
    num_threads : int, optional
        Specifies the number of threads to use for the computation. Default to 0 means
        it will choose the number of threads based on the RAYON_NUM_THREADS environment variable (if set),
        or the number of logical CPUs (otherwise)

    Returns
    -------
    ipix : `numpy.ndarray`
        A numpy array containing all the HEALPix cell indexes stored as `np.uint64`.

    Raises
    ------
    ValueError
        When the number of longitudes and latitudes given do not match.

    Examples
    --------
    >>> from cdshealpix.ring import skycoord_to_healpix
    >>> import astropy.units as u
    >>> from astropy.coordinates import SkyCoord
    >>> import numpy as np
    >>> skycoord = SkyCoord([0, 50, 25] * u.deg, [6, -12, 45] * u.deg, frame="icrs")
    >>> depth = 12
    >>> ipix = skycoord_to_healpix(skycoord, 1 << depth)
    """
    return lonlat_to_healpix(
        Longitude(skycoord.icrs.ra),
        Latitude(skycoord.icrs.dec),
        nside,
        return_offsets,
        num_threads,
    )


def healpix_to_lonlat(ipix, nside, dx=0.5, dy=0.5, num_threads=0):
    r"""Get the longitudes and latitudes of the center of some HEALPix cells at a given depth.

    This method does the opposite transformation of `lonlat_to_healpix`.
    It's wrapped around the `center <https://docs.rs/cdshealpix/0.1.5/cdshealpix/nested/struct.Layer.html#method.center>`__
    method from the `cdshealpix Rust crate <https://crates.io/crates/cdshealpix>`__.

    Parameters
    ----------
    ipix : `numpy.ndarray`
        The HEALPix cell indexes given as a `np.uint64` numpy array.
    nside : `numpy.ndarray`
        The nside of the HEALPix cells.
    dx : float, optional
        The offset position :math:`\in [0, 1[` along the X axis. By default, `dx=0.5`
    dy : float, optional
        The offset position :math:`\in [0, 1[` along the Y axis. By default, `dy=0.5`
    num_threads : int, optional
        Specifies the number of threads to use for the computation. Default to 0 means
        it will choose the number of threads based on the RAYON_NUM_THREADS environment variable (if set),
        or the number of logical CPUs (otherwise)

    Returns
    -------
    lon, lat : (`astropy.coordinates.Longitude`, `astropy.coordinates.Latitude`)
        The sky coordinates of the center of the HEALPix cells given as a longitude, latitude tuple.

    Raises
    ------
    ValueError
        When the HEALPix cell indexes given have values out of :math:`[0, 12` x :math:`N_{side} ^ 2[`.

    Examples
    --------
    >>> from cdshealpix.ring import healpix_to_lonlat
    >>> import numpy as np
    >>> ipix = np.array([42, 6, 10], dtype=np.uint64)
    >>> depth = np.array([2, 12])
    >>> nside = 2 ** depth
    >>> lon, lat = healpix_to_lonlat(ipix[:, np.newaxis], nside[np.newaxis, :])
    """
    # Check arrays
    ipix = np.atleast_1d(ipix)
    nside = np.atleast_1d(nside)
    if (nside < 1).any() or (nside > (1 << 29)).any():
        raise ValueError("nside must be in the [1, (1 << 29)[ closed range")

    if dx < 0 or dx >= 1:
        raise ValueError("dx must be between [0, 1[")

    if dy < 0 or dy >= 1:
        raise ValueError("dy must be between [0, 1[")

    _check_ipixels(data=ipix, nside=nside)

    # Broadcasting
    ipix, nside = np.broadcast_arrays(ipix, nside)

    size_skycoords = ipix.shape
    # Allocation of the array containing the resulting coordinates
    lon = np.zeros(size_skycoords)
    lat = np.zeros(size_skycoords)

    # Call the Rust extension
    nside = nside.astype(np.uint32)
    ipix = ipix.astype(np.uint64)
    num_threads = np.uint16(num_threads)

    cdshealpix.healpix_to_lonlat_ring(nside, ipix, dx, dy, lon, lat, num_threads)
    return Longitude(lon, u.rad), Latitude(lat, u.rad)


def healpix_to_skycoord(ipix, nside, dx=0.5, dy=0.5, num_threads=0):
    r"""Get the sky coordinates of the center of some HEALPix cells at a given nside.

    This method does the opposite transformation of `lonlat_to_healpix`.
    It is the equivalent of `healpix_to_lonlat` except that it returns `astropy.coordinates.SkyCoord` instead
    of `astropy.units.Quantity`.
    It's wrapped around the `center <https://docs.rs/cdshealpix/0.1.5/cdshealpix/nested/struct.Layer.html#method.center>`__
    method from the `cdshealpix Rust crate <https://crates.io/crates/cdshealpix>`__.

    Parameters
    ----------
    ipix : `numpy.ndarray`
        The HEALPix cell indexes given as a `np.uint64` numpy array.
    nside : `numpy.ndarray`
        The nside of the HEALPix cells.
    dx : float, optional
        The offset position :math:`\in [0, 1[` along the X axis. By default, `dx=0.5`
    dy : float, optional
        The offset position :math:`\in [0, 1[` along the Y axis. By default, `dy=0.5`
    num_threads : int, optional
        Specifies the number of threads to use for the computation. Default to 0 means
        it will choose the number of threads based on the RAYON_NUM_THREADS environment variable (if set),
        or the number of logical CPUs (otherwise)

    Returns
    -------
    skycoord : `astropy.coordinates.SkyCoord`
        The sky coordinates of the center of the HEALPix cells given as a `~astropy.coordinates.SkyCoord` object.

    Raises
    ------
    ValueError
        When the HEALPix cell indexes given have values out of :math:`[0, 12` x :math:`N_{side} ^ 2[`.

    Examples
    --------
    >>> from cdshealpix.ring import healpix_to_skycoord
    >>> import numpy as np
    >>> ipix = np.array([42, 6, 10])
    >>> depth = 12
    >>> nside = 2 ** depth
    >>> skycoord = healpix_to_skycoord(ipix, nside)
    """
    lon, lat = healpix_to_lonlat(ipix, nside, dx, dy, num_threads)
    return SkyCoord(ra=lon, dec=lat, frame="icrs", unit="rad")


def healpix_to_xy(ipix, nside, num_threads=0):
    r"""
    Project the center of a HEALPix cell to the xy-HEALPix plane.

    Parameters
    ----------
    ipix : `numpy.ndarray`
        The HEALPix cells which centers will be projected
    nside : `numpy.ndarray`
        The nside of the HEALPix cells
    num_threads : int, optional
        Specifies the number of threads to use for the computation. Default to 0 means
        it will choose the number of threads based on the RAYON_NUM_THREADS environment variable (if set),
        or the number of logical CPUs (otherwise)

    Returns
    -------
    x, y: (`numpy.ndarray`, `numpy.ndarray`)
        The position of the HEALPix centers in the xy-HEALPix plane.
        :math:`x \in [0, 8[` and :math:`y \in [-2, 2]`

    Examples
    --------
    >>> from cdshealpix.ring import healpix_to_xy
    >>> import astropy.units as u
    >>> import numpy as np
    >>> depth = np.array([0, 10])
    >>> nside = 2 ** depth
    >>> ipix = np.arange(12)
    >>> x, y = healpix_to_xy(ipix[:, np.newaxis], nside[np.newaxis, :])
    """
    # Check arrays
    ipix = np.atleast_1d(ipix)
    nside = np.atleast_1d(nside)

    if (nside < 1).any() or (nside > (1 << 29)).any():
        raise ValueError("nside must be in the [1, (1 << 29)[ closed range")

    _check_ipixels(data=ipix, nside=nside)

    # Broadcasting
    ipix, nside = np.broadcast_arrays(ipix, nside)

    # Allocation for the resulting arrays
    x = np.zeros(ipix.shape, dtype=np.float64)
    y = np.zeros(ipix.shape, dtype=np.float64)

    # Call the Rust extension
    ipix = ipix.astype(np.uint64)
    nside = nside.astype(np.uint32)
    num_threads = np.uint16(num_threads)
    cdshealpix.healpix_to_xy_ring(nside, ipix, x, y, num_threads)

    return x, y


def vertices(ipix, nside, step=1, num_threads=0):
    """Get the longitudes and latitudes of the vertices of some HEALPix cells at a given nside.

    This method returns the 4 vertices of each cell in `ipix`.
    This method is wrapped around the `vertices <https://docs.rs/cdshealpix/0.1.5/cdshealpix/nested/struct.Layer.html#method.vertices>`__
    method from the `cdshealpix Rust crate <https://crates.io/crates/cdshealpix>`__.

    Parameters
    ----------
    ipix : `numpy.ndarray`
        The HEALPix cell indexes given as a `np.uint64` numpy array.
    nside : int
        The nside of the HEALPix cells.
    step : int, optional
        The number of vertices returned per HEALPix side. By default it is set to 1 meaning that
        it will only return the vertices of the cell. 2 means that it will returns the vertices of
        the cell plus one more vertex per edge (the middle of it). More generally, the number
        of vertices returned is ``4 * step``.
    num_threads : int, optional
        Specifies the number of threads to use for the computation. Default to 0 means
        it will choose the number of threads based on the RAYON_NUM_THREADS environment variable (if set),
        or the number of logical CPUs (otherwise)

    Returns
    -------
    lon, lat : (`astropy.coordinates.Longitude`, `astropy.coordinates.Latitude`)
        The sky coordinates of the 4 vertices of the HEALPix cells.
        `lon` and `lat` are `~astropy.coordinates.Longitude` and `~astropy.coordinates.Latitude` instances respectively,
        containing a :math:`N` x :math:`4` numpy array where N is the number of HEALPix cell given in `ipix`.

    Warnings
    --------
    ``step`` is currently not implemented for the ring scheme. Therefore it is set to 1 by default.

    Raises
    ------
    ValueError
        When the HEALPix cell indexes given have values out of :math:`[0, 12` x :math:`N_{side} ^ 2[`.

    Examples
    --------
    >>> from cdshealpix.ring import vertices
    >>> import numpy as np
    >>> ipix = np.array([42, 6, 10])
    >>> depth = 12
    >>> lon, lat = vertices(ipix, (1 << depth))
    """
    if nside < 1 or nside > (1 << 29):
        raise ValueError("nside must be in the [1, (1 << 29)[ closed range")

    if step < 1:
        raise ValueError("The number of step must be >= 1")

    ipix = np.atleast_1d(ipix)
    _check_ipixels(data=ipix, nside=nside)
    ipix = ipix.astype(np.uint64)

    # Allocation of the array containing the resulting coordinates
    lon = np.zeros((*ipix.shape, 4 * step))
    lat = np.zeros((*ipix.shape, 4 * step))
    num_threads = np.uint16(num_threads)

    cdshealpix.vertices_ring(nside, ipix, step, lon, lat, num_threads)
    return Longitude(lon, u.rad), Latitude(lat, u.rad)


def vertices_skycoord(ipix, nside, step=1):
    """Get the sky coordinates of the vertices of some HEALPix cells at a given nside.

    This method returns the 4 vertices of each cell in `ipix`.
    This method is wrapped around the `vertices <https://docs.rs/cdshealpix/0.1.5/cdshealpix/nested/struct.Layer.html#method.vertices>`__
    method from the `cdshealpix Rust crate <https://crates.io/crates/cdshealpix>`__.

    Parameters
    ----------
    ipix : `numpy.ndarray`
        The HEALPix cell indexes given as a `np.uint64` numpy array.
    nside : int
        The nside of the HEALPix cells.
    step : int, optional
        The number of vertices returned per HEALPix side. By default it is set to 1 meaning that
        it will only return the vertices of the cell. 2 means that it will returns the vertices of
        the cell plus one more vertex per edge (the middle of it). More generally, the number
        of vertices returned is ``4 * step``.

    Returns
    -------
    vertices : `astropy.coordinates.SkyCoord`
        The sky coordinates of the 4 vertices of the HEALPix cells. `vertices` is a `~astropy.coordinates.SkyCoord` object
        containing a :math:`N` x :math:`4` numpy array where N is the number of HEALPix cells given in `ipix`.

    Raises
    ------
    ValueError
        When the HEALPix cell indexes given have values out of :math:`[0, 12` x :math:`N_{side} ^ 2[`.

    Examples
    --------
    >>> from cdshealpix.ring import vertices_skycoord
    >>> import numpy as np
    >>> ipix = np.array([42, 6, 10])
    >>> depth = 12
    >>> vertices = vertices_skycoord(ipix, 1 << depth)
    """
    lon, lat = vertices(ipix, nside, step)
    return SkyCoord(ra=lon, dec=lat, frame="icrs", unit="rad")
