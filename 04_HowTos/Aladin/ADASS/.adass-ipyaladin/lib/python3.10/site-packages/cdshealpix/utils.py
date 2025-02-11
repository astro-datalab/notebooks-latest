"""Conversions between NESTED and RING scheme."""
import functools

from astropy.coordinates import Latitude, Longitude
import numpy as np

from . import cdshealpix

__all__ = ["to_ring", "from_ring"]

# Raise a ValueError exception if the input
# HEALPix cells array contains invalid values


def _check_ipixels(data, depth):
    npix = 12 * 4 ** (depth)
    if (data >= npix).any() or (data < 0).any():
        raise ValueError(
            f"The input HEALPix cells contains value out of [0, {npix - 1}]"
        )


def _validate_lonlat(function):
    """Validate the longitude and latitudes entries of methods of the MOC class.

    Parameters
    ----------
    function : <class 'function'>
        must have the signature function(lon, lat, *args, **kwargs)

    Returns
    -------
        applies desired transformations for the `lon` and `lat` arguments and calls
        `function` with these modified arguments
    """

    @functools.wraps(function)
    def _validate_lonlat_wrap(lon, lat, *args, **kwargs):
        # be sure that lon and lat are of the same shape
        if lon.shape != lat.shape:
            raise ValueError(
                f"'lon' and 'lat' should have the same shape but are of shapes {lon.shape} and {lat.shape}",
            )
        # convert into astropy objects
        lon = lon if isinstance(lon, Longitude) else Longitude(lon)
        lat = lat if isinstance(lat, Latitude) else Latitude(lat)
        return function(lon, lat, *args, **kwargs)

    return _validate_lonlat_wrap


def _check_depth(depth):
    ravel_depth = np.ravel(np.atleast_1d(depth))
    if any(ravel_depth < 0) or any(ravel_depth > 29):
        raise ValueError("Depth must be in the [0, 29] closed range")


def to_ring(ipix, depth, num_threads=0):
    """Convert HEALPix cells from the NESTED to the RING scheme.

    Parameters
    ----------
    ipix : `numpy.ndarray`
        The HEALPix cell indexes in the NESTED scheme.
    depth : int
        The depth of the HEALPix cells.
    num_threads : int, optional
        Specifies the number of threads to use for the computation. Default to 0 means
        it will choose the number of threads based on the RAYON_NUM_THREADS environment variable (if set),
        or the number of logical CPUs (otherwise)

    Returns
    -------
    ipix_ring : `numpy.ndarray`
        The corresponding HEALPix cells in the RING scheme.

    Raises
    ------
    ValueError
        When the HEALPix cell indexes given have values out of :math:`[0, 4^{29 - depth}[`.

    Examples
    --------
    >>> from cdshealpix import to_ring
    >>> import numpy as np
    >>> ipix = np.array([42, 6, 10])
    >>> depth = 12
    >>> print(to_ring(ipix, depth))
    [100526076 100591616 100591614]
    """
    _check_depth(depth)

    ipix = np.atleast_1d(ipix)
    _check_ipixels(data=ipix, depth=depth)
    ipix = ipix.astype(np.uint64)

    # Allocation of the array containing the cells under the RING scheme
    ipix_ring = np.zeros(ipix.shape, dtype=np.uint64)

    num_threads = np.uint16(num_threads)
    cdshealpix.to_ring(depth, ipix, ipix_ring, num_threads)

    return ipix_ring


def from_ring(ipix, depth, num_threads=0):
    """Convert HEALPix cells from the RING to the NESTED scheme.

    Parameters
    ----------
    ipix : `numpy.ndarray`
        The HEALPix cell indexes in the RING scheme.
    depth : int
        The depth of the HEALPix cells.
    num_threads : int, optional
        Specifies the number of threads to use for the computation. Default to 0 means
        it will choose the number of threads based on the RAYON_NUM_THREADS environment variable (if set),
        or the number of logical CPUs (otherwise)

    Returns
    -------
    ipix_nested : `numpy.ndarray`
        The corresponding HEALPix cells in the NESTED scheme.

    Raises
    ------
    ValueError
        When the HEALPix cell indexes given have values out of :math:`[0, 4^{29 - depth}[`.

    Examples
    --------
    >>> from cdshealpix import from_ring
    >>> import numpy as np
    >>> ipix = np.array([42, 6, 10])
    >>> depth = 12
    >>> print(from_ring(ipix, depth))
    [16777203 33554430 67108862]
    """
    _check_depth(depth)

    ipix = np.atleast_1d(ipix)
    _check_ipixels(data=ipix, depth=depth)
    ipix = ipix.astype(np.uint64)

    # Allocation of the array containing the cells under the NESTED scheme
    ipix_nested = np.zeros(ipix.shape, dtype=np.uint64)

    num_threads = np.uint16(num_threads)
    cdshealpix.from_ring(depth, ipix, ipix_nested, num_threads)

    return ipix_nested
