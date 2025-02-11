# Licensed under a 3-clause BSD style license - see LICENSE.rst
import copy

import numpy as np
from astropy.coordinates import SkyCoord

__all__ = ['PixCoord']


# Define global variables for the default 'origin' and 'mode' used for
# WCS transformations.
_DEFAULT_WCS_ORIGIN = 0
_DEFAULT_WCS_MODE = 'all'


class PixCoord:
    """
    A class for pixel coordinates.

    This class can represent a scalar or an array of pixel coordinates.
    `~regions.PixCoord` objects can be added or subtracted to each
    other. They can also be compared for equality.

    The data members are either numbers or `~numpy.ndarray`
    (not `~astropy.units.Quantity` objects with unit "pixel").

    Given a `astropy.wcs.WCS` object, it can be transformed to and from
    a `~astropy.coordinates.SkyCoord` object.

    Parameters
    ----------
    x : float or array-like
        Pixel coordinate x value.
    y : float or array-like
        Pixel coordinate y value.

    Examples
    --------
    Usage examples are provided in the :ref:`getting_started-coord`
    section of the documentation.
    """

    def __init__(self, x, y):
        x, y = np.broadcast_arrays(x, y)

        if x.shape == ():
            self.x, self.y = x.item(), y.item()
        else:
            self.x, self.y = x, y

    def copy(self):
        return self.__class__(copy.deepcopy(self.x), copy.deepcopy(self.y))

    @staticmethod
    def _validate(obj, name, expected='any'):
        """
        Validate that a given object is a valid `PixCoord`.

        Parameters
        ----------
        obj : `PixCoord`
            The object to check.
        name : str
            The parameter name used for error messages.
        expected : {'any', 'scalar', 'array'}
            What kind of PixCoord to check for.

        Returns
        -------
        obj : `PixCoord`
            The input object, if valid.
        """
        if not isinstance(obj, PixCoord):
            raise TypeError(f'{name!r} must be a PixCoord')

        if expected == 'any':
            pass
        elif expected == 'scalar':
            if not obj.isscalar:
                raise ValueError(f'{name!r} must be a scalar PixCoord')
        elif expected == 'array':
            if obj.isscalar:
                raise ValueError(f'{name!r} must be a PixCoord array')
        else:
            raise ValueError(f'Invalid value for "expected": {expected!r}')

        return obj

    @property
    def isscalar(self):
        """
        Whether the instance is scalar (e.g., a single (x, y)
        coordinate).
        """
        return np.isscalar(self.x)

    def __repr__(self):
        return f'{self.__class__.__name__}(x={self.x}, y={self.y})'

    def __len__(self):
        if self.isscalar:
            raise TypeError(f'Scalar {self.__class__.__name__!r} object has '
                            'no len()')
        return len(self.x)

    def __iter__(self):
        for (x, y) in zip(self.x, self.y, strict=True):
            yield PixCoord(x=x, y=y)

    def __getitem__(self, key):
        if self.isscalar:
            raise IndexError(f'Scalar {self.__class__.__name__!r} cannot be '
                             'indexed or sliced.')

        # Let Numpy do the slicing
        x = self.x[key]
        y = self.y[key]
        return PixCoord(x=x, y=y)

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError('Can add only to another PixCoord')
        return self.__class__(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError('Can subtract only from another PixCoord')
        return self.__class__(self.x - other.x, self.y - other.y)

    def __eq__(self, other):
        """
        Check whether ``other`` is `PixCoord` object and whether their
        abscissa and ordinate values are equal using
        `np.testing.assert_allclose` with its default tolerance values.
        """
        if isinstance(other, self.__class__):
            return np.allclose([self.x, self.y], [other.x, other.y])
        return False

    def to_sky(self, wcs, origin=_DEFAULT_WCS_ORIGIN, mode=_DEFAULT_WCS_MODE):
        """
        Convert to a `~astropy.coordinates.SkyCoord`.

        Parameters
        ----------
        wcs : `~astropy.wcs.WCS`
            The WCS to use to convert pixels to world coordinates.
        origin : int, optional
            Whether to return 0 or 1-based pixel coordinates.
        mode : {'all', 'wcs'}, optional
            Whether to do the transformation including distortions
            (``'all'``) or only including only the core WCS
            transformation (``'wcs'``).

        Returns
        -------
        coord : `~astropy.coordinates.SkyCoord`
            A new object with sky coordinates corresponding to the pixel
            coordinates.
        """
        return SkyCoord.from_pixel(xp=self.x, yp=self.y, wcs=wcs,
                                   origin=origin, mode=mode)

    @classmethod
    def from_sky(cls, skycoord, wcs, origin=_DEFAULT_WCS_ORIGIN,
                 mode=_DEFAULT_WCS_MODE):
        """
        Create `PixCoord` from a `~astropy.coordinates.SkyCoord`.

        Parameters
        ----------
        skycoord : `~astropy.coordinates.SkyCoord`
            The sky coordinate.
        wcs : `~astropy.wcs.WCS`
            The WCS to use to convert pixels to world coordinates.
        origin : int, optional
            Whether to return 0 or 1-based pixel coordinates.
        mode : {'all', 'wcs'}, optional
            Whether to do the transformation including distortions
            (``'all'``) or only including only the core WCS
            transformation (``'wcs'``).

        Returns
        -------
        coord : `PixCoord`
            A new `PixCoord` object at the position of the input sky
            coordinates.
        """
        x, y = skycoord.to_pixel(wcs=wcs, origin=origin, mode=mode)
        return cls(x=x, y=y)

    def separation(self, other):
        r"""
        Calculate the separation to another pixel coordinate.

        This is the two-dimensional Cartesian separation :math:`d` where

        .. math::
            d = \sqrt{(x_1 - x_2) ^ 2 + (y_1 - y_2) ^ 2}

        Parameters
        ----------
        other : `PixCoord`
            The other pixel coordinate.

        Returns
        -------
        separation : `numpy.array`
            The separation in pixels.
        """
        dx = other.x - self.x
        dy = other.y - self.y
        return np.hypot(dx, dy)

    @property
    def xy(self):
        """
        A 2-tuple ``(x, y)`` for this coordinate.
        """
        return self.x, self.y

    def rotate(self, center, angle):
        """
        Rotate the pixel coordinate.

        Positive ``angle`` corresponds to counter-clockwise rotation.

        Parameters
        ----------
        center : `PixCoord`
            The rotation center point.
        angle : `~astropy.coordinates.Angle`
            The rotation angle.

        Returns
        -------
        coord : `PixCoord`
            The rotated coordinates (which is an independent copy).
        """
        dx = self.x - center.x
        dy = self.y - center.y
        vec = np.array([dx, dy])

        cosa, sina = np.cos(angle), np.sin(angle)
        rotation_matrix = np.array([[cosa, -sina], [sina, cosa]])

        vec = np.matmul(rotation_matrix, vec)

        return self.__class__(center.x + vec[0], center.y + vec[1])
