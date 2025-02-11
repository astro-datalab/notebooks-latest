# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
The module provides several custom descriptor classes for attribute
validation of region classes.
"""

import abc

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.units import Quantity

from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord

__all__ = []


class RegionAttribute(abc.ABC):
    """
    Base descriptor class for region attribute validation.

    Parameters
    ----------
    doc : str, optional
        The documentation string for the attribute.
    """

    def __init__(self, doc=''):
        self.__doc__ = doc

    def __set_name__(self, owner, name):
        self.name = name

    def __get__(self, instance, owner):
        if instance is None:
            return self  # pragma: no cover
        return instance.__dict__[self.name]

    def __set__(self, instance, value):
        self._validate(value)
        instance.__dict__[self.name] = value

    def __delete__(self, instance):
        raise AttributeError(f'cannot delete {self.name!r}')

    @abc.abstractmethod
    def _validate(self, value):
        """
        Validate the attribute value.

        An exception is raised if the value is invalid.
        """
        raise NotImplementedError  # pragma: no cover


class ScalarPixCoord(RegionAttribute):
    """
    Descriptor class to check that value is a scalar
    `~regions.PixCoord`.
    """

    def _validate(self, value):
        if not (isinstance(value, PixCoord) and value.isscalar):
            raise ValueError(f'{self.name!r} must be a scalar PixCoord')


class OneDPixCoord(RegionAttribute):
    """
    Descriptor class to check that value is a 1D `~regions.PixCoord`.
    """

    def _validate(self, value):
        if not (isinstance(value, PixCoord) and not value.isscalar
                and value.x.ndim == 1):
            raise ValueError(f'{self.name!r} must be a 1D PixCoord')


class PositiveScalar(RegionAttribute):
    """
    Descriptor class to check that value is a strictly positive (> 0)
    scalar float/int (not `~astropy.units.Quantity`).
    """

    def _validate(self, value):
        if isinstance(value, Quantity):
            raise ValueError(f'{self.name!r} must be a scalar integer or '
                             'float')

        if not np.isscalar(value) or value <= 0:
            raise ValueError(f'{self.name!r} must be a strictly positive '
                             'scalar')


class ScalarSkyCoord(RegionAttribute):
    """
    Descriptor class to check that value is a scalar
    `~astropy.coordinates.SkyCoord`.
    """

    def _validate(self, value):
        if not (isinstance(value, SkyCoord) and value.isscalar):
            raise ValueError(f'{self.name!r} must be a scalar SkyCoord')


class OneDSkyCoord(RegionAttribute):
    """
    Descriptor class to check that value is a 1D
    `~astropy.coordinates.SkyCoord`.
    """

    def _validate(self, value):
        if not (isinstance(value, SkyCoord) and value.ndim == 1):
            raise ValueError(f'{self.name!r} must be a 1D SkyCoord')


class ScalarAngle(RegionAttribute):
    """
    Descriptor class to check that value is a scalar angle, either an
    `~astropy.coordinates.Angle` or `~astropy.units.Quantity` with
    angular units.
    """

    def _validate(self, value):
        if isinstance(value, Quantity):
            if not value.isscalar:
                raise ValueError(f'{self.name!r} must be a scalar')

            if not value.unit.physical_type == 'angle':
                raise ValueError(f'{self.name!r} must have angular units')
        else:
            raise ValueError(f'{self.name!r} must be a scalar angle')


class PositiveScalarAngle(RegionAttribute):
    """
    Descriptor class to check that value is a strictly positive scalar
    angle, either an `~astropy.coordinates.Angle` or
    `~astropy.units.Quantity` with angular units.
    """

    def _validate(self, value):
        if isinstance(value, Quantity):
            if not value.isscalar:
                raise ValueError(f'{self.name!r} must be a scalar')

            if not value.unit.physical_type == 'angle':
                raise ValueError(f'{self.name!r} must have angular units')

            if not value > 0:
                raise ValueError(f'{self.name!r} must be strictly positive')
        else:
            raise ValueError(f'{self.name!r} must be a strictly positive '
                             'scalar angle')


class RegionType(RegionAttribute):
    """
    Descriptor class to check the region type of value.

    Parameters
    ----------
    name : str
        The name of the attribute.

    regionclass : `~regions.Region`
        The region class to check.
    """

    def __init__(self, name, regionclass):
        super().__init__(name)
        self.regionclass = regionclass

    def _validate(self, value):
        if not isinstance(value, self.regionclass):
            raise ValueError(f'{self.name!r} must be a '
                             f'{self.regionclass.__name__} object')


class RegionMetaDescr(RegionAttribute):
    """
    Descriptor class for the region meta dictionary.

    If input as a pure `dict`, it will be converted to a `RegionMeta`
    object.
    """

    def __set__(self, instance, value):
        # RegionMeta subclasses dict
        if isinstance(value, dict) and not isinstance(value, RegionMeta):
            value = RegionMeta(value)
        super().__set__(instance, value)

    def _validate(self, value):
        if not isinstance(value, RegionMeta):
            raise ValueError(f'{self.name!r} must be a dict or RegionMeta '
                             'object')


class RegionVisualDescr(RegionAttribute):
    """
    Descriptor class for the region visual dictionary.

    If input as a pure `dict`, it will be converted to a `RegionVisual`
    object.
    """

    def __set__(self, instance, value):
        # RegionVisual subclasses dict
        if isinstance(value, dict) and not isinstance(value, RegionVisual):
            value = RegionVisual(value)
        super().__set__(instance, value)

    def _validate(self, value):
        if not isinstance(value, RegionVisual):
            raise ValueError(f'{self.name!r} must be a dict or RegionVisual '
                             'object')
