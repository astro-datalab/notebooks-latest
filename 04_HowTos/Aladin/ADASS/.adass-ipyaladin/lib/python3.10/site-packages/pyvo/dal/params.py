# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains functionallity related to VOTABLE Params.
"""
import numpy as np
from collections.abc import MutableSet
import abc

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.units import Quantity, Unit
from astropy.time import Time
from astropy.io.votable.converters import (
    get_converter as get_votable_converter)

from .exceptions import DALServiceError


NUMERIC_DATATYPES = {'short', 'int', 'long', 'float', 'double'}


def find_param_by_keyword(keyword, params):
    """
    Searches for a specific param by keyword.
    This function will try to look for the keyword as-is first, and then tries
    to find the uppercase'd version of the keyword.
    """
    if keyword in params:
        return params[keyword]

    keyword = keyword.upper()
    if keyword in params:
        return params[keyword]

    raise KeyError('No param named {} defined'.format(keyword))


registry = dict()


def xtype(name):
    def decorate(cls):
        registry[name] = cls
        return cls
    return decorate


def unify_value(func):
    """
    Decorator for serialize method to do unit conversion on input value.
    The decorator converts the input value to the unit in the input param.
    """

    def wrapper(self, value):
        if self._param.unit:
            value = Quantity(value)
            if not value.unit.to_string():
                value = value * Unit(self._param.unit)
            else:
                value = value.to(self._param.unit)

        if isinstance(value, Quantity):
            value = value.value

        return func(self, value)

    return wrapper


def get_converter(param):
    if param.xtype in registry:
        return registry[param.xtype](param)

    if param.datatype in NUMERIC_DATATYPES:
        return Number(param)

    return Converter(param)


class Converter:
    """
    Base class for all converters. Each subclass handles the conversion of a
    input value based on a specific xtype.
    """

    def __init__(self, param):
        self._param = param

    def serialize(self, value):
        """
        Serialize for use in DAL Queries
        """
        if isinstance(value, list):
            # multiple values
            return [str(_) for _ in value]
        else:
            return str(value)


class Number(Converter):
    def __init__(self, param):
        if param.datatype not in {'short', 'int', 'long', 'float', 'double'}:
            pass

        super().__init__(param)

    @unify_value
    def serialize(self, value):
        """
        Serialize for use in DAL Queries
        """
        return get_votable_converter(self._param).output(
            value, np.zeros_like(value))


@xtype('timestamp')
class Timestamp(Converter):
    def __init__(self, param):
        if param.datatype != 'char':
            raise DALServiceError('Datatype is not char')

        super().__init__(param)

    def serialize(self, value):
        """
        Serialize time values for use in DAL Queries
        """
        value = Time(value)

        if value.size == 1:
            return value.isot
        else:
            raise DALServiceError('Expecting a scalar time value')


@xtype('interval')
class Interval(Number):
    def __init__(self, param):
        try:
            arraysize = int(param.arraysize)
            if arraysize % 2:
                raise DALServiceError('Arraysize is not even')
        except ValueError:
            raise DALServiceError('Arraysize is not even')

        super().__init__(param)

    @unify_value
    def serialize(self, value):
        size = np.size(value)
        if size % 2:
            raise DALServiceError('Interval size is not even')

        return super().serialize(value)


@xtype('point')
class Point(Number):
    def __init__(self, param):
        try:
            arraysize = int(param.arraysize)
            if arraysize != 2:
                raise DALServiceError('Point arraysize must be 2')
        except ValueError:
            raise DALServiceError('Point arraysize must be 2')

        super().__init__(param)

    @unify_value
    def serialize(self, value):
        size = np.size(value)
        if size != 2:
            raise DALServiceError('Point size must be 2')

        return super().serialize(value)


@xtype('circle')
class Circle(Number):
    def __init__(self, param):
        arraysize = int(param.arraysize)
        if arraysize != 3:
            raise DALServiceError('Circle arraysize must be 3')

        super().__init__(param)

    @unify_value
    def serialize(self, value):
        size = np.size(value)
        if size != 3:
            raise DALServiceError('Circle size must be 3')

        return super().serialize(value)


@xtype('polygon')
class Polygon(Number):
    def __init__(self, param):
        try:
            arraysize = int(param.arraysize)
            if arraysize % 3:
                raise DALServiceError('Arraysize is not a multiple of 3')
        except ValueError:
            if param.arraysize != '*':
                raise DALServiceError('Arraysize is not a multiple of 3')

        super().__init__(param)

    @unify_value
    def serialize(self, value):
        size = np.size(value)
        try:
            if size % 3:
                raise DALServiceError('Size is not a multiple of 3')
        except ValueError:
            raise DALServiceError('Size is not a multiple of 3')

        return super().serialize(value)


class AbstractDalQueryParam(MutableSet, metaclass=abc.ABCMeta):
    """
    Base class for a DAL parameter. In general, a DAL parameter allows
    for multiple values which are OR-ed by the service. As such, the class
    behaves like a set that stores all the values.

    When a value is added to an attribute, it is validated and formatted
    to conform to the using service (SIA2 or SODA) and value errors might be
    raised. The `dal` attribute stores the current list of formatted
    attributes.

    Subclasses must override the `dal_formatter` method that formats values
    for serialization. That includes unit conversions and string representation
    Duplicates in the set are determine based on the formatted DAL
    representation of the value.
    """

    def __init__(self, values=()):
        """
        Parameters
        ----------
        values : sequence, optional
            An initial set of values.
        """
        self.dal = []
        self._data = []
        for item in values:
            self.add(item)
        super().__init__()

    @abc.abstractmethod
    def get_dal_format(self, item):
        """
        Method to be provided by subclasses
        """
        return

    def add(self, item):
        if item in self:
            return
        self._data.append(item)
        self.dal.append(self.get_dal_format(item))

    def discard(self, item):
        # relies on the fact that both the raw and the formatted
        # attribute lists have the items in the same order. It
        # uses the formatted list (normalized units) to get the index.
        index = self.dal.index(self.get_dal_format(item))
        self._data.pop(index)
        self.dal.pop(index)

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __contains__(self, item):
        # check dal format for duplications since the quantities are known
        return self.get_dal_format(item) in self.dal


class StrQueryParam(AbstractDalQueryParam):
    """
    Representation of a unitless, single-value parameter. The formatter is
    just a str() cast
    """

    def get_dal_format(self, val):
        return str(val)


class PosQueryParam(AbstractDalQueryParam):
    """
    Representation of a position parameter. Depending on the number
    of entries, the resulting DAL format is CIRCLE, RANGE or POLYGON.
    """

    def get_dal_format(self, val):
        """
        formats the tuple values into a string to be sent to the service
        entries in values are either quantities or assumed to be degrees
        """
        self._validate_pos(val)
        if len(val) == 2 or len(val) == 3:
            shape = 'CIRCLE'
        elif len(val) == 4:
            shape = 'RANGE'
        elif len(val) > 5 and not len(val) % 2:
            shape = 'POLYGON'
        else:
            raise ValueError(
                'Invalid shape {}. Tuple with 3 (CIRCLE), 4 (RANGE) or '
                'even 6 and above (POLYGON) accepted.'.format(val))
        return '{} {}'.format(shape, ' '.join(
            [str(val.to(u.deg).value) if isinstance(val, Quantity) else
             val.transform_to('icrs').to_string() if isinstance(val, SkyCoord) else
             str((val * u.deg).value) for val in val]))

    def _validate_pos(self, pos):
        """
        validates position

        This has probably done already somewhere else
        """

        if len(pos) == 2:
            if not isinstance(pos[0], SkyCoord):
                raise ValueError
            if not isinstance(pos[1], Quantity):
                radius = pos[1] * u.deg
            else:
                radius = pos[1]
            if radius <= 0 * u.deg or radius.to(u.deg) > 90 * u.deg:
                raise ValueError('Invalid circle radius: {}'.format(radius))
        elif len(pos) == 3:
            self._validate_ra(pos[0])
            self._validate_dec(pos[1])
            if not isinstance(pos[2], Quantity):
                radius = pos[2] * u.deg
            else:
                radius = pos[2]
            if radius <= 0 * u.deg or radius.to(u.deg) > 90 * u.deg:
                raise ValueError('Invalid circle radius: {}'.format(radius))
        elif len(pos) == 4:
            ra_min = pos[0] if isinstance(pos[0], Quantity) else pos[0] * u.deg
            ra_max = pos[1] if isinstance(pos[1], Quantity) else pos[1] * u.deg
            dec_min = pos[2] if isinstance(pos[2], Quantity) \
                else pos[2] * u.deg
            dec_max = pos[3] if isinstance(pos[3], Quantity) \
                else pos[3] * u.deg
            self._validate_ra(ra_min)
            self._validate_ra(ra_max)
            if ra_max.to(u.deg) < ra_min.to(u.deg):
                raise ValueError('min > max in ra range: {} > {}'.
                                 format(ra_min, ra_max))
            self._validate_dec(dec_min)
            self._validate_dec(dec_max)
            if dec_max.to(u.deg) < dec_min.to(u.deg):
                raise ValueError('min > max in dec range: {} > {}'.
                                 format(dec_min, dec_max))
        else:
            for i, m in enumerate(pos):
                if i % 2:
                    self._validate_dec(m)
                else:
                    self._validate_ra(m)

    def _validate_ra(self, ra):
        if not isinstance(ra, Quantity):
            ra = ra * u.deg
        if ra.to(u.deg).value < 0 or ra.to(u.deg).value > 360.0:
            raise ValueError('Invalid ra: {}'.format(ra))

    def _validate_dec(self, dec):
        if not isinstance(dec, Quantity):
            dec = dec * u.deg
        if dec.to(u.deg).value < -90.0 or dec.to(u.deg).value > 90.0:
            raise ValueError('Invalid dec: {}'.format(dec))


class IntervalQueryParam(AbstractDalQueryParam):
    """
    Representation of an interval DAL parameter.
    """

    def __init__(self, unit=None, equivalencies=None):
        """
        Parameters
        ----------
        unit : `astropy.unit`
            Unit this paramter is represented in DAL format
        equivalencies: list
            List of equivalencies for unit conversion
        """
        self._unit = unit
        self._equivalencies = equivalencies
        super().__init__()

    def get_dal_format(self, val):
        if isinstance(val, (tuple, Quantity)):
            if len(val) == 1:
                high = low = val[0]
            elif len(val) == 2:
                low = val[0]
                high = val[1]
            else:
                raise ValueError('Too few/many values in interval attribute: {}'.
                                 format(val))
        else:
            high = low = val
        if isinstance(low, (int, float)) and isinstance(high, (int, float))\
                and low > high:
            raise ValueError('Invalid interval: min({}) > max({})'.format(
                low, high))
        if self._unit:
            if not isinstance(low, Quantity):
                low = low * self._unit
            low = low.to(self._unit, equivalencies=self._equivalencies).value
            if not isinstance(high, Quantity):
                high = high * self._unit
            high = high.to(self._unit, equivalencies=self._equivalencies).value

            if low > high:
                # interval could become invalid during transform (e.g. GHz->m)
                low, high = high, low

        return '{} {}'.format(low, high)


class TimeQueryParam(AbstractDalQueryParam):
    """
    Representation of a timestamp parameter.
    """

    def get_dal_format(self, val):
        if isinstance(val, tuple):
            if len(val) == 1:
                max_time = min_time = val[0]
            elif len(val) == 2:
                min_time = val[0]
                max_time = val[1]
            else:
                raise ValueError('Too few/many members in time attribute: {}'.
                                 format(val))
        else:
            max_time = min_time = val

        if not isinstance(min_time, Time):
            min_time = Time(min_time)
        if not isinstance(max_time, Time):
            max_time = Time(max_time)
        if min_time > max_time:
            raise ValueError('Invalid time interval: min({}) > max({})'.format(
                min_time, max_time
            ))
        return '{} {}'.format(min_time.mjd, max_time.mjd)


class EnumQueryParam(AbstractDalQueryParam):
    """
    Representation of an enum parameter
    """

    def __init__(self, allowed_values):
        """
        Parameters
        ----------
        allowed_values : sequence
            Sequence of allowed values for the enum
        """
        self._allowed = allowed_values
        super().__init__()

    def get_dal_format(self, val):
        if val not in self._allowed:
            raise ValueError('{} not a valid value from: {}'.
                             format(val, self._allowed))
        return str(val)
