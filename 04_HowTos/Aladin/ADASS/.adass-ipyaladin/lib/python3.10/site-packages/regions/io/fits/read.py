# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import QTable
from astropy.utils.exceptions import AstropyUserWarning

from regions.core import PixCoord, RegionMeta, Regions
from regions.core.registry import RegionsRegistry
from regions.io.fits.core import FITSParserError, shape_map

__all__ = []


@RegionsRegistry.register(Regions, 'parse', 'fits')
def _parse_fits(region_table):
    """
    Parse a FITS region table.

    Parameters
    ----------
    region_table : `~astropy.table.Table`
        The table contents of a FITS region file.

    Returns
    -------
    regions : `regions.Regions`
        A `Regions` object containing a list of `~regions.Region`
        objects.
    """
    regions = parse_table(region_table)
    return Regions(regions)


@RegionsRegistry.register(Regions, 'read', 'fits')
def _read_fits(filename, cache=False):
    """
    Read a FITS region file, converting a FITS regions table to a list
    of `~regions.Region` objects.

    Parameters
    ----------
    filename : str
        The FITS region filename. The first "REGION" FITS extension will
        be used.

    cache : bool or 'update', optional
        Whether to cache the contents of remote URLs. If 'update', check
        the remote URL for a new version but store the result in the
        cache.

    Returns
    -------
    regions : list
        A list of `~regions.Region` objects.
    """
    with fits.open(filename, cache=cache) as hdul:
        for hdu in hdul:
            # use the first 'REGION' HDU
            if hdu.name == 'REGION':
                region_table = QTable.read(hdu)
                regions = _parse_fits(region_table)
                return regions

        raise FITSParserError('An extension with name (EXTNAME) "REGION" '
                              'was not found')


def get_shape(region_row):
    include = 1
    shape_key = 'SHAPE'
    if shape_key not in region_row.colnames:
        shape = 'point'
        return shape, include

    shape = region_row[shape_key].lower()
    if shape[0] == '!':
        include = 0
        shape = shape[1:]

    supported_shapes = list(shape_map.keys())
    unsupported_shapes = ['pie', 'sector', 'diamond', 'rhombus',
                          'rotdiamond', 'rotrhombus']
    valid_shapes = supported_shapes + unsupported_shapes

    if shape not in valid_shapes:
        raise FITSParserError(f'{shape!r} is not a valid FITS region shape')
    if shape not in supported_shapes:
        warnings.warn(f'{shape!r} is not supported by the regions package, '
                      'skipping.', AstropyUserWarning)
        shape = None

    return shape, include


def get_column_values(region_row, colname):
    index = None
    if colname[-1].isdigit():
        index = int(colname[-1])
        colname = colname[:-1]

    value = np.atleast_1d(region_row[colname])
    if isinstance(value, u.Quantity) and value.unit == u.pixel:
        value = value.value  # remove pixel units

    if index is None:  # polygon uses all values
        return value

    try:
        return value[index]
    except IndexError as exc:
        raise FITSParserError(f'The {colname!r} column must have more '
                              f'than {index!r} values for the '
                              'region.') from exc


def get_shape_params(shape, region_row, shape_columns):
    values = [get_column_values(region_row, column)
              for column in shape_columns]

    if 'rectangle' in shape:
        (xmin, xmax, ymin, ymax) = values[0:4]
        xcenter = 0.5 * (xmin + xmax)
        ycenter = 0.5 * (ymin + ymax)
        xsize = xmax - xmin
        ysize = ymax - ymin

        shape_params = [PixCoord(xcenter, ycenter), xsize, ysize]
        if shape == 'rotrectangle':
            shape_params.append(values[-1])  # angle

        return shape_params

    # center (or polygon) coordinates for all other regions
    shape_params = [PixCoord(values[0], values[1])]

    # shape params
    if shape == 'ellipse':
        # FITS uses semi-axis lengths;
        # the last value is always the rotation angle
        values[2:-1] = list(np.array(values[2:-1]) * 2.)

    shape_params.extend(values[2:])

    return shape_params


def parse_row(region_row):
    shape, include = get_shape(region_row)
    if shape is None:
        return None

    region_cls, shape_columns = shape_map[shape]

    for column in shape_columns:
        if column[-1].isdigit():
            column = column[:-1]
        if column not in region_row.colnames:
            warnings.warn(f'Table columns are missing for {shape!r} shape, '
                          'skipping.', AstropyUserWarning)
            return None

    shape_params = get_shape_params(shape, region_row, shape_columns)
    region = region_cls(*shape_params)

    meta = {}
    if include == 0:
        meta = {'include': include}
        region.meta = RegionMeta(meta)

    shape_key = 'COMPONENT'
    if shape_key in region_row.colnames:
        component = int(region_row[shape_key])
        meta = {'component': component}

    if meta:
        region.meta = RegionMeta(meta)

    return region


def parse_table(region_table):
    valid_columns = ('X', 'Y', 'SHAPE', 'R', 'ROTANG', 'COMPONENT')

    for column in region_table.colnames:
        if column not in valid_columns:
            raise FITSParserError(f'{column!r} is not a valid column name')

    regions = []
    for row in region_table:
        region = parse_row(row)
        if region is not None:
            regions.append(region)

    return regions
