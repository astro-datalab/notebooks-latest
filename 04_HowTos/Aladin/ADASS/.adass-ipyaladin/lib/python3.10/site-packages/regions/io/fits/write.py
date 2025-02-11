# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings
from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import QTable
from astropy.utils.exceptions import AstropyUserWarning

from regions.core import Region, Regions, SkyRegion
from regions.core.registry import RegionsRegistry
from regions.shapes import RegularPolygonPixelRegion

__all__ = []


@RegionsRegistry.register(Region, 'serialize', 'fits')
@RegionsRegistry.register(Regions, 'serialize', 'fits')
def _serialize_fits(regions):
    region_data = []
    for region in regions:
        if isinstance(region, SkyRegion):
            warnings.warn('Sky regions cannot be serialized to the FITS '
                          'region format, skipping.', AstropyUserWarning)
            continue
        regdata = _serialize_region_fits(region)
        if regdata is not None:
            region_data.append(regdata)

    if region_data:
        fits_table = _make_table(region_data)
        return fits_table
    return QTable()


@RegionsRegistry.register(Region, 'write', 'fits')
@RegionsRegistry.register(Regions, 'write', 'fits')
def _write_fits(regions, filename, header=None, overwrite=False):
    """
    Convert a list of `~regions.Region` to a FITS region table and write
    to a file.

    Parameters
    ----------
    regions : list
        A list of `~regions.Region` objects.

    filename : str
        The filename in which the table is to be written.

    header : `~astropy.io.fits.Header`, optional
        The FITS header.

    overwrite : bool, optional
        If True, overwrite the output file if it exists. Raises an
        `OSError` if False and the output file exists. Default is False.
    """
    output = _serialize_fits(regions)

    if header is None:
        hdudoc = ('ASC-FITS-REGION-1.2: Rots, McDowell: FITS REGION '
                  'Binary Table Design')
        header = dict([('EXTNAME', 'REGION'),
                       ('EXTVER', 1),
                       ('EXTLEVEL', 1),
                       ('HDUNAME', 'REGION'),
                       ('HDUCLASS', 'ASC'),
                       ('HDUCLAS1', 'REGION'),
                       ('HDUCLAS2', 'STANDARD'),
                       ('HDUVERS', '1.2.0'),
                       ('HDUDOC', hdudoc),
                       ('CONTENT', 'REGION'),
                       ('ORIGIN', 'astropy/regions')])

    bin_table = fits.BinTableHDU(data=output, header=header)
    bin_table.writeto(filename, overwrite=overwrite)


@dataclass
class _RegionData:
    """
    Class to hold region data.
    """

    shape: str
    x: np.ndarray
    y: np.ndarray
    r: np.ndarray
    rotang: np.ndarray
    component: int


def _serialize_region_fits(region):
    if isinstance(region, RegularPolygonPixelRegion):
        region = region.to_polygon()

    region_clsname = region.__class__.__name__

    unsupported_regions = ('RectangleAnnulusPixelRegion',
                           'LinePixelRegion', 'TextPixelRegion',
                           'CompoundPixelRegion')
    if region_clsname in unsupported_regions:
        warnings.warn(f'({region_clsname} cannot be serialized using the '
                      'FITS format, skipping.', AstropyUserWarning)
        return None

    # translate region class to FITS shape name
    shape = region_clsname.lower().replace('pixelregion', '')

    if region.meta.get('include', None) == 0:
        shape = f'!{shape}'

    region_map = {'circleannulus': 'annulus',
                  'ellipseannulus': 'elliptannulus',
                  'rectangle': 'rotbox'}
    if shape in region_map:
        shape = region_map[shape]

    shape_params = []
    rotang = None
    for param in region._params:
        value = getattr(region, param)
        if param in ('center', 'vertices'):
            x, y = value.xy
        elif param == 'angle':
            rotang = value
        else:
            # ellipse region is defined by full axis lengths, but
            # FITS regions file uses semi-axis lengths
            if shape == 'ellipse':
                value /= 2.0
            shape_params.append(value)

    if not shape_params:
        shape_params = 0
    if rotang is None:
        rotang = u.Quantity(0, 'deg')

    component = region.meta.get('component', None)

    return _RegionData(shape, np.atleast_1d(x), np.atleast_1d(y),
                       np.atleast_1d(shape_params), np.atleast_1d(rotang),
                       component)


def _make_column(arrays):
    arr_sizes = [arr.size for arr in arrays]
    arr_size = np.max(arr_sizes)

    data = []
    for (arr, size) in zip(arrays, arr_sizes, strict=True):
        pad_width = arr_size - size
        if pad_width != 0:
            arr = np.pad(arr, (0, pad_width), mode='constant')
        if not isinstance(arr[0], u.Quantity):
            arr <<= u.pix
        if arr_size == 1 and pad_width == 0:
            arr = arr[0]  # use scalars instead of len-1 arrays
        data.append(arr)

    data = u.Quantity(data)

    return data


def _define_components(region_data):
    components = np.array([regdata.component for regdata in region_data])

    none_idx = [i for i, item in enumerate(components) if item is None]
    if not none_idx:
        return components

    comps = [i for i in components if i is not None]
    if comps:
        start_component = np.max(comps) + 1
        components[none_idx] = np.arange(len(none_idx)) + start_component
    else:
        # all components are set to None - do not write a COMPONENT
        # column
        components = None

    return components


def _make_table(region_data):
    tbl = QTable()
    attrs = ('shape', 'x', 'y', 'r', 'rotang')

    for attr in attrs:
        arrays = [getattr(data, attr) for data in region_data]
        if attr not in ('shape',):
            arrays = _make_column(arrays)
        tbl[attr.upper()] = arrays

    components = _define_components(region_data)
    if components is not None:
        tbl['COMPONENT'] = components

    return tbl
