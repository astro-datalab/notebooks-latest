# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import warnings
from copy import deepcopy

from astropy.coordinates import Angle, SkyCoord
from astropy.units import Quantity
from astropy.utils.exceptions import AstropyUserWarning

from regions.core import (CompoundPixelRegion, CompoundSkyRegion, PixCoord,
                          PixelRegion, Region, Regions)
from regions.core.registry import RegionsRegistry
from regions.io.ds9.core import ds9_frame_map, ds9_shape_templates
from regions.io.ds9.meta import _translate_metadata_to_ds9
from regions.shapes import RegularPolygonPixelRegion

__all__ = []


@RegionsRegistry.register(Region, 'serialize', 'ds9')
@RegionsRegistry.register(Regions, 'serialize', 'ds9')
def _serialize_ds9(regions, precision=8):
    if not regions:
        return ''

    region_data = []
    for region in regions:
        if isinstance(region, (CompoundPixelRegion, CompoundSkyRegion)):
            warnings.warn('Cannot serialize a compound region, skipping',
                          AstropyUserWarning)

        if isinstance(region, RegularPolygonPixelRegion):
            region = region.to_polygon()

        region_data.append(_serialize_region_ds9(region, precision=precision))

    # ds9 file header
    output = '# Region file format: DS9 astropy/regions\n'

    # extract common region metadata and place in the global metadata
    all_meta = []
    for region in region_data:
        region_meta = deepcopy(region['meta'])
        region_meta.pop('tag', None)  # "tag" cannot be in global metadata
        all_meta.append(region_meta)

    global_meta = dict(set.intersection(*[set(meta_dict.items())
                                          for meta_dict in all_meta]))
    if global_meta:
        output += f'global {_make_meta_str(global_meta)}\n'

    # define region frame and metadata (removing items that are
    # defined in global metadata)
    frames = []
    metadata = []
    for region in region_data:
        frames.append(region['frame'])
        meta_tmp = region['meta']
        for key in global_meta:
            meta_tmp.pop(key)
        metadata.append(meta_tmp)

    # extract the coordinate frame if it is identical for all regions
    # TODO: extract common coord frame(s) for block(s) of consecutive regions
    frames = set(frames)
    global_frame = None
    if len(frames) == 1:
        global_frame = frames.pop()
        output += f'{global_frame}\n'

    # add line for each region
    for region, region_meta in zip(region_data, metadata, strict=True):
        if global_frame is None:
            output += f'{region["frame"]}; '

        meta_str = _make_meta_str(region_meta)
        if meta_str:
            meta_str = f' # {meta_str}'

        output += f'{region["region"]}{meta_str}\n'

    return output


@RegionsRegistry.register(Region, 'write', 'ds9')
@RegionsRegistry.register(Regions, 'write', 'ds9')
def _write_ds9(regions, filename, *, precision=8, overwrite=False):
    """
    Convert a list of `~regions.Region` to a DS9 string and write to a
    file.

    Parameters
    ----------
    regions : list
        A list of `~regions.Region` objects.

    filename : str
        The output region filename.

    precision : int, optional
        The level of decimal precision given as the number of decimal
        places.

    overwrite : bool, optional
        If `True`, overwrite the output file if it exists. If `False`
        (default) and the output file exists, an `OSError` is raised.
    """
    if os.path.lexists(filename) and not overwrite:
        raise OSError(f'{filename} already exists')

    output = _serialize_ds9(regions, precision=precision)
    with open(filename, 'w') as fh:
        fh.write(output)


def _get_region_shape(region):
    shape = region.__class__.__name__.lower().replace('skyregion', '')
    shape = shape.replace('pixelregion', '')
    return shape


def _get_frame_name(region, mapping):
    if isinstance(region, PixelRegion):
        frame = 'image'
    elif 'center' in region._params:
        frame = region.center.frame.name
    elif 'vertices' in region._params:
        frame = region.vertices.frame.name
    elif 'start' in region._params:
        frame = region.start.frame.name
    else:
        raise ValueError(f'Unable to get coordinate frame for {region!r}')

    if frame not in mapping:
        warnings.warn(f'Cannot serialize region with frame={frame}, skipping',
                      AstropyUserWarning)

    return mapping[frame]


def _make_meta_str(meta):
    metalist = []
    for key, val in meta.items():
        if key == 'tag':  # can have multiple tags; value is always a list
            metalist.append(' '.join([f'tag={{{val}}}' for val in meta[key]]))
        else:
            metalist.append(f'{key}={val}')
    return ' '.join(metalist)


def _get_region_params(region, shape_template, precision=8):
    ellipse_axes = ('width', 'height', 'inner_width', 'inner_height',
                    'outer_width', 'outer_height')
    ellipse_names = ('ellipse', 'ellipseannulus')

    param = {}
    for param_name in region._params:
        if param_name in ('text',):
            continue

        value = getattr(region, param_name)

        # DS9 ellipse parameters are serialized as semi-axis lengths, but
        # ellipse region is defined by full axis lengths
        is_ellipse = (shape_template[0] in ellipse_names
                      and param_name in ellipse_axes)
        if not isinstance(value, (PixCoord, SkyCoord)) and is_ellipse:
            # deepcopy to prevent changing value in memory
            value = deepcopy(value) / 2.0  # semi-axis lengths

        if isinstance(value, PixCoord):
            # pixels; ds9's origin is (1, 1)
            if value.isscalar:
                value = (f'{value.x + 1:0.{precision}f},'
                         f'{value.y + 1:0.{precision}f}')
            else:
                value_str = ''
                for val in value:
                    value_str += (f'{val.x + 1:0.{precision}f},'
                                  f'{val.y + 1:0.{precision}f},')
                value = value_str[:-1]

        elif isinstance(value, SkyCoord):
            val = value.to_string(precision=precision)
            # polygon region has multiple SkyCoord
            value = ' '.join(val) if not value.isscalar else val
            value = value.replace(' ', ',')

        elif isinstance(value, Angle):
            value = value.to_string(unit='deg', decimal=True,
                                    precision=precision)

        elif isinstance(value, Quantity):
            # [:-4] to trim ' deg' from string end
            value = value.to_string(unit='deg', precision=precision)[:-4]

        else:
            value = f'{value:.{precision}f}'

        param[param_name] = value

    try:
        param_str = shape_template[1].format(**param)
    except KeyError as err:
        raise ValueError(
            f'Unable to get shape parameters for {region!r}') from err

    return param_str


def _serialize_region_ds9(region, precision=8):
    frame_mapping = {v: k for k, v in ds9_frame_map.items()}
    frame = _get_frame_name(region, mapping=frame_mapping)

    shape = _get_region_shape(region)
    if shape not in ds9_shape_templates:
        warnings.warn(f'Cannot serialize region shape "{shape}", '
                      'skipping', AstropyUserWarning)

    region_params = _get_region_params(region,
                                       ds9_shape_templates[shape],
                                       precision=precision)

    region_type = ds9_shape_templates[shape][0]
    region_str = f'{region_type}({region_params})'

    region_meta = _translate_metadata_to_ds9(region, shape)

    return {'frame': frame, 'region': region_str, 'meta': region_meta}
