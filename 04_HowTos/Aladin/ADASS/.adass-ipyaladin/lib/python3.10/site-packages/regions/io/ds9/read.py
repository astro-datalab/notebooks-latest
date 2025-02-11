# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re
import string
import warnings
from dataclasses import dataclass

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.utils.data import get_readable_fileobj
from astropy.utils.exceptions import AstropyUserWarning

from regions.core import PixCoord, RegionMeta, Regions, RegionVisual
from regions.core.registry import RegionsRegistry
from regions.io.ds9.core import (DS9ParserError, ds9_frame_map,
                                 ds9_params_template, ds9_shape_to_region,
                                 make_region_template)
from regions.io.ds9.meta import _split_raw_metadata, _translate_ds9_to_visual

__all__ = []


@RegionsRegistry.register(Regions, 'read', 'ds9')
def _read_ds9(filename, cache=False):
    """
    Read a DS9 region file in as a list of `~regions.Region` objects.

    Parameters
    ----------
    filename : str
        The filename of the file to access.

    cache : bool or 'update', optional
        Whether to cache the contents of remote URLs. If 'update', check
        the remote URL for a new version but store the result in the
        cache.

    Returns
    -------
    regions : `regions.Regions`
        A `Regions` object containing a list of `~regions.Region`
        objects.
    """
    with get_readable_fileobj(filename, cache=cache) as fh:
        region_string = fh.read()
        return _parse_ds9(region_string)


@RegionsRegistry.register(Regions, 'parse', 'ds9')
def _parse_ds9(region_str):
    """
    Parse a DS9 region string.

    Parameters
    ----------
    region_str : str
        The string contents of a DS9 region file.

    Returns
    -------
    regions : `regions.Regions`
        A `Regions` object containing a list of `~regions.Region`
        objects.
    """
    # first parse the input string to generate the raw region data
    region_data = _parse_raw_data(region_str)

    # now parse the raw region data into region object(s)
    regions = []
    for region_data_ in region_data:
        region = _make_region(region_data_)
        if region is not None:  # skip region if error during parsing
            regions.extend(region)
    return Regions(regions)


@dataclass
class _RegionData:
    """
    Class to hold data used to initialize a Region object.

    Data for multi-annulus regions is stored in a single object.
    """

    frame: str
    region_type: str
    shape: str
    shape_params: str
    raw_meta: dict
    region_str: str


def _split_lines(region_str):
    """
    Split a region string on newlines and all semicolons, except those
    between {} braces ({} contains ds9 text strings).

    Parameters
    ----------
    region_str : str
        The string contents of a DS9 region file.

    Returns
    -------
    lines : list of str
        A list of strings.
    """
    return [line_.strip() for line in region_str.split('\n')
            for line_ in _split_semicolon(line)]


def _parse_raw_data(region_str):
    """
    Parse the raw data to extract the region data.

    Parameters
    ----------
    region_str : str
        The string contents of a DS9 region file.

    Returns
    -------
    lines : list of `_RegionData`
        A list of `_RegionData` objects. Note that at this stage, region
        data for multi-annulus regions is stored in a single object.
    """
    global_meta = {}
    composite_meta = ''
    frame = None
    region_data = []

    supported_frames = ['image', 'icrs', 'fk5', 'j2000', 'fk4', 'b1950',
                        'galactic', 'ecliptic']
    unsupported_frames = ['linear', 'amplifier', 'detector', 'physical',
                          'tile']
    wcs_frames = ['wcs', 'wcs0'] + [f'wcs{letter}'
                                    for letter in string.ascii_lowercase]
    unsupported_frames += wcs_frames

    supported_shapes = ['circle', 'ellipse', 'box', 'annulus', 'polygon',
                        'line', 'point', 'text', 'composite']
    unsupported_shapes = ['vector', 'ruler', 'compass', 'projection',
                          'panda', 'epanda', 'bpanda']

    supported_frames_shapes = supported_frames + supported_shapes
    unsupported_frames_shapes = unsupported_frames + unsupported_shapes
    valid_frames_shapes = supported_frames_shapes + unsupported_frames_shapes

    regex_frame_or_shape = re.compile('^#? *([+-]?)([a-zA-Z0-9]+)')

    for line in _split_lines(region_str):  # split on semicolons & newlines
        # skip blank lines
        if not line:
            continue

        # skip comments
        if (line.startswith('#')
                and not line.startswith(('# text(', '# composite('))):
            continue

        original_line = line  # used to parse text and tag fields (keep case)
        line = line.lower()

        # ds9 region files can have multiple (including successive)
        # global lines
        if line.startswith('global'):
            global_meta.update(_parse_metadata(original_line[7:]))
            continue

        match = regex_frame_or_shape.search(line)
        if match is None:
            raise ValueError(f'Unable to parse line "{line}".')
        include_symbol = match.groups()[0]
        frame_or_shape = match.groups()[1]

        if frame_or_shape not in valid_frames_shapes:
            warnings.warn(f'"{frame_or_shape}" frame or shape is not a valid '
                          'frame or region shape; unable to parse line '
                          f'"{line}", skipping.', AstropyUserWarning)
            continue

        if frame_or_shape in unsupported_frames_shapes:
            warnings.warn(f'"{frame_or_shape}" frame or shape is not '
                          'supported by the regions package, skipping.',
                          AstropyUserWarning)
            if frame_or_shape in unsupported_frames:
                frame = None
            continue

        if frame_or_shape in supported_frames:
            # NOTE: frame value persists for subsequent regions until changed
            frame = frame_or_shape
            continue

        if frame_or_shape in supported_shapes:
            shape = frame_or_shape

            # a frame must be defined before any shape(s)
            if frame is None:
                warnings.warn('A coordinate frame was not found for region: '
                              f'"{line}", skipping.', AstropyUserWarning)
                continue

            if shape == 'composite':
                idx = line.find('||')
                if idx == -1:
                    raise ValueError('unable to parse line with composite '
                                     f'shape: "{line}"')
                # composite metadata applies to all regions within the
                # composite shape
                composite_meta = _parse_metadata(line[idx + 2:].strip())
                # remove "composite=1" since we split the composite
                composite_meta.pop('composite', None)

            # NOTE: include=1/0 in metadata overrides the leading
            #       "-/+" symbol; -: include=0;, + or '': include=1
            include = 0 if include_symbol == '-' else 1
            include_meta = {'include': include}

            params_str, meta_str = _parse_shape_line(shape, original_line,
                                                     match.span())

            region_meta = _parse_metadata(meta_str)
            raw_meta = _define_raw_metadata(global_meta, composite_meta,
                                            include_meta, region_meta)

            region_type = 'sky'
            if frame == 'image':
                region_type = 'pixel'

            # composite shape is used only to extract metadata
            if shape != 'composite':
                region_data.append(_RegionData(frame, region_type, shape,
                                               params_str, raw_meta, line))

            # reset composite metadata after the composite region ends
            if '||' not in line and composite_meta:
                composite_meta = {}

    return region_data


def _parse_shape_line(shape, line, span):
    """
    Parse a line of a DS9 file containing a shape.

    Parameters
    ----------
    shape : str
        The DS9 region shape (e.g., 'circle').

    line : str
        A line defining a DS9 shape.

    span : 2 tuple
        A tuple containing the (start, end) positions of the match
        defining the shape.

    Returns
    -------
    shape_params_str : str
        The region shape parameters as a string.

    meta_str : str
        The region metadata as a string.
    """
    full_line = line
    line = full_line[span[1]:]  # starts with the shape parameters

    # ds9 writes out text regions in this odd (undocumented) format
    if shape == 'text' and full_line.lower().startswith('# text'):
        idx = line.find(' ')
        if idx == -1:
            raise ValueError(f'unable to parse line "{line}"')
        meta_str = line[idx + 1:]
        shape_params_str = line[:idx]

    else:
        # split line into shape parameters and metadata
        parts = line.split('#', 1)  # color value can contain a #
        if len(parts) == 1:
            shape_params_str = line
            meta_str = ''  # no metadata found
        else:
            shape_params_str, meta_str = parts

    # strip trailing space and | chars
    shape_params_str = shape_params_str.strip(' |')
    shape_params_str = re.sub('[()]', '', shape_params_str).lower()
    meta_str = meta_str.strip()

    return shape_params_str, meta_str


def _parse_metadata(metadata_str):
    """
    Parse metadata for a single DS9 region.

    Parameters
    ----------
    metadata_str : str
        The region metadata (e.g., everything after the "#").

    Returns
    -------
    metadata : dict
        The region metadata as a dictionary.
    """
    # {.*?}    # all chars in curly braces
    # \'.*?\'  # all chars in single quotes
    # \".*?\"  # all chars in double quotes
    # [-?\d+\.?\d*\s]+\s?  # ([-/+]floats [whitespace]) incl. repeats
    #                        (e.g., dashlist=8 3, width=3, textangle=18.35)
    # [^=\s]+\s*[-?\d+\.?\d*]*   # (all chars (e.g., point=diamond) or
    #                              (all chars [whitespace] digits)
    #                                 (e.g., point=diamond 42)
    metadata_regex = re.compile(r'([a-zA-Z]+)\s*=\s*({.*?}|\'.*?\'|\".*?\"|'
                                r'[-?\d+\.?\d*\s]+\s?|'
                                r'[^=\s]+\s*[-?\d+\.?\d*]*)')

    metadata = {}
    for key, val in metadata_regex.findall(metadata_str):
        key = key.lower()
        val = val.strip().strip("'").strip('"').lstrip('{').rstrip('}')
        if key not in metadata:
            if key == 'tag':
                val = [val]  # tag value is always a list
            metadata[key] = val
        elif key == 'tag':
            metadata[key].append(val)
        else:
            warnings.warn(f'Found duplicate metadata for "{key}", '
                          'skipping', AstropyUserWarning)
    return metadata


def _define_raw_metadata(global_meta, composite_meta, include_meta,
                         region_meta):
    """
    Define the raw metadata dictionary for the region.

    The raw metadata is not modified, except that invalid metadata is
    ignored.

    Parameters
    ----------
    global_meta : dict
        The global metadata.

    composite_meta : dict
        The metadata for a composite region.

    include_meta : dict
        The include/exclude metadata.

    region_meta : dict
        The region metadata.

    Returns
    -------
    meta : tuple of dict
        The (valid) raw metadata extracted from the region file.
    """
    all_meta = global_meta.copy()
    all_meta.update(composite_meta)
    all_meta.update(include_meta)
    # region_meta must come after include_meta because include=1/0 in
    # metadata overrides the leading "-/+" include symbol
    all_meta.update(region_meta)

    # valid DS9 point symbols
    valid_points = ('circle', 'box', 'diamond', 'cross', 'x', 'arrow',
                    'boxcircle')

    # valid DS9 line values
    valid_lines = ('0 0', '0 1', '1 0', '1 1')

    # binary keys must have values of 0 or 1
    # NOTE: textrotate is not documented in the DS9 region file spec
    binary_keys = ('dash', 'select', 'highlite', 'fixed', 'edit', 'move',
                   'rotate', 'delete', 'include', 'source', 'background',
                   'fill', 'vector', 'textrotate')

    metadata = {}
    for key, value in all_meta.items():
        try:
            value = float(value)
            if value.is_integer():
                value = int(value)
        except (ValueError, TypeError):
            pass

        is_invalid = False
        # point value can either be ["symbol int"] or ["symbol"]
        if key == 'point':
            val = value.split()
            if val[0] not in valid_points:
                is_invalid = True
            if len(val) == 2 and not float(val[1]).is_integer():
                is_invalid = True

        if key == 'line' and value not in valid_lines:
            is_invalid = True

        if key in binary_keys and value not in (0, 1):
            is_invalid = True

        if is_invalid:
            warnings.warn(f'DS9 "{key}={value}" is invalid and will be '
                          'ignored', AstropyUserWarning)
            continue

        metadata[key] = value

    return metadata


def _parse_pixel_coord(param_str):
    invalid_chars = ('d', 'r', 'p', ':', 'h', 'm', 's')
    for char in invalid_chars:
        if char in param_str:
            raise DS9ParserError('Cannot parse pixel region position '
                                 'coordinates')

    if param_str[-1] == 'i':
        param_str = param_str[:-1]

    # DS9 uses 1-indexed pixels
    return float(param_str) - 1


def _parse_sky_coord(param_str, frame, index):
    invalid_chars = ('i', 'p')
    for char in invalid_chars:
        if char in param_str:
            raise DS9ParserError('Cannot parse sky region position '
                                 'coordinates')

    if param_str[-1] == 'r':
        return Angle(param_str[:-1], unit=u.radian)

    elif param_str[-1] == 'd':
        return Angle(param_str[:-1], unit=u.degree)

    elif 'd' in param_str or 'h' in param_str:
        return Angle(param_str)

    elif ':' in param_str:
        if index % 2 == 0 and frame not in ('galactic', 'ecliptic'):
            return Angle(param_str, u.hourangle)
        else:
            return Angle(param_str, u.degree)

    else:  # unit not specified
        return Angle(float(param_str), unit=u.degree)


def _parse_coord(region_type, param_str, frame, index):
    if region_type == 'pixel':
        return _parse_pixel_coord(param_str)
    else:
        return _parse_sky_coord(param_str, frame, index)


def _parse_angle(param_str):
    invalid_chars = ('p', 'i')
    for char in invalid_chars:
        if char in param_str:
            raise DS9ParserError('Cannot parse sky region angle parameter')

    unit_mapping = {'"': u.arcsec,
                    "'": u.arcmin,
                    'd': u.deg,
                    'r': u.rad}
    if param_str[-1] not in string.digits:
        unit = unit_mapping[param_str[-1]]
        return u.Quantity(float(param_str[:-1]), unit=unit)
    else:
        return u.Quantity(float(param_str), unit=u.degree)


def _parse_size(region_type, param_str):
    if region_type == 'pixel':
        invalid_chars = ('"', "'", 'd', 'r', 'p')
        for char in invalid_chars:
            if char in param_str:
                raise DS9ParserError('Cannot parse pixel region size '
                                     'parameters - must not be in angular '
                                     'units')

        if param_str[-1] == 'i':
            param_str = param_str[:-1]

        return float(param_str)

    else:
        # size in angular units
        return _parse_angle(param_str)


def _parse_shape_params(region_data):
    """
    Parse the shape parameters for a region line.

    Parameters
    ----------
    region_data : `_RegionData` instance
        A `_RegionData` instance containing the data for a region line.
        Data for multi-annulus regions is stored in a single object.

    Returns
    -------
    shape_params : list of list(s)
        The shape parameters for each region(s). ``shape_params`` is
        a list of lists of shape parameters. A separate shape list is
        returned for multiple annulus/box/ellipse regions. Otherwise,
        the output list contains only one list.
    """
    region_type = region_data.region_type
    shape = region_data.shape
    frame = region_data.frame
    params = [val for val in re.split(r'\s|\,', region_data.shape_params)
              if val]  # split values on space or comma

    nparams = len(params)
    n_annulus = 0

    if shape in ('ellipse', 'box') and nparams > 5:
        if nparams % 2 != 1:
            raise ValueError(f'incorrect number of parameters ({nparams}) '
                             f'for shape "{shape}"')
        n_annulus = ((nparams - 3) // 2) - 1

    if shape in ('annulus',):
        n_annulus = nparams - 3

    if shape in ('ellipse', 'box', 'annulus'):
        # reset the cycle iterators
        shape_template = make_region_template()
    else:
        shape_template = ds9_params_template[shape]

    shape_params = []
    # TODO: check zip strict=True
    for idx, (param_type, value) in enumerate(zip(shape_template, params)):
        if shape in ('ellipse', 'box') and idx == nparams - 1:
            param_type = 'angle'  # last parameter is always an angle

        if param_type == 'coord':
            param = _parse_coord(region_type, value, frame, idx)
        elif param_type in ('length',):
            param = _parse_size(region_type, value)
            if shape == 'ellipse':
                param *= 2.0  # ds9 uses semi-axis lengths
        elif param_type in ('angle',):
            param = _parse_angle(value)
        else:
            raise ValueError('cannot parse shape parameters')

        shape_params.append(param)

    if n_annulus > 1:
        tmp_params = []
        for i in range(n_annulus):
            if shape in ('ellipse', 'box'):
                idx = (i + 1) * 2
                params = [shape_params[0], shape_params[1], shape_params[idx],
                          shape_params[idx + 1], shape_params[idx + 2],
                          shape_params[idx + 3], shape_params[-1]]
            else:
                idx = i + 2
                params = [shape_params[0], shape_params[1], shape_params[idx],
                          shape_params[idx + 1]]

            tmp_params.append(params)
        shape_params = tmp_params
    else:
        shape_params = [shape_params]

    if n_annulus > 0:
        if shape == 'ellipse':
            shape = 'ellipse_annulus'
        elif shape == 'box':
            shape = 'rectangle_annulus'

    return shape, shape_params


def _define_coords(region_type, params, frame=None):
    if region_type == 'pixel':
        coords = PixCoord(*params)
    else:
        coords = SkyCoord(*params, frame=frame)
    return coords


def _define_region_params(region_type, shape, shape_params, frame=None):
    if frame is not None:
        frame = ds9_frame_map[frame]

    if shape == 'polygon':
        coord_params = (shape_params[0::2], shape_params[1::2])
        params = [_define_coords(region_type, coord_params, frame=frame)]

    elif shape == 'line':
        params = [_define_coords(region_type, shape_params[0:2], frame=frame),
                  _define_coords(region_type, shape_params[2:4], frame=frame)]

    elif shape in ('ellipse_annulus', 'rectangle_annulus'):
        size_params = shape_params[2:-1]
        tmp = [size_params[0::2], size_params[1::2]]
        tmp_flat = [item for sublist in tmp for item in sublist]
        params = [_define_coords(region_type, shape_params[0:2], frame=frame),
                  *tmp_flat, shape_params[-1]]

    else:
        params = ([_define_coords(region_type, shape_params[0:2],
                                  frame=frame)] + shape_params[2:])

    return params


def _make_region(region_data):
    """
    Make a region object from the region data.

    Parameters
    ----------
    region_data : `_RegionData` instance
        A `_RegionData` instance.
    """
    try:
        # NOTE: returned shape can be different from region_data.shape
        shape, shape_params_list = _parse_shape_params(region_data)
    except DS9ParserError as err:
        # raise a warning and skip the region
        msg = f'{str(err)}: {region_data.region_str}'
        warnings.warn(msg, AstropyUserWarning)
        return None

    # define the parameters to initialize a Region
    # NOTE: region_params can be longer than region_data for
    # multi-annulus regions
    region_type = region_data.region_type
    region_params = []
    for shape_params in shape_params_list:
        region_params.extend([_define_region_params(region_type, shape,
                                                    shape_params,
                                                    region_data.frame)])

    # separate the metadata and visual metadata and then translate the
    # visual metadata to valid mpl kwargs for the particular region
    meta, visual = _split_raw_metadata(region_data.raw_meta)
    visual = _translate_ds9_to_visual(shape, visual)

    regions = []
    for shape_params in region_params:
        # for Text region, we need to add meta['text'] to params and
        # remove it from meta; set to '' if the text meta value was not
        # specified
        if shape == 'text':
            shape_params.append(region_data.raw_meta.get('text', ''))
            meta.pop('text', None)

        region = ds9_shape_to_region[region_type][shape](*shape_params)

        region.meta = RegionMeta(meta)
        region.visual = RegionVisual(visual)
        region._raw_meta = region_data.raw_meta

        regions.append(region)

    return regions


def _find_text_delim_idx(region_str):
    """
    Find the indices of the DS9 text field delimiters ({}, '', or "") in
    a string.
    """
    pattern = re.compile(r'(text\s*=\s*[{\'"])')
    idx0 = []
    delim = []
    start_idx = []
    for match in pattern.finditer(region_str):
        idx0.append(match.start())
        end_delim = match.group()[-1]
        if end_delim == '{':
            end_delim = '}'
        delim.append(end_delim)
        start_idx.append(match.span()[1])

    idx1 = []
    for sidx, char in zip(start_idx, delim, strict=True):
        idx1.append(region_str.find(char, sidx))

    return idx0, idx1


def _split_semicolon(region_str):
    r"""
    Split a DS9 region string on semicolons.

    The line is not split on semicolons found in a text field.

    This turned out to be a very tricky problem (attempts with regex
    failed):

      * text fields are delimited by {}, '', or ""
      * the text delimiters do not have to be consistent within a file
      * region strings can have unpaired ' (arcmin) or " (arcsec)
      * region strings can have "#" in a color value (cannot split on #)
      * the text field can contain "text={str}", "text='str'",
        'test="str"' etc., e.g., text={text="hello"}
      * escaped delimiter characters used by the text field
        are not allowed within the text field, e.g., "text={my
        field\{test\}}" and "text='my field \'test\''" are invalid.
        However, "text={my field, 'test'}" is valid

    This code finds the text field delimiters and then finds the indices
    of the opening and closing delimiters. Text fields that contain
    "text={str}", etc. are included (as a smaller range between the
    larger text field range), but this is fine because all we need are
    index ranges where to exclude semicolons (for splitting). Semicolons
    found at indices between the open/close delimiter indices are
    excluded from splitting.
    """
    idx0, idx1 = _find_text_delim_idx(region_str)

    semi_idx = [pos for pos, char in enumerate(region_str) if char == ';']
    fidx = []
    for i in semi_idx:
        for i0, i1 in zip(idx0, idx1, strict=True):
            if i0 <= i <= i1:
                break
        else:
            fidx.append(i + 1)
    fidx.insert(0, 0)

    return [region_str[i:j].rstrip(';')
            for i, j in zip(fidx, fidx[1:] + [None], strict=True)]
