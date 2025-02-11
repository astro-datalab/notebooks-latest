# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

from astropy.utils.exceptions import AstropyUserWarning

from regions.io.ds9.core import DS9ParserError, ds9_valid_symbols

__all__ = []


def _split_raw_metadata(raw_metadata):
    """
    Split the raw metadata into meta and visual dictionaries.

    The raw metadata is not changed, except unsupported metadata is
    ignored.

    Parameters
    ----------
    raw_metadata : dict
        The raw metadata as read from the DS9 region file/str.

    Returns
    -------
    meta, visual : dict
        The raw metadata and visual metadata dictionaries.
    """
    # not currently supported by regions
    unsupported_meta = ('line', 'vector', 'ruler', 'compass')

    # TODO: include text in visual keys to support text annotations for
    # all DS9 regions?
    ds9_visual_keys = ('color', 'dash', 'dashlist', 'fill', 'font', 'point',
                       'textangle', 'textrotate', 'width')

    meta = {}
    visual = {}
    for key, value in raw_metadata.items():
        if key in unsupported_meta:
            # don't warn for "line=0 0" (no arrows) because it works,
            # but skip adding it to metadata
            if key != 'line' or '1' in value:
                warnings.warn(f'DS9 meta "{key}={value}" is unsupported and '
                              'will be ignored', AstropyUserWarning)
            continue

        if key in ds9_visual_keys:
            visual[key] = value
        else:
            meta[key] = value

    # set the default plotting style to 'ds9' when parsing ds9 data
    visual['default_style'] = 'ds9'

    return meta, visual


def _translate_ds9_to_visual(shape, visual_meta):
    """
    Translate ds9 visual metadata to dictionary of matplotlib keywords.

    Parameters
    ----------
    shape : str
        The DS9 region shape (e.g., 'circle').

    visual_meta : dict
        The visual metadata.

    Returns
    -------
    visual : dict
        The visual dictionaries of region metadata translated
        to matplotlib keywords.
    """
    meta = visual_meta.copy()

    fill = meta.pop('fill', 0)
    # fill=1 is supported in DS9 only for the circle, ellipse, box, and
    # polygon shapes
    if fill == 1 and shape in ('circle', 'ellipse', 'box', 'polygon'):
        meta['fill'] = True

    dash = meta.pop('dash', 0)
    dashlist = meta.pop('dashlist', None)
    if int(dash) == 1:
        if dashlist is not None:
            dashes = tuple(int(i) for i in dashlist.split())
            meta['linestyle'] = (0, dashes)
        else:
            meta['linestyle'] = 'dashed'

        if shape == 'point':
            warnings.warn('dashed lines are unsupported for DS9 point '
                          'regions and will be ignored', AstropyUserWarning)
            meta.pop('linestyle', None)

    # "point=symbol [size]"; size is optional, symbol is not
    point = meta.pop('point', None)
    if point is not None:
        point_ = point.split()
        if len(point_) == 1:
            ds9_marker = point_[0]
        elif len(point_) == 2:
            ds9_marker, meta['markersize'] = point_
        else:
            raise ValueError(f'invalid point data "{point}"')
        meta['marker'] = ds9_valid_symbols[ds9_marker]

    font = meta.pop('font', None)
    if font is not None:
        vals = font.split()
        if len(vals) == 4:
            (meta['fontname'], meta['fontsize'], meta['fontweight'],
             meta['fontstyle']) = vals
        if len(vals) == 3:
            meta['fontname'], meta['fontsize'], meta['fontweight'] = vals
        if len(vals) == 2:
            meta['fontname'], meta['fontsize'] = vals
        if len(vals) == 1:
            meta['fontname'] = vals[0]

        # define default font values (helvetica 10 normal roman)
        if 'fontsize' not in meta:
            meta['fontsize'] = '10'
        if 'fontweight' not in meta:
            meta['fontweight'] = 'normal'
        if 'fontstyle' not in meta:
            meta['fontstyle'] = 'normal'

        # fontsize is a string
        try:
            meta['fontsize'] = int(meta['fontsize'])
        except ValueError:
            raise DS9ParserError('font size must be an integer, got '
                                 f'{meta["fontsize"]}') from None
        meta['fontstyle'] = meta['fontstyle'].replace('roman', 'normal')

    if shape == 'point':
        width = meta.pop('width', None)
        if width is not None:
            meta['markeredgewidth'] = width
    else:
        invalid = ('marker', 'markersize')
        for key in invalid:
            meta.pop(key, None)

    if shape == 'text':
        textangle = meta.pop('textangle', None)
        if textangle is not None:
            textrotate = meta.pop('textrotate', None)
            if textrotate != 0:  # rotate if None or 1
                meta['rotation'] = textangle

        # remove invalid mpl kwargs for matplotlib.text.Text
        invalid = ('linestyle', 'linewidth', 'fill')
        for key in invalid:
            meta.pop(key, None)
    else:
        # these kwarg are valid only for matplotlib.text.Text
        invalid = ('textangle', 'textrotate')
        for key in invalid:
            meta.pop(key, None)

    if shape not in ('point', 'line', 'text'):
        color = meta.pop('color', None)
        if color is not None:
            meta['facecolor'] = color
            meta['edgecolor'] = color

    return meta


def _remove_invalid_keys(region_meta, valid_keys):
    # TODO: instead of new dict, del region_meta in-place?
    meta = {}
    for key in region_meta:
        if key in valid_keys:
            meta[key] = region_meta[key]
    return meta


def _translate_metadata_to_ds9(region, shape):
    """
    Translate region metadata to valid ds9 meta keys.
    """
    meta = {**region.meta, **region.visual}
    # special case for Text regions
    if 'text' in region._params:
        meta = {'text': region.text, **meta}

    if 'annulus' in shape:
        # ds9 does not allow fill for annulus regions
        meta.pop('fill', None)

    fill = meta.pop('fill', None)
    if fill is not None:
        meta['fill'] = int(fill)

    if 'text' in meta:
        meta['text'] = f'{{{meta["text"]}}}'

    edgecolor = meta.pop('edgecolor', None)
    facecolor = meta.pop('facecolor', None)
    color = None
    if edgecolor is not None:
        color = edgecolor
    if facecolor is not None:
        if facecolor != color:
            warnings.warn('facecolor and edgecolor are different, edgecolor '
                          'will be used', AstropyUserWarning)
        if color is None:
            color = facecolor
    if color is not None:
        meta['color'] = color

    linewidth = meta.pop('linewidth', None)
    if linewidth is not None:
        meta['width'] = linewidth

    # point region marker width
    markeredgewidth = meta.pop('markeredgewidth', None)
    if markeredgewidth is not None:
        meta['width'] = markeredgewidth

    marker = meta.pop('marker', None)
    if marker is not None:
        symbol_map = {y: x for x, y in ds9_valid_symbols.items()}
        if marker in symbol_map:
            markersize = meta.pop('markersize', None)
            msize = ''
            if markersize is not None:
                msize = f' {markersize}'
            meta['point'] = f'{symbol_map[marker]}{msize}'
        else:
            warnings.warn(f'Unable to serialize marker "{marker}"',
                          AstropyUserWarning)

    fontname = meta.pop('fontname', None)
    if fontname is not None:
        fontsize = meta.pop('fontsize', 10)  # default 10
        fontweight = meta.pop('fontweight', 'normal')  # default normal
        # default roman
        fontstyle = meta.pop('fontstyle', 'roman').replace('normal', 'roman')
        meta['font'] = f'"{fontname} {fontsize} {fontweight} {fontstyle}"'

    linestyle = meta.pop('linestyle', None)
    if linestyle is not None:
        meta['dash'] = 1
    # if linestyle in ('dashed', '--'):
    if isinstance(linestyle, tuple):
        meta['dashlist'] = f'{linestyle[1][0]} {linestyle[1][1]}'

        # dashes = meta.pop('dashes', None)
        # if dashes is not None:
        #     meta['dashlist'] = f'{dashes[0]} {dashes[1]}'

    rotation = meta.pop('rotation', None)
    if rotation is not None:
        meta['textangle'] = rotation
        # meta['textrotate'] = 1

    # ds9 meta keys
    meta_keys = ['background', 'delete', 'edit', 'fixed', 'highlite',
                 'include', 'move', 'rotate', 'select', 'source', 'tag',
                 'text']
    visual_keys = ['color', 'dash', 'dashlist', 'fill', 'font', 'point',
                   'textangle', 'textrotate', 'width']
    valid_keys = meta_keys + visual_keys
    meta = _remove_invalid_keys(meta, valid_keys)

    return meta
