# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numbers
from warnings import warn

import astropy.units as u
from astropy.coordinates import (Angle, SkyCoord, UnitSphericalRepresentation,
                                 frame_transform_graph)
from astropy.utils.exceptions import AstropyUserWarning

from regions.core.core import PixCoord, SkyRegion
from regions.core.metadata import RegionMeta, RegionVisual
from regions.io.crtf.core import CRTFRegionParserWarning
from regions.shapes import (CircleAnnulusPixelRegion, CircleAnnulusSkyRegion,
                            CirclePixelRegion, CircleSkyRegion,
                            EllipseAnnulusPixelRegion, EllipseAnnulusSkyRegion,
                            EllipsePixelRegion, EllipseSkyRegion,
                            LinePixelRegion, LineSkyRegion, PointPixelRegion,
                            PointSkyRegion, PolygonPixelRegion,
                            PolygonSkyRegion, RectangleAnnulusPixelRegion,
                            RectangleAnnulusSkyRegion, RectanglePixelRegion,
                            RectangleSkyRegion, RegularPolygonPixelRegion,
                            TextPixelRegion, TextSkyRegion)

__all__ = []

regions_attributes = dict(circle=['center', 'radius'],
                          ellipse=['center', 'width', 'height', 'angle'],
                          rectangle=['center', 'width', 'height', 'angle'],
                          polygon=['vertices'],
                          circleannulus=['center', 'inner_radius',
                                         'outer_radius'],
                          ellipseannulus=['center', 'inner_width',
                                          'inner_height', 'outer_width',
                                          'outer_height', 'angle'],
                          line=['start', 'end'],
                          point=['center'],
                          text=['center'])

regions_attributes['rectangleannulus'] = regions_attributes['ellipseannulus']

# Map the region names in the respective format to the ones available in
# this package
reg_mapping = {'CRTF': {x: x for x in regions_attributes}}
reg_mapping['CRTF']['rotbox'] = 'rectangle'
reg_mapping['CRTF']['box'] = 'rectangle'
reg_mapping['CRTF']['centerbox'] = 'rectangle'
reg_mapping['CRTF']['poly'] = 'polygon'
reg_mapping['CRTF']['symbol'] = 'point'
reg_mapping['CRTF']['text'] = 'text'
reg_mapping['CRTF']['annulus'] = 'circleannulus'

# valid astropy coordinate frames in their respective formats
valid_coordsys = {'CRTF': ['image', 'fk5', 'fk4', 'galactic',
                           'geocentrictrueecliptic', 'supergalactic', 'icrs']}

# Map astropy's coordinate frame names with their respective name in
# the file format.
coordsys_mapping = {'CRTF': {x: x.upper() for x in valid_coordsys['CRTF']}}

# Coordinate frame names must be uppercase following the CASA CRTF syntax
coordsys_mapping['CRTF']['geocentrictrueecliptic'] = 'ECLIPTIC'
coordsys_mapping['CRTF']['fk5'] = 'J2000'
coordsys_mapping['CRTF']['fk4'] = 'B1950'
coordsys_mapping['CRTF']['supergalactic'] = 'SUPERGAL'


class RegionConversionError(ValueError):
    """
    A generic error class for Shape to Region conversions.
    """


class _ShapeList(list):
    """
    A class to hold a list of `~regions.Shape` objects.
    """

    def to_regions(self):
        """
        Convert to a list of `~regions.Region` objects.
        """
        regions = []
        for shape in self:
            # Skip elliptical multi-annulus for now
            if shape.region_type == 'ellipse' and len(shape.coord) > 5:
                msg = f'Skipping elliptical annulus {shape}'
                warn(msg, AstropyUserWarning)
                continue

            if shape.region_type in ['box'] and shape.format_type == 'CRTF':
                msg = f'Skipping {shape.region_type} {shape}'
                warn(msg, CRTFRegionParserWarning)
                continue

            region = shape.to_region()
            regions.append(region)

        return regions

    def to_crtf(self, coordsys='fk5', fmt='.6f', radunit='deg'):
        """
        Convert to CRTF (CASA Region Text Format) region strings.

        Parameters
        ----------
        coordsys : str
            An Astropy coordinate system that overrides the coordinate
            system frame for all regions.

        fmt : str
            A python string format defining the output precision.
            Default is '.6f', which is accurate to 0.0036 arcseconds.

        radunit : str
            The unit of the radius.

        Returns
        -------
        region_string : str
            A CRFT region string.
        """
        crtf_strings = {
            'circle': '{0}circle[[{1:FMT}deg, {2:FMT}deg], {3:FMT}RAD]',
            'circleannulus': ('{0}annulus[[{1:FMT}deg, {2:FMT}deg], '
                              '[{3:FMT}RAD, {4:FMT}RAD]]'),
            # Make sure that width goes to minor axis and height to
            # major axis
            'ellipse': ('{0}ellipse[[{1:FMT}deg, {2:FMT}deg], [{4:FMT}RAD, '
                        '{3:FMT}RAD], {5:FMT}deg]'),
            'rectangle': ('{0}rotbox[[{1:FMT}deg, {2:FMT}deg], [{3:FMT}RAD, '
                          '{4:FMT}RAD], {5:FMT}deg]'),
            'polygon': '{0}poly[{1}]',
            'point': '{0}point[[{1:FMT}deg, {2:FMT}deg]]',
            'symbol': '{0}symbol[[{1:FMT}deg, {2:FMT}deg], {symbol}]',
            'text': '{0}text[[{1:FMT}deg, {2:FMT}deg], \'{text}\']',
            'line': ('{0}line[[{1:FMT}deg, {2:FMT}deg], [{3:FMT}deg, '
                     '{4:FMT}deg]]')}

        output = '#CRTFv0\n'

        if radunit == 'arcsec':
            # arcseconds are allowed for all but image coordinates
            if coordsys.lower() not in ('image',):
                radunitstr = '"'
            else:
                raise ValueError('Radius unit arcsec not valid for '
                                 f'coordsys {coordsys}')
        else:
            radunitstr = radunit

        for key, val in crtf_strings.items():
            crtf_strings[key] = val.replace('FMT', fmt).replace('RAD',
                                                                radunitstr)

        # CASA does not support global coordinate specification, even
        # though the documentation for the specification explicitly
        # states that it does.
        global_coord = coordsys_mapping['CRTF'][coordsys.lower()]
        output += f'global coord={global_coord}\n'

        for shape in self:
            shape.check_crtf()
            shape.meta = _to_crtf_meta(shape.meta)

            # if unspecified, include is True. Despite the
            # specification, CASA does *not* support a preceding "+". If
            # you want a region included, leave the opening character
            # blank.
            include = ''
            if shape.include in (False, '-'):
                include = '-'

            if shape.meta.get('type', 'reg') == 'ann':
                include += 'ann '

            if shape.meta.get('label', '') != '':
                shape.meta['label'] = f"'{shape.meta['label']}'"

            keylist = ('include', 'comment', 'symbol', 'coord', 'text',
                       'range', 'corr', 'type')
            meta_pairs = []
            for key, val in shape.meta.items():
                if key not in keylist:
                    meta_pairs.append(f'{key}={val}')
            meta_str = ', '.join(meta_pairs)

            # The first item should be the coordinates, since CASA
            # cannot recognize a region without an inline coordinate
            # specification. It can be, but does not need to be,
            # comma-separated at the start.
            shape_coordsys = shape.coordsys
            if shape_coordsys.lower() != coordsys.lower():
                coord = coordsys_mapping['CRTF'][coordsys.lower()]
                if meta_str.strip():
                    meta_str = f'coord={coord}, ' + meta_str
                else:
                    # if there is no metadata at all (above), the
                    # trailing comma is incorrect
                    meta_str = f'coord={coord}'

            if 'comment' in shape.meta:
                meta_str += ', ' + shape.meta['comment']

            if 'range' in shape.meta:
                shape.meta['range'] = [str(str(x).replace(' ', '')) for x in
                                       shape.meta['range']]
                meta_str += f", range={shape.meta['range']}".replace("'", '')
            if 'corr' in shape.meta:
                meta_str += f", corr={shape.meta['corr']}".replace("'", '')

            coord = []
            if coordsys not in ['image', 'physical']:
                for val in shape.coord:
                    if (isinstance(val, Angle)
                            or (radunit == '' or radunit is None)):
                        coord.append(float(val.value))
                    else:
                        coord.append(float(val.to(radunit).value))
            else:
                for val in shape.coord:
                    if isinstance(val, u.Quantity):
                        coord.append(float(val.value))
                    else:
                        coord.append(float(val))

            if (shape.region_type in ['ellipse', 'rectangle']
                    and len(shape.coord) % 2 == 1):
                coord[-1] = float(shape.coord[-1].to('deg').value)

            if shape.region_type == 'polygon':
                vals = [f'[{x:{fmt}}deg, {y:{fmt}}deg]'
                        for x, y in zip(coord[::2], coord[1::2], strict=True)]
                coord = ', '.join(vals)
                line = crtf_strings['polygon'].format(include, coord)

            elif shape.region_type == 'point':
                if 'symbol' in shape.meta:
                    line = crtf_strings['symbol'].format(
                        include, *coord, symbol=shape.meta['symbol'])
                else:
                    line = crtf_strings['point'].format(include, *coord)

            elif shape.region_type == 'ellipse':
                coord[2:] = [x / 2 for x in coord[2:]]
                if len(coord) % 2 == 1:
                    coord[-1] *= 2
                line = crtf_strings['ellipse'].format(include, *coord)

            elif shape.region_type == 'text':
                line = crtf_strings['text'].format(
                    include, *coord, text=shape.meta['text'])
            else:
                line = crtf_strings[shape.region_type].format(include, *coord)

            if meta_str.strip():
                output += f'{line}, {meta_str}\n'
            else:
                output += f'{line}\n'

        return output


class _Shape:
    """
    Helper class to represent a CRTF Region.

    This serves as intermediate step in the parsing process.

    Parameters
    ----------
    coordsys : str
        An Astropy Coordinate system frame used in the region.

    region_type : str
        The type of the region (as defined in this package).

    coord : list of `~astropy.coordinates.Angle` or `~astropy.units.Quantity`
        The region coordinates.

    meta : dict
        The meta attributes.

    composite : bool
        Flag indicting whether the region is a Composite region.

    include : bool
        Flag indicating where to include or exclude the region.
    """

    shape_to_sky_region = dict(circle=CircleSkyRegion,
                               ellipse=EllipseSkyRegion,
                               rectangle=RectangleSkyRegion,
                               polygon=PolygonSkyRegion,
                               circleannulus=CircleAnnulusSkyRegion,
                               ellipseannulus=EllipseAnnulusSkyRegion,
                               rectangleannulus=RectangleAnnulusSkyRegion,
                               line=LineSkyRegion,
                               point=PointSkyRegion,
                               text=TextSkyRegion)

    shape_to_pixel_region = dict(circle=CirclePixelRegion,
                                 ellipse=EllipsePixelRegion,
                                 rectangle=RectanglePixelRegion,
                                 polygon=PolygonPixelRegion,
                                 circleannulus=CircleAnnulusPixelRegion,
                                 ellipseannulus=EllipseAnnulusPixelRegion,
                                 rectangleannulus=RectangleAnnulusPixelRegion,
                                 line=LinePixelRegion,
                                 point=PointPixelRegion,
                                 text=TextPixelRegion)

    error = RegionConversionError

    def __init__(self, coordsys, region_type, coord, meta, composite, include):
        self._coordsys = coordsys
        self._region_type = region_type
        self.coord = coord
        self.meta = meta
        self.composite = composite
        self.include = include

    @property
    def coordsys(self):
        return self._coordsys

    @coordsys.setter
    def coordsys(self, value):
        self._coordsys = value.lower()
        self._validate()

    @property
    def region_type(self):
        return self._region_type

    @region_type.setter
    def region_type(self, value):
        self._region_type = value.lower()
        self._validate()

    def __str__(self):
        ss = self.__class__.__name__
        ss += f'\nType: {self.meta.get("type", "reg")}'
        ss += f'\nCoord sys: {self.coordsys}'
        ss += f'\nRegion type: {self.region_type}'
        if self.region_type == 'symbol':
            ss += f"\nSymbol: {self.meta['symbol']}"
        if self.region_type == 'text':
            ss += f"\nText: {self.meta['text']}"
        ss += f'\nMeta: {self.meta}'
        ss += f'\nComposite: {self.composite}'
        ss += f'\nInclude: {self.include}'
        ss += '\n'
        return ss

    def convert_coords(self):
        """
        Process a list of coordinates.

        This mainly searches for a tuple of coordinates in the
        coordinate list and creates a SkyCoord or PixCoord object from
        them if appropriate for a given region type. This involves again
        some coordinate transformation, so this step could be moved to
        the parsing process.
        """
        if self.coordsys in ['image', 'physical']:
            coords = self._convert_pix_coords()
        else:
            coords = self._convert_sky_coords()

        if self.region_type == 'line':
            coords = [coords[0][0], coords[0][1]]

        if self.region_type == 'text':
            coords.append(self.meta['text'])

        return coords

    def _convert_sky_coords(self):
        """
        Convert to sky coordinates.
        """
        parsed_angles = []
        for x, y in zip(self.coord[:-1:2], self.coord[1::2], strict=True):
            if isinstance(x, Angle) and isinstance(y, Angle):
                parsed_angles.append((x, y))

        frame = frame_transform_graph.lookup_name(self.coordsys)

        if len(parsed_angles) == 0:
            raise ValueError('error parsing region')

        lon, lat = zip(*parsed_angles, strict=True)
        if (hasattr(lon, '__len__') and hasattr(lat, '__len__')
                and len(lon) == 1 and len(lat) == 1):
            # force entries to be scalar if they are length-1
            lon, lat = u.Quantity(lon[0]), u.Quantity(lat[0])
        else:
            # otherwise, they are vector quantities
            lon, lat = u.Quantity(lon), u.Quantity(lat)
        sphcoords = UnitSphericalRepresentation(lon, lat)
        coords = [SkyCoord(frame(sphcoords))]

        if self.region_type != 'polygon':
            coords += self.coord[len(coords * 2):]

        return coords

    def _convert_pix_coords(self):
        """
        Convert to pixel coordinates.
        """
        if self.region_type in ['polygon', 'line']:
            # have to special-case polygon in the phys coord case
            # b/c can't typecheck when iterating as in sky coord case
            coords = [PixCoord(self.coord[0::2], self.coord[1::2])]
        else:
            temp = [_.value for _ in self.coord]
            coord = PixCoord(temp[0], temp[1])
            coords = [coord] + temp[2:]

        # The angle remains as a quantity object. Modulus check makes
        # sure that it works for ellipse/rectangle annulus.
        if (self.region_type in ['ellipse', 'rectangle']
                and len(coords) % 2 == 0):
            coords[-1] = self.coord[-1]

        return coords

    def to_region(self):
        """
        Convert to a `~regions.Region` object.
        """
        coords = self.convert_coords()
        viz_keywords = ['color', 'dash', 'dashlist', 'width', 'font',
                        'symsize', 'symbol', 'symsize', 'fontsize',
                        'fontstyle', 'usetex', 'labelpos', 'labeloff',
                        'labelcolor', 'linewidth', 'linestyle', 'point',
                        'textangle', 'fontweight', 'symthick', 'default_style',
                        'fill', 'textrotate', 'fontstyle']

        if isinstance(coords[0], SkyCoord):
            reg = self.shape_to_sky_region[self.region_type](*coords)
        elif isinstance(coords[0], PixCoord):
            reg = self.shape_to_pixel_region[self.region_type](*coords)
        else:
            self._raise_error('No central coordinate')

        reg.visual = RegionVisual()
        reg.meta = RegionMeta()

        # both 'text' and 'label' should be set to the same value, where
        # we default to the 'text' value since that is the one used by
        # ds9 regions
        label = self.meta.get('text', self.meta.get('label', ''))
        if label != '':
            reg.meta['label'] = label
        for key in self.meta:
            if key in viz_keywords:
                reg.visual[key] = self.meta[key]
            else:
                reg.meta[key] = self.meta[key]
        reg.meta['include'] = self.include

        return reg

    def _raise_error(self, msg):
        raise self.error(msg)

    def check_crtf(self):
        """
        Check for CRTF compatibility.
        """
        if self.region_type not in regions_attributes:
            raise ValueError(f'"{self.region_type}" is not a valid region '
                             'type')

        if self.coordsys not in valid_coordsys['CRTF']:
            raise ValueError(f'"{self.coordsys}" is not a valid coordinate '
                             'reference frame')

    def _validate(self):
        """
        Check whether all the attributes of this object is valid.
        """
        if self.region_type not in regions_attributes:
            raise ValueError(f'"{self.region_type}" is not a valid region '
                             'type')

        if self.coordsys not in valid_coordsys['CRTF']:
            raise ValueError(f'"{self.coordsys}" is not a valid coordinate '
                             'reference frame')


def _to_shape_list(region_list, coordinate_system='fk5'):
    """
    Convert a list of regions into a `~regions.ShapeList` object.

    Parameters
    ----------
    region_list : list
        A list of `regions.Region` objects.

    coordinate_system : str, optional
        The Astropy coordinate system frame in which all the coordinates
        present in the ``region_list`` will be converted. Default is
        'fk5'.

    Returns
    -------
    shape_list: `regions.ShapeList`
        A `~regions.ShapeList` object.
    """
    shape_list = _ShapeList()

    for region in region_list:
        if isinstance(region, SkyRegion):
            reg_type = region.__class__.__name__[:-9].lower()
        elif isinstance(region, RegularPolygonPixelRegion):
            reg_type = 'polygon'
        else:
            reg_type = region.__class__.__name__[:-11].lower()

        if reg_type == 'polygon':
            coord = region.vertices
        else:
            coord = []
            for val in regions_attributes[reg_type]:
                coord.append(getattr(region, val))

        if coordinate_system:
            coordsys = coordinate_system
        else:
            coordsys = (coord[0].name
                        if isinstance(region, SkyRegion) else 'image')

        new_coord = []
        for val in coord:
            if isinstance(val, Angle):
                # convert Angle to Quantity; Angle values get units
                # stripped in serialization, but Quantity gets converted
                new_coord.append(u.Quantity(val))
            elif isinstance(val, (u.Quantity, numbers.Number)):
                new_coord.append(val)
            elif isinstance(val, PixCoord):
                new_coord.append(u.Quantity(val.x, u.dimensionless_unscaled))
                new_coord.append(u.Quantity(val.y, u.dimensionless_unscaled))
            else:
                frame = frame_transform_graph.lookup_name(coordsys)
                new_coord.append(Angle(val.transform_to(frame).spherical.lon))
                new_coord.append(Angle(val.transform_to(frame).spherical.lat))

        meta = dict(region.meta)
        meta.update(region.visual)

        if reg_type == 'text':
            meta['text'] = meta.get('text', meta.pop('label', ''))

        include = region.meta.pop('include', True)

        shape_list.append(_Shape(coordsys, reg_type, new_coord, meta, False,
                                 include))

    return shape_list


def _to_crtf_meta(shape_meta):
    """
    Make the metadata CRTF compatible by filtering and mapping the valid
    keys.

    Parameters
    ----------
    shape_meta : dict
        A meta attribute of a `regions.Region` object.

    Returns
    -------
    meta : dict
        A CRTF compatible meta dictionary.
    """
    # please refer : https://casaguides.nrao.edu/index.php/CASA_Region_Format

    # meta keys allowed in CRTF
    valid_keys = ['label', 'include', 'frame', 'range', 'veltype',
                  'restfreq', 'coord', 'type', 'text', 'corr']

    # visual keys allowed in CRTF
    valid_keys += ['color', 'width', 'font', 'symthick', 'symsize',
                   'fontsize', 'fontstyle', 'usetex', 'labelpos', 'labeloff',
                   'linewidth', 'linestyle', 'symbol']

    key_mappings = {}

    return _to_io_meta(shape_meta, valid_keys, key_mappings)


def _to_io_meta(shape_meta, valid_keys, key_mappings):
    """
    Make metadata compatible with a specific IO by filtering and mapping
    to its valid keys.

    Parameters
    ----------
    shape_meta : dict
        A meta attribute of a `regions.Region` object.

    valid_keys : list
        The valid keys of a particular file format.

    key_mappings : dict
        A dictionary mapping of the actual name of the key in the
        format.

    Returns
    -------
    meta : dict
        An IO compatible meta dictionary according to ``valid_keys`` and
        ``key_mappings``.
    """
    meta = {}

    for key in shape_meta:
        if key in valid_keys:
            meta[key_mappings.get(key, key)] = shape_meta[key]

    return meta
