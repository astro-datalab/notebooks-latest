# Licensed under a 3-clause BSD style license - see LICENSE.rst

import copy
import itertools
import re
from warnings import warn

import astropy.units as u
from astropy.coordinates import Angle, frame_transform_graph
from astropy.utils.data import get_readable_fileobj

from regions.core import Regions
from regions.core.registry import RegionsRegistry
from regions.io.crtf.core import (CRTFRegionParserError,
                                  CRTFRegionParserWarning, valid_symbols)
from regions.io.crtf.io_core import _Shape, _ShapeList, reg_mapping

__all__ = []

# All CASA files start with '#CRTF' . It may also include the version
# number like '#CRTFv0' .
regex_begin = re.compile(r'^#CRTFv?[\d]?')

# Comment format
regex_comment = re.compile(r'^#.*$')

# Identifies the global attributes Format
regex_global = re.compile(r'^global\s+(?P<parameters>.*)?')

# Coordinate Format : '[x, y]'
regex_coordinate = re.compile(r'\[([\w.+-:]*?)\s*[,]\s*([\w.+-:]*?)\]')

# Single length format, e.g., helps extract the radius of a circle
regex_length = re.compile(r'(?:\[[^=\]]*\])+[,]\s*([^\[]*)\]')

# Extracts each 'parameter=value' pair
regex_meta = re.compile(r'(?:(\w+)\s*=[\s\'\"]*([^,\[\]]+?)[\'\",]+)|(?:(\w+)\s*=\s*\[(.*?)\])')  # noqa: E501

# Region format which segregates the include ('+'|'-') parameter, the
# kind of definition ('ann' for annotations or '' for regions) and region
# type.
regex_region = re.compile(r'(?P<include>[+-])?(?P<type>ann(?=\s))?\s*(?P<regiontype>[a-z]*?)\s?\[[^=]*]')  # noqa: E501

# Line format which checks the validity of the line and segregates the
# meta attributes from the region format.
regex_line = re.compile(r'(?P<region>[+-]?(?:ann(?=\s))?\s*[a-z]+?\s?\[[^=]+\])(?:\s*,?\s*(?P<parameters>.*))?')  # noqa: E501


@RegionsRegistry.register(Regions, 'read', 'crtf')
def _read_crtf(filename, errors='strict', cache=False):
    """
    Read a CRTF region file and return a list of region objects.

    Parameters
    ----------
    filename : str
        The filename of the file to access.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.CRTFRegionParserError`. 'warn' will raise a
        `~regions.CRTFRegionParserWarning`, and 'ignore' will do nothing
        (i.e., be silent).

    cache : bool or 'update', optional
        Whether to cache the contents of remote URLs. If 'update', check
        the remote URL for a new version but store the result in the
        cache.

    Returns
    -------
    regions : list
        A list of `~regions.Region` objects.
    """
    with get_readable_fileobj(filename) as fh:
        if regex_begin.search(fh.readline()):
            region_string = fh.read()
            return _parse_crtf(region_string, errors=errors)
        else:
            raise CRTFRegionParserError('Every CRTF Region must start with '
                                        '"#CRTF"')


@RegionsRegistry.register(Regions, 'parse', 'crtf')
def _parse_crtf(region_string, errors='strict'):
    parser = _CRTFParser(region_string, errors=errors)
    return Regions(parser.shapes.to_regions())


class _CRTFParser:
    """
    Parse a CRTF string.

    This class transforms a CRTF string to a
    `~regions.io.core.ShapeList`. The result is stored in the ``shapes``
    attribute.

    Parameters
    ----------
    region_string : str
        A CRTF region string.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.CRTFRegionParserError`. 'warn' will raise a
        `~regions.CRTFRegionParserWarning`, and 'ignore' will do nothing
        (i.e., be silent).
    """

    # Each line is tested for either containing a region with meta
    # attributes or global parameters. If meta attributes or global
    # parameters are found, they are stored in the ``global_meta``
    # attribute. If a region is found, then ``_CRTFRegionParser``
    # is invoked to transform the line into a `~regions.io.core.Shape`
    # object, which is stored in the ``shapes`` attribute.

    # valid definition (region or annotation) types
    valid_definition = ('box', 'centerbox', 'rotbox', 'poly', 'circle',
                        'annulus', 'ellipse', 'line', 'vector', 'text',
                        'symbol')

    # valid parameters (attributes)
    valid_global_keys = ('coord', 'frame', 'corr', 'veltype', 'restfreq',
                         'linewidth', 'linestyle', 'symsize', 'symthick',
                         'color', 'font', 'fontsize', 'fontstyle', 'usetex',
                         'labelpos', 'labelcolor', 'labeloff', 'range')

    def __init__(self, region_string, errors='strict'):
        if errors not in ('strict', 'ignore', 'warn'):
            raise ValueError('errors must be one of "strict", "ignore", or '
                             '"warn"')

        self.region_string = region_string
        self.errors = errors
        self.global_meta = {}  # global states
        self.shapes = _ShapeList()
        self.run()

    def __str__(self):
        ss = self.__class__.__name__
        ss += f'\nErrors: {self.errors}'
        ss += f'\nGlobal meta: {self.global_meta}'
        ss += f'\nShapes: {self.shapes}'
        ss += '\n'
        return ss

    def parse_line(self, line):
        """
        Parse a single line.
        """
        # Skip blanks
        if line == '':
            return

        # Skip comments
        if regex_comment.search(line):
            return

        # Special case / header: parse global parameters into metadata
        global_parameters = regex_global.search(line)
        if global_parameters:
            self.parse_global_meta(global_parameters.group('parameters'))
            return

        # Tries to check the validity of the line.
        crtf_line = regex_line.search(line)
        if crtf_line:
            # Tries to parse the line.

            # Finds info about the region.
            region = regex_region.search(crtf_line.group('region'))
            type_ = region.group('type') or 'reg'
            include = region.group('include') or '+'
            region_type = region.group('regiontype').lower()
            if region_type in self.valid_definition:
                helper = _CRTFRegionParser(
                    self.global_meta, include, type_, region_type,
                    *crtf_line.group('region', 'parameters'))

                self.shapes.append(helper.shape)
            else:
                self._raise_error('Not a valid CRTF Region type: '
                                  f'"{region_type}".')

        else:
            self._raise_error(f'Not a valid CRTF line: "{line}".')
            return

    def _raise_error(self, msg):
        if self.errors == 'warn':
            warn(msg, CRTFRegionParserWarning)
        elif self.errors == 'strict':
            raise CRTFRegionParserError(msg)

    def run(self):
        """
        Run all the steps.

        This function splits the regions into lines and calls
        ``parse_line`` for each line.
        """
        for line in self.region_string.split('\n'):
            self.parse_line(line)

    def parse_global_meta(self, global_meta_str):
        """
        Parse the line starting with global to extract all the valid
        meta key/value pair.
        """
        if global_meta_str:
            global_meta_str = regex_meta.findall(global_meta_str + ',')
        if global_meta_str:
            for par in global_meta_str:
                if par[0] != '':
                    val1 = par[0].lower()
                    val2 = par[1]
                else:
                    val1 = par[2].lower()
                    val2 = par[3]
                val1 = val1.strip()
                val2 = val2.strip()
                if val1 in self.valid_global_keys:
                    if val1 in ('range', 'corr', 'labeloff'):
                        val2 = val2.split(',')
                        val2 = [x.strip() for x in val2 if x]
                    self.global_meta[val1] = val2
                else:
                    self._raise_error(f'"{val1}" is not a valid global meta '
                                      'key')


class _CRTFRegionParser:
    """
    Parse a CRTF region string.

    This will turn a line containing a CRTF region into a
    `~regions.Shape` object.

    Parameters
    ----------
    global_meta : dict
        Global meta data of the CRTF (CASA Region Text Format) file
        which is used as default meta values for regions.

    include : {'+', '-'}
        Flag at the beginning of the line.

    type_ : {'reg', 'ann'}
        Kind of the region definition.

    region_type : str
        Region type.

    reg_str : str
        Region string to parse.

    meta_str : str
        Meta string to parse.

    errors : {'strict', 'warn', 'ignore'}, optional
        The error handling scheme to use for handling parsing
        errors. The default is 'strict', which will raise a
        `~regions.CRTFRegionParserError`. 'warn' will raise a
        `~regions.CRTFRegionParserWarning`, and 'ignore' will do nothing
        (i.e., be silent).
    """

    # List of valid coordinate system
    # TODO : There are still many reference systems to support
    coordinate_systems = ['j2000', 'icrs', 'galactic', 'supergal', 'image',
                          'ecliptic']

    # Maps CASA coordinate frame to appropriate astropy coordinate frames.
    coordsys_mapping = dict(zip(frame_transform_graph.get_names(),
                                frame_transform_graph.get_names(),
                                strict=True))
    coordsys_mapping['j2000'] = 'fk5'
    coordsys_mapping['b1950'] = 'fk4'
    coordsys_mapping['supergal'] = 'supergalactic'
    coordsys_mapping['ecliptic'] = 'geocentrictrueecliptic'

    # CRTF Format specifications defining how a region is read:
    #   'c' denotes a coordinates
    #   'l' denotes a length
    #   'pl' denotes a pair of lengths
    #   's' denotes a string (generally text or a symbol)
    language_spec = {'circle': ['c', 'l'],
                     'box': ['c', 'c'],
                     'centerbox': ['c', 'pl'],
                     'rotbox': ['c', 'pl', 'l'],
                     'poly': itertools.cycle('c'),
                     'annulus': ['c', 'pl'],
                     'ellipse': ['c', 'pl', 'l'],
                     'line': ['c', 'c'],
                     'vector': ['c', 'c'],
                     'symbol': ['c', 's'],
                     'text': ['c', 's']}

    def __init__(self, global_meta, include, type_, region_type, reg_str,
                 meta_str, errors='strict'):
        self.global_meta = global_meta
        self.reg_str = reg_str
        self.meta_str = meta_str
        self.errors = errors

        self.coord = None
        self.coordsys = None
        self.coord_str = None
        self.type_ = type_
        self.region_type = region_type
        self.meta = copy.deepcopy(global_meta)
        self.shape = None
        self.include = include or '+'

        self.parse()

    def _raise_error(self, msg):
        if self.errors == 'warn':
            warn(msg, CRTFRegionParserWarning)
        elif self.errors == 'strict':
            raise CRTFRegionParserError(msg)

    def parse(self):
        """
        Parse the CRTF region string.
        """
        self.convert_meta()
        self.coordsys = self.meta.get('coord', 'image').lower()
        self.set_coordsys()
        self.convert_coordinates()
        self.make_shape()

    def set_coordsys(self):
        """
        Map to astropy's coordinate system frame name.
        """
        # TODO: needs expert attention (most reference systems are not
        # mapped)
        if self.coordsys.lower() in self.coordsys_mapping:
            self.coordsys = self.coordsys_mapping[self.coordsys.lower()]

    def convert_coordinates(self):
        """
        Convert coordinate string to `~astropy.coordinates.Angle` or
        `~astropy.units.quantity.Quantity` objects.
        """
        coord_list_str = (regex_coordinate.findall(self.reg_str)
                          + regex_length.findall(self.reg_str))
        coord_list = []

        if self.region_type == 'poly':
            if len(coord_list_str) < 3:
                self._raise_error(f'Not in proper format: {self.reg_str} '
                                  'polygon should have >= 3 coordinates')
        else:
            expected_nparam = len(self.language_spec[self.region_type])
            if len(coord_list_str) != expected_nparam:
                self._raise_error(f'Not in proper format: "{self.reg_str}". '
                                  'Does not contain expected number of '
                                  'parameters for the region '
                                  f'"{self.region_type}"')

        # TODO: check zip strict=True
        for attr_spec, val_str in zip(self.language_spec[self.region_type],
                                      coord_list_str):
            if attr_spec == 'c':
                if len(val_str) == 2 and val_str[1] != '':
                    coord_list.append(
                        _CRTFCoordinateParser.parse_coordinate(val_str[0]))
                    coord_list.append(
                        _CRTFCoordinateParser.parse_coordinate(val_str[1]))
                else:
                    self._raise_error(f'Not in proper format: {val_str} '
                                      'should be a coordinate')

            if attr_spec == 'pl':
                if len(val_str) == 2 and val_str[1] != '':
                    coord_list.append(
                        _CRTFCoordinateParser.parse_angular_length_quantity(
                            val_str[0]))
                    coord_list.append(
                        _CRTFCoordinateParser.parse_angular_length_quantity(
                            val_str[1]))
                else:
                    self._raise_error(f'Not in proper format: {val_str} '
                                      'should be a pair of lengths')

            if attr_spec == 'l':
                if isinstance(val_str, str):
                    coord_list.append(
                        _CRTFCoordinateParser.parse_angular_length_quantity(
                            val_str))
                else:
                    self._raise_error(f'Not in proper format: {val_str} '
                                      'should be a single length')

            if attr_spec == 's':
                if self.region_type == 'symbol':
                    if val_str in valid_symbols:
                        self.meta['symbol'] = val_str
                    else:
                        self._raise_error(f'Not in proper format: "{val_str}" '
                                          'should be a symbol')
                elif self.region_type == 'text':
                    self.meta['text'] = val_str[1:-1]

        self.coord = coord_list

    def convert_meta(self):
        """
        Parse the meta_str to a dictionary and store in the ``meta``
        attribute.
        """
        if self.meta_str:
            self.meta_str = regex_meta.findall(self.meta_str + ',')
        if self.meta_str:
            for par in self.meta_str:
                if par[0] != '':
                    val1 = par[0]
                    val2 = par[1]
                else:
                    val1 = par[2]
                    val2 = par[3]
                val1 = val1.strip()
                val2 = val2.strip()
                if val1 in _CRTFParser.valid_global_keys or val1 == 'label':
                    if val1 in ('range', 'corr', 'labeloff'):
                        val2 = val2.split(',')
                        val2 = [x.strip() for x in val2]
                    self.meta[val1] = val2
                else:
                    self._raise_error(f'"{val1}" is not a valid meta key')

        self.meta['include'] = self.include != '-'
        self.include = self.meta['include']

        if 'range' in self.meta:
            self.meta['range'] = [u.Quantity(x) for x in self.meta['range']]

        self.meta['type'] = self.type_

    def make_shape(self):
        """
        Make a `~regions.Shape` object.
        """
        if self.region_type == 'ellipse':
            self.coord[2:] = [x * 2 for x in self.coord[2:]]
            # Map major and minor axis to height and width respectively
            self.coord[2], self.coord[3] = self.coord[3], self.coord[2]
            if len(self.coord) % 2 == 1:  # check if angle is present
                self.coord[-1] /= 2

        if self.region_type == 'box':
            x = (self.coord[0] + self.coord[2]) / 2
            y = (self.coord[1] + self.coord[3]) / 2
            w = u.Quantity(self.coord[0] - self.coord[2])
            h = u.Quantity(self.coord[1] - self.coord[3])
            self.coord = [x, y, abs(w), abs(h)]

        self.meta.pop('coord', None)

        self.shape = _Shape(coordsys=self.coordsys,
                            region_type=reg_mapping['CRTF'][self.region_type],
                            coord=self.coord,
                            meta=self.meta,
                            composite=False,
                            include=self.include)


class _CRTFCoordinateParser:
    """
    Helper class to structure coordinate parser.
    """

    @staticmethod
    def parse_coordinate(string_rep):
        """
        Parse a single coordinate.
        """
        # Any CRTF coordinate representation (sexagesimal or degrees)
        if 'pix' in string_rep:
            return u.Quantity(string_rep[:-3], u.dimensionless_unscaled)
        if 'h' in string_rep or 'rad' in string_rep:
            return Angle(string_rep)

        # CRTF format apparently implicitly defines
        # 18:20:30.12
        # strings as hour units, and
        # 10.11.54.69
        # as degrees.
        # I can't find this definition in documentation anywhere, though.
        unit = u.deg
        if len(string_rep.split('.')) >= 3:
            string_rep = string_rep.replace('.', ':', 2)
        elif string_rep.count(':') == 2:
            unit = u.hour

        return Angle(string_rep, unit)

    @staticmethod
    def parse_angular_length_quantity(string_rep):
        """
        Parse a string into a Quantity object.

        Given a string that is a number and a unit, return a Quantity of
        that string. An error is raised if there is no unit, e.g.:

        * 50" -> 50 * u.arcsec
        * 50 -> CRTFRegionParserError : Units must be specified for 50
        """
        unit_mapping = {'deg': u.deg,
                        'rad': u.rad,
                        'arcmin': u.arcmin,
                        'arcsec': u.arcsec,
                        'pix': u.dimensionless_unscaled,
                        '"': u.arcsec,
                        "'": u.arcmin}

        regex_str = re.compile(r'([0-9+,-.]*)(.*)')
        restr = regex_str.search(string_rep)
        unit = restr.group(2)
        if unit:
            if unit in unit_mapping:
                return u.Quantity(restr.group(1), unit=unit_mapping[unit])
            return u.Quantity(restr.group(1))
        else:
            raise CRTFRegionParserError('Units must be specified for '
                                        f'{string_rep}')
