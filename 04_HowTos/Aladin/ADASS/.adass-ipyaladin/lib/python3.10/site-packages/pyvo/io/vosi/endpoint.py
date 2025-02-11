# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains a contains the high-level functions to read the various
VOSI Endpoints.
"""

from astropy.utils.xml import iterparser
from astropy.utils.collections import HomogeneousList
from astropy.io.votable.exceptions import vo_raise, vo_warn
from astropy.io.votable.util import version_compare

from ...utils.xml.elements import xmlattribute, xmlelement, Element
from . import voresource as vr
from . import vodataservice as vs
from . import availability as av
from .exceptions import W15, W16, E07, E10

__all__ = [
    "parse_tables", "parse_capabilities", "parse_availability",
    "TablesFile", "CapabilitiesFile", "AvailabilityFile"]


def _pedantic_settings(pedantic):
    """
    Controls the pedantic parser settings.  Based on the bool
    passed in to pedantic, create a config to be passed to
    astropy parsing to raise exceptions or ignore them on
    pedantic errors.

    Parameters
    ----------
    pedantic : bool
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        the standard Python mechanisms.  See the `warnings`
        module in the Python standard library for more information.

    Returns
    -------
    A dict containing 'verify' configuration settings.
    """
    if pedantic:
        return {'verify': 'exception'}
    else:
        return {'verify': 'warn'}


def parse_tables(source, *, pedantic=None, filename=None,
                 _debug_python_based_parser=False):
    """
    Parses a tableset xml file (or file-like object), and returns a
    `~pyvo.io.vosi.endpoint.TablesFile` object.

    Parameters
    ----------
    source : str or readable file-like object
        Path or file object containing a tableset xml file.
    pedantic : bool, optional
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        the standard Python mechanisms.  See the `warnings`
        module in the Python standard library for more information.
        Defaults to False.
    filename : str, optional
        A filename, URL or other identifier to use in error messages.
        If *filename* is None and *source* is a string (i.e. a path),
        then *source* will be used as a filename for error messages.
        Therefore, *filename* is only required when source is a
        file-like object.

    Returns
    -------
    tables_file : `~pyvo.io.vosi.endpoint.TablesFile` object

    See also
    --------
    pyvo.io.vosi.exceptions : The exceptions this function may raise.
    """
    config = _pedantic_settings(pedantic)

    if filename is None and isinstance(source, str):
        config['filename'] = source
    else:
        config['filename'] = filename

    with iterparser.get_xml_iterator(
        source,
        _debug_python_based_parser=_debug_python_based_parser
    ) as iterator:
        return TablesFile(
            config=config, pos=(1, 1)).parse(iterator, config)


def parse_capabilities(source, *, pedantic=None, filename=None,
                       _debug_python_based_parser=False):
    """
    Parses a capabilities xml file (or file-like object), and returns a
    `~pyvo.io.vosi.endpoint.CapabilitiesFile` object.

    Parameters
    ----------
    source : str or readable file-like object
        Path or file object containing a capabilities xml file.
    pedantic : bool, optional
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        the standard Python mechanisms.  See the `warnings`
        module in the Python standard library for more information.
        Defaults to False.
    filename : str, optional
        A filename, URL or other identifier to use in error messages.
        If *filename* is None and *source* is a string (i.e. a path),
        then *source* will be used as a filename for error messages.
        Therefore, *filename* is only required when source is a
        file-like object.

    Returns
    -------
    capabilities_file : `~pyvo.io.vosi.endpoint.CapabilitiesFile` object

    See also
    --------
    pyvo.io.vosi.exceptions : The exceptions this function may raise.
    """
    config = _pedantic_settings(pedantic)

    if filename is None and isinstance(source, str):
        config['filename'] = source
    else:
        config['filename'] = filename

    with iterparser.get_xml_iterator(
        source,
        _debug_python_based_parser=_debug_python_based_parser
    ) as iterator:
        return CapabilitiesFile(
            config=config, pos=(1, 1)).parse(iterator, config)


def parse_availability(source, *, pedantic=None, filename=None,
                       _debug_python_based_parser=False):
    """
    Parses a availability xml file (or file-like object), and returns a
    `~pyvo.io.vosi.endpoint.AvailabilityFile` object.

    Parameters
    ----------
    source : str or readable file-like object
        Path or file object containing a availability xml file.
    pedantic : bool, optional
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        the standard Python mechanisms.  See the `warnings`
        module in the Python standard library for more information.
        Defaults to False.
    filename : str, optional
        A filename, URL or other identifier to use in error messages.
        If *filename* is None and *source* is a string (i.e. a path),
        then *source* will be used as a filename for error messages.
        Therefore, *filename* is only required when source is a
        file-like object.

    Returns
    -------
    availability_file : `~pyvo.io.vosi.endpoint.AvailabilityFile` object

    See also
    --------
    pyvo.io.vosi.exceptions : The exceptions this function may raise.
    """
    config = _pedantic_settings(pedantic)

    if filename is None and isinstance(source, str):
        config['filename'] = source
    else:
        config['filename'] = filename

    with iterparser.get_xml_iterator(
        source,
        _debug_python_based_parser=_debug_python_based_parser
    ) as iterator:
        return AvailabilityFile(
            config=config, pos=(1, 1)).parse(iterator, config)


class TablesFile(Element):
    """
    TABLESET/TABLE element: represents an entire file.
    The keyword arguments correspond to setting members of the same
    name, documented below.
    """

    def __init__(self, *, config=None, pos=None, version="1.1"):
        Element.__init__(self, config, pos)

        self._tableset = None
        self._table = None

        self._ntables = None

        version = str(version)
        if version not in ("1.0", "1.1"):
            raise ValueError("'version' should be one of '1.0' or '1.1'")

        config['version'] = version
        self._version = version

    def __repr__(self):
        if self.table:
            return repr(self.table)
        elif self.tableset:
            return repr(self.tableset)
        else:
            return super().__repr__()

    @xmlattribute
    def version(self):
        """
        The version of the TableSet specification that the file uses.
        """
        return self._version

    @version.setter
    def version(self, version):
        version = str(version)
        if version not in ('1.0', '1.1'):
            raise ValueError(
                "pyvo.io.vosi.tables only supports VOSI versions 1.0 and 1.1")
        self._version = version

    @xmlelement
    def tableset(self):
        """
        The tableset. Must be a `TableSet` object.
        """
        return self._tableset

    @tableset.setter
    def tableset(self, tableset):
        self._tableset = tableset

    @tableset.adder
    def tableset(self, iterator, tag, data, config, pos):
        tableset = vs.TableSet(config, pos, 'tableset', **data)
        tableset.parse(iterator, config)
        self._tableset = tableset

    @xmlelement(cls=vs.VODataServiceTable)
    def table(self):
        """
        The `VODataServiceTable` root element if present.
        """
        return self._table

    @table.setter
    def table(self, table):
        self._table = table

    @property
    def ntables(self):
        """
        The number of tables in the file.
        """
        return self._ntables

    def parse(self, iterator, config):
        super().parse(iterator, config)

        if self.tableset is None and self.table is None:
            vo_raise(E07, config=config, pos=self._pos)

        self._version = config['version']
        if config['version'] not in ('1.0', '1.1'):
            vo_warn(W15, config=config, pos=self._pos)

        if self.table:
            if version_compare(config['version'], '1.1') < 0:
                vo_warn(W16, config=config, pos=self._pos)
            self._ntables = 1
        else:
            self._ntables = sum(
                len(schema.tables) for schema in self.tableset.schemas)

        return self

    def iter_tables(self):
        """
        Iterates over all tables in the VOSITables file in a "flat" way,
        ignoring the schemas.
        """
        if self.table:
            yield self.table
        else:
            for schema in self.tableset.schemas:
                yield from schema.tables

    def get_first_table(self):
        """
        When you parse table metadata for a single table here is only one table
        in the file, and that's all you need.
        This method returns that first table.
        """
        for table in self.iter_tables():
            return table
        raise IndexError("No table found in VOSITables file.")

    def get_table_by_name(self, name):
        """
        Looks up a table element by the given name.
        """
        for table in self.iter_tables():
            if table.name == name:
                return table
        raise KeyError("No table with name {} found".format(name))


class CapabilitiesFile(Element, HomogeneousList):
    """
    capabilities element: represents an entire file.
    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    def __init__(self, *, config=None, pos=None, _name='capabilities', **kwargs):
        Element.__init__(self, config=config, pos=pos, **kwargs)
        HomogeneousList.__init__(self, vr.Capability)

    @xmlelement(name='capability')
    def capabilities(self):
        """List of `~pyvo.io.vosi.voresource.Capability` objects"""
        return self

    @capabilities.adder
    def capabilities(self, iterator, tag, data, config, pos):
        capability = vr.Capability(config, pos, 'capability', **data)
        capability.parse(iterator, config)
        self.append(capability)

    def parse(self, iterator, config):
        for start, tag, data, pos in iterator:
            if start:
                if tag == "xml":
                    pass
                elif tag == "capabilities":
                    break
            else:
                vo_raise(E10, config=config, pos=pos)

        super().parse(iterator, config)

        return self


class AvailabilityFile(av.Availability):
    """
    availability element: represents an entire file.
    The keyword arguments correspond to setting members of the same
    name, documented below.
    """
    def parse(self, iterator, config):
        for start, tag, data, pos in iterator:
            if start:
                if tag == 'xml':
                    pass
                elif tag == 'availability':
                    break

        super().parse(iterator, config)

        return self
