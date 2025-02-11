# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module provides a Regions class.
"""

from regions.core.core import Region
from regions.core.registry import RegionsRegistry

__all__ = ['Regions']
__doctest_skip__ = ['Regions.read', 'Regions.write', 'Regions.parse',
                    'Regions.serialize']


class Regions:
    """
    Class to hold a list of `~regions.Region` objects.

    This class provides a unified I/O interface that supports reading,
    writing, parsing, and serializing many region data formats.

    Parameters
    ----------
    regions : list of `~regions.Region`
        The list of region objects.
    """

    def __init__(self, regions=(), /):
        if regions == ():
            regions = []
        for item in regions:
            if not isinstance(item, Region):
                raise TypeError('Input regions must be a list of Region '
                                'objects')
        self.regions = regions

    def __getitem__(self, index):
        newregions = self.regions[index]
        if isinstance(newregions, Region):  # one item
            return newregions
        else:
            newcls = object.__new__(self.__class__)
            newcls.regions = newregions
            return newcls

    def __repr__(self):
        cls_name = self.__class__.__name__
        return f'<{cls_name}({repr(self.regions)})>'

    def __str__(self):
        return str(self.regions)

    def __len__(self):
        return len(self.regions)

    def append(self, region):
        """
        Append the region to the end of the list of regions.

        Parameters
        ----------
        region : `~regions.Region`
            The region to append.
        """
        if not isinstance(region, Region):
            raise TypeError('Input region must be a Region object')
        self.regions.append(region)

    def extend(self, regions):
        """
        Extend the list of regions by appending elements from the input
        regions.

        Parameters
        ----------
        regions : `~regions.Regions` or list of `~regions.Region`
            A `~regions.Regions` object or a list of regions to include.
        """
        if isinstance(regions, Regions):
            self.regions.extend(regions.regions)
        else:
            for item in regions:
                if not isinstance(item, Region):
                    raise TypeError('Input regions must be a list of Region '
                                    'objects')
            self.regions.extend(regions)

    def insert(self, index, region):
        """
        Insert the region before index.

        Parameters
        ----------
        index : int
            The list index.
        region : `~regions.Region`
            The region to insert.
        """
        self.regions.insert(index, region)

    def reverse(self):
        """
        Reverse the list of regions in place.
        """
        self.regions.reverse()

    def pop(self, index=-1):
        """
        Remove and return the region at index.

        Parameters
        ----------
        index : int, optional
            The index of the region to remove.

        Returns
        -------
        result : `~regions.Region`
            The removed region.
        """
        return self.regions.pop(index)

    def copy(self):
        """
        Return a shallow copy of this object.
        """
        newcls = object.__new__(self.__class__)
        newcls.regions = self.regions.copy()
        return newcls

    @classmethod
    def get_formats(cls):
        """
        Get the registered I/O formats as a Table.
        """
        return RegionsRegistry.get_formats(cls)

    @classmethod
    def read(cls, filename, format=None, cache=False, **kwargs):
        """
        Read and parse a region file and return as a Regions object.

        This method allows reading a file in many supported data
        formats, e.g.,::

            >>> from regions import Regions
            >>> reg1 = Regions.read('regions.reg', format='ds9')
            >>> reg2 = Regions.read('regions.crtf', format='crtf')
            >>> reg3 = Regions.read('regions.fits', format='fits')

        A list of the available formats for `~regions.Regions` is
        available using::

            >>> Regions.get_formats()

        Parameters
        ----------
        filename : str
            The filename or URL of the file to read.

        format : str, optional
            The file format specifier.

        cache : bool or 'update', optional
            Whether to cache the contents of remote URLs. If 'update',
            check the remote URL for a new version but store the result
            in the cache.

        **kwargs : dict, optional
            Keyword arguments passed to the data reader.

        Returns
        -------
        result : `~regions.Regions`
            A `~regions.Regions` object containing the file contents.
        """
        return RegionsRegistry.read(filename, cls, format=format,
                                    cache=cache, **kwargs)

    @classmethod
    def parse(cls, data, format=None, **kwargs):
        """
        Parse a region string or table and return as a Regions object.

        This method allows parsing region data in many supported data
        formats, e.g.,::

            >>> from regions import Regions
            >>> reg1 = Regions.parse(regions_str, format='ds9')
            >>> reg2 = Regions.parse(regions_str, format='crtf')
            >>> reg3 = Regions.parse(regions_tbl, format='fits')

        A list of the available formats for `~regions.Regions` is
        available using::

            >>> Regions.get_formats()

        Parameters
        ----------
        data : str or `~astropy.table.Table`
            The region data to parse.

        format : str, optional
            The file format specifier.

        **kwargs : dict, optional
            Keyword arguments passed to the data parser.

        Returns
        -------
        result : `~regions.Regions`
            A `~regions.Regions` object containing the data contents.
        """
        return RegionsRegistry.parse(data, cls, format=format,
                                     **kwargs)

    def write(self, filename, format=None, overwrite=False, **kwargs):
        """
        Write the regions to a region file in the specified format.

        This method allows writing a file in many supported data
        formats, e.g.,::

            >>> from regions import Regions
            >>> reg = Regions.read('regions.reg', format='ds9')
            >>> reg.write('new_regions.reg', format='ds9')
            >>> reg.write('new_regions.crtf', format='crtf')
            >>> reg.write('new_regions.fits', format='fits')

        A list of the available formats for `~regions.Regions` is
        available using::

            >>> Regions.get_formats()

        Parameters
        ----------
        filename : str
            The filename or URL of the file to write.

        format : str, optional
            The file format specifier.

        overwrite : bool, optional
            If True, overwrite the output file if it exists. Raises an
            `OSError` if False and the output file exists. Default is
            False.

        **kwargs : dict, optional
            Keyword arguments passed to the data writer.
        """
        return RegionsRegistry.write(self.regions, filename,
                                     self.__class__, format=format,
                                     overwrite=overwrite, **kwargs)

    def serialize(self, format=None, **kwargs):
        """
        Serialize the regions to a region string or table.

        This method allows serializing regions in many supported data
        formats, e.g.,::

            >>> from regions import Regions
            >>> reg = Regions.read('regions.reg', format='ds9')
            >>> reg1_str = reg.serialize(format='ds9')
            >>> reg2_str = reg.serialize(format='crtf')
            >>> reg3_tbl = reg.serialize(format='fits')

        A list of the available formats for `~regions.Regions` is
        available using::

            >>> Regions.get_formats()

        Parameters
        ----------
        format : str, optional
            The file format specifier.

        **kwargs : dict, optional
            Keyword arguments passed to the data serializer.
        """
        return RegionsRegistry.serialize(self.regions, self.__class__,
                                         format=format, **kwargs)
