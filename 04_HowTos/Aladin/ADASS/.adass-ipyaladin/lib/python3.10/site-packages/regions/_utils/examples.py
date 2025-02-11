# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.table import vstack as table_vstack
from astropy.utils import lazyproperty
from astropy.wcs import WCS

from regions.core.pixcoord import PixCoord

__all__ = ['make_example_dataset']


def make_example_dataset(data='simulated', config=None):
    """
    Make an example dataset.

    This is a factory function for ``ExampleDataset`` objects.

    The following config options are available (default values shown):

    * ``crval = 0, 0``
    * ``crpix = 180, 90``
    * ``cdelt = -1, 1``
    * ``shape = 180, 360``
    * ``ctype = 'GLON-AIT', 'GLAT-AIT'``

    Parameters
    ----------
    data : {'simulated', 'fermi'}
        The dataset to use.

    config : dict or None
        Configuration options.

    Returns
    -------
    dataset : ``ExampleDataset``
        An example dataset object.

    Examples
    --------
    Make an example dataset:

    >>> from regions import make_example_dataset
    >>> config = dict(crpix=(18, 9), cdelt=(-10, 10), shape=(18, 36))
    >>> dataset = make_example_dataset(data='simulated', config=config)

    Access properties of the ``dataset`` object:

    >>> dataset.source_table   # doctest: +IGNORE_OUTPUT
    >>> dataset.event_table   # doctest: +IGNORE_OUTPUT
    >>> dataset.wcs   # doctest: +IGNORE_OUTPUT
    >>> dataset.image   # doctest: +IGNORE_OUTPUT
    >>> dataset.hdu_list   # doctest: +IGNORE_OUTPUT
    """
    if data == 'simulated':
        return ExampleDatasetSimulated(config=config)
    elif data == 'fermi':
        return ExampleDatasetFermi(config=config)
    else:
        raise ValueError(f'Invalid selection data: {data}')


class ExampleDataset:
    """
    Base class for an example dataset.

    Parameters
    ----------
    config : dict or None
        Configuration options.
    """

    def __init__(self, config=None):
        config_default = dict()
        # These parameters are FITS WCS order: (x, y)
        config_default['crval'] = 0, 0
        config_default['crpix'] = 180, 90
        config_default['cdelt'] = -1, 1

        # This parameter is numpy order: (y, x)
        config_default['shape'] = 180, 360
        config_default['ctype'] = 'GLON-AIT', 'GLAT-AIT'

        # Callers can update the default config parameters they want.
        if config:
            config_default.update(config)

        self.config = config_default

    @lazyproperty
    def wcs(self):
        """
        World coordinate system (`~astropy.wcs.WCS`).
        """
        wcs = WCS(naxis=2)
        wcs.wcs.crval = self.config['crval']
        wcs.wcs.crpix = self.config['crpix']
        wcs.wcs.cdelt = self.config['cdelt']
        wcs.wcs.ctype = self.config['ctype']

        return wcs

    @lazyproperty
    def image(self):
        """
        Return a "counts" image (`~astropy.io.fits.ImageHDU`).
        """
        events = self.event_table
        skycoord = SkyCoord(events['GLON'], events['GLAT'], unit='deg',
                            frame='galactic')
        pixcoord = PixCoord.from_sky(skycoord=skycoord, wcs=self.wcs)

        shape = self.config['shape']
        bins = [np.arange(shape[0] + 1), np.arange(shape[1] + 1)]
        sample = np.vstack((pixcoord.y, pixcoord.x)).T
        data, _ = np.histogramdd(sample=sample, bins=bins)
        data = data.astype('float32')

        header = self.wcs.to_header()
        return fits.ImageHDU(data=data, header=header, name='image')

    @lazyproperty
    def hdu_list(self):
        """
        HDU list (`~astropy.io.fits.HDUList`).

        Different pieces collected together in a HDU list.

        This method makes it easy to write the example dataset to a FITS
        file with multiple HDUs.
        """
        hdu_list = fits.HDUList()

        hdu = _table_to_bintable(self.source_table)
        hdu.name = 'sources'
        hdu_list.append(hdu)

        hdu = _table_to_bintable(self.event_table)
        hdu.name = 'events'
        hdu_list.append(hdu)

        hdu = self.image
        hdu.name = 'image'
        hdu_list.append(hdu)

        return hdu_list


class ExampleDatasetSimulated(ExampleDataset):
    """
    Example simulated dataset.

    Similar to `ExampleDatasetFermi`, but simulated, not requiring any
    data files.
    """

    @lazyproperty
    def source_table(self):
        """
        Source table  (`~astropy.table.Table`).

        Columns: GLON, GLAT, COUNTS
        """
        table = Table()
        table['GLON'] = np.array([0, 45, 45], dtype='float32')
        table['GLAT'] = np.array([0, 0, 45], dtype='float32')
        table['COUNTS'] = np.array([100, 100, 100], dtype='int32')
        return table

    @lazyproperty
    def event_table(self):
        """
        Event table (`~astropy.table.Table`).

        Columns: GLON, GLAT, SOURCE_IDX
        """
        # Create event list table for each source
        tables = []
        for source in self.source_table:
            lon = source['GLON'] * np.ones(source['COUNTS'])
            lat = source['GLAT'] * np.ones(source['COUNTS'])

            # TODO: scatter positions assuming Gaussian PSF on the sky
            # using SkyOffsetFrame.
            # coord = SkyCoord(lon, lat, unit='deg', frame='galactic')

            table = Table()
            table['GLON'] = lon
            table['GLAT'] = lat
            table['SOURCE_IDX'] = source.index

            tables.append(table)

        # Stack all tables together
        table = table_vstack(tables)

        return table


class ExampleDatasetFermi(ExampleDataset):
    """
    Example real dataset using Fermi-LAT 2FHL source catalog and event
    list.

    When saving the HDU list to a FITS file, the file size is 748K, with
    the 65k EVENTS taking up most of the space.
    """

    @lazyproperty
    def source_table(self):
        """
        Source table  (`~astropy.table.Table`).

        Columns: GLON, GLAT, COUNTS
        """
        url = ('https://github.com/gammapy/gammapy-extra/raw/master/datasets/'
               'fermi_2fhl/gll_psch_v08.fit.gz')
        table = Table.read(url, hdu='2FHL Source Catalog')
        table.rename_column('Npred', 'COUNTS')
        table.keep_columns(['GLON', 'GLAT', 'COUNTS'])
        table.meta.clear()
        return table

    @lazyproperty
    def event_table(self):
        """
        Event table (`~astropy.table.Table`).

        Columns: GLON, GLAT
        """
        url = ('https://github.com/gammapy/gammapy-extra/raw/master/datasets/'
               'fermi_2fhl/2fhl_events.fits.gz')
        table = Table.read(url, hdu='EVENTS')
        table.rename_column('L', 'GLON')
        table.rename_column('B', 'GLAT')
        table.keep_columns(['GLON', 'GLAT'])
        table.meta.clear()
        return table


def _table_to_bintable(table):
    """
    Convert a `~astropy.table.Table` to a `astropy.io.fits.BinTable`.
    """
    data = table.as_array()
    header = fits.Header()
    header.update(table.meta)
    name = table.meta.pop('name', None)
    return fits.BinTableHDU(data, header, name=name)


if __name__ == '__main__':
    dataset = make_example_dataset(data='simulated')
    dataset.hdu_list.writeto('example-dataset-simulated.fits', overwrite=True)

    dataset = make_example_dataset(data='fermi')
    dataset.hdu_list.writeto('example-dataset-fermi.fits', overwrite=True)
