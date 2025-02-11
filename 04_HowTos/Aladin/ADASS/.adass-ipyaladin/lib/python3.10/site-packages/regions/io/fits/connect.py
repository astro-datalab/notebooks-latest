# Licensed under a 3-clause BSD style license - see LICENSE.rst

import warnings

from astropy.io import fits

from regions.core import Region, Regions
from regions.core.registry import RegionsRegistry

__all__ = []


@RegionsRegistry.register(Region, 'identify', 'fits')
@RegionsRegistry.register(Regions, 'identify', 'fits')
def is_fits(methodname, filepath):
    """
    Identify a FITS region file.

    Parameters
    ----------
    methodname : {'read', 'write'}
        The method name called that needs auto-identification.

    filepath : str
        The path to the file.

    Returns
    -------
    result : bool
        Returns `True` if the given file is a FITS region file.
    """
    all_exten = ('.fits', '.fit', '.fts', '.fits.gz', '.fit.gz', '.fts.gz')
    exten = {'read': all_exten, 'write': all_exten[0:3]}

    if methodname == 'write':
        return filepath.lower().endswith(exten[methodname])

    elif methodname == 'read':
        if (isinstance(filepath, str)
                and filepath.lower().endswith(exten[methodname])):
            return True
        else:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    with fits.open(filepath):
                        return True
            except OSError:
                return False

    else:
        return False
