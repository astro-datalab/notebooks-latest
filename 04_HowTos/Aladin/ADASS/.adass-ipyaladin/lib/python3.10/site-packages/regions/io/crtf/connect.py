# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.data import get_readable_fileobj

from regions.core import Region, Regions
from regions.core.registry import RegionsRegistry

__all__ = []


@RegionsRegistry.register(Region, 'identify', 'crtf')
@RegionsRegistry.register(Regions, 'identify', 'crtf')
def is_crtf(methodname, filepath):
    """
    Identify a CRTF region file.

    Parameters
    ----------
    methodname : {'read', 'write'}
        The method name called that needs auto-identification.

    filepath : str
        The path to the file.

    Returns
    -------
    result : bool
        Returns `True` if the given file is a CRTF region file.
    """
    all_exten = ('.crtf', '.crtf.gz')
    exten = {'read': all_exten, 'write': all_exten[0]}

    if methodname == 'write':
        return filepath.lower().endswith(exten[methodname])

    elif methodname == 'read':
        if (isinstance(filepath, str)
                and filepath.lower().endswith(exten[methodname])):
            return True
        else:
            with get_readable_fileobj(filepath, encoding='binary') as fileobj:
                signature = '#CRTF'
                pos = fileobj.tell()
                sig = fileobj.read(len(signature))
                fileobj.seek(pos)
                return sig == signature or sig == signature.encode()

    else:
        return False
