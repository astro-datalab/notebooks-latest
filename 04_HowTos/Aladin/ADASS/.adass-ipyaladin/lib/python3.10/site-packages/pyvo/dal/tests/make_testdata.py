from pathlib import Path

import numpy as np

from astropy.table import Table
from astropy.io.votable.tree import VOTableFile, Info
from astropy.io.fits import ImageHDU


def _votablefile():
    table = Table([
        [23, 42, 1337],
        [b'Illuminatus', b"Don't panic, and always carry a towel", b'Elite']
    ], names=('1', '2'))

    table['1'].meta['ucd'] = 'foo;bar'

    table['2'].meta['utype'] = 'foobar'

    votable_file = VOTableFile.from_table(table)

    info = Info(name='QUERY_STATUS', value='OK')
    info.content = 'OK'
    votable_file.resources[0].infos.append(info)

    return votable_file


def votablefile():
    votable_file = _votablefile()
    return votable_file


def votablefile_errorstatus():
    votable_file = _votablefile()

    info = Info(name='QUERY_STATUS', value='ERROR')
    info.content = 'ERROR'
    votable_file.resources[0].infos[0] = info

    return votable_file


def votablefile_overflowstatus():
    votable_file = _votablefile()

    info_ok = Info(name='QUERY_STATUS', value='OK')
    info_overflow = Info(name='QUERY_STATUS', value='OVERFLOW')
    votable_file.resources[0].infos[0] = info_ok
    votable_file.resources[0].infos.append(info_overflow)

    return votable_file


def votablefile_missingtable():
    votable_file = _votablefile()
    del votable_file.resources[0].tables[0]
    return votable_file


def votablefile_missingresource():
    votable_file = _votablefile()
    del votable_file.resources[0]
    return votable_file


def votablefile_missingcolumns():
    votable_file = _votablefile()
    del votable_file.resources[0].tables[0].fields[:]
    return votable_file


def votablefile_firstresource():
    votable_file = _votablefile()
    votable_file.resources[0]._type = 'results'
    return votable_file


def votablefile_tableinfo():
    votable_file = _votablefile()
    votable_file.resources[0].tables[0].infos[:] = (
        votable_file.resources[0].infos[:])

    del votable_file.resources[0].infos[:]

    return votable_file


def votablefile_rootinfo():
    votable_file = _votablefile()
    votable_file.infos[:] = (
        votable_file.resources[0].infos[:])

    del votable_file.resources[0].infos[:]

    return votable_file


def votablefile_dataset():
    table = Table([
        [
            'image/fits',
            'application/x-votable+xml',
            'application/x-votable+xml;content=datalink'
        ],
        [
            b'http://example.com/querydata/image.fits',
            b'http://example.com/querydata/votable.xml',
            b'http://example.com/querydata/votable-datalink.xml'
        ]
    ], names=('dataformat', 'dataurl'))

    table['dataformat'].meta['ucd'] = 'meta.code.mime'
    table['dataurl'].meta['utype'] = 'Access.Reference'
    table['dataurl'].meta['ucd'] = 'meta.dataset;meta.ref.url'

    votable_file = VOTableFile.from_table(table)

    info = Info(name='QUERY_STATUS', value='OK')
    info.content = 'OK'
    votable_file.resources[0].infos.append(info)

    return votable_file


def dataset_fits():
    hdu = ImageHDU(np.random.random((256, 256)))
    return hdu


def main():
    dirname = Path(__file__).parent / 'data'

    votablefile().to_xml(
        str(dirname / 'query/basic.xml'), tabledata_format='tabledata')

    votablefile_errorstatus().to_xml(
        str(dirname / 'query/errorstatus.xml'), tabledata_format='tabledata')

    votablefile_overflowstatus().to_xml(
        str(dirname / 'query/overflowstatus.xml'), tabledata_format='tabledata')

    votablefile_missingtable().to_xml(
        str(dirname / 'query/missingtable.xml'), tabledata_format='tabledata')

    votablefile_missingresource().to_xml(
        str(dirname / 'query/missingresource.xml'),
        tabledata_format='tabledata')

    votablefile_missingcolumns().to_xml(
        str(dirname / 'query/missingcolumns.xml'),
        tabledata_format='tabledata')

    votablefile_firstresource().to_xml(
        str(dirname / 'query/firstresource.xml'), tabledata_format='tabledata')

    votablefile_tableinfo().to_xml(
        str(dirname / 'query/tableinfo.xml'), tabledata_format='tabledata')

    votablefile_rootinfo().to_xml(
        str(dirname / 'query/rootinfo.xml'), tabledata_format='tabledata')

    votablefile_dataset().to_xml(
        str(dirname / 'query/dataset.xml'), tabledata_format='tabledata')

    dataset_fits().writeto(
        str(dirname / 'querydata/image.fits'), overwrite=True)


if __name__ == '__main__':
    main()
