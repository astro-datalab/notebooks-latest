"""
Miscellenaneous utilities for writing tests.
"""

from astropy.io.votable import tree
from pyvo.dal import query as dalquery


try:
    TABLE_ELEMENT = tree.TableElement
except AttributeError:
    TABLE_ELEMENT = tree.Table


def create_votable(field_descs, records):
    """returns a VOTableFile with a a single table containing records,
    described by field_descs.
    """
    votable = tree.VOTableFile()
    resource = tree.Resource(type="results")
    votable.resources.append(resource)
    table = TABLE_ELEMENT(votable)
    resource.tables.append(table)
    table.fields.extend(
        tree.Field(votable, **desc) for desc in field_descs)
    table.create_arrays(len(records))

    for index, rec in enumerate(records):
        table.array[index] = rec

    return votable


def create_dalresults(
        field_descs,
        records,
        *,
        resultsClass=dalquery.DALResults):
    """returns a DALResults instance for a query returning records
    described by field_descs.

    The arguments are as for create_votable.
    """

    return resultsClass(
        create_votable(field_descs, records),
        url="http://testing.pyvo/test-url")
