# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Class that provides multiple getters on VOTable RESOURCE elements.
"""
from pyvo.mivot.utils.vocabulary import Constant
from pyvo.utils.prototype import prototype_feature


@prototype_feature('MIVOT')
class ResourceSeeker:
    """
    This class provides multiple getters on resource tables.
    Some methods are simple wrappers for external tools in order to have all
    the search functions on RESOURCE gathered in within a single namespace.
    """
    def __init__(self, resource):
        """
        Constructor
        Parameters
        ----------
        resource (astropy.votable.Resource): The resource object to be queried.
        """
        self._resource = resource

    def get_table_ids(self):
        """
        Return the list of table ids.
        Only resource children are considered.
        The @ID is first searched and then the @name, and finally 'AnonymousTable' is taken.
        Returns
        -------
        list of str: table ids.
        """
        ids_found = []
        for table in self._resource.tables:
            if table.ID is not None:
                ids_found.append(table.ID)
            elif table.name is not None:
                ids_found.append(table.name)
            else:
                ids_found.append(Constant.ANONYMOUS_TABLE)
        return ids_found

    def get_table(self, table_name_or_id):
        """
        Return the table matching table_name first by ID and then by name.
        Parameters
        ----------
        table_name_or_id (str): Name or id of the table to get.
        Returns
        -------
        ~astropy.votable.table: table matching the table_name.
        """
        if table_name_or_id == Constant.FIRST_TABLE:
            return self._resource.tables[0]
        for table in self._resource.tables:
            if (table_name_or_id is None or table.name == table_name_or_id
                    or table.ID == table_name_or_id):
                return table
        return None

    def get_params(self):
        """
        Return the VOTable PARAMS.
        Returns
        -------
        ~astropy.votable.Resource.params: VOTable PARAMS.
        """
        return self._resource.params

    def get_id_index_mapping(self, table_name):
        """
        Build an index binding column number with field id.
        Parameters
        ----------
        table_name (str): Name of the table.
        Returns
        -------
        dict: dictionary mapping field id to column number: {name: {ID, ref, indx}...}
        """
        column_index = {}
        table = self.get_table(table_name)
        indx = 0
        for field in table.fields:
            field_desc = {}
            if field.ID is not None:
                field_desc["ID"] = field.ID
            if field.ref is not None:
                field_desc["ref"] = field.ref
            field_desc["indx"] = indx
            if "ID" not in field_desc:
                field_desc["ID"] = field.name
            column_index[field.name] = field_desc
            indx += 1

        return column_index

    def get_id_unit_mapping(self, table_name):
        """
        Build an index binding field unit with field id.
        Parameters
        ----------
        table_name (str): Name of the table.
        Returns
        -------
        dict: A dictionary mapping field id to field unit {ID1: unit, name1: unit, ref1: unit ...}.
        """
        unit_index = {}
        table = self.get_table(table_name)
        for field in table.fields:
            unit = field.unit
            if field.ID is not None:
                unit_index[field.ID] = unit
            elif field.name is not None:
                unit_index[field.name] = unit
            elif field.ref is not None:
                unit_index[field.ref] = unit
        return unit_index
