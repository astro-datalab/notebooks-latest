# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
XMLViewer provides several getters on XML instances built by
 `pyvo.mivot.viewer.mivot_viewer`.
"""
from pyvo.mivot.utils.exceptions import MivotError
from pyvo.mivot.utils.xpath_utils import XPath
from pyvo.utils.prototype import prototype_feature


@prototype_feature('MIVOT')
class XMLViewer:
    """
    The XMLViewer  is used by `~pyvo.mivot.viewer.mivot_viewer`
    to extract from the XML serialization of the model,
    elements that will be used to build the dictionary from which
    the Python class holding the mapped model will be generated.
    """
    def __init__(self, xml_view):
        self._xml_view = xml_view

    @property
    def view(self):
        """
        getter returning the XML model view

        returns
        -------
            XML model view to be parsed
            by different methods
        """
        return self._xml_view

    def get_instance_by_role(self, dmrole, all_instances=False):
        """
        If all_instances is False, return the first INSTANCE matching with @dmrole.
        If all_instances is True, return a list of all instances matching with @dmrole.
        Parameters
        ----------
        dmrole : str
            The @dmrole to look for.
        all_instances : bool, optional
            If True, returns a list of all instances, otherwise returns the first instance.
            Default is False.
        Returns
        -------
        Union[`xml.etree.ElementTree.Element`, List[`xml.etree.ElementTree.Element`], None]
            If all_instances is False, returns the instance matching with @dmrole.
            If all_instances is True, returns a list of all instances matching with @dmrole.
            If no matching instance is found, returns None.
        Raises
        ------
        MivotElementNotFound
            If dmrole is not found.
        """
        instances = XPath.select_elements_by_atttribute(
            self._xml_view,
            "INSTANCE",
            "dmrole",
            dmrole)

        if len(instances) == 0:
            raise MivotError(
                f"Cannot find dmrole {dmrole} in any instances of the VOTable")

        if all_instances is False:
            return instances[0]
        else:
            return instances

    def get_instance_by_type(self, dmtype, all_instances=False):
        """
        Return the instance matching with @dmtype.
        If all_instances is False, returns the first INSTANCE matching with @dmtype.
        If all_instances is True, returns a list of all instances matching with @dmtype.
        Parameters
        ----------
        dmtype : str
            The @dmtype to look for.
        all : bool, optional
            If True, returns a list of all instances, otherwise returns the first instance.
            Default is False.
        Returns
        -------
        Union[~`xml.etree.ElementTree.Element`, List[~`xml.etree.ElementTree.Element`], None]
            If all_instances is False, returns the instance matching with @dmtype.
            If all_instances is True, returns a list of all instances matching with @dmtype.
            If no matching instance is found, returns None.
        Raises
        ------
        MivotElementNotFound
            If dmtype is not found.
        """
        instances = XPath.select_elements_by_atttribute(
            self._xml_view,
            "INSTANCE",
            "dmtype",
            dmtype)

        if len(instances) == 0:
            raise MivotError(
                f"Cannot find dmtype {dmtype} in any instances of the VOTable")

        if all_instances is False:
            return instances[0]
        else:
            return instances

    def get_collection_by_role(self, dmrole, all_instances=False):
        """
        Return the collection matching with @dmrole.
        If all_instances is False, returns the first COLLECTION matching with @dmrole.
        If all_instances is True, returns a list of all COLLECTION matching with @dmrole.
        Parameters
        ----------
        dmrole : str
            The @dmrole to look for.
        all_instances : bool, optional
            If True, returns a list of all COLLECTION, otherwise returns the first COLLECTION.
            Default is False.
        Returns
        -------
        Union[~`xml.etree.ElementTree.Element`, List[~`xml.etree.ElementTree.Element`], None]
            If all_instances is False, returns the collection matching with @dmrole.
            If all_instances is True, returns a list of all collections matching with @dmrole.
            If no matching collection is found, returns None.
        Raises
        ------
        MivotElementNotFound
            If dmrole is not found.
        """
        collections = XPath.select_elements_by_atttribute(
            self._xml_view,
            "COLLECTION",
            "dmrole",
            dmrole)

        if len(collections) == 0:
            raise MivotError(
                f"Cannot find dmrole {dmrole} in any collections of the VOTable")

        if all_instances is False:
            return collections[0]
        else:
            return collections
