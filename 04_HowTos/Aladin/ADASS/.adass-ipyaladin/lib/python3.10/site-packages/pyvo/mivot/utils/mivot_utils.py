'''
Some utilities making easier the transformation of Mivot elements into dictionary components.
These dictionaries are used to generate ``MivotInstance`` objects
'''
import numpy


class MivotUtils(object):
    @staticmethod
    def xml_to_dict(element):
        """
        Recursively create a nested dictionary from the XML tree structure, preserving the hierarchy.
        Each object in the dictionary is represented by a new dictionary with dmrole: {}.
        The processing of elements depends on the tag:
         - For INSTANCE, a new dictionary is created.
         - For COLLECTION, a list is created.
         - For ATTRIBUTE, a leaf is created in the tree structure with dmtype, dmrole, value, unit, and ref.
        Parameters
        ----------
        element (~`xml.etree.ElementTree.Element`) : The XML element to convert to a dictionary.
        Returns
        -------
        dict: The nested dictionary representing the XML tree structure.
        """
        dict_result = {}
        for key, value in element.attrib.items():
            dict_result[key] = value
        for child in element:
            dmrole = child.get("dmrole")
            if child.tag == "ATTRIBUTE":
                dict_result[dmrole] = MivotUtils.attribute_to_dict(child)
            elif child.tag == "INSTANCE":  # INSTANCE is recursively well managed by the function _to_dict
                dict_result[dmrole] = MivotUtils.xml_to_dict(child)
            elif child.tag == "COLLECTION":
                dict_result[dmrole] = MivotUtils.collection_to_dict(child)
        return dict_result

    @staticmethod
    def attribute_to_dict(child):
        """
        Convert an ATTRIBUTE element to a dictionary.
        ATTRIBUTE is always a leaf, so it is not recursive.
        Parameters
        ----------
        child (~`xml.etree.ElementTree.Element`): ATTRIBUTE XML element to convert.
        Returns
        -------
        dict: A dictionary representing the ATTRIBUTE element with keys:
              'dmtype', 'dmrole', 'value', 'unit', and 'ref'.
        """
        attribute = {}
        if child.get('dmtype') is not None:
            attribute['dmtype'] = child.get("dmtype")
        if child.get("value") is not None:
            attribute['value'] = MivotUtils.cast_type_value(child.get("value"), child.get("dmtype"))
        else:
            attribute['value'] = None
        if child.get("unit") is not None:
            attribute['unit'] = child.get("unit")
        else:
            attribute['unit'] = None
        if child.get("ref") is not None:
            attribute['ref'] = child.get("ref")
        else:
            attribute['ref'] = None
        return attribute

    @staticmethod
    def collection_to_dict(child):
        """
        Convert a COLLECTION element to a list of dictionaries.
        COLLECTION is always represented as a list, and each element of the COLLECTION is added to the list.
        Parameters
        ----------
        child (`~`xml.etree.ElementTree.Element``):  COLLECTION XML element to convert.
        Returns
        -------
        list: list of dictionaries representing the elements of the COLLECTION.
        """
        collection_items = []
        for child_coll in child:
            collection_items.append(MivotUtils.xml_to_dict(child_coll))
        return collection_items

    @staticmethod
    def cast_type_value(value, dmtype):
        """
        Cast the value of an ATTRIBUTE based on its dmtype.
        As the type of ATTRIBUTE values returned in the dictionary is string by default,
        this function is used to cast them based on their dmtype.
        Parameters
        ----------
        value (str): value of the ATTRIBUTE.
        dmtype (str): dmtype of the ATTRIBUTE.
        Returns
        -------
        Union[bool, float, str, None]
            The cast value based on the dmtype.
        """
        if type(value) is numpy.float32 or type(value) is numpy.float64:
            return float(value)
        lower_dmtype = dmtype.lower()
        if isinstance(value, str):
            lower_value = value.lower()
        else:
            lower_value = value
        if "bool" in lower_dmtype:
            if value == "1" or lower_value == "true" or lower_value:
                return True
            else:
                return False
        elif lower_value in ('notset', 'noset', 'null', 'none', 'nan') or value is None:
            return None
        elif (isinstance(value, numpy.ndarray) or isinstance(value, numpy.ma.core.MaskedConstant)
              or value == '--'):
            return None
        elif "real" in lower_dmtype or "double" in lower_dmtype or "float" in lower_dmtype:
            return float(value)
        else:
            return value
