"""
Utility class to process XML.
"""
from pyvo.mivot.utils.xpath_utils import XPath
import xml.etree.ElementTree as ET
from pyvo.mivot.utils.vocabulary import Constant
from pyvo.mivot.utils.vocabulary import Att
from pyvo.mivot.utils.exceptions import MivotError


class XmlUtils:
    """
    Static class implementing convenient operations on XML
    """
    @staticmethod
    def pretty_print(xmltree):
        """
        Pretty print an XML tree.
        Parameters
        ----------
        xmltree (~`xml.etree.ElementTree.Element`): XML tree to pretty print.
        """
        print(XmlUtils.pretty_string(xmltree))

    @staticmethod
    def pretty_string(xmltree):
        """
        Return a pretty string representation of an XML tree.
        Parameters
        ----------
        xmltree (~`xml.etree.ElementTree.Element`): XML tree to convert to a pretty string.
        Returns
        -------
        str: The pretty string representation of the XML tree.
        """
        if hasattr(xmltree, 'getroot'):
            XmlUtils.indent(xmltree.getroot())
            new_xml = ET.tostring(xmltree.getroot(), encoding='unicode')
        else:
            XmlUtils.indent(xmltree)
            new_xml = ET.tostring(xmltree, encoding='unicode')
        return new_xml.replace("ns0:", "")

    @staticmethod
    def indent(elem, level=0):
        """
        Indent an XML tree.
        Parameters
        ----------
        elem (~`xml.etree.ElementTree.Element`): XML tree to indent.
        level (int): level of indentation.
        Returns
        -------
        ~`xml.etree.ElementTree.Element`
            The indented XML tree.
        """
        i = "\n" + level * "  "
        j = "\n" + (level - 1) * "  "
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = i + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = i
            for subelem in elem:
                XmlUtils.indent(subelem, level + 1)
            if not elem.tail or not elem.tail.strip():
                elem.tail = j
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = j
        return elem

    @staticmethod
    def add_column_indices(mapping_block, index_map):
        """
        Add column ranks to attributes having a ref.
        Using ranks allows identifying columns even when NumPy arrays have been serialized as [].
        Parameters
        ----------
        mapping_block : ~`xml.etree.ElementTree.Element`
            The XML mapping block.
        index_map : dict
            A dictionary mapping ref values to column indices.
        """
        for ele in XPath.x_path(mapping_block, ".//ATTRIBUTE"):
            attr_ref = ele.get(Att.ref)
            if attr_ref is not None and attr_ref != Constant.NOT_SET:
                field_desc = None
                if attr_ref in index_map:
                    field_desc = index_map[attr_ref]
                else:
                    for _, value in index_map.items():
                        if value["ID"] == attr_ref:
                            field_desc = value
                            break
                if not field_desc:
                    if not ele.get(Att.value):
                        raise MivotError(
                            f"Attribute {ele.get(Att.dmrole)} can not be set:"
                            f" references a non existing column: {attr_ref} "
                            f"and has no default value")
                    else:
                        ele.attrib.pop(Att.ref, None)
                if field_desc:
                    ele.attrib[Constant.COL_INDEX] = str(field_desc["indx"])
                    if field_desc["ID"] != attr_ref:
                        ele.set(Att.ref, field_desc["ID"])

    @staticmethod
    def add_column_units(mapping_block, unit_map):
        """
        Add field units to attributes having a ref.
        Used for performing unit conversions.
        Parameters
        ----------
        mapping_block : ~`xml.etree.ElementTree.Element`
            The XML mapping block.
        unit_map : dict
            A dictionary mapping ref values to units.
        """
        for ele in XPath.x_path(mapping_block, ".//ATTRIBUTE"):
            ref = ele.get(Att.ref)
            if ref is not None and ref != Constant.NOT_SET:
                unit = None
                if ref in unit_map:
                    unit = unit_map[ref].__str__()
                else:
                    unit = ""
                ele.attrib[Constant.FIELD_UNIT] = unit
