"""
Utility class to check that generated XML match reference elements.
Only used by the tests
"""
from pyvo.utils import activate_features

# Activate MIVOT for all tests
activate_features('MIVOT')

try:
    from defusedxml import ElementTree as etree
except ImportError:
    from xml.etree import ElementTree as etree
from pyvo.mivot.utils.xml_utils import XmlUtils

class XMLOutputChecker:
    """
    This class is used to compare XML outputs, ignoring whitespace differences.
    """
    def check_output(self, want, got):
        """
        Compare two XML outputs, ignoring whitespace differences.
        Parameters
        ----------
        want : str
            The expected XML output.
        got : str
            The actual XML output.
        Returns
        -------
        bool
            True if the two XML outputs are equal, False otherwise.
        """
        return self._format_xml(want.strip()) == self._format_xml(got.strip())

    def output_difference(self, want, got):
        """
        Return a string describing the differences between two XML outputs.
        Parameters
        ----------
        want : str
            The expected XML output.
        got : str
            The actual XML output.
        Returns
        -------
        str
            A string describing the differences between the two XML outputs.
        """
        return f"Diff:\n{self._format_xml(want)}\nvs.\n{self._format_xml(got)}"

    def _format_xml(self, xml_str):
        """
        Format an XML string.
        Parameters
        ----------
        xml_str : str
            The XML string to format.
        Returns
        -------
        str
            The formatted XML string.
        """
        return "\n".join(line.strip() for line in xml_str.splitlines())
    
    @staticmethod
    def xmltree_to_file(xmltree, file_path):
        """
        Write an XML tree to a file.
        Parameters
        ----------
        xmltree : ~`xml.etree.ElementTree.Element`
            The XML tree to write to the file.
        file_path : str
            The path to the output file.
        """
        with open(file_path, 'w') as output:
            output.write(XmlUtils.pretty_string(xmltree))


    @staticmethod
    def xmltree_from_file(file_path):
        """
        Parse an XML tree from a file.
        Parameters
        ----------
        file_path : str
            The path to the XML file.
        Returns
        -------
        ~`xml.etree.ElementTree.Element`
            The parsed XML tree.
        """
        return etree.parse(file_path)
    
    @staticmethod
    def assertXmltreeEquals(xmltree1, xmltree2):
        """
        Assert that two XML trees are equal.
        Parameters
        ----------
        xmltree1 : ~`xml.etree.ElementTree.Element`
            The first XML tree for comparison.
        xmltree2 : ~`xml.etree.ElementTree.Element`
            The second XML tree for comparison.
        """
        xml_str1 = etree.tostring(xmltree1).decode("utf-8")
        xml_str2 = etree.tostring(xmltree2).decode("utf-8")
        checker = XMLOutputChecker()
        assert checker.check_output(xml_str1, xml_str2, 0), f"XML trees differ:\n{xml_str1}\n---\n{xml_str2}"
        
    @staticmethod
    def assertXmltreeEqualsFile(xmltree1, xmltree2_file):
        """
        Assert that an XML tree is equal to the content of a file.
        Parameters
        ----------
        xmltree1 : ~`xml.etree.ElementTree.Element`
            The XML tree for comparison.
        xmltree2_file : str
            The path to the file containing the second XML tree.
        """
        xmltree2 = XMLOutputChecker.xmltree_from_file(xmltree2_file).getroot()
        xml_str1 = etree.tostring(xmltree1).decode("utf-8")
        xml_str2 = etree.tostring(xmltree2).decode("utf-8")
        checker = XMLOutputChecker()
        assert checker.check_output(xml_str1, xml_str2), f"XML trees differ:\n{xml_str1}\n---\n{xml_str2}"
