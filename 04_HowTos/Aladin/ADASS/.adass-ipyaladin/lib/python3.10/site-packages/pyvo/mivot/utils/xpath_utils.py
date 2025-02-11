"""
Utility class performing XPath queries on XML trees.
"""


class XPath:
    """
    Static class use to perform XPath queries on XML trees.
    """
    @staticmethod
    def x_path(etree, path):
        """
        Return all the elements of the XML tree that match the given XPath query.
        Parameters
        ----------
        etree : ~`xml.etree.ElementTree.Element`
            The XML tree to query.
        path : str
            The XPath query to perform.
        Returns
        -------
        list
            The list of all the elements of the XML tree that match the given XPath query.
        """
        return etree.findall(path)

    @staticmethod
    def x_path_contains(etree, path, key, value):
        """
        Return all the elements of the XML tree that match the given
        XPath query with a given attribute containing a given value.
        Example of a path: ".//INSTANCE[contains(@dmtype,'dmtype_pattern')]"
        Parameters
        ----------
        etree : ~`xml.etree.ElementTree.Element`
            The XML tree to query.
        path : str
            The XPath query to perform.
        key : str
            The attribute to look for.
        value : str
            The value to look for.
        Returns
        -------
        list
            The list of all the elements of the XML tree that match the
            given XPath query with a given attribute containing a given value.
        """
        result = []
        xset = etree.findall(path)
        for ele in xset:
            if value in ele.get(key):
                result.append(ele)
        return result

    @staticmethod
    def x_path_startwith(etree, path):
        """
        Return all the elements of the XML tree that match the given
        XPath query with a given Tag starting with a given value.
        Example of a path: ".//*[starts-with(name(), 'REFERENCE_')]"
        This function is only used in the static reference resolver to find all the REFERENCEs and JOINs.
        It adds a counter to the path to find all the REFERENCEs and JOINs.
        Parameters
        ----------
        etree : `xml.etree.ElementTree.Element`
            The XML tree to query.
        path : str
            The XPath query to perform.
        Returns
        -------
        list
            The list of all the elements of the XML tree that match the given XPath query.
        """
        cnt = 1
        result = []
        run = True
        while run:
            if etree.find(path + str(cnt)) is not None:
                result.append(etree.find(path + str(cnt)))
                cnt += 1
            else:
                run = False
        return result

    @staticmethod
    def select_elements_by_atttribute(etree, element, attribute, attribute_value):
        """
        Select all xml elements named 'element' having @attribute=attribute_value
        Parameters
        ----------
        etree : `xml.etree.ElementTree.Element`
            The XML tree to query.
        element: str
            name of the searched elements
        attribute: str
            name of the element attribute used for the selection
        attribute_value: str
            attribute value of the selected elements
        Returns
        -------
        list
            The list of all the elements of the XML tree that match the given XPath query.
        """
        return XPath.x_path(etree,
                            f'.//{element}[@{attribute}="{attribute_value}"]'
                            )
