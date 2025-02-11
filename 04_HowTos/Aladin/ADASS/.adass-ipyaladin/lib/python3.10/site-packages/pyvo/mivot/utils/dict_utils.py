"""
Utility class to process dictionary.
"""
import json
import logging
from pyvo.mivot.utils.exceptions import MivotError
from pyvo.mivot.utils.json_encoder import MivotJsonEncoder


class DictUtils:
    """
    Static class implementing convenient operations on dictionaries.
    """

    @staticmethod
    def read_dict_from_file(filename, fatal=False):
        """
        Read a dictionary from a file and raise an exception if something goes wrong.
        Parameters:
        - filename (str): The filename. Any file that can be processed by json.load is accepted
        - fatal (bool): Triggers a system exit if True.
        Returns:
        - dict: The dictionary extracted from the file.
        Raises:
        - DataFormatException: If the file has an incorrect format.
        """
        try:
            logging.debug("Reading json from %s", filename)
            from collections import OrderedDict
            with open(filename, 'r') as file:
                return json.load(file, object_pairs_hook=OrderedDict)
        except Exception as exception:
            if fatal:
                raise MivotError("reading {}".format(filename))
            else:
                logging.error("{} reading {}".format(exception, filename))

    @staticmethod
    def _get_pretty_json(dictionary):
        """
        Return a pretty string representation of the dictionary.
        Parameters:
        - dictionary (dict): The dictionary.
        Returns:
        - str: A pretty string representation of the dictionary.
        """
        return json.dumps(dictionary,
                          indent=2,
                          cls=MivotJsonEncoder)

    @staticmethod
    def print_pretty_json(dictionary):
        """
        Print out a pretty string representation of the dictionary.
        Parameters:
        - dictionary (dict): The dictionary.
        """
        print(DictUtils._get_pretty_json(dictionary))
