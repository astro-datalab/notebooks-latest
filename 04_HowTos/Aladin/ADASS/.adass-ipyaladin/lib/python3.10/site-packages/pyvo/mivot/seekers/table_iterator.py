# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Iterator for table rows.
"""
from pyvo.utils.prototype import prototype_feature


@prototype_feature('MIVOT')
class TableIterator:
    """
    Simple wrapper iterating over table rows.
    Some methods are simple wrappers for external tools in order to have all
    the search functions on TABLE gathered in within a single namespace.
    """
    def __init__(self, name, data_table):
        """
        Constructor of the TableIterator class.
        Parameters
        ----------
        name (str): Table name (not really used).
        data_table (~numpy.ndarray): Numpy table returned by `~astropy.votable`.
        """
        self.name = name
        self.data_table = data_table
        self.last_row = None
        self.iter = None
        # not used yet
        self.row_filter = None

    def get_next_row(self):
        """
        Return the next Numpy row or None.
        The end of table exception usually returned by Numpy is trapped.
        """
        # The iterator is set at the first iteration
        if self.iter is None:
            self.iter = iter(self.data_table)
        try:
            while True:
                row = next(self.iter)
                if row is not None:
                    if self.row_filter is None or self.row_filter.row_match(row):
                        self.last_row = row
                        return row
                else:
                    return None
        except StopIteration:
            return None

    def rewind(self):
        """
        Set the pointer on the table-top, destroys the iterator actually.
        """
        self.iter = None
