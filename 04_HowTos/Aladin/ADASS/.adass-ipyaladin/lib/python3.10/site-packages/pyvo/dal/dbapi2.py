# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
An implementation of the Database API v2.0 interface to DAL VOTable responses.
This only supports read-only access.
"""
from .query import Iter

apilevel = "2.0"
threadsafety = 2
paramstyle = "n/a"

__all__ = "STRING BINARY NUMBER DATETIME ROWID".split()


class Error(Exception):
    """
    DB-API base exception
    """
    pass


class Warning(Exception):
    """
    DB-API warning
    """
    pass


class InterfaceError(Error):
    """
    DB-API exception indicating an error related to the database interface
    rather than the database itself.
    """
    pass


class DatabaseError(Error):
    """
    DB-API exception indicating an error related to the database.
    """
    pass


class DataError(DatabaseError):
    """
    DB-API exception indicating an error while processing data (e.g. divide
    by zero, numeric value out-of-range, etc.)
    """
    pass


class OperationalError(DatabaseError):
    """
    DB-API exception indicating an error related to the database's operation
    and not necessarily under the control of the programmer.
    """
    pass


class IntegrityError(DatabaseError):
    """
    DB-API exception indicating an inconsistancy in the database integrity.
    """
    pass


class InternalError(DatabaseError):
    """
    DB-API exception indicating an internal error that might indicate that
    a connection or cursor is no longer valid.
    """
    pass


class ProgrammingError(DatabaseError):
    """
    DB-API exception indicating an erroneous request (e.g. column not found)
    """
    pass


class NotSupportedError(DatabaseError):
    """
    DB-API exception indicating a request is not supported
    """
    pass


class TypeObject:
    def __init__(self, *values):
        self._values = values

    @property
    def id(self):
        return self._values[0]

    def __eq__(self, other):
        if not isinstance(other, TypeObject):
            return False
        if other.id in self._values:
            return True
        return self.id in other._values

    def __ne__(self, other):
        return not self.__eq__(other)


STRING = TypeObject(0)
BINARY = TypeObject(1)
NUMBER = TypeObject(2)
DATETIME = TypeObject(3, STRING.id)
ROWID = TypeObject(4, NUMBER.id)


def connect(source):
    raise NotSupportedError("Connection objects not supported")


class Cursor(Iter):
    """
    A class used to walk through a query response table row by row,
    accessing the contents of each record (row) of the table.  This class
    implements the Python Database API.
    """

    def __init__(self, results):
        """Create a cursor instance.  The constructor is not typically called
        by directly applications; rather an instance is obtained from calling a
        DalQuery's execute().
        """
        super().__init__(results)
        self._description = self._mkdesc()
        self._rowcount = len(self.resultset)
        self._arraysize = 1

    def _mkdesc(self):
        flds = self.resultset.fieldnames
        out = []
        for name in flds:
            fld = self.resultset.getdesc(name)
            typ = STRING
            if fld.datatype in ("short", "int", "long", "float", "double",
                                "floatComplex", "doubleComplex", "boolean"):
                typ = NUMBER
            elif fld.datatype in "char unicodeChar unsignedByte".split():
                typ = STRING

            out.append((name, typ))

        return tuple(out)

    @property
    def description(self):
        """
        a read-only sequence of 2-item seqences.  Each seqence describes
        a column in the results, giving its name and type_code.
        """
        return self._description

    @property
    def rowcount(self):
        """
        the number of rows in the result (read-only)
        """
        return self._rowcount

    @property
    def arraysize(self):
        """
        the number of rows that will be returned by returned by a call to
        fetchmany().  This defaults to 1, but can be changed.
        """
        return self._arraysize

    @arraysize.setter
    def arraysize(self, value):
        if not value:
            value = 1
        self._arraysize = value

    def infos(self):
        """Return any INFO elements in the VOTable as a dictionary.

        Returns
        -------
        dict :
            A dictionary with each element corresponding to a single INFO,
            representing the INFO as a name:value pair.
        """
        return self.resultset._infos

    def fetchone(self):
        """Return the next row of the query response table.

        Returns
        -------
        tuple :
            The response is a tuple wherein each element is the value of the
            corresponding table field.
        """
        try:
            rec = self.next()
            out = []
            for name in self.resultset.fieldnames:
                out.append(rec[name])
            return out
        except StopIteration:
            return None

    def fetchmany(self, size=None):
        """Fetch the next block of rows from the query result.

        Parameters
        ----------
        size : int
            The number of rows to return (default: cursor.arraysize).

        Returns
        -------
        list of tuples :
            A list of tuples, one per row.  An empty sequence is returned when
            no more rows are available.  If a DictCursor is used then the
            output consists of a list of dictionaries, one per row.
        """
        if not size:
            size = self.arraysize
        out = []
        for _ in range(size):
            out.append(self.fetchone())
        return out

    def fetchall(self):
        """Fetch all remaining rows from the result set.

        Returns
        -------
        list of tuples :
            A list of tuples, one per row.  An empty sequence is returned when
            no more rows are available.  If a DictCursor is used then the
            output consists of a list of dictionaries, one per row.
        """
        out = []
        for _ in range(self._rowcount - self.pos):
            out.append(self.fetchone())
        return out

    def scroll(self, value, mode="relative"):
        """Move the row cursor.

        Parameters
        ----------
        value : str
            The number of rows to skip or the row number to position to.
        mode : str
            Either "relative" for a relative skip (default), or "absolute"
            to position to a row by its absolute index within the result set
            (zero-indexed).
        """
        if mode == "absolute":
            if value > 0:
                self.pos = value
            else:
                raise DataError("row number not valid:" + str(value))
        elif mode == "relative":
            self.pos += value

    def close(self):
        """Close the cursor object and free all resources.  This implementation
        does nothing.  It is provided for compliance with the Python Database
        API.
        """
        # this can remain implemented as "pass"
        pass
