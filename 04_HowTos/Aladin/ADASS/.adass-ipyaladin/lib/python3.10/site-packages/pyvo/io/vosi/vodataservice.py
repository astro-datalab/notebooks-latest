# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains xml element classes as defined in the VODataService standard

There are different ways of handling the various xml tags.

* Elements with complex content
* Elements with simple content and attributes
* Elements with simple content without attributes

Elements with complex content are parsed with objects inherited from `~pyvo.utils.xml.elements.Element`.

Elements with simple content are parsed with objects inherited from `~pyvo.utils.xml.elements.Element`
defining a ``value`` property.
"""

import re

from astropy.utils import deprecated
from astropy.utils.collections import HomogeneousList
from textwrap import indent
from astropy.utils.xml import check as xml_check
from astropy.io.votable.exceptions import vo_raise, vo_warn, warn_or_raise

from ...utils.xml.elements import (
    xmlattribute, xmlelement, Element, ElementWithXSIType, ContentMixin)

from . import voresource as vr
from .exceptions import (
    W01, W02, W03, W04, W05, W06, W07, W08, W09, W10, W11, W12, W13, W14, W17,
    W18, W36, W37,
    E01, E02, E03, E04, E05, E06)

__all__ = [
    "TableSet", "TableSchema", "ParamHTTP", "VODataServiceTable", "BaseParam", "TableParam",
    "InputParam", "DataType", "SimpleDataType", "TableDataType", "VOTableType",
    "TAPDataType", "TAPType", "FKColumn", "ForeignKey"]


######################################################################
# FACTORY FUNCTIONS
def _convert_boolean(value, default=None):
    return {
        'false': False,
        '0': False,
        'true': True,
        '1': True
    }.get(value, default)


######################################################################
# ATTRIBUTE CHECKERS
def check_anyuri(uri, config=None, pos=None):
    """
    Raises a `~pyvo.io.vosi.tables.exceptions.VOSITablesWarning` if
    *uri* is not a valid URI.
    As defined in RFC 2396.
    """
    if uri is not None and not xml_check.check_anyuri(uri):
        warn_or_raise(W01, W01, uri, config=config, pos=pos)
        return False
    return True


def check_datatype_flag(data, config=None, pos=None):
    """
    Checks if the datatype flag is valid
    """
    if data not in ('indexed', 'primary', 'nullable'):
        warn_or_raise(W04, W04, data, config=config, pos=pos)
        return False
    return True


######################################################################
# ELEMENT CLASSES
class TableSet(Element, HomogeneousList):
    """
    TableSet element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    The set of tables hosted by a resource.
    """

    def __init__(
        self, config=None, pos=None, _name='tableset', version='1.1', **kwargs
    ):
        HomogeneousList.__init__(self, TableSchema)
        Element.__init__(self, config, pos, _name, **kwargs)

        self._version = version

    def __repr__(self):
        return '<TableSet>... {} schemas ...</TableSet>'.format(
            len(self))

    @xmlattribute
    def version(self):
        """The version of the standard"""
        return self._version

    @version.setter
    def version(self, version):
        self._config['version'] = version
        self._version = version

    @xmlelement(name='schema')
    def schemas(self):
        """
        A list of schemas. Must contain only `Schema` objects.

        A named description of a set of logically related tables.

        The name given by the "name" child element must be unique within this
        TableSet instance.  If there is only one schema in this set and/or
        there's no locally appropriate name to provide, the name can be set to
        "default".

        This aggregation does not need to map to an actual database, catalog,
        or schema, though the publisher may choose to aggregate along such
        designations, or particular service protocol may recommend it.
        """
        return self

    @schemas.adder
    def schemas(self, iterator, tag, data, config, pos):
        schema = TableSchema(config, pos, 'schema', **data)
        schema.parse(iterator, config)
        self.append(schema)

    def parse(self, iterator, config):
        super().parse(iterator, config)

        if not self.schemas:
            warn_or_raise(W14, W14, config=config, pos=self._pos)


class TableSchema(Element, HomogeneousList):
    """
    TableSchema element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    A detailed description of a logically-related set of tables.
    """

    def __init__(self, config=None, pos=None, _name='schema', **kwargs):
        HomogeneousList.__init__(self, VODataServiceTable)
        Element.__init__(self, config, pos, _name, **kwargs)

        self._name = None
        self._title = None
        self._description = None
        self._utype = None

    def __repr__(self):
        return '<TableSchema name="{}">... {} tables ...</TableSchema>'.format(
            self.name, len(self.tables))

    @xmlelement(plain=True, multiple_exc=W05)
    def name(self):
        """
        A name for the set of tables.

        This is used to uniquely identify the table set among
        several table sets.  If a title is not present, this
        name can be used for display purposes.

        If there is no appropriate logical name associated with
        this set, the name should be explicitly set to
        "default".
        """
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @xmlelement(plain=True, multiple_exc=W13)
    def title(self):
        """
        a descriptive, human-interpretable name for the table set.

        This is used for display purposes.  There is no requirement
        regarding uniqueness.  It is useful when there are
        multiple schemas in the context (e.g. within a
        tableset; otherwise, the resource title could be
        used instead).
        """
        return self._title

    @title.setter
    def title(self, title):
        self._title = title

    @xmlelement(plain=True, multiple_exc=W06)
    def description(self):
        """
        A free text description of the tableset that should
        explain in general how all of the tables are related.
        """
        return self._description

    @description.setter
    def description(self, description):
        self._description = description

    @xmlelement(plain=True, multiple_exc=W09)
    def utype(self):
        """
        an identifier for a concept in a data model that
        the data in this schema as a whole represent.

        The format defined in the VOTable standard is strongly
        recommended.
        """
        return self._utype

    @utype.setter
    def utype(self, utype):
        self._utype = utype

    @xmlelement(name='table')
    def tables(self):
        """
        A list of tables in the schema.  Must contain only `Table` objects.

        A description of one of the tables that makes up the set.
        The table names for the table should be unique.
        """
        return self

    @tables.adder
    def tables(self, iterator, tag, data, config, pos):
        table = VODataServiceTable(config, pos, 'table', **data)
        table.parse(iterator, config)
        self.append(table)

    def parse(self, iterator, config):
        super().parse(iterator, config)

        if not self.name:
            vo_raise(E06, self._Element__name, config=config, pos=self._pos)


@vr.Interface.register_xsi_type('vs:ParamHTTP')
class ParamHTTP(vr.Interface):
    """
    ParamHTTP element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    A service invoked via an HTTP Query (either Get or Post)
    with a set of arguments consisting of keyword name-value pairs.

    Note that the URL for help with this service can be put into
    the Service/ReferenceURL element.
    """

    def __init__(self, config=None, pos=None, _name='', **kwargs):
        super().__init__(
            config=config, pos=pos, _name=_name, **kwargs)

        self._querytypes = HomogeneousList(str)
        self._resulttype = None

    @xmlelement(name='queryType', multiple_exc=W17)
    def querytypes(self):
        """
        The type of HTTP request, either GET or POST.

        The service may indicate support for both GET
        and POST by providing 2 queryType elements, one
        with GET and one with POST.
        """
        return self._querytypes

    @xmlelement(name='resultType', multiple_exc=W36)
    def resulttype(self):
        """The MIME type of a document returned in the HTTP response."""
        return self._resulttype

    @resulttype.setter
    def resulttype(self, resulttype):
        self._resulttype = resulttype

    def parse(self, iterator, config):
        super().parse(iterator, config)

        if len(self.querytypes) > 2:
            warn_or_raise(W18, W18, config=config, pos=self._pos)


class VODataServiceTable(Element):
    """
    Table element as described in
    http://www.ivoa.net/xml/VODataService/v1.1
    """

    def __init__(
        self, config=None, pos=None, _name='table', version='1.1', **kwargs
    ):
        super().__init__(config, pos, _name, **kwargs)

        self._name = None
        self._title = None
        self._description = None
        self._utype = None
        self._type = kwargs.get("type")
        self._version = version
        self._nrows = None

        self._columns = HomogeneousList(TableParam)
        self._foreignkeys = HomogeneousList(ForeignKey)

    def __repr__(self):
        return '<VODataServiceTable name="{}">... {} columns ...</VODataServiceTable>'.format(
            self.name, len(self.columns))

    def describe(self):
        print(self.name)
        if self.description is not None:
            print(indent(self.description, 4 * " "))
        else:
            print('No description')

        print()

    @xmlelement(plain=True, multiple_exc=W05)
    def name(self):
        """
        the fully qualified name of the table.  This name
        should include all catalog or schema prefixes
        needed to sufficiently uniquely distinguish it in a
        query.

        In general, the format of the qualified name may
        depend on the context; however, when the
        table is intended to be queryable via ADQL, then the
        catalog and schema qualifiers are delimited from the
        table name with dots (.).
        """
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @xmlelement(plain=True, multiple_exc=W13)
    def title(self):
        """
        a descriptive, human-interpretable name for the table.

        This is used for display purposes.  There is no requirement
        regarding uniqueness.
        """
        return self._title

    @title.setter
    def title(self, title):
        self._title = title

    @xmlelement(plain=True, multiple_exc=W06)
    def description(self):
        """
        a free-text description of the table's contents
        """
        return self._description

    @description.setter
    def description(self, description):
        self._description = description

    @xmlelement(plain=True, multiple_exc=W09)
    def utype(self):
        """
        an identifier for a concept in a data model that
        the data in this table represent.

        The format defined in the VOTable standard is highly
        recommended.
        """
        return self._utype

    @utype.setter
    def utype(self, utype):
        self._utype = utype

    @xmlelement(plain=True, multiple_exc=W09)
    def nrows(self):
        """
        the approximate number of rows in the table.

        This is None if the data provider failed to provide this
        information.
        """
        return self._nrows

    @nrows.setter
    def nrows(self, nrows):
        self._nrows = int(nrows)

    @xmlattribute
    def type(self):
        """
        a name for the role this table plays.  Recognized
        values include "output", indicating this table is output
        from a query; "base_table", indicating a table
        whose records represent the main subjects of its
        schema; and "view", indicating that the table represents
        a useful combination or subset of other tables.  Other
        values are allowed.
        """
        return self._type

    @type.setter
    def type(self, type_):
        self._type = type_

    @xmlattribute
    def version(self):
        """The version of the standard"""
        return self._version

    @version.setter
    def version(self, version):
        self._config['version'] = version
        self._version = version

    @xmlelement(name='column')
    def columns(self):
        """
        A list of columns in the table.
        Must contain only `TableParams` objects.

        A description of a table column.
        """
        return self._columns

    @columns.adder
    def columns(self, iterator, tag, data, config, pos):
        column = TableParam(config, pos, 'column', **data)
        column.parse(iterator, config)
        self.columns.append(column)

    @xmlelement(name='foreignKey')
    def foreignkeys(self):
        """
        A list of columns in the table.  Must contain only `ForeignKey` objects

        a description of a foreign keys, one or more columns
        from the current table that can be used to join with
        another table.
        """
        return self._foreignkeys

    @foreignkeys.adder
    def foreignkeys(self, iterator, tag, data, config, pos):
        foreignkey = ForeignKey(config, pos, 'foreignKey', **data)
        foreignkey.parse(iterator, config)
        self.foreignkeys.append(foreignkey)

    def parse(self, iterator, config):
        super().parse(iterator, config)

        if not self.name:
            vo_raise(E06, self._Element__name, config=config, pos=self._pos)


@deprecated("1.5", alternative="VODataServiceTable")
class Table(VODataServiceTable):
    pass


class BaseParam(Element):
    """
    BaseParam element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    a description of a parameter that places no restriction on the parameter's
    data type. As the parameter's data type is usually important, schemas
    normally employ a sub-class of this type (e.g. Param), rather than this
    type directly.
    """

    def __init__(self, config=None, pos=None, _name='', **kwargs):
        super().__init__(
            config=config, pos=pos, _name=_name, **kwargs)

        self._name = None
        self._description = None
        self._unit = None
        self._ucd = None
        self._utype = None

    def __repr__(self):
        return '<BaseParam name="{}"/>'.format(self.name)

    @xmlelement(plain=True, multiple_exc=W05)
    def name(self):
        """the name of the element"""
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @xmlelement(plain=True, multiple_exc=W06)
    def description(self):
        """
        a free-text description of the element's contents
        """
        return self._description

    @description.setter
    def description(self, description):
        self._description = description

    @xmlelement(plain=True, multiple_exc=W07)
    def unit(self):
        """the unit associated with all values in the element"""
        return self._unit

    @unit.setter
    def unit(self, unit):
        self._unit = unit

    @xmlelement(plain=True, multiple_exc=W08)
    def ucd(self):
        """
        the name of a unified content descriptor that
        describes the scientific content of the element.

        There are no requirements for compliance with any
        particular UCD standard.  The format of the UCD can
        be used to distinguish between UCD1, UCD1+, and
        SIA-UCD.  See
        http://www.ivoa.net/Documents/latest/UCDlist.html
        for the latest IVOA standard set.
        """
        return self._ucd

    @ucd.setter
    def ucd(self, ucd):
        self._ucd = ucd

    @xmlelement(plain=True, multiple_exc=W09)
    def utype(self):
        """
        an identifier for a concept in a data model that
        the data in this element represent.

        The format defined in the VOTable standard is highly recommended.
        """
        return self._utype

    @utype.setter
    def utype(self, utype):
        self._utype = utype


class TableParam(BaseParam):
    """
    TableParam element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    A description of a table parameter having a fixed data type.
    The allowed data type names match those supported by VOTable.
    """
    @classmethod
    def from_field(cls, field):
        """
        Create a instance from a `~astropy.io.votable.tree.Field` instance.
        """
        instance = cls()

        instance.name = field.name
        instance.description = field.description
        instance.unit = field.unit
        instance.ucd = field.ucd
        instance.utype = field.utype

        datatype = VOTableType(arraysize=field.arraysize)
        datatype.value = field.datatype

        instance.datatype = datatype

        return instance

    def __init__(self, config=None, pos=None, _name='', std=None, **kwargs):
        super().__init__(
            config=config, pos=pos, _name=_name, **kwargs)

        self._datatype = None
        self._flags = HomogeneousList(str)
        self._std = _convert_boolean(std)

    @xmlelement(name='dataType')
    def datatype(self):
        """The type of data contained in the element"""
        return self._datatype

    @datatype.setter
    def datatype(self, datatype):
        if datatype is not None and not isinstance(datatype, TableDataType):
            raise ValueError("datatype must be an TableDataType object")
        self._datatype = datatype

    @datatype.adder
    def datatype(self, iterator, tag, data, config, pos):
        datatype = TableDataType(config, pos, 'dataType', **data)
        datatype.parse(iterator, config)

        if self.datatype:
            warn_or_raise(
                W37, args=self._Element__name, config=config, pos=pos)
        self.datatype = datatype

    @xmlelement(name='flag')
    def flags(self):
        """
        A list of flags. Must contain only `str` objects.

        a keyword representing traits of the column. Recognized values include
        "indexed", "primary", and "nullable".
        """
        return self._flags

    @xmlattribute
    def std(self):
        """
        If true, the meaning and use of this parameter is
        reserved and defined by a standard model.  If false,
        it represents a database-specific parameter
        that effectively extends beyond the standard.  If
        not provided, then the value is unknown.
        """
        return self._std

    @std.setter
    def std(self, std):
        self._std = std

    def parse(self, iterator, config):
        super().parse(iterator, config)

        if not self.name:
            vo_raise(E06, self._Element__name, config=config, pos=self._pos)


class InputParam(BaseParam):
    """
    InputParam element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    A description of a service or function parameter having a fixed data type.
    """

    def __init__(
            self, config=None, pos=None, _name='', use="optional", std="1",
            **kwargs):
        BaseParam.__init__(self, config, pos, _name, **kwargs)

        self._datatype = None
        self._use = use
        self._std = _convert_boolean(std, True)

    @xmlelement(name='dataType')
    def datatype(self):
        """The type of data contained in the element"""
        return self._datatype

    @datatype.setter
    def datatype(self, datatype):
        if datatype is not None and not isinstance(datatype, SimpleDataType):
            raise ValueError("datatype must be an SimpleDataType object")
        self._datatype = datatype

    @xmlattribute
    def use(self):
        """
        An indication of whether this parameter is required to be provided for
        the application or service to work properly.
        Allowed values are "required" and "optional".
        """
        return self._use

    @use.setter
    def use(self, use):
        self._use = use

    @xmlattribute
    def std(self):
        """
        If true, the meaning and behavior of this parameter is
        reserved and defined by a standard interface.  If
        false, it represents an implementation-specific
        parameter that effectively extends the behavior of the
        service or application.
        """
        return self._std

    @std.setter
    def std(self, std):
        self._std = std

    def parse(self, iterator, config):
        super().parse(iterator, config)

        if not self.name:
            vo_raise(E06, self._Element__name, config=config, pos=self._pos)


class DataType(ContentMixin, ElementWithXSIType):
    """
    DataType element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    A type (in the computer language sense) associated with a parameter with an
    arbitrary name.

    This XML type is used as a parent for defining data types with a restricted
    set of names.
    """

    def __init__(
            self, config=None, pos=None, _name='dataType',
            arraysize=None, delim=None, extendedType=None, extendedSchema=None,
            **kwargs
    ):
        super().__init__(
            config=config, pos=pos, _name=_name, **kwargs)

        if arraysize is None:
            arraysize = "1"

        if delim is None:
            delim = " "

        self.arraysize = arraysize
        self._delim = delim
        self._extendedtype = extendedType
        self.extendedschema = extendedSchema

    def __repr__(self):
        return '<DataType arraysize={}>{}</DataType>'.format(
            self.arraysize, self.content)

    @xmlattribute
    def arraysize(self):
        """Specifies the size of the dataType"""
        return self._arraysize

    @arraysize.setter
    def arraysize(self, arraysize):
        if all((
                arraysize is not None,
                not re.match(r"^([0-9]+x)*[0-9]*[*]?(s\W)?$", arraysize)
        )):
            vo_raise(E01, arraysize, self._config, self._pos)
        self._arraysize = arraysize

    @xmlattribute
    def delim(self):
        """
        the string that is used to delimit elements of an array
        value when arraysize is not "1".

        Unless specifically disallowed by the context,
        applications should allow optional spaces to
        appear in an actual data value before and after
        the delimiter (e.g. "1, 5" when delim=",").

        the default is " "; i.e. the values are delimited by spaces.
        """
        return self._delim

    @delim.setter
    def delim(self, delim):
        self._delim = delim

    @xmlattribute(name='extendedType')
    def extendedtype(self):
        """
        The data value represented by this type can be
        interpreted as of a custom type identified by
        the value of this attribute.

        If an application does not recognize this
        extendedType, it should attempt to handle value
        assuming the type given by the element's value.
        string is a recommended default type.

        This element may make use of the extendedSchema
        attribute and/or any arbitrary (qualified)
        attribute to refine the identification of the type.
        """

    @extendedtype.setter
    def extendedtype(self, extendedtype):
        self._extendedtype = extendedtype

    @xmlattribute(name='extendedSchema')
    def extendedschema(self):
        """
        An identifier for the schema that the value given
        by the extended attribute is drawn from.

        This attribute is normally ignored if the
        extendedType attribute is not present.
        """
        return self._extendedschema

    @extendedschema.setter
    def extendedschema(self, extendedschema):
        if extendedschema is not None:
            check_anyuri(extendedschema, self._config, self._pos)
        self._extendedschema = extendedschema


class SimpleDataType(DataType):
    """
    SimpleDataType element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    A data type restricted to a small set of names which is imprecise as to the
    format of the individual values.

    This set is intended for describing simple input parameters to a service or
    function.
    """

    def _content_check(self, value):
        if value is not None:
            valid_values = {
                'integer', 'real', 'complex', 'boolean', 'char', 'string'}
            if value not in valid_values:
                vo_warn(W02, value, self._config, self._pos)


class TableDataType(DataType):
    """
    TableDataType element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    an abstract parent for a class of data types that can be
    used to specify the data type of a table column.

    Subtypes must be decorated with ``register_xsi_type('ns:name')``.
    """


@TableDataType.register_xsi_type('vs:VOTable')
@TableDataType.register_xsi_type('vs:VOTableType')
class VOTableType(TableDataType):
    """
    VOTableType element as described in
    http://www.ivoa.net/xml/VODataService/v1.1
    """

    def _content_check(self, value):
        if value is not None:
            valid_values = (
                'boolean', 'bit', 'unsignedByte', 'short', 'int', 'long',
                'char', 'unicodeChar', 'float', 'double',
                'floatComplex', 'doubleComplex')
            if value not in valid_values:
                vo_warn(W02, value, self._config, self._pos)


class TAPDataType(TableDataType):
    """
    TAPDataType element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    an abstract parent for the specific data types supported by the
    Table Access Protocol.
    """

    def __init__(
        self, config=None, pos=None, _name='dataType', size=None, **kwargs
    ):
        super().__init__(
            config=config, pos=pos, _name=_name, **kwargs)

        self.size = size

    @xmlattribute
    def size(self):
        """
        the length of the fixed-length value.

        This corresponds to the size Column attribute in
        the TAP_SCHEMA and can be used with data types
        that are defined with a length (CHAR, BINARY).
        """
        return self._size

    @size.setter
    def size(self, size):
        if size is not None and int(size) < 0:
            size = 0
            warn_or_raise(W03, W03, config=self._config, pos=self._pos)
        self._size = size


@TableDataType.register_xsi_type('vs:TAP')
@TableDataType.register_xsi_type('vs:TAPType')
class TAPType(TAPDataType):
    """
    TAPType element as described in
    http://www.ivoa.net/xml/VODataService/v1.1

    a data type supported explicitly by the Table Access Protocol (v1.0).
    """

    def _content_check(self, value):
        if value is not None:
            valid_values = (
                'BOOLEAN', 'SMALLINT', 'INTEGER', 'BIGINT', 'REAL', 'DOUBLE',
                'TIMESTAMP', 'CHAR', 'VARCHAR', 'BINARY', 'VARBINARY',
                'POINT', 'REGION', 'CLOB', 'BLOB')
            if value not in valid_values:
                vo_warn(W02, value, self._config, self._pos)


class FKColumn(Element):
    """
    FKColumn element as described in
    http://www.ivoa.net/xml/VODataService/v1.1
    """

    def __init__(self, config=None, pos=None, _name='fkColumn', **kwargs):
        super().__init__(
            config=config, pos=pos, _name=_name, **kwargs)

        self._fromcolumn = None
        self._targetcolumn = None

    def __repr__(self):
        return '<FKColumn fromColumn={} targetColumn={}>...</FKColumn>'.format(
            self.fromcolumn, self.targetcolumn)

    @xmlelement(name='fromColumn', plain=True, multiple_exc=W10)
    def fromcolumn(self):
        """
        The unqualified name of the column from the current table.
        """
        return self._fromcolumn

    @fromcolumn.setter
    def fromcolumn(self, fromcolumn):
        self._fromcolumn = fromcolumn

    @xmlelement(name='targetColumn', plain=True, multiple_exc=W11)
    def targetcolumn(self):
        """
        The unqualified name of the column from the target table.
        """
        return self._targetcolumn

    @targetcolumn.setter
    def targetcolumn(self, targetcolumn):
        self._targetcolumn = targetcolumn

    def parse(self, iterator, config):
        super().parse(iterator, config)

        if self.fromcolumn is None:
            vo_raise(E02, config=config, pos=self._pos)
        if self.targetcolumn is None:
            vo_raise(E03, config=config, pos=self._pos)


class ForeignKey(Element):
    """
    ForeignKey element as described in
    http://www.ivoa.net/xml/VODataService/v1.1
    """

    def __init__(self, config=None, pos=None, _name='foreignKey', **kwargs):
        Element.__init__(self, config, pos, _name, **kwargs)

        self._targettable = None
        self._fkcolumns = HomogeneousList(FKColumn)
        self._description = None
        self._utype = None

    def __repr__(self):
        return '<ForeignKey targetTable={}>...</ForeignKey>'.format(
            self.targettable)

    @xmlelement(name='targetTable', plain=True, multiple_exc=W12)
    def targettable(self):
        """
        the fully-qualified name (including catalog and schema, as
        applicable) of the table that can be joined with the
        table containing this foreign key.
        """
        return self._targettable

    @targettable.setter
    def targettable(self, targettable):
        self._targettable = targettable

    @xmlelement(name='fkColumn')
    def fkcolumns(self):
        """
        A list of foreign key columns. Must contain only `FKColumn` objects.

        a pair of column names, one from this table and one
        from the target table that should be used to join the
        tables in a query.
        """
        return self._fkcolumns

    @fkcolumns.adder
    def fkcolumns(self, iterator, tag, data, config, pos):
        fkcolumn = FKColumn(config, pos, 'fkColumn', **data)
        fkcolumn.parse(iterator, config)
        self.fkcolumns.append(fkcolumn)

    @xmlelement(plain=True, multiple_exc=W06)
    def description(self):
        """
        a free-text description of what this key points to
        and what the relationship means.
        """
        return self._description

    @description.setter
    def description(self, description):
        self._description = description

    @xmlelement(plain=True, multiple_exc=W09)
    def utype(self):
        """
        an identifier for a concept in a data model that
        the association enabled by this key represents.

        The format defined in the VOTable standard is highly
        recommended.
        """
        return self._utype

    @utype.setter
    def utype(self, utype):
        self._utype = utype

    def parse(self, iterator, config):
        super().parse(iterator, config)

        if not self.targettable:
            vo_raise(E04, config=config, pos=self._pos)
        if not self.fkcolumns:
            vo_raise(E05, config=config, pos=self._pos)
