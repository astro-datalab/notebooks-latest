# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
.. _warnings:
Warnings
--------
.. note::
    Most of the following warnings indicate violations of the VOSI
    specification.  They should be reported to the authors of the
    tools that produced the VOSI file.
    To control the warnings emitted, use the standard Python
    :mod:`warnings` module.  Most of these are of the type
    `VOSISpecWarning`.
{warnings}
.. _exceptions:
Exceptions
----------
.. note::
    This is a list of many of the fatal exceptions emitted by vosi.endpoint
    when the file does not conform to spec.  Other exceptions may be
    raised due to unforeseen cases or bugs in vosi.endpoint itself.
{exceptions}
"""

from astropy.utils.exceptions import AstropyWarning
from ...utils.xml.exceptions import XMLWarning

__all__ = ["VOSIWarning"]
__all__ += ["W{:0>2}".format(i) for i in range(1, 36)]
__all__ += ["E{:0>2}".format(i) for i in range(1, 10)]


class VOSIWarning(AstropyWarning):
    """
    The base class of all VOSI warnings and exceptions.
    Handles the formatting of the message with a warning or exception
    code, filename, line and column number.
    """


class W01(VOSIWarning, XMLWarning):
    """
    The attribute must be a valid URI as defined in `RFC 2396
    <http://www.ietf.org/rfc/rfc2396.txt>`_.
    """

    message_template = "'{}' is not a valid URI"
    default_args = ('x',)


class W02(VOSIWarning, XMLWarning):
    """
    The attribute must be any of the accepted types in the VOSI spec.
    """
    message_template = (
        "'{}' is not a valid datatype according to the VOSI spec")
    default_args = ('x',)


class W03(VOSIWarning, XMLWarning):
    """
    The attribute must be an positive integer.
    """
    message_template = "Size must be positive"


class W04(VOSIWarning, XMLWarning):
    """
    The attribute must have one of the recognized values
    'indexed', 'primary', 'nullable'.
    """
    message_template = "'{}' is not a recognized flag"
    default_args = ('x',)


class W05(VOSIWarning, XMLWarning):
    """
    A ``name`` element can only appear once within its parent element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = "{} element contains more than one name element"
    default_args = ('x',)


class W06(VOSIWarning, XMLWarning):
    """
    A ``description`` element can only appear once within its parent element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = "{} element contains more than one description element"
    default_args = ('x',)


class W07(VOSIWarning, XMLWarning):
    """
    A ``unit`` element can only appear once within its parent element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = "{} element contains more than one unit element"
    default_args = ('x',)


class W08(VOSIWarning, XMLWarning):
    """
    A ``ucd`` element can only appear once within its parent element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = "{} element contains more than one ucd element"
    default_args = ('x',)


class W09(VOSIWarning, XMLWarning):
    """
    A ``utype`` element can only appear once within its parent element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = "{} element contains more than one utype element"
    default_args = ('x',)


class W10(VOSIWarning, XMLWarning):
    """
    A ``fromColumn`` element can only appear once within its parent element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = "{} element contains more than one fromColumn element"
    default_args = ('x',)


class W11(VOSIWarning, XMLWarning):
    """
    A ``targetColumn`` element can only appear once within its parent element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = "{} element contains more than one targetColumn element"
    default_args = ('x',)


class W12(VOSIWarning, XMLWarning):
    """
    A ``targetTable`` element can only appear once within its parent element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = "{} element contains more than one targetTable element"
    default_args = ('x',)


class W13(VOSIWarning, XMLWarning):
    """
    A ``title`` element can only appear once within its parent element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = "{} element contains more than one title element"
    default_args = ('x',)


class W14(VOSIWarning, XMLWarning):
    """
    The tableset element must contain at least one schema element.
    """

    message_template = (
        "tableset element must contain at least one schema element.")


class W15(VOSIWarning, XMLWarning):
    """
    Unknown issues may arise using ``dal`` with VOSITables files
    from a version other than 1.0 or 1.1
    """

    message_template = (
        'pyvo.dal is designed for VOSITables version 1.0, and 1.1, but '
        + 'this file is {}')
    default_args = ('x',)


class W16(VOSIWarning, XMLWarning):
    """
    The table element is not a valid root element in VOSI before version 1.1
    """
    message_template = (
        "The element table is not a valid root element in VOSI below v1.1")


class W17(VOSIWarning, XMLWarning):
    """
    A ``queryType`` element can only appear once within the ParamHTTP element.
    According to the schema, it may only occur once (`1.1
    <http://www.ivoa.net/documents/VODataService/20101202/REC-VODataService-1.1-20101202.html#appA>`__,
    """

    message_template = (
        "ParamHTTP element contains more than one ParamHTTP element")


class W18(VOSIWarning, XMLWarning):
    """
    The QueryType element must not occur more than two times.
    """
    message_template = (
        "The QueryType element must not occur more than two times.")


class W19(VOSIWarning, XMLWarning):
    """
    TAP Capabilities must not have an ivo-id other than ivo://ivoa.net/std/TAP
    """
    message_template = (
        "TAP Capabilities must not have an ivo-id other than "
        "ivo://ivoa.net/std/TAP"
    )


class W20(VOSIWarning, XMLWarning):
    """
    TAP Capabilties must have at least one `language` element.
    """
    message_template = (
        "TAP Capabilties must have at least one `language` element.")


class W21(VOSIWarning, XMLWarning):
    """
    TAP Capabilties must have at least one outputFormat element.
    """
    message_template = (
        "TAP Capabilties must have at least one `outputFormat` element.")


class W22(VOSIWarning, XMLWarning):
    """
    The `retentionPeriod` element must not occur more than once.
    """
    message_template = (
        "The retentionPeriod element must not occur more than once")


class W23(VOSIWarning, XMLWarning):
    """
    The `executionDuration` element must not occur more than once.
    """
    message_template = (
        "The executionDuration element must not occur more than once")


class W24(VOSIWarning, XMLWarning):
    """
    The `outputLimit` element must not occur more than once.
    """
    message_template = (
        "The outputLimit element must not occur more than once")


class W25(VOSIWarning, XMLWarning):
    """
    The `uploadLimit` element must not occur more than once.
    """
    message_template = (
        "The uploadLimit element must not occur more than once")


class W26(VOSIWarning, XMLWarning):
    """
    The ivo-id attribute is mandatory.
    """
    message_template = "The ivo-id attribute is mandatory"


class W27(VOSIWarning, XMLWarning):
    """
    The `form` element must not occur more than once.
    """
    message_template = "The form element must not occur more than once"


class W28(VOSIWarning, XMLWarning):
    """
    The `mime` element must not occur more than once.
    """
    message_template = "The mime element must not occur more than once"


class W29(VOSIWarning, XMLWarning):
    """
    The `default` element must not occur more than once.
    """
    message_template = "The default element must not occur more than once"


class W30(VOSIWarning, XMLWarning):
    """
    The `hard` element must not occur more than once.
    """
    message_template = "The hard element must not occur more than once"


class W31(VOSIWarning, XMLWarning):
    """
    The content of the `DataLimit` element must be byte or row
    """
    message_template = (
        "The content of the DataLimit element must be byte or row")


class W32(VOSIWarning, XMLWarning):
    """
    The `available` element must not occur more than once.
    """
    message_template = "The available element must not occur more than once"


class W33(VOSIWarning, XMLWarning):
    """
    The `upSince` element must not occur more than once.
    """
    message_template = "The upSince element must not occur more than once"


class W34(VOSIWarning, XMLWarning):
    """
    The `downAt` element must not occur more than once.
    """
    message_template = "The downAt element must not occur more than once"


class W35(VOSIWarning, XMLWarning):
    """
    The `backAt` element must not occur more than once.
    """
    message_template = "The backAt element must not occur more than once"


class W36(VOSIWarning, XMLWarning):
    """
    The `resultType` element must not occur more than once.
    """
    message_template = "The resultType element must not occur more than once"


class W37(VOSIWarning, XMLWarning):
    """
    The `dataType` element must not occur more than once.
    """
    message_template = "The dataType element must not occur more than once"


class E01(VOSIWarning, XMLWarning, ValueError):
    r"""
    The attribute must be a valid arraysize according to the VOTable standard.
    From the VOTable 1.2 spec:
        A table cell can contain an array of a given primitive type,
        with a fixed or variable number of elements; the array may
        even be multidimensional. For instance, the position of a
        point in a 3D space can be defined by the following::
            <FIELD ID="point_3D" datatype="double" arraysize="3"/>
        and each cell corresponding to that definition must contain
        exactly 3 numbers. An asterisk (\*) may be appended to
        indicate a variable number of elements in the array, as in::
            <FIELD ID="values" datatype="int" arraysize="100*"/>
        where it is specified that each cell corresponding to that
        definition contains 0 to 100 integer numbers. The number may
        be omitted to specify an unbounded array (in practice up to
        =~2×10⁹ elements).
        A table cell can also contain a multidimensional array of a
        given primitive type. This is specified by a sequence of
        dimensions separated by the ``x`` character, with the first
        dimension changing fastest; as in the case of a simple array,
        the last dimension may be variable in length. As an example,
        the following definition declares a table cell which may
        contain a set of up to 10 images, each of 64×64 bytes::
            <FIELD ID="thumbs" datatype="unsignedByte" arraysize="64×64×10*"/>
    **References**: `1.1
    <http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#sec:dim>`__,
    `1.2
    <http://www.ivoa.net/Documents/VOTable/20091130/REC-VOTable-1.2.html#sec:dim>`__
    """

    message_template = "Invalid arraysize attribute '{}'"
    default_args = ('x',)


class E02(VOSIWarning, XMLWarning, ValueError):
    """
    The `FKColumn` element must have a `fromColumn`.
    """
    message_template = "fkColumn element is missing a fromColumn"


class E03(VOSIWarning, XMLWarning, ValueError):
    """
    The element must have a `targetColumn`.
    """
    message_template = "The element is missing a targetColumn"


class E04(VOSIWarning, XMLWarning, ValueError):
    """
    The element must have a `targetTable`.
    """
    message_template = "The element is missing a targetTable"


class E05(VOSIWarning, XMLWarning, ValueError):
    """
    The element must contain at least one `fkColumn`.
    """
    message_template = "The element contains no `fkColumn`"


class E06(VOSIWarning, XMLWarning, ValueError):
    """
    The element must have a ``name`` element.
    """
    message_template = "The {} element must have a name element"
    default_args = ('x',)


class E07(VOSIWarning, XMLWarning, ValueError):
    """
    Raised either when the file doesn't appear to be XML, or the root
    element is not tableset or table.
    """
    message_template = "File does not appear to be a VOSITables file"


class E08(VOSIWarning, XMLWarning, ValueError):
    """
    The element must have a ``version`` element.
    """
    message_template = "The {} element must have a version element"
    default_args = ('x',)


class E09(VOSIWarning, XMLWarning, ValueError):
    """
    The element must have a ``form`` element.
    """
    message_template = "The {} element must have a form element"
    default_args = ('x',)


class E10(VOSIWarning, XMLWarning, ValueError):
    """
    Raised when then file doesn't appear to be valid capabilities xml
    """
    message_template = "File does not appear to be a VOSICapabilities file"


class VOSIError(Exception):
    """
    Raised for non-XML VOSI errors
    """
    pass
