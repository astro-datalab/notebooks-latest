# Licensed under a 3-clause BSD style license - see LICENSE.rst

from inspect import getmembers
from functools import partial
import warnings

from astropy.utils.xml import iterparser
from astropy.io.votable.exceptions import warn_or_raise
from pyvo.utils.xml.exceptions import UnknownElementWarning

__all__ = [
    "xmlattribute", "xmlelement",
    "make_add_complexcontent", "make_add_simplecontent",
    "Element", "ElementWithXSIType", "ContentMixin", "parse_for_object"]


def parse_for_object(
    source, object_type, pedantic=None, filename=None,
        _debug_python_based_parser=False
):
    """
    Parses an xml file (or file-like object), and returns a
    object of specified object_type. object_type must be a subtype of
    `~pyvo.utils.xml.elements.Element` type

    Parameters
    ----------
    source : str or readable file-like object
        Path or file object containing a tableset xml file.
    object : object type to return (subtype `~pyvo.utils.xml.elements.Element`)
    pedantic : bool, optional
        When `True`, raise an error when the file violates the spec,
        otherwise issue a warning.  Warnings may be controlled using
        the standard Python mechanisms.  See the `warnings`
        module in the Python standard library for more information.
        Defaults to False.
    filename : str, optional
        A filename, URL or other identifier to use in error messages.
        If *filename* is None and *source* is a string (i.e. a path),
        then *source* will be used as a filename for error messages.
        Therefore, *filename* is only required when source is a
        file-like object.

    Returns
    -------
    object : `~pyvo.utils.xml.elements.Element` object or None

    See also
    --------
    pyvo.io.vosi.exceptions : The exceptions this function may raise.
    """
    config = {
        'pedantic': pedantic,
        'filename': filename
    }

    if filename is None and isinstance(source, str):
        config['filename'] = source

    with iterparser.get_xml_iterator(
            source,
            _debug_python_based_parser=_debug_python_based_parser
    ) as iterator:
        return object_type(
            config=config, pos=(1, 1)).parse(iterator, config)


class xmlattribute(property):
    def __init__(self, fget=None, fset=None, fdel=None, doc=None, name=None):
        super().__init__(fget, fset, fdel, doc)
        if name:
            self.name = name
        elif fget is not None:
            self.name = fget.__name__
        else:
            raise ValueError(
                "xmlattribute either needs a getter or a element name or both")

    def __call__(self, fget):
        return self.__class__(fget, name=self.name)

    def getter(self, fget):
        return self.__class__(
            fget, self.fset, self.fdel, self.__doc__, self.name)

    def setter(self, fset):
        return self.__class__(
            self.fget, fset, self.fdel, self.__doc__, self.name)

    def deleter(self, fdel):
        return self.__class__(
            self.fget, self.fset, fdel, self.__doc__, self.name)


class xmlelement(property):
    """

    """

    def __init__(
        self, fget=None, fset=None, fdel=None, fadd=None, fformat=None,
        doc=None, name=None, ns=None, plain=False, cls=None, multiple_exc=None
    ):
        super().__init__(fget, fset, fdel, doc)

        if name:
            self.name = name
        elif fget is not None:
            self.name = fget.__name__
        else:
            self.name = None

        self.ns = ns

        self.plain = plain
        self.cls = cls
        self.multiple_exc = multiple_exc

        self.fadd = fadd
        self.fformat = fformat

    def __call__(self, fget):
        return self.__class__(
            fget, name=self.name or fget.__name__, ns=self.ns,
            plain=self.plain, cls=self.cls, multiple_exc=self.multiple_exc
        )

    def __get__(self, obj, owner=None):
        if obj is not None:
            val = super().__get__(obj, owner)
            if self.plain:
                return val
            elif not isinstance(val, (Element, list)):
                element = ContentMixin(_name=self.name, _ns=self.ns)
                element.content = val
                return element
            else:
                return val
        else:
            return super().__get__(obj, owner)

    def getter(self, fget):
        return self.__class__(
            fget, self.fset, self.fdel, self.fadd, self.fformat, self.__doc__,
            self.name, self.ns, self.plain, self.cls, self.multiple_exc)

    def setter(self, fset):
        return self.__class__(
            self.fget, fset, self.fdel, self.fadd, self.fformat, self.__doc__,
            self.name, self.ns, self.plain, self.cls, self.multiple_exc)

    def deleter(self, fdel):
        return type(self)(
            self.fget, self.fset, fdel, self.fadd, self.fformat, self.__doc__,
            self.name, self.ns, self.plain, self.cls, self.multiple_exc)

    def adder(self, fadd):
        if self.cls:
            raise RuntimeError(
                'xmlelement cls parameter has no effect when adder is'
                ' defined')

        if self.multiple_exc:
            raise RuntimeError(
                'xmlelement multiple_exc parameter has no effect when'
                ' adder is defined')

        return self.__class__(
            self.fget, self.fset, self.fdel, fadd, self.fformat, self.__doc__,
            self.name, self.ns, self.plain, self.cls, self.multiple_exc)

    def formatter(self, fformat):
        return self.__class__(
            self.fget, self.fset, self.fdel, self.fadd, fformat, self.__doc__,
            self.name, self.ns, self.plain, self.cls, self.multiple_exc)


def object_attrs(obj):
    objtype = type(obj)
    attrs = {
        getattr(objtype, name).name: value for name, value in getmembers(obj)
        if isinstance(getattr(objtype, name, None), xmlattribute)}
    return attrs


def object_children(obj):
    objtype = type(obj)

    try:
        for child in obj:
            if isinstance(child, Element):
                yield (child._Element__name, None, child)
    except TypeError:
        for name, child in getmembers(obj):
            if child is None:
                continue
            descr = getattr(objtype, name, None)

            if isinstance(descr, xmlelement):
                element_name = descr.name

                if descr.fformat:
                    fformat = partial(descr.fformat, obj)
                else:
                    fformat = None
                yield (element_name, fformat, child)
            elif isinstance(child, Element):
                yield (child._Element__name, None, child)


def object_mapping(obj):
    objtype = type(obj)

    for name, val in getmembers(obj):
        descr = getattr(objtype, name, None)

        if isinstance(descr, xmlelement):
            if descr.fadd is None:
                if descr.cls is None:
                    fadd = make_add_simplecontent(
                        obj, descr.name, name, descr.multiple_exc)
                else:
                    fadd = make_add_complexcontent(
                        obj, descr.name, name, descr.cls, descr.multiple_exc)
            else:
                fadd = partial(descr.fadd, obj)

            yield descr.name, fadd


def make_add_complexcontent(
        self, element_name, attr_name, cls_, exc_class=None):
    """
    Factory for generating add functions for elements with complex content.
    """
    def add_complexcontent(iterator, tag, data, config, pos):
        attr = getattr(self, attr_name)

        element = cls_(
            config=config, pos=pos, _name=element_name, **data)

        if attr and exc_class is not None:
            warn_or_raise(
                exc_class, args=element_name,
                config=config, pos=pos)

        if isinstance(getattr(self, attr_name, None), list):
            getattr(self, attr_name).append(element)
        else:
            setattr(self, attr_name, element)

        element.parse(iterator, config)

    return add_complexcontent


def make_add_simplecontent(
        self, element_name, attr_name, exc_class=None, check_func=None,
        data_func=None):
    """
    Factory for generating add functions for elements with simple content.
    This means elements with no child elements.
    If exc_class is given, warn or raise if element was already set.
    """
    def add_simplecontent(iterator, tag_ignored, data_ignored, config, pos_ignored):
        # Ignored parameters are kept in the API signature to be compatible
        # with other functions.
        for start, tag, data, pos in iterator:
            if not start and tag == element_name:
                attr = getattr(self, attr_name)

                if attr and exc_class:
                    warn_or_raise(
                        exc_class, args=self._Element__name,
                        config=config, pos=pos)
                if check_func:
                    check_func(data, config, pos)
                if data_func:
                    data = data_func(data)

                if isinstance(getattr(self, attr_name), list):
                    getattr(self, attr_name).append(data)
                else:
                    setattr(self, attr_name, data or None)
                break

    return add_simplecontent


class Element:
    """
    A base class for all classes that represent XML elements.

    Subclasses and Mixins must initialize their independent attributes after
    calling ``super().__init__``.
    """

    def __init__(self, config=None, pos=None, _name='', _ns='', **kwargs):
        if config is None:
            config = {}
        self._config = config
        self._pos = pos

        self.__name = _name
        self.__ns = _ns

        self._tag_mapping = {}

    def _add_unknown_tag(self, iterator, tag, data, config, pos):
        if tag != 'xml':
            warn_or_raise(
                UnknownElementWarning, UnknownElementWarning, tag, config, pos)

    def _end_tag(self, tag, data, pos):
        pass

    def _ignore_add(self, iterator, tag, data, config, pos):
        pass

    def parse(self, iterator, config):
        """
        For internal use. Parse the XML content of the children of the element.
        Override this method and do after-parse checks after calling
        ``super().parse``, if you need to.

        Parameters
        ----------
        iterator : xml iterator
            An iterator over XML elements as returned by
            `~astropy.utils.xml.iterparser.get_xml_iterator`.
        config : dict
            The configuration dictionary that affects how certain
            elements are read.
        """
        tag_mapping = dict(object_mapping(self))

        for start, tag, data, pos in iterator:
            if start:
                tag_mapping.get(tag, self._add_unknown_tag)(
                    iterator, tag, data, config, pos)
            else:
                if tag == self._Element__name:
                    self._end_tag(tag, data, pos)
                    break
        return self

    def to_xml(self, w, **kwargs):
        if self._Element__ns:
            name = ':'.join((self._Element__ns, self._Element__name))
        else:
            name = self._Element__name

        with w.tag(name, attrib=object_attrs(self)):
            for name, formatter, child in object_children(self):
                if isinstance(child, Element):
                    child.to_xml(w, formatter=formatter)
                else:
                    if formatter:
                        child = formatter()

                    if not child:
                        child = ''

                    w.element(name, str(child))


class ElementWithXSIType(Element):
    """
    An XML element that supports type dispatch through xsi:type.

    When a class A is derived from this, it gains a decorator
    register_xsi_type, which classes derived from A can use to say
    "construct me rather than A when xsi:type has the value I'm putting
    in.

    At this point we are doing *no* namespace processing in our XML
    parsing.  Hence, we discard any prefixes both when registering
    and when matching.

    We probably should do namespaces one day; astropy.utils.xml will
    presumably learn them when they add VO-DML support.  Let's revisit
    this when it's there safely.  Meanwhere, use canonical Registry
    prefixes (cf. RegTAP 1.1, sect. 5) everywhere in your code; these
    will continue to be safe no matter what.
    """
    _xsi_type_mapping = {}

    @classmethod
    def register_xsi_type(cls, typename):
        """Decorator factory for registering subtypes."""
        def register(class_):
            """Decorator for registering subtypes"""
            cls._xsi_type_mapping[typename.split(":")[-1]] = class_
            return class_
        return register

    def __new__(cls, *args, **kwargs):
        xsi_type = None
        # Another namespace trouble: people can bind the xsi URI
        # to anything, *and* it's not unlikely they have another
        # type attribute, too.  Wiggle out of it by preferring a
        # literal xsi:type and otherwise hope for the best.  This
        # really needs to be fixed when we switch to namespace-aware
        # parsing.
        for name, val in kwargs.items():
            if name == "xsi:type":
                xsi_type = val
                break
            elif name.split(":")[-1] == "type":
                xsi_type = val

        if xsi_type is None:
            dtype = cls
        else:
            try:
                dtype = cls._xsi_type_mapping[xsi_type.split(":")[-1]]
            except KeyError:
                warnings.warn(f"Unknown xsi:type {xsi_type} ignored")
                dtype = cls

        obj = Element.__new__(dtype)
        obj.__init__(*args, **kwargs)
        return obj


class ContentMixin(Element):
    """
    Mixin class for elements with inner content.
    """

    def __init__(self, config=None, pos=None, _name=None, _ns=None, **kwargs):
        super().__init__(config, pos, _name, _ns, **kwargs)
        self._content = None

    def __bool__(self):
        return bool(self.content)

    __nonzero__ = __bool__

    def _end_tag(self, tag, data, pos):
        self.content = data

    def _content_check(self, content):
        pass

    def _content_parse(self, content):
        return content

    @property
    def content(self):
        """The inner content of the element."""
        return self._content

    @content.setter
    def content(self, content):
        self._content_check(content)
        self._content = self._content_parse(content)

    def to_xml(self, w, **kwargs):
        if self._Element__ns:
            name = ':'.join((self._Element__ns, self._Element__name))
        else:
            name = self._Element__name

        try:
            content = kwargs['formatter']()
        except (KeyError, TypeError):
            content = self.content

        if content is not None:
            w.element(name, str(content), attrib=object_attrs(self))
