# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains xml element classes as defined in the VOResource standard.

There are different ways of handling the various xml tags.

* Elements with complex content
* Elements with simple content and attributes
* Elements with simple content without attributes

Elements with complex content are parsed with objects inherited from `~pyvo.utils.xml.elements.Element`.

Elements with simple content are parsed with objects inherited from `~pyvo.utils.xml.elements.Element`
defining a ``value`` property.
"""
from astropy.utils.collections import HomogeneousList
from textwrap import indent

from ...utils.xml.elements import (
    Element, ElementWithXSIType, ContentMixin, xmlattribute, xmlelement)
from .exceptions import W06

__all__ = [
    "ValidationLevel", "Capability", "Interface", "AccessURL",
    "SecurityMethod", "WebBrowser", "WebService", "MirrorURL"]


######################################################################
# ELEMENT CLASSES
class ValidationLevel(ContentMixin, Element):
    """
    ValidationLevel element as described in
    http://www.ivoa.net/xml/VOResource/v1.0

    the allowed values for describing the resource descriptions and interfaces.
    See the RM (v1.1, section 4) for more guidance on the use of these values.

    Possible values:

    0:
        The resource has a description that is stored in a
        registry. This level does not imply a compliant
        description.
    1:
        In addition to meeting the level 0 definition, the
        resource description conforms syntactically to this
        standard and to the encoding scheme used.
    2:
        In addition to meeting the level 1 definition, the
        resource description refers to an existing resource that
        has demonstrated to be functionally compliant.

        When the resource is a service, it is consider to exist
        and functionally compliant if use of the
        service accessURL responds without error when used as
        intended by the resource. If the service is a standard
        one, it must also demonstrate the response is syntactically
        compliant with the service standard in order to be
        considered functionally compliant. If the resource is
        not a service, then the ReferenceURL must be shown to
        return a document without error.
    3:
        In addition to meeting the level 2 definition, the
        resource description has been inspected by a human and
        judged to comply semantically to this standard as well
        as meeting any additional minimum quality criteria (e.g.,
        providing values for important but non-required
        metadata) set by the human inspector.
    4:
        In addition to meeting the level 3 definition, the
        resource description meets additional quality criteria
        set by the human inspector and is therefore considered
        an excellent description of the resource. Consequently,
        the resource is expected to be operate well as part of a
        VO application or research study.
    """
    def __init__(
        self, config=None, pos=None, _name='validationLevel', validatedBy=None,
        **kwargs
    ):
        super().__init__(config, pos, _name, **kwargs)

        self._validatedby = validatedBy

    def __repr__(self):
        return '<ValidationLevel validatedBy={}>{}</ValidationLevel>'.format(
            self.validatedby, self.content)

    @xmlattribute
    def validatedby(self):
        """
        The IVOA ID of the registry or organisation that
        assigned the validation level.
        """
        return self._validatedby

    @validatedby.setter
    def validatedby(self, validatedby):
        self._validatedby = validatedby


class AccessURL(ContentMixin, Element):
    """
    AccessURL element as described in
    http://www.ivoa.net/xml/VOResource/v1.0

    The URL (or base URL) that a client uses to access the
    service.  How this URL is to be interpreted and used
    depends on the specific Interface subclass
    """
    def __init__(
        self, config=None, pos=None, _name='accessURL', use=None, **kwargs
    ):
        super().__init__(config, pos, _name, **kwargs)

        self._use = use

    def __repr__(self):
        return '<AccessURL use={}>{}</AccessURL>'.format(
            self.use, self.content)

    @xmlattribute
    def use(self):
        """
        A flag indicating whether this should be interpreted as a base
        URL, a full URL, or a URL to a directory that will produce a
        listing of files.

        Possible values:

        full:
            Assume a full URL--that is, one that can be invoked directly
            without alteration.  This usually returns a single document or
            file.
        base:
            Assume a base URL--that is, one requiring an extra portion to be
            appended before being invoked.
        dir:
            Assume URL points to a directory that will return a listing of
            files.
        """
        return self._use

    @use.setter
    def use(self, use):
        self._use = use


class MirrorURL(ContentMixin, Element):
    """
    A URL available as a mirror of an access URL.

    These come with a human-readable title intended to aid in mirror
    selection.
    """
    def __init__(
        self, config=None, pos=None, _name='accessURL', title=None, **kwargs
    ):
        super().__init__(config, pos, _name, **kwargs)
        self._title = title

    @xmlattribute
    def title(self):
        """
        A human-readable title for the mirror.
        """
        return self._title


class SecurityMethod(ContentMixin, Element):
    """
    SecurityMethod element as described in
    http://www.ivoa.net/xml/VOResource/v1.0

    A description of a security mechanism.

    this type only allows one to refer to the mechanism via a URI.
    Derived types would allow for more metadata.
    """
    def __init__(
        self, config=None, pos=None, _name='securityMethod', standardID=None,
        **kwargs
    ):
        super().__init__(config, pos, _name, **kwargs)

        self._standardid = standardID

    def __repr__(self):
        return '<SecurityMethod standardID={}>{}</SecurityMethod>'.format(
            self.standardid, self.content)

    @xmlattribute(name='standardID')
    def standardid(self):
        """
        A URI identifier for a standard security mechanism.
        """
        return self._standardid

    @standardid.setter
    def standardid(self, standardid):
        self._standardid = standardid


class Interface(ElementWithXSIType):
    """
    Interface element as described in
    http://www.ivoa.net/xml/VOResource/v1.0

    A description of a service interface.

    Since this type is abstract, one must use an Interface subclassto describe
    an actual interface.

    Additional interface subtypes (beyond WebService and WebBrowser) are
    defined in the VODataService schema.
    """
    def __init__(
        self, config=None, pos=None, _name='interface', version='1.0',
        role=None, **kwargs
    ):
        super().__init__(config, pos, _name, **kwargs)

        self._xsi_type = kwargs.get('xsi:type')

        self._version = version
        self._role = role
        self._resulttype = None
        self._testquerystring = None

        self._accessurls = HomogeneousList(AccessURL)
        self._securitymethods = HomogeneousList(SecurityMethod)
        self._mirrorurls = HomogeneousList(MirrorURL)

    def __repr__(self):
        return '<Interface role={}>...</Interface>'.format(
            self.role)

    def describe(self):
        """
        Prints out a human readable description
        """
        print('Interface {}'.format(self._xsi_type))

        accessurls = '\n'.join(
            accessurl.content for accessurl in self.accessurls)

        print(indent(accessurls, 4 * " "))

        print()

    @xmlattribute
    def version(self):
        """
        The version of a standard interface specification that this
        interface complies with.  When the interface is
        provided in the context of a Capability element, then
        the standard being refered to is the one identified by
        the Capability's standardID element.  If the standardID
        is not provided, the meaning of this attribute is
        undefined.
        """
        return self._version

    @version.setter
    def version(self, version):
        self._version = version

    @xmlattribute
    def role(self):
        """
        A tag name the identifies the role the interface plays
        in the particular capability.  If the value is equal to
        "std" or begins with "std:", then the interface refers
        to a standard interface defined by the standard
        referred to by the capability's standardID attribute.

        For an interface complying with some registered
        standard (i.e. has a legal standardID), the role can be
        match against interface roles enumerated in standard
        resource record.  The interface descriptions in
        the standard record can provide default descriptions
        so that such details need not be repeated here.
        """
        return self._role

    @role.setter
    def role(self, role):
        self._role = role

    @xmlelement(name='accessURL', cls=AccessURL)
    def accessurls(self):
        """
        A list of access urls in the interface.  Must contain only `AccessURL`
        objects.
        """
        return self._accessurls

    @xmlelement(name='mirrorURL', cls=MirrorURL)
    def mirrorurls(self):
        """
        mirror(s) for this access URL.
        """
        return self._mirrorurls

    @xmlelement(name='securityMethod', cls=SecurityMethod)
    def securitymethods(self):
        """
        the mechanism the client must employ to gain secure
        access to the service.

        when more than one method is listed, each one must
        be employed to gain access.
        """
        return self._securitymethods

    @xmlelement(name='testQueryString')
    def testquerystring(self):
        """
        a string to be used in an interface-specific way to obtain a
        non-empty result from the service.
        """
        return self._testquerystring

    @testquerystring.setter
    def testquerystring(self, testquerystring):
        self._testquerystring = testquerystring

    @xmlelement
    def resulttype(self):
        """
        The MIME type of a document returned in the HTTP response.
        """
        return self._resulttype

    @resulttype.setter
    def resulttype(self, resulttype):
        self._resulttype = resulttype


class Capability(ElementWithXSIType):
    """
    Capability element as described in
    http://www.ivoa.net/xml/VOResource/v1.0

    a description of what the service does
    (in terms of context-specific behavior), and how to use it
    (in terms of an interface)
    """
    def __init__(
        self, config=None, pos=None, _name='capability', standardID=None,
        **kwargs
    ):
        super().__init__(config, pos, _name, **kwargs)

        self._description = None
        self._standardid = standardID

        self._validationlevels = HomogeneousList(ValidationLevel)
        self._interfaces = HomogeneousList(Interface)

    def __repr__(self):
        return (
            '<Capability standardID={}>'
            '... {} validationLevels, {} interfaces ...'
            '</Capability>'
        ).format(
            self.standardid, len(self.validationlevels), len(self.interfaces))

    def describe(self):
        """
        Prints out a human readable description
        """
        print("Capability {}".format(self.standardid))
        print()

        if self.description:
            print(self.description)
            print()

        for interface in self.interfaces:
            interface.describe()

    @xmlelement(plain=True, multiple_exc=W06)
    def description(self):
        """
        A human-readable description of what this capability
        provides as part of the over-all service

        Use of this optional element is especially encouraged when
        this capability is non-standard and is one of several
        capabilities listed.
        """
        return self._description

    @description.setter
    def description(self, description):
        self._description = description

    @xmlelement(name='validationLevel', cls=ValidationLevel)
    def validationlevels(self):
        """
        A numeric grade describing the quality of the
        capability description and interface, when applicable,
        to be used to indicate the confidence an end-user
        can put in the resource as part of a VO application
        or research study.
        """
        return self._validationlevels

    @xmlelement(name='interface', cls=Interface)
    def interfaces(self):
        """
        a description of how to call the service to access this capability

        Since the Interface type is abstract, one must describe
        the interface using a subclass of Interface, denoting
        it via xsi:type.

        Multiple occurances can describe different interfaces to
        the logically same capability--i.e. data or functionality.
        That is, the inputs accepted and the output provides should
        be logically the same.  For example, a WebBrowser interface
        given in addition to a WebService interface would simply
        provide an interactive, human-targeted interface to the
        underlying WebService interface.
        """
        return self._interfaces

    @xmlattribute(name='standardID')
    def standardid(self):
        """
        A URI identifier for a standard service.

        This provides a unique way to refer to a service
        specification standard, such as a Simple Image Access service.
        The use of an IVOA identifier here implies that a
        VOResource description of the standard is registered and
        accessible.
        """
        return self._standardid

    @standardid.setter
    def standardid(self, standardid):
        self._standardid = standardid


@Interface.register_xsi_type('vr:WebBrowser')
class WebBrowser(Interface):
    """
    WebBrowser element as described in
    http://www.ivoa.net/xml/VOResource/v1.0

    A (form-based) interface intended to be accesed interactively by a user via
    a web browser. The accessURL represents the URL of the web form itself.
    """


@Interface.register_xsi_type('vr:WebService')
class WebService(Interface):
    """
    WebService element as described in
    http://www.ivoa.net/xml/VOResource/v1.0

    A Web Service that is describable by a WSDL document.
    The accessURL element gives the Web Service's endpoint URL.
    """
    def __init__(self, config=None, pos=None, _name='interface', **kwargs):
        super().__init__(config, pos, _name, **kwargs)

        self._wsdlurls = HomogeneousList(str)

    @xmlelement(name='wsdlURL')
    def wsdlurls(self):
        """
        The location of the WSDL that describes this Web Service.
        If not provided, the location is assumed to be the accessURL with
        "?wsdl" appended.

        Multiple occurances should represent mirror copies of the same WSDL
        file.
        """
        return self._wsdlurls
