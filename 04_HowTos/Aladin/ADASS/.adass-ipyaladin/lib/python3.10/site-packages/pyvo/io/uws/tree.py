# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This file contains xml element classes as defined in the VOResource standard.
"""
from functools import partial

from astropy.utils.collections import HomogeneousList
from astropy.time import Time, TimeDelta

from ...utils.xml.elements import (
    xmlattribute, xmlelement, Element, ContentMixin)

uwselement = partial(xmlelement, ns='uws')


def XSInDate(val):
    if not val:
        return None

    try:
        return Time(val, format='iso')
    except ValueError:
        pass

    try:
        return Time(val, format='isot')
    except ValueError:
        pass

    raise ValueError('Cannot parse datetime {}'.format(val))


InDuration = partial(TimeDelta, format='sec')
XSOutDate = partial(Time, out_subfmt='date')


__all__ = [
    'UWSElement', 'Reference', 'JobSummary', 'Parameters', 'Parameter',
    'Results', 'Result']


def _convert_boolean(value, default=None):
    return {
        'false': False,
        '0': False,
        'true': True,
        '1': True
    }.get(value, default)


class UWSElement(Element):
    def __init__(self, config=None, pos=None, _name='', _ns='uws', **kwargs):
        super().__init__(config, pos, _name, 'uws', **kwargs)


class Reference(UWSElement):
    """standard xlink references"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.type = kwargs.get('xlink:type')
        self.href = kwargs.get('xlink:href')

    @xmlattribute(name='xlink:type')
    def type(self):
        """the type of the result"""
        return self._type

    @type.setter
    def type(self, type_):
        self._type = type_

    @xmlattribute(name='xlink:href')
    def href(self):
        """the url the result can be retrieved"""
        return self._href

    @href.setter
    def href(self, href):
        self._href = href


class JobSummary(Element):
    def __init__(self, config=None, pos=None, _name='job', **kwargs):
        super().__init__(config, pos, _name, **kwargs)
        self.jobid = kwargs.get('id')
        self._runid = None
        self._ownerid = None
        self._phase = None
        self._quote = None
        self._creationtime = None
        self._starttime = None
        self._endtime = None
        self._executionduration = None
        self._destruction = None
        self._parameters = Parameters()
        self._results = Results()
        self._errorsummary = None
        self._message = None

    @uwselement(name='jobId', plain=True)
    def jobid(self):
        """
        The identifier for the job
        """
        return self._jobid

    @jobid.setter
    def jobid(self, jobid):
        self._jobid = jobid

    @uwselement(name='runId', plain=True)
    def runid(self):
        """client supplied identifier"""
        return self._runid

    @runid.setter
    def runid(self, runid):
        self._runid = runid

    @uwselement(name='ownerId', plain=True)
    def ownerid(self):
        """the owner (creator) of the job"""
        return self._ownerid

    @ownerid.setter
    def ownerid(self, ownerid):
        self._ownerid = ownerid

    @uwselement(plain=True)
    def phase(self):
        """the execution phase"""
        return self._phase

    @phase.setter
    def phase(self, phase):
        self._phase = phase

    @uwselement(plain=True)
    def quote(self):
        """estimated completion time"""
        return self._quote

    @quote.setter
    def quote(self, quote):
        self._quote = XSInDate(quote)

    @quote.formatter
    def quote(self):
        try:
            return str(XSOutDate(self._quote))
        except ValueError:
            return None

    @uwselement(name='creationTime', plain=True)
    def creationtime(self):
        """The instant at which the job was created."""
        return self._creationtime

    @creationtime.setter
    def creationtime(self, creationtime):
        self._creationtime = XSInDate(creationtime)

    @creationtime.formatter
    def creationtime(self):
        try:
            return str(XSOutDate(self._creationtime))
        except ValueError:
            return None

    @uwselement(name='startTime', plain=True)
    def starttime(self):
        """The instant at which the job started execution."""
        return self._starttime

    @starttime.setter
    def starttime(self, starttime):
        self._starttime = XSInDate(starttime)

    @starttime.formatter
    def starttime(self):
        try:
            return str(XSOutDate(self._starttime))
        except ValueError:
            return None

    @uwselement(name='endTime', plain=True)
    def endtime(self):
        """The instant at which the job finished execution"""
        return self._endtime

    @endtime.setter
    def endtime(self, endtime):
        self._endtime = XSInDate(endtime)

    @endtime.formatter
    def endtime(self):
        try:
            return str(XSOutDate(self._endtime))
        except ValueError:
            return None

    @uwselement(name='executionDuration', plain=True)
    def executionduration(self):
        """
        The duration (in seconds) for which the job should be allowed to run -
        a value of 0 is intended to mean unlimited
        """
        return self._executionduration

    @executionduration.setter
    def executionduration(self, executionduration):
        if not isinstance(executionduration, TimeDelta):
            executionduration = InDuration(float(executionduration))

        self._executionduration = executionduration

    @executionduration.formatter
    def executionduration(self):
        if self.executionduration:
            return str(int(self._executionduration.value))

    @uwselement(plain=True)
    def destruction(self):
        """The time at which the whole job will be destroyed"""
        return self._destruction

    @destruction.setter
    def destruction(self, destruction):
        self._destruction = XSInDate(destruction)

    @destruction.formatter
    def destruction(self):
        try:
            return str(XSOutDate(self._destruction))
        except ValueError:
            return None

    @uwselement
    def parameters(self):
        """The parameters to the job"""
        return self._parameters

    @parameters.adder
    def parameters(self, iterator, tag, data, config, pos):
        parameters = Parameters(config, pos, 'parameters', **data)
        parameters.parse(iterator, config)
        self._parameters = parameters

    @uwselement
    def results(self):
        """The results for the job"""
        return self._results

    @results.adder
    def results(self, iterator, tag, data, config, pos):
        results = Results(config, pos, 'results', **data)
        results.parse(iterator, config)
        self._results = results

    @uwselement(name='errorSummary', plain=True)
    def errorsummary(self):
        """The error summary of the job."""
        return self._errorsummary

    @errorsummary.adder
    def errorsummary(self, iterator, tag, data, config, pos):
        res = ErrorSummary(config, pos, 'errorSummary', **data)
        res.parse(iterator, config)
        self._errorsummary = res


class Jobs(HomogeneousList, UWSElement):
    """A parsed representation of the joblist endpoint.
    """
    def __init__(self, config=None, pos=None, _name='jobs', **kwargs):
        HomogeneousList.__init__(self, JobSummary)
        UWSElement.__init__(self, config, pos, _name, **kwargs)

    @uwselement
    def jobs(self):
        return self

    @jobs.adder
    def jobs(self, iterator, tag, data, config, pos):
        return

    @uwselement(name='jobref')
    def joblist(self):
        return self

    @joblist.adder
    def joblist(self, iterator, tag, data, config, pos):
        job = JobSummary(config, pos, 'jobref', **data)
        job.parse(iterator, config)
        self.append(job)


class Parameters(UWSElement, HomogeneousList):
    """
    Parameters element of a job
    """
    def __init__(self, config=None, pos=None, _name='parameters', **kwargs):
        """ """
        # Note: Above is a load-bearing empty comment.
        # Do not remove, or else the Sphinx build may fail (see PR #193).
        HomogeneousList.__init__(self, Parameter)
        UWSElement.__init__(self, config, pos, _name, **kwargs)

    @uwselement(name='parameter')
    def parameters(self):
        return self

    @parameters.adder
    def parameters(self, iterator, tag, data, config, pos):
        parameter = Parameter(config, pos, 'parameter', **data)
        parameter.parse(iterator, config)
        self.append(parameter)


class Parameter(ContentMixin, UWSElement):
    def __init__(self, config=None, pos=None, _name='parameter', **kwargs):
        super().__init__(config, pos, _name, **kwargs)

        self.byreference = _convert_boolean(kwargs.get('byReference'))
        self.id_ = kwargs.get('id')

    @xmlattribute
    def byreference(self):
        """
        if this attribute is true then the content of the parameter represents
        a URL to retrieve the actual parameter value.
        """
        return self._byreference

    @byreference.setter
    def byreference(self, byreference):
        self._byreference = byreference

    @xmlattribute(name='id')
    def id_(self):
        """the identifier for the parameter"""
        return self._id

    @id_.setter
    def id_(self, id_):
        self._id = id_


class Results(UWSElement, HomogeneousList):
    """ """
    def __init__(self, config=None, pos=None, _name='results', **kwargs):
        HomogeneousList.__init__(self, Result)
        UWSElement.__init__(self, config, pos, _name, **kwargs)

    @uwselement(name='result')
    def results(self):
        return self

    @results.adder
    def results(self, iterator, tag, data, config, pos):
        result = Result(config, pos, 'result', **data)
        result.parse(iterator, config)
        self.append(result)


class Result(Reference, UWSElement):
    """A reference to a UWS result."""
    def __init__(self, config=None, pos=None, _name='result', **kwargs):
        super().__init__(config, pos, _name, **kwargs)

        self.id_ = kwargs.get('id')
        self.size = int(kwargs.get('size') or 0)
        self.mimetype = kwargs.get('mime-type')

    @xmlattribute(name='id')
    def id_(self):
        """the identifier for the result"""
        return self._id

    @id_.setter
    def id_(self, id_):
        self._id = id_

    @xmlattribute
    def size(self):
        """the size of the result"""
        return self._size

    @size.setter
    def size(self, size):
        self._size = size

    @xmlattribute
    def mimetype(self):
        """the mimetype of the result"""
        return self._mimetype

    @mimetype.setter
    def mimetype(self, mimetype):
        self._mimetype = mimetype


class ErrorSummary(UWSElement):
    """A UWS Error summary."""
    def __init__(self, config=None, pos=None, _name='errorSummary', **kwargs):
        super().__init__(config, pos, _name, **kwargs)

        self.type_ = kwargs.get('type')
        self.has_detail = _convert_boolean(kwargs.get('hasDetail'))
        self.message = None

    @xmlattribute(name='type')
    def type_(self):
        """the type of the error"""
        return self._type

    @type_.setter
    def type_(self, type_):
        self._type = type_

    @xmlattribute
    def has_detail(self):
        """whether error has details"""
        return self._has_detail

    @has_detail.setter
    def has_detail(self, has_detail):
        self._has_detail = has_detail

    @uwselement(name='message')
    def message(self):
        """The error message"""
        return self._message

    @message.setter
    def message(self, message):
        self._message = message


class Message(ContentMixin, UWSElement):
    """The actual UWS Error message."""
    def __init__(self, config=None, pos=None, _name='message', **kwargs):
        super().__init__(config, pos, _name, **kwargs)
