# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.exceptions import AstropyWarning

__all__ = ['XMLWarning', 'UnknownElementWarning']


def _format_message(message, name, config=None, pos=None):
    if config is None:
        config = {}
    if pos is None:
        pos = ('?', '?')
    filename = config.get('filename', '?')

    return '{}:{}:{}: {}: {}'.format(filename, pos[0], pos[1], name, message)


class XMLWarning(AstropyWarning):
    """
    Base warning for violations of XML specifications
    """
    def __init__(self, args, config=None, pos=None):
        if config is None:
            config = {}
        if not isinstance(args, tuple):
            args = (args, )
        msg = self.message_template.format(*args)

        self.formatted_message = _format_message(
            msg, self.__class__.__name__, config, pos)
        Warning.__init__(self, self.formatted_message)


class UnknownElementWarning(XMLWarning):
    """
    Warning for missing xml elements
    """
    message_template = "Unknown element {}"
    default_args = ('x',)
