# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.utils.exceptions import AstropyUserWarning

__all__ = ['CRTFRegionParserWarning', 'CRTFRegionParserError']


class CRTFRegionParserWarning(AstropyUserWarning):
    """
    A generic warning class for CRTF region parsing.
    """


class CRTFRegionParserError(ValueError):
    """
    A generic error class for CRTF region parsing.
    """


# Valid symbols in CRTF
valid_symbols = {'.': 'point',
                 ',': 'pixel',
                 'o': 'circle',
                 'v': 'triangle_down',
                 '^': 'triangle_up',
                 '<': 'triangle_left',
                 '>': 'triangle_right',
                 '1': 'tri_down',
                 '2': 'tri_up',
                 '3': 'tri_left',
                 '4': 'tri_right',
                 's': 'square',
                 'p': 'pentagon',
                 '*': 'star',
                 'h': 'hexagon1',
                 'H': 'hexagon2',
                 '+': 'plus',
                 'x': 'x',
                 'D': 'diamond',
                 'd': 'thin_diamond',
                 '|': 'vline',
                 '_': 'hline'}
