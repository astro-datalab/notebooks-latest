# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
A collection of routines to format metadata
"""

import re
from itertools import chain
import textwrap


_parasp = re.compile(r"(?:[ \t\r\f\v]*\n){2,}[ \t\r\f\v]*")
_ptag = re.compile(r"\s*(?:<p\s*/?>)|(?:\\para(?:\\ )*)\s*")


def para_format_desc(text, width=78):
    """
    format description text into paragraphs suitable for display in the
    shell.  That is, the output will be one or more plain text paragraphs
    of the prescribed width (78 characters, the default).  The text will
    be split into separate paragraphs where there occurs (1) a two or more
    consecutive carriage return, (2) an HTML paragraph tag, or (2)
    a LaTeX paragraph control sequence.  It will attempt other substitutions
    of HTML and LaTeX markup that sometimes find their way into resource
    descriptions.
    """
    paras = _parasp.split(text)
    paras = filter(
        bool, chain.from_iterable(_ptag.split(para) for para in paras))
    paras = ("\n".join(
        map(lambda ll: ll.strip(), para.splitlines())
    ) for para in paras)
    paras = map(deref_markup, paras)

    return "\n\n".join(textwrap.fill(para, width) for para in paras)


_musubs = [
    (re.compile(r"&lt;"), "<"), (re.compile(r"&gt;"), ">"),
    (re.compile(r"&amp;"), "&"), (re.compile(r"<br\s*/?>"), ''),
    (re.compile(r"</p>"), ''), (re.compile(r"&#176;"), " deg"),
    (re.compile(r"\$((?:[^\$]*[\*\+=/^_~><\\][^\$]*)|(?:\w+))\$"), r'\1'),
    (re.compile(r"\\deg"), " deg"),
]

_alink = re.compile(r'''<a .*href=(["])([^\1]*)(?:\1).*>\s*(\S.*\S)\s*</a>''')


def deref_markup(text):
    """
    perform some substitutions of common markup suitable for text display.
    This includes HTML escape sequence
    """
    for pat, repl in _musubs:
        text = pat.sub(repl, text)
    text = _alink.sub(r"\3 <\2>", text)
    return text
