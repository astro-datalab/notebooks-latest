# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Tests for the metadata module.
"""
import pytest

from regions.core.metadata import RegionMeta, RegionVisual


def test_region_meta():
    meta_dict = {'text': 'hello world', 'tag': ['Tag1', 'Tag2']}
    meta = RegionMeta(meta_dict)

    for key, val in meta_dict.items():
        assert val == meta[key]

    meta2 = meta.copy()
    assert isinstance(meta2, RegionMeta)

    text = 'new'
    meta['text'] = text
    assert meta2['text'] != text

    with pytest.raises(KeyError):
        RegionMeta({'invalid': 1})
    with pytest.raises(KeyError):
        meta.update({'invalid': 1})
    with pytest.raises(KeyError):
        meta.update(invalid=1)
    with pytest.raises(KeyError):
        meta.setdefault('invalid', 1)


def test_region_visual():
    meta_dict = {'color': 'blue', 'fontsize': 12}
    meta = RegionVisual(meta_dict)

    for key, val in meta_dict.items():
        assert val == meta[key]

    meta2 = meta.copy()
    assert isinstance(meta2, RegionVisual)

    color = 'green'
    meta['color'] = color
    assert meta2['color'] != color

    with pytest.raises(KeyError):
        RegionVisual({'invalid': 1})
    with pytest.raises(KeyError):
        meta.update({'invalid': 1})
    with pytest.raises(KeyError):
        meta.update(invalid=1)
    with pytest.raises(KeyError):
        meta.setdefault('invalid', 1)


def test_region_visual_mpl_kwargs():
    meta_dict = {'color': 'blue'}
    meta = RegionVisual(meta_dict)

    kwargs = meta.define_mpl_kwargs('Patch')
    expected = {'edgecolor': 'blue', 'fill': False}
    assert kwargs == expected

    kwargs = meta.define_mpl_kwargs('Line2D')
    expected = {'markeredgecolor': 'blue', 'fillstyle': 'none',
                'marker': 'o'}
    assert kwargs == expected

    kwargs = meta.define_mpl_kwargs('Text')
    expected = {'color': 'blue'}
    assert kwargs == expected

    with pytest.raises(ValueError):
        meta.define_mpl_kwargs('invalid')
