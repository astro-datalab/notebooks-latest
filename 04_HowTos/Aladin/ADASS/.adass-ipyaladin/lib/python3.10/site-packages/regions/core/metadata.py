# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module proves classes to handle region metadata.
"""
from copy import deepcopy

__all__ = ['Meta', 'RegionMeta', 'RegionVisual']


class Meta(dict):
    """
    A base class for region metadata.

    Parameters
    ----------
    seq : dict-like, optional
        A dictionary or other mapping object to initialize the metadata.

    **kwargs
        Additional keyword arguments to initialize the metadata.
    """

    valid_keys = []
    key_mapping = {}

    def __init__(self, seq=None, **kwargs):
        super().__init__()

        if seq:
            if isinstance(seq, dict):
                for key, val in seq.items():
                    self.__setitem__(key, val)
            else:
                for key, val in seq:
                    self.__setitem__(key, val)

        if len(kwargs) > 0:
            for key, val in kwargs.items():
                self.__setitem__(key, val)

    def __setitem__(self, key, value):
        key = self.key_mapping.get(key, key)
        if key in self.valid_keys:
            super().__setitem__(key, value)
        else:
            raise KeyError(f'{key} is not a valid key for this class.')

    def __getitem__(self, item):
        item = self.key_mapping.get(item, item)
        return super().__getitem__(item)

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise ValueError('Only one argument can be input')
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]

    def copy(self):
        """
        Make a deep copy of this object.
        """
        return deepcopy(self)


class RegionMeta(Meta):
    """
    A dictionary subclass that holds the meta attributes of the region.
    """

    valid_keys = ['background', 'comment', 'component', 'composite', 'corr',
                  'delete', 'edit', 'fixed', 'frame', 'highlite', 'include',
                  'label', 'line', 'move', 'name', 'range', 'restfreq',
                  'rotate', 'select', 'source', 'tag', 'text', 'textrotate',
                  'type', 'veltype']

    key_mapping = {}


class RegionVisual(Meta):
    """
    A dictionary subclass which holds the visual attributes of the
    region.
    """

    valid_keys = ['color', 'dash', 'dashlist', 'fill', 'font',
                  'fontname', 'fontsize', 'fontstyle', 'fontweight',
                  'labeloff', 'labelpos', 'labelcolor', 'line', 'linestyle',
                  'linewidth', 'marker', 'markersize', 'symbol', 'symsize',
                  'symthick', 'textangle', 'textrotate', 'usetex',
                  'default_style', 'dashes', 'markeredgewidth', 'rotation',
                  'facecolor', 'edgecolor']

    key_mapping = {'point': 'symbol', 'width': 'linewidth'}

    def _define_default_mpl_kwargs(self, artist):
        """
        Define the default matplotlib kwargs for the specified artist.

        The kwargs depend on the value of self.visual['default_style'],
        which can be set when reading region files. If this keywords is
        not set or set to 'mpl' or `None`, then the matplotlib defaults
        will be used, with the exception fill is turned off for Patch
        and Line2D objects.

        Parameters
        ----------
        artist : {'Text', 'Line2D', 'Patch'}
            The matplotlib artist type.

        Returns
        -------
        result : dict
            A dictionary of matplotlib keyword arguments.
        """
        kwargs = {}

        default_style = self.get('default_style', None)
        if default_style is None or default_style == 'mpl':
            # do not fill by default, which is the only change from
            # matplotlib defaults
            if artist == 'Patch':
                kwargs['fill'] = False
            elif artist == 'Line2D':
                kwargs['fillstyle'] = 'none'
                kwargs['marker'] = 'o'
            return kwargs

        # 'ds9' style is set when reading from ds9 region files
        elif default_style == 'ds9':
            kwargs['color'] = '#00ff00'  # green
            if artist == 'Text':
                kwargs['ha'] = 'center'  # text horizontal alignment
                kwargs['va'] = 'center'  # text vertical alignment
            elif artist == 'Line2D':
                from regions.io.ds9.core import ds9_valid_symbols

                kwargs['marker'] = ds9_valid_symbols['boxcircle']
                kwargs['markersize'] = 11
                kwargs['markeredgecolor'] = kwargs.pop('color')
                kwargs['fillstyle'] = 'none'
            elif artist == 'Patch':
                kwargs['edgecolor'] = kwargs.pop('color')
                kwargs['fill'] = False

            else:
                raise ValueError('invalid visual["default"] value')

        return kwargs

    def _to_mpl_kwargs(self, artist):
        """
        Convert the visual metadata to a dictionary of matplotlib
        keyword arguments for the given artist.

        Parameters
        ----------
        artist : {'Text', 'Line2D', 'Patch'}
            The matplotlib artist type.

        Returns
        -------
        result : dict
            A dictionary of matplotlib keyword arguments.
        """
        if artist == 'Text':
            keymap = {'font': 'family',
                      'fontstyle': 'style',
                      'fontweight': 'weight',
                      'fontsize': 'size',
                      'textangle': 'rotation'}

        elif artist == 'Line2D':
            keymap = {'symsize': 'markersize',
                      'color': 'markeredgecolor',
                      'linewidth': 'markeredgewidth',
                      'fill': 'fillstyle'}

        elif artist == 'Patch':
            keymap = {'color': 'edgecolor',
                      'fill': 'fill'}

        else:
            raise ValueError('invalid artist type')

        kwargs = {}
        for name, val in self.items():
            if name in keymap:
                # NOTE: this will override existing mpl kwargs
                kwargs[keymap[name]] = val
            else:
                kwargs[name] = val

        default_style = kwargs.pop('default_style', None)
        if default_style == 'ds9':
            for key, val in kwargs.items():
                if val == 'green':
                    # X11/mpl green is #008000, ds9 uses #00ff00
                    kwargs[key] = '#00ff00'

        return kwargs

    def define_mpl_kwargs(self, artist):
        """
        Define a dictionary of matplotlib keywords for the input
        ``artist`` from the region's ``visual`` properties.

        Parameters
        ----------
        artist : {'Text', 'Line2D', 'Patch'}
            The matplotlib artist type.

        Returns
        -------
        result : dict
            A dictionary of matplotlib keyword arguments.
        """
        if artist not in ('Patch', 'Line2D', 'Text'):
            raise ValueError(f'artist "{artist}" is not supported')

        kwargs = self._define_default_mpl_kwargs(artist)
        kwargs.update(self._to_mpl_kwargs(artist))

        remove_keys = []
        if artist != 'Text':
            remove_keys.extend(['fontname', 'fontsize', 'fontweight',
                                'fontstyle'])
        else:
            remove_keys.extend(['linewidth'])

        for key in remove_keys:
            kwargs.pop(key, None)

        if 'fontstyle' in kwargs:
            kwargs['fontstyle'].replace('roman', 'normal')

        return kwargs
