# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines text regions in both pixel and sky coordinates.
"""

from regions._utils.wcs_helpers import pixel_scale_angle_at_skycoord
from regions.core.attributes import (RegionMetaDescr, RegionVisualDescr,
                                     ScalarPixCoord, ScalarSkyCoord)
from regions.shapes.point import PointPixelRegion, PointSkyRegion

__all__ = ['TextSkyRegion', 'TextPixelRegion']


class TextPixelRegion(PointPixelRegion):
    """
    A text string in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The leftmost point of the text string before rotation.
    text : str
        The text string.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import PixCoord, TextPixelRegion, RegionVisual
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        center = PixCoord(x=15, y=10)
        visual = RegionVisual({'textangle': 30})
        reg = TextPixelRegion(center=center, text="Hello World!",
                              visual=visual)
        reg.plot(ax=ax)

        ax.set_xlim(10, 30)
        ax.set_ylim(2.5, 20)
        ax.set_aspect('equal')
    """

    _params = ('center', 'text')
    _mpl_artist = 'Text'
    center = ScalarPixCoord('The leftmost pixel position (before rotation) '
                            'as a |PixCoord|.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, text, meta=None, visual=None):
        super().__init__(center, meta, visual)
        self.text = text

    def to_sky(self, wcs):
        center = wcs.pixel_to_world(self.center.x, self.center.y)

        # rotation value is relative to the coordinate system axes;
        # convert from counterclockwise angle from the positive x axis
        # to angle relative to WCS longitude axis
        visual = self.visual
        if 'rotation' in self.visual:
            _, _, angle = pixel_scale_angle_at_skycoord(center, wcs)
            visual = visual.copy()
            visual['rotation'] -= angle.to('deg').value - 90.

        return TextSkyRegion(center, self.text, meta=self.meta.copy(),
                             visual=visual.copy())

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib Text object for this region
        (`matplotlib.text.Text`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : dict
            Any keyword arguments accepted by `~matplotlib.text.Text`.
            These keywords will override any visual meta attributes of
            this region.

        Returns
        -------
        artist : `~matplotlib.text.Text`
            A matplotlib Text object.
        """
        from matplotlib.text import Text

        mpl_kwargs = self.visual.define_mpl_kwargs(self._mpl_artist)
        mpl_kwargs.update(kwargs)

        return Text(self.center.x - origin[0], self.center.y - origin[1],
                    self.text, **mpl_kwargs)


class TextSkyRegion(PointSkyRegion):
    """
    A text string in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The leftmost position of the text string before rotation.
    text : str
        The text string.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center', 'text')
    center = ScalarSkyCoord('The leftmost position (before rotation) as a '
                            '|SkyCoord|.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, text, meta=None, visual=None):
        super().__init__(center, meta, visual)
        self.text = text

    def to_pixel(self, wcs):
        center, _, angle = pixel_scale_angle_at_skycoord(self.center, wcs)

        # rotation value is relative to the WCS longitude axis;
        # convert to counterclockwise angle from the positive x axis
        visual = self.visual
        if 'rotation' in self.visual:
            visual = visual.copy()
            visual['rotation'] += angle.to('deg').value - 90.

        return TextPixelRegion(center, self.text, meta=self.meta.copy(),
                               visual=visual.copy())
