# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines line regions in both pixel and sky coordinates.
"""

import numpy as np

from regions.core.attributes import (RegionMetaDescr, RegionVisualDescr,
                                     ScalarPixCoord, ScalarSkyCoord)
from regions.core.bounding_box import RegionBoundingBox
from regions.core.core import PixelRegion, SkyRegion
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord

__all__ = ['LinePixelRegion', 'LineSkyRegion']


class LinePixelRegion(PixelRegion):
    """
    A line in pixel coordinates.

    Parameters
    ----------
    start : `~regions.PixCoord`
        The start position.
    end : `~regions.PixCoord`
        The end position.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import PixCoord, LinePixelRegion
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        start = PixCoord(x=15, y=10)
        end = PixCoord(x=20, y=25)
        reg = LinePixelRegion(start=start, end=end)
        patch = reg.plot(ax=ax, edgecolor='red', lw=2, label='Line')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(0, 30)
        ax.set_ylim(0, 30)
        ax.set_aspect('equal')
    """

    _params = ('start', 'end')
    _mpl_artist = 'Patch'
    start = ScalarPixCoord('The start pixel position as a |PixCoord|.')
    end = ScalarPixCoord('The end pixel position as a |PixCoord|.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, start, end, meta=None, visual=None):
        self.start = start
        self.end = end
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    @property
    def area(self):
        return 0

    def contains(self, pixcoord):
        in_reg = (False if pixcoord.isscalar
                  else np.zeros(pixcoord.x.shape, dtype=bool))

        if self.meta.get('include', True):
            return in_reg
        else:
            return np.logical_not(in_reg)

    def to_sky(self, wcs):
        start = wcs.pixel_to_world(self.start.x, self.start.y)
        end = wcs.pixel_to_world(self.end.x, self.end.y)
        return LineSkyRegion(start, end, meta=self.meta.copy(),
                             visual=self.visual.copy())

    @property
    def bounding_box(self):
        xmin = min(self.start.x, self.end.x)
        xmax = max(self.start.x, self.end.x)
        ymin = min(self.start.y, self.end.y)
        ymax = max(self.start.y, self.end.y)

        return RegionBoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=5):
        # TODO: needs to be implemented
        raise NotImplementedError

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib patch object for this region
        (`matplotlib.patches.Arrow`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : dict
            Any keyword arguments accepted by
            `~matplotlib.patches.Arrow`. These keywords will override
            any visual meta attributes of this region.

        Returns
        -------
        artist : `~matplotlib.patches.Arrow`
            A matplotlib line patch.
        """
        # Note: Long term we want to support DS9 lines with arrow heads.
        # We may want to use Line2D instead of arrow for lines because the
        # width of the arrow is non-scalable in patches.
        from matplotlib.patches import Arrow

        x = self.start.x - origin[0]
        y = self.start.y - origin[1]
        dx = self.end.x - self.start.x
        dy = self.end.y - self.start.y
        kwargs.setdefault('width', 0.1)

        mpl_kwargs = self.visual.define_mpl_kwargs(self._mpl_artist)
        mpl_kwargs.update(kwargs)

        return Arrow(x, y, dx, dy, **mpl_kwargs)

    def rotate(self, center, angle):
        """
        Rotate the region.

        Positive ``angle`` corresponds to counter-clockwise rotation.

        Parameters
        ----------
        center : `~regions.PixCoord`
            The rotation center point.
        angle : `~astropy.coordinates.Angle`
            The rotation angle.

        Returns
        -------
        region : `LinePixelRegion`
            The rotated region (which is an independent copy).
        """
        start = self.start.rotate(center, angle)
        end = self.end.rotate(center, angle)
        return self.copy(start=start, end=end)


class LineSkyRegion(SkyRegion):
    """
    A line in sky coordinates.

    Parameters
    ----------
    start : `~astropy.coordinates.SkyCoord`
        The start position.
    end : `~astropy.coordinates.SkyCoord`
        The end position.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('start', 'end')
    start = ScalarSkyCoord('The start position as a |SkyCoord|.')
    end = ScalarSkyCoord('The end position as a |SkyCoord|.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, start, end, meta=None, visual=None):
        self.start = start
        self.end = end
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def contains(self, skycoord, wcs):  # pylint: disable=unused-argument
        # lines never contain anything
        return not self.meta.get('include', True)

    def to_pixel(self, wcs):
        start_x, start_y = wcs.world_to_pixel(self.start)
        start = PixCoord(start_x, start_y)
        end_x, end_y = wcs.world_to_pixel(self.end)
        end = PixCoord(end_x, end_y)
        return LinePixelRegion(start, end, meta=self.meta.copy(),
                               visual=self.visual.copy())
