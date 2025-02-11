# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines circular regions in both pixel and sky coordinates.
"""

import math

import astropy.units as u
import numpy as np
from astropy.coordinates import Angle

from regions._geometry import circular_overlap_grid
from regions._utils.wcs_helpers import pixel_scale_angle_at_skycoord
from regions.core.attributes import (PositiveScalar, PositiveScalarAngle,
                                     RegionMetaDescr, RegionVisualDescr,
                                     ScalarPixCoord, ScalarSkyCoord)
from regions.core.bounding_box import RegionBoundingBox
from regions.core.core import PixelRegion, SkyRegion
from regions.core.mask import RegionMask
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord

__all__ = ['CirclePixelRegion', 'CircleSkyRegion']


class CirclePixelRegion(PixelRegion):
    """
    A circle defined using pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The center position.
    radius : float
        The radius in pixels.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import PixCoord, CirclePixelRegion
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        reg = CirclePixelRegion(PixCoord(x=8, y=7), radius=3.5)
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Circle')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(0, 15)
        ax.set_ylim(0, 15)
        ax.set_aspect('equal')
    """

    _params = ('center', 'radius')
    _mpl_artist = 'Patch'
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    radius = PositiveScalar('The radius in pixels as a float.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, radius, meta=None, visual=None):
        self.center = center
        self.radius = radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    @property
    def area(self):
        return math.pi * self.radius ** 2

    def contains(self, pixcoord):
        pixcoord = PixCoord._validate(pixcoord, name='pixcoord')
        in_circle = self.center.separation(pixcoord) < self.radius
        if self.meta.get('include', True):
            return in_circle
        else:
            return np.logical_not(in_circle)

    def to_sky(self, wcs):
        # TODO: write a pixel_to_skycoord_scale_angle
        center = wcs.pixel_to_world(self.center.x, self.center.y)
        _, pixscale, _ = pixel_scale_angle_at_skycoord(center, wcs)
        radius = Angle(self.radius * u.pix * pixscale, 'arcsec')
        return CircleSkyRegion(center, radius, meta=self.meta.copy(),
                               visual=self.visual.copy())

    @property
    def bounding_box(self):
        """
        Bounding box (`~regions.RegionBoundingBox`).
        """
        xmin = self.center.x - self.radius
        xmax = self.center.x + self.radius
        ymin = self.center.y - self.radius
        ymax = self.center.y + self.radius

        return RegionBoundingBox.from_float(xmin, xmax, ymin, ymax)

    def to_mask(self, mode='center', subpixels=1):
        self._validate_mode(mode, subpixels)

        if mode == 'center':
            mode = 'subpixels'
            subpixels = 1

        # Find bounding box and mask size
        bbox = self.bounding_box
        ny, nx = bbox.shape

        # Find position of pixel edges and recenter so that circle is at
        # origin
        xmin = float(bbox.ixmin) - 0.5 - self.center.x
        xmax = float(bbox.ixmax) - 0.5 - self.center.x
        ymin = float(bbox.iymin) - 0.5 - self.center.y
        ymax = float(bbox.iymax) - 0.5 - self.center.y

        use_exact = 0 if mode == 'subpixels' else 1

        fraction = circular_overlap_grid(xmin, xmax, ymin, ymax, nx, ny,
                                         self.radius, use_exact, subpixels)

        return RegionMask(fraction, bbox=bbox)

    def as_artist(self, origin=(0, 0), **kwargs):
        """
        Return a matplotlib patch object for this region
        (`matplotlib.patches.Circle`).

        Parameters
        ----------
        origin : array_like, optional
            The ``(x, y)`` pixel position of the origin of the displayed
            image.

        **kwargs : dict
            Any keyword arguments accepted by
            `~matplotlib.patches.Circle`. These keywords will override
            any visual meta attributes of this region.

        Returns
        -------
        artist : `~matplotlib.patches.Circle`
            A matplotlib circle patch.
        """
        from matplotlib.patches import Circle

        xy = self.center.x - origin[0], self.center.y - origin[1]
        radius = self.radius
        mpl_kwargs = self.visual.define_mpl_kwargs(self._mpl_artist)
        mpl_kwargs.update(kwargs)

        return Circle(xy=xy, radius=radius, **mpl_kwargs)

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
        region : `CirclePixelRegion`
            The rotated region (which is an independent copy).
        """
        center = self.center.rotate(center, angle)
        return self.copy(center=center)


class CircleSkyRegion(SkyRegion):
    """
    A circle defined using sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The center position.
    radius : `~astropy.units.Quantity`
        The radius in angular units.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center', 'radius')
    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    radius = PositiveScalarAngle('The radius as a |Quantity| angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, radius, meta=None, visual=None):
        self.center = center
        self.radius = radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

    def to_pixel(self, wcs):
        center, pixscale, _ = pixel_scale_angle_at_skycoord(self.center, wcs)
        radius = (self.radius / pixscale).to(u.pix).value
        return CirclePixelRegion(center, radius, meta=self.meta.copy(),
                                 visual=self.visual.copy())
