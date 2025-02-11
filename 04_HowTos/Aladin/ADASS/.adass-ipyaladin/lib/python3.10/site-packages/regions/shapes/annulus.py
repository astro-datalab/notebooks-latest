# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This module defines annulus regions in both pixel and sky coordinates.
"""

import abc
import operator

import astropy.units as u

from regions._utils.wcs_helpers import pixel_scale_angle_at_skycoord
from regions.core.attributes import (PositiveScalar, PositiveScalarAngle,
                                     RegionMetaDescr, RegionVisualDescr,
                                     ScalarAngle, ScalarPixCoord,
                                     ScalarSkyCoord)
from regions.core.compound import CompoundPixelRegion
from regions.core.core import PixelRegion, SkyRegion
from regions.core.metadata import RegionMeta, RegionVisual
from regions.core.pixcoord import PixCoord
from regions.shapes.circle import CirclePixelRegion
from regions.shapes.ellipse import EllipsePixelRegion, EllipseSkyRegion
from regions.shapes.rectangle import RectanglePixelRegion, RectangleSkyRegion

__all__ = ['AnnulusPixelRegion', 'AsymmetricAnnulusPixelRegion',
           'AsymmetricAnnulusSkyRegion',
           'CircleAnnulusPixelRegion', 'CircleAnnulusSkyRegion',
           'EllipseAnnulusPixelRegion', 'EllipseAnnulusSkyRegion',
           'RectangleAnnulusPixelRegion', 'RectangleAnnulusSkyRegion']


class AnnulusPixelRegion(PixelRegion, abc.ABC):
    """
    A base class for annulus pixel regions.
    """

    @property
    def _compound_region(self):
        return CompoundPixelRegion(self._inner_region, self._outer_region,
                                   operator.xor, self.meta, self.visual)

    @property
    def area(self):
        return self._outer_region.area - self._inner_region.area

    @property
    def bounding_box(self):
        return self._outer_region.bounding_box

    def contains(self, pixcoord):
        return self._compound_region.contains(pixcoord)

    def as_artist(self, origin=(0, 0), **kwargs):
        return self._compound_region.as_artist(origin, **kwargs)

    def to_mask(self, mode='center', subpixels=5):
        return self._compound_region.to_mask(mode, subpixels)

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
        region : `~regions.AnnulusPixelRegion`
            The rotated region (which is an independent copy).
        """
        changes = {}
        changes['center'] = self.center.rotate(center, angle)
        if hasattr(self, 'angle'):
            changes['angle'] = self.angle + angle
        return self.copy(**changes)


class CircleAnnulusPixelRegion(AnnulusPixelRegion):
    """
    A circular annulus in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the annulus.
    inner_radius : float
        The inner radius of the annulus in pixels.
    outer_radius : float
        The outer radius of the annulus in pixels.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from regions import PixCoord, CircleAnnulusPixelRegion
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        reg = CircleAnnulusPixelRegion(PixCoord(x=6, y=6),
                                       inner_radius=5.5,
                                       outer_radius=8.0)
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Circle Annulus')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(-5, 20)
        ax.set_ylim(-5, 20)
        ax.set_aspect('equal')
    """

    _component_class = CirclePixelRegion
    _params = ('center', 'inner_radius', 'outer_radius')
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    inner_radius = PositiveScalar('The inner radius in pixels as a float.')
    outer_radius = PositiveScalar('The outer radius in pixels as a float.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, inner_radius, outer_radius, meta=None,
                 visual=None):
        self.center = center
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

        if inner_radius >= outer_radius:
            raise ValueError('outer_radius must be greater than inner_radius')

    @property
    def _inner_region(self):
        return self._component_class(self.center, self.inner_radius,
                                     self.meta, self.visual)

    @property
    def _outer_region(self):
        return self._component_class(self.center, self.outer_radius,
                                     self.meta, self.visual)

    def to_sky(self, wcs):
        center = wcs.pixel_to_world(self.center.x, self.center.y)
        _, pixscale, _ = pixel_scale_angle_at_skycoord(center, wcs)
        inner_radius = self.inner_radius * u.pix * pixscale
        outer_radius = self.outer_radius * u.pix * pixscale
        return CircleAnnulusSkyRegion(center, inner_radius, outer_radius,
                                      self.meta.copy(), self.visual.copy())


class CircleAnnulusSkyRegion(SkyRegion):
    """
    A circular annulus in sky coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the annulus.
    inner_radius : `~astropy.units.Quantity`
        The inner radius of the annulus in angular units.
    outer_radius : `~astropy.units.Quantity`
        The outer radius of the annulus in angular units.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center', 'inner_radius', 'outer_radius')
    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    inner_radius = PositiveScalarAngle('The inner radius as a |Quantity| '
                                       'angle.')
    outer_radius = PositiveScalarAngle('The outer radius as a |Quantity| '
                                       'angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, inner_radius, outer_radius, meta=None,
                 visual=None):
        self.center = center
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

        if inner_radius >= outer_radius:
            raise ValueError('outer_radius must be greater than inner_radius')

    def to_pixel(self, wcs):
        center, pixscale, _ = pixel_scale_angle_at_skycoord(self.center, wcs)
        inner_radius = (self.inner_radius / pixscale).to(u.pix).value
        outer_radius = (self.outer_radius / pixscale).to(u.pix).value
        return CircleAnnulusPixelRegion(center, inner_radius, outer_radius,
                                        meta=self.meta.copy(),
                                        visual=self.visual.copy())


class AsymmetricAnnulusPixelRegion(AnnulusPixelRegion):
    """
    Helper class for asymmetric annuli sky regions.

    Used for ellipse and rectangle pixel annulus regions.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the annulus.
    center : `~regions.PixCoord`
        The position of the center of the annulus.
    inner_width : float
        The inner width of the annulus (before rotation) in pixels.
    outer_width : float
        The outer width of the annulus (before rotation) in pixels.
    inner_height : float
        The inner height of the annulus (before rotation) in pixels.
    outer_height : float
        The outer height of the annulus (before rotation) in pixels.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the annulus, measured anti-clockwise. If
        set to zero (the default), the width axis is lined up with the x
        axis.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center', 'inner_width', 'outer_width', 'inner_height',
               'outer_height', 'angle')
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    inner_width = PositiveScalar('The inner width (before rotation) in '
                                 'pixels as a float.')
    outer_width = PositiveScalar('The outer width (before rotation) in '
                                 'pixels as a float.')
    inner_height = PositiveScalar('The inner height (before rotation) in '
                                  'pixels as a float.')
    outer_height = PositiveScalar('The outer height (before rotation) in '
                                  'pixels as a float.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        self.center = center
        self.inner_width = inner_width
        self.outer_width = outer_width
        self.inner_height = inner_height
        self.outer_height = outer_height
        self.angle = angle
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

        if inner_width >= outer_width:
            raise ValueError('outer_width must be greater than inner_width')
        if inner_height >= outer_height:
            raise ValueError('outer_height must be greater than inner_height')

    @property
    def _inner_region(self):
        return self._component_class(self.center, self.inner_width,
                                     self.inner_height, self.angle,
                                     self.meta, self.visual)

    @property
    def _outer_region(self):
        return self._component_class(self.center, self.outer_width,
                                     self.outer_height, self.angle,
                                     self.meta, self.visual)

    def to_sky_args(self, wcs):
        center = wcs.pixel_to_world(self.center.x, self.center.y)
        _, pixscale, north_angle = pixel_scale_angle_at_skycoord(center, wcs)
        inner_width = (self.inner_width * u.pix * pixscale).to(u.arcsec)
        outer_width = (self.outer_width * u.pix * pixscale).to(u.arcsec)
        inner_height = (self.inner_height * u.pix * pixscale).to(u.arcsec)
        outer_height = (self.outer_height * u.pix * pixscale).to(u.arcsec)
        # region sky angles are defined relative to the WCS longitude axis;
        # photutils aperture sky angles are defined as the PA of the
        # semimajor axis (i.e., relative to the WCS latitude axis)
        angle = self.angle - (north_angle - 90 * u.deg)

        return (center, inner_width, outer_width, inner_height, outer_height,
                angle)


class AsymmetricAnnulusSkyRegion(SkyRegion):
    """
    Helper class for asymmetric annuli sky regions.

    Used for ellipse and rectangle sky annulus regions.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the annulus.
    inner_width : `~astropy.units.Quantity`
        The inner width of the annulus (before rotation) as an angle.
    outer_width : `~astropy.units.Quantity`
        The outer width of the annulus (before rotation) as an angle.
    inner_height : `~astropy.units.Quantity`
        The inner height of the annulus (before rotation) as an angle.
    outer_height : `~astropy.units.Quantity`
        The outer height of the annulus (before rotation) as an angle.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the annulus, measured anti-clockwise. If
        set to zero (the default), the width axis is lined up with the
        longitude axis of the celestial coordinates.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    _params = ('center', 'inner_width', 'outer_width', 'inner_height',
               'outer_height', 'angle')

    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    inner_width = PositiveScalarAngle('The inner width (before rotation) as '
                                      'a |Quantity| angle.')
    outer_width = PositiveScalarAngle('The outer width (before rotation) as '
                                      'a |Quantity| angle.')
    inner_height = PositiveScalarAngle('The inner height (before rotation) '
                                       'as a |Quantity| angle.')
    outer_height = PositiveScalarAngle('The outer height (before rotation) '
                                       'as a |Quantity| angle.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a'
                        '|Quantity| angle.')
    meta = RegionMetaDescr('The meta attributes as a |RegionMeta|')
    visual = RegionVisualDescr('The visual attributes as a |RegionVisual|.')

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        self.center = center
        self.inner_width = inner_width
        self.outer_width = outer_width
        self.inner_height = inner_height
        self.outer_height = outer_height
        self.angle = angle
        self.meta = meta or RegionMeta()
        self.visual = visual or RegionVisual()

        if inner_width >= outer_width:
            raise ValueError('outer_width must be greater than inner_width')
        if inner_height >= outer_height:
            raise ValueError('outer_height must be greater than inner_height')

    def to_pixel_args(self, wcs):
        center, pixscale, north_angle = pixel_scale_angle_at_skycoord(
            self.center, wcs)
        center = PixCoord(center.x, center.y)
        inner_width = (self.inner_width / pixscale).to(u.pix).value
        outer_width = (self.outer_width / pixscale).to(u.pix).value
        inner_height = (self.inner_height / pixscale).to(u.pix).value
        outer_height = (self.outer_height / pixscale).to(u.pix).value
        # region sky angles are defined relative to the WCS longitude axis;
        # photutils aperture sky angles are defined as the PA of the
        # semimajor axis (i.e., relative to the WCS latitude axis)
        angle = self.angle + (north_angle - 90 * u.deg)

        return (center, inner_width, outer_width, inner_height, outer_height,
                angle)


class EllipseAnnulusPixelRegion(AsymmetricAnnulusPixelRegion):
    """
    A elliptical annulus in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the elliptical annulus.
    inner_width : float
        The inner width of the elliptical annulus (before rotation) in
        pixels.
    outer_width : float
        The outer width of the elliptical annulus (before rotation) in
        pixels.
    inner_height : float
        The inner height of the elliptical annulus (before rotation) in
        pixels.
    outer_height : float
        The outer height of the elliptical annulus (before rotation) in
        pixels.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the elliptical annulus, measured
        anti-clockwise. If set to zero (the default), the width axis is
        lined up with the x axis.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from astropy.coordinates import Angle
        from regions import PixCoord, EllipseAnnulusPixelRegion
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        reg = EllipseAnnulusPixelRegion(PixCoord(6, 6),
                                        inner_width=5.5,
                                        outer_width=8.5,
                                        inner_height=3.5,
                                        outer_height=6.5,
                                        angle=Angle('45deg'))
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Ellipse Annulus')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(-5, 20)
        ax.set_ylim(-5, 20)
        ax.set_aspect('equal')
    """

    # duplicated from AsymmetricAnnulusPixelRegion because otherwise Sphinx
    # ignores the docstrings in the parent class
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    inner_width = PositiveScalar('The inner width (before rotation) in '
                                 'pixels as a float.')
    outer_width = PositiveScalar('The outer width (before rotation) in '
                                 'pixels as a float.')
    inner_height = PositiveScalar('The inner height (before rotation) in '
                                  'pixels as a float.')
    outer_height = PositiveScalar('The outer height (before rotation) in '
                                  'pixels as a float.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')

    _component_class = EllipsePixelRegion

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        super().__init__(center, inner_width, outer_width, inner_height,
                         outer_height, angle, meta, visual)

    def to_sky(self, wcs):
        return EllipseAnnulusSkyRegion(*self.to_sky_args(wcs),
                                       meta=self.meta.copy(),
                                       visual=self.visual.copy())


class EllipseAnnulusSkyRegion(AsymmetricAnnulusSkyRegion):
    """
    A elliptical annulus in `~astropy.coordinates.SkyCoord` coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the elliptical annulus.
    inner_width : `~astropy.units.Quantity`
        The inner width of the elliptical annulus (before rotation) as
        an angle.
    outer_width : `~astropy.units.Quantity`
        The outer width of the elliptical annulus (before rotation) as
        an angle.
    inner_height : `~astropy.units.Quantity`
        The inner height of the elliptical annulus (before rotation) as
        an angle.
    outer_height : `~astropy.units.Quantity`
        The outer height of the elliptical annulus (before rotation) as
        an angle.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the elliptical annulus, measured
        anti-clockwise. If set to zero (the default), the width axis is
        lined up with the longitude axis of the celestial coordinates.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    # duplicated from AsymmetricAnnulusSkyRegion because otherwise Sphinx
    # ignores the docstrings in the parent class
    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    inner_width = PositiveScalarAngle('The inner width (before rotation) as '
                                      'a |Quantity| angle.')
    outer_width = PositiveScalarAngle('The outer width (before rotation) as '
                                      'a |Quantity| angle.')
    inner_height = PositiveScalarAngle('The inner height (before rotation) '
                                       'as a |Quantity| angle.')
    outer_height = PositiveScalarAngle('The outer height (before rotation) '
                                       'as a |Quantity| angle.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')

    _component_class = EllipseSkyRegion

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        super().__init__(center, inner_width, outer_width, inner_height,
                         outer_height, angle, meta, visual)

    def to_pixel(self, wcs):
        return EllipseAnnulusPixelRegion(*self.to_pixel_args(wcs),
                                         meta=self.meta.copy(),
                                         visual=self.visual.copy())


class RectangleAnnulusPixelRegion(AsymmetricAnnulusPixelRegion):
    """
    A rectangular annulus in pixel coordinates.

    Parameters
    ----------
    center : `~regions.PixCoord`
        The position of the center of the rectangular annulus.
    inner_width : float
        The inner width of the rectangular annulus (before rotation) in
        pixels.
    outer_width : float
        The outer width of the rectangular annulus (before rotation) in
        pixels.
    inner_height : float
        The inner height of the rectangular annulus (before rotation) in
        pixels.
    outer_height : float
        The outer height of the rectangular annulus (before rotation) in
        pixels.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the rectangular annulus, measured
        anti-clockwise. If set to zero (the default), the width axis is
        lined up with the x axis.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.

    Examples
    --------
    .. plot::
        :include-source:

        from astropy.coordinates import Angle
        from regions import PixCoord, RectangleAnnulusPixelRegion
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1)

        reg = RectangleAnnulusPixelRegion(PixCoord(x=6, y=6),
                                          inner_width=5.5,
                                          outer_width=8.5,
                                          inner_height=3.5,
                                          outer_height=6.5,
                                          angle=Angle('45deg'))
        patch = reg.plot(ax=ax, facecolor='none', edgecolor='red', lw=2,
                         label='Rectangle Annulus')

        ax.legend(handles=(patch,), loc='upper center')
        ax.set_xlim(-5, 20)
        ax.set_ylim(-5, 20)
        ax.set_aspect('equal')
    """

    # duplicated from AsymmetricAnnulusPixelRegion because otherwise Sphinx
    # ignores the docstrings in the parent class
    center = ScalarPixCoord('The center pixel position as a |PixCoord|.')
    inner_width = PositiveScalar('The inner width (before rotation) in '
                                 'pixels as a float.')
    outer_width = PositiveScalar('The outer width (before rotation) in '
                                 'pixels as a float.')
    inner_height = PositiveScalar('The inner height (before rotation) in '
                                  'pixels as a float.')
    outer_height = PositiveScalar('The outer height (before rotation) in '
                                  'pixels as a float.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')

    _component_class = RectanglePixelRegion

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        super().__init__(center, inner_width, outer_width, inner_height,
                         outer_height, angle, meta, visual)

    def to_sky(self, wcs):
        return RectangleAnnulusSkyRegion(*self.to_sky_args(wcs),
                                         meta=self.meta.copy(),
                                         visual=self.visual.copy())


class RectangleAnnulusSkyRegion(AsymmetricAnnulusSkyRegion):
    """
    A rectangular annulus in `~astropy.coordinates.SkyCoord`
    coordinates.

    Parameters
    ----------
    center : `~astropy.coordinates.SkyCoord`
        The position of the center of the rectangular annulus.
    inner_width : `~astropy.units.Quantity`
        The inner width of the rectangular annulus (before rotation) as
        an angle.
    outer_width : `~astropy.units.Quantity`
        The outer width of the rectangular annulus (before rotation) as
        an angle.
    inner_height : `~astropy.units.Quantity`
        The inner height of the rectangular annulus (before rotation) as
        an angle.
    outer_height : `~astropy.units.Quantity`
        The outer height of the rectangular annulus (before rotation) as
        an angle.
    angle : `~astropy.units.Quantity`, optional
        The rotation angle of the rectangular annulus, measured
        anti-clockwise. If set to zero (the default), the width axis is
        lined up with the longitude axis of the celestial coordinates.
    meta : `~regions.RegionMeta` or `dict`, optional
        A dictionary that stores the meta attributes of the region.
    visual : `~regions.RegionVisual` or `dict`, optional
        A dictionary that stores the visual meta attributes of the
        region.
    """

    # duplicated from AsymmetricAnnulusSkyRegion because otherwise Sphinx
    # ignores the docstrings in the parent class
    center = ScalarSkyCoord('The center position as a |SkyCoord|.')
    inner_width = PositiveScalarAngle('The inner width (before rotation) as '
                                      'a |Quantity| angle.')
    outer_width = PositiveScalarAngle('The outer width (before rotation) as '
                                      'a |Quantity| angle.')
    inner_height = PositiveScalarAngle('The inner height (before rotation) '
                                       'as a |Quantity| angle.')
    outer_height = PositiveScalarAngle('The outer height (before rotation) '
                                       'as a |Quantity| angle.')
    angle = ScalarAngle('The rotation angle measured anti-clockwise as a '
                        '|Quantity| angle.')

    _component_class = RectangleSkyRegion

    def __init__(self, center, inner_width, outer_width, inner_height,
                 outer_height, angle=0 * u.deg, meta=None, visual=None):
        super().__init__(center, inner_width, outer_width, inner_height,
                         outer_height, angle, meta, visual)

    def to_pixel(self, wcs):
        return RectangleAnnulusPixelRegion(*self.to_pixel_args(wcs),
                                           meta=self.meta.copy(),
                                           visual=self.visual.copy())
