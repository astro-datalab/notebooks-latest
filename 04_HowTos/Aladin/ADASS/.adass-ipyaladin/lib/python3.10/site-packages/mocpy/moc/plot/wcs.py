from collections.abc import Sequence

import astropy.units as u
import numpy as np
from astropy import coordinates, wcs


class WCS:
    """
    Create a WCS for visualizing a MOC in a matplotlib axis.

    This method is a quick astropy WCS instantiation and does not fit every use case.
    It is especially not designed for all-sky plots.

    Parameters
    ----------
    fig : `~matplotlib.pyplot.figure`
        The matplotlib figure used for plotting the MOC.
    fov : `~astropy.units.Quantity` or Sequence[`~astropy.units.Quantity`, `~astropy.units.Quantity`]
        Size of the field of view.
        If it is a sequence, it must be of length 2 and represent the longitudinal and latitudinal
        field of views.
    center : `~astropy.coordinates.SkyCoord`, optional
        World coordinates matching with the center of the plot. Default to (0 deg, 0 deg) (in ICRS frame).
    coordsys : str, optional
        Coordinate system. Default to "icrs". Must be in ["icrs", "galactic"].
    projection : str, optional
        World base -> Image base projection type. See http://docs.astropy.org/en/stable/wcs/#supported-projections for
        the projections currently supported in astropy. Default to Aitoff.
    rotation : `~astropy.coordinates.Angle`, optional
        The angle of rotation. Default to no rotation.

    Returns
    -------
    wcs : `~astropy.wcs.WCS`
        The WCS that can be passed to mocpy.MOC.fill/border.

    Examples
    --------
    >>> from mocpy import MOC, WCS
    >>> from astropy.coordinates import Angle, SkyCoord
    >>> import matplotlib.pyplot as plt
    >>> import astropy.units as u
    >>> # Create a MOC
    >>> moc = MOC.from_elliptical_cone(lon = 0 * u.rad,
    ...                                lat = 0 * u.rad,
    ...                                a = 50 * u.deg,
    ...                                b = 20 * u.deg,
    ...                                pa = 90 * u.deg,
    ...                                max_depth = 13,
    ...                                )
    >>> fig = plt.figure(figsize=(10, 10))
    >>> with WCS(fig,
    ...           fov=200 * u.deg,
    ...           center=SkyCoord(0, 20, unit='deg', frame='icrs'),
    ...           coordsys="icrs",
    ...           rotation=Angle(0, u.degree),
    ...           projection="AIT") as wcs: wcs = wcs
    >>> ax = fig.add_subplot(projection=wcs)
    >>> # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    >>> moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="green")
    >>> moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    """

    def __init__(
        self,
        fig,
        fov,
        center=coordinates.SkyCoord(0, 0, unit="deg", frame="icrs"),
        coordsys="icrs",
        projection="AIT",
        rotation=coordinates.Angle(0, u.radian),
    ):
        self.w = wcs.WCS(naxis=2)

        width_px, height_px = fig.get_size_inches() * float(fig.dpi)

        if isinstance(fov, Sequence):
            cdelt_x = fov[0].to_value("deg") / float(width_px)
            cdelt_y = fov[1].to_value("deg") / float(height_px)
            self.w.wcs.cdelt = [-cdelt_x, cdelt_y]
        else:
            cdelt_x = fov.to_value("deg") / float(width_px)
            cdelt_y = fov.to_value("deg") / float(height_px)
            cdelt = max(cdelt_x, cdelt_y)
            self.w.wcs.cdelt = [-cdelt, cdelt]

        self.w.wcs.crpix = [width_px / 2.0 + 0.5, height_px / 2.0 + 0.5]

        if coordsys == "icrs":
            self.w.wcs.crval = [center.icrs.ra.deg, center.icrs.dec.deg]
            self.w.wcs.ctype = ["RA---" + projection, "DEC--" + projection]
        elif coordsys == "galactic":
            self.w.wcs.crval = [center.galactic.l.deg, center.galactic.b.deg]
            self.w.wcs.ctype = ["GLON-" + projection, "GLAT-" + projection]

        theta = rotation.radian
        self.w.wcs.pc = [
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ]

    def __enter__(self):
        return self.w

    def __exit__(self, exception_type, exception_value, traceback):
        pass
