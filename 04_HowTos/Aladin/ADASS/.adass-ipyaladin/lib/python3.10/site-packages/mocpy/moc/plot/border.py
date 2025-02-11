import cdshealpix
import numpy as np
from astropy.coordinates import ICRS, SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from matplotlib.patches import PathPatch
from matplotlib.path import Path


def border(moc, ax, wcs, **kw_mpl_pathpatch):
    """Highlight the borders of a MOC in a plot."""
    from .utils import _set_wcs, build_plotting_moc

    moc_to_plot = build_plotting_moc(moc, wcs)

    if moc_to_plot.empty():
        return

    max_order = moc_to_plot.max_order
    ipixels_open = moc_to_plot.flatten()
    num_ipixels = 3 << (2 * (max_order + 1))

    # Take the complement if the MOC covers more than half of the sky
    sky_fraction = moc_to_plot.sky_fraction

    if sky_fraction > 0.5:
        ipixels_all = np.arange(num_ipixels)
        ipixels_open = np.setdiff1d(ipixels_all, ipixels_open, assume_unique=True)

    neighbors = cdshealpix.neighbours(ipixels_open, max_order).T
    # Select the direct neighbors (i.e. those in WEST, NORTH, EAST and SOUTH directions)
    neighbors = neighbors[[3, 7, 5, 1], :]
    ipix_moc = np.isin(neighbors, ipixels_open)

    west_edge = ipix_moc[0, :]
    south_edge = ipix_moc[1, :]
    east_edge = ipix_moc[2, :]
    north_edge = ipix_moc[3, :]

    num_ipix_moc = ipix_moc.sum(axis=0)

    ipixels_border_id = num_ipix_moc < 4
    # The border of each HEALPix cells is drawn one at a time
    path_vertices_l = []
    codes = []

    west_border = west_edge[ipixels_border_id]
    south_border = south_edge[ipixels_border_id]
    east_border = east_edge[ipixels_border_id]
    north_border = north_edge[ipixels_border_id]
    ipixels_border = ipixels_open[ipixels_border_id]

    ipix_lon_boundaries, ipix_lat_boundaries = cdshealpix.vertices(
        ipixels_border,
        max_order,
    )
    ipix_boundaries = SkyCoord(ipix_lon_boundaries, ipix_lat_boundaries, frame=ICRS())

    # Projection on the given WCS
    xp, yp = skycoord_to_pixel(coords=ipix_boundaries, wcs=wcs)
    from . import culling_backfacing_cells

    xp, yp, frontface_id = culling_backfacing_cells.backface_culling(xp, yp)

    west_border = west_border[frontface_id]
    south_border = south_border[frontface_id]
    east_border = east_border[frontface_id]
    north_border = north_border[frontface_id]

    for i in range(xp.shape[0]):
        vx = xp[i]
        vy = yp[i]
        if not north_border[i]:
            path_vertices_l += [(vx[0], vy[0]), (vx[1], vy[1]), (0, 0)]
            codes += [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]

        if not east_border[i]:
            path_vertices_l += [(vx[1], vy[1]), (vx[2], vy[2]), (0, 0)]
            codes += [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]

        if not south_border[i]:
            path_vertices_l += [(vx[2], vy[2]), (vx[3], vy[3]), (0, 0)]
            codes += [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]

        if not west_border[i]:
            path_vertices_l += [(vx[3], vy[3]), (vx[0], vy[0]), (0, 0)]
            codes += [Path.MOVETO, Path.LINETO, Path.CLOSEPOLY]

    if len(path_vertices_l) > 0:
        path = Path(path_vertices_l, codes)
        perimeter_patch = PathPatch(path, **kw_mpl_pathpatch)

        ax.add_patch(perimeter_patch)
        _set_wcs(ax, wcs)
