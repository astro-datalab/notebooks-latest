import cdshealpix
import numpy as np
from astropy.coordinates import ICRS, SkyCoord
from astropy.wcs.utils import skycoord_to_pixel


def backface_culling(xp, yp):
    """Remove cells crossing the MOC after projection.

    The remaining HEALPix cells are used for computing the patch of the MOC
    """
    vx = xp
    vy = yp

    def compute_vector_at_index(X, Y, i):
        cur = i
        next_index = (i + 1) % 4

        x = X[:, next_index] - X[:, cur]
        y = Y[:, next_index] - Y[:, cur]
        z = np.zeros(x.shape)

        return np.vstack((x, y, z)).T

    # A ----- B
    #  \      |
    #   D-----C
    # Compute the cross product between AB and BC
    # and the cross product between BC and CD
    AB = compute_vector_at_index(vx, vy, 0)
    BC = compute_vector_at_index(vx, vy, 1)
    CD = compute_vector_at_index(vx, vy, 2)
    DA = compute_vector_at_index(vx, vy, 3)

    ABC = np.cross(AB, BC)
    BCD = np.cross(BC, CD)
    CDA = np.cross(CD, DA)
    DAB = np.cross(DA, AB)

    tol = -0.1
    frontface_cells = (
        (ABC[:, 2] <= tol)
        & (CDA[:, 2] <= tol)
        & (BCD[:, 2] <= tol)
        & (DAB[:, 2] <= tol)
    )

    frontface_cells = np.asarray(frontface_cells)
    vx = vx[frontface_cells]
    vy = vy[frontface_cells]

    return vx, vy, frontface_cells


def from_moc(depth_ipix_d, wcs):
    """Create a new MOC that do not contain the HEALPix cells that are backfacing the projection."""
    depths = [int(depth_str) for depth_str in depth_ipix_d]
    min_depth = min(depths)
    max_depth = max(depths)
    ipixels = np.asarray(depth_ipix_d[str(min_depth)])

    # Split the cells located at the border of the projection
    # until at least the depth 7
    max_split_depth = max(7, max_depth)

    ipix_d = {}
    for depth in range(min_depth, max_split_depth + 1):
        if depth < 3:
            too_large_ipix = ipixels
        else:
            ipix_lon, ipix_lat = cdshealpix.vertices(ipixels, depth)

            ipix_lon = ipix_lon[:, [2, 3, 0, 1]]
            ipix_lat = ipix_lat[:, [2, 3, 0, 1]]
            ipix_vertices = SkyCoord(ipix_lon, ipix_lat, frame=ICRS())

            # Projection on the given WCS
            xp, yp = skycoord_to_pixel(coords=ipix_vertices, wcs=wcs)
            _, _, frontface_id = backface_culling(xp, yp)

            # Get the pixels which are backfacing the projection
            backfacing_ipix = ipixels[~frontface_id]  # backfacing
            frontface_ipix = ipixels[frontface_id]

            depth_str = str(depth)
            ipix_d.update({depth_str: frontface_ipix})

            too_large_ipix = backfacing_ipix

        next_depth = str(depth + 1)
        ipixels = []

        if next_depth in depth_ipix_d:
            ipixels = depth_ipix_d[next_depth]

        for ipix in too_large_ipix:
            child_ipix = ipix << 2
            ipixels.extend([child_ipix, child_ipix + 1, child_ipix + 2, child_ipix + 3])

        ipixels = np.asarray(ipixels)

    return ipix_d
