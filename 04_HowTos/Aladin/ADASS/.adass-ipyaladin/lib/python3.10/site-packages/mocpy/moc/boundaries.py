import networkx as nx
import numpy as np

try:  # noqa: SIM105
    from astropy_healpix import HEALPix
except ImportError:
    pass

import astropy.units as u
from astropy.coordinates import ICRS, SkyCoord


class Boundaries:  # noqa: D101
    @staticmethod
    def get(moc, order):  # noqa: D102
        boundaries_l = []

        # Get the ipixels of the MOC at the deepest order
        max_order, ipixels = Boundaries._compute_HEALPix_indices(moc, order)

        # Split the global MOC graph into all its non connected subgraphs.
        graph_boundaries = Boundaries._compute_graph_HEALPix_boundaries(
            max_order,
            ipixels,
        )
        boundaries_l.extend(Boundaries._retrieve_skycoords(graph_boundaries))

        return boundaries_l

    @staticmethod
    def _compute_HEALPix_indices(m, order):
        moc = m
        if order and m.max_order > order:
            moc = m.degrade_to_order(order)

        max_order = moc.max_order
        ipixels = moc.flatten()

        # Take the complement if the MOC covers more than half of the sky => the perimeter(MOC) = perimeter(complement(MOC))
        # but we process less hpx cells
        num_ipixels = 3 << (2 * (max_order + 1))
        ipixels.shape[0] / float(num_ipixels)

        # if sky_fraction > 0.5:
        #    ipixels_all = np.arange(num_ipixels)
        #    ipixels = np.setdiff1d(ipixels_all, ipixels, assume_unique=True)

        return max_order, ipixels

    @staticmethod
    def _compute_graph_HEALPix_boundaries(depth, ipixels):
        def insert_edge(G, l1, l2, p1_lon, p1_lat, p2_lon, p2_lat):
            # Nodes are indexed by str(skycoord). When getting ordered nodes, one can retrieve back the skycoord instance
            # by accessing the python dict `pts_d`.
            try:
                # Avoid the special case where holes are touching to each other
                # 'x' belongs to the MOC
                # ' ' is part of the holes in the MOC
                #    |xxx
                #    |xxx
                # ---A---
                # xxx|
                # xxx|
                #
                # If this case occurs we split the node A into 2. One is attached to the bottom left graph and the other to the
                # top right one. When computing the MST (minimal spanning tree) from a graph, we need our graphs to have
                # only nodes of degrees 1 or 2 (i.e. to be lines).
                if G.degree[l1] >= 2:
                    l1 += "_"
            # TODO: Fix bare except
            except:  # noqa: S110, E722
                pass

            try:
                if G.degree[l2] >= 2:
                    l2 += "_"
            # TODO: fix bare except
            except:  # noqa: S110, E722
                pass
            # Set the skycoord instance as an attribute of the nodes
            G.add_node(l1, ra=p1_lon, dec=p1_lat)
            G.add_node(l2, ra=p2_lon, dec=p2_lat)
            G.add_edge(l1, l2)

        # Phase 1: Retrieve the ipixels located at the border of
        # this connexe MOC component
        ipixels = ipixels.astype(int)
        hp = HEALPix(nside=(1 << depth), order="nested", frame=ICRS())
        neighbours = hp.neighbours(ipixels)[[0, 2, 4, 6], :]
        ipixels = ipixels.astype(np.uint64)

        neighbours = neighbours.astype(np.uint64)

        isin = np.isin(neighbours, ipixels)
        border = isin.sum(axis=0) < 4

        ipixels_border = ipixels[border]
        isin_border = isin[:, border]

        # Phase 2: Build the graph from the positions of the ipixels boundaries
        # import cdshealpix as healpix
        # ipix_lon, ipix_lat = healpix.vertices(ipixels_border, depth)
        # ipix_lon[:, [1, 3]] = ipix_lon[:, [3, 1]]
        # ipix_lat[:, [1, 3]] = ipix_lat[:, [3, 1]]
        ipix_lon, ipix_lat = hp.boundaries_lonlat(ipixels_border, step=1)

        ipix_lon_deg = ipix_lon.to_value(u.deg)
        ipix_lat_deg = ipix_lat.to_value(u.deg)

        ipix_lon_repr = np.around(
            np.asarray(ipix_lon.reshape((1, -1))[0]),
            decimals=3,
        ).tolist()
        ipix_lat_repr = np.around(
            np.asarray(ipix_lat.reshape((1, -1))[0]),
            decimals=3,
        ).tolist()

        west_border = ~isin_border[0, :]
        south_border = ~isin_border[1, :]
        east_border = ~isin_border[2, :]
        north_border = ~isin_border[3, :]

        E = nx.Graph()

        for i in range(ipixels_border.shape[0]):
            lon_deg = ipix_lon_deg[i]
            lat_deg = ipix_lat_deg[i]

            p0_lon = lon_deg[0]
            p1_lon = lon_deg[1]
            p2_lon = lon_deg[2]
            p3_lon = lon_deg[3]

            p0_lat = lat_deg[0]
            p1_lat = lat_deg[1]
            p2_lat = lat_deg[2]
            p3_lat = lat_deg[3]

            off = 4 * i
            off2 = 4 * (i + 1)
            repr_lon = ipix_lon_repr[off:off2]
            repr_lat = ipix_lat_repr[off:off2]

            s0 = str(repr_lon[0]) + " " + str(repr_lat[0])
            s1 = str(repr_lon[1]) + " " + str(repr_lat[1])
            s2 = str(repr_lon[2]) + " " + str(repr_lat[2])
            s3 = str(repr_lon[3]) + " " + str(repr_lat[3])

            # WEST border
            if west_border[i]:
                insert_edge(E, s1, s2, p1_lon, p1_lat, p2_lon, p2_lat)

            # NORTH border
            if north_border[i]:
                insert_edge(E, s2, s3, p2_lon, p2_lat, p3_lon, p3_lat)

            # EAST border
            if east_border[i]:
                insert_edge(E, s3, s0, p3_lon, p3_lat, p0_lon, p0_lat)

            # SOUTH border
            if south_border[i]:
                insert_edge(E, s0, s1, p0_lon, p0_lat, p1_lon, p1_lat)

        return E

    @staticmethod
    def _retrieve_skycoords(V):
        coords_l = []
        # Accessing the borders one by one. At this step, V_subgraphs contains a list of cycles
        # (i.e. one describing the external border of the MOC component and several describing the holes
        # found in the MOC component).
        # See https://stackoverflow.com/questions/61154740/attributeerror-module-networkx-has-no-attribute-connected-component-subgraph
        #     https://networkx.org/documentation/networkx-2.1/reference/algorithms/generated/networkx.algorithms.components.connected_component_subgraphs.html
        # V_subgraphs = nx.connected_component_subgraphs(V)
        V_subgraphs = (V.subgraph(c) for c in nx.connected_components(V))
        for v in V_subgraphs:
            # Compute the MST for each cycle
            v = nx.convert_node_labels_to_integers(v)
            mst = nx.minimum_spanning_tree(v)
            # Get one end of the span tree by looping over its node and checking if the degree is one
            src = None
            for node, deg in mst.degree():
                if deg == 1:
                    src = node
                    break

            # Get the unordered lon and lat
            ra = np.asarray(list(nx.get_node_attributes(v, "ra").values()))
            dec = np.asarray(list(nx.get_node_attributes(v, "dec").values()))
            coords = np.vstack((ra, dec)).T
            # Get the ordering from the MST
            ordering = np.asarray(list(nx.dfs_preorder_nodes(mst, src)))
            # Order the coords
            coords = coords[ordering]
            # Get a skycoord containing N coordinates computed from the Nx2 `coords` array
            coords = SkyCoord(coords, unit="deg")
            coords_l.append(coords)

        return coords_l
