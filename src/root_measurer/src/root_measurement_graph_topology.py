"""Distinguish adjacency cycles from H1 of closed cubical skeleton voxels.

The full 26-neighbor graph is the intersection graph of unit cubes. A triple
of these axis-aligned cubes has a common intersection exactly when each pair
intersects. Its triangle is filled, not a material tunnel. H1 uses the rank of
these triangle boundary relations over GF(2); no coordinates or edges are
deleted. Higher intersection simplices do not change the image of d2.

This measures the supplied voxel complex. It does not certify thinning,
anatomical identity, observed boundary confidence, or physical centerlines.
"""
import numpy as np
from root_reconstruction_measurement import _indices, neighbor_graph


def skeleton_voxel_h1(cells, max_nodes=100000, max_edges=1000000,
                      max_triangles=2000000):
    cells = _indices(cells)
    if len(cells) > max_nodes:
        raise ValueError('TOPOLOGY_NODE_GUARD')
    pairs, labels = neighbor_graph(cells)
    if len(pairs) > max_edges:
        raise ValueError('TOPOLOGY_EDGE_GUARD')
    adjacency = [set() for _ in cells]
    edge_id = {}
    for index, (a, b) in enumerate(pairs):
        a, b = int(a), int(b)
        adjacency[a].add(b); adjacency[b].add(a)
        edge_id[(a,b)] = index
    pivots = {}
    triangles = 0
    for a, neighbors in enumerate(adjacency):
        for b in sorted(n for n in neighbors if n > a):
            for c in sorted(n for n in neighbors.intersection(adjacency[b]) if n > b):
                triangles += 1
                if triangles > max_triangles:
                    raise ValueError('TOPOLOGY_TRIANGLE_GUARD')
                # Pairwise closed intervals on each coordinate intersect jointly.
                if np.any(np.ptp(cells[[a,b,c]], axis=0) > 1):
                    raise AssertionError('TRIANGLE_WITHOUT_CUBE_INTERSECTION')
                column = (1 << edge_id[(a,b)]) ^ (1 << edge_id[(a,c)]) ^ (1 << edge_id[(b,c)])
                while column:
                    pivot = column.bit_length()-1
                    if pivot in pivots:
                        column ^= pivots[pivot]
                    else:
                        pivots[pivot] = column
                        break
    components = len(np.unique(labels))
    graph_rank = len(pairs)-len(cells)+components
    h1 = graph_rank-len(pivots)
    if h1 < 0:
        raise AssertionError('NEGATIVE_FIRST_HOMOLOGY')
    return dict(status='COMPLETE_CLOSED_SKELETON_VOXEL_H1', nodes=len(cells), edges=len(pairs),
                foreground_components=components, adjacency_cycle_rank=graph_rank,
                filled_intersection_triangles=triangles, triangle_boundary_rank_GF2=len(pivots),
                skeleton_voxel_H1=h1, H2=None,
                coordinates_unchanged=True, edges_deleted=0, occupied_material_changed=False,
                biological_centerline_qualified=False, thinning_topology_qualified=False,
                definition='H1 of closed unit cubes at supplied skeleton indices; filled triangle boundary relations over GF(2)')
