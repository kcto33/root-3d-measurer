"""Indexed surface topology only; never repair geometry or infer a material union."""
from collections import Counter
import numpy as np


def _components(vertices, edges):
    parent = {int(v): int(v) for v in vertices}
    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for a, b in edges:
        a, b = find(int(a)), find(int(b))
        if a != b:
            parent[b] = a
    groups = {}
    for v in parent:
        groups.setdefault(find(v), []).append(v)
    return list(groups.values())


def indexed_surface_topology(mesh):
    """Certify combinatorial 2-manifold links, independently of geometric embedding.

    Betti numbers, when available, are over GF(2) for the indexed surface.
    They are neither solid Betti numbers nor biological branching measurements.
    """
    f = np.asarray(mesh.faces)
    if f.ndim != 2 or f.shape[1] != 3 or not np.issubdtype(f.dtype, np.integer):
        raise ValueError('INTEGER_TRIANGLE_FACES_REQUIRED')
    if len(f) == 0 or f.min() < 0 or f.max() >= len(mesh.vertices):
        raise ValueError('NONEMPTY_VALID_INDEXED_SURFACE_REQUIRED')
    repeated = int(np.count_nonzero((f[:, 0] == f[:, 1]) |
                                   (f[:, 0] == f[:, 2]) | (f[:, 1] == f[:, 2])))
    duplicate = len(f) - len(np.unique(np.sort(f, axis=1), axis=0))
    raw_edges = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    edges, counts = np.unique(np.sort(raw_edges, axis=1), axis=0, return_counts=True)
    used = np.unique(f)
    groups = _components(used, edges)
    boundary = edges[counts == 1]
    boundary_degree = Counter(int(x) for x in boundary.ravel())
    boundary_groups = _components(list(boundary_degree), boundary)
    boundary_types = Counter()
    for group in boundary_groups:
        degrees = [boundary_degree[v] for v in group]
        if all(d == 2 for d in degrees):
            boundary_types['cycle'] += 1
        elif degrees.count(1) == 2 and all(d in (1, 2) for d in degrees):
            boundary_types['path'] += 1
        else:
            boundary_types['branched_or_singular'] += 1
    center = f.ravel()
    opposite = np.stack([f[:, [1, 2]], f[:, [2, 0]], f[:, [0, 1]]], axis=1).reshape(-1, 2)
    order = np.argsort(center, kind='stable')
    sorted_centers = center[order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_centers)) + 1, len(center)]
    link_types = Counter()
    bad_vertices = []
    for left, right in zip(starts[:-1], starts[1:]):
        vertex = int(sorted_centers[left])
        pairs = opposite[order[left:right]]
        adjacency = {}
        degree = Counter()
        for a, b in pairs:
            a, b = int(a), int(b)
            adjacency.setdefault(a, set()).add(b)
            adjacency.setdefault(b, set()).add(a)
            degree[a] += 1
            degree[b] += 1
        seen = set()
        components = 0
        for seed in adjacency:
            if seed in seen:
                continue
            components += 1
            stack = [seed]
            seen.add(seed)
            while stack:
                for neighbor in adjacency[stack.pop()]:
                    if neighbor not in seen:
                        seen.add(neighbor)
                        stack.append(neighbor)
        simple = (len(pairs) == len({tuple(sorted(map(int, p))) for p in pairs})
                  and all(a != b for a, b in pairs))
        degrees = list(degree.values())
        if simple and components == 1 and all(d == 2 for d in degrees) and vertex not in boundary_degree:
            kind = 'interior_cycle'
        elif (simple and components == 1 and degrees.count(1) == 2
              and all(d in (1, 2) for d in degrees) and boundary_degree.get(vertex) == 2):
            kind = 'boundary_path'
        else:
            kind = 'invalid_link'
            bad_vertices.append(vertex)
        link_types[kind] += 1
    qualified = not (repeated or duplicate or bad_vertices or np.any(counts > 2))
    chi = int(len(used) - len(edges) + len(f))
    closed_components = sum(not any(v in boundary_degree for v in g) for g in groups)
    betti = [len(groups), int(len(groups) + closed_components - chi), closed_components] if qualified else None
    if qualified:
        assert betti[1] >= 0 and boundary_types.get('path', 0) == 0 and boundary_types.get('branched_or_singular', 0) == 0
    return dict(status='INDEXED_SURFACE_TOPOLOGY_DIAGNOSTICS',
                referenced_vertices=int(len(used)), unreferenced_vertices=int(len(mesh.vertices)-len(used)),
                edges=int(len(edges)), faces=int(len(f)), connected_components=len(groups),
                euler=chi, repeated_index_faces=repeated, duplicate_faces=int(duplicate),
                boundary_edges=int(len(boundary)), nonmanifold_edges=int(np.count_nonzero(counts > 2)),
                boundary_components=len(boundary_groups), boundary_component_types=dict(boundary_types),
                vertex_link_types=dict(link_types), invalid_link_vertex_indices=bad_vertices,
                indexed_two_manifold_with_boundary=qualified,
                closed_indexed_two_manifold=bool(qualified and len(boundary) == 0),
                surface_betti_GF2=betti, geometry_repaired=False,
                embedding_self_intersections_checked=False, material_union_topology=None,
                physical_volume_qualified=False, biological_topology_qualified=False,
                limitation='Indexed surface complex only; coincident indices are not welded; no material or biological topology inferred')
