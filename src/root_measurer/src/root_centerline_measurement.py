"""Topology-aware measurement graphs and a separate display-curve ledger.

Only isolated, filled equal-diagonal triangle faces are rewritten as a star
at their common cube corner. This preserves the local contractible face and
does not connect components or add occupied material. Ambiguous complexes
remain unsupported. The resulting graph is a geometric representation, not
an anatomical centerline or a physical-tip qualification.
"""
from collections import Counter
import math
import numpy as np
from root_reconstruction_measurement import _indices
from root_measurement_graph_topology import skeleton_voxel_h1
from root_volume_topology import cubical_euler
from root_skeleton_topology import branches_from_adjacency
from root_polyline_geometry import fit_polylines


def _validate_graph(points, adjacency):
    points = np.asarray(points, float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError('FINITE_GRAPH_POINTS_REQUIRED')
    if len(points) != len(adjacency) or len(np.unique(points, axis=0)) != len(points):
        raise ValueError('GRAPH_NODE_IDENTITY')
    out = [set(map(int, ns)) for ns in adjacency]
    for i, ns in enumerate(out):
        if i in ns or any(j < 0 or j >= len(out) for j in ns):
            raise ValueError('GRAPH_SELF_OR_INVALID_EDGE')
        if any(i not in out[j] for j in ns):
            raise ValueError('GRAPH_ASYMMETRIC_EDGE')
    return points, out


def _components(adjacency):
    labels = np.full(len(adjacency), -1, int)
    count = 0
    for start in range(len(adjacency)):
        if labels[start] >= 0:
            continue
        labels[start] = count; todo = [start]
        while todo:
            i = todo.pop()
            for j in adjacency[i]:
                if labels[j] < 0:
                    labels[j] = count; todo.append(j)
        count += 1
    return count, labels


def canonicalize_corner_junctions(cells, adjacency, pitch):
    """Represent isolated corner-intersection triangular faces by stars.

    Existing node coordinates and all non-triangle edges are retained.
    New corner nodes have no invented radius and keep contact ambiguity.
    Unsupported faces/cavities cannot be silently promoted to centerlines.
    """
    cells = _indices(cells)
    if not math.isfinite(pitch) or pitch <= 0:
        raise ValueError('FINITE_POSITIVE_PITCH_REQUIRED')
    points, out = _validate_graph((cells+.5)*pitch, adjacency)
    for i, ns in enumerate(out):
        for j in ns:
            if np.max(np.abs(cells[i]-cells[j])) > 1:
                raise ValueError('NON_NEIGHBOR_VOXEL_EDGE')
    count_before, labels_before = _components(out)
    triangles = [(i,j,k) for i, ns in enumerate(out) for j in sorted(ns) if j>i
                 for k in sorted(ns.intersection(out[j])) if k>j]
    multiplicity = Counter(tuple(sorted(e)) for tri in triangles for e in
                           ((tri[0],tri[1]),(tri[0],tri[2]),(tri[1],tri[2])))
    candidates, unsupported = [], []
    for tri in triangles:
        edges = [(tri[0],tri[1]),(tri[0],tri[2]),(tri[1],tri[2])]
        equal_diagonal = all(int(np.sum((cells[i]-cells[j])**2)) == 2 for i,j in edges)
        lo = cells[list(tri)].max(axis=0)
        hi = (cells[list(tri)]+1).min(axis=0)
        if equal_diagonal and np.array_equal(lo,hi) and all(multiplicity[e] == 1 for e in edges):
            candidates.append((tuple(lo), tri, edges))
        else:
            unsupported.append(dict(nodes=list(tri), reason='NONISOLATED_OR_UNSUPPORTED_INTERSECTION_FACE'))
    corner_counts = Counter(corner for corner,_,_ in candidates)
    corrected = []
    positions = points.tolist()
    for corner, tri, edges in sorted(candidates):
        if corner_counts[corner] != 1:
            unsupported.append(dict(nodes=list(tri), reason='AMBIGUOUS_MULTIPLE_FACES_AT_CORNER'))
            continue
        center = np.array(corner, float)*pitch
        if any(np.array_equal(center, p) for p in positions):
            unsupported.append(dict(nodes=list(tri), reason='EXISTING_NODE_AT_CORNER'))
            continue
        for i in tri:
            if not np.all((center >= cells[i]*pitch-1e-14) & (center <= (cells[i]+1)*pitch+1e-14)):
                raise AssertionError('NEW_CONTACT_OUTSIDE_SOURCE_CUBE')
        node = len(positions); positions.append(center.tolist()); out.append(set(tri))
        for i,j in edges:
            out[i].remove(j); out[j].remove(i)
        for i in tri:
            out[i].add(node)
        corrected.append(dict(source_triangle=list(tri), new_contact_node=node,
                              coordinates=center.tolist(), new_radius=None,
                              contact_ambiguity_preserved=True,
                              spokes_contained_in_original_closed_cubes=True))
    positions, out = _validate_graph(np.array(positions, float).reshape(-1,3), out)
    count_after, labels_after = _components(out)
    if count_before != count_after:
        raise AssertionError('CORNER_REPRESENTATION_COMPONENT_CHANGE')
    # Compare partitions, not arbitrary component label numbers.
    for component in range(count_before):
        if len(np.unique(labels_after[:len(cells)][labels_before == component])) != 1:
            raise AssertionError('SOURCE_COMPONENT_PARTITION_CHANGED')
    edge_count = sum(map(len,out))//2
    graph_rank = edge_count-len(out)+count_after
    homology = skeleton_voxel_h1(cells)
    voxel_euler = cubical_euler(cells)['euler']
    voxel_h2 = voxel_euler-homology['foreground_components']+homology['skeleton_voxel_H1']
    passed = (count_after == homology['foreground_components'] and
              graph_rank == homology['skeleton_voxel_H1'] and voxel_h2 == 0 and not unsupported)
    audit = dict(status='PASS_SUPPORTED_GEOMETRIC_GRAPH_TOPOLOGY' if passed else
                 'UNRESOLVED_VOXEL_COMPLEX_NOT_ADMITTED_AS_CENTERLINE_GRAPH',
                 original_nodes=len(cells), added_contact_nodes=len(corrected),
                 corrected_faces=corrected, unsupported_faces=unsupported,
                 source_component_count=count_before, graph_component_count=count_after,
                 graph_cycle_rank=graph_rank, skeleton_voxel_H1=homology['skeleton_voxel_H1'],
                 skeleton_voxel_H2=voxel_h2, original_positions_unchanged=True,
                 occupied_material_added=0, unknown_cells_filled=0, components_joined=0,
                 genuine_cycles_forced_to_tree=False, anatomical_connectivity_inferred=False,
                 physical_centerline_qualified=False)
    return positions, out, audit


def measurement_ledger(points, adjacency, radii, pitch):
    """Keep measurement-graph arclength distinct from display fairing length."""
    points, adjacency = _validate_graph(points, adjacency)
    if not math.isfinite(pitch) or pitch <= 0:
        raise ValueError('FINITE_POSITIVE_PITCH_REQUIRED')
    if len(radii) != len(points):
        raise ValueError('RADIUS_NODE_COUNT')
    radii = np.asarray([np.nan if r is None else float(r) for r in radii])
    if np.isinf(radii).any() or np.any(radii[np.isfinite(radii)] < 0):
        raise ValueError('INVALID_RADIUS')
    radii[radii == 0] = np.nan
    branches, audit = branches_from_adjacency(points, adjacency, radii)
    for b in branches:
        b['radii'] = [float(r) if np.isfinite(r) else None for r in b['radii']]
    edges = [(i,j) for i, ns in enumerate(adjacency) for j in ns if i<j]
    length = math.fsum(float(np.linalg.norm(points[i]-points[j])) for i,j in edges)
    branch_sum = math.fsum(b['arc_length'] for b in branches)
    if not math.isclose(length, branch_sum, rel_tol=1e-12, abs_tol=1e-12):
        raise AssertionError('BRANCH_LENGTH_NOT_CONSERVED')
    degree = np.array([len(ns) for ns in adjacency])
    isolated = np.flatnonzero(degree == 0)
    decoded = dict(status='OK', branches=branches, audit=dict(pitch=pitch),
                   total_length=length, component_count=audit['components'],
                   endpoints=points[degree==1].tolist(), branchpoints=points[degree>=3].tolist())
    display = None; display_error = None; shift_bound = None
    try:
        if branches:
            display = fit_polylines(decoded)
            key = lambda p: tuple(np.round(p,12))
            mapping = {key(p):p for p in points}
            for b in display['branches']:
                for old,new in zip(b['native_nodes'],b['nodes']):
                    mapping[key(old)] = np.asarray(new)
            displacement = np.array([np.linalg.norm(mapping[key(p)]-p) for p in points])
            shift_bound = math.fsum(float((displacement[i]+displacement[j])) for i,j in edges)
            if abs(display['total_length']-length) > shift_bound+1e-11:
                raise AssertionError('DISPLAY_LENGTH_CHANGE_EXCEEDS_PROVED_BOUND')
    except ValueError as error:
        display_error = str(error)
    return dict(status='COMPLETE_SEPARATE_GRAPH_AND_DISPLAY_METRIC_LEDGER',
                measurement_graph_length=length, serialized_branch_length=branch_sum,
                display_curve_length=None if display is None else display['total_length'],
                display_total_length_displacement_bound=shift_bound,
                display_error=display_error,
                display_never_overwrites_measurement_metric=True,
                biological_total_root_length=None, physical_endpoint_length=None,
                graph_components=audit['components'], physical_edges=len(edges),
                all_edges_once=True, isolated_nodes=isolated.tolist(),
                isolated_coordinates=points[isolated].tolist(),
                endpoints=decoded['endpoints'], branchpoints=decoded['branchpoints'],
                branches=branches, display_branches=None if display is None else display['branches'],
                radius_definition='supplied geometric radius at original node; no positive radius fabricated at contact nodes',
                units='normalized_model_coordinates', scientific_accuracy_qualified=False)
