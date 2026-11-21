"""Local interpretation of bridges; every cyclic edge stays in the source ledger.

Bridge orientation relative to an explicit root is unique even when the original
component has cycles. No path through a cyclic core is chosen or invented.
"""
from collections import deque
import math
import numpy as np


def graph_bridges(adjacency):
    """Iterative Tarjan low-link test, including disconnected and deep graphs."""
    discovery = {}; low = {}; parent = {}; bridges = set(); clock = 0
    for start in range(len(adjacency)):
        if start in discovery:
            continue
        discovery[start] = low[start] = clock; clock += 1; parent[start] = None
        stack = [(start, iter(adjacency[start]))]
        while stack:
            node, neighbors = stack[-1]
            other = next(neighbors, None)
            if other is None:
                stack.pop(); up = parent[node]
                if up is not None:
                    if low[node] > discovery[up]:
                        bridges.add(tuple(sorted((up, node))))
                    low[up] = min(low[up], low[node])
            elif other == parent[node]:
                continue
            elif other in discovery:
                low[node] = min(low[node], discovery[other])
            else:
                parent[other] = node; discovery[other] = low[other] = clock; clock += 1
                stack.append((other, iter(adjacency[other])))
    return bridges


def edge_regions(edges):
    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(a, []).append(b); adjacency.setdefault(b, []).append(a)
    seen = set(); regions = []
    for start in sorted(adjacency):
        if start in seen:
            continue
        queue = [start]; seen.add(start); nodes = []
        while queue:
            node = queue.pop(); nodes.append(node)
            for other in adjacency[node]:
                if other not in seen:
                    seen.add(other); queue.append(other)
        nodes.sort()
        owned = sorted((a, b) for a in nodes for b in adjacency[a] if a < b)
        regions.append((nodes, owned))
    return regions


def build_with_cycle_regions(graph, *, annotations, scale_mm_per_model_unit,
                             radius_samples, inference_options, tree_builder):
    from root_graph_measurement_ledger import graph_measurement_ledger
    from root_axis_architecture import _annotations, _axis_id
    if not isinstance(graph, dict) or 'points' not in graph or 'adjacency' not in graph:
        return tree_builder(graph, annotations=annotations, scale_mm_per_model_unit=scale_mm_per_model_unit,
                            radius_samples=radius_samples, inference_options=inference_options)
    try:
        points = np.asarray(graph['points'], dtype=float)
    except (TypeError, ValueError):
        return tree_builder(graph, annotations=annotations, scale_mm_per_model_unit=scale_mm_per_model_unit,
                            radius_samples=radius_samples, inference_options=inference_options)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all() or not len(points):
        return tree_builder(graph, annotations=annotations, scale_mm_per_model_unit=scale_mm_per_model_unit,
                            radius_samples=radius_samples, inference_options=inference_options)
    ledger = graph_measurement_ledger(points, graph['adjacency'], coordinate_units='model_units')
    adjacency = graph['adjacency']; components = ledger['component_node_ids']
    bases, overrides, endpoints, types, confirmed, confirmed_components, normalized = _annotations(
        annotations, len(points), components, adjacency)
    bridges = graph_bridges(adjacency) if ledger['source_cycle_rank'] else set()
    targets = {}
    owner = {node: cid for cid, nodes in enumerate(components) for node in nodes}
    component_edges = [[] for _ in components]
    for edge in ledger['source_edge_ids']:
        component_edges[owner[edge[0]]].append(tuple(edge))
    for cid, nodes in enumerate(components):
        root = [n for n in nodes if n in bases]
        edges = component_edges[cid]
        if len(root) == 1 and len(edges)-len(nodes)+1 > 0 and any(e in bridges for e in edges):
            targets[cid] = (root[0], edges)
    kwargs = dict(scale_mm_per_model_unit=scale_mm_per_model_unit,
                  radius_samples=radius_samples, inference_options=inference_options)
    if not targets:
        return tree_builder(graph, annotations=annotations, **kwargs)
    if confirmed_components.intersection(targets):
        raise ValueError('CYCLIC_COMPONENT_GLOBAL_CONFIRMATION_UNAVAILABLE')
    safe = dict(normalized, root_types=[], confirmed_axes=[],
        junction_continuations=[r for r in normalized['junction_continuations'] if owner[r['node_id']] not in targets])
    result = tree_builder(graph, annotations=safe, **kwargs)
    unaffected = {a['axis_id'] for a in result['axes'] if a['component_id'] not in targets}
    if unaffected.intersection(set(types) | confirmed):
        restored = dict(safe, root_types=[r for r in normalized['root_types'] if r['axis_id'] in unaffected],
                        confirmed_axes=[n for n in confirmed if n in unaffected])
        result = tree_builder(graph, annotations=restored, **kwargs)
    axes = [a for a in result['axes'] if a['component_id'] not in targets]
    orders = [e for e in result['rooted_edge_orders'] if owner[e['edge_nodes'][0]] not in targets]
    path_ids = {tuple(e): f"path{b['branch_id']:05d}" for b in ledger['branches'] for e in b['source_edge_ids']}
    samples = {}
    for sample in radius_samples or []:
        if sample.get('edge_nodes') is not None:
            samples.setdefault(tuple(sorted(sample['edge_nodes'])), []).append(sample)
    all_regions = []; consumed = {n for n in overrides if owner[n] not in targets}
    for cid, (root, source_edges) in targets.items():
        distance = {root: 0}; queue = deque([root])
        while queue:
            node = queue.popleft()
            for other in adjacency[node]:
                if other not in distance:
                    distance[other] = distance[node]+1; queue.append(other)
        tree_edges = [e for e in source_edges if e in bridges]
        core_edges = [e for e in source_edges if e not in bridges]
        core_nodes = {n for e in core_edges for n in e}
        row = result['components'][cid]
        row.update(status='PARTIAL_ROOTED_CYCLIC_COMPONENT', axis_ids=[],
            unresolved_junction_node_ids=[], auto_inferred_junction_node_ids=[], review_recommended_junction_node_ids=[],
            cyclic_core_edge_count=len(core_edges), analysis_region_ids=[],
            global_hierarchy_status='UNAVAILABLE_ACROSS_CYCLIC_CORE')
        total = math.fsum(math.dist(points[a], points[b]) for a, b in source_edges)
        core_length = math.fsum(math.dist(points[a], points[b]) for a, b in core_edges)
        row['cyclic_core_length_model_units'] = core_length
        row['rooted_fragment_arc_fraction'] = None if not total else (total-core_length)/total
        chunks = [(False, ns, es) for ns, es in edge_regions(tree_edges)]
        chunks += [(True, ns, es) for ns, es in edge_regions(core_edges)]
        for index, (core, nodes, edges) in enumerate(chunks):
            mapping = {n: i for i, n in enumerate(nodes)}
            sub_adj = [[] for _ in nodes]
            for a, b in edges:
                sub_adj[mapping[a]].append(mapping[b]); sub_adj[mapping[b]].append(mapping[a])
            entry = None if core else min(nodes, key=lambda n: (distance[n], n))
            region_id = f'component{cid:05d}_region{index:05d}'
            entry_is_actual = entry == root
            region = dict(region_id=region_id, component_id=cid, kind='CYCLIC_CORE' if core else 'ROOTED_BRIDGE_FRAGMENT',
                node_ids=nodes, edge_ids=[list(e) for e in edges], entry_node_id=entry,
                entry_is_explicit_root_base=entry_is_actual,
                boundary_node_ids=sorted(set(nodes).intersection(core_nodes)),
                length_model_units=math.fsum(math.dist(points[a], points[b]) for a, b in edges),
                global_hierarchy_available=entry_is_actual and not core,
                crossing_route_selected=False, source_edges_removed=False)
            all_regions.append(region); row['analysis_region_ids'].append(region_id)
            sub_annotations = dict(root_bases=[] if core else [dict(node_id=mapping[entry])],
                endpoint_types=[dict(node_id=mapping[n], type=t) for n,t in endpoints.items()
                    if n in mapping and t != 'root_base' and n != entry], junction_continuations=[])
            if entry_is_actual:
                sub_annotations['root_bases'][0]['root_type'] = bases[root]
            edge_set = set(edges)
            for n, (incoming, outgoing) in overrides.items():
                if n not in mapping or core:
                    continue
                if tuple(sorted((n, incoming))) not in edge_set or (outgoing is not None and tuple(sorted((n,outgoing))) not in edge_set):
                    continue
                sub_annotations['junction_continuations'].append(dict(node_id=mapping[n],incoming_node_id=mapping[incoming],
                    outgoing_node_id=None if outgoing is None else mapping[outgoing]))
                consumed.add(n)
            sub_samples = None if radius_samples is None else []
            if sub_samples is not None:
                for edge in edges:
                    for sample in samples.get(edge, []):
                        sub_samples.append(dict(sample, edge_nodes=[mapping[n] for n in sample['edge_nodes']],
                            node_id=None if sample.get('node_id') is None else mapping[sample['node_id']]))
            sub = tree_builder(dict(points=points[nodes].tolist(), adjacency=sub_adj), annotations=sub_annotations,
                scale_mm_per_model_unit=scale_mm_per_model_unit, radius_samples=sub_samples, inference_options=inference_options)
            identities = {a['axis_id']: _axis_id([nodes[n] for n in a['node_ids']], not core) for a in sub['axes']}
            for axis in sub['axes']:
                axis['axis_id'] = identities[axis['axis_id']]
                axis['node_ids'] = [nodes[n] for n in axis['node_ids']]
                axis['edge_ids'] = [sorted([nodes[a], nodes[b]]) for a,b in axis['edge_ids']]
                axis['component_id'] = cid; axis['analysis_region_id'] = region_id
                axis['axis_scope'] = region['kind']
                axis['cyclic_boundary_node_ids'] = sorted(set(axis['node_ids']).intersection(core_nodes))
                for key in ('base_node_id','tip_node_id','attachment_node_id'):
                    if axis.get(key) is not None:
                        axis[key] = nodes[axis[key]]
                axis['parent_axis_id'] = identities.get(axis['parent_axis_id'])
                axis['children_axis_ids'] = [identities[n] for n in axis['children_axis_ids']]
                axis['source_path_ids'] = sorted({path_ids[tuple(e)] for e in axis['edge_ids']})
                axis['local_candidate_developmental_order'] = axis['candidate_developmental_order']
                axis['local_candidate_strahler_orders'] = axis['candidate_strahler_orders']
                axis['global_hierarchy_status'] = 'ROOT_BASE_CONNECTED_CANDIDATE' if entry_is_actual and not core else 'UNAVAILABLE_ACROSS_CYCLIC_CORE'
                if not entry_is_actual:
                    axis['candidate_developmental_order'] = None
                axis['developmental_order'] = None
                axis['candidate_strahler_orders'] = axis['strahler_orders'] = None
                axis['candidate_strahler_order_base'] = axis['strahler_order_base'] = None
                if core:
                    axis['relationship_status'] = 'CYCLIC_CORE_RETAINED_GEOMETRY'
                else:
                    if axis['confirmation_status'] != 'UNRESOLVED':
                        axis['confirmation_status'] = 'CANDIDATE'
                    axis['relationship_status'] = 'LOCAL_FRAGMENT_PARENT_RELATION' if axis['parent_axis_id'] else 'ROOT_RELATIVE_FRAGMENT_START'
                    if not entry_is_actual and axis['node_ids'][0] == entry:
                        axis['basal_endpoint_type'] = 'unknown'
                axes.append(axis); row['axis_ids'].append(axis['axis_id'])
            for decision in sub['junction_decisions']:
                decision.update(component_id=cid, analysis_region_id=region_id,
                    node_id=nodes[decision['node_id']], incoming_node_id=nodes[decision['incoming_node_id']],
                    outgoing_node_id=None if decision['outgoing_node_id'] is None else nodes[decision['outgoing_node_id']])
                for alternative in decision['alternatives']:
                    alternative['outgoing_node_id'] = nodes[alternative['outgoing_node_id']]
                result['junction_decisions'].append(decision)
            for key in ('unresolved_junction_node_ids','auto_inferred_junction_node_ids','review_recommended_junction_node_ids'):
                row[key].extend(nodes[n] for c in sub['components'] for n in c[key])
            if not core:
                orders.extend(dict(edge_nodes=[nodes[n] for n in e['edge_nodes']], candidate_strahler_order=None,
                    strahler_order=None,local_candidate_strahler_order=e['candidate_strahler_order'],analysis_region_id=region_id)
                    for e in sub['rooted_edge_orders'])
    if set(overrides) != consumed:
        raise ValueError('CONTINUATION_CROSSES_UNRESOLVED_CYCLIC_CORE')
    by_id = {a['axis_id']: a for a in axes}
    if (set(types) | confirmed) - set(by_id):
        raise ValueError('UNKNOWN_AXIS_ID_OR_CHANGED_AXIS_PARTITION')
    if any(by_id[n]['component_id'] in targets for n in confirmed):
        raise ValueError('CYCLIC_FRAGMENT_GLOBAL_CONFIRMATION_UNAVAILABLE')
    for axis in axes:
        if axis['axis_id'] in types:
            axis.update(root_type=types[axis['axis_id']],root_type_status='AXIS_ANNOTATION')
        if axis['axis_id'] in confirmed and axis['confirmation_status'] != 'CONFIRMED':
            if axis['confirmation_status'] in ('UNRESOLVED','GEOMETRY_ONLY'):
                raise ValueError('AXIS_CONFIRMATION_REQUIRES_ORIENTED_ROOT_AXIS')
            axis['confirmation_status'] = 'CONFIRMED'
            parent_id = axis.get('parent_axis_id')
            if parent_id is None or by_id[parent_id]['developmental_order'] is not None:
                axis['developmental_order'] = axis['candidate_developmental_order']
    actual = [tuple(e) for a in axes for e in a['edge_ids']]
    expected = {tuple(e) for e in ledger['source_edge_ids']}
    if len(actual) != len(set(actual)) or set(actual) != expected:
        raise AssertionError('CYCLIC_REGION_EDGE_PARTITION')
    if not math.isclose(math.fsum(a['length_model_units'] for a in axes), ledger['measurement_graph_length'],rel_tol=1e-12):
        raise AssertionError('CYCLIC_REGION_LENGTH_PARTITION')
    result.update(axes=axes, rooted_edge_orders=orders, annotations=normalized,
        status='PARTIAL_ROOTED_ARCHITECTURE', analysis_regions=all_regions,
        cycle_analysis_method='rooted-bridge-fragments-v001', global_cyclic_hierarchy_available=False)
    result['qualifications'].append('Bridge fragments are uniquely root-relative; cyclic cores retain all original edges. No cross-core route, complete biological axis or global cyclic hierarchy is inferred.')
    return result
