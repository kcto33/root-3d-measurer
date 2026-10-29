"""Explicitly rooted, reviewable root-axis interpretations of an unchanged graph.

Geometric paths are not biological roots.  A root base or crown must be annotated
before a tree component can be oriented.  Radius-qualified, windowed tangent, downstream length and hierarchy evidence
can identify high-confidence *geometric* continuations. Low-confidence cases
remain reviewable. Cycles and competing bases are never cut into a tree.

All coordinates and arc lengths in this module use model units.  The annotation
document must be bound to the source measurement by the caller.  Axis identifiers
include their edge set and oriented base, so an old confirmation cannot silently
move to a changed axis or reversed orientation.
"""
import hashlib
import json
import math
from collections import deque

import numpy as np

from root_graph_measurement_ledger import graph_measurement_ledger


METHOD = 'root-axis-architecture-v001'
ANNOTATION_SCHEMA = 'root-axis-annotations-v001'
ENDPOINT_TYPES = {'root_tip', 'cut_end', 'unknown', 'root_base'}
_ANNOTATION_KEYS = {'schema_version', 'root_bases', 'junction_continuations',
                    'endpoint_types', 'root_types', 'confirmed_axes',
                    'confirmed_components'}


def _edge(a, b):
    return (a, b) if a < b else (b, a)


def _identifier(value, name, upper=None):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(name + '_INTEGER_REQUIRED')
    value = int(value)
    if value < 0 or (upper is not None and value >= upper):
        raise ValueError(name + '_UNKNOWN_ID')
    return value


def _records(value, name):
    if not isinstance(value, list):
        raise ValueError(name + '_LIST_REQUIRED')
    return value


def _record(value, allowed, required, name):
    if not isinstance(value, dict) or set(value) - allowed or required - set(value):
        raise ValueError(name + '_INVALID_FIELDS')
    return value


def _text(value, name):
    if not isinstance(value, str) or not value.strip() or len(value) > 200:
        raise ValueError(name + '_NONEMPTY_TEXT_REQUIRED')
    return value.strip()


def _annotations(value, count, components, adjacency):
    value = {} if value is None else value
    if not isinstance(value, dict) or set(value) - _ANNOTATION_KEYS:
        raise ValueError('UNKNOWN_ARCHITECTURE_ANNOTATION_FIELD')
    if value.get('schema_version', ANNOTATION_SCHEMA) != ANNOTATION_SCHEMA:
        raise ValueError('UNKNOWN_ARCHITECTURE_ANNOTATION_SCHEMA')
    bases, continuation, endpoints, types = {}, {}, {}, {}
    for item in _records(value.get('root_bases', []), 'ROOT_BASES'):
        _record(item, {'node_id', 'root_type'}, {'node_id'}, 'ROOT_BASE')
        node = _identifier(item['node_id'], 'ROOT_BASE_NODE', count)
        if node in bases:
            raise ValueError('DUPLICATE_ROOT_BASE_NODE')
        bases[node] = _text(item.get('root_type', 'unknown'), 'ROOT_TYPE')
    for item in _records(value.get('junction_continuations', []), 'CONTINUATIONS'):
        _record(item, {'node_id', 'incoming_node_id', 'outgoing_node_id'},
                {'node_id', 'incoming_node_id', 'outgoing_node_id'}, 'CONTINUATION')
        node = _identifier(item['node_id'], 'CONTINUATION_NODE', count)
        incoming = _identifier(item['incoming_node_id'], 'INCOMING_NODE', count)
        outgoing = item['outgoing_node_id']
        if outgoing is not None:
            outgoing = _identifier(outgoing, 'OUTGOING_NODE', count)
        if node in continuation:
            raise ValueError('DUPLICATE_CONTINUATION_NODE')
        if incoming not in adjacency[node] or (outgoing is not None and
                (outgoing not in adjacency[node] or outgoing == incoming)):
            raise ValueError('CONTINUATION_INCIDENT_DISTINCT_EDGES_REQUIRED')
        continuation[node] = (incoming, outgoing)
    for item in _records(value.get('endpoint_types', []), 'ENDPOINT_TYPES'):
        _record(item, {'node_id', 'type'}, {'node_id', 'type'}, 'ENDPOINT_TYPE')
        node = _identifier(item['node_id'], 'ENDPOINT_NODE', count)
        kind = item['type']
        if not isinstance(kind, str) or kind not in ENDPOINT_TYPES:
            raise ValueError('UNKNOWN_ENDPOINT_TYPE')
        if node in endpoints:
            raise ValueError('DUPLICATE_ENDPOINT_TYPE_NODE')
        if kind in {'root_tip', 'cut_end'} and (len(adjacency[node]) != 1 or node in bases):
            raise ValueError('ROOT_TIP_OR_CUT_REQUIRES_NONBASE_DEGREE_ONE_NODE')
        if kind == 'root_base' and node not in bases:
            raise ValueError('ROOT_BASE_LABEL_REQUIRES_ROOT_BASE_ANNOTATION')
        endpoints[node] = kind
    for node in bases:
        if endpoints.get(node, 'root_base') not in {'unknown', 'root_base'}:
            raise ValueError('CONFLICTING_ROOT_BASE_ENDPOINT_TYPE')
        endpoints[node] = 'root_base'
    for item in _records(value.get('root_types', []), 'ROOT_TYPES'):
        _record(item, {'axis_id', 'type'}, {'axis_id', 'type'}, 'AXIS_ROOT_TYPE')
        axis_id = _text(item['axis_id'], 'AXIS_ID')
        if axis_id in types:
            raise ValueError('DUPLICATE_AXIS_ROOT_TYPE')
        types[axis_id] = _text(item['type'], 'ROOT_TYPE')
    confirmed_axes = [_text(i, 'AXIS_ID') for i in
                      _records(value.get('confirmed_axes', []), 'CONFIRMED_AXES')]
    confirmed_components = [_identifier(i, 'COMPONENT', len(components)) for i in
                            _records(value.get('confirmed_components', []), 'CONFIRMED_COMPONENTS')]
    if len(set(confirmed_axes)) != len(confirmed_axes) or len(set(confirmed_components)) != len(confirmed_components):
        raise ValueError('DUPLICATE_CONFIRMATION_ID')
    normalized = dict(schema_version=ANNOTATION_SCHEMA,
        root_bases=[dict(node_id=i, root_type=t) for i, t in sorted(bases.items())],
        junction_continuations=[dict(node_id=i, incoming_node_id=p, outgoing_node_id=c)
                                for i, (p, c) in sorted(continuation.items())],
        endpoint_types=[dict(node_id=i, type=t) for i, t in sorted(endpoints.items())],
        root_types=[dict(axis_id=i, type=t) for i, t in sorted(types.items())],
        confirmed_axes=sorted(confirmed_axes), confirmed_components=sorted(confirmed_components))
    return bases, continuation, endpoints, types, set(confirmed_axes), set(confirmed_components), normalized


def _axis_id(nodes, oriented):
    payload = dict(edges=sorted(_edge(a, b) for a, b in zip(nodes, nodes[1:])),
                   oriented_base_node_id=nodes[0] if oriented else None)
    return 'axis_' + hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:16]


def _axis(nodes, component_id, points, lengths, path_ids, endpoints, oriented):
    edges = [_edge(a, b) for a, b in zip(nodes, nodes[1:])]
    arc = [0.]
    for edge in edges:
        arc.append(arc[-1] + lengths[edge])
    return dict(axis_id=_axis_id(nodes, oriented), component_id=component_id,
        node_ids=list(nodes), edge_ids=[list(e) for e in edges],
        points=[points[i].tolist() for i in nodes], arc_positions_model_units=arc,
        length_model_units=math.fsum(lengths[e] for e in edges),
        source_path_ids=sorted({path_ids[e] for e in edges}), closed=nodes[0] == nodes[-1],
        base_node_id=nodes[0] if oriented else None, tip_node_id=nodes[-1] if oriented else None,
        parent_axis_id=None, attachment_node_id=None, attachment_arc_model_units=None,
        children_axis_ids=[], orientation_status='ROOTED_BASE_TO_APEX' if oriented else 'UNORIENTED',
        confirmation_status='CANDIDATE' if oriented else 'GEOMETRY_ONLY',
        relationship_status='BASAL_AXIS_EXPLICIT_ROOT' if oriented else 'ROOT_BASE_REQUIRED',
        root_type='unknown', root_type_status='UNANNOTATED',
        basal_endpoint_type=endpoints.get(nodes[0], 'unknown'),
        apical_endpoint_type=endpoints.get(nodes[-1], 'unknown'),
        candidate_developmental_order=None, developmental_order=None,
        candidate_strahler_orders=None, strahler_orders=None,
        candidate_strahler_order_base=None, strahler_order_base=None)


def build_axis_architecture(graph, *, annotations=None, scale_mm_per_model_unit=None,
                            radius_samples=None, inference_options=None):
    """Partition every graph edge into candidate, automatically inferred or reviewed axes.

    ``annotations`` accepts only ``root_bases`` (node_id, optional root_type),
    ``junction_continuations`` (node_id, incoming_node_id, outgoing_node_id or
    null for explicit termination), ``endpoint_types`` (node_id, type),
    ``root_types`` (axis_id, type), ``confirmed_axes`` and ``confirmed_components``.
    Caller-owned measurement-hash binding is required before importing annotations.

    A crown at an internal node creates one basal axis per outgoing root.  Multiple
    bases in different components are supported.  Competing bases in one component
    remain unoriented.  Developmental order starts at zero for basal axes and rises
    at lateral attachments.  Strahler order is computed per rooted graph edge:
    terminal edges are one; two or more equally highest child orders add one.
    """
    if not isinstance(graph, dict) or 'points' not in graph or 'adjacency' not in graph:
        raise ValueError('POINTS_AND_ADJACENCY_REQUIRED')
    try:
        points = np.asarray(graph['points'], dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError('FINITE_N_BY_3_POINTS_REQUIRED') from exc
    if points.size == 0 and points.ndim == 1:
        points = points.reshape(0, 3)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError('FINITE_N_BY_3_POINTS_REQUIRED')
    if scale_mm_per_model_unit is not None and (isinstance(scale_mm_per_model_unit, bool) or
            not isinstance(scale_mm_per_model_unit, (int, float, np.integer, np.floating)) or
            not math.isfinite(scale_mm_per_model_unit) or scale_mm_per_model_unit <= 0):
        raise ValueError('POSITIVE_FINITE_SCALE_REQUIRED')
    ledger = graph_measurement_ledger(points, graph['adjacency'], coordinate_units='model_units')
    adjacency = [sorted(map(int, ns)) for ns in graph['adjacency']]
    components = ledger['component_node_ids']
    edges = [tuple(e) for e in ledger['source_edge_ids']]
    lengths = {e: math.dist(points[e[0]], points[e[1]]) for e in edges}
    path_ids = {tuple(e): f"path{b['branch_id']:05d}" for b in ledger['branches'] for e in b['source_edge_ids']}
    bases, overrides, endpoints, root_types, confirm_axes, confirm_components, normalized = _annotations(
        annotations, len(points), components, adjacency)
    multi_cue = radius_samples is not None
    if multi_cue:
        from root_axis_auto_inference import radius_support_by_edge, infer_continuation, METHOD as INFERENCE_METHOD
        from root_axis_traits import axis_trait_settings
        inference_config = axis_trait_settings(inference_options)
        support = radius_support_by_edge(radius_samples)
        factor = 1.0 if scale_mm_per_model_unit is None else float(scale_mm_per_model_unit)
    elif inference_options is not None:
        raise ValueError('INFERENCE_OPTIONS_REQUIRE_RADIUS_SAMPLES')
    source_payload = dict(points=points.tolist(), edges=[list(e) for e in edges])
    graph_sha = hashlib.sha256(json.dumps(source_payload, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    axes, component_rows, decisions = [], [], []
    edge_orders = {}
    consumed_overrides = set()
    branches_by_component = {i: [] for i in range(len(components))}
    for branch in ledger['branches']:
        branches_by_component[branch['component_id']].append(branch)
    for component_id, nodes in enumerate(components):
        marked_bases = [node for node in nodes if node in bases]
        edge_count = sum(len(adjacency[i]) for i in nodes) // 2
        cycle_rank = edge_count - len(nodes) + 1
        reason = ('NO_EDGES' if not edge_count else 'BLOCKED_CYCLIC_COMPONENT' if cycle_rank else
                  'ROOT_BASE_REQUIRED' if not marked_bases else
                  'BLOCKED_MULTIPLE_BASES_IN_COMPONENT' if len(marked_bases) != 1 else None)
        row = dict(component_id=component_id, node_count=len(nodes), edge_count=edge_count,
                   cycle_rank=cycle_rank, root_base_node_ids=marked_bases,
                   candidate_root_base_node_ids=[i for i in nodes if len(adjacency[i]) == 1],
                   status=reason or 'ROOTED_CANDIDATE', unresolved_junction_node_ids=[],
                   auto_inferred_junction_node_ids=[], review_recommended_junction_node_ids=[], axis_ids=[])
        if reason:
            for branch in branches_by_component[component_id]:
                axis = _axis(branch['source_node_ids'], component_id, points, lengths, path_ids, endpoints, False)
                axis['relationship_status'] = reason
                axes.append(axis); row['axis_ids'].append(axis['axis_id'])
            component_rows.append(row)
            continue
        base = marked_bases[0]
        parent = {base: None}; children = {}; traversal = []
        queue = deque([base])
        while queue:
            node = queue.popleft(); traversal.append(node)
            child_nodes = [i for i in adjacency[node] if i != parent[node]]
            children[node] = child_nodes
            for child in child_nodes:
                if child in parent:
                    raise AssertionError('TREE_ORIENTATION_REVISITED_NODE')
                parent[child] = node; queue.append(child)
        if len(parent) != len(nodes):
            raise AssertionError('ROOTED_COMPONENT_COVERAGE')
        node_orders = {}
        longest = {}
        for node in reversed(traversal):
            child_orders = [node_orders[c] for c in children[node]]
            maximum = max(child_orders, default=1)
            node_orders[node] = maximum + int(bool(child_orders) and child_orders.count(maximum) >= 2)
            longest[node] = max((lengths[_edge(node, child)] + longest[child] for child in children[node]), default=0.0)
            if parent[node] is not None:
                edge_orders[_edge(parent[node], node)] = node_orders[node]
        continuation = {}
        for node in traversal:
            child_nodes = children[node]
            if parent[node] is None or len(child_nodes) < 2:
                if node in overrides:
                    raise ValueError('CONTINUATION_REQUIRES_ORIENTED_NONBASE_BRANCH_NODE')
                continuation[node] = child_nodes[0] if child_nodes else None
                continue
            incoming = parent[node]
            if multi_cue:
                proposed, proposed_status, alternatives, confidence, evidence = infer_continuation(
                    node, incoming, child_nodes, points=points, parent=parent, children=children,
                    lengths=lengths, node_orders=node_orders, longest=longest,
                    radius_support=support, config=inference_config, factor=factor)
            else:
                vector = points[node] - points[incoming]
                norm = math.hypot(*vector)
                alternatives = []
                for child in child_nodes:
                    other = points[child] - points[node]
                    other_norm = math.hypot(*other)
                    score = None if not norm or not other_norm else float(np.clip(np.dot(vector / norm, other / other_norm), -1, 1))
                    alternatives.append(dict(outgoing_node_id=child, alignment_cosine=score))
                confidence, evidence = 'NOT_ASSESSED', []
            if node in overrides:
                annotated_in, selected = overrides[node]
                if annotated_in != incoming or (selected is not None and selected not in child_nodes):
                    raise ValueError('CONTINUATION_CONFLICTS_WITH_ROOTED_DIRECTION')
                consumed_overrides.add(node)
                status = 'EXPLICIT_TERMINATION' if selected is None else 'EXPLICIT_CONTINUATION'
                confidence = 'USER_SPECIFIED'
            elif multi_cue:
                selected, status = proposed, proposed_status
            else:
                finite = sorted((a for a in alternatives if a['alignment_cosine'] is not None),
                                key=lambda a: (-a['alignment_cosine'], a['outgoing_node_id']))
                tied = (not finite or len(finite) != len(alternatives) or
                        (len(finite) > 1 and abs(finite[0]['alignment_cosine'] - finite[1]['alignment_cosine']) <= 1e-9))
                if tied:
                    selected = None; status = 'UNRESOLVED_CONTINUATION'
                else:
                    selected = finite[0]['outgoing_node_id']; status = 'INFERRED_STRAIGHTEST_CANDIDATE'
            if status == 'UNRESOLVED_CONTINUATION':
                row['unresolved_junction_node_ids'].append(node)
            elif status == 'AUTO_INFERRED_CONTINUATION':
                row['auto_inferred_junction_node_ids'].append(node)
            elif status == 'INFERRED_MULTI_CUE_CANDIDATE':
                row['review_recommended_junction_node_ids'].append(node)
            continuation[node] = selected
            decisions.append(dict(component_id=component_id, node_id=node, incoming_node_id=incoming,
                outgoing_node_id=selected, status=status, alternatives=alternatives,
                confidence_tier=confidence, supporting_evidence=evidence,
                method=INFERENCE_METHOD if multi_cue else 'incident-edge cosine; review candidate',
                biological_identity_validated=False))
        pending = [(base, child, None, None, 0) for child in reversed(children[base])]
        while pending:
            start, child, parent_axis, attach_arc, order = pending.pop()
            axis_nodes = [start, child]; side_roots = []; running_arc = lengths[_edge(start, child)]
            while True:
                node = axis_nodes[-1]
                following = continuation[node]
                for side in children[node]:
                    if side != following:
                        side_roots.append((node, side, running_arc))
                if following is None:
                    break
                axis_nodes.append(following)
                running_arc += lengths[_edge(node, following)]
            axis = _axis(axis_nodes, component_id, points, lengths, path_ids, endpoints, True)
            axis.update(parent_axis_id=parent_axis, attachment_node_id=start if parent_axis else None,
                        attachment_arc_model_units=attach_arc, candidate_developmental_order=order,
                        candidate_strahler_orders=[edge_orders[tuple(e)] for e in axis['edge_ids']])
            axis['candidate_strahler_order_base'] = axis['candidate_strahler_orders'][0]
            axis['relationship_status'] = 'CANDIDATE_PARENT_RELATION' if parent_axis else 'BASAL_AXIS_EXPLICIT_ROOT'
            if parent_axis is None:
                axis.update(root_type=bases[base], root_type_status='ROOT_BASE_ANNOTATION' if bases[base] != 'unknown' else 'UNANNOTATED')
            if row['unresolved_junction_node_ids']:
                axis['confirmation_status'] = 'UNRESOLVED'
                # Stopping at a tied junction only preserves the edge partition;
                # it is not a choice that the mother root biologically ends here.
                axis['candidate_developmental_order'] = None
                axis['relationship_status'] = 'UNRESOLVED_CONTINUATION'
            axes.append(axis); row['axis_ids'].append(axis['axis_id'])
            for node, side, location in reversed(side_roots):
                pending.append((node, side, axis['axis_id'], location, order + 1))
        if row['unresolved_junction_node_ids']:
            row['status'] = 'ROOTED_UNRESOLVED_CONTINUATION'
        elif multi_cue:
            row['status'] = ('ROOTED_REVIEW_RECOMMENDED' if row['review_recommended_junction_node_ids']
                             else 'ROOTED_AUTO_INFERRED')
        component_rows.append(row)
    if set(overrides) != consumed_overrides:
        raise ValueError('CONTINUATION_REQUIRES_UNIQUE_ROOTED_TREE_COMPONENT')
    by_id = {a['axis_id']: a for a in axes}
    if len(by_id) != len(axes):
        raise AssertionError('AXIS_ID_COLLISION')
    if (set(root_types) | confirm_axes) - set(by_id):
        raise ValueError('UNKNOWN_AXIS_ID_OR_CHANGED_AXIS_PARTITION')
    for component_id in confirm_components:
        row = component_rows[component_id]
        if row['status'] == 'ROOTED_UNRESOLVED_CONTINUATION':
            raise ValueError('UNRESOLVED_CONTINUATION_REQUIRES_EXPLICIT_OVERRIDE')
        if row['status'] not in ('ROOTED_CANDIDATE', 'ROOTED_AUTO_INFERRED', 'ROOTED_REVIEW_RECOMMENDED'):
            raise ValueError('COMPONENT_CONFIRMATION_REQUIRES_UNIQUE_ROOTED_TREE')
        confirm_axes.update(row['axis_ids'])
    for axis_id in confirm_axes:
        axis = by_id[axis_id]
        if axis['confirmation_status'] == 'UNRESOLVED':
            raise ValueError('UNRESOLVED_CONTINUATION_REQUIRES_EXPLICIT_OVERRIDE')
        if axis['orientation_status'] != 'ROOTED_BASE_TO_APEX':
            raise ValueError('AXIS_CONFIRMATION_REQUIRES_ORIENTED_ROOT_AXIS')
    for axis in axes:
        axis_id = axis['axis_id']; parent_id = axis['parent_axis_id']
        if axis_id in root_types:
            axis.update(root_type=root_types[axis_id], root_type_status='AXIS_ANNOTATION')
        if parent_id:
            by_id[parent_id]['children_axis_ids'].append(axis_id)
        if axis['orientation_status'] == 'ROOTED_BASE_TO_APEX' and component_rows[axis['component_id']]['status'] == 'ROOTED_AUTO_INFERRED':
            axis['confirmation_status'] = 'AUTO_INFERRED'
            if parent_id is None or by_id[parent_id]['developmental_order'] is not None:
                axis['developmental_order'] = axis['candidate_developmental_order']
            if parent_id:
                axis['relationship_status'] = 'AUTO_INFERRED_PARENT_RELATION'
        if axis_id in confirm_axes:
            axis['confirmation_status'] = 'CONFIRMED'
            if parent_id is None or by_id[parent_id]['developmental_order'] is not None:
                axis['developmental_order'] = axis['candidate_developmental_order']
            if parent_id:
                axis['relationship_status'] = ('CONFIRMED_PARENT_RELATION' if parent_id in confirm_axes
                                               else 'CANDIDATE_PARENT_RELATION')
    for row in component_rows:
        component_axes = [by_id[i] for i in row['axis_ids']]
        confirmed = [a for a in component_axes if a['confirmation_status'] == 'CONFIRMED']
        if confirmed and len(confirmed) == len(component_axes):
            row['status'] = 'ROOTED_CONFIRMED'
            for axis in component_axes:
                axis['strahler_orders'] = list(axis['candidate_strahler_orders'])
                axis['strahler_order_base'] = axis['candidate_strahler_order_base']
        elif confirmed:
            row['status'] = 'ROOTED_PARTIALLY_CONFIRMED'
        elif component_axes and all(a['confirmation_status'] == 'AUTO_INFERRED' for a in component_axes):
            row['status'] = 'ROOTED_AUTO_INFERRED'
            for axis in component_axes:
                axis['strahler_orders'] = list(axis['candidate_strahler_orders'])
                axis['strahler_order_base'] = axis['candidate_strahler_order_base']
    assigned = [tuple(e) for a in axes for e in a['edge_ids']]
    if len(assigned) != len(edges) or len(set(assigned)) != len(edges) or set(assigned) != set(edges):
        raise AssertionError('ROOT_AXES_MUST_PARTITION_ALL_ORIGINAL_EDGES_EXACTLY_ONCE')
    total = math.fsum(a['length_model_units'] for a in axes)
    if not math.isclose(total, ledger['measurement_graph_length'], rel_tol=1e-12, abs_tol=1e-12):
        raise AssertionError('ROOT_AXIS_PARTITION_CHANGED_LENGTH')
    rooted = [a for a in axes if a['orientation_status'] == 'ROOTED_BASE_TO_APEX']
    confirmed_count = sum(a['confirmation_status'] == 'CONFIRMED' for a in axes)
    auto_count = sum(a['confirmation_status'] == 'AUTO_INFERRED' for a in axes)
    status = ('GEOMETRY_ONLY' if not rooted else
              'CONFIRMED_ARCHITECTURE' if confirmed_count == len(axes) and not ledger['isolated_node_ids'] else
              'AUTO_INFERRED_ARCHITECTURE' if auto_count and auto_count + confirmed_count == len(axes) and not ledger['isolated_node_ids'] else
              'UNRESOLVED_ARCHITECTURE' if any(a['confirmation_status'] == 'UNRESOLVED' for a in axes) and len(rooted) == len(axes) else
              'CANDIDATE_ARCHITECTURE' if len(rooted) == len(axes) and not ledger['isolated_node_ids'] else
              'PARTIAL_ARCHITECTURE')
    edge_components = {tuple(e): a['component_id'] for a in axes for e in a['edge_ids']}
    return dict(schema_version=METHOD, status=status, source_graph_sha256=graph_sha,
        coordinate_units='model_units', scale_mm_per_model_unit=None if scale_mm_per_model_unit is None else float(scale_mm_per_model_unit),
        axes=axes, components=component_rows, junction_decisions=decisions,
        rooted_edge_orders=[dict(edge_nodes=list(e), candidate_strahler_order=o,
            strahler_order=o if component_rows[edge_components[e]]['status'] in ('ROOTED_CONFIRMED', 'ROOTED_AUTO_INFERRED') else None)
            for e, o in sorted(edge_orders.items())],
        annotations=normalized,
        edge_partition=dict(source_edge_count=len(edges), assigned_edge_count=len(assigned),
            all_source_edges_owned_exactly_once=True, total_length_model_units=total),
        isolated_node_ids=list(ledger['isolated_node_ids']), preserved_cycle_rank=ledger['source_cycle_rank'],
        input_graph_unchanged=True, automatic_base_selection=False,
        inference_method=INFERENCE_METHOD if multi_cue else 'incident-edge cosine',
        inference_uses_qualified_radii=multi_cue, confidence_is_calibrated_probability=False,
        developmental_order_definition='basal axes = 0; each lateral-axis attachment increments by one',
        strahler_order_definition='per rooted graph edge; terminal edge = 1; at least two maximum-order children increment the maximum',
        qualifications=['Multi-cue high-confidence choices are geometric AUTO_INFERRED, never biological proof; uncertain choices require review.' if multi_cue else 'Automatic straightness choices are candidates until explicitly reviewed.',
            'A graph endpoint is not a biological root tip without annotation.',
            'Cycles and competing root bases are preserved without tree repair.',
            'A single continuous axis may contain several Strahler edge orders.',
            'Caller must bind imported annotations to the measurement hash.'])
