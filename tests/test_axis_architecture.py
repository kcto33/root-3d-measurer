"""Root identity, topology, coverage and annotation invariants on explicit graphs."""
import copy
import json
import math

import pytest

from root_axis_architecture import build_axis_architecture


def graph(points, edges):
    adjacency = [[] for _ in points]
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    return dict(points=points, adjacency=adjacency)


def branched():
    # A continuous three-unit main axis with two lateral roots of known origins.
    return graph([[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 0, 3],
                  [1, 0, 1.25], [2, 0, 1.5], [0, 1, 2.2]],
                 [(0, 1), (1, 2), (2, 3), (1, 4), (4, 5), (2, 6)])


def annotation(**extra):
    return dict(root_bases=[dict(node_id=0, root_type='primary')], **extra)


def assert_partition(result, source):
    expected = sorted((a, b) for a, ns in enumerate(source['adjacency']) for b in ns if a < b)
    actual = [tuple(e) for axis in result['axes'] for e in axis['edge_ids']]
    assert sorted(actual) == expected
    assert len(set(actual)) == len(actual)
    expected_length = math.fsum(math.dist(source['points'][a], source['points'][b]) for a, b in expected)
    assert sum(a['length_model_units'] for a in result['axes']) == pytest.approx(expected_length)
    assert result['edge_partition']['all_source_edges_owned_exactly_once']


def test_no_base_keeps_geometric_paths_without_inventing_biology():
    source = branched(); before = copy.deepcopy(source)
    result = build_axis_architecture(source)
    assert result['status'] == 'GEOMETRY_ONLY'
    assert result['components'][0]['candidate_root_base_node_ids'] == [0, 3, 5, 6]
    assert all(a['parent_axis_id'] is None and a['developmental_order'] is None for a in result['axes'])
    assert all(a['confirmation_status'] == 'GEOMETRY_ONLY' for a in result['axes'])
    assert all(a['apical_endpoint_type'] == 'unknown' for a in result['axes'])
    assert source == before
    assert_partition(result, source)


def test_continuous_root_axes_cross_junctions_but_remain_candidates():
    source = branched()
    result = build_axis_architecture(source, annotations=annotation())
    main = next(a for a in result['axes'] if a['parent_axis_id'] is None)
    assert main['node_ids'] == [0, 1, 2, 3]
    assert main['length_model_units'] == 3
    assert main['candidate_developmental_order'] == 0
    assert main['developmental_order'] is None
    assert main['candidate_strahler_orders'] == [2, 2, 1]
    assert main['strahler_orders'] is None
    children = [a for a in result['axes'] if a['parent_axis_id'] == main['axis_id']]
    assert [a['attachment_arc_model_units'] for a in children] == [1, 2]
    assert [a['node_ids'] for a in children] == [[1, 4, 5], [2, 6]]
    assert all(a['confirmation_status'] == 'CANDIDATE' for a in result['axes'])
    assert all(j['status'] == 'INFERRED_STRAIGHTEST_CANDIDATE' for j in result['junction_decisions'])
    assert_partition(result, source)


def test_confirmation_does_not_infer_tips_and_strahler_is_per_edge():
    source = branched()
    result = build_axis_architecture(source, annotations=annotation(confirmed_components=[0],
        endpoint_types=[dict(node_id=3, type='root_tip'), dict(node_id=5, type='cut_end')]))
    assert result['status'] == 'CONFIRMED_ARCHITECTURE'
    assert [a['developmental_order'] for a in result['axes']] == [0, 1, 1]
    main = result['axes'][0]
    assert main['strahler_orders'] == [2, 2, 1]
    assert main['strahler_order_base'] == 2
    assert main['apical_endpoint_type'] == 'root_tip'
    assert result['axes'][1]['apical_endpoint_type'] == 'cut_end'
    assert result['axes'][2]['apical_endpoint_type'] == 'unknown'
    assert all(r['strahler_order'] is not None for r in result['rooted_edge_orders'])
    assert_partition(result, source)
    json.dumps(result, allow_nan=False)


def test_manual_continuation_changes_identity_and_developmental_order_only():
    source = branched()
    natural = build_axis_architecture(source, annotations=annotation())
    result = build_axis_architecture(source, annotations=annotation(confirmed_components=[0],
        junction_continuations=[dict(node_id=1, incoming_node_id=0, outgoing_node_id=4)]))
    assert [a['node_ids'] for a in result['axes']] == [[0, 1, 4, 5], [1, 2, 3], [2, 6]]
    assert [a['developmental_order'] for a in result['axes']] == [0, 1, 2]
    assert result['axes'][0]['axis_id'] != natural['axes'][0]['axis_id']
    # Strahler belongs to the unchanged rooted graph, not to a main-axis choice.
    assert [(r['edge_nodes'], r['candidate_strahler_order']) for r in result['rooted_edge_orders']] == [
        (r['edge_nodes'], r['candidate_strahler_order']) for r in natural['rooted_edge_orders']]
    assert_partition(result, source)


def test_symmetric_branch_requires_a_choice_before_confirmation():
    source = graph([[0, 0, 0], [0, 0, 1], [1, 0, 2], [-1, 0, 2]], [(0, 1), (1, 2), (1, 3)])
    candidate = build_axis_architecture(source, annotations=annotation())
    assert candidate['status'] == 'UNRESOLVED_ARCHITECTURE'
    assert candidate['components'][0]['unresolved_junction_node_ids'] == [1]
    assert all(a['confirmation_status'] == 'UNRESOLVED' for a in candidate['axes'])
    assert all(a['candidate_developmental_order'] is None for a in candidate['axes'])
    assert all(a['relationship_status'] == 'UNRESOLVED_CONTINUATION' for a in candidate['axes'])
    assert candidate['rooted_edge_orders']  # Graph Strahler does not need a continuation choice.
    assert_partition(candidate, source)
    with pytest.raises(ValueError, match='UNRESOLVED_CONTINUATION'):
        build_axis_architecture(source, annotations=annotation(confirmed_components=[0]))
    confirmed = build_axis_architecture(source, annotations=annotation(confirmed_components=[0],
        junction_continuations=[dict(node_id=1, incoming_node_id=0, outgoing_node_id=None)]))
    assert confirmed['axes'][0]['node_ids'] == [0, 1]
    assert confirmed['axes'][0]['apical_endpoint_type'] == 'unknown'
    assert [a['developmental_order'] for a in confirmed['axes']] == [0, 1, 1]


def test_internal_crown_creates_multiple_basal_axes_without_single_main_root():
    source = graph([[0, 0, 0], [0, 0, 1], [-1, 0, 0], [1, 0, 0]], [(1, 0), (1, 2), (1, 3)])
    result = build_axis_architecture(source, annotations=dict(root_bases=[dict(node_id=1, root_type='crown')],
                                                             confirmed_components=[0]))
    assert len(result['axes']) == 3
    assert all(a['parent_axis_id'] is None and a['developmental_order'] == 0 for a in result['axes'])
    assert all(a['basal_endpoint_type'] == 'root_base' and a['node_ids'][0] == 1 for a in result['axes'])
    assert not result['junction_decisions']
    assert_partition(result, source)


def test_multiple_bases_supported_across_components_and_ambiguous_within_one():
    source = graph([[0, 0, 0], [0, 0, 1], [2, 0, 0], [2, 0, 1]], [(0, 1), (2, 3)])
    result = build_axis_architecture(source, annotations=dict(root_bases=[dict(node_id=0), dict(node_id=2)],
                                                             confirmed_components=[0, 1]))
    assert result['status'] == 'CONFIRMED_ARCHITECTURE'
    assert [a['developmental_order'] for a in result['axes']] == [0, 0]
    joined = graph(source['points'], [(0, 1), (1, 2), (2, 3)])
    ambiguous = build_axis_architecture(joined, annotations=dict(root_bases=[dict(node_id=0), dict(node_id=3)]))
    assert ambiguous['components'][0]['status'] == 'BLOCKED_MULTIPLE_BASES_IN_COMPONENT'
    assert all(a['orientation_status'] == 'UNORIENTED' for a in ambiguous['axes'])
    assert_partition(ambiguous, joined)


def test_cycle_is_preserved_and_never_turned_into_a_biological_tree():
    source = graph([[0, 0, 0], [1, 0, 0], [0, 1, 0]], [(0, 1), (1, 2), (2, 0)])
    result = build_axis_architecture(source, annotations=annotation())
    assert result['components'][0]['status'] == 'BLOCKED_CYCLIC_COMPONENT'
    assert result['preserved_cycle_rank'] == 1
    assert result['axes'][0]['closed']
    assert result['axes'][0]['candidate_developmental_order'] is None
    assert_partition(result, source)
    with pytest.raises(ValueError, match='UNIQUE_ROOTED_TREE'):
        build_axis_architecture(source, annotations=annotation(confirmed_components=[0]))


def test_individual_axis_confirmation_does_not_confirm_unreviewed_ancestors():
    source = branched()
    candidate = build_axis_architecture(source, annotations=annotation())
    child_id = candidate['axes'][1]['axis_id']
    result = build_axis_architecture(source, annotations=annotation(confirmed_axes=[child_id],
        root_types=[dict(axis_id=child_id, type='lateral')]))
    child = next(a for a in result['axes'] if a['axis_id'] == child_id)
    assert child['confirmation_status'] == 'CONFIRMED'
    assert child['root_type'] == 'lateral'
    assert child['developmental_order'] is None
    assert child['candidate_developmental_order'] == 1
    assert child['relationship_status'] == 'CANDIDATE_PARENT_RELATION'
    assert child['strahler_orders'] is None


def test_old_confirmation_cannot_follow_reversed_or_changed_root_axis():
    source = graph([[0, 0, 0], [0, 0, 1], [0, 0, 2]], [(0, 1), (1, 2)])
    old = build_axis_architecture(source, annotations=annotation())['axes'][0]['axis_id']
    with pytest.raises(ValueError, match='UNKNOWN_AXIS_ID'):
        build_axis_architecture(source, annotations=dict(root_bases=[dict(node_id=2)], confirmed_axes=[old]))


@pytest.mark.parametrize('bad', [
    {'invented': True},
    {'root_bases': [{'node_id': True}]},
    {'root_bases': [{'node_id': 100}]},
    {'root_bases': [{'node_id': 0}, {'node_id': 0}]},
    {'root_bases': [{'node_id': 0, 'unexpected': 1}]},
    {'endpoint_types': [{'node_id': 1, 'type': 'root_tip'}]},
    {'endpoint_types': [{'node_id': 3, 'type': 'live'}]},
    {'confirmed_components': [True]},
    {'confirmed_axes': ['not-an-axis']},
    {'root_bases': [{'node_id': 0}], 'junction_continuations': [
        {'node_id': 1, 'incoming_node_id': 2, 'outgoing_node_id': 4}]},
    {'root_bases': [{'node_id': 0}], 'junction_continuations': [
        {'node_id': 4, 'incoming_node_id': 1, 'outgoing_node_id': 5}]},
])
def test_invalid_or_stale_annotations_fail_explicitly(bad):
    with pytest.raises(ValueError):
        build_axis_architecture(branched(), annotations=bad)


def test_long_chain_is_iterative_and_units_do_not_change_partition():
    count = 8000
    source = graph([[0, 0, i] for i in range(count)], [(i, i + 1) for i in range(count - 1)])
    result = build_axis_architecture(source, annotations=annotation(confirmed_components=[0]), scale_mm_per_model_unit=10)
    assert len(result['axes']) == 1
    assert result['axes'][0]['length_model_units'] == count - 1
    assert result['axes'][0]['developmental_order'] == 0
    assert set(result['axes'][0]['strahler_orders']) == {1}
    assert result['scale_mm_per_model_unit'] == 10
    assert_partition(result, source)


def test_empty_and_isolated_graphs_do_not_claim_confirmed_roots():
    empty = build_axis_architecture(graph([], []))
    assert empty['axes'] == [] and empty['status'] == 'GEOMETRY_ONLY'
    isolated = build_axis_architecture(graph([[0, 0, 0]], []), annotations=annotation())
    assert isolated['isolated_node_ids'] == [0]
    assert isolated['components'][0]['status'] == 'NO_EDGES'
    assert isolated['axes'] == []


@pytest.mark.parametrize('scale', [0, -1, True, float('nan'), float('inf'), '10', 1j])
def test_invalid_scales_are_rejected(scale):
    with pytest.raises(ValueError, match='POSITIVE_FINITE_SCALE'):
        build_axis_architecture(branched(), scale_mm_per_model_unit=scale)
