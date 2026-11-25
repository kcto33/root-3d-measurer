"""Local search distance and boundary stability are separate requirements."""
import itertools
import numpy as np
from root_measurer.runtime import activate
activate()
from root_reconstruction_measurement import ReconstructedSolid
from root_branch_section_measurement import BranchOwnership
from root_native_endpoint_retry import retry_short_boundary


def probe(right, *, foreign=False):
    solid = ReconstructedSolid(np.array(list(itertools.product(
        range(-2, right), range(-3, 4), range(-3, 4))), np.int64), 1., {'units': 'model'})
    points = [[.5, .5, .5], [-.5, .5, .5]]
    adjacency = [[1], [0]]
    if foreign:
        points.extend([[2.2, .5, .5], [2.8, .5, .5]])
        adjacency.extend([[3], [2]])
    record = dict(node_id=0, status='endpoint_unresolved', reason='INITIAL_SECTION_REFUSAL')
    ownership = BranchOwnership(points, adjacency)
    graph, _, records = retry_short_boundary(solid, points, adjacency, [record], ownership)
    return graph, records[0]


def test_stable_local_boundary_can_be_beyond_two_pitches():
    graph, record = probe(3)
    assert record['status'] == 'geometric_surface_estimate'
    assert graph[0] == [3., .5, .5]
    evidence = record['short_boundary_evidence']
    assert evidence['spread_pass']
    assert evidence['thresholds']['spread_pitches'] == 1
    assert all(p['entire_segment_source_supported'] for p in evidence['probes'])
    assert not record['physical_tip_qualified']


def test_nonlocal_boundary_still_refuses():
    _, record = probe(5)
    assert record['status'] == 'endpoint_unresolved'
    assert record['short_boundary_evidence']['refusal'] == 'SHORT_ENDPOINT_EXIT_NOT_LOCAL'


def test_source_supported_but_unstable_boundary_still_refuses():
    _, record = probe(4)
    assert record['status'] == 'endpoint_unresolved'
    assert record['short_boundary_evidence']['refusal'] == 'SHORT_ENDPOINT_DIRECTION_CHALLENGE_UNSTABLE'


def test_boundary_route_crossing_foreign_branch_still_refuses():
    _, record = probe(3, foreign=True)
    assert record['status'] == 'endpoint_unresolved'
    assert record['short_boundary_evidence']['refusal'] == 'SHORT_ENDPOINT_FOREIGN_BRANCH_CHALLENGE'
