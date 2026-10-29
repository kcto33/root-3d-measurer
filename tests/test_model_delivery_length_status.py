"""Partial values never silently become complete finite-end length passes."""
import sys
from pathlib import Path
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from root_model_delivery_validation import complete_length_state


def test_missing_endpoint_records_block_complete_length():
    state=complete_length_state(dict(candidate_centerline=dict(adjacency=[[1],[0]])))
    assert state['complete_length_status']=='BLOCKED_UNRESOLVED_GEOMETRIC_ENDPOINTS'
    assert state['unresolved_geometric_endpoint_ids']==[0,1] and state['full_length_acceptance']=='BLOCKED'


def test_estimates_do_not_certify_length_accuracy():
    state=complete_length_state(dict(candidate_centerline=dict(adjacency=[[1],[0]]),
        geometric_endpoints=[dict(node_id=i,status='geometric_surface_estimate') for i in [0,1]]))
    assert state['complete_length_status']=='GEOMETRIC_ENDPOINT_ESTIMATES_AVAILABLE'
    assert state['full_length_acceptance']=='NOT_SCORED'


def test_isolated_material_is_retained_as_partial():
    state=complete_length_state(dict(candidate_centerline=dict(adjacency=[[1],[0],[]]),
        geometric_endpoints=[dict(node_id=i,status='geometric_surface_estimate') for i in [0,1]]))
    assert state['isolated_component_node_ids']==[2]
    assert state['complete_length_status']=='PARTIAL_COMPONENTS_WITHOUT_UNIQUE_AXIS'


def test_closed_cycle_requires_no_invented_endpoint():
    state=complete_length_state(dict(candidate_centerline=dict(adjacency=[[1,2],[0,2],[0,1]])))
    assert state['complete_length_status']=='CLOSED_AXES_NO_GEOMETRIC_ENDPOINTS'
    assert state['full_length_acceptance']=='NOT_SCORED'


def test_point_only_or_absent_graph_is_not_a_measurable_axis():
    assert complete_length_state({})['complete_length_status']=='BLOCKED_NO_MEASURABLE_AXIS'
    assert complete_length_state(dict(candidate_centerline=dict(adjacency=[[]])))['full_length_acceptance']=='BLOCKED'


def test_duplicate_or_conflicting_endpoint_records_cannot_certify_length():
    state=complete_length_state(dict(candidate_centerline=dict(adjacency=[[1],[0]]),
        geometric_endpoints=[dict(node_id=i,status='geometric_surface_estimate') for i in [0,0,1]]))
    assert state['unresolved_geometric_endpoint_ids']==[0] and state['full_length_acceptance']=='BLOCKED'
