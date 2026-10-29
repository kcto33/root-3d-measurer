"""Only a complete, same-input/topology/signature candidate may replace a failed axis."""
import copy,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from root_native_axis_order import endpoint_candidate_admission

def example():
    return dict(candidate_centerline=dict(points=[[0,0,0],[1,0,0]],adjacency=[[1],[0]]),
        geometric_endpoints=[dict(node_id=0,status='geometric_surface_estimate'),dict(node_id=1,status='geometric_surface_estimate')],
        entity_binding=dict(cells_sha256='immutable_input',stage='regularized_negative_indices',pitch=1.),
        stages=dict(skeleton=dict(thinning=dict(full_source_topology_check=dict(B0_B1_B2_Euler_equal=True)))),failures=[])

def test_fully_qualified_candidate_preserves_complete_graph_signature():
    primary=example();primary['geometric_endpoints'][1]['status']='endpoint_unresolved'
    assert endpoint_candidate_admission(primary,example())['qualified']

def test_fewer_unresolved_is_insufficient_for_replacement():
    a=example();b=example();b['geometric_endpoints'][1]['status']='endpoint_unresolved'
    assert not endpoint_candidate_admission(a,b)['qualified']

def test_lost_branch_or_input_mismatch_is_rejected_even_if_endpoints_are_resolved():
    a=example();b=example();b['entity_binding']['cells_sha256']='different'
    assert not endpoint_candidate_admission(a,b)['qualified']
    b=example();b['candidate_centerline']=dict(points=[[0,0,0],[1,0,0],[0,1,0]],adjacency=[[1,2],[0],[0]])
    b['geometric_endpoints']=[dict(node_id=i,status='geometric_surface_estimate') for i in [1,2]]
    # The three-node chain still has two ports and one path, so add a genuine third arm.
    b['candidate_centerline']['points'].append([0,0,1]);b['candidate_centerline']['adjacency'].append([0]);b['candidate_centerline']['adjacency'][0].append(3)
    b['geometric_endpoints'].append(dict(node_id=3,status='geometric_surface_estimate'))
    assert not endpoint_candidate_admission(a,b)['qualified']

def test_duplicate_missing_records_and_source_topology_failure_cannot_pass():
    a=example();b=example();b['geometric_endpoints']=[b['geometric_endpoints'][0]]*2
    assert not endpoint_candidate_admission(a,b)['qualified']
    b=example();b['stages']['skeleton']['thinning']['full_source_topology_check']['B0_B1_B2_Euler_equal']=False
    assert not endpoint_candidate_admission(a,b)['qualified']
