"""A ring is never cut; bridge-local traits remain available and scoped."""
import copy
import math
import numpy as np
import pytest
from root_axis_architecture import build_axis_architecture
from root_axis_cycle_regions import graph_bridges
from root_axis_traits import build_axis_traits


def graph():
    # Two branched trees connected through a real four-edge square ring.
    p=[[0,0,-4],[0,0,-2],[0,0,0],[2,0,0],[2,0,2],[0,0,2],
       [0,0,4],[0,0,6],[2,0,-2],[2,0,4]]
    edges=[(0,1),(1,2),(1,8),(2,3),(3,4),(4,5),(5,2),(5,6),(6,7),(6,9)]
    a=[[] for _ in p]
    for i,j in edges:a[i].append(j);a[j].append(i)
    return dict(points=p,adjacency=a)


def radii(g):
    return [dict(edge_nodes=[i,j],node_id=None,sampling='edge_midpoint',radius=.1,
                 available=True,status='CLOSED_DISK_CONTROL',arc_weight=math.dist(g['points'][i],g['points'][j]))
            for i,ns in enumerate(g['adjacency']) for j in ns if i<j]


def test_low_link_bridges_match_edge_removal_connectivity():
    g=graph();a=g['adjacency'];found=graph_bridges(a)
    expected=set()
    for i,ns in enumerate(a):
        for j in ns:
            if i>=j:continue
            seen={i};stack=[i]
            while stack:
                node=stack.pop()
                for n in a[node]:
                    if {node,n}=={i,j}:continue
                    if n not in seen:seen.add(n);stack.append(n)
            if j not in seen:expected.add((i,j))
    assert found==expected and len(found)==6


def test_deep_chain_has_no_recursive_stack_limit():
    a=[[] for _ in range(3000)]
    for i in range(2999):a[i].append(i+1);a[i+1].append(i)
    assert len(graph_bridges(a))==2999


def test_ring_edges_retained_once_and_local_traits_recovered():
    g=graph();before=copy.deepcopy(g)
    architecture=build_axis_architecture(g,annotations=dict(root_bases=[dict(node_id=0)],
        junction_continuations=[dict(node_id=1,incoming_node_id=0,outgoing_node_id=2),dict(node_id=6,incoming_node_id=5,outgoing_node_id=7)]))
    assert architecture['status']=='PARTIAL_ROOTED_ARCHITECTURE'
    assert architecture['preserved_cycle_rank']==1
    edges=[tuple(e) for a in architecture['axes'] for e in a['edge_ids']]
    assert len(edges)==len(set(edges))==10
    assert sum(a['length_model_units'] for a in architecture['axes'])==pytest.approx(20)
    assert g==before
    cores=[r for r in architecture['analysis_regions'] if r['kind']=='CYCLIC_CORE']
    assert len(cores)==1 and len(cores[0]['edge_ids'])==4
    downstream=[a for a in architecture['axes'] if a.get('analysis_region_id') and 7 in a['node_ids']][0]
    assert downstream['candidate_developmental_order'] is None
    assert downstream['local_candidate_developmental_order']==0
    assert all(a.get('candidate_strahler_orders') is None for a in architecture['axes'])
    traits=build_axis_traits(g,architecture,radii(g),options=dict(tangent_window=1,diameter_offset=.1,diameter_window=.5))
    assert len(traits['branch_relations'])==2
    assert all(r['local_insertion_angle_deg']['value']==pytest.approx(90) for r in traits['branch_relations'])
    trunk=[a for a in traits['axes'] if a.get('cyclic_boundary_node_ids') and a['orientation_status']=='ROOTED_BASE_TO_APEX'][0]
    assert trunk['branching']['lrbd_whole_axis']['value'] is None
    assert trunk['branching']['lrbd_observed_fragment']['value'] is not None


def test_annotation_replay_and_root_direction_are_deterministic():
    g=graph();s=radii(g);settings=dict(tangent_window=1,diameter_offset=.1,diameter_window=.5)
    first=build_axis_architecture(g,annotations=dict(root_bases=[dict(node_id=0)]),radius_samples=s,inference_options=settings)
    assert build_axis_architecture(g,annotations=first['annotations'],radius_samples=s,inference_options=settings)==first
    reverse=build_axis_architecture(g,annotations=dict(root_bases=[dict(node_id=7)]),radius_samples=s,inference_options=settings)
    assert {tuple(e) for a in reverse['axes'] for e in a['edge_ids']}=={tuple(e) for a in first['axes'] for e in a['edge_ids']}
    upstream=[r for r in reverse['analysis_regions'] if r['entry_is_explicit_root_base']]
    assert upstream[0]['entry_node_id']==7


def test_no_base_no_speculative_orientation_and_global_confirmation_rejected():
    g=graph()
    assert build_axis_architecture(g)['status']=='GEOMETRY_ONLY'
    with pytest.raises(ValueError,match='CYCLIC_COMPONENT_GLOBAL_CONFIRMATION'):
        build_axis_architecture(g,annotations=dict(root_bases=[dict(node_id=0)],confirmed_components=[0]))


def test_pure_ring_stays_unoriented():
    g=graph();ids=[2,3,4,5];p=[g['points'][i] for i in ids];a=[[1,3],[0,2],[1,3],[0,2]]
    result=build_axis_architecture(dict(points=p,adjacency=a),annotations=dict(root_bases=[dict(node_id=0)]))
    assert result['status']=='GEOMETRY_ONLY' and result['preserved_cycle_rank']==1


def test_unaffected_component_confirmation_and_types_still_work():
    g=graph();g['points'] += [[10,0,0],[10,0,2]];g['adjacency'] += [[11],[10]]
    annotations=dict(root_bases=[dict(node_id=0),dict(node_id=10)],confirmed_components=[1])
    first=build_axis_architecture(g,annotations=annotations)
    other=[a for a in first['axes'] if a['component_id']==1][0]
    assert other['confirmation_status']=='CONFIRMED' and other['developmental_order']==0
    annotations=dict(first['annotations'],root_types=[dict(axis_id=other['axis_id'],type='known_type')],confirmed_axes=[other['axis_id']])
    result=build_axis_architecture(g,annotations=annotations)
    other=[a for a in result['axes'] if a['component_id']==1][0]
    assert other['root_type']=='known_type' and other['confirmation_status']=='CONFIRMED'


def test_cross_core_override_and_fragment_confirmation_cannot_fabricate_hierarchy():
    g=graph()
    with pytest.raises(ValueError,match='CROSSES_UNRESOLVED_CYCLIC_CORE'):
        build_axis_architecture(g,annotations=dict(root_bases=[dict(node_id=0)],
            junction_continuations=[dict(node_id=2,incoming_node_id=1,outgoing_node_id=3)]))
    first=build_axis_architecture(g,annotations=dict(root_bases=[dict(node_id=0)]))
    axis=next(a for a in first['axes'] if a.get('axis_scope')=='ROOTED_BRIDGE_FRAGMENT')
    with pytest.raises(ValueError,match='CYCLIC_FRAGMENT_GLOBAL_CONFIRMATION'):
        build_axis_architecture(g,annotations=dict(first['annotations'],confirmed_axes=[axis['axis_id']]))


def test_bridge_regions_have_identical_partition_under_uniform_scale():
    g=graph();a=build_axis_architecture(g,annotations=dict(root_bases=[dict(node_id=0)]))
    enlarged=dict(g,points=(np.array(g['points'])*13.7).tolist())
    b=build_axis_architecture(enlarged,annotations=dict(root_bases=[dict(node_id=0)]))
    assert [(r['kind'],r['edge_ids'],r['entry_node_id']) for r in a['analysis_regions']]==[(r['kind'],r['edge_ids'],r['entry_node_id']) for r in b['analysis_regions']]
    assert sum(x['length_model_units'] for x in b['axes'])==pytest.approx(sum(x['length_model_units'] for x in a['axes'])*13.7)


@pytest.mark.parametrize('points',[None,1,[[0,0]],[[0,0,float('nan')]],'invalid'])
def test_invalid_point_shapes_use_existing_public_validation(points):
    with pytest.raises(ValueError,match='FINITE_N_BY_3_POINTS_REQUIRED'):
        build_axis_architecture(dict(points=points,adjacency=[]))
