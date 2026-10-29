"""Analytic root-axis controls: direction, biological semantics and measurement gaps."""
import copy
import json
import math

import numpy as np
import pytest

from root_axis_architecture import build_axis_architecture
from root_axis_traits import axis_trait_settings, build_axis_traits
from root_model_delivery_exports import measurement_tables


def fixture(points, edges, *, base=0, confirmed=True, tips=(), radius=1.0, halfedge=False):
    adjacency=[[] for _ in points]
    for a,b in edges:
        adjacency[a].append(b);adjacency[b].append(a)
    graph=dict(points=points,adjacency=adjacency,coordinate_units='model_units')
    if halfedge:
        rows=[dict(node_id=i,radius=radius(i,p) if callable(radius) else radius,
                   available=True,status='CLOSED_DISK_CONTROL') for i,p in enumerate(points)]
    else:
        rows=[dict(node_ids=[a,b],radius=radius(a,b) if callable(radius) else radius,
                   available=True,status='CLOSED_DISK_CONTROL') for a,b in edges]
    report=dict(candidate_centerline=graph,coordinate_units='model_units',
                radius_curve=dict(sampling='node_halfedge' if halfedge else 'edge_midpoint',rows=rows))
    samples=measurement_tables(report)[1]
    annotations=dict(root_bases=[] if base is None else [dict(node_id=base,root_type='primary')],
                     endpoint_types=[dict(node_id=i,type='root_tip') for i in tips],
                     confirmed_components=[0] if confirmed and base is not None else [])
    architecture=build_axis_architecture(graph,annotations=annotations)
    return graph,architecture,samples


def tree(**kwargs):
    return fixture([[0,0,0],[0,0,2],[0,0,4],[0,0,6],[2,0,2],[4,0,2],[2,0,4]],
                   [(0,1),(1,2),(2,3),(1,4),(4,5),(2,6)],tips=(3,5,6),**kwargs)


def root(result):
    return next(a for a in result['axes'] if a.get('parent_axis_id') is None and a['orientation_status']=='ROOTED_BASE_TO_APEX')


def test_mother_axis_density_denominators_intervals_and_unbranched_zones():
    g,a,s=tree()
    before=copy.deepcopy((g,a,s))
    result=build_axis_traits(g,a,s)
    axis=root(result);branching=axis['branching']
    assert axis['length']==6
    assert branching['lateral_count']['value']==2
    assert branching['lrbd_whole_axis']['value']==pytest.approx(2/6)
    assert branching['lrbd_whole_axis']['denominator']==6
    assert branching['lrbd_branched_zone']['value']==1
    assert branching['lrbd_branched_zone']['denominator']==2
    assert branching['attachment_arc_positions']==[2,4]
    assert branching['interbranch_distances']==[2]
    assert branching['interbranch_distance_statistics']['cv']==0
    assert branching['basal_unbranched_length']['value']==2
    assert branching['apical_unbranched_length']['value']==2
    assert branching['lateral_count']['status']=='ANNOTATION_CONFIRMED'
    assert (g,a,s)==before
    json.dumps(result,allow_nan=False)


def test_no_child_is_zero_density_but_first_and_last_branch_are_undefined():
    g,a,s=fixture([[0,0,0],[0,0,4]],[(0,1)],tips=(1,))
    branch=root(build_axis_traits(g,a,s))['branching']
    assert branch['lateral_count']['value']==0
    assert branch['lrbd_whole_axis']['value']==0
    assert branch['lrbd_branched_zone']['value'] is None
    assert branch['basal_unbranched_length']['value'] is None
    assert branch['apical_unbranched_length']['status']=='NO_LATERAL_ROOTS'
    assert branch['interbranch_distance_statistics']['mean'] is None


def test_one_lateral_has_no_interbranch_interval_or_positive_branching_zone():
    g,a,s=fixture([[0,0,0],[0,0,2],[0,0,4],[2,0,2]],[(0,1),(1,2),(1,3)],tips=(2,3))
    branch=root(build_axis_traits(g,a,s))['branching']
    assert branch['lateral_count']['value']==1
    assert branch['interbranch_distances']==[]
    assert branch['interbranch_distance_statistics']['median'] is None
    assert branch['lrbd_branched_zone']['denominator']==0
    assert branch['lrbd_branched_zone']['value'] is None


def test_multiple_laterals_at_same_attachment_keep_zero_interval_without_infinite_density():
    g,a,s=fixture([[0,0,0],[0,0,2],[0,0,4],[2,0,2],[-2,0,2]],
                  [(0,1),(1,2),(1,3),(1,4)],tips=(2,3,4))
    branch=root(build_axis_traits(g,a,s))['branching']
    assert branch['interbranch_distances']==[0]
    assert branch['interbranch_distance_statistics']['cv'] is None
    assert branch['lrbd_branched_zone']['value'] is None
    json.dumps(branch,allow_nan=False)


@pytest.mark.parametrize('endpoint',['cut_end','unknown'])
def test_truncated_or_unlabelled_tip_only_has_observed_apical_interval(endpoint):
    g,a,s=tree()
    mother=next(axis for axis in a['axes'] if axis['parent_axis_id'] is None)
    mother['apical_endpoint_type']=endpoint
    branch=root(build_axis_traits(g,a,s))['branching']
    assert branch['observed_apical_unbranched_length']['value']==2
    assert branch['apical_unbranched_length']['value'] is None
    assert branch['apical_unbranched_length']['status']=='TRUE_ROOT_TIP_NOT_IDENTIFIED'


def test_unrooted_geometry_never_becomes_zero_biological_branching_or_direction():
    g,a,s=tree(base=None,confirmed=False)
    result=build_axis_traits(g,a,s)
    assert result['summaries']['by_developmental_order']==[]
    for axis in result['axes']:
        assert axis['branching']['lateral_count']['value'] is None
        assert axis['branching']['lateral_count']['unit']=='count'
        assert axis['branching']['lrbd_whole_axis']['unit']=='model_units^-1'
        assert axis['whole_axis_direction']['growth_angle_from_down_deg'] is None
        assert axis['taper']['slope'] is None


def test_closed_component_keeps_geometric_diameter_but_no_taper_or_root_hierarchy():
    g,a,s=fixture([[0,0,0],[2,0,0],[2,2,0]],[(0,1),(1,2),(2,0)],base=None,confirmed=False)
    result=build_axis_traits(g,a,s)
    assert len(result['axes'])==1
    axis=result['axes'][0]
    assert axis['closed']
    assert axis['diameter_distribution']['mean']==2
    assert axis['taper']['status']=='ORIENTED_OPEN_AXIS_REQUIRED'
    assert result['branch_relations']==[]


def test_candidate_orders_and_annotation_confirmations_remain_separate():
    g,a,s=tree(confirmed=False)
    result=build_axis_traits(g,a,s)
    assert all(row['order_status']=='CANDIDATE' for row in result['summaries']['by_developmental_order'])
    assert root(result)['branching']['lateral_count']['status']=='CANDIDATE_AXIS_RELATIONSHIPS'
    assert all(not row['biological_accuracy_validated'] for row in result['axes'])
    assert all(relation['confirmation_status']=='CANDIDATE' for relation in result['branch_relations'])


def test_unresolved_symmetric_junction_does_not_invent_two_biological_laterals():
    g,a,s=fixture([[0,0,0],[0,0,1],[-1,0,2],[1,0,2]],
                  [(0,1),(1,2),(1,3)],confirmed=False,tips=(2,3))
    result=build_axis_traits(g,a,s)
    assert any(axis['confirmation_status']=='UNRESOLVED' for axis in result['axes'])
    assert result['summaries']['by_developmental_order']==[]
    assert result['summaries']['by_root_type']==[]
    assert result['summaries']['by_strahler_edge_order']
    for axis in result['axes']:
        assert axis['branching']['lateral_count']['value'] is None
        assert axis['branching']['lateral_count']['status']=='UNRESOLVED_AXIS_CONTINUATION'
        assert axis['taper']['slope'] is None
        assert axis['whole_axis_direction']['vector'] is None
    for relation in result['branch_relations']:
        assert relation['confirmation_status']=='UNRESOLVED'
        assert relation['local_insertion_angle_deg']['value'] is None


@pytest.mark.parametrize('halfedge',[False,True])
def test_zero_length_id_edges_are_retained_without_blocking_other_geometry(halfedge):
    g,a,s=fixture([[0,0,0],[0,0,0],[0,0,1],[0,0,1]],
                  [(0,1),(1,2),(2,3)],tips=(3,),halfedge=halfedge)
    result=build_axis_traits(g,a,s)
    axis=root(result)
    assert axis['node_ids']==[0,1,2,3]
    assert axis['zero_length_segment_count']==2
    assert axis['length']==1
    assert axis['diameter_distribution']['mean']==2
    assert axis['whole_axis_direction']['growth_angle_from_down_deg']==pytest.approx(180)
    assert sum(row['arc_weight'] for row in axis['diameter_profile'])==1
    json.dumps(result,allow_nan=False)


def test_all_coincident_axis_returns_na_for_length_normalized_traits():
    g,a,s=fixture([[0,0,0],[0,0,0]],[(0,1)],tips=(1,))
    axis=root(build_axis_traits(g,a,s))
    assert axis['length']==0
    assert axis['branching']['lrbd_whole_axis']['value'] is None
    assert axis['whole_axis_direction']['vector'] is None
    assert axis['diameter_distribution']['mean'] is None
    json.dumps(axis,allow_nan=False)


def test_developmental_orders_count_whole_axes_and_strahler_counts_only_ordered_segments():
    g,a,s=tree()
    result=build_axis_traits(g,a,s)
    groups={row['order']:row for row in result['summaries']['by_developmental_order']}
    assert groups[0]['axis_count']==1 and groups[0]['total_length']==6
    assert groups[1]['axis_count']==2 and groups[1]['total_length']==6
    strahler=result['summaries']['by_strahler_edge_order']
    assert sum(row['total_length'] for row in strahler)==12
    assert all('axis_count' not in row for row in strahler)
    assert len(set(root(result)['strahler_orders']))>1


def test_local_insertion_is_not_global_pca_or_vertical_angle():
    g,a,s=fixture([[0,0,0],[0,0,2],[2,0,4],[2,0,2]],
                  [(0,1),(1,2),(1,3)],tips=(2,3))
    result=build_axis_traits(g,a,s,options=dict(tangent_window=1,down_sign=1))
    relation=result['branch_relations'][0]
    assert relation['local_insertion_angle_deg']['value']==pytest.approx(90)
    assert relation['whole_axis_pca_angle_deg']['value']<90
    assert relation['parent_tangent']['effective_window']==1
    assert relation['child_tangent']['effective_window']==1
    assert 0<root(result)['whole_axis_direction']['growth_angle_from_down_deg']<90


def test_growth_angle_is_directed_against_down_with_upward_axis_at_180():
    g,a,s=fixture([[0,0,0],[0,0,4]],[(0,1)],tips=(1,))
    up=build_axis_traits(g,a,s,options=dict(down_sign=-1))
    down=build_axis_traits(g,a,s,options=dict(down_sign=1))
    assert root(up)['whole_axis_direction']['growth_angle_from_down_deg']==pytest.approx(180)
    assert root(down)['whole_axis_direction']['growth_angle_from_down_deg']==pytest.approx(0)


@pytest.mark.parametrize('length,available',[(.75,True),(.25,False)])
def test_local_physical_window_reports_shortening_and_rejects_too_short_axis(length,available):
    g,a,s=fixture([[0,0,0],[0,0,2],[0,0,4],[length,0,2]],
                  [(0,1),(1,2),(1,3)],tips=(2,3))
    result=build_axis_traits(g,a,s,options=dict(tangent_window=1))
    relation=result['branch_relations'][0]
    assert relation['child_tangent']['effective_window']==pytest.approx(length)
    assert relation['child_tangent']['window_shortened']
    assert (relation['local_insertion_angle_deg']['value'] is not None)==available


def test_continuous_pca_is_invariant_to_nonuniform_vertex_subdivision():
    p1=[[0,0,0],[0,0,2],[2,0,4]]
    p2=[[0,0,0],[0,0,.1],[0,0,.2],[0,0,.3],[0,0,1.9],[0,0,2],[2,0,4]]
    g,a,s=fixture(p1,[(0,1),(1,2)],tips=(2,))
    h,b,t=fixture(p2,[(i,i+1) for i in range(len(p2)-1)],tips=(6,))
    first=root(build_axis_traits(g,a,s))['whole_axis_direction']
    second=root(build_axis_traits(h,b,t))['whole_axis_direction']
    assert second['vector']==pytest.approx(first['vector'])
    assert second['eigenvalues']==pytest.approx(first['eigenvalues'])


def test_diameter_statistics_are_length_weighted_and_never_use_connected_union_radius():
    g,a,s=fixture([[0,0,0],[0,0,1],[0,0,4]],[(0,1),(1,2)],tips=(2,),
                  radius=lambda x,y:1 if x==0 else 2)
    for sample in s:
        sample['connected_section_radius']=1000
    axis=root(build_axis_traits(g,a,s))
    assert axis['diameter_distribution']['mean']==pytest.approx(3.5)
    assert axis['diameter_distribution']['median']==4
    s[1]['connected_section_classification']='mixed_branch'
    axis=root(build_axis_traits(g,a,s))
    assert axis['diameter_distribution']['mean']==2
    assert axis['diameter_distribution']['valid_arc_fraction']==pytest.approx(.25)
    assert axis['diameter_profile'][1]['diameter'] is None


def test_offset_caliber_windows_exclude_junction_mixed_sections_and_report_ratio():
    points=[[0,0,z] for z in range(6)]+[[1,0,4],[2,0,4]]
    edges=[(i,i+1) for i in range(5)]+[(4,6),(6,7)]
    g,a,s=fixture(points,edges,tips=(5,7),radius=lambda x,y:2 if y<=5 else .5)
    mixed=next(row for row in s if set(row['edge_nodes'])=={3,4})
    mixed.update(status='MISSING_MIXED_BRANCH_SECTION',radius=1000,available=True)
    result=build_axis_traits(g,a,s,options=dict(diameter_offset=1,diameter_window=1))
    relation=result['branch_relations'][0]
    assert relation['parent_diameter']['arc_interval']==[2,3]
    assert relation['parent_diameter']['value']==4
    assert relation['child_diameter']['arc_interval']==[1,2]
    assert relation['child_diameter']['value']==1
    assert relation['child_parent_diameter_ratio']['value']==pytest.approx(.25)


def test_missing_radius_support_is_not_interpolated_for_window_or_profile():
    g,a,s=tree()
    child=next(axis for axis in a['axes'] if axis['node_ids']==[1,4,5])
    for row in s:
        if set(row['edge_nodes'])=={1,4}:
            row.update(available=False,radius=None,connected_section_radius=999)
    result=build_axis_traits(g,a,s,options=dict(diameter_offset=.25,diameter_window=.5))
    relation=next(row for row in result['branch_relations'] if row['child_axis_id']==child['axis_id'])
    assert relation['child_diameter']['value'] is None
    assert relation['child_parent_diameter_ratio']['value'] is None
    axis=next(row for row in result['axes'] if row['axis_id']==child['axis_id'])
    assert axis['diameter_distribution']['valid_arc_fraction']==.5
    assert not axis['diameter_profile'][0]['available']
    assert not axis['diameter_distribution']['missing_interpolated']


def test_linear_diameter_taper_uses_qualified_stations_and_reports_partial_coverage():
    points=[[0,0,z] for z in range(5)]
    g,a,s=fixture(points,[(i,i+1) for i in range(4)],tips=(4,),
                  radius=lambda i,j:(5-.5*(i+j)/2)/2)
    result=build_axis_traits(g,a,s)
    taper=root(result)['taper']
    assert taper['slope']==pytest.approx(-.5)
    assert taper['intercept']==pytest.approx(5)
    assert taper['taper_rate']==pytest.approx(.5)
    assert taper['r_squared']==pytest.approx(1)
    s[1].update(available=False,radius=None)
    taper=root(build_axis_traits(g,a,s))['taper']
    assert taper['slope']==pytest.approx(-.5)
    assert taper['valid_arc_fraction']==pytest.approx(.75)
    assert taper['status']=='PARTIAL_VALID_ARC'


def test_taper_needs_three_distinct_stations_and_constant_radius_has_no_fake_r_squared():
    g,a,s=fixture([[0,0,0],[0,0,1],[0,0,2]],[(0,1),(1,2)],tips=(2,))
    assert root(build_axis_traits(g,a,s))['taper']['slope'] is None
    g,a,s=fixture([[0,0,z] for z in range(4)],[(0,1),(1,2),(2,3)],tips=(3,))
    taper=root(build_axis_traits(g,a,s))['taper']
    assert taper['slope']==0 and taper['r_squared'] is None


def test_node_halfedge_radius_mapping_follows_reversed_root_direction():
    g,a,s=fixture([[0,0,0],[0,0,2]],[(0,1)],base=1,tips=(0,),halfedge=True,
                  radius=lambda i,p:1 if i==0 else 3)
    axis=root(build_axis_traits(g,a,s))
    assert [r['s'] for r in axis['diameter_profile']]==[0,2]
    assert [r['diameter'] for r in axis['diameter_profile']]==[6,2]
    assert axis['diameter_distribution']['mean']==4
    assert sum(r['arc_weight'] for r in axis['diameter_profile'])==2


def test_scaling_changes_units_lengths_and_density_but_not_angles_or_ratios():
    g,a,s=tree()
    first=build_axis_traits(g,a,s,options=dict(diameter_offset=.25,diameter_window=.5))
    second=build_axis_traits(g,a,s,scale_mm_per_model_unit=10,
                             options=dict(diameter_offset=2.5,diameter_window=5,tangent_window=10))
    x,y=root(first),root(second)
    assert second['unit']=='mm'
    assert y['length']==pytest.approx(x['length']*10)
    assert y['diameter_distribution']['mean']==pytest.approx(x['diameter_distribution']['mean']*10)
    assert y['branching']['lrbd_whole_axis']['value']==pytest.approx(x['branching']['lrbd_whole_axis']['value']/10)
    assert y['whole_axis_direction']['growth_angle_from_down_deg']==pytest.approx(x['whole_axis_direction']['growth_angle_from_down_deg'])
    for p,q in zip(first['branch_relations'],second['branch_relations']):
        assert q['child_parent_diameter_ratio']['value']==pytest.approx(p['child_parent_diameter_ratio']['value'])
        assert q['local_insertion_angle_deg']['value']==pytest.approx(p['local_insertion_angle_deg']['value'])


def test_missing_edges_are_explicit_profile_gaps_and_duplicate_support_is_rejected():
    g,a,s=fixture([[0,0,0],[0,0,1],[0,0,2]],[(0,1),(1,2)],tips=(2,))
    axis=root(build_axis_traits(g,a,s[:1]))
    assert axis['diameter_profile'][1]['status']=='NOT_MEASURED'
    assert axis['diameter_distribution']['valid_arc_fraction']==.5
    with pytest.raises(ValueError,match='OVERLAPPING'):
        build_axis_traits(g,a,s+[s[0]])


@pytest.mark.parametrize('option',[{'tangent_window':0},{'tangent_window':True},
    {'diameter_offset':-1},{'diameter_window':float('nan')},
    {'minimum_radius_coverage':1.1},{'minimum_window_fraction':0},
    {'vertical_axis':'a'},{'down_sign':True},{'unknown':1}])
def test_invalid_options_are_rejected(option):
    with pytest.raises(ValueError):
        axis_trait_settings(option)


def test_empty_graph_is_serializable_without_fabricated_root_statistics():
    result=build_axis_traits(dict(points=[],adjacency=[]),dict(axes=[]),[])
    assert result['axes']==[] and result['branch_relations']==[]
    assert all(rows==[] for rows in result['summaries'].values())
    json.dumps(result,allow_nan=False)
