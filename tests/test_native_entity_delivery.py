"""Native package integrity, denominator distinctions and unchanged finite geometry."""
import hashlib
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT/'src'))
from root_model_delivery_exports import measurement_tables, native_surface_buffers, export_tables
from root_native_entity_measurement import load_entity, inherit_bound_measurement, axis_aligned_section, NativeSectionIndex,geometric_endpoints,subdivide_measurement_graph,_tangent
from root_reconstruction_measurement import ReconstructedSolid


@pytest.fixture
def case_dir(request):
    base=Path(os.environ['NATIVE_ENTITY_CONTROL_OUTPUT'])
    path=base/hashlib.sha256(request.node.name.encode()).hexdigest()[:8]
    path.mkdir(parents=True,exist_ok=False)
    with (path/'case.json').open('x',encoding='utf-8') as stream:
        json.dump(dict(test=request.node.name),stream)
    return path


def report(points,adjacency,rows,sampling='edge_midpoint'):
    return dict(candidate_centerline=dict(points=points,adjacency=adjacency),
                radius_curve=dict(rows=rows,sampling=sampling),
                coordinate_units='model_units',scale_mm_per_model_unit=None,
                provenance=dict(mode='fresh'),entity_binding=dict(stage='regularized_negative_indices'))


def test_nonuniform_edges_have_arc_weight_not_station_fraction():
    r=report([[0,0,0],[1,0,0],[4,0,0]],[[1],[0,2],[1]],
             [dict(node_ids=[0,1],radius=.2,available=True,status='VALID'),
              dict(node_ids=[1,2],radius=None,available=False,status='MISSING')])
    paths,samples,coverage=measurement_tables(r)
    assert len(paths)==1 and coverage['output_total_length']==4
    assert coverage['V_radius_out']==.25
    assert coverage['available_sample_count']/coverage['sample_count']==.5
    assert coverage['C_axis_ref'] is None and coverage['C_radius_ref'] is None
    assert samples[1]['radius'] is None


def test_junction_and_halfedge_weights_do_not_repeat_shared_edges():
    r=report([[0,0,0],[1,0,0],[0,2,0],[0,0,3]],[[1,2,3],[0],[0],[0]],
             [dict(node_id=i,radius=None if i==0 else .1,available=i!=0,status='S')
              for i in range(4)],'node_halfedge')
    paths,samples,c=measurement_tables(r)
    assert len(paths)==3 and sum(p['length'] for p in paths)==6
    assert sum(s['arc_weight'] for s in samples)==6
    assert c['V_radius_out']==.5


def test_cycle_is_closed_and_each_edge_counted_once():
    r=report([[0,0,0],[1,0,0],[0,1,0]],[[1,2],[0,2],[0,1]],[])
    paths,_,c=measurement_tables(r)
    assert len(paths)==1 and paths[0]['closed']
    assert math.isclose(c['output_total_length'],2+math.sqrt(2))
    assert c['V_radius_out']==0


def test_cube_display_has_all_six_faces_and_does_not_change_source(case_dir):
    s=ReconstructedSolid(np.array([[0,0,0]],np.int64),1.,{})
    before=s.cells.tobytes()
    folder,data=native_surface_buffers(s,case_dir)
    assert data['quads']==6 and data['triangles']==12
    positions=np.fromfile(folder/data['buffers']['positions']['path'],dtype='<f4').reshape(-1,3)
    assert np.array_equal(positions.min(axis=0),[0,0,0])
    assert np.array_equal(positions.max(axis=0),[1,1,1])
    assert s.cells.tobytes()==before and not s.cells.flags.writeable
    assert native_surface_buffers(s,case_dir)[1]==data


def write_package(folder):
    crop=folder/'entity_crop0000.npz'
    np.savez_compressed(crop,source_negative_indices=np.array([[0,0,0]],np.int64),
                        regularized_negative_indices=np.array([[0,0,0],[1,0,0]],np.int64),
                        pitch=np.array(.25),origin=np.array([.125,.125,.125]),
                        shape=np.array([3,2,2],np.int64))
    sha=hashlib.sha256(crop.read_bytes()).hexdigest()
    meta=dict(fields=[dict(cache=dict(path='/remote/project/entity_crop0000.npz',
                                     sha256=sha,bytes=crop.stat().st_size))],
              duplicate_occupied_cells=False,all_source_negative_cells_accounted=1)
    path=folder/'entity_prediction.json'
    with path.open('x',encoding='utf-8') as f:json.dump(meta,f)
    return path,crop


def test_native_package_binds_stage_coordinates_and_original_hash(case_dir):
    metadata,crop=write_package(case_dir)
    before=hashlib.sha256(crop.read_bytes()).hexdigest()
    solid,bound,_=load_entity(metadata)
    assert len(solid.cells)==2 and solid.pitch==.25
    assert bound['stage']=='regularized_negative_indices'
    assert bound['source_cell_count']==1 and not bound['ply_revoxelization_performed']
    source,_,_=load_entity(metadata,stage='source_negative_indices')
    assert len(source.cells)==1
    assert hashlib.sha256(crop.read_bytes()).hexdigest()==before


def test_wrong_native_crop_hash_is_rejected(case_dir):
    metadata,crop=write_package(case_dir)
    with crop.open('ab') as f:f.write(b'corrupt')
    with pytest.raises(ValueError,match='SHA_BYTES_MISMATCH'):
        load_entity(metadata)


def test_inherited_regularized_values_cannot_be_used_for_raw_stage(case_dir):
    metadata,_=write_package(case_dir)
    _,bound,_=load_entity(metadata,stage='source_negative_indices')
    witness=case_dir/'witness.json'
    with witness.open('x',encoding='utf-8') as f:
        json.dump(dict(case=dict(source_frame={})),f)
    with pytest.raises(ValueError,match='MATCHING_REGULARIZED_ENTITY_STAGE'):
        inherit_bound_measurement(witness,bound)


def test_csv_missing_values_and_scale_are_explicit(case_dir):
    r=report([[0,0,0],[2,0,0]],[[1],[0]],
             [dict(node_ids=[0,1],radius=None,available=False,status='MISSING')])
    r['scale_mm_per_model_unit']=10
    refs,paths,c=export_tables(r,case_dir/'measurement')
    import csv
    with (case_dir/refs['samples']['path']).open(encoding='utf-8-sig',newline='') as f:
        row=next(csv.DictReader(f))
    assert row['radius']=='' and row['radius_mm']=='' and row['status']=='MISSING'
    assert paths[0]['length_mm']==20 and c['C_radius_ref'] is None


def test_grid_face_section_deduplicates_squares_without_moving_plane():
    solid=ReconstructedSolid(np.array([[0,0,0],[1,0,0],[1,1,0]],np.int64),.25,{})
    at_face=axis_aligned_section(solid,[.25,.125,.125],[1,0,0])
    just_before=axis_aligned_section(solid,[np.nextafter(.25,0),.125,.125],[1,0,0])
    assert at_face['area']==2*.25**2
    assert just_before['area']==.25**2
    assert math.isclose(at_face['radius'],math.sqrt(.125/math.pi))
    assert not at_face['plane_was_moved']


def test_axis_section_does_not_fill_holes_or_admit_corner_contact():
    cells=np.array([[0,y,z] for y in range(3) for z in range(3) if (y,z)!=(1,1)],np.int64)
    result=axis_aligned_section(ReconstructedSolid(cells,1.,{}),[.5,.5,.5],[1,0,0])
    assert result['radius'] is None and result['status']=='NONDISK_SECTION_HOLES_NOT_FILLED'
    cells=np.array([[0,0,0],[0,1,1]],np.int64)
    result=axis_aligned_section(ReconstructedSolid(cells,1.,{}),[.5,.5,.5],[1,0,0])
    assert result['radius'] is None and result['status']=='NONMANIFOLD_OR_CONTACT_SECTION'


def test_axis_section_requires_strict_query_interior_and_preserves_oblique_route():
    solid=ReconstructedSolid(np.array([[0,0,0]],np.int64),1.,{})
    result=axis_aligned_section(solid,[.5,0,.5],[1,0,0])
    assert result['radius'] is None and result['status']=='QUERY_ON_SECTION_OUTER_BOUNDARY'
    assert axis_aligned_section(solid,[.5,.5,.5],[1,.1,0]) is None


def test_query_index_matches_full_kernel_and_preserves_far_coordinate_tolerance():
    from root_cubical_disk_section import interior_boundary_section_class
    local=[[x,y,z] for x in range(-2,3) for y in range(-2,3) for z in range(-4,5)]
    cells=np.array(local+[[1000,1000,1000]],np.int64)
    solid=ReconstructedSolid(cells,.125,{})
    kind=interior_boundary_section_class();full=kind(solid);indexed=NativeSectionIndex(solid,kind)
    before=solid.cells.tobytes()
    for normal in [[1,.2,.1],[.2,1,.3],[.3,.2,1],[1,1,1]]:
        point=[.0625,.0625,.0625]
        a=full.section(point,normal);b=indexed.section(point,normal)
        assert a['status']==b['status']
        assert a['radius']==pytest.approx(b['radius'],rel=1e-13) if a['radius'] is not None else b['radius'] is None
        assert a['area']==pytest.approx(b['area'],rel=1e-13)
        if a.get('section_centroid_evidence'):
            assert np.allclose(a['section_centroid_evidence']['centroid'],
                               b['section_centroid_evidence']['centroid'],rtol=0,atol=1e-13)
    assert solid.cells.tobytes()==before
    assert indexed.last_stats['selected_cells']<len(solid.cells)


def test_query_index_retains_hole_and_contact_refusals():
    from root_cubical_disk_section import interior_boundary_section_class
    kind=interior_boundary_section_class()
    shapes=[[[x,y,z] for x in range(-3,4) for y in range(3) for z in range(3) if (y,z)!=(1,1)],
            [[0,0,0],[0,1,1]]]
    for cells in shapes:
        solid=ReconstructedSolid(np.array(cells,np.int64),1.,{})
        a=kind(solid).section([.5,.5,.5],[1,.01,.02])
        b=NativeSectionIndex(solid,kind).section([.5,.5,.5],[1,.01,.02])
        assert a['status']==b['status'] and a['radius']==b['radius']


def test_finite_geometric_endpoints_recover_flat_model_extent_without_tip_claim():
    solid=ReconstructedSolid(np.array([[x,0,0] for x in range(10)],np.int64),.1,{})
    source=dict(points=[[.05,.05,.05],[.45,.05,.05],[.95,.05,.05]],adjacency=[[1],[0,2],[1]])
    before=solid.cells.tobytes()
    graph,records=geometric_endpoints(solid,source)
    assert np.allclose(graph['points'][0],[0,.05,.05],atol=1e-14)
    assert np.allclose(graph['points'][2],[1,.05,.05],atol=1e-14)
    assert [i for i,ns in enumerate(graph['adjacency']) if len(ns)==1]==[0,2]
    assert all(r['status']=='geometric_surface_estimate' for r in records)
    assert all(not r['physical_tip_qualified'] for r in records)
    assert solid.cells.tobytes()==before and source['points'][0]==[.05,.05,.05]


def test_tangent_roundoff_rule_is_defined_before_section_outcomes():
    path=np.array([[0,1e-18,0],[.1,-1e-18,0],[.2,1e-18,0],[.3,0,0]])
    normal=_tangent(path,np.array([0,.1,.2,.3]),.15,.01,False)
    assert np.array_equal(normal,[1,0,0])


def test_subdivision_preserves_unique_length_and_cycle_topology():
    graph=dict(points=[[0,0,0],[1,0,0],[0,1,0]],adjacency=[[1,2],[0,2],[0,1]])
    result=subdivide_measurement_graph(graph,.1)
    from root_graph_measurement_ledger import graph_measurement_ledger
    ledger=graph_measurement_ledger(result['points'],result['adjacency'],coordinate_units='model_units')
    assert ledger['source_cycle_rank']==1
    assert ledger['measurement_graph_length']==pytest.approx(2+math.sqrt(2))
    assert all(math.dist(result['points'][a],result['points'][b])<=.1+1e-12 for a,b in ledger['source_edge_ids'])


def test_reference_query_selection_is_independent_of_radius_success():
    import copy
    sys.path.insert(0,str(PROJECT/'scripts'))
    from score_native_model_delivery import prepare_contract
    pred=report([[0,0,0],[1,0,0]],[[1],[0]],
                [dict(node_ids=[0,1],position=[.5,0,0],normal=[1,0,0],radius=.1,available=True,status='VALID')])
    ref=dict(points=[[0,0,0],[1,0,0]],adjacency=[[1],[0]],edges=[[0,1]],
             construction_segment_radii=[.1],max_radius=.1,pitch=.1)
    before=prepare_contract(pred,ref)
    changed=copy.deepcopy(pred);changed['radius_curve']['rows'][0].update(radius=None,available=False,status='MISSING')
    assert prepare_contract(changed,ref)==before
