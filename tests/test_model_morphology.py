"""Analytic geometry controls, independent formulas and invariance checks."""
import itertools
import json
import math
import numpy as np
import pytest

from root_model_morphology import (box_count, build_morphology, depth_profile,
    fractal_dimension, morphology_csv, settings, solid_hull)
from root_reconstruction_measurement import ReconstructedSolid


def report(points, edges, scale=None):
    adjacency=[[] for _ in points]
    for i,j in edges:
        adjacency[i].append(j);adjacency[j].append(i)
    return dict(status='NATIVE_OUTPUT_WITH_EXPLICIT_GAPS',coordinate_units='model_units',
        scale_mm_per_model_unit=scale,candidate_centerline=dict(points=np.asarray(points).tolist(),adjacency=adjacency),
        geometric_endpoints=[dict(node_id=i,status='geometric_surface_estimate') for i,a in enumerate(adjacency) if len(a)==1])


def values(result):
    return {m['parameter']:m['value'] for m in result['metrics']}


def cuboid():
    return ReconstructedSolid(np.array(list(itertools.product(range(2),range(3),range(4))),dtype=np.int64),.5,{})


@pytest.mark.parametrize('cells',[
    [[0,0,0]], [[0,0,0],[3,0,0]], [[0,0,0],[1,0,0],[0,2,0]],
    list(itertools.product(range(2),range(3),range(4)))])
def test_cube_union_convex_hull_matches_all_original_corners(cells):
    from scipy.spatial import ConvexHull
    p=np.array(cells,dtype=np.int64)
    corners,hull=solid_hull(p)
    all_corners=(p[:,None,:]+np.array(list(itertools.product((0,1),repeat=3)))).reshape(-1,3)
    reference=ConvexHull(all_corners)
    assert hull.volume==pytest.approx(reference.volume)
    assert hull.area==pytest.approx(reference.area)


def test_exact_cuboid_dimensions_density_projection_sphericity():
    r=report([[.5,.5,0],[.5,.5,2]],[[0,1]])
    result=build_morphology(r,solid=cuboid());m=values(result)
    assert m['root_depth']==2
    assert m['horizontal_spread']==pytest.approx(math.hypot(1,1.5))
    assert m['convex_hull_volume']==pytest.approx(3)
    assert m['volume_solidity']==pytest.approx(1)
    assert m['projected_area']==pytest.approx(1.5)
    assert m['projected_solidity']==pytest.approx(1)
    assert m['hull_length_density']==pytest.approx(2/3)
    assert m['inclination_mean']==0 and m['tortuosity_mean']==1
    surface=2*(1*1.5+1*2+1.5*2)
    assert m['surface_volume_ratio']==pytest.approx(surface/3)
    assert m['sphericity']==pytest.approx(math.pi**(1/3)*18**(2/3)/surface)
    assert m['axis_depth_d50']==pytest.approx(1)
    assert m['axis_depth_d95']==pytest.approx(1.9)
    assert sum(row['solid_volume'] for row in result['depth_profile'])==pytest.approx(3)
    assert sum(row['axis_length'] for row in result['depth_profile'])==pytest.approx(2)
    assert result['settings']['vertical_axis']=='z' and result['settings']['down_sign']==-1
    json.dumps(result,allow_nan=False)


def test_physical_scaling_and_orientation_do_not_change_dimensionless_traits():
    r=report([[.5,.5,0],[.5,.5,2]],[[0,1]])
    a=values(build_morphology(r,solid=cuboid()))
    r['scale_mm_per_model_unit']=10
    b=values(build_morphology(r,solid=cuboid()))
    for key,power in [('root_depth',1),('convex_hull_volume',3),('projected_area',2),('hull_length_density',-2),('surface_volume_ratio',-1),('volume_solidity',0)]:
        assert b[key]==pytest.approx(a[key]*10**power)
    changed=build_morphology(r,solid=cuboid(),options=dict(vertical_axis='x',down_sign=1))
    assert values(changed)['root_depth']==10 and values(changed)['inclination_mean']==90


def test_depth_integrates_edges_not_vertex_counts_and_handles_direction():
    pairs=np.array([[[0,0,0],[0,0,10]],[[0,0,9],[3,0,9]]],float)
    rows,d50,d95,outside=depth_profile(pairs,np.array([10.,3.]),10,2,-1,10)
    assert sum(r['axis_length'] for r in rows)==pytest.approx(13)
    assert d50==pytest.approx(3.5) and d95==pytest.approx(9.35)
    assert outside==0
    _,d50,_,_=depth_profile(pairs,np.array([10.,3.]),0,2,1,10)
    assert d50==pytest.approx(6.5)


def test_missing_axis_keeps_solid_traits_and_volume_profile_without_fake_length():
    result=build_morphology(report([],[]),solid=cuboid());m=values(result)
    assert result['partial_axis'] and m['root_depth']==2
    assert m['hull_length_density'] is None and m['fractal_dimension_3d'] is None
    assert m['tortuosity_mean'] is None
    assert sum(r['solid_volume'] for r in result['depth_profile'])==pytest.approx(3)


def test_blocked_axis_never_claims_full_hull_length_density_and_closed_path_excluded():
    r=report([[0,0,0],[1,0,0],[1,1,0]],[[0,1],[1,2],[2,0]])
    result=build_morphology(r,solid=cuboid())
    assert values(result)['tortuosity_mean'] is None
    r['status']='NATIVE_MEASUREMENT_BLOCKED'
    result=build_morphology(r,solid=cuboid())
    assert result['partial_axis'] and values(result)['hull_length_density'] is None


@pytest.mark.parametrize('phase',[(0,0,0),(.5,.5,.5),(0,.5,0)])
def test_box_counts_continuous_line_invariant_to_subdivision_translation_scale(phase):
    a=np.array([[0.,0.,0.],[64.,0.,0.]])
    b=np.array([[i,0,0] for i in range(65)],float)
    for size in (1.,2.,4.,8.):
        expected=65 if size==1 else int(64/size)+1
        assert box_count(a,[[0,1]],size,phase)==expected
        assert box_count(b,[[i,i+1] for i in range(64)],size,phase)==expected
        assert box_count(b*3+[104,-90,5],[[i,i+1] for i in range(64)],size*3,phase)==expected


def test_diagonal_grid_crossings_long_and_short_segments_agree():
    a=np.array([[0.,0.,0.],[16.,16.,16.]])
    b=np.array([[i,i,i] for i in range(17)],float)
    for phase in itertools.product((0.,.5),repeat=3):
        assert box_count(a,[[0,1]],1,phase)==box_count(b,[[i,i+1] for i in range(16)],1,phase)


def test_fractal_uses_all_fixed_scales_reports_phase_and_no_cherry_pick():
    p=np.array([[0.,0.,0.],[64.,0.,0.]])
    fd=fractal_dimension(p,[[0,1]],1,settings())
    assert [s['box_size_model_units'] for s in fd['scales']]==[1,2,4,8,16]
    x=np.log(64/np.array([1,2,4,8,16]));y=np.log([65,33,17,9,5])
    assert fd['dimension']==pytest.approx(np.polyfit(x,y,1)[0])
    assert len(fd['phase_dimensions'])==8
    assert fractal_dimension(p,[[0,1]],8,settings())['dimension'] is None
    assert fractal_dimension(p,[[0,1]],None,settings())['status']=='INPUT_RESOLUTION_REQUIRED'
    with pytest.raises(ValueError,match='BUDGET'):
        box_count(p,[[0,1]],.001,(0,0,0),visit_budget=100)


def test_mesh_original_volume_remains_unqualified():
    from types import SimpleNamespace
    mesh=SimpleNamespace(vertices=np.array(list(itertools.product((0.,1.),repeat=3))))
    result=build_morphology(report([[0,0,0],[0,0,1]],[[0,1]]),mesh=mesh)
    assert values(result)['convex_hull_volume']==pytest.approx(1)
    assert values(result)['volume_solidity'] is None
    assert values(result)['projected_area'] is None


@pytest.mark.parametrize('option',[{'vertical_axis':'a'},{'down_sign':True},{'down_sign':0},{'minimum_scales':2}])
def test_invalid_settings_rejected(option):
    with pytest.raises(ValueError):settings(option)


def test_csv_has_units_status_and_default_orientation():
    result=build_morphology(report([],[]),solid=cuboid())
    assert b'fractal_dimension_3d' in morphology_csv(result)
    assert b'vertical_axis,down_sign' in morphology_csv(result,profile=True)
    import csv,io
    result['measurement']={'sha256':'bound-measurement-sha'}
    row=next(csv.DictReader(io.StringIO(morphology_csv(result).decode('utf-8-sig'))))
    assert row['method']=='root-morphology-v001' and row['measurement_sha256']=='bound-measurement-sha'
    assert row['vertical_axis']=='z' and row['down_sign']=='-1'


def test_depth_outside_bounds_is_separate_from_in_bounds_profile():
    pairs=np.array([[[0,0,-1],[0,0,3]]],float)
    rows,_,_,outside=depth_profile(pairs,np.array([4.]),2,2,-1,2)
    assert sum(r['axis_length'] for r in rows)==pytest.approx(2)
    assert outside==pytest.approx(2)


def test_sierpinski_control_separates_branched_space_filling_from_line():
    # Finite recursive triangular gasket, independent construction and vertices.
    triangles=[np.array([[0.,0.,0.],[256.,0.,0.],[128.,128*math.sqrt(3),0.]])]
    for _ in range(7):
        triangles=[(t+v)/2 for t in triangles for v in t]
    points=[];edges=[]
    for t in triangles:
        offset=len(points);points.extend(t)
        edges.extend([(offset,offset+1),(offset+1,offset+2),(offset+2,offset)])
    fd=fractal_dimension(np.array(points),edges,2,settings())
    assert fd['dimension']==pytest.approx(math.log(3)/math.log(2),abs=.13)
    assert fd['dimension']>1.4
