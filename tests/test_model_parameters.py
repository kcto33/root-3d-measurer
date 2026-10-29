"""Complete parameter output must retain units, missing arcs and blockers."""
import csv,io,math,sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from root_model_parameters import build_parameters,parameter_csv,radius_statistics
from root_reconstruction_measurement import ReconstructedSolid

def report(scale=None,missing=False):
    return dict(status='NATIVE_OUTPUT_WITH_EXPLICIT_GAPS',coordinate_units='model_units',scale_mm_per_model_unit=scale,
        case_label='contract_control',candidate_centerline=dict(points=[[0,0,0],[1,0,0],[4,0,0]],adjacency=[[1],[0,2],[1]]),
        radius_curve=dict(sampling='edge_midpoint',rows=[
            dict(node_ids=[0,1],position=[.5,0,0],normal=[1,0,0],radius=1.,available=True,status='VALID'),
            dict(node_ids=[1,2],position=[2.5,0,0],normal=[1,0,0],radius=None if missing else 3.,available=not missing,status='MISSING' if missing else 'VALID')]),
        geometric_endpoints=[dict(node_id=i,status='geometric_surface_estimate') for i in [0,2]],
        model_bounds=dict(minimum=[0,0,0],maximum=[1,.5,.5]))

def metrics(value):return {r['parameter']:r for r in value['metrics']}

def test_arc_weighted_radius_is_not_unweighted_sample_mean():
    p=build_parameters(report());m=metrics(p)
    assert m['model_length']['value']==4
    assert m['radius_mean']['value']==2.5 and m['radius_median']['value']==3
    assert m['diameter_mean']['value']==5
    assert p['paths'][0]['radius_mean']==2.5

def test_missing_radius_stays_in_output_denominator_without_imputation():
    p=build_parameters(report(missing=True));m=metrics(p)
    assert m['radius_mean']['value']==1 and m['radius_mean']['covered_arc_fraction']==.25
    assert p['radius_model_units']['missing_output_arc']==3
    assert m['radius_mean']['status']=='CONDITIONAL_VALID_ARC'

def test_native_volume_surface_and_mm_powers_are_exact_for_two_cubes():
    solid=ReconstructedSolid(np.array([[0,0,0],[1,0,0]],np.int64),.5,{})
    p=build_parameters(report(scale=10),solid=solid);m=metrics(p)
    assert p['geometry']['volume']==.25 and p['geometry']['surface_area']==2.5
    assert m['volume']['value']==250 and m['volume']['unit']=='mm^3'
    assert m['surface_area']['value']==250 and m['surface_area']['unit']=='mm^2'
    assert m['radius_mean']['value']==25 and m['bbox_x']['value']==10

def test_unresolved_endpoint_keeps_full_length_blank_in_json_and_csv():
    r=report();r['geometric_endpoints'][0]['status']='endpoint_unresolved'
    p=build_parameters(r);m=metrics(p)
    assert m['model_length']['value'] is None and m['model_length']['status']=='BLOCKED'
    assert m['extracted_axis_length']['value']==4
    values={row['parameter']:row for row in csv.DictReader(io.StringIO(parameter_csv(p).decode('utf-8-sig')))}
    assert values['model_length']['value']=='' and values['extracted_axis_length']['value']=='4.0'

def test_empty_valid_radius_is_null_not_zero():
    r=report(missing=True);r['radius_curve']['rows'][0].update(radius=None,available=False)
    p=build_parameters(r)
    assert metrics(p)['radius_mean']['value'] is None
    assert p['radius_model_units']['valid_arc_fraction']==0

def test_isolated_component_prevents_complete_length_claim():
    r=report();r['candidate_centerline']['points'].append([10,10,10]);r['candidate_centerline']['adjacency'].append([])
    p=build_parameters(r)
    assert metrics(p)['model_length']['value'] is None
    assert p['length_state']['complete_length_status']=='PARTIAL_COMPONENTS_WITHOUT_UNIQUE_AXIS'

def test_weighted_statistics_reject_nonfinite_measured_radius():
    with pytest.raises(ValueError):radius_statistics([dict(arc_weight=1,available=True,radius=float('nan'))])
