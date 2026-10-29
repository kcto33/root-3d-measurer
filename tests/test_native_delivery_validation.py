"""Delivery evidence binding and a missing-whole-branch denominator control."""
import json
import os
import sys
from pathlib import Path

import pytest

PROJECT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PROJECT/'src'))
sys.path.insert(0,str(PROJECT/'scripts'))
from root_model_delivery_exports import file_reference, measurement_tables, write_json
from root_model_delivery_validation import register_score, validation_for
from score_native_model_delivery import prepare_contract


def test_missing_whole_branch_does_not_disappear_when_output_validity_is_one():
    pred=dict(candidate_centerline=dict(points=[[0,0,0],[1,0,0]],adjacency=[[1],[0]]),
              radius_curve=dict(sampling='edge_midpoint',rows=[dict(node_ids=[0,1],
                position=[.5,0,0],normal=[1,0,0],radius=.1,available=True,status='VALID')]),
              coordinate_units='model_units',scale_mm_per_model_unit=None,
              provenance=dict(mode='fresh'),entity_binding=dict(stage='regularized_negative_indices'))
    ref=dict(points=[[0,0,0],[1,0,0],[0,1,0],[1,1,0]],adjacency=[[1],[0],[3],[2]],
             edges=[[0,1],[2,3]],construction_segment_radii=[.1,.1],max_radius=.1,pitch=.1)
    contract=prepare_contract(pred,ref)
    covered=sum(b['weight'] for b in contract['bins'] if b['axis_covered'])
    assert contract['reference_total_length']==2
    assert covered/contract['reference_total_length']==pytest.approx(.5)
    assert not contract['topology']['all_paths_matched']
    assert any(p['relative_error'] is None for p in contract['per_path'])
    assert measurement_tables(pred)[2]['V_radius_out']==1


def test_validation_is_bound_to_measurement_and_score_bytes():
    base=Path(os.environ['NATIVE_ENTITY_CONTROL_OUTPUT'])/'score_binding'
    base.mkdir(parents=True,exist_ok=False)
    experiment=next(p for p in base.parents if (p/'EXPERIMENT_LOG.md').is_file())
    prefix=base/'control';measurement=base/'control_measurement.json'
    write_json(measurement,dict(control=True))
    contract=base/'contract.json';write_json(contract,dict(measurement=file_reference(measurement)))
    score=base/'score.json'
    write_json(score,dict(contract=file_reference(contract),precision_pass=False,
                          checks=dict(topology=False),metrics=dict(C_axis_ref=.5)))
    registered=register_score(prefix,score)
    assert registered['precision_status']=='FAILED_FIXED_GATES'
    assert validation_for(prefix)['metrics']['C_axis_ref']==.5
    wrong=base/'wrong_measurement.json';write_json(wrong,dict(other=True))
    with pytest.raises(ValueError,match='DIFFERENT_MEASUREMENT'):
        register_score(base/'wrong',score)
    # Only mutate this test's own tiny score, never an experimental score.
    with score.open('a',encoding='utf-8') as stream:stream.write('\n')
    with pytest.raises(ValueError,match='SCORE_CHANGED'):
        validation_for(prefix)


def test_partial_unadmitted_graph_is_explicitly_blocked_even_with_output_files():
    base=Path(os.environ['NATIVE_ENTITY_CONTROL_OUTPUT'])/'blocked_control'
    base.mkdir(parents=True,exist_ok=False)
    write_json(base/'blocked_measurement.json',dict(status='NATIVE_OUTPUT_WITH_EXPLICIT_GAPS',
               failures=[dict(code='UNSUPPORTED_SOURCE_GRAPH_REPRESENTATION')]))
    status=validation_for(base/'blocked')
    assert status['runtime_status']=='BLOCKED_SOURCE_GEOMETRY_OR_GRAPH'
    assert status['precision_status']=='NOT_SCORED'
