"""Downloaded coverage is the bound score; underlying geometry and units stay exact."""
import csv
import hashlib
import io
import json
import os
import sys
from pathlib import Path

import pytest

P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from root_model_delivery_exports import export_tables,write_json,file_reference
from root_model_delivery_validation import register_score
from root_model_delivery_downloads import export_download,display_validation


@pytest.fixture
def result_prefix(request):
    folder=Path(os.environ['NATIVE_ENTITY_CONTROL_OUTPUT'])/hashlib.sha256(request.node.name.encode()).hexdigest()[:8]
    folder.mkdir(parents=True,exist_ok=False);prefix=folder/'fixture'
    report=dict(candidate_centerline=dict(points=[[0,0,0],[1,0,0]],adjacency=[[1],[0]]),
        radius_curve=dict(sampling='edge_midpoint',rows=[dict(node_ids=[0,1],position=[.5,0,0],normal=[1,0,0],radius=.1,available=True,status='VALID')]),
        coordinate_units='model_units',scale_mm_per_model_unit=2.,provenance=dict(mode='fresh'),
        entity_binding=dict(stage='regularized_negative_indices'),failures=[],status='NATIVE_OUTPUT_WITH_EXPLICIT_GAPS')
    measure=folder/'fixture_measurement.json';write_json(measure,report)
    exports,paths,coverage=export_tables(report,prefix)
    write_json(folder/'fixture_delivery.json',dict(measurement=file_reference(measure,folder),exports=exports,summary=dict(coverage=coverage)))
    return prefix


def bind_score(prefix,topology=True):
    measure=prefix.with_name(prefix.name+'_measurement.json');contract=prefix.parent/'contract.json'
    write_json(contract,dict(measurement=file_reference(measure)))
    score=prefix.parent/'score.json';write_json(score,dict(contract=file_reference(contract),
        precision_pass=topology,checks=dict(topology=topology,reference_comparable=True),
        metrics=dict(C_axis_ref=.97,C_radius_ref=.94,V_radius_out=1.)))
    return register_score(prefix,score)


def test_scored_download_has_reference_rates_without_changing_geometry_or_mm(result_prefix):
    prefix=result_prefix;bound=bind_score(prefix);measurement=prefix.with_name(prefix.name+'_measurement.json')
    before=measurement.read_bytes();raw,_,_=export_download(prefix,'paths')
    row=next(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))))
    assert float(row['length'])==1 and float(row['length_mm'])==2
    assert float(row['global_C_axis_ref'])==.97 and float(row['global_C_radius_ref'])==.94
    assert float(row['global_V_radius_out'])==1 and row['score_sha256']==bound['score']['sha256']
    graph=json.loads(export_download(prefix,'graph')[0])
    assert graph['candidate_centerline']==json.loads(before)['candidate_centerline']
    assert graph['coverage']['C_axis_ref']==.97 and measurement.read_bytes()==before


def test_unscored_download_keeps_unknown_reference_blank(result_prefix):
    row=next(csv.DictReader(io.StringIO(export_download(result_prefix,'samples')[0].decode('utf-8-sig'))))
    assert row['global_C_axis_ref']==row['global_C_radius_ref']==''
    assert float(row['global_V_radius_out'])==1 and row['precision_status']=='NOT_SCORED'


def test_failed_path_correspondence_is_not_exported_as_certified_coverage(result_prefix):
    bind_score(result_prefix,topology=False)
    value=display_validation(result_prefix)
    assert value['metrics']['C_axis_ref'] is None and value['metrics']['candidate_C_axis_ref']==.97
    row=next(csv.DictReader(io.StringIO(export_download(result_prefix,'paths')[0].decode('utf-8-sig'))))
    assert row['global_C_axis_ref']==row['global_C_radius_ref']==''


def test_changed_measurement_is_rejected_even_without_a_score(result_prefix):
    file=result_prefix.with_name(result_prefix.name+'_measurement.json')
    with file.open('a',encoding='utf-8') as stream:stream.write('\n')
    with pytest.raises(ValueError,match='SOURCE_MEASUREMENT_CHANGED'):export_download(result_prefix,'graph')
