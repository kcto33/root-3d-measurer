"""Separate, SHA-bound presentation of real reconstruction+measurement errors."""
import csv
import io
import json
from pathlib import Path
from root_model_delivery_exports import file_reference,write_json

PROJECT=Path(__file__).resolve().parents[1];WORKSPACE=PROJECT.parents[1]
SCOPE='RECONSTRUCTION_PLUS_MEASUREMENT_ON_STORED_SOURCE_RENDERER_REPLAY'

def verified(path,ref):
    actual=file_reference(path)
    if any(actual[k]!=ref[k] for k in ['sha256','bytes']):raise ValueError('RECONSTRUCTION_EVIDENCE_CHANGED')

def register_reconstruction_comparison(prefix,result_path):
    prefix=Path(prefix).resolve();result_path=Path(result_path).resolve()
    experiment=next(p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file())
    if not result_path.is_relative_to(experiment):raise ValueError('RECONSTRUCTION_EVIDENCE_OUTSIDE_EXPERIMENT')
    result=json.loads(result_path.read_text(encoding='utf-8'))
    query_ref=result['query_contract'];query_path=(WORKSPACE/query_ref['path']).resolve()
    if not query_path.is_relative_to(experiment):raise ValueError('RECONSTRUCTION_EVIDENCE_OUTSIDE_EXPERIMENT')
    verified(query_path,query_ref)
    query=json.loads(query_path.read_text(encoding='utf-8'))
    measurement=prefix.with_name(prefix.name+'_measurement.json');verified(measurement,query['measurement'])
    if result['scope']!=SCOPE or result.get('meter_precision_certified_by_this_comparison') is not False:
        raise ValueError('RECONSTRUCTION_SCOPE_MUST_BE_SEPARATE')
    target=prefix.with_name(prefix.name+'_reconstruction_comparison.json')
    write_json(target,dict(schema_version='native-reconstruction-binding-v001',
        measurement=file_reference(measurement,prefix.parent),result=file_reference(result_path,prefix.parent)))
    return target

def reconstruction_comparison(prefix):
    prefix=Path(prefix).resolve();binding=prefix.with_name(prefix.name+'_reconstruction_comparison.json')
    if not binding.is_file():return dict(available=False)
    record=json.loads(binding.read_text(encoding='utf-8'))
    experiment=next(p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file())
    measurement=prefix.with_name(prefix.name+'_measurement.json');verified(measurement,record['measurement'])
    path=(prefix.parent/record['result']['path']).resolve()
    if not path.is_relative_to(experiment):raise ValueError('RECONSTRUCTION_EVIDENCE_OUTSIDE_EXPERIMENT')
    verified(path,record['result']);result=json.loads(path.read_text(encoding='utf-8'))
    if result['scope']!=SCOPE or result.get('meter_precision_certified_by_this_comparison') is not False:
        raise ValueError('RECONSTRUCTION_SCOPE_MUST_BE_SEPARATE')
    return dict(available=True,scope=SCOPE,root=result['root'],metrics=result['metrics'],topology=result['topology'],
        historical_render_surface_byte_identity_certified=False,meter_precision_certified_by_this_comparison=False,
        measurement_sha256=record['measurement']['sha256'],comparison_sha256=record['result']['sha256'],
        reference_status_counts=result['reference_status_counts'])

def reconstruction_csv(prefix):
    value=reconstruction_comparison(prefix)
    if not value['available']:raise ValueError('NO_RECONSTRUCTION_COMPARISON')
    row=dict(root=value['root'],scope=value['scope'],**value['metrics'],
        measurement_sha256=value['measurement_sha256'],comparison_sha256=value['comparison_sha256'])
    stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=list(row),lineterminator='\r\n')
    writer.writeheader();writer.writerow(row)
    return stream.getvalue().encode('utf-8-sig')
