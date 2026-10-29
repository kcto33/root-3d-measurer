"""Read immutable, SHA-bound independent scores without rewriting predictions."""
import json
from collections import Counter
from pathlib import Path

from root_model_delivery_exports import file_reference, write_json,local_date


def complete_length_state(pred):
    adjacency=(pred.get('candidate_centerline') or {}).get('adjacency',[])
    endpoints={i for i,row in enumerate(adjacency) if len(row)==1}
    isolated=[i for i,row in enumerate(adjacency) if not row]
    records=pred.get('geometric_endpoints',[]);counts=Counter(r['node_id'] for r in records)
    resolved={r['node_id'] for r in records if r['status']=='geometric_surface_estimate' and counts[r['node_id']]==1}
    missing=sorted(endpoints-resolved)
    status=('BLOCKED_NO_MEASURABLE_AXIS' if not any(adjacency) else
            'BLOCKED_UNRESOLVED_GEOMETRIC_ENDPOINTS' if missing else
            'PARTIAL_COMPONENTS_WITHOUT_UNIQUE_AXIS' if isolated else
            'GEOMETRIC_ENDPOINT_ESTIMATES_AVAILABLE' if endpoints else 'CLOSED_AXES_NO_GEOMETRIC_ENDPOINTS')
    return dict(complete_length_status=status,unresolved_geometric_endpoint_ids=missing,
                geometric_endpoint_count=len(endpoints),estimated_geometric_endpoint_count=len(endpoints&resolved),
                isolated_component_node_ids=isolated,full_length_acceptance='BLOCKED' if status.startswith('BLOCKED') else 'NOT_SCORED')


def validation_for(prefix):
    prefix = Path(prefix)
    measurement = prefix.with_name(prefix.name + '_measurement.json')
    pred = json.loads(measurement.read_text(encoding='utf-8'))
    blockers = [f for f in pred.get('failures', []) if f['code'] in
                ['UNSUPPORTED_SOURCE_GRAPH_REPRESENTATION', 'NATIVE_PIPELINE_BLOCKED','MESH_PIPELINE_BLOCKED']]
    runtime = 'BLOCKED_SOURCE_GEOMETRY_OR_GRAPH' if blockers or pred.get('status') in ['NATIVE_MEASUREMENT_BLOCKED','MESH_MEASUREMENT_BLOCKED'] else 'OUTPUT_WITH_EXPLICIT_GAPS'
    path = prefix.with_name(prefix.name + '_validation.json')
    length=complete_length_state(pred)
    if runtime.startswith('BLOCKED'):length['full_length_acceptance']='BLOCKED'
    if not path.is_file():
        return dict(precision_status='NOT_SCORED', metrics=None, runtime_status=runtime, blockers=blockers,**length)
    record = json.loads(path.read_text(encoding='utf-8'))
    if file_reference(measurement)['sha256'] != record['measurement']['sha256']:
        raise ValueError('SCORE_MEASUREMENT_CHANGED')
    score_path = (prefix.parent / record['score']['path']).resolve()
    experiment = next(p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file())
    if not score_path.is_relative_to(experiment):
        raise ValueError('SCORE_OUTSIDE_EXPERIMENT')
    actual = file_reference(score_path)
    if any(actual[k] != record['score'][k] for k in ['sha256', 'bytes']):
        raise ValueError('SCORE_CHANGED')
    checks=record.get('checks',{});required=['reference_comparable','topology','reference_endpoints','predicted_geometric_endpoints','total_length','per_path_length']
    if not runtime.startswith('BLOCKED') and not length['complete_length_status'].startswith('BLOCKED') and not length['isolated_component_node_ids']:
        if all(checks.get(k) is True for k in required):length['full_length_acceptance']='PASSED_FIXED_LENGTH_GATES'
        elif any(checks.get(k) is False for k in required):length['full_length_acceptance']='FAILED_FIXED_LENGTH_GATES'
    return dict(record, runtime_status=runtime, blockers=blockers,**length)


def register_score(prefix, score_path):
    prefix, score_path = Path(prefix).resolve(), Path(score_path).resolve()
    experiment = next(p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file())
    if not score_path.is_relative_to(experiment):
        raise ValueError('SCORE_OUTSIDE_EXPERIMENT')
    measurement = prefix.with_name(prefix.name + '_measurement.json')
    score = json.loads(score_path.read_text(encoding='utf-8'))
    contract_path = Path(score['contract']['path'])
    actual = file_reference(contract_path)
    if any(actual[k] != score['contract'][k] for k in ['sha256', 'bytes']):
        raise ValueError('SCORE_CONTRACT_CHANGED')
    contract = json.loads(contract_path.read_text(encoding='utf-8'))
    bound = file_reference(measurement)
    if any(contract['measurement'][k] != bound[k] for k in ['sha256', 'bytes']):
        raise ValueError('SCORE_BINDS_DIFFERENT_MEASUREMENT')
    passed = score['precision_pass'] is True and all(score['checks'].values())
    record = dict(schema_version='native-delivery-validation-v001',
                  measurement=file_reference(measurement, prefix.parent),
                  score=file_reference(score_path, prefix.parent),
                  precision_status='REFERENCE_NOT_COMPARABLE' if score['checks'].get('reference_comparable') is False else
                  'PASSED_FIXED_GEOMETRY_GATES' if passed else 'FAILED_FIXED_GATES',
                  metrics=score['metrics'], checks=score['checks'],
                  scope='this finite entity and declared reference; real reconstruction accuracy is separate')
    target = prefix.with_name(prefix.name + '_validation.json')
    with (experiment/'EXPERIMENT_LOG.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'+local_date()+' PRELOG register immutable score: '+target.relative_to(experiment).as_posix()+'\n')
    write_json(target, record)
    with (experiment/'EXPERIMENT_LOG.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'+local_date()+' RESULT score sidecar: '+target.relative_to(experiment).as_posix()+' '+record['precision_status']+'\n')
    return record
