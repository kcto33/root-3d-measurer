"""Deterministic download view of immutable measurements plus bound validation."""
import csv
import io
import json
from pathlib import Path

from root_model_delivery_exports import file_reference
from root_model_delivery_validation import validation_for


def connected_validation_for(prefix):
    prefix=Path(prefix);sidecar=prefix.with_name(prefix.name+'_connected_section_validation.json')
    if not sidecar.exists():return dict(available=False)
    value=json.loads(sidecar.read_text(encoding='utf-8'))
    for key in ['measurement','score']:
        ref=value[key];actual=file_reference(Path(ref['path']))
        if any(actual[k]!=ref[k] for k in ['sha256','bytes']):raise ValueError('CONNECTED_VALIDATION_BINDING_CHANGED')
    expected=file_reference(prefix.with_name(prefix.name+'_measurement.json'))
    if expected['sha256']!=value['measurement']['sha256']:raise ValueError('CONNECTED_MEASUREMENT_MISMATCH')
    score=json.loads(Path(value['score']['path']).read_text(encoding='utf-8'))
    if value['radius_quantity']!='whole_original_connected_section_area_equivalent' or any(value[k]!=score[k] for k in ['precision_pass','metrics','radius_quantity']):raise ValueError('CONNECTED_SCORE_MISMATCH')
    audit_path=prefix.with_name(prefix.name+'_output_section_validation.json')
    if audit_path.exists():
        audit=json.loads(audit_path.read_text(encoding='utf-8'));ref=audit['audit'];actual=file_reference(ref['path'])
        if actual['sha256']!=ref['sha256'] or audit['measurement']['sha256']!=expected['sha256']:raise ValueError('OUTPUT_SECTION_AUDIT_CHANGED')
        result=json.loads(Path(ref['path']).read_text(encoding='utf-8'))
        if result['measurement']['sha256']!=expected['sha256'] or result['summary']!=audit['summary']:raise ValueError('OUTPUT_SECTION_AUDIT_MISMATCH')
        value['output_section_audit']=audit
    return dict(value,available=True)


def connected_columns(value):
    metrics=value.get('metrics',{})
    return dict(connected_precision_pass=value.get('precision_pass'),
        connected_C_axis_ref=metrics.get('C_axis_ref'),connected_C_radius_ref=metrics.get('C_radius_ref'),
        connected_V_radius_out=metrics.get('V_radius_out'),connected_score_sha256=value.get('score',{}).get('sha256',''))


def parameter_download(prefix,kind='json'):
    """Bound parameter view with current independent validation kept separate."""
    if kind not in ['json','csv']:raise ValueError('UNKNOWN_PARAMETER_EXPORT')
    prefix=Path(prefix);delivery=json.loads(prefix.with_name(prefix.name+'_delivery.json').read_text(encoding='utf-8'))
    measurement=delivery['measurement'];actual=file_reference(prefix.parent/measurement['path'])
    if any(actual[k]!=measurement[k] for k in ['sha256','bytes']):raise ValueError('SOURCE_MEASUREMENT_CHANGED')
    ref=delivery.get('parameters',{}).get('json')
    if ref is None:raise ValueError('PARAMETERS_NOT_PRESENT_IN_HISTORICAL_RESULT')
    path=(prefix.parent/ref['path']).resolve();actual=file_reference(path)
    if any(actual[k]!=ref[k] for k in ['sha256','bytes']):raise ValueError('SOURCE_PARAMETERS_CHANGED')
    value=json.loads(path.read_text(encoding='utf-8'));validation=display_validation(prefix)
    value['external_validation']=validation;value['measurement_sha256']=measurement['sha256']
    value['parameter_source_sha256']=ref['sha256']
    connected=connected_validation_for(prefix);value['connected_section_validation']=connected
    if kind=='json':return (json.dumps(value,ensure_ascii=False,allow_nan=False)+'\n').encode('utf-8'),'application/json','model_parameters.json'
    output=io.StringIO(newline='');fields=list(value['metrics'][0])+['precision_status','complete_length_status','measurement_sha256']+list(connected_columns(connected))
    writer=csv.DictWriter(output,fieldnames=fields,lineterminator='\r\n');writer.writeheader()
    for row in value['metrics']:
        writer.writerow(dict(row,precision_status=validation['precision_status'],
            complete_length_status=validation['complete_length_status'],measurement_sha256=measurement['sha256'],**connected_columns(connected)))
    return output.getvalue().encode('utf-8-sig'),'text/csv; charset=utf-8','model_parameters.csv'


def display_validation(prefix):
    prefix=Path(prefix)
    value=validation_for(prefix)
    metrics=value.get('metrics')
    if metrics is not None and metrics.get('V_radius_out') is not None:
        delivery=json.loads(prefix.with_name(prefix.name+'_delivery.json').read_text(encoding='utf-8'))
        measured=delivery['summary']['coverage'].get('V_radius_out')
        if measured is None or abs(measured-metrics['V_radius_out'])>1e-12:
            raise ValueError('SCORE_OUTPUT_VALIDITY_MISMATCH')
    checks=value.get('checks',{})
    if metrics is not None and (checks.get('topology') is False or checks.get('reference_comparable') is False):
        metrics=dict(metrics)
        for name in ['C_axis_ref','C_radius_ref']:
            if metrics.get(name) is not None:metrics['candidate_'+name]=metrics[name]
            metrics[name]=None
        value=dict(value,metrics=metrics,coverage_correspondence_qualified=False)
    return value


def export_download(prefix,kind):
    prefix=Path(prefix)
    if kind not in ['paths','samples','graph']:raise ValueError('UNKNOWN_NATIVE_EXPORT')
    delivery=json.loads(prefix.with_name(prefix.name+'_delivery.json').read_text(encoding='utf-8'))
    measurement=delivery['measurement'];actual_measurement=file_reference(prefix.parent/measurement['path'])
    if any(actual_measurement[k]!=measurement[k] for k in ['sha256','bytes']):raise ValueError('SOURCE_MEASUREMENT_CHANGED')
    reference=delivery['exports'][kind];path=(prefix.parent/reference['path']).resolve()
    actual=file_reference(path)
    if any(actual[k]!=reference[k] for k in ['sha256','bytes']):raise ValueError('SOURCE_EXPORT_CHANGED')
    validation=display_validation(prefix);connected=connected_validation_for(prefix)
    coverage=dict(delivery['summary']['coverage'])
    for name in ['C_axis_ref','C_radius_ref']:
        if validation.get('metrics') is not None and name in validation['metrics']:coverage[name]=validation['metrics'][name]
    measured=delivery['measurement']['sha256'];score=validation.get('score',{}).get('sha256','')
    if kind=='graph':
        value=json.loads(path.read_text(encoding='utf-8'));value['coverage']=coverage
        value['external_validation']=validation
        value['connected_section_validation']=connected
        value['measurement_sha256']=measured
        raw=(json.dumps(value,ensure_ascii=False,allow_nan=False,separators=(',',':'))+'\n').encode('utf-8')
        return raw,'application/json','measurement_graph.json'
    reader=csv.DictReader(io.StringIO(path.read_text(encoding='utf-8-sig')))
    extra=['global_C_axis_ref','global_C_radius_ref','global_V_radius_out','precision_status','measurement_sha256','score_sha256','complete_length_status','full_length_acceptance']
    extra+=list(connected_columns(connected));fields=list(reader.fieldnames);section_audits={}
    if kind=='samples' and connected.get('output_section_audit'):
        audit=connected['output_section_audit'];checked=json.loads(Path(audit['audit']['path']).read_text(encoding='utf-8'))
        original=json.loads((prefix.parent/measurement['path']).read_text(encoding='utf-8'))
        section_audits={tuple(sorted(original['radius_curve']['rows'][r['station']]['node_ids'])):r for r in checked['rows']}
        extra+=['independent_section_qualified','independent_section_status','independent_section_radius','output_section_audit_sha256']
    output=io.StringIO(newline='');writer=csv.DictWriter(output,fieldnames=fields+[name for name in extra if name not in fields],lineterminator='\r\n')
    writer.writeheader()
    for row in reader:
        row.update(global_C_axis_ref=coverage.get('C_axis_ref'),global_C_radius_ref=coverage.get('C_radius_ref'),
                   global_V_radius_out=coverage.get('V_radius_out'),precision_status=validation['precision_status'],
                   measurement_sha256=measured,score_sha256=score,
                   complete_length_status=validation['complete_length_status'],full_length_acceptance=validation['full_length_acceptance'])
        row.update(connected_columns(connected))
        if section_audits:
            result=section_audits[tuple(sorted([int(row['edge_start']),int(row['edge_end'])]))]
            row.update(independent_section_qualified=result['reference'].get('radius') is not None,
                independent_section_status=result['reference']['status'],independent_section_radius=result['reference'].get('radius'),
                output_section_audit_sha256=connected['output_section_audit']['audit']['sha256'])
        writer.writerow(row)
    return output.getvalue().encode('utf-8-sig'),'text/csv; charset=utf-8','measurement_'+kind+'.csv'
