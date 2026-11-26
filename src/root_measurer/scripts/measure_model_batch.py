"""Sequential mixed-input batches with per-sample settings and complete denominators."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / 'src'))
from root_model_delivery_exports import file_reference, write_json, local_date
from root_measurement_inputs import DEFAULTS, measurement_command


def _entry(raw, base, defaults, index):
    item = {'path': raw} if isinstance(raw, (str, Path)) else dict(raw)
    path = item.pop('path')
    label = str(item.pop('id', f'case{index:03d}'))
    if not isinstance(path, (str, Path)) or not str(path).strip():
        raise ValueError('BATCH_INPUT_PATH_REQUIRED')
    config = dict(defaults, **item)
    if config.get('voxel_metadata') is not None:
        config['voxel_metadata'] = str((base / config['voxel_metadata']).resolve())
    return label, (base / path).resolve(), config


def _cell(value):
    if isinstance(value, (list, dict)):
        return json.dumps(value, ensure_ascii=False, default=str)
    return "'" + value if isinstance(value, str) and value.lstrip().startswith(('=', '+', '-', '@')) else value


def run(inputs, folder, scale=None, max_memory_mb=8192, cancel_file=None, vertical_axis='z', down_sign=-1,
        section_workers=1, *, defaults=None, manifest_base=None):
    if not isinstance(inputs, list) or not 1 <= len(inputs) <= 1000:
        raise ValueError('BATCH_SIZE_1_TO_1000_REQUIRED')
    global_options = dict(DEFAULTS, scale_mm=scale, max_memory_mb=max_memory_mb,
                          vertical_axis=vertical_axis, down_sign=down_sign, section_workers=section_workers)
    global_options.update(defaults or {})
    if set(global_options) - DEFAULTS.keys():
        raise ValueError('UNKNOWN_BATCH_DEFAULT_OPTION')
    base = Path(manifest_base or Path.cwd()).resolve()
    folder = Path(folder).resolve()
    if folder.parent.name != 'model_delivery_local' or not (folder.parents[2] / 'ROOT_MEASURER_WORKSPACE.json').is_file():
        raise ValueError('EXPERIMENT_MODEL_DELIVERY_RUN_REQUIRED')
    if folder.exists() and any(folder.iterdir()):
        raise ValueError('BATCH_OUTPUT_ALREADY_EXISTS')
    folder.mkdir(exist_ok=True)
    experiment = folder.parents[2]
    started = time.monotonic()
    batch_id = hashlib.sha256(folder.relative_to(experiment).as_posix().encode()).hexdigest()[:16]
    write_json(folder / 'batch_request.json', dict(schema='root-measurer-batch-v002', id=batch_id,
               total=len(inputs), inputs=inputs, defaults=global_options, manifest_base=str(base), started_utc=time.time()))
    with (experiment / 'EXPERIMENT_LOG.md').open('a', encoding='utf-8') as stream:
        stream.write('\n' + local_date() + ' PRELOG mixed batch ' + folder.relative_to(experiment).as_posix() +
                     '; per-sample settings and complete parameter/status rows retained.\n')
    records, parameter_rows = [], []
    for index, raw in enumerate(inputs):
        prefix = folder / f'case{index:03d}'
        row = dict(index=index, sample_id=f'case{index:03d}', input=str(raw),
                   prefix=prefix.relative_to(experiment).as_posix(), exit_code=None, complete_length=False,
                   radius_valid_arc_fraction=None, status='FAILED', run_status='FAILED', error=None)
        event = dict(index=index, total=len(inputs), phase='measurement')
        with (folder / 'batch_progress.jsonl').open('a', encoding='utf-8') as stream:
            stream.write(json.dumps(event) + '\n')
        print(json.dumps(event), flush=True)
        if cancel_file is not None and Path(cancel_file).exists():
            row.update(status='CANCELLED_BEFORE_START', run_status='CANCELLED')
        else:
            try:
                label, source, config = _entry(raw, base, global_options, index)
                row.update(sample_id=label, input=str(source), requested_options=config)
                command, config, source = measurement_command(source, prefix, config, cancel_file)
                row['input_kind'] = 'volume' if '--volume' in command else ('mesh' if '--mesh' in command else 'native_entity')
                with prefix.with_name(prefix.name + '_launcher.log').open('x', encoding='utf-8') as log:
                    completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                        env=dict(os.environ, PYTHONUTF8='1', PYTHONDONTWRITEBYTECODE='1'),
                        creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0)
                row['exit_code'] = completed.returncode
                supervision = prefix.with_name(prefix.name + '_supervision.json')
                if supervision.is_file():
                    outcome = json.loads(supervision.read_text(encoding='utf-8'))
                    row['error'] = outcome.get('stop_reason')
                    if row['error'] == 'USER_CANCELLED':
                        row.update(run_status='CANCELLED', status='CANCELLED')
                measured = prefix.with_name(prefix.name + '_measurement.json')
                params = prefix.with_name(prefix.name + '_parameters.json')
                if params.is_file():
                    value = json.loads(params.read_text(encoding='utf-8'))
                    metrics = value['metrics']
                    row.update({m['parameter']: m['value'] for m in metrics})
                    row.update(unit=value['unit'], vertical_axis=config['vertical_axis'], down_sign=config['down_sign'],
                               input_preparation=value.get('input_preparation'),
                               complete_length=value['complete_model_length_available'],
                               radius_valid_arc_fraction=value['coverage']['V_radius_out'],
                               connected_radius_valid_arc_fraction=value.get('connected_section_radius_model_units', {}).get('valid_arc_fraction'),
                               parameters=file_reference(params, experiment), measurement=file_reference(measured, experiment))
                    if row['run_status'] != 'CANCELLED':
                        row['status'] = 'COMPLETE_AXIS_WITH_EXPLICIT_RADIUS_COVERAGE' if row['complete_length'] else 'INCOMPLETE_AXIS'
                        row['run_status'] = 'COMPLETE' if completed.returncode == 0 else 'PARTIAL'
                    for metric in metrics:
                        parameter_rows.append(dict(index=index, sample_id=label, input=str(source),
                                                   run_status=row['run_status'], **metric))
                if completed.returncode and not row['error']:
                    row['error'] = 'MEASUREMENT_FAILED_OR_PARTIAL_SEE_LAUNCHER_LOG'
            except (ValueError, TypeError, KeyError, OSError) as exc:
                row['error'] = str(exc)
        if not any(r['index'] == index for r in parameter_rows):
            parameter_rows.append(dict(index=index, sample_id=row['sample_id'], input=row['input'],
                                       run_status=row['run_status'], parameter=None, value=None, status=row['error'] or row['status']))
        records.append(row)
        with (folder / 'batch_progress.jsonl').open('a', encoding='utf-8') as stream:
            stream.write(json.dumps(dict(index=index, total=len(inputs), phase='case_complete', record=row), ensure_ascii=False) + '\n')
        print(json.dumps(row, ensure_ascii=False), flush=True)
    elapsed = time.monotonic() - started
    summary = dict(schema='root-measurer-batch-v002', id=batch_id, total=len(inputs), measurements=records,
                   elapsed_seconds=elapsed, all_complete_lengths=all(r['complete_length'] for r in records),
                   complete_radius_claim=False, missing_intervals_retained=True,
                   failed_count=sum(r['run_status'] == 'FAILED' for r in records),
                   cancelled_count=sum(r['run_status'] == 'CANCELLED' for r in records),
                   partial_count=sum(r['run_status'] == 'PARTIAL' for r in records))
    summary['status'] = 'COMPLETE' if all(r['run_status'] == 'COMPLETE' for r in records) else 'COMPLETE_WITH_EXPLICIT_GAPS'
    for name, rows in [('batch_summary.csv', records), ('batch_parameters.csv', parameter_rows)]:
        fields = list(dict.fromkeys(k for row in rows for k in row))
        with (folder / name).open('x', encoding='utf-8-sig', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows({k: _cell(v) for k, v in row.items()} for row in rows)
    write_json(folder / 'batch_summary.json', summary)
    with (experiment / 'EXPERIMENT_LOG.md').open('a', encoding='utf-8') as stream:
        stream.write('\n' + local_date() + ' RESULT mixed batch ' + folder.relative_to(experiment).as_posix() +
                     ': ' + str(len(records)) + ' input rows retained; parameters and missing statuses exported.\n')
    return records


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    ap.add_argument('--defaults-json', default='{}')
    ap.add_argument('--manifest-base', type=Path)
    ap.add_argument('--scale-mm', type=float)
    ap.add_argument('--vertical-axis', choices=['x', 'y', 'z'], default='z')
    ap.add_argument('--down-sign', type=int, choices=[-1, 1], default=-1)
    ap.add_argument('--max-memory-mb', type=int, default=8192)
    ap.add_argument('--section-workers', type=int, default=1)
    ap.add_argument('--cancel-file', type=Path)
    a = ap.parse_args()
    value = json.loads(a.manifest.read_text(encoding='utf-8-sig'))
    if not isinstance(value, dict) or not isinstance(value.get('defaults', {}), dict):
        raise ValueError('BATCH_MANIFEST_OBJECT_REQUIRED')
    defaults = dict(json.loads(a.defaults_json), **value.get('defaults', {}))
    run(value['inputs'], a.output_dir, a.scale_mm, a.max_memory_mb, a.cancel_file,
        a.vertical_axis, a.down_sign, a.section_workers, defaults=defaults, manifest_base=a.manifest_base or a.manifest.resolve().parent)
