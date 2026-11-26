"""Mixed formats, per-case calibration, explicit failures and restored history."""
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from root_measurer import batch
from root_measurer.runtime import initialize_workspace
from root_model_delivery_service import batch_history
from measure_model_batch import run


def test_real_batch_continues_and_calibrates_each_volume(tmp_path):
    raw = np.zeros((5, 5, 21), np.uint8)
    raw[1:4, 1:4, :] = 1
    raw[0, 0, 0] = 2
    np.save(tmp_path / 'first.npy', raw)
    np.savez(tmp_path / 'second.npz', mask=raw.transpose(2, 1, 0),
             spacing=np.array([.02, .02, .02]), axis_order=np.array('zyx'), unit=np.array('cm'))
    manifest = tmp_path / 'inputs.json'
    manifest.write_text(json.dumps({'inputs': [
        {'id': 'missing', 'path': 'absent.npy', 'voxel_pitch': .1},
        {'id': 'first', 'path': 'first.npy', 'voxel_pitch': .1, 'voxel_unit': 'mm', 'voxel_labels': [1]},
        {'id': 'second', 'path': 'second.npz', 'voxel_labels': [1], 'vertical_axis': 'y', 'down_sign': 1},
    ]}))
    result = batch(manifest, workspace=tmp_path / 'ws')
    summary = result.summary
    assert summary['total'] == 3 and summary['failed_count'] == 1
    missing, first, second = summary['measurements']
    assert missing['run_status'] == 'FAILED' and missing['error']
    assert first['run_status'] == second['run_status'] == 'COMPLETE'
    assert first['unit'] == second['unit'] == 'mm'
    assert first['volume'] == pytest.approx(189 * .1**3)
    assert second['volume'] == pytest.approx(189 * .2**3)
    assert second['vertical_axis'] == 'y' and second['down_sign'] == 1
    assert first['input_preparation']['source']['labels'] == [1]
    assert (result.folder / 'case001_imported_voxel/entity_prediction.json').is_file()
    assert (result.folder / 'case002_imported_voxel/entity_prediction.json').is_file()
    rows = list(csv.DictReader((result.folder / 'batch_parameters.csv').open(encoding='utf-8-sig')))
    assert {r['sample_id'] for r in rows} == {'missing', 'first', 'second'}
    assert any(r['parameter'] == 'volume' and r['unit'] == 'mm^3' for r in rows)
    assert all('status' in r and 'covered_arc_fraction' in r for r in rows)
    history = batch_history(tmp_path / 'ws')
    assert summary['id'] in history and history[summary['id']]['history']


def test_mixed_commands_and_partial_metrics_keep_denominators(tmp_path, monkeypatch):
    import measure_model_batch as runner
    for name in ('root.ply', 'entity.json'):
        (tmp_path / name).write_text('ply\n' if name.endswith('ply') else '{}')
    np.save(tmp_path / 'root.npy', np.ones((3, 3, 9), bool))
    workspace = initialize_workspace(tmp_path / 'ws')
    commands = []

    def execute(command, **kwargs):
        commands.append(command)
        prefix = Path(command[command.index('--output-prefix') + 1])
        prefix.with_name(prefix.name + '_measurement.json').write_text('{}')
        prefix.with_name(prefix.name + '_parameters.json').write_text(json.dumps(dict(
            unit='mm', complete_model_length_available=False, coverage={'V_radius_out': .25},
            metrics=[dict(parameter='extracted_axis_length', value=4, unit='mm', status='PARTIAL', covered_arc_fraction=None),
                     dict(parameter='radius_mean', value=None, unit='mm', status='MISSING', covered_arc_fraction=.25)])))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(runner.subprocess, 'run', execute)
    rows = run([
        {'path': 'root.ply', 'mesh_target_pitch': .1, 'mesh_target_reason': 'source feature scale'},
        {'path': 'root.npy', 'voxel_pitch': .1, 'voxel_labels': [1]},
        {'path': 'root.npy', 'voxel_pitch': -1},
        {'path': 'entity.json', 'scale_mm': 10},
    ], workspace / 'results/model_delivery_local/run001', manifest_base=tmp_path)
    assert len(rows) == 4 and len(commands) == 3
    assert '--mesh-workflow' in commands[0] and '--volume' in commands[1] and '--entity' in commands[2]
    assert commands[2][commands[2].index('--scale-mm') + 1] == '10'
    assert rows[2]['run_status'] == 'FAILED'
    assert rows[0]['status'] == 'INCOMPLETE_AXIS' and rows[0]['extracted_axis_length'] == 4
    assert rows[0]['radius_mean'] is None and rows[0]['radius_valid_arc_fraction'] == .25
    with pytest.raises(ValueError, match='ALREADY_EXISTS'):
        run(['entity.json'], workspace / 'results/model_delivery_local/run001', manifest_base=tmp_path)


def test_cancelled_inputs_and_interrupted_history_are_preserved(tmp_path):
    workspace = initialize_workspace(tmp_path / 'ws')
    cancel = tmp_path / 'cancel.json'
    cancel.write_text('{}')
    folder = workspace / 'results/model_delivery_local/run001'
    rows = run(['a', 'b', 'c'], folder, cancel_file=cancel)
    assert len(rows) == 3 and all(r['run_status'] == 'CANCELLED' for r in rows)
    interrupted = workspace / 'results/model_delivery_local/run002'
    interrupted.mkdir()
    (interrupted / 'batch_request.json').write_text(json.dumps(dict(total=2)))
    records = batch_history(workspace)
    assert len(records) == 2
    assert all(record['process'] is None for record in records.values())


def test_http_restores_history_downloads_and_accepts_mixed_manifest(tmp_path, monkeypatch):
    import threading
    from http.server import ThreadingHTTPServer
    from unittest.mock import Mock
    from urllib.request import Request, urlopen
    from urllib.error import HTTPError
    import root_model_delivery_service as service
    workspace = initialize_workspace(tmp_path / 'ws')
    folder = workspace / 'results/model_delivery_local/run001'
    folder.mkdir()
    (folder / 'batch_request.json').write_text(json.dumps(dict(total=2)))
    (folder / 'batch_summary.json').write_text(json.dumps(dict(total=2, measurements=[], schema='root-measurer-batch-v002')))
    (folder / 'batch_parameters.csv').write_text('parameter,value,status\nvolume,1,AVAILABLE\n')
    servers = []

    def factory(address, handler):
        server = ThreadingHTTPServer(address, handler)
        servers.append(server)
        server.serve_forever = lambda: None
        return server

    monkeypatch.setattr(service, 'ThreadingHTTPServer', factory)
    service.serve_workspace(workspace)
    server = servers[0]
    thread = threading.Thread(target=lambda: ThreadingHTTPServer.serve_forever(server, poll_interval=.01), daemon=True)
    thread.start()
    url = f'http://127.0.0.1:{server.server_address[1]}'

    def request(path, data=None, token=None):
        headers = {'Content-Type': 'application/json'}
        if token:
            headers['X-Root-Local-Token'] = token
        with urlopen(Request(url + path, data=None if data is None else json.dumps(data).encode(), headers=headers), timeout=10) as response:
            raw = response.read()
            return json.loads(raw) if response.headers.get('Content-Type', '').startswith('application/json') else raw

    try:
        history = request('/api/state')['batches'][0]
        assert history['history'] and not history['running']
        for _ in range(5):
            with pytest.raises(HTTPError) as denied:
                request('/api/jobs', {'entity_id': 'unauthorized-control'})
            assert denied.value.code == 403
        assert b'volume,1,AVAILABLE' in request('/api/batches/' + history['id'] + '/parameters')
        manifest = tmp_path / 'batch.json'
        manifest.write_text(json.dumps(dict(defaults={'section_workers': 2}, inputs=[
            {'path': 'missing.npy', 'voxel_pitch': .1},
            {'path': 'root.ply', 'mesh_target_pitch': .2, 'mesh_target_reason': 'source scale'},
        ])))
        process = Mock()
        process.poll.return_value = None
        launcher = Mock(return_value=process)
        monkeypatch.setattr(service.subprocess, 'Popen', launcher)
        token = request('/api/interaction-session')['token']
        created = request('/api/batches', dict(manifest_path=str(manifest)), token)
        assert created['batch_id']
        command = launcher.call_args.args[0]
        assert command[command.index('--manifest-base') + 1] == str(tmp_path)
        recorded = json.loads(Path(command[command.index('--manifest') + 1]).read_text())
        assert recorded['inputs'][0]['path'] == 'missing.npy'
        assert recorded['defaults']['section_workers'] == 2
    finally:
        server.shutdown();server.server_close();thread.join(timeout=5)
