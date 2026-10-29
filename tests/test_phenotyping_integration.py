"""Public API, immutable provenance, CLI snapshots and the real loopback routes."""
import copy
import csv
import io
import json
import threading
from http.server import ThreadingHTTPServer
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pytest

from root_measurer import measure, phenotype
from root_measurer.cli import main
from root_model_delivery_exports import file_reference
from root_model_phenotyping import (ANNOTATION_SCHEMA, export_phenotyping,
                                    latest_annotation_profile, phenotyping_csv)


@pytest.fixture
def measurement(tmp_path):
    points = [[0, 0, 0], [0, 0, 2], [0, 0, 4], [0, 0, 6],
              [2, 0, 2], [4, 0, 2], [2, 0, 4]]
    edges = [(0, 1), (1, 2), (2, 3), (1, 4), (4, 5), (2, 6)]
    adjacency = [[] for _ in points]
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    report = dict(status='NATIVE_OUTPUT_WITH_EXPLICIT_GAPS', coordinate_units='model_units',
        scale_mm_per_model_unit=10,
        candidate_centerline=dict(points=points, adjacency=adjacency),
        model_bounds=dict(minimum=[0, 0, 0], maximum=[4, 1, 6]),
        geometric_endpoints=[dict(node_id=i, status='geometric_surface_estimate')
                             for i, row in enumerate(adjacency) if len(row) == 1],
        radius_curve=dict(sampling='edge_midpoint', rows=[
            dict(node_ids=list(edge), radius=.5 if i < 3 else .25,
                 available=True, status='CLOSED_DISK_CONTROL') for i, edge in enumerate(edges)]))
    path = tmp_path / 'tree_measurement.json'
    path.write_text(json.dumps(report), encoding='utf-8')
    return path


def rooted_profile(value, *, confirmed=False):
    profile = copy.deepcopy(value['annotation_profile'])
    profile['annotations'] = dict(schema_version='root-axis-annotations-v001',
        root_bases=[dict(node_id=0, root_type='primary')],
        endpoint_types=[dict(node_id=n, type='root_tip') for n in (3, 5, 6)],
        confirmed_components=[0] if confirmed else [])
    return profile


def test_public_geometry_candidate_and_confirmed_outputs_preserve_input(measurement):
    original = measurement.read_bytes()
    geometry = phenotype(measurement)
    assert geometry['architecture']['status'] == 'GEOMETRY_ONLY'
    assert geometry['distributions']['total_axis_length'] == 120
    assert geometry['unit'] == geometry['axis_traits']['unit'] == 'mm'
    assert geometry['measurement_coverage']['length_unit'] == 'model_units'
    assert geometry['measurement_coverage']['output_total_length'] == 12
    assert geometry['distributions']['distributions']['radial']['status'].startswith('UNAVAILABLE')
    candidate = phenotype(measurement, annotation_profile=rooted_profile(geometry))
    assert all(a['confirmation_status'] in ('CANDIDATE', 'AUTO_INFERRED') for a in candidate['architecture']['axes'])
    result = phenotype(measurement, annotation_profile=rooted_profile(candidate, confirmed=True))
    axes = result['axis_traits']['axes']
    mother = next(a for a in axes if a['parent_axis_id'] is None)
    assert len(axes) == 3 and mother['length'] == 60
    assert mother['branching']['lateral_count']['value'] == 2
    assert mother['branching']['lrbd_whole_axis']['value'] == pytest.approx(2 / 60)
    assert mother['developmental_order'] == 0
    assert all(a['confirmation_status'] == 'CONFIRMED' for a in axes)
    assert not result['confirmation_is_accuracy_validation']
    assert measurement.read_bytes() == original
    assert phenotype(measurement, annotation_profile=result['annotation_profile']) == result


def test_unscaled_declared_units_survive_all_public_layers(measurement):
    report = json.loads(measurement.read_text())
    report['coordinate_units'] = 'cm'
    report['scale_mm_per_model_unit'] = None
    measurement.write_text(json.dumps(report), encoding='utf-8')
    value = phenotype(measurement)
    assert value['unit'] == value['axis_traits']['unit'] == value['distributions']['unit'] == 'cm'
    assert value['distributions']['total_axis_length'] == 12


def test_profiles_are_bound_to_exact_measurement_bytes(measurement):
    value = phenotype(measurement)
    profile = rooted_profile(value, confirmed=True)
    measurement.write_bytes(measurement.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='ANNOTATION_MEASUREMENT_SHA256_MISMATCH'):
        phenotype(measurement, annotation_profile=profile)


def test_original_source_change_rejected_without_rewriting_measurement(measurement):
    source = measurement.parent / 'source.ply'
    source.write_bytes(b'original material control')
    report = json.loads(measurement.read_text())
    report['input_file_binding'] = file_reference(source)
    measurement.write_text(json.dumps(report), encoding='utf-8')
    assert phenotype(measurement)['source_status'] == 'ORIGINAL_MESH_HASH_VERIFIED'
    before = measurement.read_bytes()
    source.write_bytes(b'changed material control')
    with pytest.raises(ValueError, match='ORIGINAL_MESH_CHANGED'):
        phenotype(measurement)
    assert measurement.read_bytes() == before


@pytest.mark.parametrize('bad', [[], False, {'unknown': 1}, {'diameter_edges': [1, .5]}])
def test_public_invalid_settings_rejected(measurement, bad):
    with pytest.raises(ValueError):
        phenotype(measurement, options=bad)


def test_cli_snapshot_provenance_roundtrip_and_no_overwrite(measurement, capsys):
    prefix = measurement.parent / 'exports' / 'first'
    argv = ['phenotype', str(measurement), '--output', str(prefix), '--root-base', '0',
            '--confirm-component', '0', '--origin', '0', '0', '0', '--depth-bins', '3']
    before = file_reference(measurement)
    assert main(argv) == 0
    response = json.loads(capsys.readouterr().out)
    assert set(response['files']) == {'phenotyping', 'annotations', 'distributions',
                                     'axes', 'branches', 'profiles', 'summaries'}
    paths = {k: Path(p) for k, p in response['files'].items()}
    original_exports = {k: p.read_bytes() for k, p in paths.items()}
    result = phenotype(measurement, annotation_profile=paths['annotations'])
    assert result == json.loads(paths['phenotyping'].read_text(encoding='utf-8'))
    for kind in ('distributions', 'axes', 'branches', 'profiles', 'summaries'):
        rows = list(csv.DictReader(io.StringIO(paths[kind].read_text(encoding='utf-8-sig'))))
        assert rows and all(r['measurement_sha256'] == before['sha256'] for r in rows)
        assert all(r['output_unit'] == 'mm' and r['settings_sha256'] == result['settings_sha256']
                   for r in rows)
    assert main(argv) == 1
    assert 'PHENOTYPE_OUTPUT_EXISTS' in capsys.readouterr().err
    assert {k: p.read_bytes() for k, p in paths.items()} == original_exports
    assert file_reference(measurement) == before


def test_missing_values_remain_empty_and_root_labels_are_not_csv_formulas(measurement):
    report = json.loads(measurement.read_text())
    report['radius_curve']['rows'] = []
    measurement.write_text(json.dumps(report), encoding='utf-8')
    value = phenotype(measurement)
    profile = rooted_profile(value, confirmed=True)
    profile['annotations']['root_bases'][0]['root_type'] = '=custom_label'
    result = phenotype(measurement, annotation_profile=profile)
    rows = list(csv.DictReader(io.StringIO(phenotyping_csv(result, 'axes').decode('utf-8-sig'))))
    mother = next(r for r in rows if r['root_type'].endswith('custom_label'))
    assert mother['root_type'] == "'=custom_label"
    assert mother['diameter_distribution_mean'] == ''
    assert mother['branching_lateral_count_value'] == '2'
    assert result['distributions']['sample_coverage']['qualified_radius_arc_fraction'] == 0


def test_cli_can_change_the_base_of_an_exported_profile(measurement, capsys):
    first = phenotype(measurement)
    rooted = phenotype(measurement, annotation_profile=rooted_profile(first, confirmed=True))
    profile_path = Path(export_phenotyping(rooted, measurement.parent / 'old')['annotations'])
    assert main(['phenotype', str(measurement), '--annotations', str(profile_path),
                 '--root-base', '3', '--output', str(measurement.parent / 'reversed')]) == 0
    response = json.loads(capsys.readouterr().out)
    profile = json.loads(Path(response['files']['annotations']).read_text(encoding='utf-8'))
    assert profile['annotations']['root_bases'] == [dict(node_id=3, root_type='unknown')]
    assert not profile['annotations']['confirmed_axes']
    assert not profile['annotations']['confirmed_components']
    assert not any(e['node_id'] == 0 and e['type'] == 'root_base'
                   for e in profile['annotations']['endpoint_types'])


def test_latest_saved_profile_skips_corrupt_and_foreign_snapshots(measurement):
    value = phenotype(measurement)
    profile = rooted_profile(value)
    result = phenotype(measurement, annotation_profile=profile)
    folder = measurement.parent / 'phenotyping'
    export_phenotyping(result, folder / 'run00001' / 'analysis')
    for name, contents in [('run00002', '[]'), ('run00003', '{'),
                           ('run00004', json.dumps(dict(profile, measurement_sha256='0' * 64)))]:
        target = folder / name
        target.mkdir()
        (target / 'analysis_annotations.json').write_text(contents, encoding='utf-8')
    assert latest_annotation_profile(measurement) == result['annotation_profile']


def test_empty_axis_still_has_a_serializable_report(tmp_path):
    path = tmp_path / 'empty_measurement.json'
    path.write_text(json.dumps(dict(coordinate_units='model_units',
        candidate_centerline=dict(points=[], adjacency=[]))), encoding='utf-8')
    result = phenotype(path)
    assert result['axis_traits']['axes'] == []
    assert result['distributions']['total_axis_length'] is None
    json.dumps(result, allow_nan=False)


@pytest.fixture
def local_app(tmp_path, monkeypatch):
    """Exercise the production Handler through a real TCP server, without spawning a browser."""
    from build_native_measurement_acceptance import build
    import root_model_delivery_service as service
    inputs = tmp_path / 'inputs'
    build('straight', inputs)
    workspace = tmp_path / 'workspace'
    result = measure(inputs / 'straight', workspace=workspace, scale_mm=100,
                     vertical_axis='x', down_sign=1)
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
    url = 'http://127.0.0.1:%d' % server.server_address[1]

    def request(path, data=None, token=None):
        headers = {'Content-Type': 'application/json'}
        if token:
            headers['X-Root-Local-Token'] = token
        req = Request(url + path, data=None if data is None else json.dumps(data).encode(), headers=headers)
        with urlopen(req, timeout=10) as response:
            raw = response.read()
            return json.loads(raw) if 'application/json' in response.headers.get('Content-Type', '') else raw

    try:
        yield workspace, result, request
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_browser_routes_token_save_restore_export_and_source_binding(local_app):
    workspace, measured, request = local_app
    original = measured.measurement_path.read_bytes()
    result_id = request('/api/state')['deliveries'][0]['id']
    route = '/result/' + result_id + '/phenotyping'
    value = request(route + '.json')
    assert value['settings']['vertical_axis'] == 'x'
    assert value['settings']['down_sign'] == 1
    assert b'installPhenotyping' in request('/model_delivery_phenotyping.js')
    with pytest.raises(HTTPError) as denied:
        request(route, {'action': 'save'})
    assert denied.value.code == 403
    token = request('/api/interaction-session')['token']
    profile = copy.deepcopy(value['annotation_profile'])
    base = value['architecture']['components'][0]['candidate_root_base_node_ids'][0]
    profile['annotations'] = dict(schema_version='root-axis-annotations-v001',
        root_bases=[dict(node_id=base)], confirmed_components=[0])
    saved = request(route, dict(action='save', annotation_profile=profile,
                               options={'diameter_edges': [1, 3, 10]}), token)
    assert len(saved['saved_files']) == 7
    assert all((workspace / ref['path']).is_file() for ref in saved['saved_files'].values())
    restored = request(route + '.json')
    assert restored['annotation_profile'] == saved['annotation_profile']
    raw = request(route, dict(action='export', kind='distributions',
                             annotation_profile=restored['annotation_profile']), token)
    assert raw == phenotyping_csv(restored, 'distributions')
    bad = copy.deepcopy(profile)
    bad['measurement_sha256'] = '0' * 64
    with pytest.raises(HTTPError) as mismatch:
        request(route, dict(action='save', annotation_profile=bad), token)
    assert mismatch.value.code == 400
    assert len(list(measured.measurement_path.parent.glob('phenotyping/run*'))) == 1
    assert measured.measurement_path.read_bytes() == original
    measured.measurement_path.write_bytes(original + b'\n')
    try:
        with pytest.raises(HTTPError) as stale:
            request(route + '.json')
        assert stale.value.code == 409
    finally:
        measured.measurement_path.write_bytes(original)
