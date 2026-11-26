"""Tests for universal 3D binary/numeric voxel volume import (.npy/.npz/.tif)."""
import math
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock
import numpy as np
import pytest

from root_measurer import cli, import_volume, measure, morphology, query_radius
from root_measurer.runtime import activate


def _cylinder_volume(radius=0.1, length=1.0, pitch=0.02):
    xs = np.arange(-6, 7) * pitch
    ys = np.arange(-6, 7) * pitch
    zs = np.arange(0, int(round(length / pitch)) + 1) * pitch
    grid_x, grid_y, _ = np.meshgrid(xs, ys, zs, indexing="ij")
    return (grid_x * grid_x + grid_y * grid_y) <= (radius * radius)


def test_import_npy_and_end_to_end_measure_with_workers(tmp_path):
    vol = _cylinder_volume()
    npy_path = tmp_path / "ct_root.npy"
    np.save(npy_path, vol)

    result = measure(
        npy_path,
        workspace=tmp_path / "ws",
        voxel_pitch=0.02,
        voxel_axis_order="xyz",
        scale_mm=100,
        section_workers=2,
    )
    params = result.parameters
    assert params["complete_model_length_available"]
    metrics = {m["parameter"]: m["value"] for m in params["metrics"]}
    assert 95.0 <= metrics["model_length"] <= 106.0
    assert params["geometry"]["kind"] == "EXACT_NATIVE_CUBE_UNION"

    sec = query_radius(result.measurement_path, "path00000", 0.5)
    assert sec["radius"] is not None
    assert not sec["stored_samples_interpolated"]
    assert math.isclose(sec["radius"], 10.0, rel_tol=0.12)

    morph = morphology(result.measurement_path, vertical_axis="z", down_sign=-1)
    assert morph["settings"]["vertical_axis"] == "z"


def test_import_npz_and_zyx_axis_order_and_tiff(tmp_path):
    vol_xyz = np.zeros((5, 7, 9), dtype=np.uint8)
    vol_xyz[1:3, 2:4, 3:7] = 255
    vol_zyx = np.transpose(vol_xyz, (2, 1, 0))

    npz_path = tmp_path / "stack.npz"
    np.savez_compressed(npz_path, mask=vol_zyx)
    out_npz = tmp_path / "out_npz"
    info_npz = import_volume(npz_path, out_npz, pitch=0.05, axis_order="zyx", threshold=127)
    assert info_npz["occupied_voxels"] == int((vol_xyz > 0).sum())
    assert info_npz["axis_order"] == "zyx"

    from skimage.io import imsave
    tif_path = tmp_path / "stack.tif"
    imsave(str(tif_path), vol_zyx, check_contrast=False)
    out_tif = tmp_path / "out_tif"
    info_tif = import_volume(tif_path, out_tif, pitch=0.05, axis_order="zyx", threshold=127)

    activate()
    from root_native_entity_measurement import load_entity
    s_npz, _, _ = load_entity(info_npz["entity_path"])
    s_tif, _, _ = load_entity(info_tif["entity_path"])
    assert np.array_equal(s_npz.cells, s_tif.cells)
    assert s_npz.pitch == pytest.approx(0.05)


def test_import_validation_and_cli(tmp_path):
    bad_2d = tmp_path / "slice.npy"
    np.save(bad_2d, np.ones((4, 4), dtype=bool))
    with pytest.raises(ValueError, match="VOXEL_VOLUME_3D_REQUIRED"):
        import_volume(bad_2d, tmp_path / "o1", pitch=0.1)

    empty = tmp_path / "empty.npy"
    np.save(empty, np.zeros((3, 3, 3), dtype=bool))
    with pytest.raises(ValueError, match="EMPTY_VOXEL_VOLUME"):
        import_volume(empty, tmp_path / "o2", pitch=0.1)

    multi = tmp_path / "multi.npz"
    np.savez(multi, a=np.ones((2, 2, 2), dtype=bool), b=np.ones((2, 2, 2), dtype=bool))
    with pytest.raises(ValueError, match="NPZ_ARRAY_KEY_REQUIRED"):
        import_volume(multi, tmp_path / "o3", pitch=0.1)
    ok = import_volume(multi, tmp_path / "o3_ok", pitch=0.1, array_key="b")
    assert ok["occupied_voxels"] == 8

    good = tmp_path / "cube.npy"
    np.save(good, np.ones((3, 3, 3), dtype=bool))
    with pytest.raises(ValueError, match="POSITIVE_FINITE_VOXEL_PITCH_REQUIRED"):
        measure(good, workspace=tmp_path / "ws_bad")
    assert not (tmp_path / "ws_bad").exists()

    rc = cli.main(["import-volume", str(good), "--output", str(tmp_path / "cli_out"), "--pitch", "0.05", "--axis-order", "zyx"])
    assert rc == 0
    assert (tmp_path / "cli_out" / "entity_prediction.json").is_file()


@pytest.mark.parametrize("shape", [(3, 7, 9), (4, 7, 9), (9, 7, 3), (9, 7, 4)])
@pytest.mark.parametrize("axis_order", ["xyz", "zyx"])
def test_tiff_preserves_grayscale_volume_axes(tmp_path, shape, axis_order):
    from tifffile import imwrite
    from root_native_entity_measurement import load_entity
    from root_voxel_import import load_3d_voxel_array

    volume = np.zeros(shape, np.uint8)
    volume[0, 1, 1] = 255
    volume[-1, -2, -2] = 255
    path = tmp_path / "grayscale_stack.tif"
    imwrite(path, volume, photometric="minisblack")
    loaded, _ = load_3d_voxel_array(path)
    assert loaded.shape == shape
    assert np.array_equal(loaded, volume)
    imported = import_volume(path, tmp_path / "solid", pitch=0.05, axis_order=axis_order, threshold=127)
    solid, _, metadata = load_entity(imported["entity_path"])
    expected = np.argwhere(volume)
    if axis_order == "zyx":
        expected = expected[:, [2, 1, 0]]
    assert sorted(map(tuple, solid.cells)) == sorted(map(tuple, expected))
    assert metadata["source_volume"]["raw_shape"] == list(shape)


def test_api_dispatches_volume_without_importing_in_caller(tmp_path, monkeypatch):
    from root_measurer import api
    import root_voxel_import
    source = tmp_path / "volume.npy"
    np.save(source, np.ones((3, 3, 9), dtype=bool))
    importer = Mock(side_effect=AssertionError("import must run in supervised worker"))
    monkeypatch.setattr(root_voxel_import, "import_voxel_volume", importer)
    runner = Mock(return_value=SimpleNamespace(returncode=0))
    monkeypatch.setattr(api.subprocess, "run", runner)
    result = measure(source, workspace=tmp_path / "workspace", voxel_pitch=.05,
                     voxel_axis_order="zyx", voxel_threshold=127, voxel_array_key="mask", max_memory_mb=256)
    command = runner.call_args.args[0]
    assert command[command.index("--volume") + 1] == str(source.resolve())
    assert command[command.index("--voxel-array-key") + 1] == "mask"
    assert command[command.index("--max-memory-mb") + 1] == "256"
    assert not (result.prefix.parent / "imported_voxel").exists()
    importer.assert_not_called()


def test_volume_measurement_low_budget_retains_supervision(tmp_path):
    from root_measurer.api import MeasurementError
    source = tmp_path / "volume.npy"
    np.save(source, _cylinder_volume())
    with pytest.raises(MeasurementError) as failed:
        measure(source, workspace=tmp_path / "workspace", voxel_pitch=.02, max_memory_mb=1)
    prefix = failed.value.prefix
    report = json.loads(prefix.with_name(prefix.name + "_supervision.json").read_text())
    assert report["stop_reason"] == "WORKER_MEMORY_BUDGET_EXCEEDED"
    assert report["memory_scope"] == "worker_process_tree"


def test_existing_supervisor_counts_children_and_stops_them(monkeypatch):
    import psutil
    from measure_model_delivery import _worker_memory, _stop_worker
    parent = Mock()
    child = Mock()
    vanished = Mock()
    parent.memory_info.return_value = SimpleNamespace(rss=32, private=24, peak_wset=40)
    child.memory_info.return_value = SimpleNamespace(rss=2048, private=1800, peak_wset=2100)
    vanished.memory_info.side_effect = psutil.NoSuchProcess(123)
    parent.children.return_value = [child, vanished]
    monkeypatch.setattr(psutil, "Process", Mock(return_value=parent))
    process = Mock(pid=123)
    process.poll.return_value = None
    memory = _worker_memory(process)
    assert memory == dict(working_set=2080, private_bytes=1824, peak_working_set=2140)
    waiter = Mock(return_value=([vanished], [child]))
    monkeypatch.setattr(psutil, "wait_procs", waiter)
    _stop_worker(process)
    process.terminate.assert_called_once()
    child.terminate.assert_called_once()
    child.kill.assert_called_once()
    process.wait.assert_called_once()


def test_web_dispatches_volume_without_blocking_server_for_import(tmp_path, monkeypatch):
    import threading
    from http.server import ThreadingHTTPServer
    from urllib.request import Request, urlopen
    from root_measurer.runtime import initialize_workspace
    import root_model_delivery_service as service
    import root_voxel_import

    source = tmp_path / "volume.npy"
    np.save(source, np.ones((3, 3, 9), dtype=bool))
    workspace = initialize_workspace(tmp_path / "workspace")
    servers = []

    def factory(address, handler):
        server = ThreadingHTTPServer(address, handler)
        servers.append(server)
        server.serve_forever = lambda: None
        return server

    monkeypatch.setattr(service, "ThreadingHTTPServer", factory)
    importer = Mock(side_effect=AssertionError("HTTP handler must not import volumes"))
    monkeypatch.setattr(root_voxel_import, "import_voxel_volume", importer)
    service.serve_workspace(workspace, input_paths=[source])
    server = servers[0]
    thread = threading.Thread(target=lambda: ThreadingHTTPServer.serve_forever(server, poll_interval=.01), daemon=True)
    thread.start()
    url = f"http://127.0.0.1:{server.server_address[1]}"

    def request(path, data=None, token=None):
        headers = {"Content-Type": "application/json"}
        if token:
            headers["X-Root-Local-Token"] = token
        req = Request(url + path, data=None if data is None else json.dumps(data).encode(), headers=headers)
        with urlopen(req, timeout=10) as response:
            return json.loads(response.read())

    process = Mock()
    process.poll.return_value = None
    launcher = Mock(return_value=process)
    try:
        entity = request("/api/state")["entities"][0]
        token = request("/api/interaction-session")["token"]
        monkeypatch.setattr(service.subprocess, "Popen", launcher)
        result = request("/api/jobs", dict(entity_id=entity["id"], voxel_pitch=.05,
                         voxel_axis_order="zyx", voxel_threshold=127, voxel_array_key="mask", section_workers=2), token)
        assert result["job_id"]
        command = launcher.call_args.args[0]
        assert command[command.index("--volume") + 1] == str(source.resolve())
        assert command[command.index("--voxel-array-key") + 1] == "mask"
        assert command[command.index("--section-workers") + 1] == "2"
        assert not list(workspace.glob("results/model_delivery_local/run*/imported_voxel"))
        importer.assert_not_called()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def test_label_selection_and_explicit_resampling_preserve_source(tmp_path):
    import json
    from root_native_entity_measurement import load_entity
    raw = np.zeros((3, 4, 5), np.uint16)
    raw[1, 1:3, 1:4] = 7
    raw[0, 0, 0] = 2
    source = tmp_path / 'labels.npy'
    np.save(source, raw)
    original = source.read_bytes()
    metadata = tmp_path / 'acquisition.json'
    metadata.write_text(json.dumps(dict(spacing=[1, 2, 3], axis_order='xyz', unit='mm', origin=[10, 20, 30])))
    with pytest.raises(ValueError, match='ANISOTROPIC_VOXELS_REQUIRE_EXPLICIT'):
        import_volume(source, tmp_path / 'rejected', labels=[7], metadata=metadata)
    assert not (tmp_path / 'rejected').exists()
    info = import_volume(source, tmp_path / 'resampled', labels=[7], metadata=metadata, resample_pitch=1)
    solid, _, meta = load_entity(info['entity_path'])
    assert len(solid.cells) == 6 * 2 * 3
    assert meta['source_volume']['labels'] == [7]
    assert meta['source_volume']['source_origin_xyz'] == [10, 20, 30]
    assert meta['resampling']['applied']
    assert meta['resampling']['original_volume'] == meta['resampling']['output_volume'] == 36
    assert source.read_bytes() == original


def test_npz_metadata_sets_axis_order_spacing_and_unit(tmp_path):
    from root_native_entity_measurement import load_entity
    raw = np.zeros((9, 7, 5), np.uint8)
    raw[1:7, 2:5, 1:4] = 3
    source = tmp_path / 'metadata.npz'
    np.savez(source, segmentation=raw, spacing=np.array([.1, .1, .1]), axis_order=np.array('zyx'), unit=np.array('cm'))
    info = import_volume(source, tmp_path / 'solid', labels=[3])
    solid, _, _ = load_entity(info['entity_path'])
    assert sorted(map(tuple, solid.cells)) == sorted(map(tuple, np.argwhere(raw)[:, [2, 1, 0]]))
    assert info['settings']['scale_mm_per_model_unit'] == 10


def test_ome_tiff_reads_physical_spacing_without_reordering_array(tmp_path):
    from tifffile import imwrite
    from root_voxel_import import volume_metadata
    raw = np.ones((4, 7, 9), np.uint8)
    source = tmp_path / 'root.ome.tif'
    imwrite(source, raw, ome=True, photometric='minisblack', metadata=dict(axes='ZYX',
            PhysicalSizeX=50, PhysicalSizeY=50, PhysicalSizeZ=200,
            PhysicalSizeXUnit='µm', PhysicalSizeYUnit='µm', PhysicalSizeZUnit='µm'))
    meta = volume_metadata(source)
    assert meta['spacing'] == [50, 50, 200]
    assert meta['axis_order'] == 'zyx'
    with pytest.raises(ValueError, match='ANISOTROPIC'):
        import_volume(source, tmp_path / 'rejected')
    info = import_volume(source, tmp_path / 'solid', resample_pitch=50)
    assert info['occupied_voxels'] == raw.size * 4
    assert info['settings']['scale_mm_per_model_unit'] == .001


def test_imported_native_package_keeps_acquisition_calibration(tmp_path):
    raw = np.ones((3, 3, 21), bool)
    source = tmp_path / 'root.npy'
    np.save(source, raw)
    imported = import_volume(source, tmp_path / 'solid', pitch=.02, unit='um')
    result = measure(imported['entity_path'], workspace=tmp_path / 'ws')
    params = result.parameters
    assert params['unit'] == 'mm' and params['scale_mm_per_model_unit'] == .001
    volume = next(m for m in params['metrics'] if m['parameter'] == 'volume')
    assert volume['value'] == pytest.approx(raw.size * (.02 * .001)**3)
