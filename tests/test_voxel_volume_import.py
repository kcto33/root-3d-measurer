"""Tests for universal 3D binary/numeric voxel volume import (.npy/.npz/.tif)."""
import math
from pathlib import Path
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
