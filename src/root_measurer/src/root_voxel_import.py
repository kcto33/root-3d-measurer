"""Import standard 3D binary/numeric voxel volumes (.npy, .npz, 3D .tif/.tiff)
into the canonical SHA-bound native solid package consumed by Root 3D Measurer.
"""
from __future__ import annotations

import math
from pathlib import Path
import numpy as np

from root_model_delivery_exports import file_reference, write_json


VOXEL_VOLUME_SUFFIXES = (".npy", ".npz", ".tif", ".tiff")


def is_voxel_volume_path(path: Path | str) -> bool:
    return Path(path).suffix.lower() in VOXEL_VOLUME_SUFFIXES


def _validate_pitch(pitch: float) -> float:
    if isinstance(pitch, bool) or not isinstance(pitch, (int, float)) or not math.isfinite(float(pitch)) or float(pitch) <= 0.0:
        raise ValueError("POSITIVE_FINITE_VOXEL_PITCH_REQUIRED")
    return float(pitch)


def _validate_axis_order(axis_order: str) -> str:
    value = str(axis_order).strip().lower()
    if value not in ("xyz", "zyx"):
        raise ValueError("AXIS_ORDER_XYZ_OR_ZYX_REQUIRED")
    return value


def load_3d_voxel_array(source: Path, *, array_key: str | None = None) -> tuple[np.ndarray, str | None]:
    source = Path(source).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    suffix = source.suffix.lower()
    used_key: str | None = None
    if suffix == ".npy":
        arr = np.load(source, allow_pickle=False)
    elif suffix == ".npz":
        with np.load(source, allow_pickle=False) as archive:
            if array_key is not None:
                if array_key not in archive.files:
                    raise ValueError(f"NPZ_ARRAY_KEY_NOT_FOUND: {array_key}")
                used_key = array_key
                arr = archive[array_key]
            else:
                candidates = [k for k in archive.files if archive[k].ndim == 3]
                if len(candidates) == 1:
                    used_key = candidates[0]
                    arr = archive[used_key]
                elif not candidates:
                    raise ValueError("VOXEL_VOLUME_3D_REQUIRED")
                else:
                    raise ValueError(f"NPZ_ARRAY_KEY_REQUIRED: {', '.join(candidates)}")
    elif suffix in (".tif", ".tiff"):
        from skimage.io import imread

        arr = np.asarray(imread(str(source)))
    else:
        raise ValueError(f"UNSUPPORTED_VOXEL_VOLUME_FORMAT: {suffix}")

    arr = np.asarray(arr)
    if arr.ndim != 3:
        raise ValueError("VOXEL_VOLUME_3D_REQUIRED")
    return arr, used_key


def import_voxel_volume(
    source: Path | str,
    output_dir: Path | str,
    *,
    pitch: float,
    axis_order: str = "xyz",
    threshold: float = 0.0,
    array_key: str | None = None,
    label: str | None = None,
) -> dict:
    """Convert a 3D binary/numeric volume into a SHA-bound native entity package."""
    pitch_val = _validate_pitch(pitch)
    order = _validate_axis_order(axis_order)
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(float(threshold)):
        raise ValueError("FINITE_VOXEL_THRESHOLD_REQUIRED")
    threshold_val = float(threshold)

    src = Path(source).resolve()
    folder = Path(output_dir).resolve()
    if folder.exists() and any(folder.iterdir()):
        raise ValueError("IMPORT_OUTPUT_EXISTS")

    arr, used_key = load_3d_voxel_array(src, array_key=array_key)
    mask = arr.astype(bool) if arr.dtype == np.bool_ else (arr > threshold_val)
    coords = np.argwhere(mask)
    if not len(coords):
        raise ValueError("EMPTY_VOXEL_VOLUME")
    if order == "zyx":
        coords = coords[:, [2, 1, 0]]
    cells = np.ascontiguousarray(coords, dtype=np.int64)

    lo = cells.min(axis=0) - 1
    shape = (cells.max(axis=0) - lo + 2).astype(np.int64)
    local = np.ascontiguousarray(cells - lo, dtype=np.int64)
    origin = (lo.astype(np.float64) + 0.5) * pitch_val

    folder.mkdir(parents=True, exist_ok=True)
    crop = folder / "entity_crop0000.npz"
    np.savez_compressed(
        crop,
        source_negative_indices=local,
        regularized_negative_indices=local,
        pitch=np.array(pitch_val, dtype=np.float64),
        origin=origin,
        shape=shape,
    )
    crop_ref = file_reference(crop, folder)
    source_ref = file_reference(src)
    entity_label = (label or src.stem).strip() or "voxel_volume"
    occupied_shape = (cells.max(axis=0) - cells.min(axis=0) + 1).tolist()

    metadata = dict(
        schema_version="imported-voxel-entity-v001",
        root=entity_label,
        status="COMPLETE_IMPORTED_VOXEL_ENTITY",
        fields=[dict(cache=crop_ref)],
        duplicate_occupied_cells=False,
        all_source_negative_cells_accounted=int(len(cells)),
        storage=dict(
            pitch=pitch_val,
            occupied_voxels=int(len(cells)),
            global_occupied_shape=occupied_shape,
        ),
        source_volume=dict(
            file=source_ref,
            format=src.suffix.lower().lstrip("."),
            array_key=used_key,
            axis_order=order,
            threshold=threshold_val,
            raw_shape=[int(x) for x in arr.shape],
        ),
        geometry="imported 3D binary voxel cube union; coordinates=(integer_index+0.5)*pitch",
        no_closing_or_postprocessing=True,
    )
    entity_json = folder / "entity_prediction.json"
    write_json(entity_json, metadata)
    return dict(
        entity_path=str(entity_json),
        crop_path=str(crop),
        root=entity_label,
        pitch=pitch_val,
        occupied_voxels=int(len(cells)),
        axis_order=order,
        source_volume=source_ref,
    )
