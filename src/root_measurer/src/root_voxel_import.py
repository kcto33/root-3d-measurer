"""Import standard 3D binary/numeric voxel volumes (.npy, .npz, 3D .tif/.tiff)
into the canonical SHA-bound native solid package consumed by Root 3D Measurer.
"""
from __future__ import annotations

import math
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

from root_model_delivery_exports import file_reference, write_json


VOXEL_VOLUME_SUFFIXES = (".npy", ".npz", ".tif", ".tiff")
UNIT_TO_MM = {"mm": 1.0, "cm": 10.0, "m": 1000.0, "um": .001, "µm": .001, "μm": .001, "nm": .000001}


def is_voxel_volume_path(path: Path | str) -> bool:
    return Path(path).suffix.lower() in VOXEL_VOLUME_SUFFIXES


def _validate_pitch(pitch: float) -> float:
    if isinstance(pitch, bool) or not isinstance(pitch, (int, float)) or not math.isfinite(float(pitch)) or float(pitch) <= 0.0:
        raise ValueError("POSITIVE_FINITE_VOXEL_PITCH_REQUIRED")
    return float(pitch)


def _validate_axis_order(axis_order: str) -> str:
    value = str(axis_order).strip().lower()
    if value not in ("auto", "xyz", "zyx"):
        raise ValueError("AXIS_ORDER_AUTO_XYZ_OR_ZYX_REQUIRED")
    return value


def _vector(value, name, *, positive=False):
    if not isinstance(value, (list, tuple, np.ndarray)) or len(value) != 3:
        raise ValueError(name + "_XYZ_TRIPLE_REQUIRED")
    if any(isinstance(v, (bool, np.bool_)) or not math.isfinite(float(v)) or (positive and float(v) <= 0) for v in value):
        raise ValueError(name + "_FINITE_VALUES_REQUIRED")
    return [float(v) for v in value]


def volume_metadata(source, metadata=None):
    """Read acquisition metadata only; spacing and origin are always XYZ."""
    source = Path(source)
    value = {}
    if source.suffix.lower() == ".npz":
        with np.load(source, allow_pickle=False) as archive:
            for name in ("spacing", "origin", "axis_order", "unit"):
                if name in archive.files:
                    value[name] = archive[name].tolist()
    elif source.suffix.lower() in (".tif", ".tiff"):
        from tifffile import TiffFile
        with TiffFile(source) as tif:
            axes = tif.series[0].axes.lower()
            if axes in ("xyz", "zyx"):
                value["axis_order"] = axes
            elif axes in ("qyx", "iyx"):
                value["axis_order"] = "zyx"
            if tif.ome_metadata:
                root = ET.fromstring(tif.ome_metadata)
                pixels = next((n for n in root.iter() if n.tag.rsplit("}", 1)[-1] == "Pixels"), None)
                if pixels is not None and all("PhysicalSize" + a in pixels.attrib for a in "XYZ"):
                    spacing = [float(pixels.attrib["PhysicalSize" + a]) for a in "XYZ"]
                    units = [pixels.attrib.get("PhysicalSize" + a + "Unit") for a in "XYZ"]
                    if len(set(units)) == 1:
                        value.update(spacing=spacing, unit=units[0])
                    elif all(u in UNIT_TO_MM for u in units):
                        value.update(spacing=[v * UNIT_TO_MM[u] for v, u in zip(spacing, units)], unit="mm")
                    else:
                        raise ValueError("INCONSISTENT_OME_SPACING_UNITS")
    if metadata is not None:
        extra = json.loads(Path(metadata).read_text(encoding="utf-8-sig")) if isinstance(metadata, (str, Path)) else metadata
        if not isinstance(extra, dict):
            raise ValueError("VOLUME_METADATA_OBJECT_REQUIRED")
        value.update({k: extra[k] for k in ("spacing", "origin", "axis_order", "unit") if k in extra})
    return value


def resolve_volume_settings(source, *, pitch=None, spacing=None, axis_order="auto", unit=None, metadata=None, resample_pitch=None):
    meta = volume_metadata(source, metadata)
    if pitch is not None and spacing is not None:
        raise ValueError("CHOOSE_VOXEL_PITCH_OR_SPACING")
    raw_spacing = spacing if spacing is not None else ([_validate_pitch(pitch)] * 3 if pitch is not None else meta.get("spacing"))
    if raw_spacing is None:
        raise ValueError("POSITIVE_FINITE_VOXEL_PITCH_REQUIRED: supply pitch, XYZ spacing or acquisition metadata")
    original = _vector(raw_spacing, "SPACING", positive=True)
    order = _validate_axis_order(axis_order)
    if order == "auto":
        order = _validate_axis_order(meta.get("axis_order", "xyz"))
        if order == "auto":
            order = "xyz"
    anisotropic = not np.allclose(original, original[0], rtol=1e-12, atol=0)
    if anisotropic and resample_pitch is None:
        raise ValueError("ANISOTROPIC_VOXELS_REQUIRE_EXPLICIT_RESAMPLE_PITCH")
    target = _validate_pitch(resample_pitch) if resample_pitch is not None else original[0]
    chosen_unit = unit if unit is not None else meta.get("unit")
    if chosen_unit is not None:
        chosen_unit = str(chosen_unit).strip().lower()
        if chosen_unit not in UNIT_TO_MM:
            raise ValueError("VOLUME_UNIT_MM_CM_M_UM_NM_REQUIRED")
    return dict(spacing_xyz=original, axis_order=order, pitch=target, unit=chosen_unit,
                scale_mm_per_model_unit=UNIT_TO_MM.get(chosen_unit),
                source_origin_xyz=_vector(meta.get("origin", [0, 0, 0]), "ORIGIN"),
                resampled=not np.allclose(original, target, rtol=1e-12, atol=0))


def load_3d_voxel_array(source: Path, *, array_key: str | None = None) -> tuple[np.ndarray, str | None]:
    source = Path(source).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    suffix = source.suffix.lower()
    used_key: str | None = None
    if suffix == ".npy":
        arr = np.load(source, allow_pickle=False, mmap_mode="r")
    elif suffix == ".npz":
        with np.load(source, allow_pickle=False) as archive:
            if array_key is not None:
                if array_key not in archive.files:
                    raise ValueError(f"NPZ_ARRAY_KEY_NOT_FOUND: {array_key}")
                used_key = array_key
                arr = archive[array_key]
            else:
                candidates = [k for k in archive.files if k not in ("spacing", "origin", "axis_order", "unit") and archive[k].ndim == 3]
                if len(candidates) == 1:
                    used_key = candidates[0]
                    arr = archive[used_key]
                elif not candidates:
                    raise ValueError("VOXEL_VOLUME_3D_REQUIRED")
                else:
                    raise ValueError(f"NPZ_ARRAY_KEY_REQUIRED: {', '.join(candidates)}")
    elif suffix in (".tif", ".tiff"):
        from tifffile import imread

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
    pitch: float | None = None,
    axis_order: str = "auto",
    threshold: float = 0.0,
    array_key: str | None = None,
    label: str | None = None,
    spacing=None,
    labels=None,
    unit=None,
    metadata=None,
    resample_pitch=None,
) -> dict:
    """Convert a 3D binary/numeric volume into a SHA-bound native entity package."""
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(float(threshold)):
        raise ValueError("FINITE_VOXEL_THRESHOLD_REQUIRED")
    threshold_val = float(threshold)

    src = Path(source).resolve()
    folder = Path(output_dir).resolve()
    if folder.exists() and any(folder.iterdir()):
        raise ValueError("IMPORT_OUTPUT_EXISTS")

    settings = resolve_volume_settings(src, pitch=pitch, spacing=spacing, axis_order=axis_order,
                                       unit=unit, metadata=metadata, resample_pitch=resample_pitch)
    pitch_val, order = settings["pitch"], settings["axis_order"]
    if labels is not None:
        if not isinstance(labels, (list, tuple)) or not labels or any(isinstance(v, bool) or not isinstance(v, int) for v in labels):
            raise ValueError("NONEMPTY_INTEGER_VOXEL_LABELS_REQUIRED")
        labels = sorted(set(labels))
    arr, used_key = load_3d_voxel_array(src, array_key=array_key)
    mask = np.isin(arr, labels) if labels is not None else (arr if arr.dtype == np.bool_ else (arr > threshold_val))
    source_occupied = int(np.count_nonzero(mask))
    if order == "zyx":
        mask = mask.transpose(2, 1, 0)
    if settings["resampled"]:
        from scipy.ndimage import affine_transform
        ratio = pitch_val / np.asarray(settings["spacing_xyz"])
        shape_xyz = tuple(int(math.ceil(n * s / pitch_val)) for n, s in zip(mask.shape, settings["spacing_xyz"]))
        mask = affine_transform(mask, np.diag(ratio), offset=(ratio - 1) / 2,
                                output_shape=shape_xyz, order=0, mode="grid-constant", cval=0, prefilter=False)
    coords = np.argwhere(mask)
    if not len(coords):
        raise ValueError("EMPTY_VOXEL_VOLUME")
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
            threshold=threshold_val,
            labels=labels,
            raw_shape=[int(x) for x in arr.shape],
            source_occupied_voxels=source_occupied,
            **settings,
        ),
        geometry="imported 3D binary voxel cube union; coordinates=(integer_index+0.5)*pitch",
        no_closing_or_postprocessing=not settings["resampled"],
        resampling=dict(applied=settings["resampled"], method="nearest_neighbor_cell_centers" if settings["resampled"] else None,
                        original_volume=source_occupied * math.prod(settings["spacing_xyz"]),
                        output_volume=len(cells) * pitch_val**3,
                        source_data_unchanged=True),
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
        settings=settings,
    )
