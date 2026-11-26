"""Shared public measurement options for single, batch and HTTP inputs."""
import math
import sys
from pathlib import Path

DEFAULTS = dict(scale_mm=None, max_memory_mb=8192, section_workers=1,
                mesh_target_pitch=None, mesh_target_reason="", vertical_axis="z", down_sign=-1,
                no_axis_candidate=False, voxel_pitch=None, voxel_axis_order="auto", voxel_threshold=0.0,
                voxel_array_key=None, voxel_spacing=None, voxel_labels=None, voxel_unit=None,
                voxel_metadata=None, voxel_resample_pitch=None)


def add_volume_options(parser, *, short=False):
    prefix = "" if short else "voxel-"
    parser.add_argument("--" + prefix + "pitch", type=float, help="Isotropic source pitch; optional when metadata supplies spacing")
    parser.add_argument("--" + prefix + "axis-order", choices=["auto", "xyz", "zyx"], default="auto")
    parser.add_argument("--" + prefix + "threshold", type=float, default=0.0)
    parser.add_argument("--" + prefix + "array-key")
    parser.add_argument("--" + prefix + "spacing", type=float, nargs=3, metavar=("SX", "SY", "SZ"), help="Source spacing in XYZ order")
    parser.add_argument("--" + prefix + "labels", type=int, nargs="+", help="Measure only these segmentation values")
    parser.add_argument("--" + prefix + "unit", choices=["mm", "cm", "m", "um", "nm"])
    parser.add_argument("--" + prefix + "metadata", type=Path, help="Acquisition metadata JSON")
    parser.add_argument("--" + prefix + "resample-pitch", type=float, help="Explicit nearest-neighbor isotropic target pitch")


def measurement_command(source, prefix, options=None, cancel_file=None):
    from root_model_morphology import settings
    from root_voxel_import import _validate_pitch, _validate_axis_order, _vector, is_voxel_volume_path
    options = dict(options or {})
    if set(options) - DEFAULTS.keys():
        raise ValueError("UNKNOWN_MEASUREMENT_OPTION: " + ", ".join(sorted(set(options) - DEFAULTS.keys())))
    config = dict(DEFAULTS, **options)
    settings({k: config[k] for k in ("vertical_axis", "down_sign")})
    for key in ("max_memory_mb", "section_workers"):
        if isinstance(config[key], bool) or not isinstance(config[key], int) or config[key] < 1:
            raise ValueError("POSITIVE_" + ("WORKER_MEMORY_BUDGET" if key == "max_memory_mb" else "SECTION_WORKERS") + "_REQUIRED")
    scale = config["scale_mm"]
    if scale is not None and (isinstance(scale, bool) or not math.isfinite(float(scale)) or float(scale) <= 0):
        raise ValueError("POSITIVE_FINITE_MM_SCALE_REQUIRED")
    source = Path(source).expanduser().resolve()
    if source.is_dir():
        source /= "entity_prediction.json"
    if not source.is_file():
        raise FileNotFoundError(source)
    volume = is_voxel_volume_path(source)
    mesh = source.suffix.lower() == ".ply"
    if not volume and not mesh and source.suffix.lower() != ".json":
        raise ValueError("UNSUPPORTED_MEASUREMENT_INPUT")
    command = [sys.executable, "-B", str(Path(__file__).resolve().parents[1] / "scripts/measure_model_delivery.py"),
               "--volume" if volume else ("--mesh" if mesh else "--entity"), str(source), "--output-prefix", str(prefix),
               "--max-memory-mb", str(config["max_memory_mb"]), "--section-workers", str(config["section_workers"]),
               "--vertical-axis", config["vertical_axis"], "--down-sign", str(config["down_sign"])]
    if scale is not None:
        command += ["--scale-mm", str(scale)]
    if config["no_axis_candidate"]:
        command += ["--no-axis-candidate"]
    if mesh:
        command += ["--mesh-workflow", "--mesh-target-reason", config["mesh_target_reason"]]
        if config["mesh_target_pitch"] is not None:
            _validate_pitch(config["mesh_target_pitch"])
            command += ["--mesh-target-pitch", str(config["mesh_target_pitch"])]
    if volume:
        _validate_axis_order(config["voxel_axis_order"])
        if config["voxel_pitch"] is not None:
            _validate_pitch(config["voxel_pitch"])
        if config["voxel_spacing"] is not None:
            _vector(config["voxel_spacing"], "SPACING", positive=True)
        if config["voxel_pitch"] is not None and config["voxel_spacing"] is not None:
            raise ValueError("CHOOSE_VOXEL_PITCH_OR_SPACING")
        if source.suffix.lower() == ".npy" and all(config[k] is None for k in ("voxel_pitch", "voxel_spacing", "voxel_metadata")):
            raise ValueError("POSITIVE_FINITE_VOXEL_PITCH_REQUIRED")
        threshold = config["voxel_threshold"]
        if isinstance(threshold, bool) or not math.isfinite(float(threshold)):
            raise ValueError("FINITE_VOXEL_THRESHOLD_REQUIRED")
        labels = config["voxel_labels"]
        if labels is not None and (not isinstance(labels, (list, tuple)) or not labels or any(isinstance(v, bool) or not isinstance(v, int) for v in labels)):
            raise ValueError("NONEMPTY_INTEGER_VOXEL_LABELS_REQUIRED")
        if config["voxel_resample_pitch"] is not None:
            _validate_pitch(config["voxel_resample_pitch"])
        for key in DEFAULTS:
            if not key.startswith("voxel_") or config[key] is None:
                continue
            value = config[key]
            command += ["--" + key.replace("_", "-")]
            command += [str(v) for v in value] if isinstance(value, (list, tuple)) else [str(value)]
    if cancel_file is not None:
        command += ["--cancel-file", str(cancel_file)]
    return command, config, source
