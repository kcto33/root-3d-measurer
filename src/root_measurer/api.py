"""Stable API around supervised measurement and explicit user-selected intervals."""
import json
import math
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from .runtime import activate, initialize_workspace, script, unused_run


class MeasurementError(RuntimeError):
    """A failed/blocked run retains its prefix and all available evidence."""
    def __init__(self, prefix, exit_code):
        self.prefix = Path(prefix)
        self.exit_code = exit_code
        super().__init__(f"Measurement failed or blocked (exit {exit_code}); inspect {self.prefix}_supervision.json and process log")


@dataclass(frozen=True)
class Measurement:
    prefix: Path

    @property
    def measurement_path(self):
        return self.prefix.with_name(self.prefix.name + "_measurement.json")

    @property
    def parameters(self):
        path = self.prefix.with_name(self.prefix.name + "_parameters.json")
        return json.loads(path.read_text(encoding="utf-8"))


@dataclass(frozen=True)
class BatchMeasurement:
    folder: Path

    @property
    def summary_path(self):
        return self.folder / "batch_summary.json"

    @property
    def summary(self):
        return json.loads(self.summary_path.read_text(encoding="utf-8"))


def batch(inputs, *, workspace="runs/default", defaults=None):
    """Measure a mixed input list or a JSON manifest; preserve every failure row."""
    activate()
    from root_measurement_inputs import DEFAULTS
    base = Path.cwd()
    if isinstance(inputs, (str, Path)):
        manifest = Path(inputs).expanduser().resolve()
        value = json.loads(manifest.read_text(encoding="utf-8-sig"))
        base = manifest.parent
    else:
        value = inputs if isinstance(inputs, dict) else dict(inputs=inputs)
    if not isinstance(value, dict) or not isinstance(value.get("inputs"), list) or not 1 <= len(value["inputs"]) <= 1000:
        raise ValueError("BATCH_SIZE_1_TO_1000_REQUIRED")
    options = dict(defaults or {}, **value.get("defaults", {}))
    if set(options) - DEFAULTS.keys():
        raise ValueError("UNKNOWN_BATCH_DEFAULT_OPTION")
    root = initialize_workspace(workspace)
    folder = unused_run(root)
    requests = root / "snapshot/batch_requests"
    requests.mkdir(exist_ok=True)
    request = requests / (folder.name + ".json")
    with request.open("x", encoding="utf-8") as stream:
        json.dump(dict(inputs=value["inputs"], defaults=options), stream, ensure_ascii=False, default=str)
    command = [sys.executable, "-B", str(script("measure_model_batch.py")), "--manifest", str(request),
               "--manifest-base", str(base), "--output-dir", str(folder)]
    completed = subprocess.run(command, env=dict(os.environ, PYTHONUTF8="1", PYTHONDONTWRITEBYTECODE="1"))
    if completed.returncode:
        raise RuntimeError("BATCH_RUN_FAILED_SEE_LOGS: " + str(folder))
    return BatchMeasurement(folder)


def import_volume(source, output_dir, *, pitch=None, axis_order="auto", threshold=0.0, array_key=None, label=None,
                  spacing=None, labels=None, unit=None, metadata=None, resample_pitch=None):
    """Import a 3D binary/numeric voxel volume (.npy/.npz/.tif/.tiff) as a SHA-bound native solid."""
    activate()
    from root_voxel_import import import_voxel_volume
    return import_voxel_volume(
        source,
        output_dir,
        pitch=pitch,
        axis_order=axis_order,
        threshold=threshold,
        array_key=array_key,
        label=label,
        spacing=spacing, labels=labels, unit=unit, metadata=metadata, resample_pitch=resample_pitch,
    )


def measure(model, *, workspace="runs/default", scale_mm=None, max_memory_mb=8192,
            mesh_target_pitch=None, mesh_target_reason="", vertical_axis="z", down_sign=-1,
            no_axis_candidate=False, section_workers=1, voxel_pitch=None, voxel_axis_order="auto",
            voxel_threshold=0.0, voxel_array_key=None, voxel_spacing=None, voxel_labels=None,
            voxel_unit=None, voxel_metadata=None, voxel_resample_pitch=None):
    activate()
    from root_measurement_inputs import DEFAULTS, measurement_command
    options = {k: v for k, v in locals().copy().items() if k in DEFAULTS}
    # Validate before allocating output; acquisition metadata and arrays are read by the worker.
    args, _, _ = measurement_command(model, "__pending__", options)
    root = initialize_workspace(workspace)
    run_dir = unused_run(root)
    prefix = run_dir / "model"
    args[args.index("--output-prefix") + 1] = str(prefix)
    result = subprocess.run(args, env=dict(os.environ, PYTHONUTF8="1", PYTHONDONTWRITEBYTECODE="1"))
    if result.returncode:
        raise MeasurementError(prefix, result.returncode)
    return Measurement(prefix)


def _interaction(measurement):
    activate()
    from root_model_interaction import ModelInteraction
    from root_native_entity_measurement import check_entity_unchanged
    report = json.loads(Path(measurement).read_text(encoding="utf-8"))
    if report.get("entity_binding"):
        check_entity_unchanged(report["entity_binding"])
    return ModelInteraction(report)


def query_length(measurement, intervals):
    """Intervals contain explicit path_id/start/end in model-unit arc distances."""
    result = _interaction(measurement).length(intervals)
    from root_model_delivery_exports import file_reference
    return dict(result, measurement=file_reference(measurement))


def query_radius(measurement, path_id, arc, *, normal=None):
    """Fresh original-solid section; ambiguous/contact sections remain unavailable."""
    request = {"path_id": path_id, "arc": arc}
    if normal is not None:
        request["normal"] = normal
    result = _interaction(measurement).radius(request)
    from root_model_delivery_exports import file_reference
    return dict(result, measurement=file_reference(measurement))


def morphology(measurement, *, vertical_axis="z", down_sign=-1):
    """Additive analysis: verify bound original input; never rewrite measurements."""
    activate()
    from root_model_morphology import build_morphology, settings
    from root_model_delivery_exports import file_reference
    from root_native_entity_measurement import load_entity, check_entity_unchanged
    options = dict(vertical_axis=vertical_axis, down_sign=down_sign)
    settings(options)
    measured = Path(measurement).resolve()
    before = file_reference(measured)
    report = json.loads(measured.read_text(encoding="utf-8"))
    solid = mesh = None
    geometry = {}
    if report.get("entity_binding"):
        binding = report["entity_binding"]
        check_entity_unchanged(binding)
        solid, actual, _ = load_entity(binding["metadata"]["path"], stage=binding["stage"])
        if actual["cells_sha256"] != binding["cells_sha256"]:
            raise ValueError("ORIGINAL_ENTITY_CHANGED")
        from root_reconstruction_measurement import boundary_faces
        geometry["surface_area"] = len(boundary_faces(solid)[0])*solid.pitch**2
    elif report.get("input_file_binding"):
        from root_triangle_measurement_input import read_triangle_ply_bytes
        binding = report["input_file_binding"]
        source = Path(binding["path"])
        actual = file_reference(source)
        if any(actual[k] != binding[k] for k in ("sha256", "bytes")):
            raise ValueError("ORIGINAL_MESH_CHANGED")
        mesh = read_triangle_ply_bytes(source.read_bytes())
    value = build_morphology(report, solid=solid, mesh=mesh, options=options, geometry=geometry)
    if report.get("entity_binding"):
        check_entity_unchanged(report["entity_binding"])
    elif mesh is not None and file_reference(source) != actual:
        raise ValueError("ORIGINAL_MESH_CHANGED")
    if file_reference(measured) != before:
        raise ValueError("MEASUREMENT_CHANGED")
    return dict(value, measurement=before, source_binding=report.get("entity_binding") or report.get("input_file_binding"))


def phenotype(measurement, *, options=None, annotation_profile=None):
    """Analyze distributions/annotated root axes without changing source geometry.

    Linear settings use mm when the measurement is calibrated, otherwise model
    units. Imported annotation profiles must match the measurement SHA256.
    """
    activate()
    from root_model_phenotyping import analyze_measurement
    return analyze_measurement(measurement, options=options, annotation_profile=annotation_profile)
