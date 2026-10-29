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


def measure(model, *, workspace="runs/default", scale_mm=None, max_memory_mb=8192,
            mesh_target_pitch=None, mesh_target_reason="", vertical_axis="z", down_sign=-1):
    activate()
    from root_model_morphology import settings
    settings(dict(vertical_axis=vertical_axis, down_sign=down_sign))
    source = Path(model).expanduser().resolve()
    if source.is_dir():
        source /= "entity_prediction.json"
    if not source.is_file():
        raise FileNotFoundError(source)
    if scale_mm is not None and (isinstance(scale_mm, bool) or not math.isfinite(scale_mm) or scale_mm <= 0):
        raise ValueError("POSITIVE_FINITE_MM_SCALE_REQUIRED")
    root = initialize_workspace(workspace)
    prefix = unused_run(root) / "model"
    mesh = source.suffix.lower() == ".ply"
    args = [sys.executable, "-B", str(script("measure_model_delivery.py")),
            "--mesh" if mesh else "--entity", str(source), "--output-prefix", str(prefix),
            "--max-memory-mb", str(max_memory_mb), "--vertical-axis", vertical_axis, "--down-sign", str(down_sign)]
    if scale_mm is not None:
        args += ["--scale-mm", str(scale_mm)]
    if mesh:
        args += ["--mesh-workflow", "--mesh-target-reason", mesh_target_reason]
        if mesh_target_pitch is not None:
            args += ["--mesh-target-pitch", str(mesh_target_pitch)]
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
