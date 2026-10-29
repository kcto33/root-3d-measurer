"""A real finite entity traverses the installed supervisor, worker and queries."""
import math
from root_measurer import measure, query_length, query_radius


def test_generated_entity_complete_measurement_and_fresh_section(control_workspace):
    from build_native_measurement_acceptance import build
    inputs = control_workspace / "snapshot/d"
    build("straight", inputs)
    result = measure(inputs / "straight", workspace=control_workspace / "results/api", scale_mm=100)
    parameters = result.parameters
    values = {r["parameter"]: r["value"] for r in parameters["metrics"]}
    assert parameters["complete_model_length_available"]
    assert math.isclose(values["model_length"], 120, rel_tol=0, abs_tol=1e-9)
    assert parameters["geometry"]["kind"] == "EXACT_NATIVE_CUBE_UNION"
    segment = query_length(result.measurement_path, [{"path_id":"path00000", "start":.12345, "end":.78901}])
    assert math.isclose(segment["length"], 66.556)
    section = query_radius(result.measurement_path, "path00000", .53217)
    assert math.isclose(section["area"], 320)
    assert math.isclose(section["radius"], math.sqrt(320 / math.pi))
    assert not section["stored_samples_interpolated"]
