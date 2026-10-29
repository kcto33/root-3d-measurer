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
    from root_measurer import morphology
    from root_model_delivery_exports import file_reference
    before=file_reference(result.measurement_path)
    changed=morphology(result.measurement_path,vertical_axis='x',down_sign=1)
    assert changed['settings']['vertical_axis']=='x'
    assert changed['measurement']==before==file_reference(result.measurement_path)
    assert changed['source_binding']['cells_sha256']
    assert parameters['morphology']['settings']['vertical_axis']=='z'
    assert parameters['morphology']['settings']['down_sign']==-1
    traits={m['parameter']:m['value'] for m in changed['metrics']}
    assert math.isclose(traits['root_depth'],20)
    assert traits['fractal_dimension_3d'] is not None
