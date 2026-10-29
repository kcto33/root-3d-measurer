import json
from pathlib import Path
import pytest

from root_measurer import measure, query_length
from root_measurer.runtime import initialize_workspace


def test_workspace_refuses_to_adopt_an_unrelated_folder(tmp_path):
    original = tmp_path / "original.txt"
    original.write_text("preserve", encoding="utf-8")
    with pytest.raises(ValueError, match="EMPTY_FOLDER"):
        initialize_workspace(tmp_path)
    assert original.read_text(encoding="utf-8") == "preserve"


def test_workspace_can_be_reopened_and_keeps_existing_result(tmp_path):
    root = initialize_workspace(tmp_path / "workspace")
    result = root / "results/retained.json"
    result.write_text('{"value":7}', encoding="utf-8")
    assert initialize_workspace(root) == root
    assert json.loads(result.read_text(encoding="utf-8")) == {"value": 7}


def test_public_interval_keeps_three_dimensional_arc_and_source_binding(tmp_path):
    path = tmp_path / "curve_measurement.json"
    path.write_text(json.dumps(dict(candidate_centerline=dict(
        points=[[0,0,0],[1,0,0],[1,1,0]], adjacency=[[1],[0,2],[1]]))), encoding="utf-8")
    result = query_length(path, [{"path_id":"path00000", "start":.1, "end":1.7}])
    assert result["length"] == pytest.approx(1.6)
    assert not result["graph_route_searched"]
    assert result["measurement"]["bytes"] == path.stat().st_size


def test_invalid_scale_does_not_start_or_allocate_a_run(tmp_path):
    source = tmp_path / "entity_prediction.json"
    source.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="POSITIVE_FINITE"):
        measure(source, workspace=tmp_path / "workspace", scale_mm=-1)
    assert not (tmp_path / "workspace").exists()
