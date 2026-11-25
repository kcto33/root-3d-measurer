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


def test_section_workers_cli_and_validation(tmp_path, monkeypatch):
    from root_measurer import cli
    from root_measurer.runtime import activate

    args = cli.parser().parse_args(["measure", str(tmp_path / "m"), "--section-workers", "4"])
    assert args.section_workers == 4
    batch_args = cli.parser().parse_args(["batch", str(tmp_path / "b.json"), "--section-workers", "3"])
    assert batch_args.section_workers == 3
    assert cli.parser().parse_args(["measure", str(tmp_path / "m")]).section_workers == 1
    assert cli.parser().parse_args(["batch", str(tmp_path / "b.json")]).section_workers == 1

    source = tmp_path / "entity_prediction.json"
    source.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="POSITIVE_SECTION_WORKERS_REQUIRED"):
        measure(source, workspace=tmp_path / "workspace", section_workers=0)
    with pytest.raises(ValueError, match="POSITIVE_SECTION_WORKERS_REQUIRED"):
        measure(source, workspace=tmp_path / "workspace", section_workers=True)
    assert cli.main(["batch", str(tmp_path / "b.json"), "--workspace", str(tmp_path / "ws"), "--section-workers", "0"]) == 1

    captured = {}
    def fake_measure(*a, **kw):
        captured["section_workers"] = kw.get("section_workers")
        class _Dummy:
            prefix = tmp_path / "out/model"
            parameters = {}
        return _Dummy()
    monkeypatch.setattr(cli, "measure", fake_measure)
    assert cli.main(["measure", str(source), "--section-workers", "2"]) == 0
    assert captured["section_workers"] == 2

    activate()
    from root_model_delivery_service import _parse_section_workers
    assert _parse_section_workers(None) == 1
    assert _parse_section_workers(4) == 4
    for bad in (0, -2, True, "2", 1.5):
        with pytest.raises(ValueError, match="POSITIVE_SECTION_WORKERS_REQUIRED"):
            _parse_section_workers(bad)
