"""Portable local CLI; no historical experiment or cloud account is required."""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from . import __version__
from .api import MeasurementError, measure, query_length, query_radius
from .runtime import activate, initialize_workspace, script


def parser():
    ap = argparse.ArgumentParser(prog="root-measurer", description="根系三维测量器 / Root 3D Measurer")
    ap.add_argument("--version", action="version", version=__version__)
    commands = ap.add_subparsers(dest="command", required=True)
    gui = commands.add_parser("gui", help="启动本地浏览器工作台")
    gui.add_argument("--workspace", type=Path, default=Path("runs/default"))
    gui.add_argument("--input", type=Path, action="append", default=[])
    gui.add_argument("--port", type=int, default=0)
    gui.add_argument("--open-browser", action="store_true")
    one = commands.add_parser("measure", help="测量一个原生实体或三角面 PLY")
    one.add_argument("model", type=Path)
    one.add_argument("--workspace", type=Path, default=Path("runs/default"))
    one.add_argument("--scale-mm", type=float)
    one.add_argument("--max-memory-mb", type=int, default=8192)
    one.add_argument("--mesh-target-pitch", type=float)
    one.add_argument("--mesh-target-reason", default="")
    batch = commands.add_parser("batch", help="按清单顺序批量测量原生实体")
    batch.add_argument("manifest", type=Path)
    batch.add_argument("--workspace", type=Path, default=Path("runs/default"))
    batch.add_argument("--scale-mm", type=float)
    batch.add_argument("--max-memory-mb", type=int, default=8192)
    demo = commands.add_parser("demo", help="现场生成小型合成实体；无需下载数据")
    demo.add_argument("--output", type=Path, required=True)
    demo.add_argument("--case", choices=["straight", "curve", "branching_y", "loop"], default="straight")
    length = commands.add_parser("length", help="显式根段区间的三维弧长")
    length.add_argument("measurement", type=Path)
    length.add_argument("--path", required=True)
    length.add_argument("--start", type=float, required=True)
    length.add_argument("--end", type=float, required=True)
    radius = commands.add_parser("radius", help="连续轴位置的新原实体截面半径")
    radius.add_argument("measurement", type=Path)
    radius.add_argument("--path", required=True)
    radius.add_argument("--arc", type=float, required=True)
    radius.add_argument("--normal", type=float, nargs=3)
    score = commands.add_parser("score", help="独立参考评分；与数值预测分开")
    score.add_argument("measurement", type=Path)
    score.add_argument("--reference", type=Path, required=True)
    score.add_argument("--output", type=Path, required=True)
    return ap


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        activate()
        if args.command == "gui":
            from root_model_delivery_service import serve_workspace
            root = initialize_workspace(args.workspace)
            serve_workspace(root, args.port, args.open_browser, input_paths=args.input)
        elif args.command == "measure":
            result = measure(args.model, workspace=args.workspace, scale_mm=args.scale_mm,
                max_memory_mb=args.max_memory_mb, mesh_target_pitch=args.mesh_target_pitch,
                mesh_target_reason=args.mesh_target_reason)
            print(json.dumps({"prefix": str(result.prefix), "parameters": result.parameters}, ensure_ascii=False))
        elif args.command == "batch":
            root = initialize_workspace(args.workspace)
            base = root / "results/model_delivery_local"
            number = 1
            while (base / f"run{number:03d}").exists():
                number += 1
            cmd = [sys.executable, "-B", str(script("measure_model_batch.py")), "--manifest",
                   str(args.manifest.resolve()), "--output-dir", str(base / f"run{number:03d}"),
                   "--max-memory-mb", str(args.max_memory_mb)]
            if args.scale_mm is not None:
                cmd += ["--scale-mm", str(args.scale_mm)]
            return subprocess.call(cmd, env=dict(os.environ, PYTHONUTF8="1", PYTHONDONTWRITEBYTECODE="1"))
        elif args.command == "demo":
            from build_native_measurement_acceptance import build
            print(json.dumps(build(args.case, args.output.resolve()), ensure_ascii=False))
        elif args.command == "length":
            print(json.dumps(query_length(args.measurement, [{"path_id": args.path, "start": args.start, "end": args.end}]), ensure_ascii=False))
        elif args.command == "radius":
            print(json.dumps(query_radius(args.measurement, args.path, args.arc, normal=args.normal), ensure_ascii=False))
        elif args.command == "score":
            from score_native_model_delivery import score
            from root_model_delivery_exports import file_reference, write_json
            output = args.output.resolve()
            output.parent.mkdir(parents=True, exist_ok=True)
            value = score(args.measurement.resolve(), args.reference.resolve(), output,
                          reference_engine="native_plane_geos_connected_v004")
            measured = args.measurement.resolve()
            if not measured.name.endswith("_measurement.json"):
                raise ValueError("MEASUREMENT_FILENAME_REQUIRED")
            sidecar = measured.with_name(measured.name.removesuffix("_measurement.json") + "_connected_section_validation.json")
            write_json(sidecar, dict(measurement=file_reference(measured), score=file_reference(output),
                radius_quantity=value["radius_quantity"], precision_pass=value["precision_pass"], metrics=value["metrics"]))
        return 0
    except (MeasurementError, ValueError, OSError, KeyError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
