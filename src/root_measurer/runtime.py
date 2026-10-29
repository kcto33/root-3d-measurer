"""Installed source resources and user-owned writable workspaces."""
import json
import sys
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent
MARKER = "ROOT_MEASURER_WORKSPACE.json"


def activate():
    # The verified geometry modules retain their import and inspectable source
    # identities. Internal imports are not part of the supported public API.
    for folder in [PACKAGE / "scripts", PACKAGE / "src"]:
        name = str(folder)
        if name not in sys.path:
            sys.path.insert(0, name)


def initialize_workspace(folder):
    root = Path(folder).expanduser().resolve()
    if root == PACKAGE or root.is_relative_to(PACKAGE):
        raise ValueError("WORKSPACE_MUST_BE_OUTSIDE_INSTALLED_PACKAGE")
    root.mkdir(parents=True, exist_ok=True)
    marker = root / MARKER
    if marker.exists():
        value = json.loads(marker.read_text(encoding="utf-8"))
        if value.get("schema") != "root-measurer-workspace-v001":
            raise ValueError("UNRECOGNIZED_WORKSPACE_MARKER")
    else:
        if any(root.iterdir()):
            raise ValueError("USE_EMPTY_FOLDER_OR_EXISTING_MEASURER_WORKSPACE")
        with marker.open("x", encoding="utf-8") as stream:
            from . import __version__
            json.dump({"schema": "root-measurer-workspace-v001", "version": __version__}, stream)
    for relative in ["results/model_delivery_local", "analysis", "snapshot"]:
        (root / relative).mkdir(parents=True, exist_ok=True)
    log = root / "EXPERIMENT_LOG.md"
    if not log.exists():
        with log.open("x", encoding="utf-8") as stream:
            stream.write("# Root 3D Measurer workspace\n\nGenerated inputs and measurements retain source hashes; output prefixes are never overwritten.\n")
    return root


def unused_run(root):
    base = Path(root) / "results/model_delivery_local"
    number = 1
    while True:
        folder = base / f"run{number:03d}"
        try:
            folder.mkdir()
            return folder
        except FileExistsError:
            number += 1


def script(name):
    return PACKAGE / "scripts" / name
