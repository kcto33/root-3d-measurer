"""Run geometry tests against packaged sources with isolated, disposable test inputs."""
import os
from pathlib import Path
import pytest
from root_measurer.runtime import activate, initialize_workspace

activate()


@pytest.fixture(scope="session", autouse=True)
def control_workspace(tmp_path_factory):
    root = initialize_workspace(tmp_path_factory.mktemp("w"))
    controls = root / "snapshot/c"
    controls.mkdir()
    old = os.environ.get("NATIVE_ENTITY_CONTROL_OUTPUT")
    os.environ["NATIVE_ENTITY_CONTROL_OUTPUT"] = str(controls)
    yield root
    if old is None:
        os.environ.pop("NATIVE_ENTITY_CONTROL_OUTPUT", None)
    else:
        os.environ["NATIVE_ENTITY_CONTROL_OUTPUT"] = old
