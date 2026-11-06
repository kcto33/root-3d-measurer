"""Windows desktop shell: a native window around the local workspace server.

Run with the installed interpreter: ``pythonw -m root_measurer.desktop``.
The measurement pipeline is untouched; this module only starts the same
loopback server used by ``root-measurer gui`` inside a pywebview window and
bridges native file dialogs, save-as and folder opening to the page.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.request
from pathlib import Path

APP_NAME = "Root3DMeasurer"
WINDOW_TITLE = "根系三维测量器"
DEMO_CASES = ("straight", "curve", "branching_y", "loop")
RECENT_LIMIT = 20


def app_dir():
    base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
    return Path(base) / APP_NAME


def setup_logging():
    # Under pythonw.exe sys.stdout/stderr are None; every print would raise.
    folder = app_dir()
    folder.mkdir(parents=True, exist_ok=True)
    stream = (folder / "app.log").open("a", encoding="utf-8", buffering=1)
    sys.stdout = stream
    sys.stderr = stream


def resolve_workspace(explicit):
    if explicit:
        return Path(explicit).expanduser().resolve()
    settings = app_dir() / "settings.json"
    if settings.is_file():
        try:
            value = json.loads(settings.read_text(encoding="utf-8-sig")).get("workspace")
            if value:
                return Path(value).expanduser().resolve()
        except (ValueError, OSError):
            pass
    return app_dir() / "workspace"


def load_recent():
    path = app_dir() / "recent_entities.json"
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
        return [str(item) for item in value] if isinstance(value, list) else []
    except (ValueError, OSError):
        return []


def save_recent(paths):
    path = app_dir() / "recent_entities.json"
    try:
        path.write_text(json.dumps(paths, ensure_ascii=False, indent=1), encoding="utf-8")
    except OSError:
        pass


class DesktopBridge:
    """Exposed to the desktop page as ``window.pywebview.api.*``."""

    def __init__(self, base_url, workspace):
        self.base_url = base_url
        self.workspace = Path(workspace)

    def _window(self):
        import webview

        return webview.windows[0]

    def pick_file(self):
        import webview

        result = self._window().create_file_dialog(
            webview.OPEN_DIALOG, allow_multiple=False,
            file_types=("实体或三角网格 (*.json;*.ply)", "所有文件 (*.*)"))
        if not result:
            return None
        return str(result[0] if isinstance(result, (list, tuple)) else result)

    def save_file(self, url, suggested_name="download"):
        # Only same-server paths of the form "/result/..." may be saved.
        if not isinstance(url, str) or not url.startswith("/") or "//" in url:
            raise ValueError("BAD_DOWNLOAD_URL")
        with urllib.request.urlopen(self.base_url.rstrip("/") + url, timeout=300) as response:
            data = response.read()
            matched = re.search(r'filename="([^"]+)"', response.headers.get("Content-Disposition", ""))
            name = matched[1] if matched else (suggested_name or "download")
        target = self._save_dialog(name)
        if not target:
            return None
        Path(target).write_bytes(data)
        print("saved", target, flush=True)
        return str(target)

    def save_text(self, name, text):
        target = self._save_dialog(name)
        if not target:
            return None
        Path(target).write_text(str(text), encoding="utf-8")
        print("saved", target, flush=True)
        return str(target)

    def _save_dialog(self, name):
        import webview

        result = self._window().create_file_dialog(
            webview.SAVE_DIALOG, directory=str(Path.home() / "Downloads"), save_filename=name)
        if not result:
            return None
        return str(result if isinstance(result, str) else result[0])

    def open_workspace(self):
        os.startfile(str(self.workspace))  # noqa: S606 - intended shell open
        return True

    def recent_entities(self):
        return load_recent()

    def remember_entity(self, path):
        try:
            value = str(Path(str(path)).resolve())
        except OSError:
            return False
        items = [item for item in load_recent() if item != value]
        items.insert(0, value)
        save_recent(items[:RECENT_LIMIT])
        return True

    def forget_entity(self, path):
        value = str(Path(str(path)).resolve())
        save_recent([item for item in load_recent() if item != value])
        return True

    def generate_demo(self, case="straight"):
        if case not in DEMO_CASES:
            raise ValueError("UNKNOWN_DEMO_CASE")
        output = self.workspace / "demo_inputs"
        command = [sys.executable, "-m", "root_measurer", "demo", "--output", str(output), "--case", case]
        completed = subprocess.run(
            command, capture_output=True, text=True, encoding="utf-8", errors="replace",
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
        entity = output / case / "entity_prediction.json"
        if completed.returncode or not entity.is_file():
            detail = (completed.stderr or completed.stdout or "").strip().splitlines()
            raise RuntimeError(detail[-1] if detail else "DEMO_GENERATION_FAILED")
        return {"entity_path": str(entity.parent)}


def cancel_running_jobs(base_url, token):
    request = urllib.request.Request(
        base_url + "api/cancel-all", data=b"{}", method="POST",
        headers={"Content-Type": "application/json", "X-Root-Local-Token": token})
    with urllib.request.urlopen(request, timeout=5) as response:
        print("cancel-all:", response.read().decode("utf-8", "replace"), flush=True)


def main(argv=None):
    setup_logging()
    parser = argparse.ArgumentParser(description="根系三维测量器桌面版")
    parser.add_argument("--workspace", default=None, help="工作区目录，默认 %%LOCALAPPDATA%%\\Root3DMeasurer\\workspace")
    parser.add_argument("--port", type=int, default=0, help="服务端口，默认随机空闲端口")
    args = parser.parse_args(argv)

    from . import __version__
    print(time.strftime("%Y-%m-%d %H:%M:%S"), "desktop starting, version", __version__, flush=True)

    workspace = resolve_workspace(args.workspace)
    from .runtime import activate, initialize_workspace

    experiment = initialize_workspace(workspace)
    activate()
    from root_model_delivery_service import start_workspace_server

    server, url, token = start_workspace_server(experiment, args.port)
    threading.Thread(target=server.serve_forever, daemon=True, name="workspace-server").start()
    print("serving", url, flush=True)

    try:
        import webview
    except ImportError:
        print("pywebview not installed; opening the default browser instead", flush=True)
        import webbrowser

        webbrowser.open(url + "desktop")
        try:
            while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            pass
        finally:
            server.shutdown()
        return

    webview.settings["ALLOW_DOWNLOADS"] = True
    bridge = DesktopBridge(url, experiment)
    webview.create_window(
        WINDOW_TITLE, url + "desktop", js_api=bridge,
        width=1280, height=840, min_size=(1000, 640), background_color="#f4f6f8")
    try:
        webview.start()
    except Exception as exc:  # noqa: BLE001 - any shell failure must not kill the server silently
        print("webview failed:", exc, flush=True)
        import webbrowser

        webbrowser.open(url + "desktop")
        try:
            while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            pass
    finally:
        print("window closed", flush=True)
        try:
            cancel_running_jobs(url, token)
        except OSError as exc:
            print("cancel-all failed:", exc, flush=True)
        server.shutdown()
        print("desktop stopped", flush=True)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except BaseException:
        import traceback

        traceback.print_exc(file=sys.stderr)
        raise SystemExit(1)
