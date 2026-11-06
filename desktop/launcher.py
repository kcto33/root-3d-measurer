"""Tiny windowed launcher for the portable desktop build.

Compiled by PyInstaller into ``Root3DMeasurer.exe``; its only job is to start
``runtime\\pythonw.exe -m root_measurer.desktop`` next to itself and exit.
"""
import subprocess
import sys
from pathlib import Path


def main():
    here = Path(sys.executable).resolve().parent if getattr(sys, "frozen", False) else Path(__file__).resolve().parent
    pythonw = here / "runtime" / "pythonw.exe"
    if not pythonw.is_file():
        import ctypes

        ctypes.windll.user32.MessageBoxW(
            None, "未找到 runtime\\pythonw.exe。\n请保持 Root3DMeasurer 整个文件夹完整，不要单独移动 exe。",
            "根系三维测量器", 0x10)
        return 1
    subprocess.Popen(
        [str(pythonw), "-m", "root_measurer.desktop"],
        cwd=str(here),
        creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess, "CREATE_NO_WINDOW") else 0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
