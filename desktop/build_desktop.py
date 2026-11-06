"""Assemble the portable Windows desktop build into ``dist/Root3DMeasurer``.

Usage (any Python >= 3.11 with internet access)::

    python desktop/build_desktop.py [--zip]

The result is a self-contained folder: an embedded official CPython runtime
(``runtime\\``) with the project wheel and pywebview installed via pip, plus a
small ``Root3DMeasurer.exe`` launcher. Freezing the measurement pipeline is
deliberately avoided: measurement jobs spawn ``sys.executable -B <script>``
subprocesses, which keep working unchanged because the embedded runtime ships
a real python.exe/pythonw.exe pair.
"""
import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build" / "desktop"
DIST = ROOT / "dist" / "Root3DMeasurer"
BUILD_ENV = ROOT / ".build-desktop"

PYTHON_EMBED_VERSION = os.environ.get("ROOT_MEASURER_PYTHON_EMBED", "3.12.10")
EMBED_URL = "https://www.python.org/ftp/python/{version}/python-{version}-embed-amd64.zip"
GET_PIP_URL = "https://bootstrap.pypa.io/get-pip.py"
PIP_INDEX = os.environ.get("ROOT_MEASURER_PIP_INDEX", "https://mirrors.aliyun.com/pypi/simple/")

BUILD_ENV_PACKAGES = ["build>=1,<2", "pyinstaller>=6,<7", "pillow>=10"]


def pip_command(python_exe, *arguments, **kwargs):
    return run([python_exe, "-m", "pip", "install", "--quiet", "-i", PIP_INDEX, *arguments], **kwargs)

LAUNCHER_NAME = "Root3DMeasurer.exe"
RUNTIME = DIST / "runtime"


def run(command, **kwargs):
    print("+", " ".join(str(part) for part in command), flush=True)
    completed = subprocess.run([str(part) for part in command], **kwargs)
    if completed.returncode:
        raise SystemExit(f"command failed ({completed.returncode}): {command}")
    return completed


def download(url, target):
    if target.is_file():
        print("cached:", target)
        return
    print("downloading:", url, flush=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + ".part")
    with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as stream:
        shutil.copyfileobj(response, stream)
    partial.rename(target)


def ensure_build_env():
    if not (BUILD_ENV / "Scripts" / "python.exe").is_file():
        run([sys.executable, "-m", "venv", BUILD_ENV])
    packages = BUILD_ENV / "Scripts" / "python.exe"
    run([packages, "-m", "pip", "install", "--quiet", "--upgrade", "pip", "-i", PIP_INDEX])
    pip_command(packages, *BUILD_ENV_PACKAGES)
    return packages


def build_wheel(build_python):
    wheels = BUILD / "wheels"
    run([build_python, "-m", "build", "--wheel", "--outdir", wheels, ROOT])
    candidates = sorted(wheels.glob("root_3d_measurer-*.whl"))
    if not candidates:
        raise SystemExit("project wheel was not produced")
    return candidates[-1]


def prepare_runtime(build_python):
    DIST.mkdir(parents=True, exist_ok=True)
    if RUNTIME.is_dir():
        print("removing previous runtime")
        shutil.rmtree(RUNTIME)
    archive = BUILD / f"python-{PYTHON_EMBED_VERSION}-embed-amd64.zip"
    download(EMBED_URL.format(version=PYTHON_EMBED_VERSION), archive)
    with zipfile.ZipFile(archive) as bundle:
        bundle.extractall(RUNTIME)
    tag = PYTHON_EMBED_VERSION.rpartition(".")[0].replace(".", "")
    pth = RUNTIME / f"python{tag}._pth"
    lines = pth.read_text(encoding="utf-8").splitlines()
    if "import site" not in lines:
        lines = [line for line in lines if line.strip() != "# import site"]
        lines += ["Lib/site-packages", "import site"]
        pth.write_text("\n".join(lines) + "\n", encoding="utf-8")
    get_pip = BUILD / "get-pip.py"
    download(GET_PIP_URL, get_pip)
    python_exe = RUNTIME / "python.exe"
    run([python_exe, get_pip, "--no-warn-script-location", "-i", PIP_INDEX], cwd=RUNTIME)
    run([python_exe, "-m", "pip", "install", "--quiet", "--upgrade", "pip", "-i", PIP_INDEX], cwd=RUNTIME)
    # The embeddable runtime cannot create pip build-isolation environments;
    # ship setuptools/wheel and build sdist-only dependencies in place.
    pip_command(python_exe, "setuptools", "wheel", cwd=RUNTIME)
    wheel = build_wheel(build_python)
    pip_command(python_exe, "--no-build-isolation", "--prefer-binary",
                str(wheel), "pywebview>=5", cwd=RUNTIME)
    run([python_exe, "-c", "import root_measurer,webview,numpy,scipy,skimage,shapely,psutil;print('runtime check OK')"], cwd=RUNTIME)


def make_icon():
    icon = Path(__file__).resolve().parent / "app.ico"
    if icon.is_file():
        return icon
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        print("pillow unavailable; launcher keeps the default icon")
        return None
    size = 256
    image = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle([8, 8, size - 8, size - 8], radius=52, fill=(15, 118, 110, 255))
    ink = (240, 253, 250, 255)
    draw.line([128, 52, 128, 196], fill=ink, width=14)          # taproot
    draw.line([128, 100, 76, 148], fill=ink, width=11)           # left branch
    draw.line([128, 128, 180, 168], fill=ink, width=11)          # right branch
    draw.line([76, 148, 60, 196], fill=ink, width=8)
    draw.line([180, 168, 196, 208], fill=ink, width=8)
    draw.line([128, 52, 108, 30], fill=ink, width=8)             # stem tips
    draw.line([128, 52, 152, 30], fill=ink, width=8)
    draw.ellipse([116, 40, 140, 64], fill=(240, 253, 250, 255))
    image.save(icon, sizes=[(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    return icon


def build_launcher(build_python):
    icon = make_icon()
    command = [build_python, "-m", "PyInstaller", "--onefile", "--noconsole", "--clean", "--noconfirm",
               "--distpath", DIST, "--workpath", BUILD / "pyinstaller", "--specpath", BUILD,
               "--name", "Root3DMeasurerLauncher"]
    if icon:
        command += ["--icon", str(icon)]
    command.append(str(Path(__file__).resolve().parent / "launcher.py"))
    run(command)
    target = DIST / LAUNCHER_NAME
    if target.exists():
        target.unlink()
    (DIST / "Root3DMeasurerLauncher.exe").rename(target)


def write_readme():
    source = Path(__file__).resolve().parent / "使用说明.txt"
    if source.is_file():
        shutil.copyfile(source, DIST / "使用说明.txt")


def maybe_zip(enable):
    if not enable:
        return
    print("packing zip …", flush=True)
    archive = shutil.make_archive(str(ROOT / "dist" / f"Root3DMeasurer-win64-x64"), "zip", root_dir=str(ROOT / "dist"), base_dir="Root3DMeasurer")
    print("zip:", archive)


def main():
    parser = argparse.ArgumentParser(description="构建免安装 Windows 桌面版")
    parser.add_argument("--zip", action="store_true", help="额外打包 dist\\Root3DMeasurer-win64-x64.zip")
    args = parser.parse_args()
    build_python = ensure_build_env()
    prepare_runtime(build_python)
    build_launcher(build_python)
    write_readme()
    print("\nDONE:", DIST)
    maybe_zip(args.zip)


if __name__ == "__main__":
    main()
