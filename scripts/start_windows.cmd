@echo off
setlocal
cd /d "%~dp0.."
if exist ".venv\Scripts\python.exe" (
  ".venv\Scripts\python.exe" -m root_measurer gui --workspace runs\default --open-browser
) else (
  python -m root_measurer gui --workspace runs\default --open-browser
)
if errorlevel 1 (
  echo Install this tool first: python -m pip install -e .
  pause
)
