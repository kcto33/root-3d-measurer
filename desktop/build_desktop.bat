@echo off
rem 一键构建免安装桌面版（需要已安装 Python 3.11+ 并可访问网络）。
setlocal
cd /d "%~dp0.."
python desktop\build_desktop.py %*
if errorlevel 1 (
  echo.
  echo 构建失败，请查看上方错误信息。
  pause
  exit /b 1
)
echo.
echo 构建完成：dist\Root3DMeasurer
pause
