@echo off
rem Windows Application Control blocks .venv\Scripts\python.exe here, so use the uv base
rem interpreter the venv was built from, with the venv's site-packages on PYTHONPATH.
setlocal
set "APPDIR=%~dp0"
set "PY=%APPDATA%\uv\python\cpython-3.10-windows-x86_64-none\python.exe"
if not exist "%PY%" set "PY=%APPDIR%.venv\Scripts\python.exe"
set "PYTHONPATH=%APPDIR%.venv\Lib\site-packages;%APPDIR%src;%APPDIR%."
set "PYTHONUNBUFFERED=1"
cd /d "%APPDIR%"
"%PY%" %*
exit /b %ERRORLEVEL%
