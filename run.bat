@echo off
REM ============================================================================
REM  RepurposeIQ / NeuroRepurpose - launcher
REM
REM  Runs the Flask app with the Python 3.12 virtual environment (.venv312).
REM  This MUST be used instead of the system Python 3.14 install: Windows
REM  Smart App Control blocks RDKit's DLLs on 3.14, which breaks docking and
REM  3D structure generation. The 3.12 venv's RDKit wheel is SAC-trusted.
REM ============================================================================

cd /d "%~dp0"

set "VENV_PY=%~dp0.venv312\Scripts\python.exe"

if not exist "%VENV_PY%" (
    echo [ERROR] Virtual environment not found at .venv312
    echo         Recreate it with:  py -3.12 -m venv .venv312
    echo         then install deps into it. See .interface-design/ and project notes.
    pause
    exit /b 1
)

echo Starting RepurposeIQ on http://127.0.0.1:5000  (Ctrl+C to stop)
"%VENV_PY%" flask_app.py
