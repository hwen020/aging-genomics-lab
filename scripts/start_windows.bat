@echo off
setlocal
cd /d "%~dp0.."
if not exist .venv\Scripts\python.exe (
  echo [1/3] Creating Python virtual environment...
  py -3 -c "import sys; sys.exit(sys.version_info < (3,10))" >nul 2>nul
  if errorlevel 1 (
    python -c "import sys; sys.exit(sys.version_info < (3,10))" >nul 2>nul
    if errorlevel 1 goto :fail
    python -m venv .venv
  ) else (
    py -3 -m venv .venv
  )
  if errorlevel 1 goto :fail
)
.venv\Scripts\python.exe -c "import aginggenomics, streamlit, plotly, scipy, pytest" >nul 2>nul
if errorlevel 1 (
  echo [2/3] Installing package and dependencies...
  .venv\Scripts\python.exe -m pip install -e ".[dev]"
  if errorlevel 1 goto :fail
)
echo [3/3] Opening Aging Genomics Lab at http://localhost:8501
echo Press Ctrl+C to stop.
.venv\Scripts\python.exe -m streamlit run app.py --server.address 127.0.0.1 --server.port 8501
goto :eof
:fail
echo Setup failed. Verify Python 3.10+ and internet connection for initial installation.
pause
exit /b 1
