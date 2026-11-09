@echo off
REM Double-click to start SeqGuide on Windows (needs Docker Desktop).
cd /d "%~dp0"
set "SEQGUIDE_HOME=%USERPROFILE%\SeqGuide"
if not exist "%SEQGUIDE_HOME%\input" mkdir "%SEQGUIDE_HOME%\input"
if not exist "%SEQGUIDE_HOME%\runs" mkdir "%SEQGUIDE_HOME%\runs"

docker info >nul 2>&1
if errorlevel 1 (
  echo Docker Desktop is not running. Start Docker Desktop, wait until it says "running", then try again.
  pause
  exit /b 1
)

echo Starting SeqGuide...
echo (The very first start downloads and installs all tools: 20-40 minutes. Later starts take seconds.)
docker compose up -d --build
if errorlevel 1 (
  echo Something went wrong, see the messages above.
  pause
  exit /b 1
)
timeout /t 10 /nobreak >nul
start "" http://localhost:8501
echo.
echo SeqGuide is running at http://localhost:8501
echo Your files: %SEQGUIDE_HOME%\input    Results: %SEQGUIDE_HOME%\runs
echo To stop it, double-click stop-windows.bat
pause
