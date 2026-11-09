@echo off
REM End-to-end validation on Windows: builds SeqGuide, checks all tools, runs every sample in samples.tsv,
REM copies the results next to this script, then starts the GUI at http://localhost:8501
setlocal
cd /d "%~dp0"
set "HERE=%CD%"
cd ..
set "APP=%CD%"
set "SEQGUIDE_HOME=%USERPROFILE%\SeqGuide_validation"
set "LOG=%HERE%\validation.log"

docker info >nul 2>&1
if errorlevel 1 (
  echo Docker Desktop is not running. Start it, wait until it says "Engine running", then try again.
  pause
  exit /b 1
)
if not exist "%SEQGUIDE_HOME%\input" mkdir "%SEQGUIDE_HOME%\input"
if not exist "%SEQGUIDE_HOME%\runs" mkdir "%SEQGUIDE_HOME%\runs"
xcopy /y /q "%HERE%\input\*" "%SEQGUIDE_HOME%\input\" >nul
copy /y "%HERE%\samples.tsv" "%SEQGUIDE_HOME%\samples.tsv" >nul
if exist "%SEQGUIDE_HOME%\runs\validation" rmdir /s /q "%SEQGUIDE_HOME%\runs\validation"
if exist "%LOG%" del "%LOG%"
cd /d "%APP%"

call :log "=== SeqGuide validation %DATE% %TIME% ==="
call :run cmd /c ver
call :run docker version
call :log "--- build ---"
call :run docker compose build
call :log "--- tool check ---"
call :run docker compose run --rm seqguide python -m seqguide check
call :log "--- batch ---"
call :run docker compose run --rm seqguide python -m seqguide batch /data/samples.tsv --out /data/runs/validation
call :log "--- copy results back ---"
if exist "%HERE%\results" rmdir /s /q "%HERE%\results"
xcopy /e /i /y /q "%SEQGUIDE_HOME%\runs\validation" "%HERE%\results" >nul
call :log "--- start GUI ---"
call :run docker compose up -d
call :log "=== finished %DATE% %TIME% ==="
call :log "VALIDATION SCRIPT DONE"
start "" http://localhost:8501
echo.
echo Done. Please send validation.log and the results folder (both next to this script).
pause
exit /b 0

:log
echo %~1
>>"%LOG%" echo %~1
exit /b 0

:run
powershell -NoProfile -Command "& { $ErrorActionPreference='Continue'; & %* 2>&1 | ForEach-Object { \"$_\" } | Tee-Object -FilePath $env:LOG -Append }"
exit /b 0
