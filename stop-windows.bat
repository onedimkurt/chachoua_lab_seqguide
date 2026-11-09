@echo off
cd /d "%~dp0"
set "SEQGUIDE_HOME=%USERPROFILE%\SeqGuide"
docker compose down
echo SeqGuide stopped.
pause
