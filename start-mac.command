#!/bin/bash
# Double-click to start SeqGuide on a Mac (needs Docker Desktop).
cd "$(dirname "$0")"
export TZ="${TZ:-$(readlink /etc/localtime | sed 's#.*/zoneinfo/##')}"
export PATH="$HOME/.docker/bin:/usr/local/bin:/opt/homebrew/bin:$PATH"
export SEQGUIDE_HOME="${SEQGUIDE_HOME:-$HOME/SeqGuide}"
mkdir -p "$SEQGUIDE_HOME/input" "$SEQGUIDE_HOME/runs"

if ! docker info >/dev/null 2>&1; then
  echo "Starting Docker Desktop..."
  open -a Docker
  for _ in $(seq 1 90); do docker info >/dev/null 2>&1 && break; sleep 2; done
fi
if ! docker info >/dev/null 2>&1; then
  echo "Docker Desktop is not running. Open Docker from Applications, wait until it says 'Engine running', and try again."
  read -r -p "Press Enter to close."
  exit 1
fi

echo "Starting SeqGuide..."
echo "(The very first start downloads and installs all tools: 10-30 minutes. Later starts take seconds.)"
if ! docker compose up -d --build; then
  read -r -p "Something went wrong (see messages above). Press Enter to close."
  exit 1
fi
for _ in $(seq 1 60); do curl -s http://localhost:8501 >/dev/null && break; sleep 2; done
open "http://localhost:8501"
echo
echo "SeqGuide is running at http://localhost:8501"
echo "Your files: $SEQGUIDE_HOME/input    Results: $SEQGUIDE_HOME/runs"
echo "To stop it, double-click stop-mac.command."
