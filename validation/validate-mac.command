#!/bin/bash
# End-to-end validation: builds SeqGuide, checks all tools, runs every sample in samples.tsv,
# then starts the GUI on the same data.
# macOS does not let Docker read Desktop/Documents/Downloads, so the data is copied to
# ~/SeqGuide_validation, analysed there, and the results are copied back next to this script.
HERE="$(cd "$(dirname "$0")" && pwd)"
APP="$(cd "$HERE/.." && pwd)"
export TZ="${TZ:-$(readlink /etc/localtime | sed 's#.*/zoneinfo/##')}"
export PATH="$HOME/.docker/bin:/usr/local/bin:/opt/homebrew/bin:$PATH"
export SEQGUIDE_HOME="$HOME/SeqGuide_validation"
LOG="$HERE/validation.log"
cd "$APP"
{
  echo "=== SeqGuide validation $(date) ==="
  uname -a
  docker version --format 'Docker client {{.Client.Version}}, server {{.Server.Version}} ({{.Server.Arch}})' || { echo "DOCKER NOT RUNNING"; exit 1; }
  echo "--- copy data to $SEQGUIDE_HOME ---"
  mkdir -p "$SEQGUIDE_HOME/input" "$SEQGUIDE_HOME/runs"
  cp -n "$HERE"/input/* "$SEQGUIDE_HOME/input/"
  cp "$HERE/samples.tsv" "$SEQGUIDE_HOME/samples.tsv"
  echo "--- build ---"
  time docker compose build
  echo "--- tool check ---"
  docker compose run --rm seqguide python -m seqguide check
  echo "--- batch ---"
  rm -rf "$SEQGUIDE_HOME/runs/validation"
  time docker compose run --rm seqguide python -m seqguide batch /data/samples.tsv --out /data/runs/validation
  echo "--- copy results back ---"
  rm -rf "$HERE/results"; cp -R "$SEQGUIDE_HOME/runs/validation" "$HERE/results"
  echo "--- start GUI ---"
  docker compose up -d
  echo "=== finished $(date) ==="
} 2>&1 | tee "$LOG"
echo "VALIDATION SCRIPT DONE" >> "$LOG"
