#!/bin/bash
cd "$(dirname "$0")"
export TZ="${TZ:-$(readlink /etc/localtime | sed 's#.*/zoneinfo/##')}"
export PATH="$HOME/.docker/bin:/usr/local/bin:/opt/homebrew/bin:$PATH"
export SEQGUIDE_HOME="${SEQGUIDE_HOME:-$HOME/SeqGuide}"
docker compose down
echo "SeqGuide stopped."
