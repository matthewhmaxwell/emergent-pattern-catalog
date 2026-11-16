#!/bin/bash
# Repro v2 group A, queue 1 of 2: P23 minority game (full protocol). One process at a time.
set -euo pipefail
cd "$(dirname "$0")/../../.."
PY=/home/matthewhmaxwell/epc-venv/bin/python
S=analysis/reproductions/v2/p23_minority_game.py
export PYTHONPATH=.
echo "queue1 start $(date -Is) commit $(git rev-parse HEAD)"
nice -n 10 $PY $S t1 --tie random --seeds 0-19
nice -n 10 $PY $S stationarity --tie random --seeds 0-19 --rerun
nice -n 10 $PY $S t1 --tie first --seeds 0-9
nice -n 10 $PY $S stationarity --tie first --seeds 0-9 --rerun
nice -n 10 $PY $S t2 --payoff linear --seeds 0-9
nice -n 10 $PY $S t2 --payoff sign --seeds 0-9
echo "queue1 done $(date -Is)"
