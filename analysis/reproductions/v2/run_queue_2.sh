#!/bin/bash
# Repro v2 group A, queue 2 of 2: P18 voter, then P27 Nowak-May (full protocols). One process at a time.
set -euo pipefail
cd "$(dirname "$0")/../../.."
PY=/home/matthewhmaxwell/epc-venv/bin/python
V=analysis/reproductions/v2/p18_voter.py
N=analysis/reproductions/v2/p27_nowak_may.py
export PYTHONPATH=.
echo "queue2 start $(date -Is) commit $(git rev-parse HEAD)"
nice -n 10 $PY $V theory
nice -n 10 $PY $V run --L 512 --nb von_neumann --seeds 0-19
nice -n 10 $PY $V precision
nice -n 10 $PY $V run --L 256 --nb von_neumann --seeds 0-9
nice -n 10 $PY $V run --L 256 --nb moore --seeds 0-9
nice -n 10 $PY $N invariance
nice -n 10 $PY $N run --L 1000 --b 1.9 --seeds 0-9
nice -n 10 $PY $N run --L 500 --b 1.9 --seeds 0-9
nice -n 10 $PY $N run --L 200 --b 1.9 --seeds 0-9
nice -n 10 $PY $N run --L 1000 --b 1.8 --seeds 0-9
echo "queue2 done $(date -Is)"
