#!/usr/bin/env bash
# P3 MISS -> novelty protocol (osf.io/up7m9, PREREG §4). Step 1: replicate the MISS condition (rho=1, full + fair
# memoryless baseline) with seeds 4-6 on a byte-identical copy of the frozen code. Step 2: cue-landmark conflict test +
# resample ablation (p3_conflict.py) over every P3 model. Niced, resumable. Run inside tmux as matthewhmaxwell.
set -u
cd "$(dirname "$0")"
export PY=/home/matthewhmaxwell/epc-venv/bin/python PYTHONDONTWRITEBYTECODE=1 PREREG_THREADS=2
mkdir -p replication_P3/prereg_runs/logs
[ -f replication_P3/prereg_envs.py ] || cp prereg_envs.py replication_P3/prereg_envs.py
cmp -s prereg_envs.py replication_P3/prereg_envs.py || { echo "FROZEN CODE MISMATCH - abort"; exit 1; }

run_one() {
  v=$1; s=$2; tag="P3_${v}_rho1.0_s$s"
  [ -f "replication_P3/prereg_runs/$tag.json" ] && { echo "skip $tag"; return 0; }
  if (cd replication_P3 && nice -n 10 "$PY" prereg_envs.py train P3 "$v" "$s" --rho 1.0 > "prereg_runs/logs/$tag.log" 2>&1)
  then echo "$(date +%H:%M) done $tag"; else echo "$(date +%H:%M) FAIL $tag"; fi
}
export -f run_one

echo "$(date) START P3 replication (6 runs, PAR=3)"
printf "full 4\nfull 5\nfull 6\nmemoryless 4\nmemoryless 5\nmemoryless 6\n" | xargs -P 3 -L 1 bash -c 'run_one "$0" "$1"'
echo "$(date) START conflict test"
nice -n 10 "$PY" p3_conflict.py 2>&1 | tee P3_conflict.log
echo "$(date) P3 PROTOCOL DONE"
