#!/usr/bin/env bash
# EPC Ring-3 pre-registration v1 (osf.io/up7m9): launch all 42 registered training runs (niced, resumable:
# finished runs are skipped), then P4 cross-play, then the registered verdicts. Run inside tmux as matthewhmaxwell.
set -u
cd "$(dirname "$0")"
export PY=/home/matthewhmaxwell/epc-venv/bin/python PYTHONDONTWRITEBYTECODE=1
export PREREG_THREADS=${PREREG_THREADS:-2}
PAR=${PAR:-3}
mkdir -p prereg_runs/logs

joblist() {
  for s in 1 2 3; do
    echo "P1 full $s"; echo "P1 nochannel $s"
    echo "P4 full $s"; echo "P4 nochannel $s"
    echo "P2 full $s"; echo "P2 nochannel $s"; echo "P2 noobs $s"
    for r in 1.0 0.75 0.5 0.25 0.0; do echo "P3 full $s --rho $r"; done
    for r in 1.0 0.0; do echo "P3 memoryless $s --rho $r"; done
  done
}

run_one() {
  task=$1; var=$2; seed=$3; rho=""; [ "${4:-}" = "--rho" ] && rho=$5
  tag="${task}_${var}${rho:+_rho$rho}_s${seed}"
  [ -f "prereg_runs/$tag.json" ] && { echo "skip $tag"; return 0; }
  if nice -n 10 "$PY" prereg_envs.py train "$task" "$var" "$seed" ${rho:+--rho $rho} > "prereg_runs/logs/$tag.log" 2>&1
  then echo "$(date +%H:%M) done $tag"; else echo "$(date +%H:%M) FAIL $tag (see prereg_runs/logs/$tag.log)"; fi
}
export -f run_one

echo "$(date) START: $(joblist | wc -l) jobs, PAR=$PAR, threads/job=$PREREG_THREADS"
joblist | xargs -P "$PAR" -I{} bash -c 'run_one {}'
"$PY" prereg_envs.py crossplay > prereg_runs/logs/crossplay.log 2>&1 && echo "crossplay done"
"$PY" prereg_analysis.py | tee prereg_runs/logs/verdicts.log
echo "$(date) ALLDONE"
