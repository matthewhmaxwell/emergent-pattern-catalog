#!/bin/bash
# Repro v2 group C scheduler (not a compute process). When P28 (slot B) has finished, stop the single-process
# P12 L = 100 primary run and continue all remaining P12 work (primary L = 100, then D2 L = 40, then D2 L = 100)
# with two queue workers, so at most two compute processes run at any time.
cd /home/matthewhmaxwell/epc-repro-c || exit 1
LOG=analysis/outputs/repro_v2/logs
PY=/home/matthewhmaxwell/epc-venv/bin/python
ENVS="PYTHONPATH=. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1"
echo "supervisor started $(date -u)" >> $LOG/supervisor.log
until grep -q P28_EXIT $LOG/slotB.log 2>/dev/null; do sleep 60; done
echo "P28 finished: $(grep P28_EXIT $LOG/slotB.log) at $(date -u)" >> $LOG/supervisor.log
pkill -f "^$PY analysis/repro_v2/p12_extinction.py run --variant primary --L 100"
sleep 10
rm -rf analysis/outputs/repro_v2/raw/p12/claims
for k in 1 2; do
  tmux new-session -d -s repro_c_q$k "cd $PWD && env $ENVS nice -n 10 $PY analysis/repro_v2/p12_extinction.py queue --jobs primary:100,D2:40,D2:100 >> $LOG/p12_queue_$k.log 2>&1; echo Q${k}_EXIT \$? >> $LOG/supervisor.log"
  sleep 5
done
echo "queue workers started $(date -u)" >> $LOG/supervisor.log
