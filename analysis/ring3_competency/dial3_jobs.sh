#!/bin/bash
# dial3 sweep: 6 alphas x 2 conditions x 10 seeds, 6 in parallel, resumable (finished runs are skipped).
cd /home/matthewhmaxwell/emergent-pattern-catalog/analysis/ring3_competency
export DIAL_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
for s in $(seq 0 9); do for a in 0.0 0.1 0.25 0.5 0.75 1.0; do for c in channel nochannel; do echo "$a $s $c"; done; done; done |
  xargs -P 6 -L 1 bash -c "nice -n 10 /home/matthewhmaxwell/epc-venv/bin/python dial3_ppo.py run --alpha \$0 --seed \$1 --cond \$2 >> dial3_runs.log 2>&1"
/home/matthewhmaxwell/epc-venv/bin/python dial3_ppo.py summarize > dial3_summary.txt 2>&1
echo DIAL3_DONE >> dial3_runs.log
