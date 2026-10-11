#!/bin/bash
# EPC Emergence Census v1 — the registered run (OSF e5hkf; plan docs/census/PREREG_census_v1.md; code 52cee0f).
# Order: G1 census <= 22 bits -> G2 census <= 21 bits -> novelty arm. Runner stages are resumable: re-running this
# script skips finished programs. Run inside tmux as matthewhmaxwell:  tmux new -s census_run 'bash census/run_census.sh'
cd /home/matthewhmaxwell/epc-census
export PYTHONDONTWRITEBYTECODE=1
PY="nice -n 10 /home/matthewhmaxwell/epc-venv/bin/python -u"
R=census/runs; mkdir -p $R; L=$R/census_v1.log
echo "START $(date '+%F %T') commit $(git rev-parse --short HEAD)" >> $L
echo "G1 <= 22 start $(date '+%F %T')" >> $L
$PY -m census.runner --out $R/census_v1_g1 --max-bits 22 --workers 5 --grammar g1 >> $R/census_v1_g1.log 2>&1
echo "G1 done $(date '+%F %T')" >> $L
echo "G2 <= 21 start $(date '+%F %T')" >> $L
$PY -m census.runner --out $R/census_v1_g2 --max-bits 21 --workers 5 --grammar g2 >> $R/census_v1_g2.log 2>&1
echo "G2 done $(date '+%F %T')" >> $L
if [ ! -f $R/census_v1_novelty/elites.json ] || ! grep -q "NOVELTY DONE" $L; then
  echo "novelty arm start $(date '+%F %T')" >> $L
  $PY -m census.novelty --out $R/census_v1_novelty --evals 150000 --patience 20000 --workers 5 >> $R/census_v1_novelty.log 2>&1
  echo "NOVELTY DONE $(date '+%F %T')" >> $L
fi
echo "ALL DONE $(date '+%F %T')" >> $L
