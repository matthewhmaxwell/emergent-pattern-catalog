#!/bin/bash
# Post-run registered analyses (prereg S3, S4.5, S5, S10). Unattended; survives the owner's laptop disconnecting.
#   tmux new -d -s postrun 'bash census/postrun.sh'
cd /home/matthewhmaxwell/epc-census
set -a; . ~/.config/epc-litcheck.env; set +a
export PYTHONPATH=.
PY="nice -n 10 /home/matthewhmaxwell/epc-venv/bin/python -u"
O=census/runs/postrun; mkdir -p $O
echo "START $(date '+%F %T') commit $(git rev-parse --short HEAD)" > $O/postrun.log
$PY -m census.postrun ffaudit > $O/ffaudit.log 2>&1 &
FF=$!
for s in label minimal robust review; do
  echo "--- $s $(date '+%F %T')" >> $O/postrun.log
  $PY -m census.postrun $s > $O/$s.log 2>&1 || echo "FAIL $s" >> $O/postrun.log
done
echo "--- waiting for ffaudit $(date '+%F %T')" >> $O/postrun.log
wait $FF || echo "FAIL ffaudit" >> $O/postrun.log
$PY -m census.postrun summary > $O/summary.log 2>&1 || echo "FAIL summary" >> $O/postrun.log
echo "ALL DONE $(date '+%F %T')" >> $O/postrun.log
