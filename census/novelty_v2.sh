#!/bin/bash
# DEVIATION 2: corrected quality-diversity arm (search, then naming + grouping + literature check). Unattended, resumable:
#   tmux new -d -s qd2 'bash census/novelty_v2.sh'     (re-run the same command after a crash/reboot to resume)
cd /home/matthewhmaxwell/epc-census
set -a; . ~/.config/epc-litcheck.env; set +a
export PYTHONPATH=.
PY="nice -n 10 /home/matthewhmaxwell/epc-venv/bin/python -u"
O=census/runs/census_v1_novelty_v2; mkdir -p $O
echo "START/RESUME $(date '+%F %T') commit $(git rev-parse --short HEAD)" >> $O/qd2.log
[ -f $O/SEARCH_DONE.json ] || $PY -m census.novelty_v2 search --out $O --workers 5 >> $O/search.log 2>&1 || echo "FAIL search" >> $O/qd2.log
[ -f $O/SEARCH_DONE.json ] && { $PY -m census.novelty_v2 downstream --out $O >> $O/downstream.log 2>&1 || echo "FAIL downstream" >> $O/qd2.log; }
echo "ALL DONE $(date '+%F %T')" >> $O/qd2.log
