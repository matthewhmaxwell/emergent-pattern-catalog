#!/bin/bash
# Daily census check (cron, matthewhmaxwell): triage + digest of every run stage, plus health -> census/runs/DAILY.md
cd /home/matthewhmaxwell/epc-census
export PYTHONDONTWRITEBYTECODE=1
PY="nice -n 15 /home/matthewhmaxwell/epc-venv/bin/python"
R=census/runs; LIB=census/reflib/census_v3/library.json; O=$R/DAILY.md
{
echo "# Census daily check — $(date '+%F %H:%M')"; echo
echo "## Health"
# Liveness = the actual run process + whether results advanced recently (tmux has-session is unreliable from cron/ssh:
# the socket often isn't reachable and it falsely reports dead even while the run is writing — 2026-10-14).
NW=$(pgrep -u matthewhmaxwell -f 'census\.(runner|novelty) --out census/runs/census_v1' | wc -l)
NEW=$(find $R/census_v1_*/results_*.jsonl -mmin -30 2>/dev/null | wc -l)
if [ "$NW" -gt 0 ] && [ "$NEW" -gt 0 ]; then ALIVE="yes (writing)"; elif [ "$NW" -gt 0 ]; then ALIVE="process up, NO new results in 30 min — CHECK"; else ALIVE="NO — run process not found, CHECK"; fi
echo "- run alive: $ALIVE"
echo "- run processes: $NW ; result files updated in last 30 min: $NEW"
echo "- disk free: $(df -h / | awk 'NR==2{print $4" ("$5" used)"}')"
echo "- last log lines: $(tail -2 $R/census_v1.log | tr '\n' ' ')"
for s in g1 g2 novelty; do
  d=$R/census_v1_$s; [ -d $d ] || continue
  if [ "$s" = novelty ]; then echo "- novelty: $(tail -1 $R/census_v1_novelty.log 2>/dev/null)"; fi
  $PY -m census.triage $d --lib $LIB > $d/triage.log 2>&1
  $PY -m census.digest $d > $d/digest.log 2>&1
  echo; echo "## Stage $s"; sed -n '3,4p' $d/TRIAGE.md 2>/dev/null; sed -n '3p' $d/digest.log 2>/dev/null
  echo; awk '/^## Behaviour map/{f=1;next} /^## Routes/{f=0} f' $d/TRIAGE.md 2>/dev/null | head -30
done
} > $O.tmp 2>&1 && mv $O.tmp $O
