#!/bin/bash
# Daily census check (cron, matthewhmaxwell): triage + digest of every run stage, plus health -> census/runs/DAILY.md
cd /home/matthewhmaxwell/epc-census
export PYTHONDONTWRITEBYTECODE=1
PY="nice -n 15 /home/matthewhmaxwell/epc-venv/bin/python"
R=census/runs; LIB=census/reflib/census_v3/library.json; O=$R/DAILY.md
{
echo "# Census daily check — $(date '+%F %H:%M')"; echo
echo "## Health"
echo "- run session alive: $(tmux has-session -t census_run 2>/dev/null && echo yes || echo NO)"
echo "- worker processes: $(pgrep -u matthewhmaxwell -f 'census.runner|census.novelty' | wc -l)"
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
