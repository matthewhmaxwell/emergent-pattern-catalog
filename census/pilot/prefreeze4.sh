#!/bin/bash
# Pre-freeze re-check after the phase-view fix + continuous-evidence fix (2026-10-02). Everything on fresh seeds.
cd /home/matthewhmaxwell/epc-census
export PYTHONDONTWRITEBYTECODE=1
PY="nice -n 10 /home/matthewhmaxwell/epc-venv/bin/python -u"
L=census/pilot/prefreeze4.log
echo "START $(date '+%F %T')" > $L
$PY -m census.t_phase_noise 10 2>&1 | grep -v -i warn >> $L
echo "--- library $(date '+%T')" >> $L
$PY -m census.reflib build --out census/reflib/census_v3 --workers 6 --seedsets 2 --seedset-start 30 > census/reflib_census_v3.log 2>&1
tail -3 census/reflib_census_v3.log >> $L
for s in a b; do
  echo "--- placebo $s $(date '+%T')" >> $L
  $PY -m census.placebo --from-run census/pilot/speed22 census/pilot/rehearsal_g1_14 census/pilot/rehearsal_g2_14 --out census/validation/placebo2_$s --workers 6 --salt fix2_$s > census/validation_placebo2_$s.log 2>&1
  sed -n 3,14p census/validation/placebo2_$s/PLACEBO.md >> $L
done
i=0
for off in 9000011 10000019 11000027 12000017; do
  i=$((i+1)); echo "--- gate $i (seed offset $off) $(date '+%T')" >> $L
  $PY -m census.validate4 --lib census/reflib/census_v3/library.json --workers 6 --seed-offset $off --free-start 2800 --extended --out census/validation/gate2_$i > census/validation_gate2_$i.log 2>&1
  sed -n '/^## Criteria/,/reported) P4/p' census/validation/gate2_$i/VALIDATION.md | cut -c1-200 >> $L
  grep -A6 "Negatives flagged" census/validation/gate2_$i/VALIDATION.md >> $L
  grep "GATE:" census/validation/gate2_$i/VALIDATION.md >> $L
done
for g in g1 g2; do
  echo "--- rehearsal $g $(date '+%T')" >> $L
  rm -rf census/pilot/rehearsal4_${g}_14
  $PY -m census.runner --out census/pilot/rehearsal4_${g}_14 --max-bits 14 --workers 6 --grammar $g > census/pilot/rehearsal4_${g}_14.log 2>&1
  $PY -m census.triage census/pilot/rehearsal4_${g}_14 --lib census/reflib/census_v3/library.json > census/pilot/rehearsal4_${g}_14_triage.log 2>&1
  $PY -m census.digest census/pilot/rehearsal4_${g}_14 > census/pilot/rehearsal4_${g}_14_digest.log 2>&1
  sed -n 3,4p census/pilot/rehearsal4_${g}_14/TRIAGE.md >> $L
done
echo "ALL DONE $(date '+%F %T')" >> $L
