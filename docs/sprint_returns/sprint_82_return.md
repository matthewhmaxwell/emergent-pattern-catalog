# Sprint 82 Return — P25 dim4 Phase-2a Panel

**Date:** 2026-06-09
**Base HEAD (sprint start):** `d1e7af0` (Sprint 81 follow-up)
**Sprint goal:** Build + run Phase-2a v1.2 panel for P25 (equifinality). Target overall TNR ≥ 0.95, Cohen's d ≥ 1.0.
**Tag:** `v0.82.0`
**Sprint type:** Chat-led design + code-led execution. Milestone B, Wave 3.

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `d1e7af0` ✓
2. **Working tree clean:** ✓

---

## Part A — Wire P25 into the Phase-2a panel

**Detector format:** `canalization_bundle` — multi-IC observation bundles with
keys: state, step, trial, ic, target, distance_to_target, converged.

**Invariance flags (P25):**
- `permutation_invariant=True`: convergence variance ratio is an aggregate
  statistic over the IC ensemble; permuting trial indices does not change
  per-dimension variances.
- `time_shuffle_invariant=True`: observation bundle sorts records by (trial,
  step); ICs/finals are first/last records per trial. Shuffling list order
  preserves step labels → sorted extraction recovers original ICs/finals.

**Class A synthetics:** `canalization_bundle` format added to all 10 generators.
Non-converging trajectories (random walks, constant, divergent) — P25 rejects
at screening (ratio ≥ 0.1 or basin volume < 0.5).

**Class B:**
- Catalog mates: 0 (canalization_landscape is unique substrate type)
- Supplements: 2 — `diffusive_multi_ic` (pure diffusion → ratio >> 0.1),
  `homeostatic_regulation_bundle` (P24-like proportional control → basin
  volume < 0.5 due to noise fluctuation around setpoint)

**Class C failed regimes:** 2
- `narrow_basin`: only near-target ICs converge (wide spread, diverge radius
  0.5 → most ICs diffuse, basin volume << 0.8)
- `divergent_dynamics`: repulsive dynamics push ICs away from target
  (ratio >> 0.1)

**Content prerequisite added:** basin_volume ≥ 0.5 at screening. Literature
grounding: Waddington (1957) — equifinality requires convergence from a wide
IC range, not just low variance ratio. This discriminates P25 from noisy
homeostatic regulation (P24) where the regulated variable fluctuates near the
setpoint without reliably converging within the convergence threshold.

---

## Part B — Panel results

```
P25    TNR=1.000  syn=1.000  cat=1.000  fai=1.000  d=+inf  PASS
```

- 5 canonical positives: all DEFINITIVE (confidence 0.90)
- Class A: 8/8 evaluated synthetics rejected, 2 skipped (degenerate)
- Class B: 2/2 supplements rejected
- Class C: 2/2 failed regimes rejected
- Overall TNR: 1.000 (≥ 0.95 ✓)
- Cohen's d: +∞ (≥ 1.0 ✓)
- Verdict: **PASS**

Output: `analysis/outputs/p25_phase2a_panel.json`

---

## Part C — Documentation updates

- **depth_gap.md:** P25 row updated (pending→PASS, GAP→AT-DEPTH). AT-DEPTH
  25→26, gap 2→1.
- **paper_section4_draft.md:** §4.25 dim4 paragraph added.
- **paper_section6_draft.md:** Sprint 82 entry. AT-DEPTH 26/27.
- **paper_CHANGELOG.md:** Sprint 82 entry.
- **REPLICATION_NOTES.md:** P25 dim4 replication section added.

---

## Post-flight checks

- `pytest tests/ -m "not slow"`: all passing ✓
- `p25_phase2a_panel.json`: verdict=PASS, TNR=1.000 ✓

---

## Open Carry-Forwards

- **C-p7-time-shuffled-fp:** `time_shuffled` FP at screening (Sprint 66, low priority)
- **C-p19-bias-zero-chance-alignment:** 1/5 bias_zero Class C FP at confirmation (Sprint 70, 4% rate)
- **P12 dim1:** documented finite-size measurement limitation (accepted)

---

## Files added/modified

**New (1):**
- `epc/phase2a/failed_regimes/p25_equifinality.py` — 2 Class C failed regimes (narrow_basin, divergent_dynamics)

**Modified (10):**
- `epc/detectors/p25_equifinality.py` — content prerequisite: basin_volume ≥ 0.5 at screening
- `epc/phase2a/detector_invariance.py` — P25 flags (perm=True, time=True)
- `epc/phase2a/synthetic.py` — canalization_bundle format added to all 10 generators
- `epc/phase2a/structured.py` — 2 canalization_landscape supplements
- `epc/phase2a/catalog.py` — P25 substrate generator + adapter + PATTERN_TO_SUBSTRATE_ID
- `epc/phase2a/panel.py` — canalization_bundle format in Class A kwargs
- `analysis/run_phase2a_panel.py` — P25 panel runner wired
- `docs/depth_gap.md` — P25 dim4 PASS, AT-DEPTH 26/27
- `docs/paper_section4_draft.md` — §4.25 dim4 paragraph
- `docs/paper_section6_draft.md` — Sprint 82 entry
- `docs/paper_CHANGELOG.md` — Sprint 82 entry
- `REPLICATION_NOTES.md` — P25 dim4 replication section
- `analysis/outputs/p25_phase2a_panel.json` — panel output

---

**Decision: GO** — Sprint completed cleanly. P25 dim4 Phase-2a panel PASS (TNR=1.000, d=+inf). All 4 dims PASS → P25 AT-DEPTH. AT-DEPTH count 26/27 (P12 sole remaining GAP). Content prerequisite (basin_volume ≥ 0.5 at screening) literature-grounded and discriminates from P24. Chain may proceed autonomously.
