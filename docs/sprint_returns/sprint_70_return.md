# Sprint 70 Return — P19 dim4 Phase-2a Panel (Completes Wave 1)

**Date:** 2026-06-08
**Base HEAD (sprint start):** `e4f2fa0` (Sprint 69 follow-up)
**Sprint goal:** Build + run Phase-2a v1.2 panel for P19 (emergent leadership). Target TNR ≥ 0.95, Cohen's d ≥ 1.0.
**Tag:** `v0.70.0`
**Sprint type:** Chat-led design + code-led execution (Milestone B, Wave 1 completion).

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `e4f2fa0` ✓
2. **Working tree clean:** ✓

---

## Part A — Wire P19 into Phase-2a panel

### Class A (synthetic)
continuous_2d particle substrates (10 synthetic nulls). P19 rejects all Class A
substrates because:
- Most lack `informed_mask` → prerequisite rejection
- `time_shuffled` and `permutation_shuffled` preserve `informed_mask` from the
  canonical positive, but the early-window leadership content prerequisite
  (Sprint 70) blocks false positives from scrambled temporal order

### Class B (catalog-mate discrimination)
5 continuous_2d catalog mates: P2_abp, P5_vicsek, P6_dorsogna, P7_lane_formation,
P17_collective_sensing. All lack `informed_mask` in history → prerequisite
rejection (detected=False). P5 (Vicsek) was the key risk per brief: directional
flocking with high φ but no informed minority → correctly rejected.

### Class C (failed regimes)
10 regimes in `epc/phase2a/failed_regimes/p19_informed_minority.py`:
- **5 rho_zero** (ρ=0, no informed agents): no leadership possible → all rejected
- **5 bias_zero** (ρ=0.1, ω=0, informed agents have no preference): 4/5 rejected;
  1 FP (seed=410) at confirmation due to chance alignment toward θ_pref

### Invariance flags
- `permutation_invariant=False`: influence asymmetry depends on informed/naive
  partition identity (Couzin 2005)
- `time_shuffle_invariant=False`: TE-based influence measurement requires
  temporal order

Both flags unchanged from Sprint 69 values in `epc/phase2a/detector_invariance.py`.

---

## Part B — Run + fix

### Initial run
Overall TNR = 0.880, verdict = PARTIAL. Three FPs:
1. `time_shuffled` at screening (0.500): each frame of the shuffled canonical
   positive still has high accuracy
2. `bias_zero_seed=410` at confirmation (0.700): chance alignment toward θ_pref
3. `bias_zero_seed=412` at screening (0.500): chance accuracy > 0.3

### Content prerequisite (Sprint 30 rule: literature-grounded, no threshold relaxation)
**Added:** Early-window informed→naive leadership gap (Couzin 2005).

In the convergence phase (10–40% of trajectory), informed agents' mean alignment
with θ_pref must exceed naive agents' alignment. This verifies that the informed
minority *leads* the alignment process — the defining P19 mechanism.

Grounding: Couzin et al. (2005) describe leadership as a process where the
informed minority has a "persistent directional preference" that steers the group.
Leaderless coordination (Vicsek alone) or time-shuffled trajectories show no
such early-window informed→naive gap.

### Final run

| Class | Size | TN | FP | TNR |
|-------|------|----|----|-----|
| A (synthetic) | 10 | 10 | 0 | 1.000 |
| B (catalog) | 5 | 5 | 0 | 1.000 |
| C (failed) | 10 | 9 | 1 | 0.900 |
| **Overall** | **25** | **24** | **1** | **0.960** |

**Cohen's d:** 5.418
**Verdict:** PASS (TNR ≥ 0.95, d ≥ 1.0)

Positives: 4/5 DEFINITIVE (0.900); seed 1 at SCREENING (0.500).

Output: `analysis/outputs/p19_phase2a_panel.json`.

---

## Part C — depth_gap + REPLICATION_NOTES + paper

- **depth_gap.md:** P19 dim4 pending→PASS, grade GAP→AT-DEPTH. AT-DEPTH count
  20→21. Gap count 2→1. Sprint 70 finding entry added.
- **REPLICATION_NOTES.md:** Sprint 70 P19 dim4 section added with full panel
  results, content prerequisite description, per-class TNR, carry-forward.
- **paper_section4_draft.md:** §4.22 (P19 emergent leadership) added — model,
  detector, architecture decision, dim1/dim2 results, Phase-2a panel results,
  distinctness from P5/P17/P18.
- **paper_section6_draft.md:** Sprints 67–70 entries added — P17 implementation
  + dim4, P19 implementation + dim4, Wave 1 completion.
- **paper_CHANGELOG.md:** Sprint 70 entry added.

---

## Post-flight checks

- `pytest tests/test_orchestration.py tests/test_transfer_matrix_counts.py tests/test_informed_minority_p19_e2e.py tests/test_cross_detection_matrix.py -m "not slow"`: **128 passed** ✓
- P19 panel verdict: PASS (TNR=0.960, d=5.418) ✓
- No threshold relaxation (Sprint 30 rule) ✓
- Content prerequisite is literature-grounded (Couzin 2005) ✓

---

## Carry-forwards

- **C-p19-bias-zero-chance-alignment:** 1/5 bias_zero regimes (seed=410)
  reaches confirmation by chance. The flock aligns toward θ_pref=0.0 by Vicsek
  dynamics alone, and the 20 "informed" agents happen to lead during early
  convergence. This is a 1/25 FP rate (4%), within tolerance for TNR ≥ 0.95.
  Potential fix: require minimum pull_mean effect size at confirmation (not just
  significance), but this would be threshold tightening → deferred.

- **C-p19-abrupt-saturation** (from Sprint 69): accuracy saturates to ~1.0 at
  ρ=0.025, more abruptly than Couzin (2005) Fig. 2a. Root cause: η=0.1 is low.
  Not a validation failure.

- **C-p19-te-vs-pull** (from Sprint 69): TE (KSG) on mean heading collapses in
  converged regime. Documented architectural decision; label-shuffle pull used.

---

## Files added/modified

**New (1):**
- `epc/phase2a/failed_regimes/p19_informed_minority.py` — Class C failed regimes (5 rho_zero + 5 bias_zero)

**Modified (7):**
- `epc/phase2a/catalog.py` — +1 substrate (P19_informed_minority) + generator + PATTERN_TO_SUBSTRATE_ID entry
- `analysis/run_phase2a_panel.py` — P19 panel wiring (build_positives, make_detector_fn with content prereq, run_p19, main dispatch)
- `analysis/outputs/p19_phase2a_panel.json` — Phase-2a panel output
- `docs/depth_gap.md` — P19 dim4 PASS, AT-DEPTH count 20→21, gap count 2→1
- `REPLICATION_NOTES.md` — Sprint 70 P19 dim4 section
- `docs/paper_section4_draft.md` — §4.22 P19 section added
- `docs/paper_section6_draft.md` — Sprints 67–70 entries added
- `docs/paper_CHANGELOG.md` — Sprint 70 entry

---

**Decision: GO**
