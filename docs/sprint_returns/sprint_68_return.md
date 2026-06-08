# Sprint 68 Return — P17 dim4 Phase-2a Panel (PASS, AT-DEPTH)

**Date:** 2026-06-08
**Base HEAD (sprint start):** `5b211e9` (Sprint 67 follow-up)
**Sprint goal:** Build + run the Phase-2a v1.2 panel for P17 (collective gradient sensing). Target overall TNR ≥ 0.95, Cohen's d ≥ 1.0.
**Tag:** `v0.68.0`
**Sprint type:** Chat-led design + code-led execution (Milestone B, Wave 1).

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `5b211e9` ✓
2. **Working tree clean:** ✓

---

## Part A — Wire P17 into the panel

### Invariance flags (`epc/phase2a/detector_invariance.py`)
- `permutation_invariant=True`: Group CI is computed from the center-of-mass position; permuting agent indices does not change the CoM → permutation_shuffled is degenerate-by-construction.
- `time_shuffle_invariant=False`: Gradient climbing is a temporal trajectory; shuffling temporal order destroys the directional approach signal.
- Citation: Berdahl et al. (2013) — CI measures net displacement toward field peak over time.

### Class A (synthetic nulls)
9 evaluated (1 skipped: permutation_shuffled). All use "particles" format. Rejected at prerequisite 1: no `field_samples` in history.

### Class B (catalog-mate discrimination)
Substrate type: `continuous_2d`. Catalog mates: P2_abp, P5_vicsek, P6_dorsogna, P7_lane_formation.
- **P5 Vicsek flocking:** Coordinated motion WITHOUT external-field inference → no `field_samples` in history → rejected at prerequisite 1. ✓
- **P6 D'Orsogna milling:** Same — particle trajectories without field → rejected. ✓
- **P2 ABP MIPS:** Same. ✓
- **P7 Lane formation:** Same. ✓

### Class C (failed regimes — `epc/phase2a/failed_regimes/p17_collective_sensing.py`)
10 regimes in two failure modes:

1. **social_off** (5 regimes, social_strength=0.0, noise ∈ {0.5, 0.8, 1.0, 1.2, 1.5}):
   Without social coupling, agents disperse uniformly → mean distance-to-CoM / L > 0.20 → rejected at prerequisite 3 (social cohesion). ✓

2. **field_too_strong** (5 regimes, amplitude ∈ {5, 8, 10, 10, 15}, noise ∈ {0.1, 0.15, 0.1, 0.05, 0.05}):
   Individual SNR >> 3.0 → single agents can reliably climb gradient → no group advantage → rejected at prerequisite 2. ✓

### Prerequisites (Sprint 30 rule: literature-grounded content prerequisites only)

1. **Field presence** (`field_samples` in history): P17 requires an external scalar field. Berdahl 2013: "the mechanism requires environmental information."
2. **Group advantage** (individual SNR ≤ 3.0): Derived from field_samples temporal statistics. If individual signal-to-noise is high, no collective amplification needed.
3. **Social cohesion** (mean dist-to-CoM / L ≤ 0.20): Berdahl 2013: social interactions localize the group for effective noise averaging.

---

## Part B — Run

```
PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p17
```

### Panel results

| Metric | Value |
|--------|-------|
| Overall TNR | **1.000** |
| Class A (synthetic) TNR | 1.000 (9/9 TN) |
| Class B (catalog) TNR | 1.000 (4/4 TN) |
| Class C (failed regime) TNR | 1.000 (10/10 TN) |
| Cohen's d | **11.117** |
| Verdict | **PASS** |

Positive scores: 0.900, 0.500, 0.700, 0.900, 0.900 (mean=0.780; 3 definitive, 1 confirmation, 1 screening).

Output: `analysis/outputs/p17_phase2a_panel.json` (gitignored).

---

## Part C — depth_gap + REPLICATION_NOTES + paper

- **depth_gap.md:** P17 dim4 pending→PASS, grade GAP→AT-DEPTH. AT-DEPTH count 19→20, gap count 2→1.
- **REPLICATION_NOTES.md:** Sprint 68 P17 Phase-2a section added with full class breakdown.
- **paper_CHANGELOG.md:** Sprint 68 entry added.
- Remaining gap: P12 (dim1 — closed-as-documented-limitation per Sprint 63).

---

## Post-flight checks

- `pytest tests/test_orchestration.py tests/test_transfer_matrix_counts.py tests/test_cross_detection_matrix.py -m "not slow"`: 109 passed ✓
- `pytest tests/test_collective_sensing_p17_e2e.py -m "not slow"`: 10 passed ✓
- P17 panel verdict: PASS (TNR=1.000, d=11.117) ✓
- No threshold relaxation (Sprint 30 rule) ✓

---

## Carry-forwards

- **C-p17-positive-seed-variance:** 1/5 canonical positives only reaches SCREENING (CI sensitive to initial offset direction at seed=1). Not a panel issue (TNR is the target metric) but notes variance in positive detection strength.

---

## Files added/modified

**New (1):**
- `epc/phase2a/failed_regimes/p17_collective_sensing.py` — Class C failed regimes (5 social_off + 5 field_too_strong)

**Modified (6):**
- `epc/phase2a/detector_invariance.py` — P17 permutation_invariant False→True
- `epc/phase2a/catalog.py` — +1 substrate (P17_collective_sensing) + generator + PATTERN_TO_SUBSTRATE_ID
- `analysis/run_phase2a_panel.py` — P17 panel wiring (build_p17_positives, make_p17_detector_fn, run_p17, dispatch)
- `docs/depth_gap.md` — P17 AT-DEPTH, counts updated
- `REPLICATION_NOTES.md` — Sprint 68 section
- `docs/paper_CHANGELOG.md` — Sprint 68 entry

---

**Decision: GO**
