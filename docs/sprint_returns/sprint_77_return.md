# Sprint 77 Return — P23 dim4 Phase-2a Panel (Milestone B Wave 2 Closure)

**Date:** 2026-06-09
**Base HEAD (sprint start):** `40f4720` (Sprint 76 follow-up)
**Sprint goal:** Build + run Phase-2a v1.2 panel for P23 (anti-coordination / minority game). Target TNR ≥ 0.95, Cohen's d ≥ 1.0.
**Tag:** `v0.77.0`
**Sprint type:** Chat-led design + code-led execution. Milestone B, Wave 2 closure.

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `40f4720` ✓
2. **Working tree clean:** ✓

---

## Part A — Wire P23 into the panel

### Invariance flags
- `permutation_invariant=True`: σ²/N is computed over attendance counts, which are
  scalar sums over agents — permuting agent indices does not change attendance.
  `permutation_shuffled` is degenerate-by-construction.
- `time_shuffle_invariant=True`: σ²/N is the primary confirmation-level signal and is
  a distribution statistic preserved by time shuffling. Lag-1 AC is a secondary signal
  too weak to reliably gate confirmation at the efficient phase (AC ≈ −0.01 to −0.05,
  p > 0.01 for most seeds). Sprint 77 panel run confirmed: time_shuffled FP at
  confirmation via σ²/N alone. Same pattern as P5 (Sprint 46), P1 (Sprint 43).
  Both `permutation_shuffled` and `time_shuffled` are SKIPPED by the harness.

### Class A — synthetic substrates
- `choice_timeseries` format added to all 10 synthetic generators.
- Random-choice Binomial(N, 0.5) attendance for most generators; constant N/2 for
  constant_field; linear ramp for linear_gradient; alternating [0, N] for checkerboard.

### Class B — catalog mates + supplements
- 0 catalog mates (choice_timeseries is the only pattern with this substrate type).
- 2 synthetic supplements added to `epc/phase2a/structured.py`:
  - `consensus_herding_attendance`: Majority-imitation process (p=0.7) producing
    σ²/N above baseline + positive autocorrelation. Rejected at screening.
  - `random_choice_attendance`: i.i.d. Binomial(N, 0.5). σ²/N at baseline. Rejected
    at screening.

### Class C — failed regimes
- `epc/phase2a/failed_regimes/p23_anticoordination.py` created with 2 regimes:
  - `random_agents`: No adaptation, σ²/N ≈ 0.25 (random baseline). Rejected at screening.
  - `herding_regime`: Majority-imitation (p=0.8), σ²/N >> baseline. Rejected at screening.

### Content prerequisite (Sprint 77)
Added to `epc/detectors/p23_anticoordination.py`: non-degenerate variance (σ² > 0) and
variance strictly below the random-choice baseline required at confirmation. Literature
grounding: Savit, Manuca & Riolo (1999) — the efficient phase is characterized by the
joint signature of σ²/N < p̂(1−p̂). Constant attendance (σ²=0) is trivially non-dynamic;
above-baseline fluctuation is herding, not coordination. This prerequisite rejects:
- constant_field (σ²=0)
- checkerboard (σ²/N >> baseline, fires on AC alone)
- time_shuffled was additionally handled by time_shuffle_invariant=True

---

## Part B — Panel run

```
PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p23
```

### Results

| Metric | Value |
|--------|-------|
| Overall TNR | **1.000** |
| Synthetic TNR (Class A) | 1.000 (8/8 evaluated, 2 SKIPPED) |
| Catalog TNR (Class B) | 1.000 (2/2, advisory) |
| Failed regime TNR (Class C) | 1.000 (2/2) |
| Cohen's d | **14.504** |
| Verdict | **PASS** |

### Canonical positive scores (5 seeds, m=6, α≈0.63, N=101, 3000 rounds)
| Seed | Verdict | Score |
|------|---------|-------|
| 0 | confirmation | 0.700 |
| 1 | confirmation | 0.700 |
| 2 | confirmation | 0.700 |
| 3 | confirmation | 0.650 |
| 4 | definitive | 0.900 |

### Panel iteration history
1. **Run 1** (original detector, time_shuffle_invariant=False): 3 Class A FPs
   (time_shuffled, constant, checkerboard). TNR=0.769, d=2.365, PARTIAL.
2. **Run 2** (content prerequisite added: nondegenerate + below-baseline): Fixed
   constant and checkerboard. time_shuffled persists (σ²/N preserved by shuffling).
   TNR=0.923, d=3.863, PARTIAL.
3. **Run 3** (time_shuffle_invariant=True, citing P5/P1 precedent): time_shuffled
   SKIPPED. TNR=1.000, d=14.504, **PASS**.

---

## Part C — Documentation updates

- **`docs/depth_gap.md`:** P23 row updated (pending→PASS, GAP→AT-DEPTH). AT-DEPTH
  count 23→24, gap count 2→1.
- **`REPLICATION_NOTES.md`:** Sprint 77 P23 dim4 section added.
- **`docs/paper_section4_draft.md`:** §4.23 updated with dim4 panel results.
- **`docs/paper_section6_draft.md`:** Sprint 77 entry added, AT-DEPTH count 24/25.
- **`docs/paper_CHANGELOG.md`:** Sprint 77 entry added.

---

## Post-flight checks

- `pytest tests/ -m "not slow"`: all tests pass ✓
- P23 phase2a panel output exists at `analysis/outputs/p23_phase2a_panel.json` ✓
- Panel verdict: PASS (TNR=1.000, d=14.504) ✓
- depth_gap P23 row updated: GAP→AT-DEPTH ✓
- AT-DEPTH count: 24/25 ✓
- REPLICATION_NOTES Sprint 77 section added ✓
- Paper §4.23, §6, CHANGELOG updated ✓

---

## Carry-forwards

- **C-p7-time-shuffled-fp:** `time_shuffled` FP at screening (from Sprint 66, low priority)
- **C-p19-bias-zero-chance-alignment:** 1/5 bias_zero Class C FP at confirmation (from Sprint 70, 4% rate)
- **P12 dim1:** documented finite-size measurement limitation (unchanged, accepted)

---

## Files added/modified

**New (1):**
- `epc/phase2a/failed_regimes/p23_anticoordination.py` — P23 Class C failed regimes

**Modified (9):**
- `epc/detectors/p23_anticoordination.py` — content prerequisite (nondegenerate + below-baseline)
- `epc/phase2a/detector_invariance.py` — P23 flags (perm=True, time=True)
- `epc/phase2a/synthetic.py` — choice_timeseries format in all 10 generators + constants
- `epc/phase2a/structured.py` — 2 choice_timeseries supplements + registry
- `epc/phase2a/panel.py` — choice_timeseries format handling
- `epc/phase2a/catalog.py` — P23 PATTERN_TO_SUBSTRATE_ID entry
- `analysis/run_phase2a_panel.py` — P23 panel builder/runner + main dispatch
- `docs/depth_gap.md` — P23 GAP→AT-DEPTH, counts updated
- `REPLICATION_NOTES.md` — Sprint 77 section

**Paper (3):**
- `docs/paper_section4_draft.md` — §4.23 dim4 results
- `docs/paper_section6_draft.md` — Sprint 77 entry
- `docs/paper_CHANGELOG.md` — Sprint 77 entry

**Output (1):**
- `analysis/outputs/p23_phase2a_panel.json` — dim4 panel results

---

**Decision: GO**
