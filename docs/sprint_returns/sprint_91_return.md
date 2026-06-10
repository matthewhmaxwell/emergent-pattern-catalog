# Sprint 91 Return — P32 dim4 Phase-2a Panel

**Date:** 2026-06-10
**Base HEAD (sprint start):** `59f20cf` (Sprint 90 follow-up)
**Tag:** `v0.91.0`
**Sprint type:** Chat-led design + code-led execution.

## Part A — Wire P32 into Phase-2a panel

Phase-2a infrastructure extended for `task_allocation_timeseries` format:

- **detector_invariance.py:** P32 flags added (permutation_invariant=True,
  time_shuffle_invariant=False). Rationale: mean per-agent entropy is
  invariant under agent relabelling (degenerate-by-construction for
  permutation_shuffled); early-vs-late window comparison requires
  temporal ordering.

- **catalog.py:** P32_response_threshold substrate params + generator +
  PATTERN_TO_SUBSTRATE_ID entry + `_adapt_to_task_allocation` adapter
  (non-task-allocation substrates get random uniform assignments → P32
  screening rejects).

- **synthetic.py:** `_task_allocation_null_history` builder +
  task_allocation_timeseries format case in all 10 Class A generators.
  Permutation-shuffled handler: permutes agent indices (preserves per-agent
  history → degenerate, skipped by perm_inv=True).

- **structured.py:** Two Class B' supplements:
  `single_task_collapse_allocation` (all agents converge to task 0 — entropy
  declines but coverage drops → rejected by late_coverage prerequisite),
  `constant_rebalancing_allocation` (random assignments each step — no
  entropy decline → rejected at screening).

- **panel.py:** task_allocation_timeseries format case in `run_panel`.

- **failed_regimes/p32_specialization.py:** 10 Class C regimes:
  5 no-reinforcement (NoReinforcementModel, ξ=0 → no entropy decline) +
  5 high-forgetting (φ=0.2, ξ=0.005–0.02 → thresholds drift uniformly,
  no stable roles).

- **run_phase2a_panel.py:** `build_p32_positives` (5 seeds, 500 steps),
  `make_p32_detector_fn`, `run_p32` + main dispatch.

- **Invariance flags:** perm_inv=False → True (Sprint 91 correction).
  time_shuffle_inv=False (correct, no change).

## Part B — Panel run

Initial run surfaced two FPs:

1. **permutation_shuffled** (confirmation FP, score=0.700): Mean per-agent
   entropy is invariant under agent relabelling. Fix: corrected
   permutation_invariant to True → substrate skipped.

2. **single_task_collapse_allocation** (screening FP, score=0.600): Entropy
   declines when all agents collapse to one task, but that's not division
   of labor. Fix: added `late_coverage ≥ 0.5` prerequisite at screening
   (literature-grounded: Bonabeau 1996 defines DoL as allocation of
   DIFFERENT tasks).

Both fixes are content prerequisites per Sprint 30 rule (no threshold
relaxation).

Re-run result:

| Class | TNR | n | Detail |
|-------|-----|---|--------|
| Synthetic (A) | 1.000 | 9 evaluated (1 skipped) | permutation_shuffled SKIPPED; all 9 remaining rejected |
| Catalog (B) | 1.000 | 2 (advisory) | single_task_collapse rejected (low coverage); constant_rebalancing rejected (no decline) |
| Failed regimes (C) | 1.000 | 10 | 5 no-reinforcement + 5 high-forgetting, all rejected |
| **Overall** | **1.000** | **21** | **Cohen's d = 2.683** |

**Verdict: PASS.** TNR ≥ 0.95 ✓, d ≥ 1.0 ✓.

## Part C — depth_gap + REPLICATION_NOTES + paper

- `docs/depth_gap.md`: P32 row updated (dim4 pending → PASS, GAP → AT-DEPTH).
  AT-DEPTH count 29 → 30. Gap count 2 → 1 (P12 sole remaining).
- `REPLICATION_NOTES.md`: Sprint 91 P32 dim4 section added.
- `docs/paper_section4_draft.md`: §4.32 updated with dim4 panel results.
- `docs/paper_section6_draft.md`: Sprint 91 entry added.
- `docs/paper_CHANGELOG.md`: Sprint 91 entry added.

## Files changed

**New (1):**
- `epc/phase2a/failed_regimes/p32_specialization.py` — 10 Class C regimes

**Modified (12):**
- `epc/detectors/p32_specialization.py` — screening prerequisite: late_coverage ≥ 0.5
- `epc/phase2a/detector_invariance.py` — P32 flags (perm_inv=True, time_shuffle_inv=False)
- `epc/phase2a/catalog.py` — P32 substrate + generator + adapter
- `epc/phase2a/synthetic.py` — task_allocation_timeseries format in all 10 generators
- `epc/phase2a/structured.py` — 2 supplements for task_allocation_timeseries
- `epc/phase2a/panel.py` — task_allocation_timeseries format case
- `analysis/run_phase2a_panel.py` — P32 build/make/run + dispatch
- `docs/depth_gap.md` — P32 AT-DEPTH, counts updated
- `docs/paper_section4_draft.md` — §4.32 dim4 results
- `docs/paper_section6_draft.md` — Sprint 91 entry
- `docs/paper_CHANGELOG.md` — Sprint 91 entry
- `REPLICATION_NOTES.md` — Sprint 91 section

## Carry-forwards

- **C-p12-dim1**: P12 (RPS λ ∝ √M scaling) remains sole dim1 GAP pattern.
  Documented finite-size measurement limitation (Sprints 54/58/59/63).
- **C-p32-multiseed-low-confirmation**: Only 4/20 seeds reach confirmation
  at default parameters (ξ=0.05, 199 perms, 1000 steps). Sprint 90 carry-forward;
  no change.

## Summary

| Metric | Value |
|--------|-------|
| New files | 1 |
| Modified files | 12 |
| Sprint-specific tests | 25/25 PASS |
| Panel verdict | PASS (TNR=1.000, d=2.683) |
| AT-DEPTH count | 30 / 31 |
| Gap count | 1 (P12) |

**Decision: GO** — P32 dim4 Phase-2a panel PASS. TNR=1.000 across all classes,
Cohen's d=2.683. P32 advances GAP → AT-DEPTH. Two FPs resolved via
literature-grounded content prerequisites (perm_inv correction + late_coverage
screening gate). 30/31 patterns now AT-DEPTH; sole remaining gap is P12 (dim1).
Chain may proceed to Sprint 92 (P30 implement).
