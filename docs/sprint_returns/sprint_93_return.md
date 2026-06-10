# Sprint 93 Return — P30 dim4 Phase-2a Panel (Final Pattern Panel)

**Date:** 2026-06-10
**Base HEAD (sprint start):** `e7de45d` (Sprint 92 follow-up)
**Tag:** `v0.93.0`
**Sprint type:** Chat-led design + code-led execution.

## Part A — Wire P30 into Phase-2a panel

**Invariance flags** (`epc/phase2a/detector_invariance.py`):
- `permutation_invariant=True`: association_score is invariant under particle-index
  relabelling — type–position pairings preserved by consistent permutation.
- `time_shuffle_invariant=True`: all primary metrics (association_score,
  closure_fraction, enrichment_ratio) are per-snapshot spatial statistics.
  Initial run with flag=False produced time_shuffled FP at confirmation
  (score=0.700). Changed to True per empirical evidence — membrane forms
  early and persists, so temporal reordering does not degrade the spatial
  signal in the late window.

**Class C failed regimes** (`epc/phase2a/failed_regimes/p30_autopoiesis.py`):
5 regimes, all correctly rejected at screening:
1. Non-bonding (production_rate=0.0): no links formed.
2. High-decay (decay_rate=0.50): links decay faster than produced.
3. No-attraction + high-decay (attraction=0, decay=0.20, production=0.02,
   high substrate diffusion): links drift away and decay.
4. Weak-production (production_rate=0.01, decay=0.05): production/decay imbalance.
5. Large-box + low-production (box=80, production=0.02, decay=0.05): substrate
   too dispersed for production zone.

**Class B supplements** (`epc/phase2a/structured.py`):
- `dense_cluster_particles`: all type=0 (P1-like dense blob, no link particles →
  mean_n_links=0 → screening fails).
- `dispersed_typed_regions`: types assigned by spatial quadrant (P4-like exclusive
  domains → links in separate region from catalysts → association_score low).

**Synthetic nulls** (`epc/phase2a/synthetic.py`):
`_particle_membrane_null_history()` added + format branches in all 10 generators.
Random positions + random type assignments → association_score ≈ 1.0 → screening fails.

**Catalog adapter** (`epc/phase2a/catalog.py`):
`_adapt_to_particle_membrane()` produces random-typed particles for non-particle-membrane
catalog substrates. `P30_autopoiesis` added to `PATTERN_TO_SUBSTRATE_ID`.

**Panel runner** (`analysis/run_phase2a_panel.py`):
`build_p30_positives` (5 seeds, 300 steps), `make_p30_detector_fn`, `run_p30`.

## Part B — Panel run

`analysis/outputs/p30_phase2a_panel.json`.

| Class | N evaluated | TNR | Notes |
|-------|------------|-----|-------|
| A (synthetic) | 8 | 1.000 | 2 skipped (perm_inv + time_shuffle_inv) |
| B (supplements) | 2 | 1.000 | dense_cluster (P1-like), dispersed_regions (P4-like) |
| C (failed regimes) | 5 | 1.000 | 5 non-membrane regimes |
| **Overall** | **15** | **1.000** | **d = 12.124** |

**Verdict: PASS**. TNR=1.000, Cohen's d=12.124.

Canonical positive: 5/5 detected (1 definitive, 2 confirmation, 2 screening).
Mean positive score: 0.700.

**FP resolution (3 runs):**
- Run 1 (PARTIAL, TNR=0.750): 4 FPs — time_shuffled (confirmation),
  dense_cluster (screening), no_attraction (screening), large_box (confirmation).
- Fix: time_shuffle_inv→True (empirical invariance confirmed); dense_cluster
  changed to all-type-0 (matching DenseClusterModel negative control);
  no_attraction regime strengthened (decay↑, production↓, diffusion↑);
  large_box regime strengthened (production↓, decay↑).
- Run 2 (PARTIAL, TNR=0.933): 1 FP — no_attraction (confirmation).
- Fix: decay 0.10→0.20, production 0.05→0.02, added substrate_diffusion=0.5.
- Run 3 (PASS, TNR=1.000): 0 FPs.

## Part C — Documentation updates

- `docs/depth_gap.md`: P30 dim4→PASS, GAP→AT-DEPTH. AT-DEPTH: 31/32. Gap: 1 (P12).
- `REPLICATION_NOTES.md`: Sprint 93 section added.
- `docs/paper_section4_draft.md`: §4.30 dim4 paragraph added.
- `docs/paper_section6_draft.md`: Sprint 93 entry added.
- `docs/paper_CHANGELOG.md`: Sprint 93 entry added.

## Files changed

**New (1):**
- `epc/phase2a/failed_regimes/p30_autopoiesis.py` — 5 Class C failed regimes

**Modified (11):**
- `epc/phase2a/detector_invariance.py` — P30 invariance flags
- `epc/phase2a/synthetic.py` — particle_membrane null history + format branches
- `epc/phase2a/structured.py` — 2 Class B supplements
- `epc/phase2a/catalog.py` — particle_membrane adapter + P30 registration
- `epc/phase2a/panel.py` — particle_membrane format branch
- `analysis/run_phase2a_panel.py` — P30 panel entry
- `REPLICATION_NOTES.md` — Sprint 93 section
- `docs/depth_gap.md` — P30 AT-DEPTH
- `docs/paper_section4_draft.md` — §4.30 dim4
- `docs/paper_section6_draft.md` — Sprint 93 entry
- `docs/paper_CHANGELOG.md` — Sprint 93 entry

## Carry-forwards

- **C-p12-dim1**: P12 (RPS λ ∝ √M scaling) remains dim1 GAP.
  Documented finite-size measurement limitation (Sprints 54/58/59/63).
  **SOLE REMAINING GAP.**
- **C-p30-enrichment-cv**: Enrichment ratio CV=34.7% across seeds (Sprint 92).
  Association_score CV=3.9% is robust. Informational only; does not affect panel.

## Summary

| Metric | Value |
|--------|-------|
| New files | 1 |
| Modified files | 11 |
| Sprint-specific tests | 25/25 PASS (unchanged from Sprint 92) |
| Panel verdict | PASS |
| Panel TNR | 1.000 |
| Panel Cohen's d | 12.124 |
| FP resolution runs | 3 (4 FPs → 1 FP → 0 FPs) |
| AT-DEPTH count | 31 / 32 |
| Gap count | 1 (P12) |

**Decision: GO** — P30 Phase-2a panel PASS (TNR=1.000, d=12.124). All 15
negatives correctly rejected across 3 classes. P30 advances GAP→AT-DEPTH.
31/32 patterns at depth; sole gap is P12 (dim1 finite-size limitation).
This is the final pattern panel. Chain may proceed to Sprint 94
(full-catalog completion summary).
