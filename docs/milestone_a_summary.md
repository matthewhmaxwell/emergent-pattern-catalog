# Milestone A Summary — Implemented-Pattern Depth Closure

**Date:** 2026-05-30 (Sprint 64)
**Milestone A goal:** Bring all 19 implemented patterns to AT-DEPTH status per the depth rubric (`docs/depth_rubric.md`).
**Actual outcome:** **18 / 19 AT-DEPTH.** One pattern (P12) remains GAP with a documented finite-size measurement limitation.

---

## Milestone A Gap Closure Table

The Milestone A queue (Sprints 59–63) targeted 6 remaining depth gaps across 5 patterns:

| Pattern | Gap dimension | Sprint | Outcome | Detail |
|---------|--------------|--------|---------|--------|
| P2 (MIPS) | dim2 (multi-seed) | 60 | **CLOSED** | 20-seed: two_phase_score=0.1134±0.0790 (CV=69.7%), density-speed r=−0.9585±0.0196 (CV=2.1%). AT-DEPTH. |
| P21 (Polarization) | dim2 (multi-seed) | 60 | **CLOSED** | 20-seed: cluster_count=1.90±0.31 (CV=16.2%, median=2, 18/20→2 clusters). AT-DEPTH. |
| P22 (Cascade) | dim2 (multi-seed) | 60 | **CLOSED** | 20-seed: wavefront speed=0.4606±0.0163 (CV=3.5%, 19/20 valid). AT-DEPTH. |
| P1 (Aggregation) | dim4 (negative sweep) | 61 | **CLOSED** | Phase-2a panel v1.2 re-run: TNR=1.000 (↑ from 0.731). Class C calibration fixed (grid_size 32→50, threshold range corrected below critical ≈0.13). AT-DEPTH. |
| P8 (Jamming) | dim4 (negative sweep) | 62 | **CLOSED** | Phase-2a panel v1.2 re-run: TNR=1.000 (↑ from 0.714). Class C correction (density range below jamming onset ρ_c≈0.10). AT-DEPTH. |
| P12 (Cyclic dominance) | dim1 (replication) | 59, 63 | **DOCUMENTED-LIMITATION** | Four attempts (Sprints 54, 58, 59, 63) across L=100 ACF and L=200 FFT estimators; slopes 0.37, 0.107, 0.244, 0.161 — all outside [0.4, 0.6]. Accepted as finite-size measurement limitation. Grade stays GAP. |

---

## Patterns NOT at AT-DEPTH

**P12 — Cyclic dominance (spatial RPS):** dim1 PARTIAL.

The Reichenbach-Mobilia-Frey (2007) λ ∝ M^(1/2) wavelength scaling law could not be quantitatively reproduced within this project's compute and lattice-size budget. Per the Sprint 30 depth rubric, dim1 requires reproduction of a specific published figure/table with stated tolerance. Four sprint-length attempts established that the analytical formula is valid only for M/M_c ∈ [0.44, 1.00], a range providing ~2 usable data points at any lattice size tested (L=100, L=200) — insufficient for a reliable slope determination.

This is a **genuine measurement limitation**, not a force-pass candidate. The Sprint 30 rule is clear: PARTIAL scores must not be upgraded without meeting the rubric threshold. P12's validation rests on:
- Phase-2a panel PASS (TNR=1.000, Cohen's d=+∞)
- dim2 multi-seed PASS (20 seeds, CV=20%)
- dim3 methods note PASS
- dim4 negative sweep PASS
- Qualitative spiral presence confirmed across all tested M values at L=100 and L=200

The scaling-law gap is a documented finite-size measurement limitation, not a validation gap.

---

## Open Carry-Forwards After Milestone A

| ID | Description | Status | Priority |
|----|------------|--------|----------|
| C1 | Phase-2a uniformity sprint (dim4 gaps) | Effectively closed — all 19 patterns now have dim4 PASS | — |
| C2 | P12 dim1 replication gap (λ ∝ √M) | Closed-as-documented-limitation (Sprint 63) | — |
| C3 | P12 RPS wavelength scaling not replicated | Closed-as-documented-limitation (Sprint 63) | — |
| C4 | P14 BTW multi-seed dispersion | Closed (Sprint 55) | — |
| C5 | P21/P22 methods notes thin | Closed (Sprint 57) | — |
| C-p21-time-shuffled-fp | P21 `time_shuffled` FP at CONFIRMATION tier | Open (low priority) | Low |
| C-p9-constant-field | P9 `constant_field` trivial-sync Class A FP | Open (low priority) | Low |
| C-p10-perm-shuffled-fp | P10 `permutation_shuffled` FP at screening | Open (low priority) | Low |
| C-p14-class-c-borderline | P14 borderline at p_diss=0.350 | Open (low priority) | Low |

All named carry-forwards from the original audit (C1–C5) are closed. The remaining open items are minor Phase-2a panel edge cases that do not affect AT-DEPTH status.

---

## Milestone A Timeline

| Sprint | Work | AT-DEPTH delta |
|--------|------|----------------|
| 28 | Depth rubric + audit | 0/19 baselined |
| 32–49 | Phase-2a panels + invariance flags | 0→10/19 |
| 50–54 | dim1 closures (P11, P22, P2, P21, P12 attempt) | 10→11/19 |
| 55–58 | dim2 closures + methods notes + P12 attempts | 11→16/19 |
| 59 | P12 dim1 near-M_c attempt (FAIL) | 16/19 |
| 60 | P2/P21/P22 dim2 batch closure | 16→18/19 (note: Sprint 60 closed the last 3 dim2 gaps but P1/P8 dim4 were already closed implicitly by Sprint 49+61+62 sequence) |
| 61 | P1 dim4 Class C correction | → P1 AT-DEPTH |
| 62 | P8 dim4 Class C correction | → P8 AT-DEPTH, 18/19 |
| 63 | P12 dim1 final attempt — accepted limitation | 18/19 (unchanged) |
| 64 | This summary (no code changes) | 18/19 final |

---

## Recommended Next Step: Milestone B

Milestone B is the implementation of the **13 unimplemented patterns** from the original 32-pattern catalog:

P4, P7, P16, P17, P19, P20, P23, P24, P25, P26, P29, P30, P32

This would expand the catalog from 19 to 32 implemented patterns and is a substantially larger effort than Milestone A. Each new pattern requires a canonical model, detector, Phase-2a panel, and multi-dimensional validation.

**This is a recommendation awaiting human go-ahead.** The autonomous sprint chain is intentionally paused at this boundary. The human should review this summary and the current project state before authorizing Milestone B work.
