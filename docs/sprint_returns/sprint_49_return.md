# Sprint 49 Return — Invariance-flag batch update: P1, P2, P5, P6, P8, P21

**Date:** 2026-05-25
**Base HEAD (sprint start):** `6630f0b`
**Sprint goal:** Update `epc/phase2a/detector_invariance.py` with empirically-observed invariance flags for P1, P2, P5, P6, P8, P21 (each one's pattern of Class A FPs across Sprints 39–48 reveals which invariance they actually have). Re-run all 6 panels under v1.2. Close ~7 Class A degenerate-FP carry-forwards in a single batch.
**Tag:** `v0.49.0`

---

## Part A — Invariance-flag updates

**File:** `epc/phase2a/detector_invariance.py`

Flags updated per empirical FP evidence and mathematical grounding (Sprint 30 rule: panel metadata only; no detector logic changes):

| Pattern | Before (perm_inv, time_inv) | After (perm_inv, time_inv) | Carry-forward closed |
|---------|---------------------------|---------------------------|----------------------|
| P1  | (False, False) | (False, True)  | C-p1-time-shuffle-fp |
| P2  | (absent → False, False) | (True, False)  | C-p2-perm-shuffled-fp |
| P5  | (True, False)  | (True, True)   | C-p5-time-shuffle-fp |
| P6  | (False, False) | (True, True)   | C-p6-time-shuffle-fp |
| P8  | (absent → False, False) | (True, True)   | C-p8-perm-shuffled-fp, C-p8-time-shuffle-fp |
| P21 | (False, False) | (True, False)  | C-p21-perm-shuffled-fp |
| P10 | (False, False) | (False, False) | — (no change; adapter artifact, not invariance) |

**Mathematical grounding (literature-first principle):**
- **P1** `time_shuffle_inv=True`: Moran's I is per-frame; each segregated Schelling frame has high I regardless of temporal order.
- **P2** `perm_inv=True`: `two_phase_score` thresholds a coarse density grid; particle identities do not affect spatial density.
- **P5** `time_shuffle_inv=True`: mean polar order parameter φ = |⟨e^{iθ}⟩| per frame; every flocked frame has high φ.
- **P6** `perm_inv=True`, `time_shuffle_inv=True`: angular momentum |L| = |Σ r_i × v_i|/N is a sum over particles; every milled frame has high |L|.
- **P8** `perm_inv=True`, `time_shuffle_inv=True`: `stopped_fraction` is a time-average; the NaSch model is translation-invariant; the mean is preserved under any temporal reordering.
- **P21** `perm_inv=True`: Hartigan's dip test operates on the sorted value histogram; permuting agent opinion indices leaves the dip statistic unchanged.
- **P10** *no change*: C-p10-perm-shuffled-fp arises from catalog-adapter binarization preserving bimodal phase structure. This is an adapter artifact, NOT mathematical permutation invariance of `local_r`. Flagging perm_inv=True would be the wrong fix; documented in rationale comment.

**Test update:** `tests/test_phase2a_panel.py` parametrize table updated to match new flag values (P1, P2 added, P5, P6, P8 added, P21 corrected). `test_panel_runs_all_class_a_for_non_invariant_detector` updated to use `P15` (both flags False) instead of P1 (now time_shuffle_inv=True). All 157 key tests pass.

---

## Part B — Panel re-runs (v1.2)

All 6 panels re-run. P10 NOT re-run (flags unchanged per brief).

| Pattern | Overall TNR before | Overall TNR after | syn TNR before | syn TNR after | Cohen's d before | Cohen's d after | Verdict before | Verdict after |
|---------|-------------------|-------------------|----------------|---------------|-----------------|-----------------|----------------|---------------|
| P1  | 0.704 | 0.731 | 0.800 | 0.889 | 1.624 | 1.740 | PARTIAL | PARTIAL |
| P2  | 0.958 | 1.000 | 0.900 | 1.000 | 3.401 | 4.245 | PASS | PASS |
| P5  | 0.957 | 1.000 | 0.889 | 1.000 | 4.987 | +inf  | PASS-with-weakness | PASS |
| P6  | 0.958 | 1.000 | 0.900 | 1.000 | 5.087 | +inf  | PASS | PASS |
| P8  | 0.652 | 0.714 | 0.800 | 1.000 | 1.751 | 1.772 | PARTIAL | PARTIAL |
| P21 | 0.913 | 0.955 | 0.800 | 0.889 | 4.543 | 5.487 | PARTIAL | PASS-with-weakness |

**Commentary per pattern:**

**P1:** `time_shuffled` SKIPPED (time_shuffle_invariant=True). syn TNR 0.800→0.889. Residual FPs: `linear_gradient` at CONFIRMATION (C-p1-linear-gradient-fp; Moran's I responds to spatial gradient — different issue from invariance) and Class C fai=0.400 (C-p1-class-c-subthreshold-fp; threshold ∈ [0.161, 0.250] above empirical critical threshold at density=0.9 — calibration issue, out-of-scope). Verdict unchanged: PARTIAL.

**P2:** `permutation_shuffled` SKIPPED (permutation_invariant=True). syn TNR 0.900→1.000. All FPs closed. Overall TNR 0.958→1.000. Clean PASS.

**P5:** `time_shuffled` SKIPPED (time_shuffle_invariant=True, in addition to already-skipped `permutation_shuffled`). syn TNR 0.889→1.000. Overall TNR 0.957→1.000, Cohen's d 4.987→+inf. PASS-with-weakness → clean PASS.

**P6:** Both `permutation_shuffled` and `time_shuffled` SKIPPED (both flags True). syn TNR 0.900→1.000. Overall TNR 0.958→1.000, Cohen's d 5.087→+inf. PASS → cleaner PASS.

**P8:** Both `permutation_shuffled` and `time_shuffled` SKIPPED (both flags True). syn TNR 0.800→1.000. Class C fai=0.400 unchanged (C-p8-class-c-near-onset: 6/10 regimes at rho ≥ 0.1167 are at or above NS jamming onset at p=0.3; a restricted Class C sweep to rho ∈ linspace(0.01, 0.09, 10) is the fix). Overall TNR 0.652→0.714. Verdict unchanged: PARTIAL.

**P21:** `permutation_shuffled` SKIPPED (permutation_invariant=True). `time_shuffled` NOT skipped (time_shuffle_invariant=False; HK trajectory is NOT time-order-independent — pre-convergence unimodal steps differ meaningfully from post-convergence bimodal steps). `time_shuffled` FP at CONFIRMATION (0.850) persists: C-p21-time-shuffled-fp remains open (convergence-timing issue). syn TNR 0.800→0.889. Overall TNR 0.913→0.955, Cohen's d 4.543→5.487. Verdict: PASS-with-weakness (advances from PARTIAL). P21 dim4 PARTIAL → PASS.

---

## Part C — depth_gap.md + REPLICATION_NOTES updates

| Pattern | dim4 before | dim4 after | grade before | grade after | AT-DEPTH Δ |
|---------|-------------|------------|--------------|-------------|------------|
| P1  | PARTIAL | PARTIAL          | GAP | GAP | 0 |
| P2  | PASS (TNR=0.958) | PASS (TNR=1.000) | GAP | GAP | 0 |
| P5  | PASS (PASS-w-w) | PASS (clean)     | AT-DEPTH | AT-DEPTH | 0 |
| P6  | PASS (PASS-w-w) | PASS (clean)     | GAP | GAP | 0 |
| P8  | PARTIAL | PARTIAL          | GAP | GAP | 0 |
| P21 | PARTIAL | PASS (PASS-w-w)  | GAP | GAP | 0 |

**AT-DEPTH count: 10 / 19** (unchanged). P21 dim4 advances from PARTIAL to PASS, but dims 1–3 remain PARTIAL, so P21 does not reach AT-DEPTH.

`docs/depth_gap.md` rows for P1, P2, P5, P6, P8, P21 updated with Sprint 49 re-run results.
`REPLICATION_NOTES.md` Sprint 49 section appended with full flag table + panel result table.

---

## Part D — Paper sync

- `docs/paper_section3_draft.md`: §3.5 appended — Sprint 49 invariance-flag batch rationale (literature grounding, P10 exception, empirical outcomes).
- `docs/paper_section4_draft.md`: §4.27 appended — per-pattern Sprint 49 panel re-run narrative for all 6 patterns.
- `docs/paper_section6_draft.md`: Sprint 49 AT-DEPTH +0 paragraph appended.
- `docs/paper_CHANGELOG.md`: Sprint 49 entry appended.

---

## Part E — Carry-forwards

**CLOSED (7 total):**
- C-p1-time-shuffle-fp — `time_shuffled` SKIPPED via time_shuffle_invariant=True
- C-p2-perm-shuffled-fp — `permutation_shuffled` SKIPPED via permutation_invariant=True
- C-p5-time-shuffle-fp — `time_shuffled` SKIPPED via time_shuffle_invariant=True
- C-p6-time-shuffle-fp — `time_shuffled` SKIPPED via time_shuffle_invariant=True
- C-p8-perm-shuffled-fp — `permutation_shuffled` SKIPPED via permutation_invariant=True
- C-p8-time-shuffle-fp — `time_shuffled` SKIPPED via time_shuffle_invariant=True
- C-p21-perm-shuffled-fp — `permutation_shuffled` SKIPPED via permutation_invariant=True

**REMAIN OPEN (structural issues, not invariance-fixable):**

| ID | Pattern | Description | Priority | Recommended effort |
|----|---------|-------------|----------|--------------------|
| C-p1-linear-gradient-fp | P1 | `linear_gradient` fires at CONFIRMATION (0.700); Moran's I correctly responds to spatial gradient structure | Low | S (tighter spatial-gradient guard or suppress gradient substrates for Schelling-type detectors) |
| C-p1-class-c-subthreshold-fp | P1 | fai TNR=0.400; Class C low-threshold regimes above empirical critical threshold at density=0.9 | Low | M (recalibrate Class C sweep below empirical threshold) |
| C-p8-class-c-near-onset | P8 | fai TNR=0.400; NS jamming onset at rho≈0.12 (p=0.3) overlaps Class C sweep; 6/10 regimes fire | Medium | S (restrict sweep to rho ∈ linspace(0.01, 0.09, 10)) |
| C-p10-perm-shuffled-fp | P10 | Catalog-adapter binarization artifact; `permutation_shuffled` fires at SCREENING via preserved bimodal phase structure | Low | M (refine catalog adapter to avoid bimodal preservation for chimera-state input) |
| C-p21-time-shuffled-fp | P21 | `time_shuffled` FP at CONFIRMATION (0.850); HK pre-convergence unimodal steps reduce persistence in shuffled trajectory | Low | S (extend canonical positive run length so pre-convergence steps are a smaller fraction, or use post-convergence-only window) |
| C-class-a-constant-field-trivial-sync | Various | `constant_field` Class A substrate triggers trivial-sync FPs on several detectors | Low | M (separate methodology sprint) |

---

## Part F — Post-flight

```
PYTHONPATH=. python3.12 -m pytest tests/test_orchestration.py tests/test_phase2a_panel.py -x -q
157 passed in 11.58s
```

No regressions. Test count +3 (new parametrize entries for P2, P6, P8 in invariance flags table).

---

**Decision: GO**
