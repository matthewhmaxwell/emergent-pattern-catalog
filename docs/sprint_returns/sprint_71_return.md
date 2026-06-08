# Sprint 71 Return — Milestone B Wave 1 Checkpoint + PAUSE

**Date:** 2026-06-08
**Base HEAD (sprint start):** `c43932a` (Sprint 70 follow-up)
**Sprint goal:** Summarize Wave 1 (P7, P17, P19), recount AT-DEPTH + implemented, HALT for operator review.
**Tag:** `v0.71.0`
**Sprint type:** Chat-led checkpoint. Documentation-only.

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `c43932a` ✓
2. **Working tree clean:** ✓

---

## Part A — Recount

Source: `docs/depth_gap.md`

- **Implemented pattern count:** 22
- **AT-DEPTH count:** 21 / 22
- **Gap count:** 1 (P12 — dim1 documented finite-size measurement limitation)

### Wave 1 per-pattern outcome

| Pattern | Sprints | dim1 | dim2 | dim3 | dim4 | AT-DEPTH? |
|---------|---------|------|------|------|------|-----------|
| P7 (lane formation) | 65–66 | PASS | PASS | PASS | PASS (TNR=0.955, d=6.932) | **Yes** |
| P17 (collective sensing) | 67–68 | PASS | PASS | PASS | PASS (TNR=1.000, d=11.117) | **Yes** |
| P19 (emergent leadership) | 69–70 | PASS | PASS | PASS | PASS (TNR=0.960, d=5.418) | **Yes** |

All three Wave 1 patterns reached AT-DEPTH. None stalled at dim4.

---

## Part B — Wave 1 summary doc

Written to `docs/milestone_b_wave1_summary.md`. Contains:
- Full outcome table (model/detector files, dim1 reproduction results + tolerance verdicts, dim4 panel TNR/d, AT-DEPTH status)
- Content prerequisites added during panels with literature grounding
- Open carry-forwards from Wave 1
- Confirmation that each dim1 reproduction hits a real published anchor

---

## Post-flight checks

- `pytest tests/test_orchestration.py tests/test_transfer_matrix_counts.py -m "not slow"`: **85 passed** ✓
- No code changes made (documentation-only sprint) ✓
- No Wave 2 / new-pattern work performed ✓

---

## Carry-forwards (unchanged from Sprint 70)

- **C-p7-time-shuffled-fp:** `time_shuffled` FP at screening (low priority)
- **C-p19-bias-zero-chance-alignment:** 1/5 bias_zero Class C FP at confirmation (4% rate, within tolerance)
- **C-p19-abrupt-saturation:** accuracy saturation more abrupt than Couzin Fig. 2a (not a validation failure)
- **C-p19-te-vs-pull:** TE collapse in converged regime (documented architectural decision)
- **P12 dim1:** sole remaining GAP — documented finite-size measurement limitation (Sprints 54/58/59/63)

---

## Wave 1 outcome recap

**3/3** new patterns (P7, P17, P19) reached AT-DEPTH. Total: **22 implemented**, **21 AT-DEPTH**. The new-pattern pipeline (model → detector → dim1 reproduction → dim4 panel) produced literature-faithful results across all three self-propelled-particle patterns. Full details in `docs/milestone_b_wave1_summary.md`.

---

**Decision: NO-GO**

Milestone B Wave 1 complete; chain paused for operator review before Wave 2 (P24/P26/P23). Deliberate review gate, not a failure.
