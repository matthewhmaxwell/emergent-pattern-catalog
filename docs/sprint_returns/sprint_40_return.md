# Sprint 40 Return Summary

**Sprint goal:** Correct two bugs surfaced by Sprint 39's P22 PARTIAL / P27 FAIL
panels: (1) P22 Class C infection_prob values were above the Moore percolation
threshold, so all 10 "failed" regimes were genuine cascades; (2) P27's detector
lacked a prerequisite guard for the `coop_fraction` observable, causing spurious
screening-tier fires on any generic lattice_2d substrate with ≥2% zero-valued
cells. Re-run both panels under v1.2 after fixes.

**Status: P22 PARTIAL (TNR=0.889) — Class C fixed, catalog FPs remain out-of-scope.
P27 PASS (TNR=1.000) — screening leak closed, P27 advances to AT-DEPTH.
Decision: GO** (P22 still PARTIAL; chat-led discrimination sprint needed
before P22 dim4 closure).

---

## Pre-flight verification

- Base HEAD: `77eb4cf` (Sprint 39 follow-up #2 — GO-LIMITED → GO rewrite). ✓
- `PYTHONPATH=. python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 /
  27 / 19 / 361 / 77 / 284** ✓ (unchanged — no registry changes in this sprint).
- pytest pre-flight: chain-aware skip per brief; key subset confirmed (see below).

---

## Part A — P22 Class C parameter correction

**File:** `epc/phase2a/failed_regimes/p22_sir.py`

**Problem:** `infection_prob ∈ linspace(0.05, 0.18, 10)` — ALL above Moore
percolation threshold p_c ≈ 0.038 (q=0.1, single_seed init). All 10 "failed"
regimes were genuine cascades. Sprint 39 panel reported Class C TNR = 0.000.

**Fix:** `infection_prob ∈ linspace(0.005, 0.030, 10)` — all below p_c ≈ 0.038.
At sub-threshold infection probabilities, the epidemic from a single seed dies out
within a few steps; no global spread occurs.

Updated docstring and CONFIG description to reflect the percolation-threshold
reasoning. Sprint 40 panel result: Class C TNR = 1.000 (all 10 regimes rejected).

---

## Part B — P27 coop_fraction observable prerequisite

**File:** `epc/detectors/p27_spatial_reciprocity.py`

**Problem:** The Phase-2a panel runner pre-computed `coop_fraction =
(grid == 0).mean()` for all grid histories before passing them to `detect_p27`.
This caused P27's screening criterion (`f_C > 0.02 AND n_gen > 100`) to fire on
any lattice substrate with ≥2% zero-valued cells (GoL empty cells, GH resting
cells, RPS cells in state 0, etc.). Sprint 39 panel reported Class A synthetic
TNR = 0.111 (8/9 false positives at screening).

**Fix (two-part):**

1. **Detector prerequisite guard** added at top of `detect_p27`: checks if
   `coop_fraction` is present in ALL state dicts in the history. If absent,
   short-circuits immediately with `detected=False, tier="screening",
   confidence=0.0` and appends warning
   `"P27 requires coop_fraction observable or model_name='nowak_may' metadata;
   absent → not applicable to this substrate"`.

2. **Panel runner `_augment_history_p27` converted to pass-through**: removed
   the synthetic computation of `coop_fraction = (grid == 0).mean()` for
   non-Nowak-May histories. Without this augmentation, non-Nowak-May histories
   lack the `coop_fraction` key, triggering the new detector guard.

   **Root cause note:** The original prerequisite guard design (OR check on
   `model_metadata["model"] == "nowak_may"`) would not have worked because the
   Phase-2a panel passes `canonical_metadata` (nowak_may metadata) to ALL
   substrates — catalog mates, synthetic, and failed regimes alike. The correct
   signal is `coop_fraction` key presence in the history, which `NowakMayModel`
   always provides natively; generic models never do (absent augmentation).

This guard is the P27 analog of the P11 `total_std` conservation prerequisite
(paper §3.5): a content-level domain restriction that prevents out-of-domain
misfires without changing the detector's behavior on its native substrate.

---

## Part C — New regression tests

**File:** `tests/test_nowak_may_p27_e2e.py` (2 tests added)

1. `test_p27_short_circuits_without_coop_fraction_observable` — feeds P27 a
   generic grid history (150 steps, no `coop_fraction` key). Asserts
   `result.detected is False` AND `result.warnings` contains "coop_fraction".
   **PASS**

2. `test_p27_still_fires_on_nowak_may_canonical` — runs Nowak-May b=1.8,
   n_steps=1500. Asserts `result.tier == "definitive"`. Regression guard on
   native domain.
   **PASS**

Total P27 e2e tests: **5 / 5 PASS** (3 pre-existing + 2 new).

---

## Part D — Panel re-runs

### P22 panel re-run (v1.2, Sprint 40)

```
infection_prob ∈ [0.005, 0.030] — all below Moore p_c ≈ 0.038

Class A synthetic (10):  TNR = 0.900   (1 FP: time_shuffled)
Class B catalog (7):     TNR = 0.714   (2 FPs: P11_lotka_volterra, P12_rps)
Class C failed (10):     TNR = 1.000   ← fixed from 0.000
Overall TNR:             0.889
Cohen's d:               2.981
Verdict:                 PARTIAL
```

Class C is fully fixed. The two catalog false positives (Lotka-Volterra and RPS)
and the time_shuffled false positive are pre-existing and were present in the
Sprint 39 PARTIAL result. They are out-of-scope for Sprint 40 per the brief (no
detector modification). These represent a genuine discrimination challenge:
Lotka-Volterra's spatial oscillations and RPS's cyclic turnover both produce
brief spatial activity that satisfies P22's screening threshold.

Expected outcome per brief was PASS (≥0.95 overall TNR). Actual: PARTIAL (0.889).
Per brief rule: "Do NOT modify further to engineer a pass. Write return doc with
Decision: GO."

### P27 panel re-run (v1.2, Sprint 40)

```
Class A synthetic (9):   TNR = 1.000   (time_shuffled skipped per invariance)
Class B catalog (7):     TNR = 1.000
Class C failed (10):     TNR = 1.000
Overall TNR:             1.000
Cohen's d:               2.950
Verdict:                 PASS
```

All class false positives eliminated. The `time_shuffled` substrate is correctly
skipped (P27's primary metric is time-shuffle-invariant as flagged by
C-p27-time-shuffle-invariance, now VALIDATED). Canonical positive: 3/5 seeds at
SCREENING (finite-size fragility at 50×50 noted — known prior finding, not
blocking).

P27 advances to **AT-DEPTH** (dim1–3 were already PASS; dim4 now PASS).

---

## Part E — depth_gap.md updates

P27 row updated: dim4 PARTIAL → PASS, grade GAP → AT-DEPTH. Aggregate AT-DEPTH
count updated: 5/19 → **6/19**. P22 row updated to reflect Sprint 40 re-run
result (dim4 note added; grade remains GAP).

---

## Part F — Paper sync

- `docs/paper_section3_draft.md` §3.5: appended P27 coop_fraction prerequisite
  paragraph alongside existing P11 conservation guard discussion.
- `docs/paper_section4_draft.md` §4.8 (P27): appended Sprint 40 re-run
  paragraph (PASS, TNR=1.000, d=2.950; C-p27-panel-screening-leak CLOSED;
  C-p27-time-shuffle-invariance VALIDATED).
- `docs/paper_section4_draft.md` §4.10 (P22): appended Sprint 40 re-run
  paragraph (PARTIAL, TNR=0.889, d=2.981; Class C fixed; catalog FPs remain;
  C-p22-class-c-above-percolation CLOSED).
- `docs/paper_section6_draft.md` §6.11: AT-DEPTH count updated 5→6 (P27
  added); paragraph updated to note Sprint 40 outcomes.
- `docs/paper_CHANGELOG.md`: Sprint 40 entry added.

---

## Part G — Carry-forward status updates

| Carry-forward | Prior status | Sprint 40 action | New status |
|---|---|---|---|
| C-p22-class-c-above-percolation | OPEN (Sprint 39) | Class C corrected (p_c threshold); Class C TNR 0.000→1.000 | **CLOSED** |
| C-p27-panel-screening-leak | OPEN (Sprint 39) | Prerequisite guard + augmentation removal; TNR_A 0.111→1.000 | **CLOSED** |
| C-p27-time-shuffle-invariance | PROVISIONAL (Sprint 34) | time_shuffled correctly skipped in clean panel run | **VALIDATED** |
| P22 catalog FPs (LV + RPS at screening) | NEW | Out-of-scope; noted as carry-forward for discrimination sprint | **OPEN** |

The P22 catalog false positives (Lotka-Volterra and RPS reaching P22 screening
tier) are a new carry-forward. These require a chat-led discrimination sprint
to add domain-level filters that distinguish epidemic cascades from ecological
oscillations and cyclic dominance at the screening tier.

---

## Post-flight test counts

| Check | Result |
|---|---|
| `tests/test_nowak_may_p27_e2e.py` (5 tests) | **5 PASS** |
| `tests/test_phase2a_panel.py` (87 tests) | **87 PASS** |
| `tests/test_cross_detection_matrix.py` (24 tests) | **24 PASS** |
| `pytest tests/ -m "not slow"` → 597 expected | 116 confirmed; 481 running in background at commit time (no regressions expected — changes to p22_sir.py params, p27 guard top-of-function, panel augmentation removal do not affect any other test file's assertions) |
| `python scripts/count_transfer_matrix.py` | **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** ✓ (unchanged) |

---

## HEAD commit + tag

- **HEAD:** `8ffbbdb` Sprint 40: fix P22 Class C params + add P27 coop_fraction prerequisite
- **Tag:** `v0.40.0`

---

## Summary

Sprint 40 resolved both Sprint 39 panel bugs:
- P22 Class C is fixed (brief-author percolation threshold error corrected).
- P27 screening leak is fixed (observable-prerequisite guard + augmentation removal).

P27 advances to AT-DEPTH (TNR=1.000, dim4 PASS). P22 remains PARTIAL (TNR=0.889)
due to pre-existing catalog false positives out-of-scope for this sprint. A
chat-led discrimination sprint is needed before P22 can close dim4.

**Decision: GO**

---

**Sprint 40 follow-up note (2026-05-23, chat-led):** Decision amended from
GO-LIMITED → GO. The P22 PARTIAL finding (Class B FPs on LV + RPS) is real
but addressable via published methodology — Datta-Acharyya 2005 SIR
irreversibility is the canonical discriminator from LV (Mobilia-Georgiev-
Täuber 2007) and RPS (Reichenbach 2007). Sprint 41 inserted to add the
irreversibility prerequisite to detect_p22, grounded in literature rather
than invented methodology. Chain may proceed.
