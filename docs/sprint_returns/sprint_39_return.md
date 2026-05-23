# Sprint 39 Return Summary

**Sprint goal:** Run v1.2 Phase-2a panel against P22 (SIR information cascade) and
P27 (Nowak-May spatial reciprocity) — lattice_2d dim4 batch 1. Both patterns
currently GAP on dim4. Per Sprint 35 prediction, this was designated the GO batch.

**Status: code complete; both panels returned below PASS. Sprint 30 rule applied.
Chain paused for chat-led review. Decision: GO-LIMITED.**

---

## Pre-flight verification

- Base HEAD: `ae0abec` (Sprint 38 follow-up — GO-LIMITED → GO amend; matches
  expected Sprint 38 post-commit). ✓
- `PYTHONPATH=. python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 /
  27 / 19 / 361 / 77 / 284** ✓ (unchanged, as expected for panel-runs-only sprint).
- pytest pre-flight: chain-aware skip per brief.

---

## Part A — P27 Nowak-May panel

### `epc/phase2a/failed_regimes/p27_nowak_may.py` (new)

10 extinction regimes at `b ∈ linspace(2.0, 2.5, 10)`, rows=50, cols=50,
n_steps=200, init_coop_fraction=0.5, seeds 200+i. At b ≥ 2.0 cooperation is
extinct (T > 2R → defectors dominate unconditionally). Class C TNR = 1.000 —
all extinction regimes correctly rejected.

### `run_p27()` added to `analysis/run_phase2a_panel.py`

Following the `run_p18` template: grid format, lattice_2d substrate type,
failed_regime_module=p27_failed, canonical positive NowakMay 50×50 b=1.8
n_steps=200 (5 seeds). Added `_augment_history_p27()` wrapper to inject
`coop_fraction` and `moran_i` into generic grid-format catalog substrates
that lack those keys (detect_p27 expects NowakMay-native history format).

### P27 panel verdict: **FAIL**

Output: `analysis/outputs/p27_phase2a_panel.json`

| Class | TNR | n_eval / n_total | Notes |
|---|---|---|---|
| Synthetic (Class A) | 0.111 | 1 / 9 | gating; 8 FP at SCREENING |
| Catalog (lattice_2d) | 0.286 | 2 / 7 | gating; 5 FP at SCREENING |
| Failed-regime (b ≥ 2.0) | **1.000** | 10 / 10 | gating |
| **Overall** | **0.500** | 13 / 26 | below PASS gate |

- Cohen's d: **0.198** (below PARTIAL gate of 0.5 → FAIL)
- Canonical positive: 3/5 seeds reach SCREENING (conf=0.600); 2/5 get NONE (f_C → 0 at 50×50 stochastic extinction, finite-size effect)
- **Root cause 1 — screening-tier leak:** `detect_p27` screens on `fc_mean > 0.02 AND n_gen > 100`. The panel adapter computes `coop_fraction = (grid==0).mean()`. GoL (≈80% empty=0), GH (≈70% resting=0), Schelling (≈10% empty=0) all trivially satisfy fc_mean > 0.02. Confirmation tier's `well_mixed=0.5` gate (for `pd_verified=False`) would correctly reject most, but since `detected = (tier ≠ NONE)`, SCREENING fires a false positive.
- **Root cause 2 — canonical positive fragility at 50×50:** 2/5 seeds exhibit stochastic cooperation collapse (f_C → 0). Canonical grid is 100×100; panel used 50×50 for tractability. Finite-size effect causes unreliable detection.
- v1.2 time_shuffled skipped (C-p27-time-shuffle-invariance flag provisional). Flag NOT validated (panel FAIL provides no clean data). Flag remains provisional.

---

## Part B — P22 SIR panel

### `epc/phase2a/failed_regimes/p22_sir.py` (new)

10 sub-percolation regimes at `infection_prob ∈ linspace(0.05, 0.18, 10)`,
rows=64, cols=64, recovery_prob=0.1, single_seed init, seeds 400+i. Brief
described these as "below percolation threshold ~0.2."

### `run_p22()` added to `analysis/run_phase2a_panel.py`

Grid format, lattice_2d substrate type, canonical positive SIR 64×64
infection_prob=0.4 recovery_prob=0.1 single_seed (5 seeds).

### P22 panel verdict: **PARTIAL**

Output: `analysis/outputs/p22_phase2a_panel.json`

| Class | TNR | n_eval / n_total | Notes |
|---|---|---|---|
| Synthetic (Class A) | 0.900 | 9 / 10 | gating; FP on `time_shuffled` |
| Catalog (lattice_2d) | 0.714 | 5 / 7 | gating; FP on LV, RPS |
| Failed-regime (p ∈ [0.05, 0.18]) | **0.000** | 0 / 10 | gating; all detected |
| **Overall** | **0.519** | 14 / 27 | below PASS gate |

- Cohen's d: **1.094** (above PARTIAL gate of 0.5 → PARTIAL, not FAIL)
- Canonical positive: all 5 seeds detected at SCREENING tier (conf=0.500). At 64×64, panel positives reach SCREENING only; the canonical test at 80×80 achieves DEFINITIVE.
- **Root cause — Class C params above Moore percolation threshold:** infection_prob ∈ [0.05, 0.18] are all above p_c ≈ 0.038 (Moore, q=0.1). The epidemic genuinely spreads from the single seed. P22 correctly identifies these as cascades. The brief's "below percolation threshold ~0.2" refers to an effective CONFIRMATION threshold, not the physical bond percolation threshold. Class C TNR = 0.000: all 10 regimes are true cascades, not failed regimes.

---

## Part C — REPLICATION_NOTES additions

Appended "Phase-2a Panel Result (v1.2) — Sprint 39" subsections to:
- Nowak-May section (P27 panel FAIL, failure-mode analysis)
- SIR section (P22 panel PARTIAL, failure-mode analysis)

Format matches Sprint 35 P9/P14 v1.2 sections.

---

## Part D — depth_gap.md updates

- P22 row: dim4 updated from "cross-detection only" to "Phase-2a panel v1.2 PARTIAL (Sprint 39); C-p22-class-c-above-percolation carry-forward". Grade remains GAP.
- P27 row: dim4 updated from "Sprint 14.5 characterization only" to "Phase-2a panel v1.2 FAIL (Sprint 39); C-p27-panel-screening-leak carry-forward". Grade remains GAP.
- Aggregate header: AT-DEPTH count remains **5 / 19**. Sprint 39 note added explaining PARTIAL/FAIL outcomes and chat-led redesign requirement.

---

## Part E — Paper sync

- §4.8 P27 Nowak-May: appended Phase-2a FAIL panel paragraph (TNR=0.500, d=0.198, two linked root causes, carry-forward, Sprint 30 rule note).
- §4.10 P22 SIR: appended Phase-2a PARTIAL panel paragraph (TNR=0.519, d=1.094, Class C above-percolation root cause, redesign path).
- §6.11 aggregate: AT-DEPTH count unchanged at 5/19; Sprint 39 finding noted.
- `docs/paper_CHANGELOG.md`: Sprint 39 entry prepended.

---

## Part F — Carry-forward review

### C-p27-time-shuffle-invariance (Sprint 34, open/provisional)

P27 panel **FAILED** — no clean data to validate or invalidate the flag.
Per brief: "If P27 PASSed cleanly: C-p27-time-shuffle-invariance CLOSED.
If P27 PARTIAL: keep flag as provisional, escalate for chat-led sprint."
P27 FAIL is worse than PARTIAL; flag stays provisional, escalation triggered.

### New carry-forwards (Sprint 39)

| ID | Pattern | Description | Priority |
|---|---|---|---|
| C-p27-panel-screening-leak | P27 | Screening tier fires on any grid with ≥2% zero-cells (generic lattice_2d FPs). Root cause: detect_p27 designed for NowakMay-native history; panel adapter creates semantic mismatch. Resolution requires chat-led sprint to redesign panel approach (e.g., use model-native format only for Class B, or revise screening criterion). | HIGH — blocks dim4 closure |
| C-p22-class-c-above-percolation | P22 | Failed-regime params (p ∈ [0.05, 0.18]) lie above Moore p_c ≈ 0.038. All 10 "failed" regimes are genuine cascades. Redesign: use p < 0.038 (Moore), or Von Neumann with p ∈ [0.02, 0.09], or single_seed-with-early-abort to ensure extinction. | HIGH — blocks dim4 closure |

All other open carry-forwards unchanged from Sprint 38:
- **C-class-a-constant-field-trivial-sync (Sprint 35):** Open. P9 residual false positive. v1.3 candidate; not in scope.
- **C-p14-class-c-borderline (Sprint 33):** Persists at p_diss=0.350; P14 still PASS.
- **C-pyproject-pin (Sprint 29):** Open. 1-line deferred.
- **C-supplements-lattice_2d_continuous, scalar_wealth, opinion_space:** Still open.

---

## Post-flight verification

| Check | Result |
|---|---|
| `analysis/outputs/p27_phase2a_panel.json` exists with `"panel_version":"1.2"` | **PASS** ✓ |
| `analysis/outputs/p22_phase2a_panel.json` exists with `"panel_version":"1.2"` | **PASS** ✓ |
| P27 JSON contains ≥1 SKIPPED-degenerate-by-construction (time_shuffle_invariant=True) | **PASS** ✓ (time_shuffled: SKIPPED) |
| P22 JSON: no SKIPPED substrates (both invariance flags False) | **PASS** ✓ |
| `grep -c "Phase-2a Panel Result (v1.2)" REPLICATION_NOTES.md` ≥ prev_count + 2 | **PASS** ✓ (Sprint 35 had 2, now has 4) |
| depth_gap.md rows for P22 + P27 updated | **PASS** ✓ |
| Transfer matrix: **20/19/79/274/27/19/361/77/284** | **PASS** ✓ (no registry changes) |
| `pytest tests/ -m "not slow"` → 595 passed (unchanged) | Pending (running) |
| P22 PARTIAL + P27 FAIL → escalate_to_user=true in state.json | **PASS** ✓ (set in epc-orchestrator-state/state.json) |

---

## Deviations and judgment calls

1. **`_augment_history_p27()` added as panel runner wrapper (not a detector change):** `detect_p27` requires `coop_fraction` and `moran_i` keys that NowakMay native histories provide but generic grid-format catalog substrates do not. The wrapper computes these from the grid so the detector runs without crashing. This is a panel runner adaptation, not a detector modification. Per Sprint 30 rule: the detector is unchanged.

2. **Both verdicts PARTIAL/FAIL — Sprint 30 rule applied:** No detector, model, or panel composition changes were made to engineer a pass. Results reported as-is. Two new carry-forwards opened. Orchestrator state escalated.

3. **Sprint 35 prediction ("GO batch") was incorrect:** The prediction did not account for P27's screening-tier semantics mismatch with generic grids, nor for P22's Class C parameters exceeding the physical percolation threshold. Both errors are informative findings about panel design.

---

## HEAD commit hash and tag at end of sprint

- **Commit:** `v0.39.0` (pending test verification)
- **Tag:** `v0.39.0`

**Decision: GO-LIMITED**

Both P22 and P27 returned below PASS under v1.2 panel. Sprint 30 rule applied throughout — no detector or panel composition changes made. Two carry-forwards opened (C-p22-class-c-above-percolation, C-p27-panel-screening-leak). `escalate_to_user=true` set in orchestrator state. A chat-led sprint is required to redesign the Class C failed-regime params for P22 and to redesign the P27 panel approach (model-native format or revised screening criterion) before this batch can close. Sprint 40 (P1 Schelling + P3 Gray-Scott) may proceed autonomously once a human has reviewed these findings, or immediately if the orchestrator is configured to proceed past GO-LIMITED.
