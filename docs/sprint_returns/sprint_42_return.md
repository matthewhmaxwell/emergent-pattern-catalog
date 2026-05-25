# Sprint 42 Return Summary

**Sprint goal:** Run v1.2 Phase-2a panel against P1 Schelling and P3 Gray-Scott.

**Status: P1 PARTIAL (TNR=0.593, Cohen's d=1.298). P3 PAUSED (C-lattice_2d_continuous-substrate-undercount: 0 Class B mates, threshold 3). AT-DEPTH count unchanged: 6/19.**

---

## Pre-flight verification

- Base HEAD: `fd6167f` (Sprint 41 — P22 irreversibility prerequisite, v0.41.0). ✓
- `PYTHONPATH=. python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** ✓ (unchanged — no registry changes in this sprint).
- `pytest tests/ -m "not slow" --co -q`: **599 tests collected** (pre-flight). ✓

---

## Part A — P1 panel

### Class B pre-check

`class_b_for_pattern("P1")` returns:
```python
{
  "catalog_mates": ["P11_lotka_volterra", "P12_rps", "P13_greenberg_hastings",
                    "P14_btw_sandpile", "P15_gol", "P22_sir_epidemic", "P27_nowak_may"],
  "synthetic_supplements": [],
  "substrate_type": "lattice_2d"
}
```
7 lattice_2d mates ≥ 3 → no undercount condition. Panel proceeds.

### P1 invariance flags

`detector_invariance.py` for P1: `(permutation_invariant=False, time_shuffle_invariant=False)`. Both `permutation_shuffled` and `time_shuffled` Class A substrates evaluated (not skipped per v1.2 spec). Brief anticipated both would be TN; empirically `time_shuffled` was FP (see below).

### Files added / modified

**New:**
- `epc/phase2a/failed_regimes/p1_schelling.py` — 10 sub-threshold Schelling regimes: threshold ∈ linspace(0.05, 0.25, 10), grid_size=32, density=0.9, n_steps=80, seeds 100–109.

**Modified:**
- `analysis/run_phase2a_panel.py` — added `build_p1_positives()`, `make_p1_detector_fn()`, `run_p1()`; updated docstring + `main()` to handle "p1".

**Generated:**
- `analysis/outputs/p1_phase2a_panel.json` — panel output (panel_version=1.2).

### Panel result

Command: `PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p1`

| Class | TNR | n_eval | Notes |
|---|---|---|---|
| Synthetic (Class A) | **0.800** | 10 | FPs: `time_shuffled` (confirmation 0.700), `linear_gradient` (confirmation 0.700) |
| Catalog (substrate-typed: lattice_2d) | **0.571** | 7 | FPs: P11_LV (confirmation 0.700), P15_GoL (confirmation 0.700), P12_RPS (screening 0.600) |
| Failed-regime (Class C: threshold ∈ [0.05, 0.25]) | **0.400** | 10 | FPs: threshold=0.050, 0.161, 0.183, 0.206, 0.228, 0.250 |
| **Overall** | **0.593** | 27 | **PARTIAL** |

- Canonical positive: 5/5 seeds at CONFIRMATION (confidence=0.700).
- Cohen's d: **1.298** (positive mean 0.700 vs. negative mean 0.333).
- **Verdict: PARTIAL** (TNR 0.593 < 0.95; d=1.298 ≥ 0.5).

### FP root-cause analysis

**Class A false positives:**

1. `time_shuffled` — time-shuffles the canonical positive's frames. Schelling segregation is a fast-converging process: high Moran's I is established in the first 20–50 steps of a 200-step run. The time-shuffled final state is a random frame drawn from the canonical run; most frames in the run have high Moran's I. P1's final-state Moran's I primary metric cannot distinguish "final frame of time-shuffled canonical run" from "final frame of genuine run." Root cause: P1 is a spatial snapshot detector, not a temporal formation detector.

2. `linear_gradient` — a spatially smooth linear gradient has high spatial autocorrelation (adjacent cells have similar values) and therefore high Moran's I. No confirmation-tier gate distinguishes "gradient structure" from "clustering structure."

**Class B false positives:**

- P11_lotka_volterra: predator-prey lattice in coexistence generates persistent spatial clusters (spirals + patches) with Moran's I comparable to Schelling at confirmation tier.
- P15_gol: Game of Life random-dense positive maintains persistent localized structures (gliders, blinkers, stable oscillators) with elevated spatial autocorrelation.
- P12_rps: RPS spiral domains maintain high Moran's I; triggers at SCREENING tier (0.600).

**Class C false positives:**

- threshold=0.050: most tolerant regime; seed-100 initial placement at density=0.9 has enough incidental local correlation on a 32×32 grid to clear the CONFIRMATION gate even though agents never move.
- threshold=0.161–0.250: these regimes are close enough to the canonical threshold that some initial random clustering is present and agents move only minimally; the initial configuration's Moran's I alone clears the confirmation gate.

### Sprint 30 rule: PARTIAL → no detector/model changes

Per Sprint 30 rule, no detector or model modifications are applied. Four carry-forwards opened:

- `C-p1-time-shuffle-fp`: `time_shuffled` fires on P1 because Schelling's intermediate frames already exhibit segregation. Resolution path: add a temporal-formation secondary gate (e.g., Moran's I increase from step 0 to step T).
- `C-p1-linear-gradient-fp`: `linear_gradient` fires because Moran's I is not gradient-insensitive. Resolution path: Fourier-domain secondary metric to distinguish clustered from gradient spatial autocorrelation.
- `C-p1-class-b-lattice2d-fp`: P11/P15/P12 lattice_2d models with persistent spatial structure trigger Moran's I at confirmation tier. Resolution path: add a substrate-specific secondary gate or raise the confirmation threshold.
- `C-p1-class-c-subthreshold-fp`: sub-threshold regimes with accidental initial clustering fire at confirmation. Resolution path: choose initial conditions with seeded-uniform starts (no accidental clustering) for Class C, or widen sub-threshold range further from 0.375.

---

## Part B — P3 panel (paused)

### Class B pre-check: UNDERCOUNT CONDITION FIRES

`class_b_for_pattern("P3")` at Sprint 42 HEAD:
```python
{
  "catalog_mates": [],
  "synthetic_supplements": [],
  "substrate_type": "lattice_2d_continuous"
}
```

0 lattice_2d_continuous mates. Threshold is 3. Condition fires.

**Action:** P3 panel paused. No `analysis/outputs/p3_phase2a_panel.json` written.

**Carry-forward logged:** `C-lattice_2d_continuous-substrate-undercount` — P3 (Gray-Scott) is the only registered lattice_2d_continuous pattern. The Class B pool for this substrate type has 0 members, making the substrate-typed Class B approach infeasible for P3. Brief-recommended resolution: use lattice_2d catalog mates as a fallback for lattice_2d_continuous in `class_b_for_pattern()`. This is a spec call (chat-led), out of scope for Sprint 42.

**state.json update:** `in_flight=null`, `last_escalation` populated with sprint 42 P3 undercount + P1 PARTIAL context.

---

## Part C — REPLICATION_NOTES.md

Two sections appended to `REPLICATION_NOTES.md`:

1. **`## Phase-2a Panel Result (v1.2) — Sprint 42 (P1 Schelling, PARTIAL)`** — full per-class TNR table, FP breakdown, carry-forward list.
2. **`## Phase-2a Panel — Sprint 42 P3 Pause`** — documents the lattice_2d_continuous undercount condition, action taken, and carry-forward.

`grep -c "Phase-2a Panel Result (v1.2)" REPLICATION_NOTES.md` → **6** (was 5; incremented by 1 for P1 result).

---

## Part D — depth_gap.md update

- **P1 row `dim4_negative_sweep`**: remains PARTIAL. Updated notes column to document Sprint 42 panel result: overall TNR=0.593, per-class breakdown, 4 carry-forwards opened.
- **Aggregate findings**: Sprint 42 finding note added. AT-DEPTH count unchanged: **6 / 19**.

---

## Part E — Paper sync

- **§4.21** (`docs/paper_section4_draft.md`): added Sprint 42 section documenting P1 PARTIAL result and P3 pause. Covers FP root causes, Sprint 30 rule application, 4 carry-forwards.
- **§6** (`docs/paper_section6_draft.md`): appended Sprint 42 AT-DEPTH count note.
- **`docs/paper_CHANGELOG.md`**: Sprint 42 entry added.

---

## Part F — Carry-forward status

**NEW this sprint:**
- `C-p1-time-shuffle-fp` (Sprint 42, NEW): `time_shuffled` fires on P1 because Schelling intermediate frames already have high Moran's I.
- `C-p1-linear-gradient-fp` (Sprint 42, NEW): `linear_gradient` fires on P1 because Moran's I responds to gradient structure.
- `C-p1-class-b-lattice2d-fp` (Sprint 42, NEW): P11/P15/P12 lattice_2d models are false positives under P1 confirmation gate.
- `C-p1-class-c-subthreshold-fp` (Sprint 42, NEW): sub-threshold Schelling regimes with accidental initial clustering fire at confirmation.
- `C-lattice_2d_continuous-substrate-undercount` (Sprint 42, NEW): P3 lattice_2d_continuous has 0 Class B mates; panel paused pending spec call.

**Still open (unchanged):**
- C-p21-constant-field-trivial-sync (Sprint 35): deferred
- C-p14-class-c-borderline (Sprint 35): deferred
- C-p22-dim1-no-fig-reproduced: P22 dim1 still PARTIAL
- C-p22-dim2-no-basin-fraction: P22 dim2 still PARTIAL

---

## Pre-flight + post-flight test counts

| Stage | non-slow tests |
|---|---|
| Pre-flight (before Sprint 42) | 599 |
| Post-flight (after Sprint 42) | **599** |
| Delta | 0 (no new tests — panel-runs sprint, no new detector/model code) |

Matrix counts: **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** (unchanged — no registry changes in this sprint). ✓

---

## Files changed

| File | Change |
|---|---|
| `epc/phase2a/failed_regimes/p1_schelling.py` | NEW: 10 sub-threshold Schelling Class C regimes |
| `analysis/run_phase2a_panel.py` | Added `build_p1_positives()`, `make_p1_detector_fn()`, `run_p1()`; updated docstring + main() |
| `analysis/outputs/p1_phase2a_panel.json` | NEW: panel output, TNR=0.593, verdict=PARTIAL |
| `docs/depth_gap.md` | P1 dim4 notes updated; Sprint 42 aggregate finding added |
| `REPLICATION_NOTES.md` | Appended Phase-2a Panel Result (v1.2) Sprint 42 (P1) + P3 pause section |
| `docs/paper_section4_draft.md` | §4.21: Sprint 42 panel results section |
| `docs/paper_section6_draft.md` | Sprint 42 AT-DEPTH count note appended |
| `docs/paper_CHANGELOG.md` | Sprint 42 entry |

---

## Final HEAD and tag

- **Commit 1 (sprint body):** Sprint 42 main implementation
- **Commit 2 (return doc):** `docs/sprint_returns/sprint_42_return.md`
- **Tag:** `v0.42.0`

---

**Decision: GO-LIMITED**

P1 panel returned PARTIAL (TNR=0.593, Cohen's d=1.298) per Sprint 30 rule — 4 carry-forwards opened, no detector changes applied. P3 panel paused due to lattice_2d_continuous substrate-undercount (0 Class B mates). Recommend human review before Sprint 43 to:

1. Confirm Sprint 30 rule handling is correct for P1 (carry-forwards sufficient vs. immediate fix).
2. Decide on `class_b_for_pattern("P3")` fallback spec: use lattice_2d mates for lattice_2d_continuous, or add more lattice_2d_continuous models first?

Sprint 43 brief (lattice_2d dim4 batch 3: P12 RPS + P13 GH) does not depend on P3 resolution and can proceed independently.
