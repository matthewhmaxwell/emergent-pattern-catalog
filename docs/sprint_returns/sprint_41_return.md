# Sprint 41 Return Summary

**Sprint goal:** Add a content-level **irreversibility prerequisite** to
`P22CascadeDetector`, grounded in Datta & Acharyya (2021) SIR characterization.
SIR's S→I→R transitions are irreversible — once Recovered, cells stay
Recovered. LV and RPS both feature reversible/cyclic state transitions. Adding
this prerequisite is the SIR-specific analog of the P11 conservation guard and
the P27 observable guard (paper §3.5).

**Status: P22 PASS (TNR=1.000, Cohen's d=+∞). C-p22-class-b-cascade-overlap
CLOSED. P22 dim4 → PASS; dims 1–3 remain PARTIAL, so AT-DEPTH not yet reached.
AT-DEPTH count stays at 6/19.**

---

## Pre-flight verification

- Base HEAD: `e50afaa` (Sprint 40 follow-up — GO-LIMITED → GO amendment). ✓
  (Brief cited `9397193` as base; `e50afaa` is a clean follow-up commit that
  amended the sprint_40_return.md decision line. No substantive code change.)
- `PYTHONPATH=. python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 /
  27 / 19 / 361 / 77 / 284** ✓ (unchanged — no registry changes in this sprint).
- `pytest tests/ -m "not slow" --co -q`: **597 tests collected** (pre-flight). ✓

---

## Part A — REPLICATION_NOTES.md SIR section (literature confirmation)

**File:** `REPLICATION_NOTES.md`, §SIR Epidemic CA Replication Notes (Sprint 7–8)

The SIR section contains the following relevant characterizations:

1. **Implementation table (line ~1308):**
   > "Immunity | Permanent (R never → S) | Permanent ✅"

   This is the explicit documentation that R is absorbing — no R→S transition
   is ever generated. This directly grounds the irreversibility prereq.

2. **P1/SIR comparison narrative (line ~1719–1724):**
   > "SIR's infected cells recover irreversibly; once the wavefront has passed
   > a cell, that cell stays recovered forever and the final state is near-uniform."

   This sentence was written in the P1/SIR Sprint 10 analysis to explain why
   SIR's final Moran's I is near zero. It doubles as the irreversibility
   characterization required by the prereq.

The irreversibility property is **documented implicitly via these SIR characterization sentences** — it is not explicitly named as "the irreversibility prerequisite" before Sprint 41, which is acceptable per the brief. The canonical paper (Datta & Acharyya 2021) encodes this in the recovery rule "P(I→R) = q per timestep" and the immunity rule "Permanent" — there is no S→R reversal in the model.

**Literature-grounding for the three distinguishing papers:**

- **Datta, A. & Acharyya, M. (2021/2022).** "Modelling the Spread of an Epidemic
  in Presence of Vaccination using Cellular Automata." arXiv:2104.10456, Int.
  J. Mod. Phys. C 33, 2250094. — Canonical SIR reference for this project.
  The defining feature: S→I→R transitions are strictly unidirectional. Recovery
  (I→R) is irreversible; once recovered, a cell is permanently immune. The state
  encoding is S=0, I=1, R=2 (confirmed in `epc/models/sir_epidemic.py`).

- **Mobilia, M., Georgiev, I.T. & Täuber, U.C. (2007).** "Phase Transitions and
  Spatio-Temporal Fluctuations in Stochastic Lattice Lotka-Volterra Models." J.
  Stat. Phys. 128, 447–483. — Canonical LV reference for P11. Cells cycle
  through EMPTY(0), PREY(1), PREDATOR(2). Predator death (2→0) and predation of
  prey (1→0) are backward transitions under the SIR irreversibility convention.
  The LV model has no absorbing "recovered" state — every species state is
  revisitable.

- **Reichenbach, T., Mobilia, M. & Frey, E. (2007).** "Mobility promotes and
  jeopardizes biodiversity in rock-paper-scissors games." *Nature* 448, 1046–1049.
  — Canonical RPS reference for P12. Three species {A=1, B=2, C=3} with EMPTY=0.
  Selection reactions convert species→EMPTY (2→0, 1→0, 3→0). No absorbing state
  exists in the coexistence regime — all states are freely revisitable via
  cyclic dominance.

---

## Part B — Irreversibility prerequisite implementation

**File:** `epc/detectors/p22_information_cascade.py`

**Implementation:**

Added two methods to `P22CascadeDetector`:

1. **`detect()` override** — checks irreversibility before the full pipeline:
   ```python
   found_reversible, warning = self._check_irreversibility_prereq(state_history)
   if found_reversible:
       return DetectorResult(
           pattern_id=self.pattern_id, detected=False,
           tier=DetectionTier.SCREENING, confidence=0.0,
           primary_metric={}, secondary_metrics={}, effect_size={},
           null_p_value=1.0, null_type=NullType.SHUFFLE,
           warnings=[warning],
       )
   return super().detect(state_history, model_metadata, timescale)
   ```

2. **`_check_irreversibility_prereq()`** — vectorized numpy scan:
   ```python
   for t in range(1, len(state_history)):
       if np.any(curr_grid < prev_grid):
           return True, "P22 requires irreversible S→I→R transitions per ..."
   return False, ""
   ```

**State encoding confirmation:** `SIREpidemicModel.SUSCEPTIBLE=0, INFECTED=1, RECOVERED=2` — confirmed from `epc/models/sir_epidemic.py` lines 137–140. Allowed transitions are non-decreasing: {(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)}. Forbidden = any (old, new) where new < old.

**Design choice:** Implemented as a `detect()` override (not a `_validate_prerequisites` return) because the base class pipeline does not support short-circuit returns from `_validate_prerequisites`. This follows the pattern of the P27 detector, which also short-circuits at the top of its detect function.

---

## Part C — Regression tests

**File:** `tests/test_sir_p22_e2e.py` (2 tests added to `TestP22IrreversibilityPrereq`)

1. `test_p22_short_circuits_on_lv_substrate` — runs `LotkaVolterraLattice` (50×50,
   predation_rate=2.0, n_steps=100), feeds history to P22 detector.
   Asserts: `result.detected is False`, `result.tier == SCREENING`,
   `result.confidence == 0.0`, `"irreversib" in warning (case-insensitive)`.
   **PASS** (predator death transitions 2→0 trigger the guard).

2. `test_p22_still_fires_on_sir_canonical` — runs `SIREpidemicModel` (80×80,
   p=0.20, q=0.3, single_seed, n_steps=400).
   Asserts: `result.detected is True`, `result.tier >= DEFINITIVE`.
   **PASS** — no regression on native domain.

**Pre-existing test update:** `TestP22CrossDetection::test_nowak_may_rejected_by_p22`
was asserting `moran_i_time < 0.1` but the Sprint 41 guard now fires on Nowak-May
(cooperator=1 → defector=0 transitions are backward), returning `primary_metric={}`.
The default `get("moran_i_time", 1.0)` returned 1.0, which failed the assertion.
Updated the test to assert `"irreversib" in warning` instead — the new assertion
is physically correct (NM has genuinely reversible transitions) and more precisely
captures the reason for rejection. No test deleted; behavior is correct.

**Total test count:** 599 non-slow tests (597 pre-flight + 2 new). ✓

---

## Part D — P22 panel re-run (v1.2, Sprint 41)

```
PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p22
```

### Before (Sprint 40) vs After (Sprint 41)

| Class | Sprint 40 TNR | Sprint 41 TNR | Change |
|---|---|---|---|
| Synthetic (Class A) | 0.900 | **1.000** | `time_shuffled` now correctly rejected |
| Catalog (Class B, lattice_2d) | 0.714 | **1.000** | LV + RPS now correctly rejected |
| Failed-regime (Class C) | 1.000 | **1.000** | unchanged |
| **Overall** | **0.889** | **1.000** | PASS |
| Cohen's d | 2.981 | **+∞** | all negatives at 0.0 |
| Verdict | PARTIAL | **PASS** | |

### Per-substrate guard behavior (Class B)

| Substrate | Guard trigger | Why |
|---|---|---|
| P11_lotka_volterra | 2→0 (predator death, first step) | LV PREDATOR dies → EMPTY |
| P12_rps | 3→0 or 2→0 (species → empty) | RPS cyclic dominance |
| P13_greenberg_hastings | guard fires (n−1→0 reset) | GH clock: state n−1 resets to 0 |
| P14_btw_sandpile | guard skips (no `grid` key) | sandpile uses avalanche format |
| P15_gol | 1→0 (cell dies) | GoL binarized: living→dead |
| P1_schelling | 1→0 (binarized) | binarized grid has 1→0 |
| P27_nowak_may | 1→0 (cooperator → defector) | NM binary state reversible |

Note: The panel runner binarizes categorical grids (LV, RPS) via `(grid != 0)`
before passing to detectors. Even with binarization, LV/RPS have 1→0 transitions
(species → empty), triggering the guard.

**Panel JSON:** `analysis/outputs/p22_phase2a_panel.json`

---

## Part E — depth_gap.md update

**P22 row change:** `dim4_negative_sweep`: PARTIAL → **PASS**. Grade remains **GAP**
(dims 1–3 still PARTIAL; AT-DEPTH requires all 4 PASS).

**Aggregate findings updated:**
- AT-DEPTH count: **6 / 19** (unchanged)
- Sprint 41 finding note added: P22 dim4 PASS; C-p22-class-b-cascade-overlap CLOSED

---

## Part F — Paper sync

- **§3.5** (`docs/paper_section3_draft.md`): appended P22 irreversibility-prerequisite
  paragraph alongside existing P11 conservation guard + P27 observable prereq.
  Cites Datta-Acharyya (2021), Mobilia-Georgiev-Täuber (2007), Reichenbach (2007).

- **§4.10** (`docs/paper_section4_draft.md`): appended Sprint 41 panel-result paragraph
  after Sprint 40 paragraph. Documents TNR=1.000, Cohen's d=+∞, PASS verdict;
  C-p22-class-b-cascade-overlap CLOSED.

- **§6.11** (`docs/paper_section6_draft.md`): updated aggregate narrative to note P22
  dim4 PASS; AT-DEPTH count stays 6/19 (dims 1–3 still PARTIAL).

- **`docs/paper_CHANGELOG.md`**: Sprint 41 entry added.

---

## Part G — Carry-forward status

**CLOSED this sprint:**
- `C-p22-class-b-cascade-overlap` (Sprint 40, NEW): **CLOSED**. Literature-anchored
  irreversibility prereq addresses all Class B overlap with LV/RPS. LV (predator
  death = 2→0) and RPS (cyclic species elimination = any→0) both have backward
  transitions under the SIR monotone-state convention. TNR raised from 0.714 →
  1.000 for Class B.

**Still open (unchanged):**
- C-p21-constant-field-trivial-sync (Sprint 35): deferred
- C-p14-class-c-borderline (Sprint 35): deferred
- C-p22-dim1-no-fig-reproduced: P22 dim1 still PARTIAL (no specific Fig/table
  from Datta-Acharyya 2021 reproduced with named tolerance)
- C-p22-dim2-no-basin-fraction: P22 dim2 still PARTIAL

---

## Pre-flight + post-flight test counts

| Stage | non-slow tests |
|---|---|
| Pre-flight (before Sprint 41) | 597 |
| Post-flight (after Sprint 41) | **599** |
| Delta | +2 (2 new Sprint 41 tests) |

Matrix counts: **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** (unchanged — no
registry changes in this sprint). ✓

---

## Files changed

| File | Change |
|---|---|
| `epc/detectors/p22_information_cascade.py` | Added `detect()` override + `_check_irreversibility_prereq()` |
| `tests/test_sir_p22_e2e.py` | 2 new tests + 1 updated assertion |
| `analysis/outputs/p22_phase2a_panel.json` | Re-run; TNR=1.000, verdict=PASS |
| `docs/depth_gap.md` | P22 dim4 PARTIAL→PASS; aggregate notes updated |
| `REPLICATION_NOTES.md` | Appended Phase-2a Panel Result (v1.2) Sprint 41 section |
| `docs/paper_section3_draft.md` | §3.5: P22 irreversibility prereq paragraph |
| `docs/paper_section4_draft.md` | §4.10: Sprint 41 panel-result paragraph |
| `docs/paper_section6_draft.md` | §6.11: AT-DEPTH count note updated |
| `docs/paper_CHANGELOG.md` | Sprint 41 entry |

---

## Final HEAD and tag

- **Commit 1 (sprint body):** Sprint 41 main implementation
- **Commit 2 (return doc):** `docs/sprint_returns/sprint_41_return.md`
- **Tag:** `v0.41.0`

---

**Decision: GO**
