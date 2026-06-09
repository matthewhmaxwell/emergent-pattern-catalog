# Sprint 76 Return — P23 Anti-coordination / Emergent Load Balancing (Milestone B Wave 2)

**Date:** 2026-06-09
**Base HEAD (sprint start):** `a49dff1` (Sprint 75 follow-up)
**Sprint goal:** Implement P23 — anti-coordination / emergent load balancing end-to-end with OOD-readiness (T1a/T1b).
**Tag:** `v0.76.0`
**Sprint type:** Chat-led design + code-led execution. Milestone B, Wave 2.

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `a49dff1` ✓
2. **Working tree clean:** ✓

---

## Part A — Model: `epc/models/minority_game.py`

Three classes implemented:

1. **MinorityGame** (Challet & Zhang 1997): N agents (odd N), S strategies per agent, each mapping m-bit history to binary choice. Minority side wins; strategy scores updated via virtual points. Control parameter α = 2^m / N. Deterministic from seed.

2. **ElFarolBar** (Arthur 1994): N agents, capacity C, K linear predictors of next attendance from recent history. Agents attend if best predictor forecasts attendance < capacity. Predictor scores decay with absolute error. Deterministic from seed.

3. **RandomChoiceBaseline**: Negative control — i.i.d. Binomial(N, 0.5) attendance. No strategy adaptation.

---

## Part B — Detector: `epc/detectors/p23_anticoordination.py`

**Primary metrics:**
- Scaled variance σ²/N of attendance (post burn-in) vs random-choice baseline p̂(1−p̂)
- Lag-1 autocorrelation of attendance time series

**Null model:** Random-choice surrogate — i.i.d. Binomial(N, p̂) attendance series. Tests whether reduced variance and anti-persistence require strategic adaptation.

**Tier criteria:**
- Screening: Mean attendance within 20% of capacity
- Confirmation: σ²/N < baseline OR negative autocorrelation (p < 0.01 vs surrogate)
- Definitive: BOTH variance below baseline AND negative autocorrelation (p < 0.01)

**T1a adapter:** `extract_observation_bundle(history)` reads attendance/choice time series bundle (round, attendance, n_agents, capacity). Documented in `docs/observation_schema.md`.

**P23 registered in `DETECTOR_REGISTRY`** with `choice_timeseries` substrate.

---

## Part C — Tests

### `tests/test_minority_game_p23_e2e.py` (11 tests)
- Determinism from seed ✓
- DEFINITIVE on efficient regime (m=6, α≈0.63, N=101) ✓
- Variance below random baseline ✓
- Negative autocorrelation ✓
- Detector reads via T1a adapter ✓
- Confidence ≥ 0.75 at definitive ✓
- Random-agent negative control: NOT detected ✓
- Random variance near 0.25 baseline ✓
- Symmetric phase (m=1): not definitive ✓
- El Farol cross-model: detected at confirmation+ ✓
- El Farol deterministic from seed ✓

### `tests/test_cross_model.py` — T1b (2 tests)
- P23 on El Farol detected at confirmation+ ✓
- Both MG (definitive) and El Farol (detected) pass ✓

### `tests/test_orchestration.py` — Sprint 76 registrations (6 tests)
- minority_game registered (choice_timeseries, attendance, P23) ✓
- el_farol registered (choice_timeseries, attendance, P23) ✓
- P23 registered (choice_timeseries, attendance) ✓
- minority_game × P23 compatible ✓
- el_farol × P23 compatible ✓
- P23 only pairs with choice_timeseries models ✓

### Transfer matrix counts updated
- 27 models, 25 detectors, 675 cells, 108 compatible pairs, 10 substrate types

---

## Part D — dim1: Savit curve reproduction

**Anchor:** Savit, Manuca & Riolo (1999). PRL, 82(10), 2203–2206.

**Parameters:** N=101, m=1..11, 10 seeds per point, 3000 rounds, burn-in 600.

**Results:**
- Interior minimum at α ≈ 0.32 (m=5): σ²/N = 0.077
- Random baseline: 0.25
- Symmetric phase (m=1): σ²/N ≈ 1.45
- All 3 tolerance checks **PASS**

Output: `analysis/outputs/p23_savit_reproduction.json`

---

## Part E — dim2 multi-seed + dim3 methods note

### dim2: 25 seeds at α ≈ 0.63 (m=6, N=101)
- σ²/N: mean=0.075, std=0.006, CV=8.7%
- All 25 seeds detected (23 confirmation, 2 definitive)
- Output: `analysis/outputs/p23_multiseed.json`

### dim3: `docs/methods_notes/p23_methods.md`
- MG and El Farol dynamics
- Random-choice surrogate null design
- Dual-metric (variance + autocorrelation) detection
- Savit curve anchor
- Limitations: autocorrelation not consistently negative at all efficient-phase seeds

---

## Part F — Documentation updates

- **`docs/depth_gap.md`:** P23 row added (dim1–3 PASS, dim4 pending, GAP). Pattern count 24→25. Gap count 1→2.
- **`docs/observation_schema.md`:** Attendance/choice time series bundle documented.
- **`docs/paper_section4_draft.md`:** §4.23 Sprint 76 P23 section added.
- **`docs/paper_section6_draft.md`:** Sprint 76 entry added.
- **`docs/paper_CHANGELOG.md`:** Sprint 76 entry added.
- **`REPLICATION_NOTES.md`:** Sprint 76 P23 section added.
- **`epc/orchestration.py`:** 10th substrate (choice_timeseries), 2 models, 1 detector. Docstring updated (27 models × 25 detectors).

---

## Post-flight checks

- `pytest tests/test_orchestration.py tests/test_transfer_matrix_counts.py tests/test_minority_game_p23_e2e.py -m "not slow"`: **119 passed** ✓
- `pytest tests/test_cross_model.py tests/test_cross_detection_matrix.py -m "not slow"`: **24 passed** ✓ (cross-model: 12 tests incl 2 new P23; cross-detection matrix: 24 tests)
- P23 in DETECTOR_REGISTRY ✓
- minority_game and el_farol in MODEL_REGISTRY ✓
- T1a adapter reads observation bundle ✓
- T1b El Farol cross-model detected at confirmation ✓
- `analysis/outputs/p23_savit_reproduction.json` exists; passes_tolerance=True ✓
- `analysis/outputs/p23_multiseed.json` exists; all 25 seeds detected ✓
- depth_gap P23 row added, pattern count 24→25 ✓
- paper section 4 + section 6 + CHANGELOG updated ✓
- REPLICATION_NOTES Sprint 76 section added ✓

---

## Carry-forwards

- **C-p7-time-shuffled-fp:** `time_shuffled` FP at screening (from Sprint 66, low priority)
- **C-p19-bias-zero-chance-alignment:** 1/5 bias_zero Class C FP at confirmation (from Sprint 70, 4% rate)
- **P12 dim1:** documented finite-size measurement limitation (unchanged, accepted)
- **P23 dim4:** Phase-2a panel pending (Sprint 77)

---

## Files added/modified

**New (6):**
- `epc/models/minority_game.py` — MinorityGame, ElFarolBar, RandomChoiceBaseline
- `epc/detectors/p23_anticoordination.py` — P23 detector with T1a adapter
- `tests/test_minority_game_p23_e2e.py` — 11 e2e tests
- `analysis/reproductions/p23_savit.py` — dim1 Savit curve reproduction
- `analysis/reproductions/p23_multiseed.py` — dim2 multi-seed campaign
- `docs/methods_notes/p23_methods.md` — dim3 methods note

**Modified (10):**
- `epc/orchestration.py` — +2 models, +1 detector, +1 substrate type
- `tests/test_orchestration.py` — count updates (25→27 models, 24→25 detectors, 106→108 pairs, 9→10 substrates) + 8 Sprint 76 tests
- `tests/test_transfer_matrix_counts.py` — EXPECTED counts updated
- `tests/test_cross_model.py` — +TestP23CrossModel (2 tests)
- `docs/depth_gap.md` — P23 row, counts updated
- `docs/observation_schema.md` — attendance/choice time series bundle
- `docs/paper_section4_draft.md` — §4.23 P23 section
- `docs/paper_section6_draft.md` — Sprint 76 entry
- `docs/paper_CHANGELOG.md` — Sprint 76 entry
- `REPLICATION_NOTES.md` — Sprint 76 P23 section

**Output (2):**
- `analysis/outputs/p23_savit_reproduction.json` — dim1 results
- `analysis/outputs/p23_multiseed.json` — dim2 results

---

**Decision: GO**
