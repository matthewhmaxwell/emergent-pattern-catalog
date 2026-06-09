# Sprint 74 Return — P26 Stochastic Resonance (Milestone B Wave 2)

**Date:** 2026-06-09
**Base HEAD (sprint start):** `730c6cd` (Sprint 73 follow-up)
**Sprint goal:** Implement P26 (stochastic resonance) end-to-end with OOD-readiness (T1a/T1b). dim1 reproduction, dim2 multi-seed, dim3 methods note.
**Tag:** `v0.74.0`
**Sprint type:** Chat-led design + code-led execution. Milestone B, Wave 2.

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `730c6cd` ✓
2. **Working tree clean:** ✓

---

## Part A — Model

`epc/models/stochastic_resonance.py` — two SR implementations:

1. **BistableDoubleWell** (canonical): overdamped particle in V(x) = -ax^2/2 + bx^4/4, driven by subthreshold periodic signal A*sin(2*pi*f*t) + Gaussian noise sqrt(2D)*xi(t). Parameters: a=4.0, b=1.0 (barrier=4.0), A=1.0, f=0.005, dt=0.01, n_steps=20000, n_trials=20. Multi-trial design: runs n_trials independent simulations per noise level with alternating starting wells.

2. **ThresholdUnit** (T1b): binary threshold detector with subthreshold signal (amplitude=0.7, threshold=1.0). At intermediate noise, signal+noise crosses threshold preferentially when signal is near peak.

Both output history dicts with keys: time, x, signal, noise_level, noise_level_idx, step. Substrate: `noise_sweep_timeseries` (9th substrate type).

---

## Part B — Detector

`epc/detectors/p26_stochastic_resonance.py`:

- **Primary metric:** Coherent response |mean(x * signal)| — sensitive to both amplitude and phase. At zero noise: tiny oscillation -> small metric. At optimal noise: synchronized inter-well hopping -> large metric. At high noise: random hopping -> metric ~0.
- **Null model:** Time-shuffle x at peak noise level (199 permutations). Destroys temporal synchronization, preserving marginal distribution. Shuffled coherent response ~ 0 vs observed ~ 0.918.
- **T1a adapter:** `extract_observation_bundle()` returns 5 aligned arrays (time, x, signal, noise_level, noise_level_idx).
- **Tier criteria:** Screening (gain > 0.02), Confirmation (inverted-U shape + p < 0.05), Definitive (gain > 0.05 + decline > 0.02 + d > 1.0 + p <= 0.005 + subthreshold metadata gate).

### Key design decisions

1. **Coherent response over Pearson correlation:** correlation is ~1.0 at zero noise (phase-locked despite negligible amplitude) — no inverted-U.
2. **Coherent response over FFT SNR:** infinite SNR at zero noise (no noise floor) — degenerate.
3. **Multi-trial averaging:** 20 independent trials per noise level reduces variance by sqrt(20) ~ 4.5x.
4. **Time-shuffle over circular shift:** circular shifting preserves x's periodicity -> non-zero null. Time-shuffle fully destroys temporal structure.

---

## Part C — Tests

`tests/test_stochastic_resonance_p26_e2e.py` — 14 tests:
- `TestP26DoubleWellPositive` (6): detected, definitive tier, interior peak, inverted-U shape, gain > 0.05, deterministic from seed
- `TestP26NegativeControl` (2): suprathreshold signal not detected, zero-noise-only not detected
- `TestT1aObservationBundle` (3): bundle keys, values match history, detector uses bundle (manual history dicts)
- `TestP26ThresholdUnitPositive` (3): detected, definitive tier, interior peak

`tests/test_cross_model.py` — 2 T1b tests added:
- `test_p26_on_threshold_unit_detected`: threshold unit reaches confirmation+
- `test_p26_both_models_definitive`: both SR implementations reach definitive

---

## Part D — dim1 reproduction

`analysis/reproductions/p26_collins.py` — Gammaitoni (1998) / Collins (1995):
- Peak noise: D = 1.5 (interior, idx=7 of 15)
- Peak coherent response: 0.918
- Zero-noise coherent response: 0.063
- **Gain over zero: 0.855** (tolerance > 0.05) — **PASS**
- **Decline after peak: 0.811** (tolerance > 0.02) — **PASS**
- **Interior argmax:** idx=7 of 15 — **PASS**
- Detector tier: **DEFINITIVE** (confidence=0.90) — **PASS**

Output: `analysis/outputs/p26_collins_reproduction.json`

---

## Part E — dim2 + dim3

### dim2: 20-seed multi-seed campaign

`analysis/reproductions/p26_multiseed.py`:
- Peak coherent response: 0.897 ± 0.049 (CV=5.4%)
- Gain over zero: 0.833 ± 0.049 (CV=5.8%)
- Peak noise location: D = 1.27 ± 0.25
- All 20 seeds: **DEFINITIVE** (fraction=1.00)

Output: `analysis/outputs/p26_multiseed.json`

### dim3: methods note

`docs/methods_notes/p26_methods.md` — covers:
- Bistable double-well dynamics and ThresholdUnit variant
- Coherent response metric and why alternatives fail
- Time-shuffle null model and why alternatives fail
- Multi-trial averaging design
- Tier criteria
- Four key design decisions
- Reproduction results (dim1 + dim2 + T1b)

---

## Part F — Documentation updates

- **`docs/observation_schema.md`:** +1 noise-sweep-timeseries bundle schema (P26)
- **`docs/depth_gap.md`:** +P26 row (dim1-3 PASS, dim4 pending). Implemented count 23->24. Gap count 1->2 (P12 + P26).
- **`docs/paper_section4_draft.md`:** +section 4.26 P26 stochastic resonance
- **`docs/paper_section6_draft.md`:** +Sprint 74 entry
- **`docs/paper_CHANGELOG.md`:** +Sprint 74 entry
- **`REPLICATION_NOTES.md`:** +P26 section (Sprint 74)
- **`epc/orchestration.py`:** +bistable_double_well model (noise_sweep_timeseries), +P26 detector. 25 models x 24 detectors, 106 compatible pairs.
- **`tests/test_orchestration.py`:** count updates (25 models, 24 detectors, 106 pairs, 600 cells)
- **`tests/test_transfer_matrix_counts.py`:** count updates

---

## Post-flight checks

- `pytest tests/test_orchestration.py tests/test_transfer_matrix_counts.py tests/test_stochastic_resonance_p26_e2e.py tests/test_cross_model.py -m "not slow"`: **111 passed** ✓
- P26 in DETECTOR_REGISTRY ✓
- bistable_double_well in MODEL_REGISTRY ✓
- noise_sweep_timeseries substrate registered ✓
- T1a observation bundle contract tested (3 tests) ✓
- T1b cross-model generalization (threshold unit -> DEFINITIVE) ✓
- dim1 reproduction all 4 tolerance checks PASS ✓
- dim2 multiseed 20/20 DEFINITIVE, CV=5.4% ✓
- dim3 methods note authored ✓
- depth_gap P26 row added ✓
- paper section 4 + section 6 updated ✓

---

## Carry-forwards

- **P26 dim4:** Phase-2a panel pending (future sprint)
- **C-p7-time-shuffled-fp:** `time_shuffled` FP at screening (from Sprint 66, low priority)
- **C-p19-bias-zero-chance-alignment:** 1/5 bias_zero Class C FP at confirmation (from Sprint 70, 4% rate)
- **P12 dim1:** documented finite-size measurement limitation (unchanged, accepted)

---

## Files added/modified

**New (6):**
- `epc/models/stochastic_resonance.py` — BistableDoubleWell + ThresholdUnit
- `epc/detectors/p26_stochastic_resonance.py` — coherent response + time-shuffle null
- `tests/test_stochastic_resonance_p26_e2e.py` — 14 tests
- `analysis/reproductions/p26_collins.py` — dim1 reproduction
- `analysis/reproductions/p26_multiseed.py` — dim2 multi-seed
- `docs/methods_notes/p26_methods.md` — dim3 methods note

**Modified (10):**
- `epc/orchestration.py` — +1 model, +1 detector, +1 substrate
- `tests/test_cross_model.py` — +2 T1b tests
- `tests/test_orchestration.py` — count updates
- `tests/test_transfer_matrix_counts.py` — count updates
- `docs/observation_schema.md` — +P26 bundle
- `docs/depth_gap.md` — +P26 row
- `docs/paper_section4_draft.md` — +section 4.26
- `docs/paper_section6_draft.md` — +Sprint 74 entry
- `docs/paper_CHANGELOG.md` — +Sprint 74 entry
- `REPLICATION_NOTES.md` — +P26 section

**Output (2):**
- `analysis/outputs/p26_collins_reproduction.json` — dim1 results
- `analysis/outputs/p26_multiseed.json` — dim2 results

---

**Decision: GO**
