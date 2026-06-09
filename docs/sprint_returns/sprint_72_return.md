# Sprint 72 Return — P24 Homeostatic Regulation (Milestone B Wave 2)

**Date:** 2026-06-09
**Base HEAD (sprint start):** `25b3cfb` (Sprint 71 operator release)
**Sprint goal:** Implement P24 homeostatic regulation end-to-end: model + detector + tests + registry + dim1 + dim2 + dim3, with OOD-readiness (T1a observation contract, T1b cross-model test).
**Tag:** `v0.72.0`
**Sprint type:** Chat-led design + code-led execution. Milestone B, Wave 2 (first pattern).

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `25b3cfb` ✓
2. **Working tree clean:** ✓
3. **Orchestration tests pass:** 85 passed ✓
4. **Transfer matrix counts match:** ✓

---

## Part A — Model: `epc/models/homeostasis.py`

Two controller variants implemented:
- **`ProportionalHomeostat`**: dx/dt = −gain × (x − setpoint) + perturbation(t) + noise. Canonical negative-feedback regulation. Steady-state deviation = perturbation/gain.
- **`IntegralHomeostat`**: PID with P=0, I=gain, D=0. Accumulates error integral; achieves zero steady-state error for step perturbations. Independent implementation for T1b.

Both return history dicts with keys: `time`, `x`, `setpoint`, `perturbation`, `deviation`, `step`.

Perturbation schedule parameterized by onset, optional offset, and amplitude. `get_metadata()` returns model class `homeostatic_regulation` + `has_active_feedback` flag.

---

## Part B — Detector: `epc/detectors/p24_homeostasis.py`

**T1a observation contract:** `extract_observation_bundle(history)` → aligned arrays `{time, x, setpoint, perturbation}`. Detector reads exclusively through this adapter.

**Null model:** Surrogate uncontrolled trajectories — integrate dx = perturbation × dt (no feedback) with noise estimated from pre-onset residuals. Compares deviation integral of observed vs surrogates.

**Tier criteria:**
- Screening: deviation growth ratio ≤ 2.0 (bounded deviation, not linear drift)
- Confirmation: deviation integral < surrogate null at p < 0.01 AND ratio < 0.5
- Definitive: ratio < 0.3 AND Cohen's d > 1.0 AND p ≤ 0.005 AND metadata confirms active feedback

**Key ADR:** Surrogate uncontrolled null over time-shuffle null (ADR 57). A time-shuffle null is ineffective because the regulated variable's marginal distribution is narrow; the surrogate explicitly tests "what if there were no feedback?"

---

## Part C — Tests

**`tests/test_homeostasis_p24_e2e.py`** — 13 tests:
- `TestP24ProportionalPositive`: DEFINITIVE detection, deviation ratio < 0.3, bounded growth, deterministic from seed (5 tests)
- `TestP24NegativeControl`: gain=0 not detected, x drifts (2 tests)
- `TestT1aObservationBundle`: bundle keys, values match, detector uses bundle not model (3 tests)
- `TestP24PulsePerturbation`: pulse detected (1 test)
- `TestP24MetadataInteraction`: no-metadata still detects, false feedback blocks definitive (2 tests)

**`tests/test_cross_model.py`** — +2 T1b tests:
- `test_p24_on_integral_homeostat_detected`: integral controller → confirmation+ ✓
- `test_p24_integral_vs_proportional_both_definitive`: both variants DEFINITIVE ✓

---

## Part D — dim1 reproduction

`analysis/reproductions/p24_homeostasis.py` → `analysis/outputs/p24_homeostasis_reproduction.json`

- Controlled deviation integral: 149.75
- Uncontrolled deviation integral: 56,175.03
- **Deviation ratio: 0.0027** (tolerance < 0.30) — **PASS**
- Detector tier: **DEFINITIVE** (confidence=0.90) — **PASS**
- `passes_tolerance: true`

---

## Part E — dim2 + dim3

### dim2: `analysis/outputs/p24_multiseed.json`

20-seed campaign (gain=5.0, noise_std=0.5, sustained perturbation):
- Deviation integral: 149.59 ± 1.20 (CV = 0.8%)
- Deviation ratio: 0.0027 ± 0.00003 (CV = 1.2%)
- Steady-state deviation: 0.999 ± 0.015
- All 20 seeds: **DEFINITIVE**

### dim3: `docs/methods_notes/p24_methods.md`

Covers: proportional/integral controller dynamics, perturbation schedule, Euler integration, surrogate null model design, tier criteria, T1a contract, T1b cross-model, reproduction results.

---

## Part F — Registry + docs

- **`epc/orchestration.py`:** +1 model (`proportional_homeostat`, `scalar_timeseries`), +1 detector (`P24`). New substrate: `scalar_timeseries`. Counts: 24 models × 23 detectors, 105 compatible pairs.
- **`docs/depth_gap.md`:** +P24 row (dim1–3 PASS, dim4 pending). Implemented count → 23, gap count → 2 (P12 + P24).
- **`docs/observation_schema.md`:** Created — T1a scalar-regulated-variable bundle schema documented.
- **`docs/paper_section4_draft.md`:** §4.28 stub added.
- **`docs/paper_CHANGELOG.md`:** Sprint 72 entry.
- **`tests/test_orchestration.py`:** Count updates (24 models, 23 detectors, 105 pairs, 552 cells, 447 mismatches).
- **`tests/test_transfer_matrix_counts.py`:** Count updates.

---

## Post-flight checks

- `pytest tests/test_orchestration.py tests/test_transfer_matrix_counts.py tests/test_homeostasis_p24_e2e.py tests/test_cross_model.py -m "not slow"`: **108 passed** ✓
- Model + detector exist ✓
- P24 in DETECTOR_REGISTRY ✓
- Detector reads via T1a observation-bundle adapter ✓
- gain=0 negative control: NOT detected ✓
- T1b cross-model (integral controller): DEFINITIVE ✓
- `p24_homeostasis_reproduction.json`: deviation ratio 0.0027 < 0.30 ✓
- `p24_multiseed.json`: 20 seeds, all DEFINITIVE ✓
- `p24_methods.md` written ✓
- `observation_schema.md` created ✓
- depth_gap P24 row added, implemented count → 23 ✓

---

## Carry-forwards

- **C-p24-dim4-pending:** P24 dim4 Phase-2a panel deferred to Sprint 73 (per brief). Class B: passive decay / P25 equifinality must be rejected.
- **C-p7-time-shuffled-fp:** `time_shuffled` FP at screening (low priority, from Wave 1)
- **C-p19-bias-zero-chance-alignment:** 1/5 bias_zero Class C FP at confirmation (4% rate)
- **P12 dim1:** documented finite-size measurement limitation (unchanged)

---

## Files added/modified

**New (9):**
- `epc/models/homeostasis.py` — ProportionalHomeostat + IntegralHomeostat
- `epc/detectors/p24_homeostasis.py` — P24 detector + T1a adapter
- `tests/test_homeostasis_p24_e2e.py` — 13 e2e tests
- `analysis/reproductions/p24_homeostasis.py` — dim1 reproduction script
- `analysis/reproductions/p24_multiseed.py` — dim2 multi-seed script
- `analysis/outputs/p24_homeostasis_reproduction.json` — dim1 output
- `analysis/outputs/p24_multiseed.json` — dim2 output
- `docs/methods_notes/p24_methods.md` — dim3 methods note
- `docs/observation_schema.md` — T1a observation contract schema

**Modified (6):**
- `epc/orchestration.py` — +1 model, +1 detector, +1 substrate
- `tests/test_orchestration.py` — count updates
- `tests/test_transfer_matrix_counts.py` — count updates
- `tests/test_cross_model.py` — +2 T1b cross-model tests
- `docs/depth_gap.md` — +P24 row, count updates
- `docs/paper_section4_draft.md` — §4.28 stub
- `docs/paper_CHANGELOG.md` — Sprint 72 entry

---

**Decision: GO**
