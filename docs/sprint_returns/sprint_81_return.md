# Sprint 81 Return — P25 Canalized Restoration / Equifinality

**Date:** 2026-06-09
**Base HEAD (sprint start):** `bcd0592` (Sprint 80 follow-up)
**Sprint goal:** Implement P25 — Canalized restoration (equifinality) end-to-end with OOD-readiness. dim4 panel = Sprint 82.
**Tag:** `v0.81.0`
**Sprint type:** Chat-led design + code-led execution. Milestone B, Wave 3.

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `bcd0592` ✓
2. **Working tree clean:** ✓

---

## Part A — Model: `epc/models/canalization.py`

Four classes implemented:

1. **CanalizedLandscape**: Multi-dimensional gradient-flow on a combined
   linear + quartic potential: `dx/dt = -α·(x - x*) - 0.4·|x - x*|²/d²·(x - x*) + noise`.
   The quartic term creates genuine basin structure with nonlinear,
   distance-dependent relaxation — NOT trivial exponential decay.
   Parameters: n_dims=10, basin_strength=2.0, ic_spread=5.0, n_ics=20,
   n_steps=200, dt=0.05.

2. **MultiBasinGRN**: Continuous-valued gene regulatory network with
   sigmoidal activation and Hebbian weights:
   `dx_i/dt = -decay·x_i + σ(W·x + bias)`. Provides T1b cross-model.
   Bias=3.0 ensures deep dominant basin; all 20 ICs converge by step 400.

3. **DiffusiveDynamics**: Pure random walk negative control. ICs diverge;
   final-state variance ≥ IC variance. Detector correctly rejects.

4. **TrivialCollapse**: Instant constant map negative control. State set
   to collapse point at step 0 regardless of IC. Detector rejects via
   trivial-collapse gate (relaxation_time ≤ 1).

---

## Part B — Detector: `epc/detectors/p25_equifinality.py`

**Primary metric:** Convergence variance ratio = Var(finals) / Var(ICs).
Canalized systems: ratio ≪ 1.

**Null model:** IC-distribution surrogate. Generate surrogate "final states"
by sampling from the IC distribution (multivariate Gaussian with IC mean
and covariance). If NOT canalized, ratio ≈ 1.0. P-value = fraction of
surrogates ≤ observed.

**Tier criteria:**
| Tier | Criteria |
|------|----------|
| Screening | ratio < 0.1 |
| Confirmation | ratio < 0.1 AND basin_volume ≥ 0.8 AND null p < 0.01 |
| Definitive | confirmation AND p ≤ 0.005 AND d > 1.0 AND relax > 1 AND (if pert: restore ≥ 0.8) |

**Trivial-collapse gate:** Rejects systems with relaxation_time ≤ 1 step.
Prevents FPs from instant constant maps.

**T1a adapter:** `extract_observation_bundle()` groups records by trial,
extracts per-trial ICs, finals, convergence status, steps-to-convergence.

**Registered:** P25 in `DETECTOR_REGISTRY`, `canalization_landscape` substrate.

---

## Part C — Tests

**`tests/test_canalization_p25_e2e.py` (9 tests):**
- `TestP25CanonicalPositive`: 4 tests (detected at definitive, ratio < 0.01,
  basin ≥ 0.9, relaxation > 1)
- `TestP25NegativeControls`: 2 tests (diffusive NOT detected, trivial NOT detected)
- `TestP25PerturbationRecovery`: 1 test (restoration ≥ 0.8)
- `TestP25Determinism`: 1 test (same seed → identical results)
- `TestP25ObservationBundle`: 1 test (shape validation)

**`tests/test_cross_model.py` +3 T1b tests:**
- P25 on CanalizedLandscape: detected
- P25 on MultiBasinGRN: detected
- P25 on DiffusiveDynamics: NOT detected

**Transfer matrix:** +1 detector (27 total), +2 models (31 total),
+1 substrate (12 total), 112 compatible pairs, 837 cells.

---

## Part D — dim1 reproduction

`analysis/outputs/p25_canalization_reproduction.json`:
- IC variance: 20.25
- Final variance: ≈ 0 (machine precision)
- Convergence ratio: ≈ 0 (< 0.10 PASS)
- Basin volume: 1.0 (≥ 0.80 PASS)
- Detector: DEFINITIVE (confidence 0.90)
- **passes_tolerance: true**

---

## Part E — dim2 + dim3

**dim2** (`analysis/outputs/p25_multiseed.json`): 20 seeds (0–19).
- Convergence ratio: 0.000 ± 0.000 (CV ≈ 0%)
- Basin volume: 1.000 ± 0.000
- All 20 seeds DEFINITIVE (confidence 0.90)

**dim3** (`docs/methods_notes/p25_methods.md`): gradient-flow + quartic
potential, IC-distribution surrogate null, trivial-collapse gate,
CanalizedLandscape vs MultiBasinGRN, negative control design.

---

## Part F — Documentation updates

- **depth_gap.md:** +P25 row (GAP, dim4 pending). Audited 26→27. Gap 1→2.
- **observation_schema.md:** +canalization observation bundle (P25).
- **paper_section4_draft.md:** +§4.25 P25 Sprint 81.
- **paper_section6_draft.md:** Sprint 81 entry. AT-DEPTH 25/27.
- **paper_CHANGELOG.md:** Sprint 81 entry.
- **REPLICATION_NOTES.md:** P25 replication section added.

---

## Post-flight checks

- `pytest tests/test_canalization_p25_e2e.py tests/test_orchestration.py tests/test_transfer_matrix_counts.py tests/test_cross_model.py`: **125 passed** ✓
- `p25_canalization_reproduction.json`: passes_tolerance=true ✓
- `p25_multiseed.json`: 20/20 DEFINITIVE ✓

---

## Open Carry-Forwards

- **C-p7-time-shuffled-fp:** `time_shuffled` FP at screening (Sprint 66, low priority)
- **C-p19-bias-zero-chance-alignment:** 1/5 bias_zero Class C FP at confirmation (Sprint 70, 4% rate)
- **P12 dim1:** documented finite-size measurement limitation (accepted)

---

## Files added/modified

**New (7):**
- `epc/models/canalization.py` — CanalizedLandscape + MultiBasinGRN + DiffusiveDynamics + TrivialCollapse
- `epc/detectors/p25_equifinality.py` — P25 equifinality detector + T1a adapter
- `tests/test_canalization_p25_e2e.py` — 9 e2e tests
- `analysis/reproductions/p25_canalization.py` — dim1 reproduction
- `analysis/reproductions/p25_multiseed.py` — dim2 multi-seed campaign
- `analysis/outputs/p25_canalization_reproduction.json` — dim1 results
- `analysis/outputs/p25_multiseed.json` — dim2 results
- `docs/methods_notes/p25_methods.md` — dim3 methods note

**Modified (10):**
- `epc/orchestration.py` — +canalization_landscape substrate, +2 models, +P25 detector
- `tests/test_orchestration.py` — count updates + P25 canonical pair
- `tests/test_transfer_matrix_counts.py` — EXPECTED counts updated
- `tests/test_cross_model.py` — +3 T1b tests
- `docs/depth_gap.md` — +P25 row, aggregate updates
- `docs/observation_schema.md` — +canalization bundle
- `docs/paper_section4_draft.md` — +§4.25
- `docs/paper_section6_draft.md` — Sprint 81 entry
- `docs/paper_CHANGELOG.md` — Sprint 81 entry
- `REPLICATION_NOTES.md` — P25 replication section

---

**Decision: GO** — Sprint completed cleanly. P25 implemented with DEFINITIVE detection on canonical model, two negative controls correctly rejected (diffusive + trivial collapse), T1b cross-model (MultiBasinGRN) DEFINITIVE, 20/20 multi-seed DEFINITIVE. dim4 Phase-2a panel scheduled for Sprint 82. Chain may proceed autonomously.
