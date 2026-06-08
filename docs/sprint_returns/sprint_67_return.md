# Sprint 67 Return — P17 Collective Gradient Sensing (Berdahl 2013) + dim1/dim2/dim3

**Date:** 2026-06-08
**Base HEAD (sprint start):** `3e95509` (Sprint 66 follow-up)
**Sprint goal:** Implement P17 — Distributed sensing / collective gradient detection end-to-end: model + detector + tests + registry + dim1 + dim2 + dim3.
**Tag:** `v0.67.0`
**Sprint type:** Chat-led design + code-led execution (Milestone B, Wave 1).

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `3e95509` ✓
2. **Working tree clean:** ✓
3. **Transfer matrix pre-sprint:** 21 models × 20 detectors, 86 compatible pairs ✓

---

## Part A — Model: `epc/models/collective_sensing.py`

Berdahl et al. (2013) speed-modulation mechanism for collective gradient sensing:
- Agents in periodic 2D domain with Gaussian scalar field (resource peak)
- Each agent senses local field value with additive Gaussian noise (σ=0.8)
- **Speed modulation:** agents slow down in high-field regions (α=0.95)
- **Social attraction:** heading biased toward group centroid (strength=0.2)
- The mechanism: agents near peak slow down → group CoM drifts toward peak
- SNR ∝ √N: individual SNR ≈ 0.73 (unreliable), group N=50 SNR ≈ 5.2 (strong)

Canonical parameters: box_size=20, v_max=0.4, turn_noise=0.3, sensing_noise=0.8,
alpha=0.95, social_strength=0.2, field_sigma=5.0, n_steps=1000.

---

## Part B — Detector: `epc/detectors/p17_collective_sensing.py`

Group-size scaling test (the Berdahl signature):
- Runs model at N ∈ {1, 5, 10, 25, 50} with multiple seeds
- Computes chemotactic index (CI) at each N
- Null model: α=0 (no speed modulation) — preserves social structure, removes mechanism
- Screening: CI(N_max) > CI(N=1) + 0.05 AND CI(N_max) > 0.1
- Confirmation: positive CI-vs-log(N) slope + null p < 0.05
- Definitive: monotonic + Cohen's d > 2.0 + p < 0.01 + CI(N_max) > 0.2

Registered as P17 in DETECTOR_REGISTRY (continuous_2d, requires positions + headings).

---

## Part C — Tests: `tests/test_collective_sensing_p17_e2e.py`

12 tests (10 fast + 2 slow):
- **TestModelDeterminism**: seed reproducibility ✓
- **TestGroupSizeScaling**: CI(N=50) > CI(N=1) + 0.1 ✓
- **TestNegativeControl**: α=0 → no gradient climbing; N=1 → at chance ✓
- **TestDetectorCanonical** (slow): confirmation/definitive on canonical regime
- **TestRegistration**: model + detector registered, compatible ✓

All fast tests pass in ~3.6 seconds.

---

## Part D — dim1 reproduction: `analysis/reproductions/p17_berdahl2013.py`

**Target:** Berdahl et al. (2013) Fig. 1 — CI rises with group size N.

**Results (10 seeds per N):**
| N | CI mean ± std |
|---|---|
| 1 | −0.167 ± 0.431 |
| 5 | +0.255 ± 0.375 |
| 10 | +0.224 ± 0.332 |
| 25 | +0.300 ± 0.145 |
| 50 | +0.396 ± 0.122 |

**Tolerance checks:**
- Slope CI vs log(N) = 0.133 (> 0.02) ✓
- Spearman ρ = 0.90 (> 0.7, p = 0.037) ✓
- |CI(N=1)| = 0.167 (< 0.20, chance-level) ✓
- CI(N=50) = 0.396 (> 0.15) ✓

**`passes_tolerance`: True**

Output: `analysis/outputs/p17_berdahl2013_reproduction.json`

---

## Part E — dim2 multi-seed + dim3 methods note

### dim2: `analysis/outputs/p17_multiseed.json`
20-seed campaign at canonical N=50:
- CI mean: 0.394
- CI std: 0.130
- CI CV: 33.0%
- Fraction positive: 100% (20/20)

### dim3: `docs/methods_notes/p17_methods.md`
Full methods note covering: speed-modulation mechanism, SNR ∝ √N derivation,
α=0 null model rationale, parameter table, tier criteria, distinctness from P5,
ADR 57/58/59 (speed-modulation choice, null design, periodic CoM).

---

## Part F — Registry + depth_gap + paper CHANGELOG

- **Orchestration:** +1 model (`collective_sensing`, continuous_2d), +1 detector (P17)
- **Counts:** 22 models × 21 detectors, 95 compatible pairs
- **depth_gap.md:** P17 row added (dim1 PASS, dim2 PASS, dim3 PASS, dim4 pending). Patterns audited 20→21; gap count 1→2.
- **paper_CHANGELOG.md:** Sprint 67 entry added
- **test_transfer_matrix_counts.py:** EXPECTED updated
- **test_orchestration.py:** Hardcoded counts updated

---

## Post-flight checks

- `pytest tests/ -m "not slow"` (excluding cross-detection matrix): all pass ✓
- Transfer matrix counts verified: 22×21, 95 pairs ✓
- dim1 `passes_tolerance`: True ✓
- dim2 multiseed output exists with positive CI ✓
- dim3 methods note written ✓

---

## Carry-forwards

- **C-p17-dim4-pending:** P17 dim4 Phase-2a panel not yet run (Sprint 68 scope). Class B comparison: P5 flocking (motion without field-inference) must be rejected.

---

## Files added/modified

**New (8):**
- `epc/models/collective_sensing.py` — Berdahl 2013 collective sensing model
- `epc/detectors/p17_collective_sensing.py` — P17 detector (N-scaling CI test)
- `tests/test_collective_sensing_p17_e2e.py` — 12 tests
- `analysis/reproductions/p17_berdahl2013.py` — dim1 reproduction script
- `analysis/reproductions/p17_multiseed.py` — dim2 multi-seed script
- `analysis/outputs/p17_berdahl2013_reproduction.json` — dim1 output
- `analysis/outputs/p17_multiseed.json` — dim2 output
- `docs/methods_notes/p17_methods.md` — dim3 methods note

**Modified (5):**
- `epc/orchestration.py` — +1 model, +1 detector, docstring count update
- `tests/test_orchestration.py` — count updates (21→22 models, 20→21 detectors, 86→95 pairs)
- `tests/test_transfer_matrix_counts.py` — EXPECTED dict updated
- `docs/depth_gap.md` — P17 row, patterns audited 20→21
- `docs/paper_CHANGELOG.md` — Sprint 67 entry

---

**Decision: GO**
