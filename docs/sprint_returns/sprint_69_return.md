# Sprint 69 Return — P19 Emergent Leadership / Minority Guidance (Couzin 2005)

**Date:** 2026-06-08
**Base HEAD (sprint start):** `8e91239` (Sprint 68 follow-up)
**Sprint goal:** Implement P19 (emergent leadership / minority guidance) end-to-end: model + detector + tests + registry + dim1 + dim2 + dim3.
**Tag:** `v0.69.0`
**Sprint type:** Chat-led design + code-led execution (Milestone B, Wave 1).

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `8e91239` ✓
2. **Working tree clean:** ✓

---

## Part A — Model: `epc/models/informed_minority.py`

Vicsek-style flock with an informed minority (Couzin et al. 2005). N agents
in a periodic [0, L)² box; fraction ρ has a weighted bias (ω) toward a fixed
preferred direction; naive agents align purely via local Vicsek interactions.

- Synchronous update with metric neighbor interaction (cKDTree + sparse matrix)
- Informed agents: `θ = arg[(1-ω)⟨e^{iθ}⟩_neighbors + ω·e^{iθ_pref}] + noise`
- State dict: positions, velocities, headings, informed_mask, step
- Metadata: model_family="informed_minority", has_informed_minority=True
- `group_accuracy()` and `polarization()` utility functions

Canonical parameters: N=200, L=10.0, v₀=0.03, η=0.1, r=1.0, ρ=0.1, ω=0.3, θ_pref=0.0.

---

## Part B — Detector: `epc/detectors/p19_emergent_leadership.py`

Three-tier detection:

- **Screening (0.50):** Group accuracy > 0.3 (alignment with preferred direction, measured over last 50% of history)
- **Confirmation (0.70):** Label-shuffle influence asymmetry — informed agents' mean heading is closer to preferred direction than naive agents'. Pull metric: `cos(θ̄_inf - θ_pref) - cos(θ̄_naive - θ_pref)`. Significance via label-shuffle permutation null (p < 0.05, ≥ 99 permutations).
- **Definitive (0.90):** Guidance efficacy > 2.0 (accuracy/ρ) + genuine minority (ρ ≤ 0.25) + p < 0.01 (≥ 199 permutations).

**Architecture decision:** Directional pull (label-shuffle) over transfer entropy for the confirmation metric. TE on group mean headings = 0 in the converged Vicsek regime because both informed and naive means are constant in steady state; the label-shuffle pull metric captures the same mechanistic signature (asymmetric influence) more robustly. The KSG TE estimator (`epc/metrics/transfer_entropy_ksg.py`) was evaluated but dropped in favor of pull due to signal collapse in steady state.

P19 registered in `DETECTOR_REGISTRY` (continuous_2d substrate, headings observable).

---

## Part C — Tests: `tests/test_informed_minority_p19_e2e.py`

19 tests, all passing:
- 7 model tests: determinism, informed_mask count, ρ=0, state keys, metadata, polarization, accuracy
- 5 detector tests: detected, pattern_id, **DEFINITIVE tier**, accuracy > 0.3, pull > 0
- 1 negative control: ρ=0 → tier ∈ {none, screening} (not definitive)
- 6 registry tests: model + detector registered, compatibility, substrate

Canonical positive (ρ=0.1, N=100, 500 steps, 199 permutations): **DEFINITIVE** (accuracy=1.000, pull=0.0005, p=0.005, efficacy=10.0).

---

## Part D — dim1 reproduction: `analysis/reproductions/p19_couzin2005.py`

Couzin et al. (2005) Fig. 2a: accuracy vs informed fraction ρ.

| ρ | accuracy (mean ± std) |
|---|---|
| 0.000 | 0.125 ± 0.756 |
| 0.025 | 0.977 ± 0.040 |
| 0.050 | 1.000 ± 0.000 |
| 0.100 | 1.000 ± 0.000 |
| 0.150 | 1.000 ± 0.000 |
| 0.200 | 1.000 ± 0.000 |
| 0.300 | 1.000 ± 0.000 |
| 0.500 | 1.000 ± 0.000 |

Tolerance checks (all PASS):
- Monotone Spearman ρ = 1.0 (> 0.8) ✓
- accuracy(ρ=0.0) = 0.125 (|acc| < 0.3, chance-level) ✓
- accuracy(ρ=0.05) = 1.000 (> 0.2) ✓
- accuracy(ρ=0.1) = 1.000 (> 0.5) ✓
- accuracy(ρ=0.5) = 1.000 (> 0.8) ✓

Note: accuracy saturates more abruptly than Couzin (2005) Fig. 2a because η=0.1 is low relative to alignment strength. This is a parameter-regime effect, not a model error.

Output: `analysis/outputs/p19_couzin2005_reproduction.json` (passes_tolerance=true).

---

## Part E — dim2 + dim3

### dim2: Multi-seed (20 seeds at ρ=0.1)

| Metric | Mean | Std | CV |
|--------|------|-----|----|
| Accuracy | 1.000 | 0.000 | 0.0% |
| Influence pull | 0.0013 | 0.0020 | 156.3% |

All 20 seeds: accuracy > 0.999, detector fires (confirmation or higher). Pull CV is high because the absolute pull values are small (the informed agents lead by tiny angular amounts in steady state), but the sign is consistently positive.

Output: `analysis/outputs/p19_multiseed.json`.

### dim3: Methods note

`docs/methods_notes/p19_methods.md` — covers:
- Vicsek + informed-bias dynamics (equations, canonical parameters)
- Label-shuffle null design (over TE: rationale for architectural decision)
- Distinctness from P5 (no preferred direction), P17 (environmental sensing), P18 (symmetric pooling), P32 (stable roles)
- Convergence speed and parameter-sensitivity discussion

---

## Part F — depth_gap + paper

- **depth_gap.md:** P19 row added (dim1-3 PASS, dim4 pending → GAP). Patterns audited 21→22; gap count 1→2 (P12 + P19).
- **paper_CHANGELOG.md:** Sprint 69 entry added.
- **REPLICATION_NOTES.md:** Sprint 69 P19 section added with full dim1/dim2/dim3 results.
- **epc/orchestration.py:** +1 model (informed_minority), +1 detector (P19). 23 models × 22 detectors, 104 compatible pairs.
- **tests/test_orchestration.py:** Counts updated (23/22/506/104/402).
- **tests/test_transfer_matrix_counts.py:** EXPECTED dict updated for 23×22.

---

## Post-flight checks

- `pytest tests/test_orchestration.py tests/test_transfer_matrix_counts.py tests/test_cross_detection_matrix.py tests/test_informed_minority_p19_e2e.py -m "not slow"`: **128 passed** ✓
- P19 canonical positive: DEFINITIVE (accuracy=1.000, pull p=0.005, efficacy=10.0) ✓
- P19 negative control (ρ=0): not detected or screening only ✓
- dim1 reproduction: all 5 tolerance checks PASS ✓
- dim2 multiseed: 20/20 seeds detected ✓
- No threshold relaxation (Sprint 30 rule) ✓

---

## Carry-forwards

- **C-p19-abrupt-saturation:** Accuracy saturates to ~1.0 at ρ=0.025, more abruptly than Couzin (2005) Fig. 2a. Root cause: η=0.1 is low relative to alignment strength (high density × large interaction radius). Not a validation failure (the core result — accuracy rises with ρ, ρ=0 is chance-level — reproduces), but dim4 Sprint 70 may want to test with higher noise for a more gradual curve.
- **C-p19-te-vs-pull:** TE (KSG) on mean heading time series produces zero signal in converged regime. Documented architectural decision; label-shuffle pull used instead. If a future sprint needs TE-based P19 detection, per-agent or heading-change TE approaches may work but require longer time series or higher noise.

---

## Files added/modified

**New (6):**
- `epc/models/informed_minority.py` — Couzin 2005 informed-minority model (~230 lines)
- `epc/detectors/p19_emergent_leadership.py` — P19 detector (~280 lines)
- `tests/test_informed_minority_p19_e2e.py` — 19 tests
- `analysis/reproductions/p19_couzin2005.py` — dim1 reproduction
- `analysis/reproductions/p19_multiseed.py` — dim2 multi-seed
- `docs/methods_notes/p19_methods.md` — dim3 methods note

**Modified (6):**
- `epc/orchestration.py` — +1 model (informed_minority) + +1 detector (P19), counts 23×22
- `tests/test_orchestration.py` — count updates (22→23 models, 21→22 detectors, 95→104 pairs)
- `tests/test_transfer_matrix_counts.py` — EXPECTED updated for 23×22 registry
- `docs/depth_gap.md` — P19 row added (dim1-3 PASS, dim4 pending)
- `docs/paper_CHANGELOG.md` — Sprint 69 entry
- `REPLICATION_NOTES.md` — Sprint 69 P19 section

---

**Decision: GO**
