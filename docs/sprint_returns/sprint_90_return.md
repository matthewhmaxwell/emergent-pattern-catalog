# Sprint 90 Return — P32 Emergent Specialization / Division of Labor

**Date:** 2026-06-10
**Base HEAD (sprint start):** `e0e3dc2` (Sprint 89 follow-up)
**Tag:** `v0.90.0`
**Sprint type:** Chat-led design + code-led execution.

## Part A — Model: `epc/models/division_of_labor.py`

ResponseThresholdModel implementing Bonabeau, Theraulaz & Deneubourg (1996):
- N identical agents with per-task response thresholds θ_{i,j}
- Response probability: s²/(s² + θ²) per Bonabeau Eq. 1
- Threshold reinforcement: performing task j → θ_{i,j} decreases by ξ
- Threshold forgetting: not performing task j → θ_{i,j} increases by φ
- Stimulus dynamics: accumulate at base rate, decrease by worker count

NoReinforcementModel: identical dynamics but ξ=0, φ=0 (negative control).

Helper functions: `compute_per_agent_entropy`, `compute_windowed_entropy`,
`compute_efficiency`.

## Part B — Detector: `epc/detectors/p32_specialization.py`

P32SpecializationDetector with three-tier pipeline:
- **Screening:** per-agent entropy decline > 0.1 bits (early vs late window)
- **Confirmation:** null p < 0.01 AND role_diversity ≥ 0.5 AND NOT
  single_task_collapse AND late_coverage ≥ 0.8 AND late_entropy < 1.0
- **Definitive:** null p < 0.005 AND entropy_decline > 0.3 AND
  (late_efficiency > 0.3 OR efficiency_gain > 0) AND role_diversity ≥ 0.67

Null model: time-shuffle per-agent assignment series (NullType.SHUFFLE).
Exclusion: P23 (anti-coordination) excluded via late switching frequency < 0.3.

T1a adapter: `extract_observation_bundle(history)` → per-agent task-assignment
time series. Registered as P32 in DETECTOR_REGISTRY.

## Part C — Tests

`tests/test_division_of_labor_p32_e2e.py`: **25 tests, all PASS**
- Model determinism (2 tests)
- Entropy decline with/without reinforcement (2 tests)
- Threshold differentiation (1 test)
- Metadata + state keys (2 tests)
- Canonical DEFINITIVE (4 tests: tier, entropy, role diversity, P23 exclusion)
- Negative controls (3 tests: no reinforcement, single-task collapse, short run)
- T1a adapter (2 tests)
- Orchestration registration (5 tests)
- Helper functions (4 tests)

`tests/test_cross_model.py`: **3 T1b tests, all PASS**
- ResponseThresholdModel detected
- Multiplicative-threshold variant detected (OOD-readiness)
- NoReinforcementModel not detected

`tests/test_cross_detection_matrix.py`: **P32 column added (2 cells)**
- `response_threshold × P32`: detected (DEFINITIVE)
- `no_reinforcement × P32`: not_detected

Transfer matrix total: ≥175 audited cells.

## Part D — dim1 reproduction

`analysis/reproductions/p32_bonabeau1996.py` → `analysis/outputs/p32_bonabeau1996_reproduction.json`

Parameters: N=20, M=3, ξ=0.15, φ=0.01, 1000 steps, seed=42.

| Tolerance | Measured | Threshold | Result |
|-----------|----------|-----------|--------|
| Late entropy < 50% max | 0.516 | < 0.792 | **PASS** |
| Late coverage | 1.000 | ≥ 1.0 | **PASS** |
| Efficiency gain | 0.344 | > 0 | **PASS** |

`passes_tolerance: true`

## Part E — dim2 + dim3

**dim2:** `analysis/reproductions/p32_multiseed.py` → `analysis/outputs/p32_multiseed.json`

20 seeds at ξ=0.05, φ=0.01, 1000 steps:

| Metric | Mean | Std | CV |
|--------|------|-----|----|
| Final per-agent entropy | 1.210 | 0.208 | 17.2% |
| Efficiency gain | 0.038 | 0.042 | 111.7% |

Detection: 4/20 confirmation, 16/20 screening, 0 definitive at 199 perms.

**dim3:** `docs/methods_notes/p32_methods.md` — covers response-threshold
dynamics, three-tier criteria, T1a observation contract, T1b cross-model
validation, parameter justifications.

## Part F — depth_gap + paper + schema

- `docs/depth_gap.md`: P32 row added (dim1–dim3 PASS, dim4 pending → GAP).
  Implemented-count 30→31. AT-DEPTH 29/31. Gap 1→2 (P12, P32).
- `docs/observation_schema.md`: per-agent task-assignment bundle (P32) added.
- `docs/paper_section4_draft.md`: §4.32 added.
- `docs/paper_section6_draft.md`: Sprint 90 entry added.
- `docs/paper_CHANGELOG.md`: Sprint 90 entry added.
- `REPLICATION_NOTES.md`: Sprint 90 P32 section added.
- `epc/orchestration.py`: task_allocation_timeseries substrate (16th),
  response_threshold + no_reinforcement models, P32 detector registered.
  39 models × 31 detectors.

## Files changed

**New (6):**
- `epc/models/division_of_labor.py` — ResponseThresholdModel + NoReinforcementModel
- `epc/detectors/p32_specialization.py` — P32 detector + T1a adapter
- `tests/test_division_of_labor_p32_e2e.py` — 25 e2e tests
- `analysis/reproductions/p32_bonabeau1996.py` — dim1 reproduction
- `analysis/reproductions/p32_multiseed.py` — dim2 multiseed campaign
- `docs/methods_notes/p32_methods.md` — dim3 methods note

**Modified (9):**
- `epc/orchestration.py` — P32 + task_allocation_timeseries registered
- `tests/test_cross_model.py` — P32 T1b cross-model tests
- `tests/test_cross_detection_matrix.py` — P32 column (2 cells), total ≥175
- `docs/depth_gap.md` — P32 row, counts updated
- `docs/observation_schema.md` — P32 bundle section
- `docs/paper_section4_draft.md` — §4.32
- `docs/paper_section6_draft.md` — Sprint 90 entry
- `docs/paper_CHANGELOG.md` — Sprint 90 entry
- `REPLICATION_NOTES.md` — Sprint 90 section

## Carry-forwards

- **C-p12-dim1**: P12 (RPS λ ∝ √M scaling) remains sole dim1 GAP pattern.
  Documented finite-size measurement limitation (Sprints 54/58/59/63).
- **C-p32-dim4**: P32 dim4 Phase-2a panel pending (Sprint 91). P32 is the
  31st pattern; dim4 will validate detector specificity against synthetic
  nulls + failed regimes.
- **C-p32-multiseed-low-confirmation**: Only 4/20 seeds reach confirmation
  at default parameters (ξ=0.05, 199 perms, 1000 steps). Specialization
  strength is seed-dependent due to stochastic task encounters. Canonical
  positive at 499 perms reliably reaches DEFINITIVE.

## Summary

| Metric | Value |
|--------|-------|
| New files | 6 |
| Modified files | 9 |
| Sprint-specific tests | 30/30 PASS |
| dim1 reproduction | PASS (all 3 tolerances) |
| dim2 multiseed | 20 seeds, entropy CV=17.2% |
| AT-DEPTH count | 29 / 31 |
| Gap count | 2 (P12, P32) |
| Registry | 39 models × 31 detectors |

**Decision: GO** — P32 emergent specialization implemented end-to-end.
Model + detector + tests + dim1–dim3 + T1a/T1b OOD-readiness all complete.
Canonical positive reaches DEFINITIVE (p=0.002, d=8.88). Negative controls
correctly rejected. Chain may proceed to Sprint 91 (P32 dim4 Phase-2a panel).
