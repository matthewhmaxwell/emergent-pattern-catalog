# Sprint 86 Return — P4 Territoriality / Exclusion Boundaries

**Date:** 2026-06-10
**Base HEAD (sprint start):** `feee9a7` (Sprint 85 follow-up)
**Tag:** `v0.86.0`
**Sprint type:** Chat-led design + code-led execution.

## Part A — Model: `epc/models/territoriality.py`

Three model classes implemented:

1. **ScentMarkingModel** (canonical) — N agents on L×L torus with two-stage
   territorial movement: (1) hard avoidance of cells where foreign scent >
   repulsion_strength × deposition_rate, (2) own-scent-weighted softmax
   selection among allowed cells. Default: N=4, L=48, 20000 steps. Produces
   strong exclusive territories (excl=0.902, overlap=0.034, persistence=0.865).

2. **PheromoneRepulsionModel** (T1b) — Independent implementation with hard
   avoidance threshold + uniform random selection (no home-attraction bias).
   Default: N=4, L=48, dep=5.0, decay=0.005, threshold=0.3, 5000 steps.
   Exclusivity=0.866.

3. **PlainRandomWalkModel** (negative control) — Deposits scent but ignores
   it. Uniform random movement. Exclusivity=0.742 (below detection threshold).

## Part B — Detector: `epc/detectors/p4_territoriality.py`

- T1a observation bundle adapter: `extract_observation_bundle(history)`
- Primary metrics: exclusivity index, pairwise overlap, boundary persistence,
  occupancy-scent correlation
- Null model: cell-level multinomial shuffle (redistribute per-cell visit
  counts uniformly among agents)
- Tier criteria:
  - Screening: exclusivity > 0.80 AND p < 0.10
  - Confirmation: persistence > 0.6, p < 0.05, Cohen's d > 1.0
  - Definitive: occ-scent corr < −0.005, p ≤ 0.005, metadata confirms avoidance
- P4 registered in DETECTOR_REGISTRY

## Part C — Tests

- `tests/test_territoriality_p4_e2e.py`: 12 tests — DEFINITIVE on canonical,
  negative control NOT detected, persistence, bundle adapter, metadata,
  transfer matrix registration. All pass.
- `tests/test_cross_model.py::TestP4OnModels`: 3 tests — ScentMarking
  DEFINITIVE, PheromoneRepulsion detected (CONFIRMATION), RandomWalk NOT
  detected. All pass.
- Transfer matrix counts updated: 35 models × 29 detectors = 1015 cells,
  116 compatible pairs. All 22 transfer matrix tests pass.

## Part D — dim1 reproduction

`analysis/reproductions/p4_giuggioli2011.py` →
`analysis/outputs/p4_giuggioli2011_reproduction.json`:
- Exclusivity = 0.9023 (> 0.85 PASS)
- Pairwise overlap = 0.0336 (< 0.10 PASS)
- Boundary persistence = 0.8652 (> 0.70 PASS)
- Detector tier = definitive (PASS)
- Cohen's d = 157.5
- **passes_tolerance = true**

## Part E — dim2 + dim3

**dim2:** 20-seed campaign (`analysis/outputs/p4_multiseed.json`):
- Exclusivity: 0.925 ± 0.016 (CV = 1.7%)
- Overlap: 0.025 ± 0.007
- Persistence: 0.817 ± 0.067
- 20/20 detected (2 DEFINITIVE, 18 CONFIRMATION)

**dim3:** `docs/methods_notes/p4_methods.md` — covers two-stage movement,
multinomial null, tier criteria, reproduction results.

## Part F — Registry + docs

- `epc/orchestration.py`: new substrate `territorial_agent_field`, 2 model
  entries (scent_marking_territory, pheromone_repulsion_territory), 1 detector
  entry (P4). Counts: 35 models, 29 detectors, 14 substrate types.
- `docs/depth_gap.md`: P4 row added. AT-DEPTH count: 28/29.
- `docs/observation_schema.md`: territorial-agent-field bundle documented.
- `docs/paper_section4_draft.md`: §4.4 P4 section added.

## Carry-forwards

- C-p4-definitive-rate: Only 2/20 seeds reach DEFINITIVE (10%). The
  occ-scent correlation threshold (< −0.005) is strict; most seeds land at
  −0.003 to −0.006. This is because well-formed territories create a
  floor effect (minimal overlap → weak correlation). CONFIRMATION at 90%
  of seeds is sufficient. Low priority — does not affect detection.

## Summary

| Metric | Value |
|--------|-------|
| New files | 8 (model, detector, e2e test, reproduction, multiseed, methods note, 2 JSON outputs) |
| Modified files | 7 (orchestration, cross-model test, orch test, transfer matrix test, depth_gap, obs schema, paper) |
| Tests added | 15 (12 e2e + 3 cross-model) |
| Tests passing | 37/37 (P4-specific + transfer matrix) |
| Implemented patterns | 29 |
| AT-DEPTH | 28 (P12 sole GAP) |

**Decision: GO** — P4 territoriality implemented end-to-end with DEFINITIVE
detection on canonical model, T1b cross-model fires, negative control
correctly rejected, all tolerances pass. dim4 Phase-2a panel deferred to
Sprint 87.
