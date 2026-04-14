# Emergent Pattern Catalog — Claude Code Context

## Project Summary
Systematic research program to catalog, detect, and discover emergent behavioral
competencies in minimal agent-based systems. Building a "periodic table" of
emergent patterns and a computational toolkit for detecting them.

## Current Status (Sprint 3 — Collective Motion)
- Pattern catalog: 32 atomic patterns across 10 clusters, 3 layers, 11 dimensions
- Architecture: Layer 1 (behaviors), Layer 2A (math descriptors), Layer 2B
  (cognitive-analogue annotations), Layer 3 (meta-capacities)
- Models: 5 implemented (Zhang sorting, GH excitable CA, Game of Life,
  Vicsek flocking, D'Orsogna SPP milling)
- Metrics: 13 implemented (aggregation, DG, wave propagation ×4, TE ×2,
  collective motion ×5)
- Detectors: 5 implemented (P1, P13, P31, P5, P6) + P13/P15 discriminator
- See docs/pattern_catalog.md for the full catalog

## Code Layout

All real implementations live in `epc/` — ignore `models/`, `metrics/`,
`analysis/` at repo root (legacy stubs from initial commit).

```
epc/
├── detector_result.py         # DetectorResult dataclass
├── base_detector.py           # Base detect() pipeline
├── base_model.py              # setup/step/run/get_metadata/get_timescale
├── base_metric.py             # compute/required_keys/validate_history
├── models/
│   ├── cell_view_sorting.py        # Zhang 1D sorting (Sprint 1)
│   ├── cell_view_sorting_threaded.py # Threading variant (Sprint 1)
│   ├── greenberg_hastings.py       # Excitable CA (Sprint 2)
│   ├── game_of_life.py             # Conway B3/S23 (Sprint 2)
│   ├── vicsek.py                   # Standard Vicsek Model (Sprint 3)
│   └── dorsogna_spp.py             # D'Orsogna SPP + Morse (Sprint 3)
├── metrics/
│   ├── aggregation.py              # Moran's I, segregation, clusters (Sprint 1)
│   ├── delayed_gratification.py    # Zhang's DG metric (Sprint 1)
│   ├── wave_propagation.py         # Wavefront speed, spiral tips (Sprint 2)
│   ├── transfer_entropy.py         # Plug-in TE estimator (Sprint 2)
│   └── collective_motion.py        # φ, R, L, autocorrelation, nematic (Sprint 3)
└── detectors/
    ├── p1_aggregation.py           # P1 with label-shuffle null (Sprint 1)
    ├── p31_delayed_gratification.py # P31 + non-redundancy test (Sprint 1)
    ├── p13_excitable_wave.py       # P13 with wavefront/spiral (Sprint 2)
    ├── p13_p15_discriminator.py    # TE discriminator (Sprint 2)
    ├── p5_flocking.py              # P5 3-tier, shuffle null (Sprint 3)
    └── p6_milling.py               # P6 3-tier, ring density (Sprint 3)
```

## Key Design Decisions
- Models return full state history as list of dicts with numpy arrays
- Metrics operate on state histories, not live simulations
- All randomness seeded for reproducibility
- DetectorResult is universal output; confidence is tier-capped
- Continuous-space models use `positions` (N,2) + `velocities` (N,2) state keys
- All persistence thresholds use system-intrinsic timescales (T_cross, T_prop, etc.)
- See PROJECT_STATUS.md for full architecture decisions log (23 decisions)

## Statistical Power Requirements (NON-NEGOTIABLE)
- P1 detector: n_permutations=999
- P13 detector: n_null_runs=199
- P5/P6 detectors: n_permutations=199 (floor p=0.005)
- P31 non-redundancy: ≥500 runs for 10-fold CV
- TE discriminator: n_permutations=99
- Any new permutation test: ≥99 perms (floor p=0.01)
- Any new cross-validation: ≥500 runs, ≥10 folds

## Pattern Clusters (32 total)
- A: Spatial organization (P1-P4)
- B: Collective motion (P5-P8) ← Sprint 3 focus
- C: Temporal dynamics (P9-P12)
- D: Wave propagation (P13-P14)
- E: Information processing (P15-P17)
- F: Decision-making (P18-P23)
- G: Resilience (P24-P26)
- H: Competition/cooperation (P27-P28)
- I: Structure formation (P29-P30)
- J: Agent-level competencies (P31-P32)

## Coding Standards
- Type hints on all function signatures
- Docstrings on all public methods
- NumPy for array operations
- Matplotlib for plotting
- SciPy for spatial (cKDTree) and sparse operations
- No heavy frameworks (no Mesa dependency — keep it lightweight)
- Tests for each model and metric
