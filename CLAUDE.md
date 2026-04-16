# Emergent Pattern Catalog — Claude Code Context

## Project Summary
Systematic research program to catalog, detect, and discover emergent behavioral
competencies in minimal agent-based systems. Building a "periodic table" of
emergent patterns and a computational toolkit for detecting them.

Repository: https://github.com/matthewhmaxwell/emergent-pattern-catalog

## Current Status — Sprint 9 Complete

Sprint 9 added spatial Rock-Paper-Scissors (Reichenbach 2007) as the canonical
P12 positive and built the P12 cyclic-dominance detector. Core finding: P13
already rejects RPS cleanly at screening via wavefront-speed CV; the predicted
P13 false positive does NOT occur. P12 complements P13 via a neighbor-
conditional replacement ratio ρ which is ≈ 70–200 for RPS dominance edges and
exactly 1.0 for GH clock transitions — a two-order-of-magnitude separation.

### Inventory
- Pattern catalog: 32 atomic patterns across 10 clusters, 3 layers, 11 dimensions
- Models: 13 models (lattice_1d=2, lattice_2d=7, continuous_2d=2, oscillator=1,
  opinion_space=1) — all validated against published literature
- Metrics: 9 modules, 22+ metrics — all verified against canonical results
- Detectors: 12 detectors + 1 discriminator — all with 3-tier detection + null models
- Orchestration: 5 substrate types, 13 models × 12 detectors → **42 compatible pairs**
- Cross-detection matrix: **27 audited pairs** pinned in `tests/test_cross_detection_matrix.py`
- Test suite: ~215 tests across 21 files, all passing (a few slow-marked)
- Cluster coverage: A, B, C, D, E, F, H, J (missing: G Resilience, I Structure Formation)
- Paper Sections 3–7 drafted (~8,300 words); Section 4 needs Sprint 9 matrix update

### Files added/modified in Sprint 9
New (5):
- `epc/models/rps_spatial.py` — Reichenbach 2007 RPS (May-Leonard), ~440 lines
- `epc/detectors/p12_cyclic_dominance.py` — P12 detector, ~460 lines
- `tests/test_rps_replication.py` — 11 tests (10 fast + 1 slow)
- `tests/test_rps_p13_boundary.py` — 6 tests pinning "P13 rejects RPS via CV"
- `tests/test_rps_p12_e2e.py` — 11 tests: RPS detected, GH/SIR/GoL rejected

Modified (3):
- `epc/orchestration.py` — added rps_spatial + P12 registrations
- `tests/test_orchestration.py` — count updates (12→13 models, 11→12 detectors,
  32→42 pairs) + 5 new Sprint 9 registration tests (27→32 total)
- `tests/test_cross_detection_matrix.py` — 9 new EXPECTED_OUTCOMES entries
  (RPS row + P12 column) + `test_sprint_9_pairs_covered` (15→16 total)
- `docs/detector_cards.md` — P12 card rewritten for the implemented detector

## Package Layout
```
epc/                                    # Core Python package (45 .py files)
├── detector_result.py                  # DetectorResult + DetectionTier + NullType
├── base_detector.py                    # detect() pipeline with tier/bonus logic
├── base_model.py                       # setup/step/run/get_metadata/get_timescale
├── base_metric.py                      # compute/required_keys/validate_history
├── orchestration.py                    # Substrate-aware detector dispatch
├── models/                             # 13 model implementations
│   ├── cell_view_sorting.py            # Zhang 1D (chimeric)
│   ├── cell_view_sorting_threaded.py   # non-deterministic threaded variant
│   ├── greenberg_hastings.py           # GH excitable CA
│   ├── game_of_life.py                 # Conway life
│   ├── vicsek.py, dorsogna_spp.py      # collective motion
│   ├── kuramoto.py                     # phase oscillators
│   ├── schelling.py                    # segregation
│   ├── btw_sandpile.py                 # Bak-Tang-Wiesenfeld
│   ├── nowak_may.py                    # spatial PD
│   ├── hegselmann_krause.py            # opinion dynamics
│   ├── sir_epidemic.py                 # lattice SIR CA
│   └── rps_spatial.py                  # ← Sprint 9: Reichenbach RPS
├── metrics/                            # 9 metric modules (22+ metrics)
└── detectors/                          # 12 detectors + P13/P15 discriminator
    ├── p1_aggregation.py               # P1 Moran's I + segregation
    ├── p5_flocking.py, p6_milling.py
    ├── p9_synchronization.py
    ├── p12_cyclic_dominance.py         # ← Sprint 9: neighbor-conditional ρ
    ├── p13_excitable_waves.py
    ├── p13_p15_discriminator.py
    ├── p14_soc.py
    ├── p15_fidelity_fix.py             # legacy GoL-specific
    ├── p15_persistent_computation.py   # generalized via step_fn (Sprint 8)
    ├── p21_polarization.py
    ├── p22_information_cascade.py
    ├── p27_spatial_reciprocity.py
    └── p31_delayed_gratification.py
tests/                                  # 21 test files, ~215 tests
docs/                                   # Catalog, detector cards, paper drafts
```

NOTE: Legacy `models/`, `metrics/`, `analysis/` at repo root are initial stubs
superseded by `epc/`. All real code lives under `epc/`.

## Key Design Decisions
- Models return full state history as list of dicts with numpy arrays
- Metrics operate on state histories (decoupled from execution)
- All randomness seeded for reproducibility
- DetectorResult is universal output; confidence is tier-capped
  (screening ≤ 0.60, confirmation ≤ 0.85, definitive ≤ 1.00)
- Boundary-conditioned TE (not raw average) for P13/P15 discrimination
- Substrate-aware orchestration prevents invalid model–detector pairings
- P15 generalized via step_fn + multi-checkpoint reproducibility + perturbative
  IC variations (Sprint 8; old `p15_fidelity_fix.py` retained as legacy smoke test)
- P12 (Sprint 9): neighbor-conditional replacement ratio ρ(X,Y) over best cyclic
  triple; spatial-shuffle null. `model_class` for RPS is `"cyclic_competition"`
  specifically to avoid the `"ca"`/`"excitable"` substring trap in P13's
  placeholder exclusion logic.
- See PROJECT_STATUS.md for Architecture Decisions Log (30+ entries)

## Model API Reference

```python
# lattice_2d models (all expose 'grid' + 'grid_dims' in state snapshots)
GreenbergHastings(rows, cols, n_states=3, threshold=1, neighborhood="moore",
    boundary="periodic", init_mode="random", init_density=0.3, seed=42)
GameOfLife(rows, cols, init_mode="random", init_density=0.3, seed=42)
run_schelling(grid_size, threshold, density, n_steps, seed)  # function, not class
NowakMayModel(rows=100, cols=100, b=1.8, init_mode="random",
    init_coop_fraction=0.5, boundary="periodic", seed=42)
SIREpidemicModel(rows=100, cols=100, infection_prob=0.3, recovery_prob=0.1,
    neighborhood="moore", init_mode="single_seed", seed=42)
RPSSpatialModel(rows=100, cols=100, mobility=1e-4, seed=42)  # OR exchange_rate=

# lattice_2d models without grid (sandpile: avalanche-level observables)
BTWSandpile(...)  # not compatible with grid-reading detectors

# continuous_2d / oscillator / opinion_space
VicsekModel(n_particles=300, box_size=7.0, speed=0.03, noise=0.1,
    interaction_radius=1.0, dt=1.0, init_mode="random", seed=None)
DOrsognaSPPModel(...)
KuramotoModel(...)
HegselmannKrauseModel(n_agents=50, epsilon=0.2, seed=42)
```

## Detector API Reference

```python
# Core contract: detector.detect(history, metadata) → DetectorResult
# All detectors accept n_permutations (or n_null_runs) and seed.

P1AggregationDetector(n_permutations=999)
P5FlockingDetector(n_permutations=199, seed=42)
P6MillingDetector(n_permutations=199, seed=42)
P9SynchronizationDetector(n_permutations=199, seed=42)
P12CyclicDominanceDetector(n_permutations=199, neighborhood="von_neumann", seed=42)
P13ExcitableWaveDetector(n_null_runs=199)
P13P15Discriminator(n_permutations=99)
P14SOCDetector(n_permutations=199, seed=42)
P15PersistentComputationDetector(step_fn=..., n_variations=8, seed=42)  # needs step_fn
P21PolarizationDetector(n_boot=1000, seed=42)
P22CascadeDetector(n_permutations=199, seed=42)
P27SpatialReciprocityDetector(n_permutations=199, seed=42)
P31DelayedGratificationDetector(n_permutations=199, seed=42)
```

## Statistical Power Requirements

| Context | Minimum |
|---------|---------|
| P1 label-shuffle null | n_permutations=999 |
| P5/P6/P9/P22/P27 permutation tests | n_permutations=199 (floor p=0.005) |
| P12 spatial-shuffle null | n_permutations=199 for confirmation; ≥499 for definitive (p < 0.005) |
| P13 null runs | n_null_runs=199 |
| P21 dip test bootstrap | n_boot=1000 |
| P31 non-redundancy CV | ≥500 runs for 10-fold CV |
| TE discriminator | n_permutations=99 (floor p=0.01) |
| KSG TE | n_permutations=99, T≥5000 |
| Kuramoto | N≥300 for finite-N convergence |
| P14 sandpile | ≥10,000 burn-in, ≥50,000 driving events |
| P15 generalized | n_variations ≥ 8 for definitive, ≥ 12 recommended |
| Any new permutation test | minimum 99 perms; for p<0.01: ≥199 perms; for p<0.005: ≥499 |

## Coding Standards
- Type hints on all function signatures
- Docstrings on all public methods
- NumPy for array operations; no heavy frameworks (no Mesa)
- Every test must verify the ACTUAL EFFECT, not just that code runs
- Honest assessment over moving fast
- Metadata strings matter: P13's placeholder exclusion logic checks for "ca"
  or "excitable" substrings in `model_class`. New lattice CA models should
  pick a `model_class` that avoids these substrings unless they genuinely ARE
  excitable dynamics (e.g., RPS uses "cyclic_competition", not "rps_ca").
