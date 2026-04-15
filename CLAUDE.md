# Emergent Pattern Catalog — Claude Code Context

## Project Summary
Systematic research program to catalog, detect, and discover emergent behavioral
competencies in minimal agent-based systems. Building a "periodic table" of
emergent patterns and a computational toolkit for detecting them.

Repository: https://github.com/matthewhmaxwell/emergent-pattern-catalog

## Current Status (Sprint 5 Complete — 101/101 tests passing)
- Pattern catalog: 32 atomic patterns across 10 clusters, 3 layers, 11 dimensions
- Models: 11 model files, 14 model variants — all validated against published literature
- Metrics: 9 modules, 22+ metrics — all verified against canonical results
- Detectors: 10 detectors + 1 discriminator — all with 3-tier detection + null models
- Orchestration: 5 substrate types, 11 models × 10 detectors → 24 compatible pairs
- Test suite: 101/101 passing (100 full-power, 1 reduced-power)
- Cluster coverage: A, B, C, D, E, F, H, J (missing: G Resilience, I Structure Formation)
- Paper Section 4 draft: 2,465 words with consolidated transfer matrix

## Package Layout
```
epc/                           # Core Python package
├── detector_result.py         # DetectorResult dataclass (tier-capped confidence)
├── base_detector.py           # Full detect() pipeline with tier/bonus logic
├── base_model.py              # setup/step/run/get_metadata/get_timescale
├── base_metric.py             # compute/required_keys/validate_history
├── orchestration.py           # Substrate-aware detector dispatch
├── models/                    # 11 model files
├── metrics/                   # 9 metric modules
└── detectors/                 # 11 detector files
tests/                         # 13 test files, 101 tests
docs/                          # Catalog, detector cards, paper draft
```

NOTE: Legacy `models/`, `metrics/`, `analysis/` at repo root are initial stubs
superseded by `epc/`. All real code lives under `epc/`.

## Key Design Decisions
- Models return full state history as list of dicts with numpy arrays
- Metrics operate on state histories (decoupled from execution)
- All randomness seeded for reproducibility
- DetectorResult is universal output; confidence is tier-capped
- Boundary-conditioned TE (not raw average) for P13/P15 discrimination
- Substrate-aware orchestration prevents invalid model–detector pairings
- See PROJECT_STATUS.md for Architecture Decisions Log (27 entries)

## Model API Reference (Sprint 5 verified)

```python
GreenbergHastings(rows, cols, n_states=3, threshold=1, neighborhood="moore",
    boundary="periodic", init_mode="random", init_density=0.3, seed=42)

VicsekModel(n_particles=300, box_size=7.0, speed=0.03, noise=0.1,
    interaction_radius=1.0, dt=1.0, init_mode="random", seed=None)

DOrsognaSPPModel(...)

KuramotoModel(...)

NowakMayModel(rows=100, cols=100, b=1.8, init_mode="random",
    init_coop_fraction=0.5, boundary="periodic", seed=42)

HegselmannKrauseModel(n_agents=50, epsilon=0.2, seed=42)
```

## Statistical Power Requirements

| Context | Minimum |
|---------|---------|
| P1 label-shuffle null | n_permutations=999 |
| P5/P6/P9/P27 permutation tests | n_permutations=199 (floor p=0.005) |
| P13 null runs | n_null_runs=199 |
| P21 dip test bootstrap | n_boot=1000 |
| P31 non-redundancy CV | ≥500 runs for 10-fold CV |
| TE discriminator | n_permutations=99 (floor p=0.01) |
| KSG TE | n_permutations=99, T≥5000 |
| Kuramoto | N≥300 for finite-N convergence |
| P14 sandpile | ≥10,000 burn-in, ≥50,000 driving events |
| Any new permutation test | minimum 99 perms; for p<0.01: ≥199 perms |

## Coding Standards
- Type hints on all function signatures
- Docstrings on all public methods
- NumPy for array operations; no heavy frameworks (no Mesa)
- Every test must verify the ACTUAL EFFECT, not just that code runs
- Honest assessment over moving fast
