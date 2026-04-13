# Emergent Pattern Catalog — Claude Code Context

## Project Summary
Systematic research program to catalog, detect, and discover emergent behavioral
competencies in minimal agent-based systems. Building a "periodic table" of
emergent patterns and a computational toolkit for detecting them.

## Current Status (v0.4)
- Pattern catalog: 32 atomic patterns across 10 clusters, 3 layers, 11 dimensions
- Architecture: Layer 1 (behaviors), Layer 2A (math descriptors), Layer 2B
  (cognitive-analogue annotations), Layer 3 (meta-capacities)
- Models: Base classes implemented (BaseModel, BaseMetric). All model and metric
  implementations are TODO stubs.
- See docs/pattern_catalog.md for the full catalog

## Architecture
- `models/`: Each model implements BaseModel with `setup()`, `step()`, `run()`, `get_state()`
- `metrics/`: Each metric implements BaseMetric with `compute(history)` returning named measurements
- `analysis/`: Orchestration layer — runs models, collects histories, applies metric batteries, generates reports
- `experiments/`: Config files (YAML) specifying model params, metric selections, trial counts

## Key Design Decisions
- Models return full state history as list of numpy arrays (or dicts), not just final state
- Metrics operate on state histories, not live simulations — decouples detection from execution
- All randomness seeded for reproducibility
- Results stored as JSON + numpy archives

## Next Priorities (v0.5)
1. Design detector specification cards for all 32 patterns (required observables,
   metrics, null models, pass thresholds, false positives, neighbor exclusions)
2. Implement Zhang et al. cell-view sorting (Bubble, Insertion, Selection) as first model
3. Implement aggregation metric (P1) and DG metric (P31)
4. Implement P31 three-stage non-redundancy test protocol
5. Replicate Zhang et al. Figure 8A (algotype clustering) as validation
6. Replicate Zhang et al. Figure 7 (DG vs frozen cells) as validation
7. Implement P13/P15 Transfer Entropy discriminator

## Pattern Clusters (32 total)
- A: Spatial organization (P1-P4)
- B: Collective motion (P5-P8)
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
- No heavy frameworks (no Mesa dependency — keep it lightweight)
- Tests for each model and metric
