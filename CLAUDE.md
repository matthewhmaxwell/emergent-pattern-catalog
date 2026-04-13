# Emergent Pattern Catalog — Claude Code Context

## Project Summary
Systematic research program to catalog, detect, and discover emergent behavioral 
competencies in minimal agent-based systems. Building a "periodic table" of emergent 
patterns and a computational toolkit for detecting them.

## Architecture
- `models/`: Each model implements a BaseModel interface with `setup()`, `step()`, `run()`, and `get_state()`
- `metrics/`: Each metric implements a BaseMetric interface with `compute(history)` returning a dict of named measurements
- `analysis/`: Orchestration layer that runs models, collects state histories, applies metric batteries, generates reports
- `experiments/`: Config files (YAML) specifying model params, metric selections, and trial counts

## Key Design Decisions
- Models return full state history as list of numpy arrays (or dicts), not just final state
- Metrics operate on state histories, not live simulations — decouples detection from execution
- All randomness seeded for reproducibility
- Results stored as JSON + numpy archives

## Current Priorities
1. Implement base classes (BaseModel, BaseMetric)
2. Implement Zhang et al. cell-view sorting (Bubble, Insertion, Selection) as first model
3. Implement aggregation metric (P1) and DG metric (P2)
4. Replicate Zhang et al. Figure 8A (algotype clustering) as validation
5. Replicate Zhang et al. Figure 7 (DG vs frozen cells) as validation

## Coding Standards
- Type hints on all function signatures
- Docstrings on all public methods
- NumPy for array operations
- Matplotlib for plotting
- No heavy frameworks (no Mesa dependency — keep it lightweight)
- Tests for each model and metric
