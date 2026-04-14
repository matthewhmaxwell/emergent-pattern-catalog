# CLAUDE.md — Emergent Pattern Catalog

## Project Overview
Detection toolkit for emergent behavioral competencies in agent-based systems.
32 patterns across 10 clusters, organized in a three-layer architecture.

Repo: https://github.com/matthewhmaxwell/emergent-pattern-catalog

## Current State: Sprint 1 Complete

See `PROJECT_STATUS.md` for full tracker. See `REPLICATION_NOTES.md` for
Zhang et al. comparison.

## Code Structure

```
epc/                          # Python package
├── __init__.py
├── detector_result.py        # DetectorResult dataclass, tier-capped confidence
├── base_detector.py          # Abstract detector interface (all 32 detectors subclass this)
├── base_model.py             # Abstract model interface
├── base_metric.py            # Abstract metric interface
├── models/
│   ├── cell_view_sorting.py           # Sequential sorting model (fast, deterministic)
│   └── cell_view_sorting_threaded.py  # Threaded sorting model (faithful to Zhang)
├── metrics/
│   ├── aggregation.py                 # P1: Moran's I, segregation, clusters (1D + 2D)
│   └── delayed_gratification.py       # P31: Zhang's monotonicity-based DG metric
└── detectors/
    ├── p1_aggregation.py              # P1 detector with null models and exclusions
    └── p31_delayed_gratification.py   # P31 detector + NonRedundancyTest class
docs/
├── detector_cards.md          # v0.5.2 specs for all 32 detectors
├── pattern_catalog_v0_4.md    # Full catalog with 32 patterns, 10 clusters
├── ontology_v0_4.md           # 11 ontological dimensions
└── paper_outline_v0_4.md      # Paper structure
```

## Key Design Decisions
- Models return full state history as list of dicts
- Metrics decouple from execution — operate on state histories
- DetectorResult is the universal output schema (tier-capped confidence)
- All randomness seeded for reproducibility (except threaded model by design)
- Zhang sorting model: 1D array, N=100, matching their reference implementation
- DG metric verified bit-for-bit against Zhang's `avg_wandering_range()`
- Non-redundancy tests require ≥500 runs for reliable results
- No heavy frameworks — NumPy + SciPy + Matplotlib + scikit-learn

## Running the Code
```bash
pip install numpy scipy matplotlib scikit-learn
cd emergent-pattern-catalog
python -c "
from epc.models.cell_view_sorting import CellViewSorting
m = CellViewSorting(n=100, algorithm='bubble', seed=0)
h = m.run_to_completion()
print(f'Sorted in {m._swap_count} swaps, DG={m.compute_dg():.4f}')
"
```

## Key Results (Sprint 1)
- Zhang replication: swap counts within 4%, insertion DG within 4%
- P31 SURVIVES non-redundancy (ΔR²=+0.645, p<0.000001, 600 runs)
- Bubble DG ~35% higher than Zhang's — traced to Linux vs macOS threading

## What's Next (Sprint 2)
- Second model family (Boids, Greenberg-Hastings CA, or Schelling)
- 3-4 new detectors from different clusters
- Cross-model detection (run P1/P31 on new models)
- Transfer matrix construction begins
