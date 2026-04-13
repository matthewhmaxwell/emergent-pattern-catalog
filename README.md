# Emergent Pattern Catalog

A systematic catalog, detection toolkit, and discovery platform for emergent behavioral competencies in minimal agent-based systems. This project builds a "periodic table" of emergent patterns — cataloging the surprising cognitive-like behaviors that arise in systems of simple agents following local rules, and providing computational tools to detect and quantify them across diverse model families.

## Motivation

Research in Diverse Intelligence (Levin, 2022) has revealed that cognitive-like competencies — memory, goal-directedness, self-organization, collective decision-making — are not unique to neural systems but emerge across biological scales from cells to swarms. Zhang et al. (2024) demonstrated that even minimal sorting algorithms, when viewed from a cell-centric perspective, exhibit emergent behaviors including spontaneous aggregation, delayed gratification, and context-sensitive decision-making. This project extends that insight systematically: if emergent competencies can arise in sorting arrays, where else might they be hiding?

We implement a library of minimal agent-based models (sorting networks, Schelling segregation, Boids flocking, sandpile dynamics, cellular automata, Boolean gene regulatory networks) alongside a detection toolkit of quantitative metrics for each cataloged pattern. The goal is to map the landscape of emergent competencies, identify universal signatures, and discover novel patterns through cross-model transfer analysis.

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/emergent-pattern-catalog.git
cd emergent-pattern-catalog
pip install -r requirements.txt
```

## Quick Start

```python
from models.sorting.cell_view import CellViewBubbleSort
from metrics.aggregation import AggregationMetric

# Run a cell-view bubble sort
model = CellViewBubbleSort(size=100, seed=42)
model.setup()
history = model.run(max_steps=50000)

# Detect emergent aggregation (P1)
metric = AggregationMetric()
result = metric.compute(history)
print(f"Peak aggregation: {result['peak_aggregation']:.3f}")
print(f"Significant: {result['significant']}")
```

## Project Structure

```
emergent-pattern-catalog/
├── models/          # Model implementations (sorting, Schelling, Boids, etc.)
├── metrics/         # Detection toolkit (one metric per cataloged pattern)
├── analysis/        # Orchestration: run models, apply metrics, generate reports
├── experiments/     # Experiment configs and results (replication, transfer, discovery)
├── docs/            # Living catalog, ontology, paper outline
├── notebooks/       # Jupyter notebooks for exploration
├── tests/           # Unit tests
└── scripts/         # CLI utilities
```

## Pattern Catalog (In Progress)

| ID | Pattern | Canonical Model | Status |
|----|---------|-----------------|--------|
| P1 | Spontaneous Aggregation | Zhang sorting | Planned |
| P2 | Delayed Gratification | Zhang sorting | Planned |
| P3 | Self-Organized Criticality | Bak-Tang-Wiesenfeld | Planned |
| P4 | Collective Motion / Alignment | Reynolds Boids | Planned |
| P5 | Edge of Chaos Computation | Game of Life | Planned |
| P6 | Associative Memory | Hopfield / GRN | Planned |
| P7 | Homeostatic Regulation | Various | Planned |
| P8 | Stochastic Resonance | Various | Planned |
| P9 | Equifinality / Robustness | Various | Planned |

## Citation

```bibtex
@software{emergent_pattern_catalog,
  title = {Emergent Pattern Catalog: Detection Toolkit for Behavioral Competencies in Minimal Agent Systems},
  year = {2026},
  url = {https://github.com/YOUR_USERNAME/emergent-pattern-catalog}
}
```

## References

- Levin, M. (2022). Technological Approach to Mind Everywhere: An Experimentally-Grounded Framework for Understanding Diverse Bodies and Minds. *Frontiers in Systems Neuroscience*.
- Zhang, J., et al. (2024). Emergence of cognitive-like behaviors in simple sorting algorithms. *arXiv preprint*.
- Bak, P., Tang, C., & Wiesenfeld, K. (1987). Self-organized criticality. *Physical Review Letters*.
- Reynolds, C. W. (1987). Flocks, herds and schools: A distributed behavioral model. *SIGGRAPH*.
- Schelling, T. C. (1971). Dynamic models of segregation. *Journal of Mathematical Sociology*.

## License

MIT License. See [LICENSE](LICENSE) for details.
