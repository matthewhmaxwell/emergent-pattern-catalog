# Emergent Pattern Catalog

A systematic catalog, detection toolkit, and discovery platform for emergent
behavioral competencies in minimal agent-based systems. This project builds a
"periodic table" of emergent patterns — cataloging the surprising cognitive-like
behaviors that arise in systems of simple agents following local rules, and
providing computational tools to detect and quantify them across diverse model
families.

## Motivation

Research in Diverse Intelligence (Levin, 2022) has revealed that cognitive-like
competencies — memory, goal-directedness, self-organization, collective
decision-making — are not unique to neural systems but emerge across biological
scales from cells to swarms. Zhang et al. (2024) demonstrated that even minimal
sorting algorithms, when viewed from a cell-centric perspective, exhibit emergent
behaviors including spontaneous aggregation, delayed gratification, and
context-sensitive decision-making. This project extends that insight
systematically: if emergent competencies can arise in sorting arrays, where else
might they be hiding?

We implement a library of minimal agent-based models alongside a detection
toolkit of quantitative metrics for each cataloged pattern. The goal is to map
the landscape of emergent competencies, identify universal signatures, and
discover novel patterns through cross-model transfer analysis.

## Catalog Summary (v0.4)

The catalog contains **32 atomic patterns** organized into a three-layer architecture:

- **Layer 1 — Atomic patterns:** 32 detectable macro-behaviors across 10 clusters
- **Layer 2A — Mathematical descriptors:** Phase transitions, attractor dynamics, etc.
- **Layer 2B — Cognitive-analogue descriptors:** DI-framework interpretive annotations
- **Layer 3 — Meta-capacities:** Stigmergic coordination, multiscale competency, etc.

### Pattern Clusters

| Cluster | Patterns | Count |
|---------|----------|-------|
| A: Spatial organization | Aggregation, MIPS, Turing patterns, Territoriality | 4 |
| B: Collective motion | Flocking, Milling, Lane formation, Jamming | 4 |
| C: Temporal dynamics | Synchronization, Chimera states, Predator-prey, Cyclic dominance | 4 |
| D: Wave propagation | Excitable waves, Self-organized criticality | 2 |
| E: Information processing | Persistent computation, Associative memory, Distributed sensing | 3 |
| F: Decision-making | Consensus, Leadership, Quorum sensing, Polarization, Cascades, Anti-coordination | 6 |
| G: Resilience | Homeostasis, Canalized restoration, Stochastic resonance | 3 |
| H: Competition/cooperation | Spatial reciprocity, Wealth condensation | 2 |
| I: Structure formation | Trail/network formation, Autopoiesis | 2 |
| J: Agent-level competencies | Delayed gratification, Emergent specialization | 2 |

See [docs/pattern_catalog.md](docs/pattern_catalog.md) for the full catalog with
canonical models, detection metrics, dimensional profiles, and distinctness arguments.

## Installation

```bash
git clone https://github.com/matthewhmaxwell/emergent-pattern-catalog.git
cd emergent-pattern-catalog
pip install -r requirements.txt
```

## Project Structure

```
emergent-pattern-catalog/
├── models/          # Model implementations (sorting, Schelling, Boids, etc.)
├── metrics/         # Detection toolkit (one metric per cataloged pattern)
├── analysis/        # Orchestration: run models, apply metrics, generate reports
├── experiments/     # Experiment configs and results
├── docs/            # Living catalog (v0.4), ontology, paper outline
├── notebooks/       # Jupyter notebooks for exploration
├── tests/           # Unit tests
└── scripts/         # CLI utilities
```

## Ontological Dimensions (11)

Each pattern is classified along: spatial scale, temporal character, interaction
type, interaction substrate, agent homogeneity, goal structure, feedback structure,
memory, conflict structure, external driving, and update mode.

See [docs/ontology.md](docs/ontology.md) for full dimension definitions.

## References

- Levin, M. (2022). Technological Approach to Mind Everywhere. *Frontiers in Systems Neuroscience*.
- Zhang, J., et al. (2024). Classical sorting algorithms as a model of morphogenesis. *Adaptive Behavior*.
- Bak, P., Tang, C., & Wiesenfeld, K. (1987). Self-organized criticality. *Physical Review Letters*.
- Reynolds, C. W. (1987). Flocks, herds and schools. *SIGGRAPH*.
- Schelling, T. C. (1971). Dynamic models of segregation. *Journal of Mathematical Sociology*.
- Nowak, M. A. & May, R. M. (1992). Evolutionary games and spatial chaos. *Nature*.
- Kuramoto, Y. (1975). Self-entrainment of a population of coupled non-linear oscillators.
- Arthur, W. B. (1994). Inductive reasoning and bounded rationality. *American Economic Review*.

## License

MIT License. See [LICENSE](LICENSE) for details.
