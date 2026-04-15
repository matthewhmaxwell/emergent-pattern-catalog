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

## Implementation Status

### Models (11 implemented, validated against published literature)

| Model | Cluster | Primary Patterns | Reference |
|-------|---------|-----------------|-----------|
| Zhang cell-view sorting | A, J | P1, P31 | Zhang et al. 2024 |
| Schelling segregation | A | P1 | Schelling 1971 |
| Vicsek model | B | P5 | Vicsek et al. 1995 |
| D'Orsogna SPP | B | P6 | D'Orsogna et al. 2006 |
| Kuramoto oscillators | C | P9 | Kuramoto 1975 |
| Greenberg-Hastings CA | D | P13 | Greenberg & Hastings 1978 |
| BTW sandpile | D | P14 | Bak, Tang & Wiesenfeld 1987 |
| Game of Life | E | P15 | Conway / Gardner 1970 |
| Hegselmann-Krause | F | P21 | Hegselmann & Krause 2002 |
| Nowak-May spatial PD | H | P27 | Nowak & May 1992 |

### Detectors (10 + discriminator, all with 3-tier detection)

P1 Aggregation, P5 Flocking, P6 Milling, P9 Synchronization,
P13 Excitable Waves, P14 SOC, P15 Persistent Computation,
P21 Polarization, P27 Spatial Reciprocity, P31 Delayed Gratification,
plus the P13/P15 boundary-conditioned TE discriminator.

### Test Suite: 101/101 passing

## Installation

```bash
git clone https://github.com/matthewhmaxwell/emergent-pattern-catalog.git
cd emergent-pattern-catalog
pip install -r requirements.txt
```

## Project Structure

```
emergent-pattern-catalog/
├── epc/             # Core package (models, metrics, detectors, orchestration)
├── tests/           # 13 test files, 101 tests
├── docs/            # Catalog, detector cards, ontology, paper draft
├── CLAUDE.md        # Claude Code context file
├── PROJECT_STATUS.md
└── REPLICATION_NOTES.md
```

## Ontological Dimensions (11)

Each pattern is classified along: spatial scale, temporal character, interaction
type, interaction substrate, agent homogeneity, goal structure, feedback structure,
memory, conflict structure, external driving, and update mode.

## References

- Levin, M. (2022). Technological Approach to Mind Everywhere. *Frontiers in Systems Neuroscience*.
- Zhang, T., Goldstein, A. & Levin, M. (2024). Classical sorting algorithms as a model of morphogenesis. *Adaptive Behavior*, 33, 25–54.
- Bak, P., Tang, C., & Wiesenfeld, K. (1987). Self-organized criticality. *Physical Review Letters*.
- Reynolds, C. W. (1987). Flocks, herds and schools. *SIGGRAPH*.
- Schelling, T. C. (1971). Dynamic models of segregation. *Journal of Mathematical Sociology*.
- Nowak, M. A. & May, R. M. (1992). Evolutionary games and spatial chaos. *Nature*.
- Kuramoto, Y. (1975). Self-entrainment of a population of coupled non-linear oscillators.
- Vicsek, T. et al. (1995). Novel type of phase transition in a system of self-driven particles. *Physical Review Letters*.
- D'Orsogna, M. R. et al. (2006). Self-propelled particles with soft-core interactions. *Physical Review Letters*.
- Hegselmann, R. & Krause, U. (2002). Opinion dynamics and bounded confidence models, analysis, and simulation. *JASSS*.

## License

MIT License. See [LICENSE](LICENSE) for details.
