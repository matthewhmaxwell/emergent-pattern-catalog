# Project Status Tracker

Last updated: 2026-04-13

## Current Phase: Sprint 1 — Foundation + Sorting Model Vertical Slice

### Sprint 1 Objectives
Build end-to-end pipeline: scaffolding → first model → first metrics → validation.
Proves the architecture works before scaling to full catalog.

### Sprint 1 Checklist

#### Phase 1 — Scaffolding ✅
- [x] Directory structure created
- [x] Project status tracker
- [x] `DetectorResult` dataclass with tier-capped confidence
- [x] Base detector interface
- [x] `BaseModel` abstract class
- [x] `BaseMetric` abstract class
- [x] Smoke tests passing

#### Phase 2 — First Model: Zhang Cell-View Sorting ✅ (rebuilt)
- [x] `epc/models/cell_view_sorting.py` — rebuilt to match Zhang et al. reference implementation
- [x] 1D array (N=100), cell-view with randomized sequential activation
- [x] Bubble, Insertion, Selection sort matching Zhang's cell logic
- [x] Frozen cells (movable/immovable), chimeric arrays (mixed algotypes)
- [x] Zhang's exact DG metric (monotonicity-based, per-swap trace)
- [x] Zhang's aggregation metric (same-algotype neighbor fraction)
- [x] Validated against paper: see REPLICATION_NOTES.md

#### Phase 3 — First Metrics ✅
- [x] `epc/metrics/aggregation.py` — P1 (Moran's I, segregation, clusters) — supports 1D and 2D
- [x] `epc/metrics/delayed_gratification.py` — P31 (Zhang's monotonicity-based DG metric)
- [x] `epc/detectors/p1_aggregation.py` — P1 detector with peak aggregation, label-shuffle null
- [x] `epc/detectors/p31_delayed_gratification.py` — P31 detector + non-redundancy test
- [x] P31 three-stage non-redundancy test — **P31 SURVIVES** (ΔR²=+0.645, p<0.000001)

#### Phase 4 — Validation against Zhang et al. ✅
- [x] Swap counts match within 4% for all algorithms
- [x] Insertion DG matches within 4% (1.06 vs 1.1)
- [x] Chimeric aggregation above 50% baseline confirmed
- [x] Threaded implementation built for frozen-cell sensitivity
- [x] DG metric verified exact match against Zhang's code
- [x] Bubble DG magnitude ~35% higher (Linux vs macOS thread scheduling)
- [x] Frozen-cell DG trend partially reproduced (0→2 increasing)

#### Key Findings (Sprint 1)

**Replication (Zhang et al.):**
- Faithful threading implementation + sequential approximation both provided
- Swap counts, insertion DG, chimeric aggregation all replicate
- Bubble DG offset traced to OS/Python environmental difference (documented)
- See REPLICATION_NOTES.md for detailed comparison

**P31 non-redundancy result (properly powered):**
- P31 SURVIVES non-redundancy (600 runs, 10-fold CV, p < 0.000001)
- Baseline (P1 features only): R² = −0.02 — aggregation cannot distinguish algorithms
- Extended (P1 + DG features): R² = 0.63 — DG explains 63% of outcome variance
- Ablation (P1 + shuffled DG): R² = −0.03 — temporal structure is the signal
- DG captures the sorting PROCESS (backtracking structure), while P1 captures
  the spatial OUTCOME (final clustering). These are independent signals.
- Earlier underpowered tests (48-120 runs) produced false negatives — lesson
  documented for future detector validation.

**Methodological lesson:**
- Non-redundancy tests require sufficient sample size (≥500 runs for 10-fold CV)
- Outcome variable must not be dominated by condition identity
- Baseline should NOT include condition labels unless testing within-condition
  variation specifically table

---

## Sprint 2 — Second Model Families + Cross-Model Detection

### Sprint 2 Checklist

#### Phase 1 — New Models ✅
- [x] `epc/models/greenberg_hastings.py` — Greenberg-Hastings excitable CA
  - 2D periodic grid, synchronous update, Moore neighborhood
  - States: rest/excited/refractory, threshold-based excitation
  - Init modes: random sparse excitation, broken_wave (spiral formation)
  - Validated: state counts, quiescence at high threshold, sustained spiral
- [x] `epc/models/vicsek.py` — Vicsek self-propelled particles
  - N agents in periodic 2D box, constant speed, alignment + noise
  - Optional attraction_strength for milling (D'Orsogna extension)
  - cKDTree periodic neighbor queries, circular mean for headings/COM
  - Validated: box containment, constant speed, low-noise→flocking, high-noise→disorder

#### Phase 2 — New Metrics ✅
- [x] `epc/metrics/excitable_waves.py` — P13 metrics
  - WavefrontSpeed: rest→excited tracking, COM displacement, speed CV
  - SpiralTipDetector: topological charge on 2x2 plaquettes (Barkley 1991)
  - WaveSourceCount: spontaneous excitation sites, persistent source tracking
- [x] `epc/metrics/collective_motion.py` — P5/P6 metrics
  - Polarization: φ = |mean(v̂_i)| order parameter
  - AngularMomentum: L = (1/N)Σ(r̂_i × v̂_i) with periodic COM
  - GroupSpeedRatio: R = |V_cm| / ⟨|v_i|⟩

#### Phase 3 — New Detectors ✅
- [x] `epc/detectors/p13_excitable_waves.py` — P13 excitable wave detector
  - Primary: wavefront speed persistence + CV
  - Null: frame-shuffle (mechanistic nulls as future work)
  - Exclusions: P15 stub, P12 metadata check
- [x] `epc/detectors/p5_flocking.py` — P5 flocking detector
  - Primary: polarization φ
  - Null: heading-shuffle (999 perms, φ_null ≈ 1/√N)
  - Exclusions: P6 (|L|), P7 (heading bimodality), P8 stub
- [x] `epc/detectors/p6_milling.py` — P6 milling detector
  - Primary: |L| angular momentum
  - Null: heading-shuffle
  - Secondaries: hollowness, rotation coherence, group speed ratio
  - Exclusions: P5 (R > 0.5), P12 metadata check

#### Phase 4 — Cross-Model Detection ✅
- [x] Transfer matrix tests: P13/P5/P6 on GH, Vicsek, sorting
- [x] P13 on GH (broken_wave): DETECTED
- [x] P13 on Vicsek: NOT detected (no grid)
- [x] P13 on quiescent GH: NOT detected
- [x] P5 on Vicsek (low noise): DETECTED
- [x] P5 on Vicsek (high noise): NOT detected
- [x] P5 on GH: NOT detected (no velocities)
- [x] P6 on GH: NOT detected
- [x] P6 on standard Vicsek: NOT detected or screening only

#### Scope Decisions (Sprint 2)
- **P9 deferred** to Sprint 2.5/3 — needs Kuramoto oscillators
- **P15 deferred** — needs Game of Life; GH CA is the negative control
- P15 exclusion in P13 detector: stub ("not_checked")
- P8 exclusion in P5 detector: stub ("not_checked")

## Sprint 3+ Plan (not started)
- Systematic model build-out (one per cluster minimum, ~15 total)
- Remaining 28 detectors
- Orchestration layer + transfer matrix
- Discovery analysis

---

## Full Model Inventory (17 planned)

| Model | Cluster | Primary Patterns | Status |
|-------|---------|-----------------|--------|
| Zhang cell-view sorting | A, J | P1, P31 | Sprint 1 |
| Schelling segregation | A | P1 | Sprint 2-3 |
| Active Brownian particles | A | P2 | Sprint 3+ |
| Gray-Scott reaction-diffusion | A | P3 | Sprint 3+ |
| Vicsek / Boids | B | P5, P6 | ✅ Sprint 2 |
| Nagel-Schreckenberg traffic CA | B | P8 | Sprint 3+ |
| Kuramoto oscillators | C | P9, P10 | Sprint 2-3 |
| Lotka-Volterra lattice | C | P11 | Sprint 3+ |
| Spatial rock-paper-scissors | C | P12 | Sprint 3+ |
| Greenberg-Hastings CA | D | P13 | ✅ Sprint 2 |
| BTW sandpile | D | P14 | Sprint 3+ |
| Game of Life | E | P15 | Sprint 2-3 |
| Hopfield network | E | P16 | Sprint 3+ |
| Hegselmann-Krause opinion | F | P21 | Sprint 3+ |
| SIR / threshold contagion | F | P22 | Sprint 3+ |
| Minority Game | F | P23 | Sprint 3+ |
| Nowak-May spatial PD | H | P27 | Sprint 3+ |
| Yard-Sale wealth model | H | P28 | Sprint 3+ |

## Full Detector Inventory (32 total)

| ID | Pattern | Status |
|----|---------|--------|
| P1 | Similarity-driven aggregation | Sprint 1 |
| P2 | MIPS | Sprint 3+ |
| P3 | Turing pattern formation | Sprint 3+ |
| P4 | Territoriality | Sprint 3+ |
| P5 | Flocking | ✅ Sprint 2 |
| P6 | Milling | ✅ Sprint 2 |
| P7 | Lane formation | Sprint 3+ |
| P8 | Jamming | Sprint 3+ |
| P9 | Synchronization | Sprint 2-3 |
| P10 | Chimera states | Sprint 3+ |
| P11 | Predator-prey oscillations | Sprint 3+ |
| P12 | Cyclic dominance | Sprint 3+ |
| P13 | Excitable waves | ✅ Sprint 2 |
| P14 | SOC | Sprint 3+ |
| P15 | Persistent computation | Sprint 2-3 (deferred) |
| P16 | Associative memory | Sprint 3+ |
| P17 | Distributed sensing | Sprint 3+ |
| P18 | Consensus | Sprint 3+ |
| P19 | Leadership | Sprint 3+ |
| P20 | Quorum sensing | Sprint 3+ |
| P21 | Polarization | Sprint 3+ |
| P22 | Information cascade | Sprint 3+ |
| P23 | Anti-coordination | Sprint 3+ |
| P24 | Homeostatic regulation | Sprint 3+ |
| P25 | Canalized restoration | Sprint 3+ |
| P26 | Stochastic resonance | Sprint 3+ |
| P27 | Spatial reciprocity | Sprint 3+ |
| P28 | Wealth condensation | Sprint 3+ |
| P29 | Trail formation | Sprint 3+ |
| P30 | Autopoiesis | Sprint 3+ |
| P31 | Delayed gratification | Sprint 1 |
| P32 | Emergent specialization | Sprint 3+ |

---

## Architecture Decisions Log

1. Models return full state history as list of dicts with numpy arrays
2. Metrics operate on state histories (decoupled from execution)
3. All randomness seeded for reproducibility
4. Results stored as JSON + numpy archives
5. No heavy frameworks — NumPy + Matplotlib + SciPy
6. DetectorResult is the universal output schema for all detectors
7. Confidence is tier-capped (screening ≤ 0.60, confirmation ≤ 0.85, definitive ≤ 1.00)
8. All persistence thresholds use system-intrinsic timescales, not fixed timesteps
9. Zhang sorting model uses 1D array (N=100) matching the paper, not 2D grid
10. Concurrency simulated via randomized sequential activation (see REPLICATION_NOTES.md)
11. Zhang's DG metric (monotonicity-based per-swap trace) used instead of distance-based
12. Reference implementation: https://github.com/Zhangtaining/cell_research
