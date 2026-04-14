# Project Status Tracker

Last updated: 2026-04-14

## Current Phase: Sprint 3 — Collective Motion (Vicsek/D'Orsogna)

---

## Sprint 1 — Foundation + Sorting Model ✅

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

---

## Sprint 2 — Excitable Media + Computation ✅

#### Models ✅
- [x] `epc/models/greenberg_hastings.py` — Excitable CA. κ states, θ threshold,
      Moore/VN neighborhoods, periodic/fixed BC. Vectorized numpy update. 4 init modes.
- [x] `epc/models/game_of_life.py` — Conway B3/S23. 5 init modes (random, r_pentomino,
      glider_collision, lwss, custom).

#### Metrics ✅
- [x] `epc/metrics/wave_propagation.py` — WavefrontSpeedLocal, WavefrontSpeed (centroid),
      SpiralTipDetector (topological charge), WavePersistence
- [x] `epc/metrics/transfer_entropy.py` — Plug-in TE estimator (TransferEntropy,
      LocalTransferEntropy). Validated against 4 analytical ground truths.

#### Detectors ✅
- [x] `epc/detectors/p13_excitable_wave.py` — P13 with wavefront speed, spiral tips,
      shuffle null, P15/P12 exclusion hooks
- [x] `epc/detectors/p13_p15_discriminator.py` — Boundary-conditioned TE discriminator
      with GH control null and permutation test

#### Key Findings (Sprint 2)

**GH replication:** 6/6 canonical results replicated (winding number, wave speed 1.0,
no dispersion, topological charge conservation, spiral self-organization, threshold dependence).

**GoL verification:** B3/S23 cell-by-cell match against independent reference. Still lifes,
oscillators, glider c/4, LWSS c/2, R-pentomino gen 1103 pop 116 exact match with LifeWiki.

**Full-power transfer matrix (25 cells):**

| Model          | P1              | P13             | P31             | TE disc.            |
|----------------|-----------------|-----------------|-----------------|---------------------|
| GH spiral      | screen p=.001   | CONFIRM p=.005  | —               | P13 (1.0×) p=1.000  |
| GH random      | screen p=.001   | CONFIRM p=.005  | —               | P13 (2.1×) p=1.000  |
| GoL rpent      | CONFIRM p=.001  | —               | —               | P15cand (16.1×) p=.01|
| GoL random     | CONFIRM p=.001  | —               | —               | P15cand (15.1×) p=.01|
| Sort chimera   | screen p=.001   | —               | CONFIRM p=.005  | n/a                 |
| Sort bubble    | —               | —               | CONFIRM p=.005  | n/a                 |
| Sort insertion | —               | —               | CONFIRM p=.005  | n/a                 |

**P1 on GoL:** Documented CONCEPTUAL FALSE POSITIVE — genuine spatial autocorrelation
from B3/S23 rule dynamics, but no attraction mechanism.

**P13/P15 TE discriminator:** Boundary-conditioned TE cleanly separates GH (TE ≈ 0.001,
ratio 1-2×) from GoL (TE ≈ 0.010, ratio 15-16×). Novel methodological contribution.

---

## Sprint 3 — Collective Motion (In Progress)

#### Models ✅
- [x] `epc/models/vicsek.py` — Standard Vicsek Model (1995). N point particles, 2D
      periodic box, constant speed v₀, arctan2 circular mean heading update, angular
      noise η. Vectorized sparse-matrix neighbor averaging via cKDTree. 4 init modes.
      1.2 ms/step at N=300. All 6 published claims replicated.
- [x] `epc/models/dorsogna_spp.py` — D'Orsogna SPP with Morse potential (2006).
      Second-order Newtonian dynamics, Rayleigh self-propulsion (α−β|v|²)v, pairwise
      Morse attraction-repulsion. RK4 integration. Open boundary. Published milling
      parameters replicated: |L|=0.996, φ=0.008, ring density confirmed.

#### Metrics ✅
- [x] `epc/metrics/collective_motion.py` — 5 metrics for Cluster B:
      PolarizationMetric (φ), GroupSpeedRatioMetric (R), AngularMomentumMetric (L
      with periodic-aware COM), HeadingAutocorrelationMetric, HeadingDistributionMetric
      (nematic order parameter for antiparallel detection)

#### Detectors ✅
- [x] `epc/detectors/p5_flocking.py` — P5 detector with 3-tier detection,
      heading-shuffle null (199 perms), P6 exclusion (|L|>0.5), P7 exclusion
      (nematic order), confidence scoring per detector card
- [x] `epc/detectors/p6_milling.py` — P6 detector with 3-tier detection,
      heading-shuffle null, ring density profile (hollowness), P5 exclusion (R>0.5)

#### Cross-Detection Transfer Matrix ✅
- [x] P5 × Vicsek = **DEFINITIVE** (φ=0.988, p=0.005)
- [x] P5 × D'Orsogna = not detected (φ=0.008)
- [x] P6 × D'Orsogna = **DEFINITIVE** (|L|=0.996, hollow=0.000)
- [x] P6 × Vicsek = not detected (|L|=0.017)
- [x] Perfect discrimination: no false positives, no missed detections

#### Sprint 3 Remaining
- [ ] Run P1, P13, P31, TE detectors on Vicsek and D'Orsogna → extend transfer matrix
- [ ] Consider Kuramoto oscillators for Cluster C (P9 synchronization)
- [ ] Begin orchestration layer for automated model×detector sweeps
- [ ] P1 false positive mitigation on GoL (carry-forward from Sprint 2)
- [ ] Stage 2 functional test for definitive P15 (carry-forward)

---

## Full Model Inventory (17 planned)

| Model | Cluster | Primary Patterns | Status |
|-------|---------|-----------------|--------|
| Zhang cell-view sorting | A, J | P1, P31 | ✅ Sprint 1 |
| Schelling segregation | A | P1 | Sprint 3+ |
| Active Brownian particles | A | P2 | Sprint 3+ |
| Gray-Scott reaction-diffusion | A | P3 | Sprint 3+ |
| Vicsek model | B | P5 | ✅ Sprint 3 |
| D'Orsogna SPP | B | P6 | ✅ Sprint 3 |
| Nagel-Schreckenberg traffic CA | B | P8 | Sprint 3+ |
| Kuramoto oscillators | C | P9, P10 | Sprint 3 |
| Lotka-Volterra lattice | C | P11 | Sprint 3+ |
| Spatial rock-paper-scissors | C | P12 | Sprint 3+ |
| Greenberg-Hastings CA | D | P13 | ✅ Sprint 2 |
| BTW sandpile | D | P14 | Sprint 3+ |
| Game of Life | E | P15 | ✅ Sprint 2 |
| Hopfield network | E | P16 | Sprint 3+ |
| Hegselmann-Krause opinion | F | P21 | Sprint 3+ |
| SIR / threshold contagion | F | P22 | Sprint 3+ |
| Minority Game | F | P23 | Sprint 3+ |
| Nowak-May spatial PD | H | P27 | Sprint 3+ |
| Yard-Sale wealth model | H | P28 | Sprint 3+ |

## Full Detector Inventory (32 total)

| ID | Pattern | Status |
|----|---------|--------|
| P1 | Similarity-driven aggregation | ✅ Sprint 1 |
| P2 | MIPS | Sprint 3+ |
| P3 | Turing pattern formation | Sprint 3+ |
| P4 | Territoriality | Sprint 3+ |
| P5 | Flocking | ✅ Sprint 3 |
| P6 | Milling | ✅ Sprint 3 |
| P7 | Lane formation | Sprint 3+ |
| P8 | Jamming | Sprint 3+ |
| P9 | Synchronization | Sprint 3 |
| P10 | Chimera states | Sprint 3+ |
| P11 | Predator-prey oscillations | Sprint 3+ |
| P12 | Cyclic dominance | Sprint 3+ |
| P13 | Excitable waves | ✅ Sprint 2 |
| P14 | SOC | Sprint 3+ |
| P15 | Persistent computation | Sprint 2 (screening; needs Stage 2) |
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
| P31 | Delayed gratification | ✅ Sprint 1 |
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
13. GH model uses vectorized numpy update with periodic/fixed BC
14. GoL verified cell-by-cell against LifeWiki reference
15. Transfer entropy plug-in estimator validated against 4 analytical ground truths
16. Boundary-conditioned TE (not raw average) is the correct discriminator for P13/P15
17. P1 on GoL is a documented conceptual false positive (rule-driven, not attraction-driven)
18. Sprint 3: Continuous-space models use state dict with `positions` (N,2) and `velocities` (N,2)
19. Sprint 3: Vicsek uses cKDTree(boxsize=L) for periodic neighbor finding + sparse matrix averaging
20. Sprint 3: Heading-shuffle null uses random uniform headings (not permutation of existing)
    because permuting near-identical headings in an ordered flock preserves φ
21. Sprint 3: Nematic order parameter S=|⟨exp(2iθ)⟩| for antiparallel (P7 lane) detection
22. Sprint 3: D'Orsogna model uses open boundary (no periodic BC) — attraction-repulsion
    keeps particles cohesive. T_cross measured from actual group diameter, not initial spread.
23. Continuous-space TE will need KSG estimator (not plug-in) — deferred to Sprint 3+
