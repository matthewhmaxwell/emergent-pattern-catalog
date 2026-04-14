# Project Status Tracker

Last updated: 2026-04-14

## Previous Phase: Sprint 1 — Foundation + Sorting Model Vertical Slice (COMPLETE)

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

## Current Phase: Sprint 2 — Second Model + Cross-Model Detection

### Sprint 2 Checklist

#### Phase 1 — Second Model: Greenberg-Hastings CA ✅
- [x] `epc/models/greenberg_hastings.py` — full GH excitable CA
- [x] 2D grid (100×100 default), κ states (3+), threshold θ, Moore/VN neighborhood
- [x] Periodic and fixed boundary conditions
- [x] 4 init modes: random, single_seed, broken_wave, custom
- [x] Vectorized (numpy) update — 500 steps on 100×100 in 0.37s
- [x] Scalar/vectorized verified identical output
- [x] Wave speed = exactly 1.00 cells/step for VN, θ=1 (matches theory)
- [x] Validated: random init → sustained waves, single seed → expanding ring + death,
      broken wave → persistent spirals, high threshold → quiescent

#### Phase 2 — P13 Excitable Wave Detector ✅
- [x] `epc/metrics/wave_propagation.py` — 4 metrics:
  - WavefrontSpeedLocal (inter-excitation interval, CV=0.030 for spirals)
  - WavefrontSpeed (centroid-based, noisy for multi-spiral)
  - SpiralTipDetector (topological charge, correctly finds 2 tips, net charge 0)
  - WavePersistence (active fraction, streak length, died-out flag)
- [x] `epc/detectors/p13_excitable_wave.py` — full P13 detector card implementation
- [x] P13 detects GH spirals at CONFIRMATION (confidence=0.70, p=0.005)
- [x] P13 correctly rejects: high-threshold GH (quiescent), single-seed GH (transient),
      Zhang sorting (no grid data)

#### Phase 3 — Cross-Model Detection ✅
- [x] Transfer matrix populated (7 models × 3 detectors + TE discriminator)
- [x] P1 on GH: screening only (wave autocorrelation ≠ aggregation, tier system correct)
- [x] P31 on GH: not detected (no sorting process)
- [x] P13 on sorting: not detected (graceful failure on missing grid data)
- [x] P1 type-diversity fix: single-type sorting now correctly returns not-detected
- [x] 6 canonical GH results replicated against literature (see REPLICATION_NOTES.md)

#### Key Findings (Sprint 2)

**Transfer Matrix (full power: P1=999 perms, P13=199 null, P31=deterministic, TE=99 perms):**

| Model          | P1 (aggregation)        | P13 (waves)             | P31 (DG)                | TE disc.              |
|----------------|-------------------------|-------------------------|-------------------------|-----------------------|
| GH spiral      | screening(0.60) p=.001  | confirmation(0.70) p=.005| —                       | P13 (1.0×) p=1.000   |
| GH random      | screening(0.60) p=.001  | confirmation(0.70) p=.005| —                       | P13 (2.1×) p=1.000   |
| GoL rpent      | confirm††(0.70) p=.001  | —                       | —                       | P15 cand (16.1×) p=.010|
| GoL random     | confirm††(0.70) p=.001  | —                       | —                       | P15 cand (15.1×) p=.010|
| Sort chimera   | screening(0.60) p=.001  | —                       | confirmation(0.60) p=.005| n/a                   |
| Sort bubble    | — p=1.000               | —                       | confirmation(0.60) p=.005| n/a                   |
| Sort insertion | — p=1.000               | —                       | confirmation(0.60) p=.005| n/a                   |

Grid models: 60×60, GH=600 steps, GoL=300 steps. Sorting: N=100, run to completion.

†† P1 on GoL: TRUE POSITIVE technically (Moran's I = 0.54/0.34, segregation = 0.91),
but CONCEPTUAL FALSE POSITIVE for "similarity-driven aggregation". B3/S23 dynamics
produce compact still-life/oscillator clusters without any attraction mechanism.
See Phase 4 findings.

**Diagonal dominance confirmed:** P13 strongly detects GH, P31 strongly detects
sorting, with clean off-diagonal negatives. P1 at screening on GH is real spatial
autocorrelation (Moran's I = 0.78 peak) — but the confirmation criteria correctly
filter it (segregation and CV thresholds not met by wave patterns).

**GH random vs spiral:** GH with random init self-organizes into spirals (Fisch,
Gravner & Griffeath 1991) but P13 only reaches screening (not confirmation) because
the multi-spiral system has shorter rotation counts per spiral. The detector
correctly requires ≥ 50 rotations for confirmation.

**GH validation (6 of 6 replicated):**
- Winding number theorem: defect → persists, no defect → dies ✅
- Wavefront speed: exactly 1.000 cells/step (VN and Moore, θ=1) ✅
- No dispersion: speed independent of κ (confirmed Gerhardt et al. 1990) ✅
- Topological charge conservation: 0/401 violations on torus ✅
- Spiral self-organization from random IC (when density sufficient) ✅
- Threshold dependence: θ=3 critical for Moore κ=3 ✅

**Bug found and fixed:** WavefrontSpeedLocal used wavelength = κ-1 instead of κ,
producing apparent dispersion that doesn't exist in basic GH. Discovered during
replication test 3. Speed CV (which P13 detector uses) was unaffected.
See REPLICATION_NOTES.md § Bug found and fixed.

#### Resolved Issues (Sprint 2)
- [x] P13 null model optimized: 72s → 37s via trajectory subsampling
- [x] P13 effect size fixed: now compares observed CV vs null CV distribution
      (was comparing speed mean vs CV null — wrong metric pair)
- [x] P1 exclusion for P13 evaluated: NOT needed per detector card spec.
      P13 is not in P1's exclusion list; the tier system handles it correctly.
- [x] Wavelength bug fixed in WavefrontSpeedLocal (κ-1 → κ)
- [x] P1 type-diversity fix: require ≥ 2 unique types for meaningful Moran's I.
      Single-type systems (bubble, insertion, selection) now correctly return
      not-detected instead of degenerate screening.

#### Phase 4 — P13/P15 Transfer Entropy Discriminator ✅
- [x] `epc/models/game_of_life.py` — Conway's GoL (B3/S23), verified:
  - Still lifes: block, beehive, loaf (all stable for 5+ generations)
  - Oscillators: blinker (period 2, exact cells), toad (period 2), pulsar (period 3)
  - Spaceships: glider c/4 (exact shift +1,+1 per period, verified 5 periods),
    LWSS c/2 (exact shift +2,0 per period, 9 cells, verified 5 periods)
  - R-pentomino: cell-by-cell match against independent B3/S23 implementation
    for first 6 generations. Population 116 at gen 1103 on 700×700 grid (exact
    match with published LifeWiki value).
  - Init modes: random, r_pentomino, glider_collision, lwss, custom
- [x] `epc/metrics/transfer_entropy.py` — plug-in TE estimator for discrete CAs.
  Average TE (TransferEntropy) and local TE map (LocalTransferEntropy).
  Following Lizier, Prokopenko & Zomaya (2007, 2012).
  Validated against 4 analytical ground truths:
  - Independent random: TE ≈ 0.000 ✅
  - Deterministic copy (y_{t+1} = x_neighbor): TE = 0.995 ≈ 1.0 bits ✅
  - XOR rule (y_{t+1} = y_t ⊕ x_neighbor): TE = 1.000 bits ✅
  - Noisy copy (p=0.8): TE = 0.532 bits (between 0 and 1) ✅
- [x] `epc/detectors/p13_p15_discriminator.py` — boundary-conditioned TE
  discriminator with GH control null and permutation test (Stage 1).
- [x] Validated across 10 seeds per model. Non-overlapping distributions.

**P13/P15 discriminator results (full power: 99 permutations, 60×60 grid):**

| Model             | Boundary TE | Classification | Ratio vs GH | p-value |
|-------------------|-------------|----------------|-------------|---------|
| GH spiral κ=5     | 0.000628    | P13            | 1.0×        | 1.000   |
| GH random κ=5     | 0.001350    | P13            | 2.1×        | 1.000   |
| GoL R-pentomino   | 0.010131    | P15 candidate  | 16.1×       | 0.010   |
| GoL random        | 0.009511    | P15 candidate  | 15.1×       | 0.010   |

**Separation: min GoL ratio (15.1×) >> max GH ratio (2.1×). Gap = 7.2×.**
**All p-values at permutation-appropriate levels (99-perm floor = 0.01).**

Additionally validated with 10 seeds at 40×40 (without permutation test):
GH 0.00128 ± 0.00000, GoL 0.00972 ± 0.00043. Non-overlapping distributions.

**Boundary TE ordering (validated against analytical controls):**
  self-predictable toggle (0.000) < GH waves (0.001) << GoL (0.010) < non-self-predictable (0.250)

**P1 on GoL finding (false positive):** P1 fires at confirmation on GoL
(Moran's I = 0.54 for R-pentomino, 0.34 for random). This is genuine spatial
autocorrelation from B3/S23 rule dynamics — still lifes and oscillators form
compact clusters. But NOT "similarity-driven aggregation" — no attraction
mechanism. The P1 detector cannot currently distinguish rule-driven structural
clustering from attraction-driven aggregation using state-history alone.
Documented as a P1 limitation.

**Key methodological insight:** Raw average TE does NOT discriminate P13 from P15
— GH actually has higher global TE than GoL because resting cells are maximally
responsive to any single excited neighbor. The discriminating signal is in
*boundary-conditioned* TE: TE computed only at cells with heterogeneous Moore
neighborhoods. At these state boundaries, GH transitions are nearly deterministic
given any neighbor state (TE ≈ 0), while GoL transitions depend on the specific
count of live neighbors (TE > 0). Theoretical basis: Schreiber (2000) — TE
measures information transfer BEYOND what the target's own history provides.

**Stage 2 (functional test) status:** Not yet implemented. Required for definitive
P15 classification per detector card. Stage 1 alone is sufficient for "P15 candidate"
flag, which is the appropriate Sprint 2 deliverable.

### Sprint 2 Open Items
1. P1 false positive on GoL: add temporal signature test (aggregation should
   increase monotonically for true P1; GoL autocorrelation fluctuates).
   Alternatively, add model-metadata check for binary CAs.
2. Stage 2 functional test for P15: do collision outcomes vary systematically
   with input configuration? Requires collision detection + outcome cataloguing.
3. TE discriminator compute cost: 99 perms on 60×60 takes 45-155s per model.
   Boundary TE inner loop could be further vectorized for production use.

---

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
| Vicsek / Boids | B | P5, P6 | Sprint 2 |
| Nagel-Schreckenberg traffic CA | B | P8 | Sprint 3+ |
| Kuramoto oscillators | C | P9, P10 | Sprint 2-3 |
| Lotka-Volterra lattice | C | P11 | Sprint 3+ |
| Spatial rock-paper-scissors | C | P12 | Sprint 3+ |
| Greenberg-Hastings CA | D | P13 | Sprint 2 ✅ |
| BTW sandpile | D | P14 | Sprint 3+ |
| Game of Life | E | P15 | Sprint 2 ✅ |
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
| P5 | Flocking | Sprint 2 |
| P6 | Milling | Sprint 2 |
| P7 | Lane formation | Sprint 3+ |
| P8 | Jamming | Sprint 3+ |
| P9 | Synchronization | Sprint 2-3 |
| P10 | Chimera states | Sprint 3+ |
| P11 | Predator-prey oscillations | Sprint 3+ |
| P12 | Cyclic dominance | Sprint 3+ |
| P13 | Excitable waves | Sprint 2 ✅ |
| P14 | SOC | Sprint 3+ |
| P15 | Persistent computation | Sprint 2 (Stage 1 TE disc. ✅, Stage 2 pending) |
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
13. GH model uses synchronous update (all cells update simultaneously) — standard for CAs
14. GH provides both scalar (for-loop) and vectorized (numpy roll) stepping, verified identical
15. Wavefront speed measured via local inter-excitation interval, not centroid displacement
    (centroid method fails for multi-spiral systems where centroids move erratically)
16. Spiral tips detected via topological charge on 2×2 plaquettes — standard method from
    Bray (2002) and others
17. Detectors must handle missing state keys gracefully (return not-detected, not crash)
    — required for cross-model detection where state formats differ
18. P1 detection on non-aggregation models (waves, etc.) is a valid screening result —
    spatial autocorrelation is real but the tier system correctly prevents false confirmation
19. P1 requires ≥ 2 unique types for screening — Moran's I is degenerate on constant-type
    arrays (0/0 → 0.0 trivially exceeds expected I). Added type-diversity guard.
20. GoL model uses numpy 2D convolution for neighbor counting — standard B3/S23 rule.
    Glider speed c/4 = 1 diagonal cell per 4 steps. T_prop = 4 × max(rows, cols).
21. P13/P15 boundary uses BOUNDARY-CONDITIONED Transfer Entropy, not raw average TE.
    Raw TE is dominated by the update rule's inherent responsiveness (GH > GoL because
    resting cells always respond to any excited neighbor). Conditioning on heterogeneous
    neighborhoods isolates the information-theoretic signal at wavefronts/structure edges.
22. TE estimator uses plug-in (frequency counting), appropriate for discrete CAs with
    small state alphabets (2-14 states). KSG estimator would be needed for continuous
    systems (Vicsek/Boids, Kuramoto) in later sprints.
23. P13/P15 discrimination threshold set at 3× boundary TE ratio vs GH control.
    Empirically, GH gives 0.9-1.0× and GoL gives 7.9-8.5× — clean separation.
