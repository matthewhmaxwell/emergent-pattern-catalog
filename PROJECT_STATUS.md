# Project Status Tracker

Last updated: 2026-04-14

## Current Phase: Sprint 4 — Planning

Sprints 1–3 complete. All carry-forward items resolved. 33/33 tests passing.
Code on GitHub: https://github.com/matthewhmaxwell/emergent-pattern-catalog

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

P1=999 perms | P13=199 null | P31=deterministic | TE=99 perms

| Model          | P1              | P13             | P31             | TE disc.            |
|----------------|-----------------|-----------------|-----------------|---------------------|
| GH spiral      | screen p=.001   | CONFIRM p=.005  | —               | P13 (1.0×) p=1.000  |
| GH random      | screen p=.001   | CONFIRM p=.005  | —               | P13 (2.1×) p=1.000  |
| GoL rpent      | **REJECT**      | —               | —               | P15cand (16.1×) p=.01|
| GoL random     | **REJECT**      | —               | —               | P15cand (15.1×) p=.01|
| Sort chimera   | screen p=.001   | —               | CONFIRM p=.005  | n/a                 |
| Sort bubble    | —               | —               | CONFIRM p=.005  | n/a                 |
| Sort insertion | —               | —               | CONFIRM p=.005  | n/a                 |

GoL P1 originally reported as CONFIRM in Sprint 2 — changed to **REJECT** after
temporal convergence guard + type constancy check (Sprint 3 carry-forward fix).

**P13/P15 TE discriminator:** Boundary-conditioned TE cleanly separates GH (TE ≈ 0.001,
ratio 1-2×) from GoL (TE ≈ 0.010, ratio 15-16×). Raw average TE gives the WRONG
ordering (GH > GoL). Novel methodological contribution.

---

## Sprint 3 — Collective Motion, Synchronization, Carry-Forward Resolutions ✅

### Phase 1 — Collective Motion Models ✅

#### Models
- [x] `epc/models/vicsek.py` — Standard Vicsek Model (1995). N point particles, 2D
      periodic box, constant speed v₀, arctan2 circular mean heading update, angular
      noise η. Vectorized sparse-matrix neighbor averaging via cKDTree(boxsize=L).
      4 init modes. 1.2 ms/step at N=300. All 6 published claims replicated.
- [x] `epc/models/dorsogna_spp.py` — D'Orsogna SPP with Morse potential (2006).
      Second-order Newtonian dynamics, Rayleigh self-propulsion (α−β|v|²)v, pairwise
      Morse attraction-repulsion U(r) = −Cₐ exp(−r/lₐ) + Cᵣ exp(−r/lᵣ). RK4
      integration. Open boundary. 8.5 ms/step at N=100.

#### Metrics
- [x] `epc/metrics/collective_motion.py` — 5 metrics for Cluster B:
      PolarizationMetric (φ), GroupSpeedRatioMetric (R), AngularMomentumMetric (L
      with periodic-aware COM), HeadingAutocorrelationMetric, HeadingDistributionMetric
      (nematic order parameter S=|⟨exp(2iθ)⟩| for antiparallel detection)

#### Detectors
- [x] `epc/detectors/p5_flocking.py` — P5 with 3-tier detection, heading-shuffle null
      (199 perms, random uniform headings), P6 exclusion (|L|>0.5), P7 exclusion
      (nematic order S high + polar φ low), confidence scoring per detector card
- [x] `epc/detectors/p6_milling.py` — P6 with 3-tier detection, heading-shuffle null,
      ring density profile (hollowness metric), P5 exclusion (R>0.5)

#### Cross-Detection Transfer Matrix
| | Vicsek (η=0.5) | D'Orsogna (milling) |
|---|---|---|
| P5 (flocking) | ✓ DEFINITIVE (φ=0.988, p=0.005) | ✗ none (φ=0.008) |
| P6 (milling) | ✗ none (|L|=0.017) | ✓ DEFINITIVE (|L|=0.996, hollow=0.000) |

Perfect discrimination: no false positives, no missed detections.

### Phase 2 — Kuramoto Oscillators + P9 Synchronization ✅

#### Model
- [x] `epc/models/kuramoto.py` — Standard all-to-all Kuramoto coupled oscillators
      (1975). Lorentzian/Gaussian frequency distributions. Mean-field O(N)
      reformulation, RK4 integration. Phase transition scan utility.

#### Replication (6/6 published claims)
- Disordered baseline: r = 0.071 = 1/√200 (ratio 1.00 exact)
- Phase transition onset at K=1.0 for N≥300, analytical K_c=2γ=1.0
- Ordered r(K): r(1.2)=0.448 vs analytical 0.408 (error 0.039),
  r(2.0)=0.707 vs 0.707 (error 0.000), r(4.0)=0.879 vs 0.866 (error 0.013)
- Frequency entrainment: 97% locked at K=6K_c, σ_inst=0.010 (locked only)

#### Detector
- [x] `epc/detectors/p9_synchronization.py` — P9 with 3-tier detection, uncoupled null
      (199 perms, random-phase), P10 exclusion (uniform local r check).
      Definitive on Kuramoto at K=8K_c (r=0.963, 119 T_osc, p=0.005).
      None below K_c (r=0.087, p=0.185).

### Phase 3 — KSG Transfer Entropy ✅

- [x] `epc/metrics/transfer_entropy_ksg.py` — KSG (Kraskov-Stögbauer-Grassberger) TE
      estimator for continuous variables. Frenzel & Pompe CMI extension. k-NN with
      Chebyshev distance, phase embedding for circular variables. Resolves
      architecture decision #23 (plug-in TE destroys information for continuous state).

#### Validation
- MI error 0.013 nats on Gaussian (ρ=0.8)
- Independent TE = 0.002 (p=0.39) — correctly non-significant
- Coupled TE = 0.304 (p=0.01) — correctly significant
- Kuramoto sync: TE increases with coupling, p=0.02 at K=2K_c and K=6K_c

### Phase 4 — Carry-Forward Resolutions ✅

**P1 false positive on GoL — RESOLVED.** Two-layer defense:
1. Type constancy check: GoL alive/dead states change every step → not constant
   type labels → P1 structurally inapplicable.
2. Orchestration substrate check: GoL has no 'types' observable → N/A.
Verified: Schelling correctly passes (I_initial=0.005→I_late=0.415), GoL
correctly rejected, synthetic convergent correctly passes.

- [x] `epc/metrics/aggregation_convergence.py` — Temporal convergence guard for P1.
      Moran's I trajectory, Spearman trend, plateau stability, convergence gain
      (I_initial → I_late). Type constancy flag for GoL rejection.
- [x] `epc/models/schelling.py` — Minimal Schelling segregation (1971). 2D grid,
      two types, threshold-based satisfaction. Canonical P1 positive control.
      Converges from I=0.005 to I=0.415.

**P15 definitive — RESOLVED.** Two-stage discriminator:
- Stage 1 (TE): GoL ratio 15-16× vs GH 1-2× (Sprint 2 full-power).
- Stage 2 (functional): Deterministic replay shows reproducibility=1.000
  with 2 distinct coarse outcome types (annihilation + still_life).
  Glider collisions at different phases → different qualitative outcomes.

- [x] `epc/detectors/p15_functional_test.py` — Stage 2: structure tracking, collision
      detection, outcome cataloguing. 414 collisions, 410 distinct outcomes (fine-grained).
- [x] `epc/detectors/p15_fidelity_fix.py` — Fixed fidelity: deterministic replay
      (reproducibility=1.000) + coarse outcome classification (2 distinct types:
      annihilation + still_life). Meets detector card spec.

**TE discriminator compute cost — RESOLVED.**

- [x] `epc/metrics/transfer_entropy_vectorized.py` — Fully vectorized boundary-conditioned
      TE. Replaces per-cell Python loop with batch numpy. 2.7× speedup at test scale,
      exact numerical agreement. Inner loop O(n_states³) not O(N_cells).

### Phase 5 — Orchestration ✅

- [x] `epc/orchestration.py` — Substrate-aware detector dispatch. 4 substrate types
      (lattice_1d, lattice_2d, continuous_2d, oscillator). 11 models × 8 detectors
      registered. 24/88 compatible pairs, 64 substrate mismatches correctly identified.

### Sprint 3 Test Suite (33/33 pass)

| Test file | Tests | Coverage |
|---|---|---|
| test_sprint2_carryforward.py | 11/11 | P1 guard, P15 collisions, TE vectorization |
| test_kuramoto.py | 9/9 | Model replication, P9 detector, cross-detection |
| test_orchestration_ksg.py | 7/7 | Compatibility matrix, KSG MI/TE validation |
| test_caveat_resolution.py | 6/6 | K_c convergence, P15 fidelity, P1 on Schelling, GoL rejection, KSG Kuramoto scaling, TE discriminator direction |

---

## Sprint 4 — In Progress

### Completed

#### Schelling × P1 End-to-End ✅
- [x] Format adapter: Schelling `{'grid'}` → P1-compatible `{'grid', 'grid_dims'}`
- [x] Full P1 detector on Schelling (50×50, density=0.9, 999 permutations)
- [x] Result: **CONFIRMATION** tier (p=0.001 floor, confidence=0.700)
- [x] Moran's I = 0.423 (vs null mean -0.0002), Cohen's d = 49.87
- [x] Segregation index = 0.652, sustained_i_cv = 0.000 (perfect plateau)
- [x] Temporal guard passes: I_initial=0.005 → I_late=0.415, is_plateaued=True
- [x] Negative controls: random grid not confirmed (p=0.452), GoL rejected (types not constant)
- [x] Tier analysis: CONFIRMATION is correct — 999 shuffle perms floor at p=0.001,
      definitive requires p < 0.001. Would need ≥1999 perms or mechanistic null.

#### BTW Sandpile × P14 SOC ✅
- [x] `epc/models/btw_sandpile.py` — BTW sandpile (Bak, Tang & Wiesenfeld 1987).
      2D square lattice L×L, z_c=4, open boundary, vectorized parallel toppling.
      Includes dissipative variant (p_diss bulk grain loss) for subcritical null.
- [x] `epc/detectors/p14_soc.py` — P14 detector. Clauset et al. (2009) power-law
      MLE via powerlaw package, LR tests vs exponential/lognormal, duration scaling,
      3-tier detection with dissipative null comparison.
- [x] BTW physics validated: critical state (max z=3), mean height 2.098,
      4.3 decades size span, heavy-tailed (mean/median=12.1)
- [x] Result: **DEFINITIVE** tier (confidence=0.850)
- [x] τ = 1.247 (MLE) / 1.241 (log-binned) vs published τ ≈ 1.20 (error 0.047)
- [x] Power-law strongly preferred over exponential (R=80.6, p<0.001)
- [x] Duration scaling γ = 0.642 (T ~ s^0.64)
- [x] Dissipative null correctly exponential (R=-6.0, max size 68 vs 20,972)
- [x] Dissipative sandpile correctly NOT detected as SOC
- [x] Caveat: log-normal preferred over simple power-law (R=-76.2) — documented
      property of 2D BTW multifractal scaling, not a detection failure.
- [x] Caveat: 1/f noise not detected (β=-0.17) — measurement methodology needs work.
      Detector correctly uses duration scaling as alternative secondary.

#### Paper Section 4 Draft ✅
- [x] `paper_section4_draft.md` — 2,465 words covering all 7 model families
- [x] Subsections: Zhang sorting (P1/P31), GH (P13), GoL (P15), Vicsek/D'Orsogna
      (P5/P6), Kuramoto (P9), Schelling (P1), BTW sandpile (P14)
- [x] Consolidated 10-row × 8-column transfer matrix
- [x] Methodological lessons section (power, test correctness, boundary conditioning,
      intrinsic timescales)

### Remaining Candidates

1. Full-power TE benchmark at 60×60 with 99 perms — confirm vectorized speedup
   reproduces the GH 1-2× vs GoL 15-16× separation from Sprint 2.
2. KSG TE on Vicsek/D'Orsogna — now unblocked by KSG implementation.
3. Remaining TE discriminator limitation: ratio comparison (GoL 5.45 vs GH 4.77)
   at test scale (25×25) is directionally correct but does not yet show the clean
   15× vs 2× separation from Sprint 2 full-power runs. Needs 60×60 verification.

---

## Full Model Inventory (17 planned, 12 implemented)

| Model | Cluster | Primary Patterns | Status |
|-------|---------|-----------------|--------|
| Zhang cell-view sorting | A, J | P1, P31 | ✅ Sprint 1 |
| Zhang cell-view sorting (threaded) | A, J | P1, P31 | ✅ Sprint 1 |
| Schelling segregation | A | P1 | ✅ Sprint 3 |
| Active Brownian particles | A | P2 | Sprint 4+ |
| Gray-Scott reaction-diffusion | A | P3 | Sprint 4+ |
| Vicsek model | B | P5 | ✅ Sprint 3 |
| D'Orsogna SPP | B | P6 | ✅ Sprint 3 |
| Nagel-Schreckenberg traffic CA | B | P8 | Sprint 4+ |
| Kuramoto oscillators | C | P9, P10 | ✅ Sprint 3 |
| Lotka-Volterra lattice | C | P11 | Sprint 4+ |
| Spatial rock-paper-scissors | C | P12 | Sprint 4+ |
| Greenberg-Hastings CA | D | P13 | ✅ Sprint 2 |
| BTW sandpile | D | P14 | ✅ Sprint 4 |
| Game of Life | E | P15 | ✅ Sprint 2 |
| Hopfield network | E | P16 | Sprint 4+ |
| Hegselmann-Krause opinion | F | P21 | Sprint 4+ |
| SIR / threshold contagion | F | P22 | Sprint 4+ |
| Minority Game | F | P23 | Sprint 4+ |
| Nowak-May spatial PD | H | P27 | Sprint 4+ |
| Yard-Sale wealth model | H | P28 | Sprint 4+ |

## Full Detector Inventory (32 total, 9 implemented + discriminator)

| ID | Pattern | Status |
|----|---------|--------|
| P1 | Similarity-driven aggregation | ✅ Sprint 1 (+ temporal guard Sprint 3) |
| P2 | MIPS | Sprint 4+ |
| P3 | Turing pattern formation | Sprint 4+ |
| P4 | Territoriality | Sprint 4+ |
| P5 | Flocking | ✅ Sprint 3 |
| P6 | Milling | ✅ Sprint 3 |
| P7 | Lane formation | Sprint 4+ |
| P8 | Jamming | Sprint 4+ |
| P9 | Synchronization | ✅ Sprint 3 |
| P10 | Chimera states | Sprint 4+ |
| P11 | Predator-prey oscillations | Sprint 4+ |
| P12 | Cyclic dominance | Sprint 4+ |
| P13 | Excitable waves | ✅ Sprint 2 |
| P14 | SOC | ✅ Sprint 4 |
| P15 | Persistent computation | ✅ Sprint 2–3 (TE Stage 1 + fidelity Stage 2) |
| P16 | Associative memory | Sprint 4+ |
| P17 | Distributed sensing | Sprint 4+ |
| P18 | Consensus | Sprint 4+ |
| P19 | Leadership | Sprint 4+ |
| P20 | Quorum sensing | Sprint 4+ |
| P21 | Polarization | Sprint 4+ |
| P22 | Information cascade | Sprint 4+ |
| P23 | Anti-coordination | Sprint 4+ |
| P24 | Homeostatic regulation | Sprint 4+ |
| P25 | Canalized restoration | Sprint 4+ |
| P26 | Stochastic resonance | Sprint 4+ |
| P27 | Spatial reciprocity | Sprint 4+ |
| P28 | Wealth condensation | Sprint 4+ |
| P29 | Trail formation | Sprint 4+ |
| P30 | Autopoiesis | Sprint 4+ |
| P31 | Delayed gratification | ✅ Sprint 1 |
| P32 | Emergent specialization | Sprint 4+ |

## Full Metric Inventory (20 implemented across 7 modules)

| Module | Metrics | Cluster |
|--------|---------|---------|
| `aggregation.py` | Moran's I, Segregation Index, Cluster Stats (1D + 2D) | A |
| `aggregation_convergence.py` | Convergence gain, Spearman trend, plateau stability, type constancy | A (guard) |
| `delayed_gratification.py` | Zhang's monotonicity-based DG | J |
| `wave_propagation.py` | WavefrontSpeedLocal, WavefrontSpeed, SpiralTipDetector, WavePersistence | D |
| `transfer_entropy.py` | TransferEntropy (plug-in), LocalTransferEntropy | E (lattice) |
| `transfer_entropy_vectorized.py` | Vectorized boundary-conditioned TE | E (lattice) |
| `transfer_entropy_ksg.py` | KSG TE for continuous variables | E (continuous) |
| `collective_motion.py` | Polarization φ, GroupSpeedRatio R, AngularMomentum L, HeadingAutocorrelation, HeadingDistribution (nematic S) | B |

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
16. Boundary-conditioned TE (not raw average) is the correct discriminator for P13/P15.
    Raw average gives WRONG ordering (GH > GoL).
17. P1 on GoL is a documented conceptual false positive — rule-driven spatial
    autocorrelation, not attraction-driven aggregation. Resolved in Sprint 3 via
    type constancy check + temporal convergence guard.
18. Continuous-space models use state dict with `positions` (N,2) and `velocities` (N,2)
19. Vicsek uses cKDTree(boxsize=L) for periodic neighbor finding + sparse matrix averaging
20. Heading-shuffle null for P5/P6 uses random uniform headings (NOT permutation of
    existing headings). Permuting near-identical headings in an ordered flock preserves
    φ ≈ 1 and gives uninformative p-values.
21. Nematic order parameter S = |⟨exp(2iθ)⟩| for P7 (lane) exclusion in P5.
    Detects antiparallel alignment where naive approach fails.
22. D'Orsogna T_cross from measured group diameter (not init_radius — mill compacts
    from R=5.0 to diameter 3.02).
23. KSG estimator for continuous-space TE. Plug-in TE destroys information via
    discretization of continuous state. RESOLVED: Frenzel & Pompe CMI extension,
    k-NN with Chebyshev distance, phase embedding for circular variables. Validated
    on Gaussian MI (error 0.013 nats), coupled vs independent TE, Kuramoto scaling.
24. P1 temporal convergence guard with type constancy check. Two conditions:
    (a) types_are_constant=True (persistent agent labels, not dynamic state),
    (b) has_gain AND (is_monotonic OR is_plateaued) using I_initial from timestep 0.
    GoL correctly rejected (types not constant). Schelling correctly passes
    (ΔI=0.41, monotonic convergence).
25. Substrate-aware detector dispatch. 4 substrate types (lattice_1d, lattice_2d,
    continuous_2d, oscillator). Transfer matrix is block-diagonal by substrate type.
    11 models × 8 detectors → 24/88 compatible pairs, 64 substrate mismatches.
