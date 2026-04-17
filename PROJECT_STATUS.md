# Project Status Tracker

Last updated: 2026-04-15

## Current Phase: Sprint 5 Complete — All Carry-Forward Items Resolved

Sprints 1–5 complete. 101/101 tests passing (100 full-power, 1 reduced-power).
Sprint 5 verification session completed — all bugs found and fixed. Code on GitHub:
https://github.com/matthewhmaxwell/emergent-pattern-catalog

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

1. ~~Full-power TE benchmark at 60×60 with 99 perms~~ — **Done (Sprint 5)**
2. KSG TE on Vicsek/D'Orsogna — now unblocked by KSG implementation.
3. ~~Remaining TE discriminator limitation: ratio comparison at 25×25~~ — **Resolved (Sprint 5)**

---

## Sprint 5 — Complete ✅

### Completed

#### Full-Power TE Benchmark (60×60, 99 perms) ✅

- [x] `epc/metrics/transfer_entropy_global.py` — Optimized global boundary-conditioned TE.
      Vectorized boundary detection + np.add.at batch counting. Exact numerical agreement
      with P13P15Discriminator._boundary_te. 4× faster than per-cell approach.
- [x] `tests/test_te_benchmark_60x60.py` — 2 tests: numerical agreement + full benchmark.
- [x] Sprint 2 separation **CONFIRMED**: GoL/GH = 15.1× (random) and 16.1× (R-pentomino).
      GH spiral = 1.0×, GH random = 0.7×. All 4 classifications correct.
- [x] Total compute: 167s for all 4 models (60×60, 300 steps, 99 perms).

#### KSG TE on Continuous-Space Models (Vicsek + D'Orsogna) ✅

- [x] `tests/test_ksg_te_continuous.py` — 3 tests: Vicsek ordered/disordered,
      D'Orsogna milling/free, summary.
- [x] KSG TE detects coupling in ordered Vicsek (3/5 neighbor pairs sig, p=0.020)
      and correctly shows no coupling in disordered state (0/5, p=0.306).
- [x] KSG TE detects coupling in D'Orsogna milling (5/5 pairs sig, p=0.020)
      and correctly shows no coupling in free particles (0/5, p=1.000).
- [x] Confirms KSG TE extends information-transfer measurement from lattice CAs
      (plug-in TE, Sprint 2) to continuous-space SPP models.

Results:

| System | State | Median TE | Sig pairs | Median p |
|---|---|---|---|---|
| Vicsek (η=0.5) | Ordered (φ=0.986) | +0.030 | 3/5 | 0.020 |
| Vicsek (η=5.0) | Disordered (φ=0.159) | -0.035 | 0/5 | 0.306 |
| D'Orsogna | Milling (\|L\|=2.29) | -0.130 | 5/5 | 0.020 |
| D'Orsogna | Free (no interaction) | -6.318 | 0/5 | 1.000 |

Note: Negative absolute TE values reflect KSG finite-sample bias, not negative
information transfer. The permutation test is unaffected (both observed and null
share the bias). Statistical significance is the correct measure.

#### BTW 1/f Noise Measurement Fix ✅

- [x] Root cause: `compute_spectral_beta` was measuring PSD of avalanche sizes
      (approximately IID, white noise) instead of total energy E(t) = Σh(x,t).
- [x] Fix: `btw_sandpile.py` now tracks `energy_history` and `dissipation_history`
      per driving event. `detect_p14` accepts optional `energy` parameter.
- [x] `compute_spectral_beta` prefers detrended energy signal when available,
      falls back to cumulative activity method (documented as unreliable for BTW).
- [x] Result: β = 1.41 (was -0.17). Now correctly in 1/f range [0.5, 1.5].
- [x] P14 detector still DEFINITIVE (confidence 0.850, τ=1.248, γ=0.641).
- [x] Sprint 4 1/f caveat **RESOLVED**.
- [x] Updated `test_sandpile_p14_e2e.py` to pass energy signal.

#### Nowak-May Spatial PD × P27 End-to-End ✅

- [x] `epc/models/nowak_may.py` — Nowak & May (1992) spatial prisoner's dilemma.
      L×L lattice, binary C/D strategies, synchronous imitation update,
      Moore neighborhood, deterministic highest-payoff copying.
- [x] `epc/detectors/p27_spatial_reciprocity.py` — P27 detector with 3-tier
      detection, Moran's I permutation test (199 perms), PD structure verification,
      P1 exclusion via has_movement check.
- [x] `tests/test_nowak_may_p27_e2e.py` — 3 tests: physics validation (4 b-values),
      P27 detection (positive + 2 negative controls), transfer matrix row.
- [x] Model validated: b=1.0 → all C, b=1.5 → C survives (87%), b=1.8 → C in
      clusters (41%, I=0.48), b=2.0 → C extinct.
- [x] P27 DEFINITIVE at b=1.8 (f_C=0.408, I=0.497, p=0.005, PD verified).
      Negative controls: b=2.0 → none, b=1.0 → not definitive (no dilemma).
- [x] Opens Cluster H (competition/cooperation). 13th model, 10th detector.

#### Hegselmann-Krause × P21 Polarization ✅

- [x] `epc/models/hegselmann_krause.py` — Hegselmann & Krause (2002) bounded-confidence
      opinion dynamics. N agents, continuous opinions in [0,1], synchronous averaging
      within confidence bound ε. Converges in O(10) steps.
- [x] `epc/detectors/p21_polarization.py` — P21 detector with 3-tier detection,
      Hartigan's dip test (bootstrap, 1000 samples), cluster counting, persistence
      check, unimodal IC verification, P18 exclusion.
- [x] `tests/test_hk_p21_e2e.py` — 2 tests: physics (5 ε-values), P21 detection
      (positive ε=0.2 and ε=0.1, negative ε=0.5).
- [x] Model validated: ε=0.5 → consensus (1 cluster), ε=0.2 → polarization (2 clusters),
      ε=0.1 → fragmentation (4 clusters), ε=0.05 → 7 clusters.
- [x] P21 DEFINITIVE at ε=0.2 (2 clusters, dip p=0.001, from unimodal IC).
      Negative: ε=0.5 → none (consensus). ε=0.1 also DEFINITIVE (4 clusters).
- [x] Opens Cluster F (Decision-Making). 14th model, 11th detector.

Results:

| Model | Boundary TE | Ratio vs GH | Sprint 2 | Classification |
|---|---|---|---|---|
| GH spiral (κ=5) | 0.000628 | 1.0× | 1.0× | P13 ✅ |
| GH random (κ=3) | 0.000444 | 0.7× | 2.1× | P13 ✅ |
| GoL random | 0.009511 | 15.1× | 15.1× | P15_candidate ✅ |
| GoL R-pentomino | 0.010131 | 16.1× | 16.1× | P15_candidate ✅ |

This resolves Sprint 4 candidate #1 and #3. The 25×25 ratio compression (5.45 vs 4.77)
was a grid-size artifact; at 60×60 the clean 15-16× separation reproduces exactly.

### Remaining Sprint 5 Candidates

1. ~~KSG TE on Vicsek/D'Orsogna~~ — **Done**
2. ~~BTW 1/f noise measurement fix~~ — **Done**
3. ~~Nowak-May spatial PD × P27~~ — **Done**
4. ~~Hegselmann-Krause × P21~~ — **Done**
5. Additional models: Active Brownian (P2), Gray-Scott (P3), Nagel-Schreckenberg (P8),
   Lotka-Volterra (P11), spatial RPS (P12), Hopfield (P16),
   SIR (P22), Minority Game (P23), Yard-Sale (P28).
6. Paper Sections 3, 5, 6, 7 drafting.

---

## Full Model Inventory (17 planned, 14 implemented)

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
| Hegselmann-Krause opinion | F | P21 | ✅ Sprint 5 |
| SIR / threshold contagion | F | P22 | Sprint 4+ |
| Minority Game | F | P23 | Sprint 4+ |
| Nowak-May spatial PD | H | P27 | ✅ Sprint 5 |
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
| P21 | Polarization | ✅ Sprint 5 |
| P22 | Information cascade | Sprint 4+ |
| P23 | Anti-coordination | Sprint 4+ |
| P24 | Homeostatic regulation | Sprint 4+ |
| P25 | Canalized restoration | Sprint 4+ |
| P26 | Stochastic resonance | Sprint 4+ |
| P27 | Spatial reciprocity | ✅ Sprint 5 |
| P28 | Wealth condensation | Sprint 4+ |
| P29 | Trail formation | Sprint 4+ |
| P30 | Autopoiesis | Sprint 4+ |
| P31 | Delayed gratification | ✅ Sprint 1 |
| P32 | Emergent specialization | Sprint 4+ |

## Full Metric Inventory (21 implemented across 8 modules)

| Module | Metrics | Cluster |
|--------|---------|---------|
| `aggregation.py` | Moran's I, Segregation Index, Cluster Stats (1D + 2D) | A |
| `aggregation_convergence.py` | Convergence gain, Spearman trend, plateau stability, type constancy | A (guard) |
| `delayed_gratification.py` | Zhang's monotonicity-based DG | J |
| `wave_propagation.py` | WavefrontSpeedLocal, WavefrontSpeed, SpiralTipDetector, WavePersistence | D |
| `transfer_entropy.py` | TransferEntropy (plug-in), LocalTransferEntropy | E (lattice) |
| `transfer_entropy_vectorized.py` | Vectorized boundary-conditioned TE | E (lattice) |
| `transfer_entropy_global.py` | Global aggregate boundary TE (optimized) | E (lattice) |
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
25. Substrate-aware detector dispatch. 5 substrate types (lattice_1d, lattice_2d,
    continuous_2d, oscillator, opinion_space). Transfer matrix is block-diagonal
    by substrate type. 11 models × 10 detectors → 24/110 compatible pairs.
26. Global aggregate boundary TE (transfer_entropy_global.py). Accumulates frequency
    tables across all boundary cells per timestep using np.add.at, then computes
    single global TE. Matches discriminator._boundary_te exactly (verified 1e-6).
    4× faster than per-cell averaged approach at 60×60. Sprint 5 benchmark confirms
    GoL/GH separation: 15.1× and 16.1× at full power.
27. BTW 1/f noise uses total energy E(t) = Σh(x,t), not avalanche size sequence.
    Avalanche sizes are approximately IID (white noise spectrum). Energy fluctuations
    show 1/f scaling with β ≈ 1.4. compute_spectral_beta now prefers detrended
    energy signal; btw_sandpile tracks energy_history per driving event.

---

## Sprints 6–9 — Summary Block (appended retroactively in Sprint 9)

This block captures sprints 6–9 in a compact form. Sprint 5's "Current Phase"
heading near the top of this file was not updated per-sprint; the canonical
status now lives in CLAUDE.md and this block.

### Sprint 6 (commit `4ff238c`) — Cross-Detection Matrix

Added 8 cross-detection tests covering every compatible model-detector pair
that wasn't explicitly audited in Sprints 1–5. Key findings: D'Orsogna × P5
and Vicsek × P6 correctly mutually reject; GoL × P13 had been a false positive
at screening because the `n_states < 3` guard was a warning, not a hard guard
(fix: promoted to hard guard, GoL × P13 now rejects via the guard). Added
Schelling state enrichment: state dicts now include `grid_dims`, `n_states`,
`step` to support cross-detector compatibility.

Paper sections 3, 5, 6, 7 drafted.

### Sprint 7 (commit `36306e8`) — SIR Epidemic + P22 Information Cascade

Added `epc/models/sir_epidemic.py` (lattice SIR CA) and
`epc/detectors/p22_information_cascade.py` (information cascade via Moran's I
on infection-time map). Replication targets from Datta & Acharyya (2021):
bell-shaped epidemic curve, linear wavefront radius growth, phase transition
at p_c dependent on q and neighborhood. SIR × P22 → DEFINITIVE
(Moran I_time = 0.987, Cohen's d = 109.5, p = 0.005, cascade reach = 1.0).

### Sprint 8 (commit `d9b0465`) — SIR Reference Correction + P15 Generalization

- SIR primary reference corrected from Fuks & Lawniczak 2002 (LGCA with
  moving individuals — architecture mismatch) to Datta & Acharyya 2021
  (same fixed-site CA with independent-neighbor infection).
- P15 generalized from GoL-specific (`p15_fidelity_fix.py`) to a step_fn-based
  detector (`p15_persistent_computation.py`) that works on any deterministic
  lattice_2d CA via multi-checkpoint reproducibility + perturbative IC
  variations. Legacy `p15_fidelity_fix.py` retained as smoke test.
- GH and GoL replication documentation expanded in REPLICATION_NOTES.md.

### Sprint 8b (commit `1820d28`) — Consolidation

Four fixes:
1. P22 P13-exclusion logic rewritten — old fixed-last-quarter window mishandled
   late-peaking epidemics; new logic uses final state + trailing 2×T_spread
   window.
2. Cross-detection test coverage expanded (+6 tests in test_sir_p22_e2e.py).
3. `EXPECTED_OUTCOMES` table in `tests/test_cross_detection_matrix.py` grew
   from 8 → 18 audited pairs.
4. P15 test imports modernized to use the generalized detector.

### Sprint 9 — RPS (Reichenbach 2007) + P12 Cyclic Dominance

**Scientific findings:**

1. **P13 does NOT false-positive on RPS.** The Sprint 9 prompt predicted
   that P13 might incorrectly fire on RPS because RPS satisfies the `n_states
   ≥ 3` hard guard and produces persistent spiral-like wavefronts. In practice,
   P13 rejects RPS cleanly at screening across 3 mobilities × 3 seeds via
   `wavefront_speed_cv ∈ [0.59, 0.68]`, well above the 0.2 screening
   threshold. Mechanism: GH wavefronts are clock-driven (uniform speed,
   CV ≈ 0.05–0.15); RPS wavefronts are stochastic-neighbor-driven (high CV).
   The RPS/GH CV ratio exceeds 3× robustly. Pinned in
   `tests/test_rps_p13_boundary.py`.

2. **P12 complements P13 via a neighbor-conditional replacement ratio**:
   `ρ(X, Y) = P(cell→Y | had Y-neighbor) / P(cell→Y | no Y-neighbor)`.
   For RPS dominance edges, ρ ≈ 70–200 at canonical coexistence parameters.
   For GH clock transitions, ρ = 1.0 exactly. Two orders of magnitude of
   separation. Primary metric:
   `intransitivity_score = log10(max over cyclic triples (min forward ρ))`.
   Null model: global spatial shuffle of the grid at each timestep.

**Detection results (implemented in `P12CyclicDominanceDetector`):**
- RPS (L=30, M=1e-4, 80 gens, n_perm=199): **CONFIRMATION**, score=1.83,
  ρ_min=67.6, p=0.005, Cohen's d > 200, coexistence=1.0, direction stable,
  P13 and P22 both excluded.
- RPS (L=40, M=1e-5, 150 gens, n_perm=499): **DEFINITIVE**, score=2.06,
  p=0.002, confidence=0.85.
- GH n=3 and n=5: rejected, score = 0.0 exactly (ρ = 1.0 for clock-driven).
- SIR: rejected (no cyclic dominance), score = 0.0.
- GoL and Nowak-May: rejected via `n_candidate_species < 3` prerequisite.

**Bug found during P12 development:** `_best_intransitive_triple` initially
used the cached `self._all_species` set from `_validate_prerequisites`, which
made direction-stability secondaries crash on sub-trajectories. Fix: helper
now computes candidate species locally from the passed history. Without this
fix, DEFINITIVE tier was unreachable.

### Inventory totals at end of Sprint 9
- 45 epc .py files (was 43 at Sprint 8b HEAD)
- 21 test files (was 18)
- 42 compatible model×detector pairs (was 32)
- 27 audited cross-detection pairs (was 18)
- 12 detectors, 13 models, 5 substrate types

### Architecture decisions added in Sprint 9

28. (carried from Sprint 7) P22 infection-time Moran's I — preferred over
    final-state Moran because once an epidemic reaches most cells, the final
    state is near-uniform (I → 0) while the infection-time map retains the
    full wavefront structure.
29. (carried from Sprint 8) P15 generalized via step_fn + multi-checkpoint
    reproducibility + perturbative IC variations. Legacy p15_fidelity_fix.py
    retained for backward compat.
30. (Sprint 9) P12 primary metric is the neighbor-conditional replacement
    ratio ρ(X,Y) over the best cyclic 3-triple. ρ is smoothed with
    ε = 1e-6 to avoid div-by-zero. Null model is global spatial shuffle
    at each timestep (destroys neighbor-transition correlation while
    preserving species marginals).
31. (Sprint 9) `model_class` strings matter: P13's placeholder exclusion
    logic checks for "ca" or "excitable" substrings. RPS uses
    `model_class="cyclic_competition"` explicitly to avoid this trap.
    Future lattice CA models should do the same unless they genuinely ARE
    clock-driven excitable dynamics.

---

### Sprint 10 — P1 primary metric: peak → final-state Moran's I

**Motivation:** The Sprint 9 carry-forward list flagged SIR × P1 (and after
Sprint 9, RPS × P1) as ambiguous classifications. Both models fired P1 at
screening under the old `max(peak, final)` primary metric, but their final
states looked very different. Sprint 10 ran an empirical characterization
across six canonical models before touching the detector, per the
replication philosophy.

**Empirical characterization (6 models):**

| Model                 | I_peak | I_final | I_sustained | sustained_CV | seg_final |
|-----------------------|--------|---------|-------------|--------------|-----------|
| Schelling τ=0.375     | +0.414 | +0.414  | +0.414      | 0.00         | +0.650    |
| Nowak-May b=1.8       | +0.794 | +0.530  | +0.500      | 0.07         | +0.776    |
| **SIR β=0.20 γ=0.3**  | +0.892 | **+0.019** | +0.175   | **0.99**     | +0.981    |
| **RPS M=1e-4**        | +0.582 | **+0.550** | +0.562   | 0.016        | +0.668    |
| GH n=8 random         | +0.412 | +0.204  | +0.204      | 0.00         | +0.307    |
| Random grid (noise)   | +0.028 | +0.015  | −0.007      | inf          | +0.348    |

**Scientific finding:** SIR and RPS do **not** share the same dynamics.
SIR's spatial clustering is a genuinely transient wavefront (peak 0.89 →
final 0.02, sustained CV=0.99). RPS's spatial clustering is genuinely
sustained: spiral domains rotate but persist, so final ≈ peak ≈ 0.55.
A global peak→sustained swap (what the Sprint 10 prompt initially
recommended) would not have distinguished these. A peak→final swap does.

**Detector changes:**

1. `epc/detectors/p1_aggregation.py::_compute_primary` — primary
   `morans_i` is now `morans_i_final`. The peak and sustained values are
   retained as diagnostic fields.
2. `epc/detectors/p1_aggregation.py::_check_screening` — added a 0.05
   magnitude floor on top of the trivial `I > expected_I` check. This
   corrects a pre-existing bug where a pure-noise random grid passed
   screening at p ≈ 0.03 through chance sampling variance.

**Outcome on the 6-model panel:**
- Schelling, Nowak-May: CONFIRMATION (unchanged).
- RPS: screening (unchanged in tier, but now correctly justified by
  sustained clustering rather than transient).
- SIR: **REJECTED** (was screening). Transfer matrix entry flips from `S`
  to `rej`.
- GH: screening → blocked at confirmation by secondary gates (seg=0.31
  < 0.4). Unchanged.
- Random grid: **REJECTED** (was passing screening through sampling
  variance). Bugfix.

**Test updates:**
- `TestSIRP1Screening` → `TestSIRP1Rejection`: flip assertion, verify
  peak/final gap as diagnostic.
- `TestNowakMayP1CoOccurrence::test_nowak_may_p1_confirmation`: threshold
  relaxed 0.5 → 0.4 (final-I = 0.49 vs peak-I = 0.89 on 100×100 b=1.8);
  added explicit `peak > final` assertion.
- `EXPECTED_OUTCOMES[("sir_epidemic", "P1")]`: `screening` → `rejected`.
- RPS × P1 comment updated: "transient" → "sustained spiral clustering".
- **New** `TestRPSP1ScreeningLevel` (2 tests in `test_rps_p12_e2e.py`):
  - `test_rps_p1_screening_under_new_primary` — RPS still screens.
  - `test_rps_vs_sir_p1_asymmetry` — cross-model assertion that both models
    peak ≈ 0.5+ but only RPS has final ≈ peak.

### Inventory totals at end of Sprint 10
- 45 epc .py files (unchanged)
- 21 test files (unchanged)
- 27 audited cross-detection pairs (unchanged count; SIR × P1 reclassified)
- 12 detectors, 13 models, 5 substrate types (unchanged)
- Baseline 109 Sprint 9 tests + 2 new Sprint 10 tests = 111 fast tests pass

### Architecture decisions added in Sprint 10

32. (Sprint 10) **P1 primary metric is final-state Moran's I**, not peak
    over trajectory. The peak and sustained values are reported as
    diagnostic fields (`morans_i_peak`, `morans_i_sustained`). Rationale:
    the 6-model characterization showed SIR has peak ≫ final (transient
    wavefront) while RPS has peak ≈ final (sustained spiral domains). A
    peak-based primary falsely conflated these. A final-based primary
    correctly distinguishes them. See REPLICATION_NOTES.md for the full
    empirical table.

33. (Sprint 10) **P1 screening requires final I ≥ 0.05**, not merely
    `I > expected_I`. The trivial check passes pure random noise via
    sampling variance (random grid gave I=0.015 > expected=-0.0002, null
    p=0.03). The 0.05 floor rejects this without affecting any of the
    canonical positives (Schelling 0.41, NM 0.53, RPS 0.55 all clear the
    floor comfortably).
