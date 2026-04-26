# Project Status Tracker

Last updated: Sprint 20 (April 2026)

## Current state (post Sprint 20)

The catalog has progressed through 20 sprints of work. Sprint 1 built
the foundation (Zhang sorting + P1 + P31). Sprints 2–11 filled out the
core model inventory and detectors across every substrate type.
Sprint 12 was a paper catch-up sprint. Sprint 13 added the first
continuous-field model (Gray-Scott + P3) and the first continuous-field
substrate type. Sprints 14, 14.5, and 14.6 were cleanup sprints.
Sprint 15 added Nagel-Schreckenberg traffic + P8 — the second
`lattice_1d` model and the first with integer-valued velocity
observables. Sprint 16 added Active Brownian Particles + P2 — the
third `continuous_2d` model and the first detector whose DEFINITIVE
tier depends on a metadata-based mechanistic-null gate (Decision 43),
completing the three-class discrimination framework. Sprint 17 added
Yard-Sale + P28 (wealth condensation) — the first well-mixed
(non-spatial) agent population in the registry, introducing the new
`scalar_wealth` substrate. Sprint 18 added the non-local Kuramoto
ring + P10 (chimera states) — the second oscillator-substrate model,
creating the first 2×2 within-substrate block in the transfer matrix.
Sprint 19 was a cleanup sprint that closed the Sprint 11
orchestration-registration gap for Lotka-Volterra and P11, brought
§4/§5/§6/§7 of the paper draft up to date through §4.19, and added
finite-size slow tests at L = 250–1000 for NS/YS/KN. **Sprint 20**
adds the voter model + P18 (coarsening to consensus) — expanding the
lattice_2d-with-grid block to nine models and seven compatible
detectors, the catalog's most heavily populated cross-detection
block. Voter is the canonical microscopic substrate for emergent
consensus and the simplest member of a coarsening-without-surface-
tension universality class (Dornic et al. 2001). P18 uses a
permutation-null statistical test (ADR 54), an early-window wall-
density Spearman secondary metric (ADR 55), and canonical async
voter dynamics exclusively (ADR 56).

## Inventory snapshot at Sprint 20 HEAD

- **19 model families** (20 model files, Zhang sequential and threaded
  count as the same family)
- **19 registered detectors** in `epc/orchestration.py::DETECTOR_REGISTRY`
- **7 substrate types**: lattice_1d, lattice_2d, lattice_2d_continuous,
  continuous_2d, oscillator, opinion_space, scalar_wealth
- **18 of 32 patterns covered** (P1, P2, P3, P5, P6, P8, P9, P10, P11,
  P12, P13, P14, P15, P18, P21, P22, P27, P28, P31)
- **173 audited cells** in `tests/test_cross_detection_matrix.py::
  EXPECTED_OUTCOMES` (Sprint 18: 146 → Sprint 20: 173, +27 from
  voter row + P18 column)
- **Test suite (Sprint 20)**: orchestration 63 passed (was 53);
  voter+P18 e2e 45 fast + 8 slow; cross-detection-matrix 24 passed
  (with new Sprint 20 coverage method)

## Models (18)

| Model | Substrate | Primary Patterns | Sprint | Status |
|-------|-----------|-----------------|--------|--------|
| Zhang cell-view sorting | lattice_1d | P1, P31 | 1 | ✅ |
| Zhang threaded | lattice_1d | P1, P31 | 1 | ✅ |
| Schelling segregation | lattice_2d | P1 | 2 | ✅ |
| Greenberg-Hastings CA | lattice_2d | P13 | 2 | ✅ |
| Game of Life | lattice_2d | P15 | 2 | ✅ |
| Vicsek flocking | continuous_2d | P5 | 2 | ✅ |
| D'Orsogna milling | continuous_2d | P6 | 2 | ✅ |
| Kuramoto oscillators | oscillator | P9 | 3 | ✅ |
| BTW sandpile | lattice_2d | P14 | 4 | ✅ |
| Nowak-May spatial PD | lattice_2d | P27 (+ P1) | 5 | ✅ |
| Hegselmann-Krause opinion | opinion_space | P21 | 5 | ✅ |
| SIR epidemic | lattice_2d | P22 | 7 | ✅ |
| Spatial rock-paper-scissors | lattice_2d | P12 | 9 | ✅ |
| Lotka-Volterra lattice | lattice_2d | P11 (+ P1) | 11 | ✅ |
| Gray-Scott reaction-diffusion | lattice_2d_continuous | P3 | 13 | ✅ |
| Nagel-Schreckenberg traffic | lattice_1d | P8 | 15 | ✅ |
| Active Brownian Particles | continuous_2d | P2 | 16 | ✅ |
| Yard-Sale wealth exchange | scalar_wealth | P28 | 17 | ✅ |
| **Non-local Kuramoto ring** | **oscillator** | **P10** | **18** | ✅ |

## Detectors

Registered in `epc/orchestration.py::DETECTOR_REGISTRY` (17):
P1, P2, P3, P5, P6, P8, P9, **P10**, P12, P13, P14, P15, P21, P22, P27,
P28, P31.

Implemented and test-covered but not registered (1):
**P11** (Sprint 11 LV + P11 predator-prey oscillation detector — the
detector exists as `epc/detectors/p11_predator_prey_oscillation.py`
with DEFINITIVE canonical positive, but is not registered in
orchestration). Documentation-vs-registry consistency gap carried
forward; see REPLICATION_NOTES.md Sprint 15 #3.

Additional detectors: a P13/P15 discriminator lives at
`epc/detectors/p13_p15_discriminator.py` but is an auxiliary utility,
not a standalone detector.

## Sprint history

| Sprint | Theme | Key adds |
|--------|-------|----------|
| 1 | Foundation | Zhang sorting, P1, P31, BaseModel/BaseDetector/DetectorResult |
| 2 | Flocking / wave CA | Vicsek, D'Orsogna, GH, GoL, P5, P6, P13, P15 (v1) |
| 3 | Oscillators | Kuramoto, P9, Schelling |
| 4 | SOC | BTW sandpile, P14 |
| 5 | Social games | Nowak-May, HK, P27, P21 |
| 6 | Paper v0.3, architecture decisions | Detector cards, ontology |
| 7 | Epidemic cascade | SIR, P22 |
| 8 | P15 generalization | p15_persistent_computation (generalized) |
| 9 | Cyclic dominance | RPS spatial, P12 |
| 10 | SIR × P1 reframe | Primary metric: peak → final Moran |
| 11 | Predator-prey | Lotka-Volterra lattice, P11 |
| 12 | Paper catch-up | Sections 1–5 drafted; numerical-results pass |
| 13 | Continuous-field substrate | Gray-Scott, P3, `lattice_2d_continuous` substrate |
| 14 | P1 substrate robustness | Graceful-reject for P1 on continuous-field |
| 14.5 | P1 diagnostic schema | `screening_rejection_reason` + NM size characterization |
| 14.6 | Threshold lock + test split | GS spots tier decision, sandpile test split |
| 15 | Traffic jamming | Nagel-Schreckenberg, P8, 2nd lattice_1d model |
| 16 | MIPS | Active Brownian Particles, P2, metadata-mechanism gate (Decision 43) |
| 17 | Wealth condensation | Yard-Sale, P28, `scalar_wealth` substrate, first well-mixed agent population, Decisions 47–49 |
| **18** | **Chimera states** | **Non-local Kuramoto ring, P10, first 2×2 within-substrate block (oscillator), first IC-dependent canonical positive, Decisions 50–53** |

## Architecture decisions log (53 total, Sprint 18 adds 4)

Decisions 1–36 locked before Sprint 13. Sprint 13 added 37–39
(substrate-content gating; k_max_frac=1.0; n_permutations=199).
Sprint 15 added 40–42 (P8 stopped_fraction primary; lattice_1d
+ integer velocity prereqs; jam_lifetime_p95 confirmation gate).
Sprint 16 added 43–46 (P2 mechanistic-null gate, two_phase_score
primary, confirmation gates, post-burn ≥ 3·T_rot).
Sprint 17 added 47–49 (P28 Gini-not-Pareto primary, DY Boltzmann-
Gibbs surrogate null, four-flag mechanistic-null gate).
Sprint 18 adds:

- **Decision 50**: P10 primary metric is `pos_vel_ac[lag=4]` — spatial
  autocorrelation of time-averaged per-oscillator phase velocity on
  the ring — NOT the pattern-catalog-obvious bimodal `{r_w}` window
  signature. Phase 1h/1i showed that every window-based chimera
  signature (gap, persistence_corr, time_std_ratio, per-window
  coexistence) gives false positives on ordinary all-to-all Kuramoto
  near K_c, because the mean-field model's internal ω-sort produces
  persistent window structure that is window-wise indistinguishable
  from a chimera arc. `pos_vel_ac[4]` discriminates spatial (coupling)
  ordering from frequency ordering, giving a separation gap of +0.47
  with no overlap between chimera (0.93 ± 0.01) and ordinary Kuramoto
  K=1.0 (0.31 ± 0.13). Analogous to Sprint 16 ADR 44 (Hartigan dip
  failure) and Sprint 17 ADR 47 (Pareto α failure).

- **Decision 51**: Canonical positive pins β = 0.05, seed = 0, asymmetric
  Gaussian IC — NOT the paper's β = 0.18. Chimera states are bistable
  with full sync at these parameters; at β = 0.18 only half of random
  seeds land in the chimera basin, while at β = 0.05 all tested seeds
  (0, 1, 42, 200, 500) robustly produce chimeras. This is the first
  sprint where the canonical positive depends on a specific IC and
  seed; future bistable pattern detectors (likely P16 Hopfield, P26
  stochastic resonance) will need the same IC-specific pinning. Paper
  β = 0.18 retained as a slow-half replication test.

- **Decision 52**: P10 DEFINITIVE gate uses two metadata flags:
  `has_nonlocal_coupling = True` AND `has_frequency_heterogeneity !=
  True`. Ordinary all-to-all Kuramoto (heterogeneous ω, no non-local
  kernel) cannot reach DEFINITIVE via P10 under any content
  configuration. Follows the P2 (Decision 43) / P28 (Decision 49)
  template — content-level signal alone caps at CONFIRMATION without
  metadata affirmation of the mechanism. The retrofit on the existing
  `kuramoto.py` (adding `has_nonlocal_coupling = False`,
  `has_frequency_heterogeneity = True`) provides the negative-match
  side of the gate.

- **Decision 53**: P10 coexistence gate uses a drift-invariant per-frame
  measure — `per_frame_coexistence_fraction ≥ 0.90` — NOT a per-window
  persistence measure. Phase 2b discovered that Abrams-Strogatz
  chimera arcs execute slow random-walk translations along the ring
  at long times. A per-window persistence gate (window X coherent ≥
  90% of all frames) produces false negatives at T ≥ 100 because no
  single window is persistently coherent despite the chimera being
  intact. The per-frame measure requires that EACH FRAME has at
  least one coherent and at least one incoherent window (in the same
  frame); a drifting chimera satisfies this at ~100% of frames
  regardless of where the arcs are. Per-window persistence counts
  (`n_persistent_coh`, `n_persistent_incoh`) are retained as
  spatial-stationarity diagnostics in secondary metrics but not used
  as tier gates.

## Paper drafts

| Section | Status | Last sprint updated |
|---------|--------|---------------------|
| §1 Introduction | 1,602 words | 12 |
| §2 Pattern Catalog | 3,121 words | 12 |
| §3 Detection Framework | 1,985 words | 12 |
| §4 Replication Studies | ~10,900 words (§4.15 ABP+P2 ~1,400w) | 16 |
| §5 Cross-Model Transfer | ~3,400 words | 16 |
| §6 (TBD topic) | 1,202 words (Sprint 6 era) | 6 — **needs update** |
| §7 (TBD topic) | 1,393 words (Sprint 6 era) | 6 — **needs update** |
| §8 Conclusion | not drafted | — |
| References | not compiled | — |

**Sprint 17 and Sprint 18 both deferred their paper sections** (§4.18
P28 and §4.19 P10 respectively) to a future paper-catchup sprint. The
empirical/detector work is complete; Sprint 18's paper contribution is
documented in REPLICATION_NOTES.md Sprint 18 section (~440 lines with
Phase 1 characterization tables, Phase 2 detector verification, and
ADRs 50-53) and in the P10 detector card (docs/detector_cards.md
v0.6.2).

Total paper body at Sprint 18 HEAD: ~23,400 words (unchanged from
Sprint 16 — Sprints 17 and 18 both deferred paper updates).

Paper cleanup is expected in a dedicated sprint after Sprints 19-20,
covering items #24–29 (§6/§7 consistency pass, §8 Conclusion draft,
reference list compilation, §4.18 P28 prose, §4.19 P10 prose, §5
transfer matrix table 18×17 display).

## Outstanding carry-forwards

**From earlier sprints (code/numerical):**
1. Gray-Scott inner loop is pure numpy (Sprint 13). Numba would give 10-20×.
2. LV model inner loop is pure Python (Sprint 11). Numba.
3. RPS inner loop is pure Python (Sprint 9). Numba.
4. LV × P11 finite-size scaling slow test (Sprint 11 #1).
5. P15 IC sensitivity characterization (Sprint 8).
6. RPS wavelength scaling λ ∝ √M not replicated (Sprint 9).
7. RPS M_c not pinned precisely (Sprint 9).
8. P11 canonical positive requires ≥ 1200 generations (Sprint 11) —
   documented but not enforced.
9. NM size sensitivity documented but not enforced (Sprint 14.5 D.2).

**Introduced in Sprint 15:**
10. NS inner loop is pure numpy — fast but Numba would accelerate
    long-run finite-size scaling studies.
11. P8 CONFIRMATION tier requires L ≥ 1000 for canonical demonstration.
    Finite-size scaling slow test is a future addition.
12. **P11 not registered in orchestration** (pre-existing; surfaced in
    Sprint 15 paper-table work). Low-priority cleanup.

**Introduced in Sprint 16:**
13. ABP inner loop is cKDTree-per-step (N=1000: ~20ms/step). Grid-based
    density estimator would speed up large-N runs. Low priority.
14. Vicsek and D'Orsogna metadata lack the Sprint 16 P2 rule flags
    (has_alignment_rule, has_attraction_rule,
    has_density_dependent_speed). Retrofitting would let P2 emit more
    informative exclusion reasons for those models. Low priority.
15. P2 finite-size scaling slow test. Canonical primary scales with N.
    A slow-marked test at N ∈ {250, 500, 1000, 2000} analogous to
    Sprint 15 #11 (NS finite-size) would pin the minimum N.

**Introduced in Sprint 17:**
16. YS model inner loop is pure Python per transaction (~500k tx/sec
    at N=1000). Numba JIT would give ~20× and enable N=5000+ slow
    finite-size scaling studies. Low priority.
17. P28 finite-size scaling slow test. Phase 1d.4 verified N-invariance
    at 1000 sweeps across N ∈ {200, 500, 1000, 2000}; a slow-marked
    test pinning the lower N-bound for seed-robustness would
    strengthen the robustness claim. 1 session.
18. Retrofit P28 metadata flags to any future second scalar_wealth
    model. Not urgent until a second wealth model lands. Low priority.

**Introduced in Sprint 18:**
19. Non-local Kuramoto inner loop is O(N²) per RK4 substep, scalar
    numpy. At N=128 one frame ≈ 2s for T=50. Numba or a cached G-matrix
    sparse approximation would give 5-10× for larger-N studies. Low
    priority until N > 512 tests are requested.
20. Chimera arc drift quantification not characterized. Phase 2b
    observed drift qualitatively but did not fit a diffusion
    constant. Low priority.
21. Paper-convention time rescaling not applied. Implementation uses
    a convention that rescales time by N relative to the paper's PDE
    units. Documented in code and REPLICATION_NOTES.md. Low priority;
    cosmetic convention.
22. P10 finite-size robustness slow test (analogous to #15 and #17).
    Phase 1e verified pos_vel_ac[4] at N ∈ {64, 128, 256}; pinning
    the lower N bound would strengthen robustness claim. 1 session.
23. P10 β = 0.18 paper-parameter full replication (incl. lifetime
    measurements and (A, β) parameter-space boundary mapping) was not
    done. Slow-half test confirms seed=42 at β=0.18 reaches
    DEFINITIVE; full replication is 1-2 sessions.

**Paper carry-forwards:**
24. §6 and §7 consistency pass (Sprint 12 #1). Sprint 6-era content
    should reflect Sprints 7–18 findings.
25. §8 Conclusion not drafted (Sprint 12 #2).
26. Reference list not compiled (Sprint 12 #3).
27. **§4.18 P28 section** not drafted (Sprint 17 deferral).
    Empirical work complete in REPLICATION_NOTES.md.
28. **§4.19 P10 section** not drafted (Sprint 18 deferral).
    Empirical work complete in REPLICATION_NOTES.md.
29. §5 transfer matrix table update to 18×17 display (Sprint 18
    deferral, paired with #28).

## Test totals at Sprint 18 HEAD

- **Fast-half**: 235 passed + 15 deselected (~4:40)
- **Heavy-half**: 41 passed + 1 deselected (~3:58)
- **Sandpile-slow**: 3 passed + 2 deselected (~2:57)
- **NS slow**: 3 passed + 19 deselected (~0:07)
- **ABP slow** (Sprint 16): 4 passed + 19 deselected (~0:40)
- **YS slow** (Sprint 17): 5 passed + 30 deselected (~0:10)
- **P10 slow (Sprint 18, NEW)**: 3 passed + 22 deselected (~0:50)
- **Grand total**: **294 passed + 16 deselected** (was 269+16 at
  Sprint 17)

Sprint 18 adds: 22 fast + 3 slow e2e tests in
`test_kuramoto_p10_e2e.py`, 1 new canonical-pair entry and count
updates in `test_orchestration.py`, 1 new coverage test in
`test_cross_detection_matrix.py` (`test_sprint_18_kuramoto_nonlocal_p10_covered`),
34 new `EXPECTED_OUTCOMES` cells, and the audited-pair floor bump
from 112 to 146.

Sprint 18 net fast-half delta: +22 (213 → 235).
