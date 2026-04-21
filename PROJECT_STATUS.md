# Project Status Tracker

Last updated: Sprint 17 (April 2026)

## Current state (post Sprint 17)

The catalog has progressed through 17 sprints of work. Sprint 1 built
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
completing the three-class discrimination framework. **Sprint 17**
adds Yard-Sale + P28 (wealth condensation) — the first well-mixed
(non-spatial) agent population in the registry, introducing the new
`scalar_wealth` substrate. Sprint 17 is also the first sprint of the
Scenario-A catalog-completion campaign, targeting the remaining ~16
patterns that lack canonical positive-model + detector pairs.

## Inventory snapshot at Sprint 17 HEAD

- **17 model families** (18 model files, Zhang sequential and threaded
  count as the same family)
- **16 registered detectors** in `epc/orchestration.py::DETECTOR_REGISTRY`
  + 1 additional column for P11 (Sprint 11 work implemented but not
  retrofitted into the registry — pre-existing gap; see Sprint 15 #3
  and Sprint 16 #12 in REPLICATION_NOTES.md)
- **7 substrate types**: lattice_1d, lattice_2d, lattice_2d_continuous,
  continuous_2d, oscillator, opinion_space, **scalar_wealth** (new
  at Sprint 17)
- **112 audited cells** in `tests/test_cross_detection_matrix.py::
  EXPECTED_OUTCOMES` (Sprint 16: 78 → Sprint 17: 112, +34 from
  yard_sale row + P28 column)
- **Test suite**: **269 passed + 16 deselected** across fast-half
  (213 + 15 deselected), heavy-half (41 + 1), sandpile-slow
  (3 + 2), NS slow (3 + 19), ABP slow (4 + 19), and YS slow
  (5 + 30 when run alone)

## Models (17)

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
| **Yard-Sale wealth exchange** | **scalar_wealth** | **P28** | **17** | ✅ |

## Detectors

Registered in `epc/orchestration.py::DETECTOR_REGISTRY` (16):
P1, P2, P3, P5, P6, P8, P9, P12, P13, P14, P15, P21, P22, P27, **P28**, P31.

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
| **17** | **Wealth condensation** | **Yard-Sale, P28, `scalar_wealth` substrate, first well-mixed agent population, Decisions 47–49** |

## Architecture decisions log (49 total, Sprint 17 adds 3)

Decisions 1–36 locked before Sprint 13. Sprint 13 added 37–39
(substrate-content gating; k_max_frac=1.0; n_permutations=199).
Sprint 15 added 40–42 (P8 stopped_fraction primary; lattice_1d
+ integer velocity prereqs; jam_lifetime_p95 confirmation gate).
Sprint 16 added 43–46 (P2 mechanistic-null gate, two_phase_score
primary, confirmation gates, post-burn ≥ 3·T_rot).
Sprint 17 adds:

- **Decision 47**: P28 primary metric is the Gini coefficient, NOT
  the Pareto tail exponent α. The pre-existing pattern-catalog entry
  prescribed "Pareto power-law tail" as a detection signature; Phase
  1c Hill-estimator characterization showed α drifts unstably across
  timescales — at short time α > 2 (sub-Pareto), transiently α ∈ (1,
  2), at long time α < 1 (degenerate). No stable window exists in
  which a fixed-α gate discriminates condensation regimes. Gini is
  stable, monotonic, and has a clean null under the Dragulescu-
  Yakovenko equilibrium. α_hill is carried as a diagnostic secondary
  metric only. Analogous to Sprint 16's ADR 44 (Hartigan dip → two-
  phase score).

- **Decision 48**: P28 null model is the well-mixed Boltzmann-Gibbs
  distribution (Dragulescu-Yakovenko 2001). Under the null hypothesis
  "symmetric exchange of a conserved scalar resource equilibrates to
  a Boltzmann-Gibbs distribution", the Gini of a sample from Exp(⟨w⟩)
  is ≈ 0.5 in the large-N limit. We draw N samples from Exp(observed
  mean_w), compute Gini, repeat n_permutations (≥ 199 for floor
  p = 0.005). Right-tailed p. NullType = SURROGATE (not SHUFFLE)
  because we sample from a theoretical equilibrium.

- **Decision 49**: P28 mechanistic-null gate uses FOUR metadata flags
  simultaneously. DEFINITIVE requires:
  `has_conserved_resource = True`,
  `has_multiplicative_stake = True`,
  `has_saving_propensity = False`,
  `has_redistribution = False`.
  Each flag represents an independent mechanism that can block
  condensation (CC 2000 saving propensity and chi > 0 redistribution
  both produce finite-Gini fixed points even with a full multiplicative
  stake rule). The four-flag requirement is the third generation of
  the three-class discrimination framework: substrate-level (registry)
  / substrate-content (observable values, Decisions 37, 41) /
  metadata-mechanism (rule flags, Decisions 43, now 49).

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

**Sprint 17 deferred §4.18 + §5 table update** to a future paper-catchup
sprint. The empirical/detector work is complete; Sprint 17's paper
contribution is documented in REPLICATION_NOTES.md Sprint 17 section
(~445 lines with Phase 1 characterization tables and ADRs 47-49) and
in the P28 detector card (docs/detector_cards.md v0.6.1).

Total paper body at Sprint 17 HEAD: ~23,400 words (unchanged from
Sprint 16 — Sprint 17 paper update deferred as a carry-forward).

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

**Paper carry-forwards:**
19. §6 and §7 consistency pass (Sprint 12 #1). Sprint 6-era content
    should reflect Sprints 7–17 findings.
20. §8 Conclusion not drafted (Sprint 12 #2).
21. Reference list not compiled (Sprint 12 #3).
22. **§4.18 P28 section** not drafted (Sprint 17 deferral). The
    empirical work is complete in REPLICATION_NOTES.md; the paper
    prose catch-up is deferred to a paper-cleanup sprint.
23. §5 transfer matrix table update to 17×16 display (Sprint 17
    deferral, paired with #22).

## Test totals at Sprint 17 HEAD

- **Fast-half**: 213 passed + 15 deselected (~4:14)
- **Heavy-half**: 41 passed + 1 deselected (~3:58)
- **Sandpile-slow**: 3 passed + 2 deselected (~2:57)
- **NS slow**: 3 passed + 19 deselected (~0:07)
- **ABP slow** (Sprint 16): 4 passed + 19 deselected (~0:40)
- **YS slow (Sprint 17, NEW)**: 5 passed + 30 deselected (~0:10)
- **Grand total**: **269 passed + 16 deselected** (was 213+11 at
  Sprint 16)

Sprint 17 adds: 30 fast + 5 slow e2e tests in
`test_yard_sale_p28_e2e.py`, 8 registration tests in
`test_orchestration.py` (`TestSprint17Registrations` class),
1 new coverage test in `test_cross_detection_matrix.py`
(`test_sprint_17_yard_sale_p28_covered`), 34 new `EXPECTED_OUTCOMES`
cells, and 2 other small registration tests.

Sprint 17 net fast-half delta: +41 (172 → 213).
