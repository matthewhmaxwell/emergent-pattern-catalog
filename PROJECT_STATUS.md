# Project Status Tracker

Last updated: Sprint 16 (April 2026)

## Current state (post Sprint 16)

The catalog has progressed through 16 sprints of work. Sprint 1 built
the foundation (Zhang sorting + P1 + P31). Sprints 2–11 filled out the
core model inventory and detectors across every substrate type.
Sprint 12 was a paper catch-up sprint. Sprint 13 added the first
continuous-field model (Gray-Scott + P3) and the first continuous-field
substrate type. Sprints 14, 14.5, and 14.6 were cleanup sprints
(P1 substrate robustness, diagnostic-schema polish, threshold
characterization, test hygiene). Sprint 15 added Nagel-Schreckenberg
traffic + P8 — the second `lattice_1d` model and the first with
integer-valued velocity observables. **Sprint 16** added Active
Brownian Particles + P2 — the third `continuous_2d` model (alongside
Vicsek and D'Orsogna) and the first detector whose DEFINITIVE tier
depends on a metadata-based mechanistic-null gate (Decision 43). This
completes the three-class discrimination framework: substrate-type
(registry), substrate-content (observable values, Decisions 37/41),
metadata-mechanism (rule flags, Decision 43).

## Inventory snapshot at Sprint 16 HEAD

- **16 model families** (17 model files, Zhang sequential and threaded
  count as the same family)
- **15 registered detectors** in `epc/orchestration.py::DETECTOR_REGISTRY`
  + 1 additional displayed column for P11 (Sprint 11 work implemented
  but not retrofitted into the registry — pre-existing gap carried
  forward; see Sprint 15 #3 in REPLICATION_NOTES.md)
- **6 substrate types**: lattice_1d, lattice_2d, lattice_2d_continuous,
  continuous_2d, oscillator, opinion_space
- **94 audited cells** in the registry transfer matrix (+14 from
  Sprint 15's 80); **78 audited cells** in the paper-display transfer
  matrix (16×16 with P11 displayed) regression table
  `tests/test_cross_detection_matrix.py::EXPECTED_OUTCOMES`
- **Test suite**: 213 passed + 11 deselected across fast-half
  (172 + 10 deselected), heavy-half (41 + 1), sandpile-slow
  (3 + 2), and NS slow (3 + 19 deselected when run alone)

## Models (16)

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
| **Active Brownian Particles** | **continuous_2d** | **P2** | **16** | ✅ |

## Detectors

Registered in `epc/orchestration.py::DETECTOR_REGISTRY` (15):
P1, **P2**, P3, P5, P6, P8, P9, P12, P13, P14, P15, P21, P22, P27, P31.

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
| **16** | **MIPS** | **Active Brownian Particles, P2, metadata-mechanism gate (Decision 43)** |

## Architecture decisions log (46 total, Sprint 16 adds 4)

Decisions 1–36 locked before Sprint 13. Sprint 13 added 37–39
(substrate-level discrimination via content-level gates;
k_max_frac = 1.0; n_permutations = 199 for P3). Sprint 15 added
40–42 (P8 stopped_fraction primary; lattice_1d + integer velocity
prereqs; jam_lifetime_p95 confirmation gate). Sprint 16 adds:

- **Decision 43**: P2 mechanistic discrimination via metadata flags.
  Continuous_2d substrate is shared by P2/P5/P6, so substrate-level
  gating cannot separate MIPS from flocking/milling. DEFINITIVE tier
  for P2 requires metadata affirmation of three rule flags:
  `has_density_dependent_speed=True`, `has_attraction_rule=False`,
  `has_alignment_rule=False`. Missing or negative flags cap the
  detection at CONFIRMATION. This is the metadata-mechanism analogue
  of substrate-content gates (Decisions 37, 41), completing the
  three-class discrimination framework.

- **Decision 44**: P2 primary metric is `two_phase_coexistence_score
  = min(f_gas, f_liquid)`, NOT Hartigan dip on the density histogram.
  The pre-existing detector-card recipe (v0.5.5) specified Hartigan
  dip, which Phase 1c proved empirically wrong for this substrate:
  local densities are integer counts divided by constant area,
  producing discrete distributions that are trivially non-uniform
  by Hartigan's test regardless of underlying physics. Dip floored
  at bootstrap p=0.005 across every tested regime including truly
  uniform (dilute) and truly one-phase (stuck). The two-phase
  coexistence primary exploits the mechanistic signature directly:
  simultaneous presence of dilute and dense phases.

- **Decision 45**: P2 confirmation gates are three-part and
  simultaneous: -0.99 < density_speed_r < -0.30 (anti-correlation
  without the Poisson-discrete artifact), CV_v > 0.30 (rules out
  constant-speed and thermal regimes), frac_stalled < 0.98 (rules
  out fully-stuck clusters). The upper bound on r and the
  frac_stalled gate were added after Phase 1d/e revealed that
  naive anti-correlation and dynamic-range gates each had specific
  false-positive traps (dilute Poisson artifact and stuck cluster
  respectively).

- **Decision 46**: P2 post-burn measurement must be ≥ 3·T_rot to
  distinguish steady-state MIPS from transient coarsening. At high
  packing fraction (φ ≳ 0.7), short-runtime traces (< 2·T_rot)
  show transient two-phase structure en route to the one-phase
  stuck limit — which would false-positive DEFINITIVE if not gated.
  For canonical Pe = 100, T_rot = 100 time units = 2000 steps at
  dt = 0.05, so ≥ 6000-step post-burn is the conservative
  recommendation for φ ≥ 0.7. Documented in the P2 detector card.

## Paper drafts

| Section | Status | Last sprint updated |
|---------|--------|---------------------|
| §1 Introduction | 1,602 words | 12 |
| §2 Pattern Catalog | 3,121 words | 12 |
| §3 Detection Framework | 1,985 words | 12 |
| §4 Replication Studies | ~10,900 words (Sprint 16 added §4.15 ABP+P2 ~1,400w) | 16 |
| §5 Cross-Model Transfer | ~3,400 words (table + prose updated to 16×16) | 16 |
| §6 (TBD topic) | 1,202 words (Sprint 6 era) | 6 — **needs update** |
| §7 (TBD topic) | 1,393 words (Sprint 6 era) | 6 — **needs update** |
| §8 Conclusion | not drafted | — |
| References | not compiled | — |

**Total paper body at Sprint 16 HEAD**: ~23,400 words (was ~22,000 at
Sprint 15; Sprint 16 added ~1,400 words to §4 and ~150 words to §5).

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
    Sprint 15 paper-table work). Not a Sprint 15 regression, but needs
    a future cleanup.

**Introduced in Sprint 16:**
13. ABP inner loop is cKDTree-per-step (N=1000: ~20ms/step; N=5000:
    ~200ms/step). Grid-based density estimator would speed up large-N
    finite-size scaling runs. Low priority — canonical regime already
    fits in test budget.
14. Vicsek and D'Orsogna metadata lack the Sprint 16 P2 rule flags
    (`has_alignment_rule`, `has_attraction_rule`,
    `has_density_dependent_speed`). Retrofitting would let P2 emit
    more informative exclusion reasons for those models. Low priority.
15. P2 finite-size scaling slow test. Canonical primary scales with N
    (N=400 seed-metastable; N=800 primary~0.17; N=1000 primary~0.34).
    A slow-marked test at N ∈ {250, 500, 1000, 2000} analogous to
    Sprint 15 #11 (NS finite-size) would pin the minimum N.

**Paper carry-forwards:**
16. §6 and §7 consistency pass (Sprint 12 #1). Sprint 6-era content
    should reflect Sprints 7–16 findings.
17. §8 Conclusion not drafted (Sprint 12 #2).
18. Reference list not compiled (Sprint 12 #3).

## Test totals at Sprint 16 HEAD

- **Fast-half**: 172 passed + 10 deselected (~4:04)
- **Heavy-half**: 41 passed + 1 deselected (~3:58)
- **Sandpile-slow**: 3 passed + 2 deselected (~2:57)
- **NS slow**: 3 passed + 19 deselected (~0:07)
- **ABP slow (Sprint 16)**: 4 passed (when run with `-m slow`) (~0:40)
- **Grand total**: **213 passed + 11 deselected** (was 199+9 at Sprint 15)

Sprint 16 adds: 19 fast + 4 slow e2e tests in `test_abp_p2_e2e.py`,
1 new coverage test in `test_cross_detection_matrix.py`
(`test_sprint_16_abp_p2_covered`), and 34 new `EXPECTED_OUTCOMES` cells.
