# Project Status Tracker

Last updated: Sprint 15 (April 2026)

## Current state (post Sprint 15)

The catalog has progressed through 15 sprints of work. Sprint 1 built
the foundation (Zhang sorting + P1 + P31). Sprints 2–11 filled out the
core model inventory and detectors across every substrate type.
Sprint 12 was a paper catch-up sprint. Sprint 13 added the first
continuous-field model (Gray-Scott + P3) and the first continuous-field
substrate type. Sprints 14, 14.5, and 14.6 were cleanup sprints
(P1 substrate robustness, diagnostic-schema polish, threshold
characterization, test hygiene). Sprint 15 added Nagel-Schreckenberg
traffic + P8 — the second `lattice_1d` model and the first with
integer-valued velocity observables.

## Inventory snapshot at Sprint 15 HEAD

- **15 model families** (15 model files)
- **14 registered detectors** in `epc/orchestration.py::DETECTOR_REGISTRY`
  + 1 additional displayed column for P11 (Sprint 11 work implemented
  but not retrofitted into the registry — pre-existing gap carried
  forward; see Sprint 15 #3 in REPLICATION_NOTES.md)
- **6 substrate types**: lattice_1d, lattice_2d, lattice_2d_continuous,
  continuous_2d, oscillator, opinion_space
- **80 audited cells** in the registry transfer matrix; **71 audited
  cells** in the paper-display transfer matrix (15×15 with P11
  displayed); **68 cells** in the cross-detection matrix regression
  table
- **Test suite**: 199 passed + 9 deselected across fast-half
  (152 + 6 deselected), heavy-half (41 + 1), sandpile-slow (3 + 2),
  and NS slow (3 + 19 deselected when run alone)

## Models (15)

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

## Detectors

Registered in `epc/orchestration.py::DETECTOR_REGISTRY` (14):
P1, P3, P5, P6, P8, P9, P12, P13, P14, P15, P21, P22, P27, P31.

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
| **15** | **Traffic jamming** | **Nagel-Schreckenberg, P8, 2nd lattice_1d model** |

## Architecture decisions log (42 total, Sprint 15 adds 3)

Decisions 1–36 locked before Sprint 13 (see earlier project notes and
commit history). Sprint 13 added decisions 37–39 (substrate-level
discrimination via content-level gates; k_max_frac = 1.0 for
peak-to-mean; n_permutations = 199 for P3). Sprint 15 added:

- **Decision 40**: P8 primary metric is `stopped_fraction = ⟨1[v_i(t) = 0]⟩`
  (Bette-Habel-Emig-Schreckenberg 2017 order parameter).
- **Decision 41**: P8 substrate prerequisites are `lattice_1d`
  registration + 1D integer `velocities` observable + velocity range
  [0, 64] + n_cars ≥ 20 + run length ≥ 100 post-burn-in.
  Substrate-level discrimination (cf. Decision 37 for P3).
- **Decision 42**: P8 CONFIRMATION tier requires `jam_lifetime_p95 > 5`
  to discriminate emergent NS jamming from pigeonhole density
  saturation. Null model is per-car temporal shuffle of v(t) on the
  jam_lifetime_p95 statistic.

## Paper drafts

| Section | Status | Last sprint updated |
|---------|--------|---------------------|
| §1 Introduction | 1,602 words | 12 |
| §2 Pattern Catalog | 3,121 words | 12 |
| §3 Detection Framework | 1,985 words | 12 |
| §4 Replication Studies | ~9,500 words (Sprint 15 added §4.14 NS+P8 ~1,350w) | 15 |
| §5 Cross-Model Transfer | ~3,300 words (table + prose updated to 15×15) | 15 |
| §6 (TBD topic) | 1,202 words (Sprint 6 era) | 6 — **needs update** |
| §7 (TBD topic) | 1,393 words (Sprint 6 era) | 6 — **needs update** |
| §8 Conclusion | not drafted | — |
| References | not compiled | — |

**Total paper body at Sprint 15 HEAD**: ~22,000 words (was ~20,437 at
Sprint 14.6; Sprint 15 added ~1,350 words to §4 and minor edits to §5).

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
11. P8 CONFIRMATION tier requires L ≥ 1000 for canonical demonstration
    (ρ = 0.12 fluctuates below 0.05 screening at L = 500). Finite-size
    scaling slow test is a future addition.
12. **P11 not registered in orchestration** (pre-existing; surfaced in
    Sprint 15 paper-table work). Not a Sprint 15 regression, but needs
    a future cleanup.

**Paper carry-forwards:**
13. §6 and §7 consistency pass (Sprint 12 #1). Sprint 6-era content
    should reflect Sprints 7–15 findings.
14. §8 Conclusion not drafted (Sprint 12 #2).
15. Reference list not compiled (Sprint 12 #3).

## Test totals at Sprint 15 HEAD

- **Fast-half**: 152 passed + 6 deselected (~3:25)
- **Heavy-half**: 41 passed + 1 deselected (~4:07)
- **Sandpile-slow**: 3 passed + 2 deselected (~2:57)
- **NS slow**: 3 passed + 19 deselected (~0:07)
- **Grand total**: **199 passed + 9 deselected** (was 176+1 at Sprint 14.6)

Sprint 15 adds: 19 fast + 3 slow e2e tests in `test_ns_p8_e2e.py`,
1 new coverage test in `test_cross_detection_matrix.py`
(`test_sprint_15_ns_p8_covered`), and 17 new `EXPECTED_OUTCOMES` cells.
