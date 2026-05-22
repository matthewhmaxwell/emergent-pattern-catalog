# Depth-Gap Audit Matrix

**Audit date:** 2026-05-09
**HEAD commit at audit start:** `101430f69f44753804721b6b37123a3c9d7f98fb` (Sprint 26, v0.26.0)
**Rubric:** `docs/depth_rubric.md` (Sprint 28 deliverable, byte-identical to bundle)
**Patterns audited:** 19 (all entries in `epc/orchestration.py::DETECTOR_REGISTRY`)
**At-depth count:** 5
**Gap count:** 14

> **Methodology.** Read-only audit. For each pattern, the auditor located the
> detector card in `docs/detector_cards.md`, the canonical positive model in
> the codebase, the pattern's section(s) in `REPLICATION_NOTES.md`, and its
> tests in `tests/`. Each of the four rubric dimensions was scored
> PASS / PARTIAL / FAIL / N/A. When scoring required judgment (most often:
> "is the methods note substantive enough?" and "does the cross-detection
> transfer matrix count as a Phase-2a negative sweep?"), the auditor
> defaulted to the stricter score per the brief and flagged the row in
> `SPRINT_RETURN.md` for chat review.

> **Effort sizing.** S = ≤1 sprint (paragraph-add, single-seed→multi-seed
> rerun, methods-note draft, or test extension). M = 1–2 sprints (new
> multi-seed campaign, dedicated Phase-2a sweep, or extended methods
> document). L = 3+ sprints (new canonical replication, sourcing of new
> canonical paper, substantial new model work).

| pattern_id | name | dim1_replication | dim2_multiseed | dim3_methods_note | dim4_negative_sweep | grade | notes | effort_to_close |
|---|---|---|---|---|---|---|---|---|
| P1 | Similarity-driven aggregation (Schelling) | PASS | PASS | PASS | PARTIAL | GAP | dim4: cross-detection matrix only; no dedicated Phase-2a sweep beyond NM × P1 finite-size characterization (Sprint 14.5) | M |
| P2 | Activity-induced phase separation (MIPS) | PARTIAL | PARTIAL | PARTIAL | PARTIAL | GAP | dim1: Cates & Tailleur cited, no specific Fig/table reproduced quantitatively; dim2: Sprint 16 has multi-seed but ≥5-seed bar + dispersion not clearly reported per pattern; dim3: methods note covers two_phase_score primary but unit/discretization choices implicit; dim4: cross-detection only | M |
| P3 | Turing pattern formation (Gray-Scott) | PASS | PASS | PASS | PARTIAL | GAP | dim4: Sprint 13 has spots/stripes/labyrinthine + Sprint 14.6 5-seed threshold characterization, but no dedicated substrate-diverse Phase-2a sweep | M |
| P5 | Translational alignment / flocking (Vicsek) | PASS | PASS | PASS | PARTIAL | GAP | dim4: cross-detection matrix only; specificity not reported across substrate-diverse non-positives | M |
| P6 | Milling / vortex formation (D'Orsogna) | PASS | PARTIAL | PASS | PARTIAL | GAP | dim2: deterministic RK4 with parameter-dependent milling; ≥5-seed dispersion not documented; dim4: cross-detection only | M |
| P8 | Self-organized jamming (Nagel-Schreckenberg) | PASS | PASS | PASS | PARTIAL | GAP | dim4: jam_lifetime_p95 confirmation gate discriminates pigeonhole density saturation but no dedicated Phase-2a sweep across ≥3 substrate types | M |
| P9 | Temporal synchronization (Kuramoto) | PASS | PASS | PASS | PASS | AT-DEPTH | dim4: K_c=2γ analytical confirmation + Phase-2a panel v1.2 PASS-with-weakness (Sprint 35): overall TNR=0.952 (≥0.95), syn=0.875 (weak: only `constant_field` trips — degenerate trivial-sync substrate; new carry-forward C-class-a-constant-field-trivial-sync), cat=1.000 (advisory n=3), fai=1.000 (10/10 sub-K_c), Cohen's d=4.781, see `analysis/outputs/p9_phase2a_panel.json`. v1.2 SKIPs the permutation-degenerate substrates that drove the v1.1 PARTIAL. AT-DEPTH grade earned via v1.2 PASS-with-weakness. | — |
| P10 | Chimera states (non-local Kuramoto ring) | PASS | PASS | PASS | PARTIAL | GAP | dim4: Sprint 26 phase-diagram boundary closes Sprint 18 #23 (topology PASS, lifetime inconclusive at T_max=100) but no broad Phase-2a sweep against non-oscillator substrates | M |
| P11 | Predator-prey oscillations (Lotka-Volterra) | PARTIAL | PASS | PASS | PARTIAL | GAP | dim1: Mobilia-Georgiev-Täuber 2007 cited; LV phase diagram qualitative match documented but no specific Fig/table quantitative reproduction; dim4: bilateral-vs-cyclic exclusion is a content-level discriminator, not a substrate-diverse Phase-2a sweep | M |
| P12 | Cyclic dominance (spatial RPS) | PARTIAL | PARTIAL | PASS | PARTIAL | GAP | dim1: Reichenbach 2007 cited; phase-diagram topology qualitative match but λ ∝ √M wavelength scaling NOT replicated (Sprint 9 carry-forward, still open); dim2: characterization mostly single-seed across mobilities; dim4: cross-detection only | L |
| P13 | Excitable spiral and target waves (GH) | PASS | PASS | PASS | PARTIAL | GAP | dim4: Sprint 9 RPS-vs-GH wavefront-CV discrimination is a content-level test, not a substrate-diverse Phase-2a sweep | M |
| P14 | Self-organized criticality (BTW) | PASS | PARTIAL | PASS | PASS | GAP | dim2: τ=1.247 from a single 100k-event run; ≥5-seed bootstrap dispersion not reported (open dim2 gap, not closed by this sprint). dim4: Phase-2a panel v1.2 PASS (Sprint 35): overall TNR=0.960 (≥0.95), syn=1.000 (8/8 evaluated; 2 SKIPPED), cat=1.000 (7/7 lattice_2d mates), fai=0.900 (1 borderline at p_diss=0.350; carry-forward C-p14-class-c-borderline, low priority), Cohen's d=10.585, see `analysis/outputs/p14_phase2a_panel.json`. dim4 PASS earned via v1.2 invariance-flag fix + clean substrate-typed Class B. Pattern remains GAP-narrowed (dim2 only) — close dim2 to reach AT-DEPTH. | M |
| P15 | Persistent propagating computation (GoL) | PASS | N/A | PASS | PASS | AT-DEPTH | dim2 N/A: GoL is fully deterministic and the canonical positive (R-pentomino, glider collisions) is not phase-boundary-dependent; dim4: TE discriminator + multi-checkpoint reproducibility + GH/SIR/Nowak-May rejection collectively constitute a substrate-diverse Phase-2a sweep. Phase-2a panel v1.1 PASS (Sprint 33): overall TNR=1.000, syn=1.000, cat=1.000 (7/7 lattice_2d mates), Class C N/A, Cohen's d=8.282, see `analysis/outputs/p15_phase2a_panel.json`. AT-DEPTH grade now positively confirmed by both content-level discriminators AND v1.1 panel result. Sprint 33 sanity-check role validated. | — |
| P18 | Collective consensus / decision-making (voter) | PASS | PASS | PASS | PASS | AT-DEPTH | Sprint 20/21/24 produced canonical voter L=64 positive, multi-seed validation across L ∈ {64,128,256}, substantive ADR 54/55/56 methods notes, six-row discriminator rejection table including Schelling false-positive closure (Sprint 24 #20b). Phase-2a panel v1.1 PASS (Sprint 32): overall TNR=1.000, syn=1.000, cat=1.000 (advisory n=3, network: P21 + 2 supps), Class C N/A (voter has no parameter regime that suppresses consensus), Cohen's d=+inf, see `analysis/outputs/p18_phase2a_panel.json`. AT-DEPTH grade now positively confirmed by both content-level discriminators AND v1.1 panel result. | — |
| P21 | Polarization / fragmentation (Hegselmann-Krause) | PARTIAL | PARTIAL | PARTIAL | PARTIAL | GAP | dim1: HK 2002 cited; ε-dependent cluster count qualitatively matches but no specific Fig/table reproduction; dim2: ε-sweep mostly single-seed; dim3: methods note thin (≤2 grep hits in REPLICATION_NOTES); dim4: cross-detection only | L |
| P22 | Information cascade / social contagion (SIR) | PARTIAL | PARTIAL | PARTIAL | PARTIAL | GAP | dim1: Datta & Acharyya CA cited (Sprint 8 reference correction); percolation transition replicated qualitatively but no specific Fig/table quantitative match documented; dim2: per-N basin fraction not reported; dim3: methods note exists for substitution from Fuks-Lawniczak but light on integration choices; dim4: cross-detection only | L |
| P27 | Spatial reciprocity / emergent cooperation (Nowak-May) | PASS | PASS | PASS | PARTIAL | GAP | dim4: Sprint 14.5 5-seed × 5-grid finite-size characterization is content-level, not substrate-diverse Phase-2a | M |
| P28 | Wealth condensation / spontaneous inequality (Yard-Sale) | PASS | PASS | PASS | PASS | AT-DEPTH | Sprint 17 produced YS canonical positive vs Chakraborti-Boghosian 2002/2014, multi-seed Gini convergence, substantive methods note, four-flag mechanistic-null gate is the explicit Phase-2a-style discriminator | — |
| P31 | Delayed gratification (Zhang sorting) | PASS | PASS | PASS | PASS | AT-DEPTH | Zhang et al. 2024 swap counts within 4%, insertion DG within 4%, 600-run non-redundancy test (ΔR²=+0.645, p<0.000001) is the canonical Phase-2a-style negative test against shuffled-DG ablation | — |

## Aggregate findings

- **At-depth: 5 / 19** (P9, P15, P18, P28, P31). P9 added in Sprint 35 via Phase-2a panel v1.2 PASS-with-weakness.
- **Gap: 14 / 19**.
- **Most common FAIL/PARTIAL dimension: dim4 (broad negative sweep)** — 15/19 patterns scored PARTIAL on dim4. Cross-detection transfer matrix coverage is universal but the rubric explicitly calls this insufficient (a Phase-2a-style sweep specifically tests detector specificity against models *not* in the catalog). Only P15, P18, P28, P31 have content-level discriminators or mechanistic-null gates that meet the dim4 PASS bar.
- **Second most common: dim1 (replication)** — 4/19 patterns scored PARTIAL (P2, P11, P12, P21, P22). Pattern: canonical paper named, qualitative behavior matches, but no specific named figure/table reproduced with stated tolerance.
- **Effort distribution:** 0 S (none), 13 M, 2 L (P12, P21, P22). The two L-effort gaps are P21 and P22 (decision-making / cascade patterns from the more-recent Sprint 5/7 era where multiple dimensions are PARTIAL together) and P12 (RPS λ ∝ √M wavelength scaling, an open Sprint 9 carry-forward).

## Carry-forwards surfaced (do NOT act on per Sprint 28 brief)

The audit surfaced these candidate carry-forwards. They are recorded in `SPRINT_RETURN.md` for chat review.

- C1: Dimension 4 (Phase-2a negative sweep) is the catalog's single largest depth gap. A "Phase-2a uniformity sprint" defining a standard substrate-diverse negative panel and applying it to all currently-PARTIAL detectors would close 15/15 dim4 gaps in one campaign.
- C2: Dimension 1 (specific Fig/table quantitative reproduction) is missing for 5 patterns whose qualitative replication is documented but whose tolerances against named results are not. P2 (MIPS), P11 (LV), P12 (RPS), P21 (HK), P22 (SIR) all need a specific paper figure named with tolerance.
- C3: P12 RPS λ ∝ √M wavelength scaling not replicated (open since Sprint 9 carry-forward, still open at Sprint 26 HEAD).
- C4: P14 BTW τ-bootstrap multi-seed dispersion not reported (Sprint 4 single-run convention; the canonical positive bootstrap is a one-paragraph methods extension, not a new run).
- C5: P21 (HK) and P22 (SIR) methods notes are thin; both detectors would benefit from a methods-note expansion sprint (S/M effort) before further extension work.
