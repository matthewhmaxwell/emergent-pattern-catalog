# Emergent Pattern Catalog — Full Catalog Completion Summary (32/32)

**Date:** 2026-06-10 (Sprint 94)
**Status:** All 32 patterns implemented. 31/32 AT-DEPTH. Milestone B complete.

---

## Full 32-Pattern Table

| # | Pattern | Canonical Reference | AT-DEPTH | dim1 Anchor | Limitation |
|---|---------|-------------------|----------|-------------|------------|
| P1 | Similarity-driven aggregation | Schelling 1971 | Yes | Quantitative (Moran's I) | — |
| P2 | Activity-induced phase separation (MIPS) | Fily & Marchetti 2012 | Yes | Quantitative (two_phase_score ≥ 0.10, \|r\| ≥ 0.70) | — |
| P3 | Turing pattern formation | Gray-Scott (Pearson 1993) | Yes | Quantitative (p/m ratio) | — |
| P4 | Territoriality / exclusion boundaries | Giuggioli, Potts & Harris 2011 | Yes | Quantitative (exclusivity > 0.85, overlap < 0.10) | — |
| P5 | Translational alignment / flocking | Vicsek et al. 1995 | Yes | Quantitative (order parameter) | — |
| P6 | Milling / vortex formation | D'Orsogna et al. 2006 | Yes | Quantitative (\|L\| milling parameter) | — |
| P7 | Lane formation in counterflow | Helbing & Molnár 1995 | Yes | Quantitative (lane order parameter) | — |
| P8 | Self-organized jamming | Nagel-Schreckenberg 1992 | Yes | Quantitative (fundamental diagram) | — |
| P9 | Temporal synchronization | Kuramoto 1975 | Yes | Quantitative (K_c = 2γ analytical) | — |
| P10 | Chimera states | Abrams & Strogatz 2004 | Yes | Quantitative (chimera index) | — |
| P11 | Predator-prey oscillations | Mobilia-Georgiev-Täuber 2007 | Yes | Quantitative (O(1/L) amplitude scaling, exponent −0.967 vs published −1.0) | — |
| P12 | Cyclic dominance (spatial RPS) | Reichenbach-Mobilia-Frey 2007 | **No** | Qualitative (spirals present; λ ∝ √M slope = 0.161, outside [0.4, 0.6]) | **dim1 PARTIAL**: finite-size measurement limitation; 4 attempts (Sprints 54/58/59/63); formula-valid M range too narrow for reliable slope fit. Validated via panel PASS + dim2 PASS + qualitative spirals. |
| P13 | Excitable spiral and target waves | Greenberg-Hastings 1978 | Yes | Quantitative (wavefront persistence) | — |
| P14 | Self-organized criticality | Bak-Tang-Wiesenfeld 1987 | Yes | Quantitative (τ = 1.2914 ± 0.0012) | — |
| P15 | Persistent propagating computation | Conway Game of Life | Yes | Quantitative (TE + reproducibility + perturbative IC) | dim2 N/A (deterministic) |
| P16 | Associative memory / pattern completion | Hopfield 1982 / Amit-Gutfreund-Sompolinsky 1985 | Yes | Quantitative (storage capacity α transition ≈ 0.173) | — |
| P17 | Distributed sensing / collective gradient | Berdahl et al. 2013 | Yes | Quantitative (CI slope vs log(N) = 0.133) | — |
| P18 | Collective consensus / decision-making | Voter model (Holley-Liggett 1975) | Yes | Quantitative (wall Spearman + consensus time) | — |
| P19 | Emergent leadership / minority guidance | Couzin et al. 2005 | Yes | Quantitative (accuracy-vs-ρ curve, Spearman ρ = 1.0) | — |
| P20 | Quorum sensing / threshold activation | Waters & Bassler 2005 | Yes | Quantitative (step-function R² = 1.000, hysteresis width = 1.190) | — |
| P21 | Polarization / fragmentation | Hegselmann & Krause 2002 | Yes | Quantitative (cluster-count vs ε, 8 points all within tolerance) | — |
| P22 | Information cascade / social contagion | Datta & Acharyya 2021 (SIR) | Yes | Quantitative (wavefront speed 0.461 vs published 0.441, error 4.7%) | — |
| P23 | Anti-coordination / emergent load balancing | Savit, Manuca & Riolo 1999 (Minority Game) | Yes | Quantitative (σ²/N interior minimum at α ≈ 0.32) | — |
| P24 | Homeostatic regulation | Ashby 1956 (homeostat) | Yes | Quantitative (deviation ratio = 0.0027) | — |
| P25 | Canalized restoration / equifinality | Waddington 1957 | Yes | Quantitative (convergence variance ratio ≈ 0.0, basin volume = 1.0) | — |
| P26 | Stochastic resonance | Gammaitoni 1998 / Collins 1995 | Yes | Quantitative (peak coherent response = 0.918, inverted-U confirmed) | — |
| P27 | Spatial reciprocity / emergent cooperation | Nowak & May 1992 | Yes | Quantitative (cooperation fraction) | — |
| P28 | Wealth condensation / spontaneous inequality | Chakraborti-Boghosian 2002/2014 (Yard-Sale) | Yes | Quantitative (Gini convergence) | — |
| P29 | Trail / network formation | Tero et al. 2010 (Physarum) | Yes | Quantitative (network efficiency) | — |
| P30 | Spontaneous boundary formation / autopoiesis | Varela-Maturana-Uribe 1974 (SCL) | Yes | Quantitative (association_score = 2.211 ± 0.087, closure = 0.996) | — |
| P31 | Delayed gratification | Zhang et al. 2024 | Yes | Quantitative (swap counts within 4%, ΔR² = +0.645) | — |
| P32 | Emergent specialization / division of labor | Bonabeau 1996 | Yes | Quantitative (entropy decline + coverage ≥ 0.5) | — |

---

## AT-DEPTH Summary

- **AT-DEPTH: 31 / 32**
- **GAP: 1 / 32** — P12 (cyclic dominance, dim1 finite-size measurement limitation)

### Per-wave roll-up

| Wave | Patterns | Implemented | AT-DEPTH |
|------|----------|-------------|----------|
| Milestone A (original 19) | P1–P15, P18, P21, P22, P27, P28, P31 | 19/19 | 18/19 (P12 GAP) |
| Wave 1 (Sprints 66–70) | P7, P17, P19 | 3/3 | 3/3 |
| Wave 2 (Sprints 72–77) | P24, P26, P23 | 3/3 | 3/3 |
| Wave 3 (Sprints 79–84) | P16, P25, P20 | 3/3 | 3/3 |
| Wave 4 (Sprints 86–93) | P4, P29, P32, P30 | 4/4 | 4/4 |
| **Total** | **32** | **32/32** | **31/32** |

---

## Cumulative Carry-Forwards (Honest Known-Issues List)

### Open

| ID | Pattern | Issue | Severity | Notes |
|----|---------|-------|----------|-------|
| C-p12-dim1 | P12 | λ ∝ √M scaling-law reproduction fails (slope outside [0.4, 0.6] after 4 attempts) | Accepted limitation | Finite-size measurement limitation at L ≤ 200; formula-valid M range too narrow. P12 validated via panel + dim2 + qualitative spirals. |
| C-p7-time-shuffled-fp | P7 | time_shuffled FP at screening | Low | Each frame preserves lane structure independently of temporal order; cosmetic. |
| C-p9-constant-field | P9 | constant_field Class A trips trivial sync | Low | Degenerate substrate; cosmetic. |
| C-p19-bias-zero-chance-alignment | P19 | 1 Class C bias_zero FP at confirmation | Low | Chance alignment in 1/5 bias_zero regimes. |
| C-p21-time-shuffled-fp | P21 | time_shuffled FP at confirmation (0.850) | Low | Pre-convergence unimodal steps; cosmetic. |
| C-p14-class-c-borderline | P14 | 1 borderline at p_diss=0.350 | Low | Near-SOC dissipation regime; cosmetic. |
| C-p30-enrichment-cv | P30 | Enrichment ratio CV=34.7% across seeds | Informational | Association_score CV=3.9% is robust; enrichment is a secondary metric. |

### Closed (key closures)

- C-p27-panel-screening-leak — CLOSED Sprint 40 (observable-prerequisite guard)
- C-p22-class-b-cascade-overlap — CLOSED Sprint 41 (irreversibility prereq)
- C-p1-linear-gradient-fp — CLOSED Sprint 61 (multi-cluster prerequisite)
- C-p1-class-c-subthreshold-fp — CLOSED Sprint 62 (regime correction)
- C-p8-class-c-near-onset — CLOSED Sprint 62 (density range correction)
- C2/C3 (P12 λ scaling) — CLOSED-AS-DOCUMENTED-LIMITATION Sprint 63
- C4 (P14 multi-seed) — CLOSED Sprint 55
- C5 (P21/P22 methods notes) — CLOSED Sprint 57

---

## Transfer-Matrix Final Dimensions

The cross-detection transfer matrix covers the **original 19 Milestone-A detectors × models** registered in `epc/orchestration.py`. The full Phase-2a panel system covers all 32 patterns via per-pattern panels with 3 negative classes (Class A synthetic, Class B catalog/supplements, Class C failed regimes).

- **33 model files** in `epc/models/`
- **34 detector files** in `epc/detectors/`
- **19 registered detectors** in `DETECTOR_REGISTRY`
- **7 substrate types**: lattice_1d, lattice_2d, lattice_2d_continuous, continuous_2d, oscillator, opinion_space, scalar_wealth
- **32 Phase-2a panels** (one per pattern), each with TNR ≥ 0.95

---

## Paper Readiness

### Sections drafted

| Section | Status | Notes |
|---------|--------|-------|
| §1 Introduction | Draft | `docs/paper_section1_draft.md` |
| §2 Related Work | Draft | `docs/paper_section2_draft.md` |
| §3 Methods | Draft | `docs/paper_section3_draft.md` |
| §4 Results (per-pattern) | Draft | `docs/paper_section4_draft.md` — all 32 patterns covered (§4.1–§4.32) |
| §5 Cross-Detection Analysis | Draft | `docs/paper_section5_draft.md` |
| §6 Discussion | Draft | `docs/paper_section6_draft.md` |
| §7 Conclusion | Draft | `docs/paper_section7_draft.md` |
| §8 Appendices | Draft | `docs/paper_section8_draft.md` |

### Methods notes authored (14 patterns)

P2, P4, P7, P16, P17, P19, P20, P21, P22, P23, P24, P25, P26, P29, P30, P32 — in `docs/methods_notes/`.

### Remaining for PLOS ONE submission

1. **Milestone C instrument layer** (calibration, novelty/none-of-the-above, OOD validation, external demonstration, public API) — see `docs/instrument_roadmap.md`
2. Paper revision pass to incorporate final Wave 3–4 results into §4, §5, §6
3. Figures and supplementary materials
4. Abstract and author contributions

---

## Next Phase: Milestone C (Instrument / OOD Layer)

Per `docs/instrument_roadmap.md`:

- **T2a** — Calibration layer (conformal-style calibrated scores)
- **T2b** — Novelty / "none-of-the-above" abstention signal
- **T2c** — OOD validation suite (held-out systems, precision/recall)
- **T2d** — External demonstration (multi-agent LLM/AI system profiling)
- **T2e** — Public API/CLI (`epc-detect <observation-bundle>`)

**The catalog is complete. Chain halted for operator review before Milestone C.**
