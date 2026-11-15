# New literature reproductions for the catalog (plan v2, frozen before any reproduction run)

Written 2026-11-15/16 for the catalog paper (owner: "do all three"). The targets, protocols and pass rules are those in
the three spec files committed with this plan (group_A.md, group_B.md, group_C.md and their support scripts), written
by research agents that read the source papers and the model code and ran no model. This file names the reproductions
that will run, the primary pass rule of each, the patterns dropped and why, and the rules of conduct.

## Reproductions that will run (nine)

| Pattern | Source of the target | Primary pass rule (spec section 4) | Reported, not gating |
|---|---|---|---|
| P5 flocking | Vicsek et al. 1995, Fig. 2a (N = 400, rho = 4) | v_a within 0.04 of the digitized value at all 9 noise values | N = 100, 40; beta |
| P6 milling | D'Orsogna et al. 2006, Eq. (4) ring radius | T6.1: ring reached in >= 8 of 10 seeds at every N, and R within 3% of Eq. (4) at N = 100, 200, 300, 400 | T6.2, T6.3 |
| P10 chimera | Abrams & Strogatz 2004, analytic solution at A = 0.995, beta = 0.18; existence boundary | T10.1 (five observables within their bands) and T10.2 (collapse in [0.2500, 0.2675], none beyond 0.2700) | T10.3 |
| P12 cyclic dominance | Reichenbach, Mobilia & Frey 2007, Fig. 2b extinction probabilities | the three conditions of group_C.md P12 section 4 at L = 100 and 40 | diagnostics D1, D2 |
| P18 consensus | Frachebourg & Krapivsky 1996 exact interface density; Dornic et al. 2001 slope 2/pi | T2 (exact density within 3% at t = 100, 300, 1000, 3000) and T1 (slope in [0.586, 0.688]) | finite-size check |
| P23 anti-coordination | Challet, Marsili & Zecchina 2000, alpha_c = 0.3374 | T1: minimum of sigma^2/N at alpha in [0.29, 0.39] (P = 128) | T2; tie-rule flag |
| P26 stochastic resonance | Gammaitoni et al. 1998, Fig. 7 (Floquet, A = 0.1, Omega = 0.1) | the four conditions of group_C.md P26 section 4 | SNR optimum |
| P27 spatial reciprocity | Szabo & Fath 2007, Fig. 26 (cooperator fraction 0.316 for 9/5 < b < 2) | the three conditions of group_A.md P27 section 4, at b = 1.9 | 0.318 attribution |
| P28 wealth condensation | Boghosian et al. 2017, Eq. 61 steady Gini with redistribution | P28-A: within 0.02 at tau = 0.2, 0.5, 1.0 | P28-B (Boghosian 2014 slope) |

A pattern's reproduction counts as "published value" in the paper only if its primary rule passes. A fail is reported
as a fail, with the measured values, and the pattern keeps its previous grade.

## Dropped (stated in the paper)

- P1 Schelling: no published number for the catalog's rule; the available target needs a different dynamics.
- P17 collective sensing: the source's simulation values depend on model details the catalog model does not share.
- P19 emergent leadership: the catalog model is not Couzin et al.'s model; no parameter setting closes the gap.
- P32 division of labor: the catalog model matches neither the 1996 nor the 1998 source's rules.

## Rules of conduct

1. Model changes are options with the catalog default unchanged (update order for P5, random tie-break and linear
   payoff for P23, faster inner loops that are checked bit-identical against the existing loop where the spec says so).
   Where an option changes the rules (P5 update order, P23 tie-break), the result is reported as a reproduction of
   that variant, and the paper says which runs the detectors use.
2. Before a full run, only smoke tests too short to inform any pass rule.
3. No pass rule, tolerance, target, size or seed count changes after any full-run result exists. A protocol step that
   proves unworkable is logged as a deviation, with its reason, in the results file before the affected pass rule is
   evaluated.
4. Every result is written to analysis/outputs/repro_v2/<pattern>.json with the measured values, the pass/fail of each
   rule, the code commit and the runtime; a results note docs/validation_rebuild/repro_v2_results.md summarizes them.
