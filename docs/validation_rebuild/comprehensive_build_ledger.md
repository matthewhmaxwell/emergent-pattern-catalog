# Ring 2 — Comprehensive lens build-out ledger

Tracks the build-out of the lens battery to comprehensiveness BEFORE the next hunt (the
"breadth + depth first, hunt once" directive). Source plan: `ring2_lens_coverage_audit.md`
(52 axes; 16 covered + 2 admitted + 3 deferred + 31 missing). Each new lens earns its place
against a positive control + a severable null, with a SPECIFIC discriminator (never a shared
absolute), honoring the h1_max-not-h1_total lesson.

## ADMITTED this build-out (10) — earn-its-place validated, committed

| # | Lens (file) | Discriminator | Validation result |
|---|---|---|---|
| 1 | topological defect census (`defect_census.py`) | winding-number defect density, coherence-gated | active-nematic sea 0.0061 vs ordered/noise 0.0000 (noise gated by coherence 0.17) |
| 2 | pattern symmetry / lattice class (`pattern_symmetry.py`) | azimuthal harmonic m + ring concentration | stripes m2 / square m4 / hex m6; concentration lattice 5.8-9.4 vs labyrinth 2.0; gas gated |
| 3 | velocity-field order from positions (`velocity_order.py`) | connected velocity-fluctuation length xi | flock/mill xi 0.31-0.40 vs gas 0.088; drifting-gas pol 0.66 but xi 0.088 (confound caught) |
| 4 | metastability / regime-switching (`metastability.py`) | n_macrostates + dwell_cv (irregular) | telegraph n2 dwell_cv 1.71; oscillation dwell_cv 0.17 (regular, separated); unimodal/ramp n1 |
| 5 | coarsening dynamics (`coarsening.py`) | domain-growth exponent L(t)~t^n | Allen-Cahn n0.11 R2 0.90 (clean growth) vs static 0.00 vs noise R2 0.06 (gated) |
| 6 | anomalous transport (`anomalous_transport.py`) | MSD exponent beta + ergodicity-break | brownian 1.0 / ballistic 2.0 / fBm-sub 0.70; CTRW EB 0.72 vs brownian 0.10 |
| 7 | sync texture / chimera (`sync_texture.py`) | spatial spread of local order field | KB chimera spread 0.19 (global R 0.20 = looks like partial sync) vs uniform <=0.058 |
| 8 | chaos dimension / determinism (`chaos_dimension.py`) | 0-1 test K + nonlinear-prediction determinism | chaos K1/det0.99, periodic K0/det1.0, white-noise K1/det0.0 |
| 9 | algorithmic compressibility (`compressibility.py`) | LZ76 vs shuffle ratio | random 1.00 vs periodic 0.03 / Thue-Morse 0.11 / rule-90 0.55 |
| 10 | hyperuniformity / giant fluctuations (`hyperuniformity.py`) | number-variance scaling exponent | hyperuniform 1.16 < poisson 1.62 < clustered 2.26 |

## ADMITTED — MEDIUM/LOW + deferred-trio depth (5 more)

| # | Lens (file) | Discriminator | Validation result |
|---|---|---|---|
| 11 | cross-substrate coupling (`cross_substrate.py`) | cross-spectral coherence (joint, zero in marginals) | buried-coupled 0.66 (flat autos) / coupled 0.81 vs independent-patterned 0.21 floor |
| 12 | multifractal spectrum (`multifractal.py`) | Delta-alpha minus shuffled (surrogate-gated) | cascade 0.61 vs iid-heavytail 0.20 (raw 3.5) vs uniform 0.00 -- closes DENSE-fractal gap |
| 13 | critical slowing down (`critical_slowing.py`) | early-vs-late variance + AR(1) rise, linear-detrended | approaching var_rise 0.78 vs stationary/drift 0.017 (drift confound defused) |
| 14 | nonequilibrium current (`nonequilibrium_current.py`) | state-space circulation (broken detailed balance) | directed cycle 0.61 vs standing-osc 0.007 + noise 0.044 |
| 15 | extreme-value clustering (`extreme_value.py`) | extremal index theta (runs declustering) | iid 0.85 vs clustered 0.27 (same tail weight) |

## Deferred trio — RESOLVED (depth)
- multifractal -> ADMITTED (#12).
- correlation-length / criticality -> FSS HARNESS `fss_criticality.py` (2D Ising, Binder crossing +
  chi-peak scaling; CERTIFIED chi_max 8.2->37.5, peak at Tc 2.269). Criticality is a multi-L
  diagnostic, NOT a single-observation coordinate -- the honest architectural closure.
- recurrence / RQA -> SUBSUMED by `chaos_dimension` (#8): delay-embedding determinism + 0-1 test
  capture deterministic-recurrence structure; RQA adds no distinct coordinate. Retired-with-replacement.

## LOW tier — honestly DEFERRED with per-item rationale (not built)
- **long-range memory (DFA Hurst)** -- near-duplicate: H<->MSD-exponent beta (anomalous_transport #6)
  and 1/f<->spectral lens. No distinct coordinate.
- **ageing / two-time relaxation (FDT violation)** -- needs many waiting times + a conjugate-field
  RESPONSE (intervention) passive observation lacks. Fundamentally needs perturbation.
- **cross-frequency coupling (PAC / bicoherence)** -- high false-positive risk (sharp edges /
  harmonics mimic coupling); needs waveform-aware surrogates. Revisit only if a candidate system demands it.
- **hierarchical / nested multi-scale** -- separating genuinely distinct scales from one broad S(k)
  peak / smooth crossover is ill-posed without clear scale separation (finite-size-adjacent).

## Second tier — non-prioritized MISSING axes (survey synthesis left out of the build queue)
Genuinely missing but lower priority-vs-validatability per the synthesis; a future build round:
traveling-wave / space-time (k,omega) dispersion; glassy dynamic heterogeneity (chi_4);
self-replication / population growth; temporal-network dynamics; slow-fast / canards; stochastic
& coherence resonance; developmental staging / time-arrow program; open-ended drift / novelty
generation; the remaining cross-substrate combination sub-cells (one representative built in #11).

## STATUS: prioritized build queue COMPLETE
15 coordinate lenses ADMITTED + 1 FSS harness + 1 subsumed; 4 LOW honestly deferred w/ rationale;
deferred trio fully resolved. All 15 are COORDINATE lenses (NOT classifiers) -> per the integrity
protocol they cannot move tripwire leads. NOT yet wired.

## NEXT PHASE (wire + freeze + hunt)
batch-wire all 15 into `ring2_descriptor.py` (each self-guards -> None off-substrate) -> ONE
descriptor re-validation (corpus + particle_life/RD/Kuramoto/Lenia: 0 nulls trip, classified holds,
leads provably unchanged) -> refresh sweeps on substrates each lens fires on -> FREEZE battery ->
ONE comprehensive hunt. No hunting until frozen.
