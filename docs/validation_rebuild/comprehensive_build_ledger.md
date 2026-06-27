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

All 10 are COORDINATE lenses (descriptor fingerprint axes), NOT classifiers. Per the lens-addition
integrity protocol they CANNOT change tripwire leads -> at battery freeze they batch-wire into
`ring2_descriptor.py` (each self-guards -> None off-substrate) with ONE descriptor re-validation +
refresh of sweeps on substrates each fires on. NOT yet wired (wiring deferred to freeze, so it
happens once after the battery is complete, not per-lens).

## REMAINING before freeze (then hunt once)

- **Cross-substrate combination family** (queue MEDIUM): cross-spectrum S_AB(k), closed-loop
  bidirectional cross-TE, cross-susceptibility -- each provably zero in both marginals. Extends the
  validated P36/P37 Ring-1 hits. (Largest remaining item.)
- **Deferred trio — DEPTH revisit with proper methods:**
  - multifractal spectrum f(alpha) width for DENSE fractals (the deferred `fractal_dimension` gap)
  - correlation-length / criticality via finite-size scaling (needs varying L)
  - recurrence via delay-embedding (not surrogate) -- now partly subsumed by `chaos_dimension`;
    re-assess whether RQA adds anything beyond determinism + 0-1 K.
- **Critical slowing down** (spun out of #4): AR(1)+variance rise pre-tipping; needs a saddle-node
  tipping-point positive control to validate honestly.
- **LOW tier** (queue flagged near-duplicate or needs perturbation): cross-frequency coupling (PAC),
  nonequilibrium current / entropy production, hierarchical multi-scale, long-range memory (DFA),
  extreme-value (GEV/extremal index), ageing (two-time). Build the non-duplicate ones; honestly
  defer those that need an intervention channel passive observation lacks.

## Sequence
build remaining -> batch-wire ALL admitted into descriptor -> single re-validation (corpus +
RD/Kuramoto/Lenia, confirm 0 nulls trip / classified holds / leads unchanged) -> FREEZE battery ->
ONE comprehensive hunt. No hunting until the battery is frozen.
