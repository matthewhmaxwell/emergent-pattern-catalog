# Deviations from PLANTED_PIPELINE_PLAN.md (commit 85acbe1)

Each entry is logged with its reason before any affected pipeline output (screen, knock-out, name, card, gate) exists.

## P1. N19 headline moved to the point the paper's own figure shows (2026-11-03, before the code freeze)
The plan's table gave N19's headline as "the paper's parameters (card)", the values in the caption of the paper's
Fig. 2 (kappa_a = 5, kappa_n = 14, beta = 0.02). The plug-in author found that the paper's own linear-stability
result (its Eq. 2) makes the uniform state stable at those values (largest growth rate about -6e-4; threshold near
beta = 0.018 at kappa_a = 5), and census-scale runs stayed uniform. The paper says Fig. 2 corresponds to point K of its
Fig. 3A, plotted at kappa_a = 11.1 +- 0.2, beta = 0.084 +- 0.003 (read at 300 dpi by the module author and, separately,
by the reviewer: 11.10 +- 0.19, 0.0839 +- 0.003), where a run shows one aggregate with travelling bands as in Fig. 2.
The plan's general rule is a headline "inside the range where its paper reports the new behavior"; the caption values
are not. N19 therefore runs at point K (values read off Fig. 3A; the reading and its uncertainty are recorded in the
module docstring). The fidelity check runs at point K; the caption values are reported as a diagnostic.

## P2. N11 fidelity check follows the authors' own labels (2026-11-03, before the code freeze)
The plan labelled N11's headline g = 3.14e-5 "(stripes)" and fixed the check "stripes at the headline and hexagonal
spots at g = 3.27e-6", taken from the card. The authors' code README and simulation data label 3.14e-5 the gap
pattern (anisotropy index 0.169) and put stripes at 2.356e-5 (0.916); the card read the panels of the paper's Fig. 1
in the wrong order. Gaps are also part of what the paper presents as new, so the headline value is kept. The deciding
fidelity check becomes agreement with the authors' anisotropy index at each simulated g (within 0.15), with spots at
3.27e-6; the plan's original "stripes at the headline" result is reported alongside.

## P3. N11 census-scale domain (2026-11-03, before the code freeze)
The plan uses the paper's domain when it holds at most 10,000 particles, otherwise "a domain holding about 400 at close
packing". N11's individuals are points (no close packing), and its paper's domain (side 10) holds about 160,000 at
steady state. A domain holding 400 would be smaller than the competition radius (1.25), so no pattern could exist.
N11 runs on the paper's domain, with the time per census step reduced under the plan's compute cap (recorded in the
module).

## P4. Fidelity tolerances corrected after code review (2026-11-03, before the code freeze and before any full
fidelity result was read)
- N12: the module author's criterion "fitted log-decay exponent within 0.07 of 0.34" is dominated by fit noise (the
  lever arm of ln ln L between L = 128 and 512 is 0.19); it now passes if within max(0.07, 3 fit standard errors),
  fitted on sizes L >= 128 only, else it is reported as not decisive. The point-by-point comparison with the paper's
  fitted M^2(L) (tolerance 0.01) is unchanged.
- N13: the exact-distribution part of the check counted state-level maxima, which appear almost everywhere; it now uses
  the same binned-prominence rule as the simulated part.
- N01: the exact enumeration is also run through the code path the census simulation uses.

## P5. N13 shock rate scaled with population size (2026-11-03, before the code freeze)
The plan keeps "the paper's scaling" for parameters the paper defines per population size. The paper's stationary
solutions (appendix eqs. 34-35) depend on the shock rate only through N times the rate per unit time, i.e. N^2 times
the per-step probability, so the per-step shock probability is scaled as 10^-4.5 x (500/N)^2 (1.98e-4 at 200 nodes,
1.2e-5 at 800). Holding it fixed per step, the module's first choice, would weaken the middle peak at 200 nodes to 3%
of the tallest (paper's setting: 18%) and change the shape again at 800 nodes in the verification stage.

## P6. Enforcing the verification compute cap (2026-11-04, after the code freeze and the screen/name/card stages, before any verification output)
The plan's step 5 records a condition that "cannot finish within 24 CPU-hours" as not run; the frozen driver has no
limit. The flag stage showed the module authors' compute-cap timings were low for the slowest models (single-threaded
census-scale cost of the real program, 3 seeds: N04 8.0 h against an estimate of about 50 min; N03 4.4 h; N11 3.5 h).
`census/planted_verify_capped.py` runs the driver's own job list and worker (`planted_pipeline._vworker`, unchanged)
under this rule, fixed before any verification job ran:
- projected cost = flag-stage seconds x 10/3 seeds x (1 census, 4 for 4x area, 3 for 3x length); for particles at
  fixed density and for fields this is a lower bound on the scaled cost;
- a condition projected over 36 CPU-hours (1.5 x the cap) is recorded as not run without an attempt (N03, N04 and N11
  at 4x area; N03 and N04 at 3x length);
- every other job runs in its own process with a 24-CPU-hour limit (RLIMIT_CPU); a job killed by the limit is not run;
- sub-programs of a model whose 4x-area or 3x-length condition is not run are not run: the frozen `classify_one`
  returns ARTIFACT (or NOT RUNNABLE) before it reads sub-programs, so they cannot change the class.
Behaviors lost this way are reported separately as the plan requires ("lost to the compute cap"), not as artifacts of
the pipeline's judgement.
The compute-cap reductions of the census step (Notes below) were set from the module authors' timings; the census-scale
runs therefore took longer than the plan's 60 CPU-minutes for N03, N04 and N11. Longer runs only bring the census run
closer to the paper's protocol; no result is changed by this, and it is reported as a limitation.

## Operational note (2026-11-04)
The first gate launch passed the run name positionally (the CLI expects `--run`); three processes wrote to gate/None/
and completed one card (pipeline:N02:F0) before being stopped. That output is kept in
census/runs/planted_pipeline/superseded/ and is not scored; r1-r3 were relaunched with `--run`.

## Notes (within the plan, recorded for transparency)
- Compute-cap reductions of the time per census step (plan, "Census scale"): N04 0.2 time units, N11 1e-3 time units;
  values and timings in the modules.
- The 24 CPU-hour verification cap is applied per condition (census, 4x area, 3x length) for each model.
- Expected fidelity failures traced to the papers themselves, recorded before the full checks were read: N02's "absent
  at c = 0.49" leg (the paper's own Eq. 22 makes the uniform state unstable for every c > 0 on the c1 = c, c2 = -c
  branch; its Appendix A scheme appears linearised); N03's phase diagram (the printed equations bind clusters of 4 to 7
  or more disks at every screening value up to 0.5, by the reviewer's edge-current calculation).
- N04: after review the transverse-force length scale is 1.16 (the value that reproduces the paper's Table 1 odd
  moduli; the first version's 1.0 gave about a tenth of them). The implementation's passive crystal melts near T = 14
  against the paper's about 5, and no change to a stated parameter reconciles it, so N04's fidelity check is expected
  to fail and its results carry "implementation not confirmed".
- N11 blow-up rule: the plan's 50x rule is kept as frozen for the knock-out and every sub-program; the pair-only
  sub-program settles at 47-51x and may stop as a blow-up by fluctuation (known before any result).
- N14: one census step is the paper's Monte Carlo cycle of 1,296 attempted moves (about one sweep at the paper's
  close-packed number of particles), since the particle number grows from a dilute gas; no blow-up rule applies (the
  knock-out reaches the ideal-gas number). The chemical potential follows the module's convention (thermal wavelength 1,
  no orientation factor), which the paper does not state.
- Time units after code review: N08 one census step = 1/D_r, the paper's declared time unit (the first version used
  the Lennard-Jones unit); N03 one census step = tau = 1/D_r, the unit the paper reports times in, reduced by the
  compute cap.
- Rule texts edited beyond numbers where the census run differs in kind from the paper's protocol: N12 says
  "single-spin updates, one sweep of all sites per step, at beta = 1, from random angles, for 1000 sweeps" in place of
  the paper's cluster updates at several temperatures versus system size; N17 drops "and watch the dense travelling
  stripes that form" (a description of the paper's long runs, not of the 1,000-step census run) and states the census
  start (random headings along the four lattice directions).
