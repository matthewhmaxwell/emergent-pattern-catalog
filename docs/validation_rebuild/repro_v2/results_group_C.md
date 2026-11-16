# Literature reproductions v2, group C results (P26, P12, P28)

Specs: `group_C.md` (frozen at 39688e3), sections 3 (Protocol) and 4 (Tolerance and pass rule). Plan and rules of
conduct: `REPRO_V2_PLAN.md`. Results files: `analysis/outputs/repro_v2/P26.json`, `P12.json`, `P28.json`; raw
per-run records in `analysis/outputs/repro_v2/raw/`. Scripts: `analysis/repro_v2/`.

## Notes written before any full run

These were recorded and committed before the first full-run result existed. None changes a target, tolerance,
size, seed count or pass rule.

Code changes (all options; catalog defaults unchanged):
- P26: `double_well_fourier_driver` in `epc/models/stochastic_resonance.py`, a vectorised driver applying the
  update rule of `BistableDoubleWell.simulate()` to all trials at once and returning Fourier projections.
  Equivalence test (`tests/test_sr_vectorized_driver.py`): 4 trials x 2000 steps at D = 0.15 with the same
  pre-drawn normals as `simulate()`; trajectories agree to <= 1e-12 (passed).
- P12: `RPSSpatialModel(..., kernel="numba")` in `epc/models/rps_spatial.py`, a compiled inner loop applying the
  same pre-drawn arrays in the same order; default `kernel="python"`. Bit-identity test
  (`tests/test_rps_numba_kernel.py`): 3 seeds x 50 generations at L = 30, M in {1e-5, 1e-3}; grids and reaction
  counts identical at every generation (passed). Numba 0.68.1 (llvmlite 0.50.0, numpy 2.5.0) was installed into
  the shared venv for this; pytest 9.2.0 was also installed there to run the tests.
- P28: no model change. Parameters only, and the drivers call `step(R)` with R = 2500.

Spec details fixed in the scripts (listed in full in each script's `IMPLEMENTATION_CHOICES` and copied into the
results JSON): exact M grid values, floor(G_max), Wilson interval, quadratic-vertex bootstrap seed, the reading
of the dt rule and of the stationarity rule, the P28-B fit and solve.

Two spec observations, recorded before running:
1. P12 seeds. The spec's formula 120000 + 1000 index(L) + 10 index(M) + r with r = 0..499 overlaps across M values
   within one L (index(M) = 1, r = 0 and index(M) = 0, r = 10 both give 120010). The spec's remark that streams
   never overlap holds across L only. Every (L, M) point still has 500 distinct seeds; neighbouring M points share
   initial grids and random streams (common random numbers). The formula is applied as written.
2. P28-B units. In Boghosian 2014, a(t) in Eqs. (35)-(37) is the coefficient of P_c(w) = a/w normalised to N
   agents, so a = N / ln(w_max/w_min), and Fig. 4 shows a ~ 1e5-1e6 for the egalitarian start. The fitted line
   1/a(t) = 3.93264 + 2.04046e-5 t implies a(0) ~ 0.25, which under Eq. (37) would put w_min near exp(-2e5) at
   t = 0. With W = 5e6, N = 5e4 and the poorest agents near ln w ~ -150 by t ~ 1e8 (source, p. 11), Eq. (37) as
   the spec writes it gives 1/a of order 1e-3 at t = 8e7, i.e. a slope of order 1e-11 per transaction. P28-B as
   specified is therefore expected to fail by several orders of magnitude for a reason in the target, not the
   model. It is run and reported as specified; it does not gate P28.

Compute: measured on the shared VPS before the runs, the P12 L = 100 generation costs about 0.5 ms (spec: about
0.3 ms), and the P28 pure-Python transaction about 3.2 us (spec: 1.5-2 us). Both are within about 2x of the spec's
estimates, so the protocols run as written. The optional P12 L = 200 runs (30-45 CPU-h in the spec) are not run.

Compute log, 2026-11-16 01:50 UTC (P26 complete; P12 L = 40 primary complete; P12 L = 100 primary and P28 running;
no P12 or P28 rule evaluated). From the run times of the first 207 L = 100 realisations, a generation costs about
0.73 ms on the shared VPS (four other compute processes are running), so the P12 L = 100 primary projects to about
20 CPU-h (spec: about 8) and diagnostic D2 at L = 100 to about 19 CPU-h (spec: about 5 for D1 and D2 together; D1
costs nothing because it is read from the primary runs). This is slow but workable, so nothing is changed: the
order is P12 primary, then D2 at L = 40, then D2 at L = 100, within the two-process limit. Reading those run times
also showed that most early low-M realisations ran to G_max. No target, tolerance, size, seed count or rule has
been or will be changed.

## P26 stochastic resonance: PASS

Target: Gammaitoni et al. 1998 RMP Fig. 7, Floquet spectral amplification at A = 0.1, Omega = 0.1 (digitized:
eta_max = 12.23 at D = 0.150; eta(0.20) = 10.8; eta(0.30) = 6.5). Model: `BistableDoubleWell` rule with a = b = 1,
A = 0.1, f = 0.1/(2 pi), dt = 0.01, run through the vectorised driver. 15 D values x 10 seeds (1001-1010) x 400
trials; burn-in 10 periods, measurement 50 periods.

| Rule (spec section 4) | Target and band | Measured (mean +- SE over seeds) | Outcome |
|---|---|---|---|
| 1. eta_max (quadratic vertex on D in [0.12, 0.18]) | 12.23, [11.62, 12.84] | 12.416 +- 0.068 (bootstrap 95%: 12.29-12.55) | pass |
| 2. D_max | 0.150 +- 0.015 | 0.1537 +- 0.0003 (bootstrap 95%: 0.1531-0.1545) | pass |
| 3. eta(0.20) | 10.8, [9.94, 11.66] | 10.757 +- 0.049 | pass |
| 4. eta(0.30) | 6.5, [5.98, 7.02] | 6.524 +- 0.041 | pass |

Measured curve eta(D), seed means: 1.68 (0.06), 4.73 (0.08), 8.40 (0.10), 11.05 (0.12), 11.84 (0.13), 12.25 (0.14),
12.41 (0.15), 12.31 (0.16), 12.04 (0.17), 11.70 (0.18), 10.76 (0.20), 8.39 (0.25), 6.52 (0.30), 4.23 (0.40),
2.99 (0.50); SE 0.01-0.08.

Required checks:
- dt check (D = 0.10, 0.15, 0.30 at dt = 0.005, same seed numbers): relative differences 0.98%, 0.43%, 0.46%, all
  under 2%, so the gated values come from dt = 0.01 and no rerun was needed.
- Frame check (canonical a = 4, b = 1, A = 0.8, Omega = 0.4, D x 16, dt = 0.0025): 16 eta agrees with the a = b = 1
  values to 4e-15. With the same seeds the canonical Euler map is an exact rescaling of the a = b = 1 map, so this
  confirms the scaling and the code path; it is not an independent statistical replicate.

Reported, not gated: the seed-mean SNR has one interior local maximum on the grid, at D = 0.14 (3-point vertex
0.136), against the two-state value Delta V / 2 = 0.125. Down to D = 0.06 the SNR falls with decreasing D; the
intrawell divergence the source describes as D -> 0 does not appear in this range.

No deviation. Runtime 5,235 CPU-s (1.45 CPU-h; spec estimate < 1 CPU-h): 3,273 s main, 1,285 s dt check, 678 s
frame check. Code commit of all runs: 90cfbe2.

## P28 wealth condensation with redistribution: PASS (P28-A); P28-B FAIL (reported only)

### P28-A (primary)

Target: Boghosian et al. 2017 Physica A Eq. (61), G_inf = 1/(1 + 2.3682 tau^0.5943), from steady-state solutions of
the redistributive Fokker-Planck equation. Model: `YardSale` unchanged, N = 5000, f = 0.05, lambda_save = 0,
R = redistribute_every = 2500 = N/2, chi = tau f^2, equal start, advanced only by step(2500). Seeds 2801-2810 per
tau; burn-in max(20/tau, 20) tau-bar; Gini sampled after every tax over a 50 tau-bar window (20,000 samples per
seed).

| tau | Eq. (61) | Measured G (mean +- SE, 10 seeds) | Deviation | Rule: |dev| <= 0.02 |
|---|---|---|---|---|
| 0.2 | 0.5236 | 0.51213 +- 0.00019 | -0.0114 | pass |
| 0.5 | 0.3893 | 0.38286 +- 0.00011 | -0.0065 | pass |
| 1.0 | 0.2969 | 0.29781 +- 0.00006 | +0.0009 | pass |

P28-A PASS. The deviation grows as tau falls (all three are below or at the fit, largest at tau = 0.2); a factor-2
error in the tau mapping would have moved G by 0.09-0.10 (Eq. 61 at 2 tau: 0.421, 0.297, 0.219).

Required checks, applied as fixed in the spec:
- Stationarity (10-seed means of first and second half of the window must agree within 0.005): differences
  -0.00027, +0.00014, -0.00016 at tau = 0.2, 0.5, 1.0 (largest single-seed half difference 0.0015). No rerun.
- Finite-f check (tau = 0.5 at f = 0.025, chi = 3.125e-4): G = 0.38315 +- 0.00014, |G(0.05) - G(0.025)| = 0.0003,
  under the 0.01 threshold, so the gated values are the f = 0.05 values and no extrapolation runs were needed.

Reported, not gated: tau = 0.5 from `init_mode='exponential'` (start above G_inf) gives 0.38310 +- 0.00011,
against 0.38286 from the equal start, so the steady state is approached from both sides to within 0.0003.

### P28-B (reported, does not gate P28)

Target: Boghosian 2014 PRE Fig. 5 fit, slope of 1/a(t) = 2.04046e-5 per transaction (band [1.73, 2.35]e-5). Model:
N = 5e4, f = 0.25, no tax, w0 = 100, seeds 2821-2830, 16 times t = 5e6 ... 8e7; a(t) from the fitted w_min and
Eq. (37) as the spec writes it.

Mean slope 2.679e-11 +- 0.008e-11 per transaction (all 10 seeds within 2.674-2.683e-11); mean intercept 1.23e-4
(source 3.93). P28-B FAIL, by a factor of about 7.6e5. This is the outcome recorded as expected in the pre-run
notes: with a normalised to N agents as in Eqs. (35)-(37), 1/a = ln(w_max/w_min)/N, which is of order 1e-3 here.
The run itself is consistent with the source's description of the dynamics: the fitted ln w_min falls from -4.5 at
t = 5e6 to -103 at t = 8e7, and the poorest agent reaches ln w = -158 at t = 8e7 (source, p. 11: about -150 at
t = 1e8). 1/a is close to linear in t, as the source reports, but on a scale about 1e6 times smaller than the
published fit, so the published numbers are not in the units of Eq. (37).

No deviation. Runtime 26,299 CPU-s (7.3 CPU-h; spec estimate about 4 CPU-h): tau = 0.2, 0.5, 1.0 main 5,002 + 2,880
+ 2,277 s; f = 0.025 check 11,073 s; exponential start 2,688 s; P28-B 2,379 s. Code commit of all runs: 90cfbe2.
