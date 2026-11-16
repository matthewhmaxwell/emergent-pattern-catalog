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
