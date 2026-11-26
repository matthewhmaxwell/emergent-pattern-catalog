# Literature reproductions v2, group B results: P5, P6, P10

Plan: `REPRO_V2_PLAN.md` (frozen at 39688e3). Specs: `group_B.md` sections P5, P6, P10 (P19 dropped). Every
measured value, every rule outcome, the seeds, the sizes, the runtimes and the code commit are in
`analysis/outputs/repro_v2/P5.json`, `P6.json` and `P10.json`.

## Status

Code and the implementation readings below were committed before any full run. Full runs had not started at
the time of that commit. The per-pattern sections are filled in after the runs.

## How the runs were made

- Model changes are options; the catalog defaults are unchanged and were checked bit-identical to 39688e3 at
  default settings (`VicsekModel`, `KuramotoNonlocal`).
  - `VicsekModel(position_update="old_velocity")`: Vicsek 1995 Eq. (1), positions move with theta(t).
  - `VicsekModel(wrap_guard=True)`: a position that wraps to exactly L is set to 0 before the cKDTree search.
  - `KuramotoNonlocal(dt_paper=...)`: dt = dt_paper x 2 pi / N, so dt_paper = 0.025 is the paper's step.
  - `KuramotoNonlocal(init_mode="as2004")`: the paper's IC, 6 exp(-30 x^2) U[-1/2, 1/2].
  - P6 needs parameter settings only.
- Drivers: `analysis/reproductions/repro_v2/` (`p5_vicsek1995.py`, `p6_dorsogna2006.py`,
  `p10_abrams_strogatz2004.py`, queue `group_b_queue.py`). Every measurement streams its quantity from
  `step()` or from short `run()` chunks; nothing stores every step.
- Two worker processes on the VPS (`nice -n 10`, one thread each), sharing one task queue.
- Smoke tests used only tiny sizes and lengths (N <= 40, a few hundred steps) that inform no pass rule.
- P10 timing check (1,000 RK4 steps, before the full run): 2.2 ms/step at N = 128 and 8.5 ms/step at N = 256
  on the loaded machine, inside the spec's 10 ms/step budget.

## Implementation readings declared before the full runs

group_B.md leaves some details open. The readings used are listed in
`analysis/outputs/repro_v2/declarations_group_B.json` and copied into each P*.json. The main ones:

- P5: `old_velocity` is the synchronous reading of Vicsek's Eqs. (1)-(2). A "cell" in the convergence rule is
  one (N, eta, seed) run, extended in place to 2x10^5 measurement steps when its halves differ by more than 0.02.
- P6: gates are applied to window time means; the open-space angular momentum is
  `AngularMomentumMetric.compute_instant(box_size=None)`; a seed extended to T = 400 is re-evaluated on
  [350, 400]. The init_radius sensitivity named in section 2 has no protocol in section 3 and is not run.
- P10: locked set = modal 0.001 bin of Omega plus oscillators within 0.002 of it; dt check on the two
  lowest-numbered Stage A chimera seeds; Stage B holds 0.20, 0.22, 0.24 (500 each) then 0.2400 to 0.2800 by
  0.0025 (1000 each); collapse = min r over the last 200 > 0.95 and f_drift < 0.02 there; one confirmation hold
  after the first collapse, then the ramp stops. Optional N = 512 is not run.

## Deviations

Logged in `analysis/outputs/repro_v2/deviations_group_B.json` (copied into each P*.json), each before the rule it
touches is evaluated.

- **P10-D1 (logged 2026-11-16 04:05 UTC, before T10.1 and T10.2 are evaluated): the P10 dt check is reported as
  not run.** The declared reading names the two lowest-numbered Stage A chimera seeds, and no paper-IC seed formed a
  chimera (seeds 0-4 at the time of logging; final at 05:12 UTC: none of seeds 0-9). The dt check is reported only; no
  gating rule, target, tolerance, size or seed count changed. Repeating the fallback realizations at
  dt_paper = 0.0125 instead was considered and not done, because it was not what was declared. Approved by the
  coordinator under the owner's standing rule.

## Notes logged during the runs

- **P6 dt check (pre-registered): passed.** Seed 2000, R at dt = 0.01 vs 0.005: N = 100, 0.17629847 vs
  0.17629849 (change 9e-8); N = 400, 0.07803428 vs 0.07803574 (change 1.9e-5). Both are below 0.1%, so T6.1 runs
  at dt = 0.01.
- **P10 IC transcription check.** The arXiv posting (nlin/0407045v1, Fig. 1 caption) gives phi(x) = a exp(-b x^2) r
  with a = 6, b = 30, r uniform on [-1/2, 1/2], x in [-pi, pi], dt = 0.025, 200,000 iterations, A = 0.995,
  beta = 0.18, N = 256. This matches group_B.md and `init_mode="as2004"`.
- **P10 paper-IC seeds 0-9 at N = 256, beta = 0.18: all ten full sync** (basin fraction 0/10). Each has
  r_bar = 1.0000, f_drift = 0, R_min = R_max = 1.0000 and Delta = 0.9838 = cos(0.18). That is the in-phase state's
  rotation rate in paper units, which confirms the time conversion. Because fewer than 5 of 10 seeds reached a
  chimera, the pre-registered fallback runs (seeds 100, 101, ...), and T10.1 and T10.2 are evaluated on fallback
  realizations, labelled as such.


## Results

Aggregated on 2026-11-26 by `group_b_queue.py aggregate`, the evaluation code committed before the full runs
(9be14e5; the only later change to `analysis/reproductions/repro_v2/` is the read-only `gating_done.py` check and a
change in 08afe11 that copies the deviation log into each results file). All runs finished on 2026-11-17 with no error records. No target, tolerance, size,
seed count or pass rule was changed. The one deviation (P10-D1, above) concerns a reported-only check.

| Pattern | Primary rule | Outcome | Main measured values |
|---|---|---|---|
| P5 flocking | v_a within 0.04 of Fig. 2a at all 9 eta, N = 400 | FAIL | within at eta = 0.5 and 1.0; 0.054 to 0.292 above the figure at eta = 1.5 to 4.5 |
| P6 milling | T6.1: >= 8 of 10 seeds pass all gates at every N, and R within 3% of Eq. (4) | FAIL | 0 of 10 seeds valid at every N: rings at the exact radius, but not co-rotating |
| P10 chimera | T10.1 (five observables) and T10.2 (existence boundary) | FAIL | T10.1 fails on f_drift (0.752 against 0.7008 +- 0.03); T10.2 passes (median collapse at 0.2588) |

### P5 flocking (Vicsek 1995): FAIL

N = 400, L = 10, `position_update="old_velocity"` (Vicsek's Eq. 1), 10 seeds per eta. Seed means against the
digitized Fig. 2a:

| eta | 0.5 | 1.0 | 1.5 | 2.0 | 2.5 | 3.0 | 3.5 | 4.0 | 4.5 |
|---|---|---|---|---|---|---|---|---|---|
| measured v_a | 0.986 | 0.947 | 0.884 | 0.797 | 0.690 | 0.567 | 0.425 | 0.277 | 0.140 |
| Fig. 2a | 0.981 | 0.920 | 0.830 | 0.709 | 0.566 | 0.380 | 0.133 | 0.081 | 0.062 |
| difference | +0.005 | +0.027 | +0.054 | +0.088 | +0.124 | +0.187 | +0.292 | +0.196 | +0.078 |

Seed SE 0.0001 to 0.0015. The model stays more ordered than the published curve at every noise level, and the
gap is largest near the published transition (eta about 3 to 3.5).

Reported, not gating:
- N = 100 and N = 40 fail the same 0.04 rule (largest gaps 0.221 and 0.166). The ordering v(40) > v(100) > v(400)
  holds at all five eta >= 2.5.
- The catalog's own update order (`new_velocity`) at N = 400 is off by up to 0.242, so the update order does not
  explain the gap.
- Optional N = 4000 (5 seeds): 0.737, 0.601, 0.452 at eta = 2.0, 2.5, 3.0 against 0.645, 0.463, 0.142.
- beta by Vicsek's straightest-line procedure: 0.226 with eta_c = 3.00 (published 0.45 +- 0.07).
- Target check: re-digitizing Fig. 2a reproduces the frozen targets at N = 400, 100 and 4000; at N = 40 two values
  (eta = 2.5 and 3.0) differ from the frozen ones by 0.0005, a rounding difference with no effect on any rule.

Runtime 65.8 CPU-h over 375 runs.

### P6 milling (D'Orsogna et al. 2006): FAIL

alpha = 1, beta = 0.5, C_a = 1, l_a = 1, C_r = 0.5, l_r = 0.5, dt = 0.01 (the pre-registered dt check passed: changes of
9e-8 and 1.9e-5 at N = 100 and 400), 10 seeds per N.

Every seed at every N settled into a ring whose radius was identical across seeds and matched the exact ring
radius to within 0.01% (R/R_exact = 0.99993 to 0.99999) and Eq. (4) to within 1.2% (R/R_Eq4 = 1.0117, 1.0055, 1.0036,
1.0026 at N = 100, 200, 300, 400). The particles moved at the expected speed sqrt(2) with zero polarization, but they
did not all circle the same way: the normalized angular momentum |L| was 0.010 to 0.180, against the co-rotation gate
|L| > 0.95. No seed passed all gates, so the first condition of T6.1 (at least 8 of 10 valid seeds at every N) fails
at every N and T6.1 fails.

Reported, not gating:
- T6.2 (radius against the exact solution within 1%) is defined on valid seeds only and is therefore not evaluated.
  Measured on all seeds, the radius is within 0.01% of the exact value.
- T6.3 energy exponent: 2.11 for the catastrophic set (band [1.9, 2.1]) and 2.02 for the H-stable set (band
  [0.9, 1.1]).
- The init_radius sensitivity named in the spec without a protocol was not run.

Runtime 22.8 CPU-h over 82 runs.

### P10 chimera (Abrams & Strogatz 2004): FAIL

N = 256, A = 0.995, beta = 0.18, dt_paper = 0.025. The paper's initial condition gave full synchrony in all 10 seeds
(and in all 10 at N = 128), so, as pre-registered, T10.1 and T10.2 were evaluated on the 10 fallback realizations
(seeds 100 to 109), every one of which formed a chimera.

T10.1, medians over the 10 chimera realizations against the analytic solution:

| Observable | Target | Band | Measured | Within |
|---|---|---|---|---|
| r_bar | 0.6287 | 0.02 | 0.6309 | yes |
| Delta | 0.7320 | 0.010 | 0.7321 | yes |
| f_drift | 0.7008 | 0.03 | 0.7520 | no |
| R_min | 0.4732 | 0.03 | 0.4778 | yes |
| R_max | 0.8068 | 0.03 | 0.8081 | yes |

The drifting fraction was 0.738 to 0.801 across the ten realizations, all above the band. T10.1 fails.

T10.2: the chimera collapsed between beta = 0.2525 and 0.2625 in every ramp (median 0.2588, band [0.2500, 0.2675],
continuum fold 0.2635), and no ramp held a chimera at beta >= 0.27. T10.2 passes. The primary rule needs both, so
P10 fails.

Reported, not gating:
- dt check: not run (deviation P10-D1).
- T10.3 at beta = 0.05 (paper IC): 9 of 10 seeds formed chimeras. r_bar, R_min and R_max are within their bands;
  Delta splits between two groups (six seeds near 0.145, three near 0.715, against 0.716), and f_drift is 0.926
  against 0.8526.
- Optional N = 512 not run.

Runtime 34.6 CPU-h over 50 runs.
