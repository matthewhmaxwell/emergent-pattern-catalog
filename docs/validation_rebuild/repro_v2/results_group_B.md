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
