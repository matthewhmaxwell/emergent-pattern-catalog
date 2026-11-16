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

None at the time of the pre-run commit. Any later deviation is entered here, with its reason, before the
affected rule is evaluated.
