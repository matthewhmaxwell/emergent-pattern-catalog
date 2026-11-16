# Repro v2, group A results: P23 minority game, P27 Nowak-May, P18 voter

Spec: `docs/validation_rebuild/repro_v2/group_A.md` (frozen at commit 39688e3), plan `REPRO_V2_PLAN.md`.
Code: `analysis/reproductions/v2/`. Results: `analysis/outputs/repro_v2/P23.json`, `P27.json`, `P18.json`,
raw per-run records in `analysis/outputs/repro_v2/raw/`. P1 is dropped (plan) and was not run.

## Fixed before the full runs

This section was committed before any full-protocol run. It records how the code reads the parts of the
protocol that the spec leaves open. No target, tolerance, size, seed count or pass rule was changed.

### Code changes

- `epc/models/minority_game.py`: two options, `tie_break` ('first' default, 'random') and `payoff` ('sign'
  default, 'linear'). With the defaults the model is unchanged: `check_bit_identical.py` compares it with
  the frozen-commit model on the catalog setting (N = 101, m = 6, 3000 rounds, seeds 0-4) and two small cases,
  record by record. 'random' draws `rng.integers(0, 2, N)` every round and uses it where the two scores tie
  (the spec's formula). 'linear' adds `-(2 pred - 1)(2 n1 - N)`, computed in int64 because the strategy
  tables are int8 (the spec's expression, evaluated in int8, would overflow).
  `tests/test_minority_game_p23_e2e.py` gives the same result with and without this change: 10 pass,
  2 fail (`test_definitive_on_efficient_regime`, `test_confidence_at_definitive`). The two failures are
  already present at 39688e3 and are detector assertions, not model behaviour.
- No change to `nowak_may.py` or `voter.py`. P27 uses `b = 1.9`; P18 uses `neighborhood='von_neumann'`.
- Faster inner loops (`fastloops.c`, compiled at run time into a temporary directory; wrappers
  `mg_fast.py`, `nm_fast.py`, `voter_fast.py`). Each performs the model's arithmetic in the same order and
  draws its random numbers with the model's own Generator. `check_bit_identical.py` checks them against the
  model code on small cases: for the minority game, attendance, winning side, scores, history bits and the
  full rng state, all four option combinations, N = 11-511; for the voter, the grid after every sweep and
  the rng state, both neighbourhoods; for Nowak-May, the grid and f_C after every generation at
  b = 1.6, 1.8, 1.85, 1.9. All pass (`analysis/outputs/repro_v2/bit_identity_checks.json`).
- The minority-game loop draws a whole chunk's coin tosses in one call. NumPy's `integers(0, 2)` takes one
  32-bit value per coin and keeps the unused half of a 64-bit draw between calls, so the coins and the final
  rng state are identical to one call per round; the check above includes the rng state.
- Measurement is streamed: the minority game keeps per-chunk sums of A and A^2 (T1) or per-mu counts and sums
  of A plus the sign history of Delta_i (T2); Nowak-May keeps the cooperator count per generation; the voter
  keeps `wall_density` at the recording times. No grid or choice history is stored.

### Readings of the protocol

**P23.**
- Rounds are played in chunks of P. Burn-in 200P rounds, then 1000P measurement rounds (T1); 500P and 2000P
  (T2). sigma^2/N = mean over the measurement rounds of (2 n1 - N)^2 / N; the per-alpha value is the mean
  over seeds; alpha* is the grid alpha at the minimum of that curve.
- Bootstrap: 2000 resamples of the seed indices with replacement, applied jointly to all grid points of one
  P; alpha* recomputed each time; 95% CI from the 2.5th and 97.5th percentiles; rng seed 20261116.
- Stationarity check (pre-registered in the spec): at grid points with |alpha - 0.3374| < 0.03, the
  seed-mean sigma^2/N over the first 500P measurement rounds (h1) is compared with that over the second
  500P (h2). "Differ by more than 5%" is read as |h2 - h1| / min(h1, h2) > 0.05, the stricter of the
  possible denominators. Failing points are rerun once with burn-in 400P (same seeds, same 1000P
  measurement); the rerun replaces the original at those points and the check is reported again. The same
  check is applied to the tie_break='first' sensitivity grid.
- Tie-rule flag: raised if |alpha*(first, seeds 0-9) - alpha*(random, seeds 0-19)| > 0.05 at P = 128 or at
  P = 64. alpha*(random) restricted to seeds 0-9 is also reported.
- T2: Delta_i = U_i0 - U_i1 is snapshotted after the burn-in and after every P measurement rounds (2001
  snapshots). phi_hat is per run; the rule compares the seed mean (10 seeds) with phi_RS, and the seed mean
  of H_hat/N with H_RS. mu is the history index the agents used in that round, which equals the index rebuilt
  from the concatenated winning-side series including the m rounds before the window (checked).
- RS values are recomputed by `group_A_support/rs_minority_game.py` at analysis time.

**P27.**
- f_C(t) is the cooperator fraction after generation t; the per-seed estimate is the mean over t = 501..1500
  (1000 values) and the temporal SD uses ddof = 1.
- A seed is absorbing if f_C reaches 0 or 1 at any generation up to 1500. Absorbing seeds are counted and
  excluded from the mean at every L. Ensemble SE = SD/sqrt(n), ddof = 1.
- The break-point invariance check uses `NowakMayModel.step()` itself (not the fast loop) and is run before
  any gating run.

**P18.**
- rho(t) is `voter._wall_density` after t sweeps. rho_bar(t) is the seed mean; SE = SD/sqrt(n), ddof = 1.
- Precision requirement: if SE(3000) > 0.01 rho_exact(3000) with seeds 0-19, seeds 20-39 are run and every
  rule is judged on all 40 seeds (no stopping in between).
- Slope: OLS of 1/rho_bar(t) on ln t at t = 300, 500, 1000, 2000, 3000. Bootstrap over seeds, 2000
  resamples, percentile 95% CI, rng seed 20261117.
- Finite-size check: |rho_bar_256(3000) - rho_bar_512(3000)| / rho_bar_512(3000) <= 0.02.
- rho_exact is recomputed by `group_A_support/voter_exact_crw.py` at analysis time.
- Spec risk item 1 (the Sprint 20 Moore figure 0.22 at L = 128): that number came from the P18 detector's
  `_wall_density_moore`, which uses the same half-offset bond count as the model's `_wall_density` with Moore
  offsets. It is the model's quantity. The Moore run here is non-gating.

### Run plan

At most two compute processes, `nice -n 10`, tmux sessions owned by matthewhmaxwell.

- Queue 1 (`run_queue_1.sh`): P23 T1 random (800 runs), stationarity check and reruns, T1 'first' (400 runs)
  and its check, T2 linear and sign (80 runs each).
- Queue 2 (`run_queue_2.sh`): P18 theory, L = 512 von Neumann seeds 0-19, precision check (and seeds 20-39
  if needed), L = 256 von Neumann and Moore; then P27 invariance check, L = 1000, 500, 200 at b = 1.9 and
  L = 1000 at b = 1.8.
