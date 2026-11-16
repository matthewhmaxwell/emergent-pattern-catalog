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

## Run record

- All full runs used code commit `a21b04d` with a clean tracked tree (recorded in every raw record).
  The bit-identity checks were rerun at that commit and pass.
- Queue 1 ran 01:05-01:31 UTC on 2026-11-16. Queue 2 started at 01:05. A third session (`repro_a_q3`)
  ran the non-gating P27 b = 1.8 set from 01:47, after queue 1 had finished, so no more than two compute
  processes ran at any time. Queue 2 resumes per seed, so it skipped seeds that session had already written.
- No protocol deviation was needed. `deviations_group_A.json` is empty.

## P23: minority game

**Primary rule T1: FAIL.** alpha*(P = 128) = 0.3988 (N = 321). The pass window is [0.29, 0.39]; the
minimum is 0.0088 above its upper edge.

| alpha (P = 128) | 0.3395 | 0.3507 | 0.3606 | 0.3798 | **0.3988** | 0.4491 | 0.5311 |
|---|---|---|---|---|---|---|---|
| sigma^2/N, seed mean of 20 | 0.219 | 0.192 | 0.169 | 0.153 | **0.146** | 0.156 | 0.185 |
| SE | 0.025 | 0.024 | 0.017 | 0.011 | 0.003 | 0.001 | 0.002 |

- Bootstrap over seeds (2000 resamples): 95% CI for alpha* is [0.379, 0.399]. alpha* = 0.3988 in 1519
  resamples, 0.3798 in 428, 0.3606 in 44, 0.3507 in 6 and 0.4491 in 3, so 23.9% of resamples land inside
  the pass window.
- The curve is shallow on the symmetric-phase side. The SE there is 0.011-0.025 because sigma^2 varies
  strongly between seeds below and near alpha_c.

**Reported, not gating.**
- alpha*(P = 64) = 0.4476 (N = 143), bootstrap 95% CI [0.379, 0.529]. alpha* moves down from P = 64 to
  P = 128 (0.448 to 0.399). That fits a minimum drifting towards 0.337 as P grows, but two values of P do
  not test it.
- Minimum sigma^2/N: 0.146 (P = 128) and 0.168 (P = 64), against the RS value 0.114 at alpha_c, which
  is expected to be lower than any finite-N value.
- Stationarity check. 12 points were checked (6 per P). Two failed at P = 128: alpha = 0.3606 (5.5%) and
  0.3507 (5.8%). Both were rerun with burn-in 400P and then passed (3.6%, 3.0%). The rerun values are the
  ones used above.
- Tie-rule sensitivity (tie_break = 'first', seeds 0-9): alpha* = 0.3988 at P = 128 and 0.4476 at P = 64,
  the same grid points as under 'random'. **The flag is not raised** (difference 0 at both P). The 'first'
  curve is lower on the symmetric side (0.183 against 0.219 at alpha = 0.3395, P = 128), but the minimum
  does not move. Its stationarity check failed at three P = 128 points (6.0%, 6.7%, 5.2%), which passed
  after the 400P rerun (4.2%, 2.9%, 4.9%).

**T2, linear payoff, tie_break = 'random' (gating for the variant; reported only in the plan): FAIL.**

| (P, N) | alpha | phi_hat (SE) | phi_RS | abs diff (<= 0.03) | H_hat/N (SE) | H_RS | rel diff (<= 10%) |
|---|---|---|---|---|---|---|---|
| (512, 511) gating | 1.002 | 0.4354 (0.0065) | 0.4294 | 0.006 pass | 0.1288 (0.0032) | 0.1486 | -13.3% fail |
| (512, 255) gating | 2.008 | 0.2706 (0.0124) | 0.2392 | 0.031 fail | 0.2447 (0.0066) | 0.2795 | -12.5% fail |
| (256, 255) finite size | 1.004 | 0.4322 (0.0088) | 0.4288 | 0.003 | 0.1282 (0.0067) | 0.1490 | -13.9% |
| (256, 127) finite size | 2.016 | 0.2669 (0.0091) | 0.2381 | 0.029 | 0.2444 (0.0069) | 0.2802 | -12.8% |

- H_hat/N is 12.5-13.9% below the RS value at all four points. The shortfall does not shrink from N = 127
  to N = 511. Measured sigma^2/N (0.314 and 0.494 at P = 512) is also below the approximate RS column
  (0.347, 0.555), which is non-gating.
- phi_hat agrees with RS at alpha ~ 1 and is 0.03 high at alpha ~ 2 at both P.
- The half-window phi_hat values are within 0.01 of the full-window values at every point.
- Sign-payoff companion (non-gating): phi_hat = 0.444 / 0.288 at P = 512 (RS 0.429 / 0.239). H_hat/N is
  5.0% and 6.3% below RS.
- The spec lists endogenous histories (the RS theory assumes random information) as a known approximation
  for T2. These runs do not separate that effect from the others.

**Runtime.** 1,538 CPU-seconds summed over the 1,350 runs: T1 random 811 s (840 runs incl. 40 reruns),
T1 'first' 219 s (430), T2 linear 246 s, T2 sign 262 s. Wall time was 26 min.

## P18: 2-D voter model

**Primary rules: T2 PASS and T1 PASS.** von Neumann, L = 512, 20 seeds, 3000 sweeps.

Precision requirement met with 20 seeds: SE(3000) = 0.00125 <= 0.01 rho_exact(3000) = 0.00136. No extra
seeds were run.

| t (sweeps) | rho_bar (SE) | rho_exact | rel diff | T2 (<= 3%) |
|---|---|---|---|---|
| 100 | 0.19080 (0.00040) | 0.19153 | -0.38% | pass |
| 300 | 0.16934 (0.00046) | 0.16954 | -0.11% | pass |
| 500 | 0.16086 (0.00062) | 0.16089 | -0.02% | (not gated) |
| 1000 | 0.15066 (0.00071) | 0.15044 | +0.14% | pass |
| 2000 | 0.14155 (0.00095) | 0.14124 | +0.22% | (not gated) |
| 3000 | 0.13644 (0.00125) | 0.13635 | +0.06% | pass |

- T1: OLS slope of 1/rho_bar on ln t over t = 300-3000 = **0.6165**, bootstrap 95% CI [0.568, 0.663]
  (SD 0.024). The pass window is [0.586, 0.688]. The exact finite-t slope over the same times is 0.6234,
  and the asymptote is 2/pi = 0.6366.
- Non-gated early times: rho_bar(10) = 0.2605, rho_bar(30) = 0.2231.
- The exact values were regenerated by `voter_exact_crw.py` and match the frozen spec table to 5 digits.
  The script took 1,598 s on the loaded VPS (spec: about 80 s).

**Reported, not gating.**
- Finite-size check: **FAIL.** rho_bar(3000) = 0.13180 (SE 0.00263) at L = 256 against 0.13644 at
  L = 512, a 3.4% difference (rule: <= 2%). The difference is 1.6 combined SE, and the L = 256 SE alone is
  1.9% of rho, so with 10 seeds the check cannot resolve 2%. The L = 256 values fall below the exact curve
  from t = 500 on (-1.9%, -1.5%, -2.4%, -3.3% at t = 500, 1000, 2000, 3000). The mean squared magnetization
  at t = 3000 is 0.054 at L = 256 and 0.005 at L = 512. The target formula (Dornic et al. Eq. 1) carries a
  factor (1 - m^2). That factor is about 5% lower at L = 256 and 0.5% lower at L = 512, the same order as
  the observed 3.4% gap. This was not tested further.
- L = 256 von Neumann slope: 0.711, CI [0.588, 0.827].
- Moore, L = 256, 10 seeds, compared with the exact Moore column: +0.16%, -1.4%, +0.12%, -3.1%, +0.5% and
  +1.9% at t = 100, 300, 500, 1000, 2000, 3000. The slope is 0.361, CI [0.266, 0.471], against the exact
  finite-t value 0.418 and the derived asymptote 4/(3 pi) = 0.424. The CI is wide with 10 seeds.
- Spec risk item 1 (catalog Moore figure 0.22 at L = 128, t in [600, 800]): the Moore model output here
  matches the exact Moore solution at L = 256 to within about 3%. The catalog figure was not rerun.

**Runtime.** 1,213 CPU-seconds of simulation (L = 512: 940 s; L = 256 von Neumann: 130 s; Moore: 143 s)
plus 1,598 s for the exact solution.

## P27: Nowak-May spatial prisoner's dilemma

**Break-point invariance check: PASS** (run first, at 01:52, before any gating run). At L = 200, seed 0,
`NowakMayModel.step()` gave grids at b = 1.85 and b = 1.95 that are bit-identical to b = 1.9 at all
300 generations.

**Primary rule: PASS** (all three conditions). b = 1.9, random start, periodic, 1500 generations, f_C
averaged over generations 501-1500, seeds 0-9.

| Set | mean f_C (SE) | SD across seeds | mean temporal SD | absorbing seeds |
|---|---|---|---|---|
| L = 1000, b = 1.9 (gating) | 0.31559 (0.00003) | 0.00008 | 0.0017 | 0 of 10 |
| L = 500, b = 1.9 | 0.31549 (0.00005) | 0.00015 | 0.0034 | 0 of 10 |
| L = 200, b = 1.9 (reported) | 0.31535 (0.00010) | 0.00032 | 0.0084 | 0 of 10 |
| L = 1000, b = 1.8 (break point, reported) | 0.40797 (0.00002) | 0.00006 | 0.0015 | 0 of 10 |

1. |0.31559 - 0.316| = 0.0004 <= 0.010: pass.
2. |0.31549 - 0.31559| = 0.0001 <= 0.005: pass.
3. 0 of 10 seeds absorbing at L = 1000: pass.

**Reported, not gating.**
- The gating mean is 0.0022 below 12 ln 2 - 8 = 0.31777. With an SE of 0.00003, the two values are
  clearly distinct. The measured value agrees with the Fig. 26 reading (0.3159) to 0.0003, not with the
  analytic value attributed to Nowak & May 1992. The spec already flags that attribution as unverified.
- At the break point b = 1.8, f_C = 0.408 at L = 1000. This matches the catalog's 0.408 at L = 100 and is
  0.092 above the 9/5 < b < 2 value. b = 1.8 is the break point where, as the spec notes, ties are
  resolved by scan order.
- f_C barely depends on L (0.3154-0.3156 from L = 200 to 1000). The temporal SD scales roughly as 1/L.

**Runtime.** 2,784 CPU-seconds: L = 1000 b = 1.9 1,191 s, L = 500 298 s, L = 200 47 s, b = 1.8 1,236 s,
invariance check 12 s.

## Summary

| Pattern | Primary rule | Measured vs target | Reported-only results |
|---|---|---|---|
| P23 minority game | T1 **FAIL** | alpha*(P = 128) = 0.3988 vs window [0.29, 0.39] (alpha_c = 0.3374); bootstrap CI [0.379, 0.399] | alpha*(P = 64) = 0.448; tie-rule flag not raised; T2 (linear-payoff variant) FAIL: H/N 12.5-13.3% below RS, phi 0.031 off at alpha ~ 2 |
| P27 Nowak-May | **PASS** | f_C(L = 1000) = 0.3156 vs 0.316 +/- 0.010; L = 500 within 0.0001; 0 absorbing | 0.0022 below 12 ln 2 - 8; b = 1.8 gives 0.408 |
| P18 voter | T2 **PASS**, T1 **PASS** | rho within 0.4% of exact at t = 100-3000; slope 0.6165 vs [0.586, 0.688] | finite-size check FAIL (3.4% vs 2%); Moore within 3.1% of exact |

Following the plan, P23 keeps its previous grade, and P27 and P18 count as "published value"
reproductions under their primary rules.
