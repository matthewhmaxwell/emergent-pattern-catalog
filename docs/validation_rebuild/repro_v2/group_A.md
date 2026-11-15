# Group A literature-reproduction specs: P23, P27, P18, P1

Frozen 2026-11-15, before any simulation is run. No model or simulation code was executed to prepare this
document. Model files were read only:

- `sources/code/epc/models/minority_game.py`
- `nowak_may.py`
- `voter.py`
- `schelling.py`
- `sources/code/run_phase2a_panel.py`, the `build_p{1,18,23,27}_positives` functions

Two analytic predictions were evaluated numerically from the sources' own equations. Neither is a run of our
models.

- **Minority game.** The replica-symmetric (RS) saddle-point equations, from Challet, Marsili & Zecchina (CMZ)
  2000 and Marsili, Challet & Zecchina (MCZ) 2000. Solving them reproduces the published alpha_c = 0.33740 to six
  digits, which validates the transcription.
- **Voter model.** The exact coalescing-random-walk solution of Frachebourg & Krapivsky (FK) 1996, Eqs. (9)-(11).

Scripts are in `group_A_support/`: `rs_minority_game.py` (under 1 s) and `voter_exact_crw.py` (about 80 s).

**Figure readings.** Two values were read off published figures:

- Szabó & Fáth Fig. 26, by pixel centroid at 600 dpi.
- Gauvin et al. Fig. 3, by eye at 400 dpi.

Both are labelled as readings, with a reading uncertainty.

**Quotes.** Each section gives one verbatim quote, under 15 words because of the copyright policy this work runs
under. Exact locators (page, equation, figure) are given instead of longer excerpts. Page numbers refer to the
arXiv versions cited.

**Standard.** A reproduction counts only if a specific published number, or an analytic result stated in the source,
is compared with model output under the pass rule fixed here.

| Pattern | Primary frozen target | Code change needed | Verdict |
|---|---|---|---|
| P23 minority game | alpha_c = 0.3374 as the location of the sigma^2/N minimum (CMZ 2000). Exact RS frozen fraction and predictability at alpha ~ 1 and 2. | Random tie-breaking (all targets). A `payoff='linear'` option (T2 only). | **Target found, feasible.** |
| P27 Nowak-May | Cooperator frequency 0.316 +/- 0.005 for 9/5 < b < 2, random start, read from Szabó & Fáth 2007 Fig. 26. The analytic 12 ln 2 - 8 = 0.3178 is attributed to Nowak & May 1992; the primary text was not accessible. | None. Parameter changes only: **b must move off 1.8**, which is a break point. | **Target found, feasible.** The figure reading is the gating number. Verify the Nature text before citing 12 ln 2 - 8. |
| P18 voter | Exact rho(t) of the n.n. voter model at t = 100-3000 (FK 1996 Eqs. 9-11). Published asymptotic slope d(1/rho)/d ln t = 1/(2 pi D) = 2/pi (Dornic et al. 2001 Eq. 1). | None. Use `neighborhood='von_neumann'`. | **Target found, feasible.** FK's closed form, Eq. (23), is 7-11% off at these times; do not use it as a gate. The voter_exact_crw.py rerun reproduced the frozen table. |
| P1 Schelling | Gauvin et al. 2009 transitions at discrete tolerance values: frozen-to-segregated T_f(rho) (Table 1) and segregated-to-mixed between T = 5/7 and 3/4 (Fig. 3). | A new dynamics mode, because Gauvin's rules differ from ours. | **Weak.** No published number exists for the catalog's rule. Either reproduce the Gauvin variant as a variant, or drop P1's quantitative reproduction claim. |

---

## P23: Anti-coordination, the minority game

### 1. Published targets

**T1. alpha_c, the location of the sigma^2/N minimum (gating).**

**Primary source.** Challet D, Marsili M, Zecchina R (2000). *Statistical mechanics of systems with heterogeneous
agents: Minority Games.* Phys. Rev. Lett. 84:1824. doi:10.1103/PhysRevLett.84.1824. arXiv cond-mat/9904392v2.

- **p. 3, after Eq. (7):** alpha_c = 0.33740..., exact within the replica-symmetric solution for S = 2. Quote:
  "For α > αc = 0.33740 . . ., the solution leads to Q = q < 1".
- **p. 2:** sigma^2 "shows a cusp" at alpha_c ~ 0.34. Fig. 1 shows sigma^2/N vs alpha for P = 2^6, averaged over
  200 realizations.
- **Scope.** The result is for the linear payoff U_s += -a_s A/P. The same page (p. 2) says the original sign payoff
  leads "to similar conclusions".

**Same rule as ours (sign payoff, real histories).** Challet D, Marsili M (1999). *Phase transition and symmetry
breaking in the minority game.* Phys. Rev. E 60:R6271. doi:10.1103/PhysRevE.60.R6271. arXiv cond-mat/9904071v2.

- Sign payoff, Eq. (4), with endogenous histories.
- p. 2: ties are "broken by coin tossing".
- p. 3 and Fig. 1: a singularity at alpha_c ~ 0.34 in the order parameter theta and the frozen fraction phi, and
  the sigma^2 minimum. Shown for P = 32, 64 and 128. The value is numerical, with no error bar.
- Same paper, p. 3: phi ~ e^{-(0.37 +/- 0.02) alpha} at large alpha. This is not used as a gate, because the fit
  range is unstated.

**Corroboration.** Granath M, Perez-Diaz A (2016). J. Stat. Phys. 165:94-114. doi:10.1007/s10955-016-1607-8.
arXiv 1509.08236.

- Studies the sign-payoff game, with ties decided by coin toss (p. 2).
- Gives alpha_c = 0.3374 (p. 1).
- Appendix C (p. 11): frozen fractions are "quite insensitive" to sign vs linear payoff.

**Not used.** Savit, Manuca & Riolo (1999), PRL 82:2203, doi:10.1103/PhysRevLett.82.2203, which the catalog cites.

- The preprint (adap-org/9712006v1, p. 8) gives z_c ~ 0.5 from integer m.
- The long version (Manuca, Li, Riolo & Savit, Physica A 282:559 (2000); adap-org/9811005, footnote 16, PDF p. 38)
  says non-integer m points to a value "closer to 1/3".
- Neither gives a number with an uncertainty. CMZ 2000 supersedes both.

**T2. Exact RS predictions in the asymmetric phase (gating for the linear-payoff variant).**

Equations:

- CMZ 2000 Eq. (9), p. 3:
  - z = sqrt(alpha/(1+Q)).
  - Frozen weight per sign (1 - erf(z/sqrt2))/2. The total frozen fraction is therefore phi = 1 - erf(z/sqrt2).
  - The text calls the per-sign weight "the fraction of frozen agents", but the normalization of Eq. (9) makes
    the total twice that.
- MCZ 2000, Physica A 280:522, doi:10.1016/S0378-4371(99)00610-X, arXiv cond-mat/9908480v3, Eqs. (A.12)-(A.14),
  pp. 34-35:
  - H/N = G/(1+chi)^2.
  - For S = 2, G = (1+Q)/2, which gives chi = erf(z/sqrt2)/(alpha - erf(z/sqrt2)).
  - Q(z) = 1 - sqrt(2/pi) e^{-z^2/2}/z - (1 - 1/z^2) erf(z/sqrt2), and alpha = z^2 (1+Q).
- Check: setting chi -> infinity gives alpha_c = 0.337400.

Frozen values, in literature units (a_i = +/-1, A = sum a_i, random agents give sigma^2/N = 1 and H/N = 1/2):

| P | N | alpha | phi (total frozen) | H/N | Q | sigma^2/N (approx., non-gating) |
|---|---|---|---|---|---|---|
| 512 | 511 | 1.00196 | 0.4294 | 0.1486 | 0.6042 | 0.3465 |
| 512 | 255 | 2.00784 | 0.2392 | 0.2795 | 0.4493 | 0.5549 |
| 256 | 255 | 1.00392 | 0.4288 | 0.1490 | 0.6038 | 0.3471 |
| 256 | 127 | 2.01575 | 0.2381 | 0.2802 | 0.4483 | 0.5560 |

The sigma^2 column uses sigma^2/N ~ H/N + (1-Q)/2 (CMZ Eq. 5, assuming independent s_i). It is an approximation,
so it is reported only. Sanity check: at alpha = 1.5 it gives 0.468. That lies inside the 0.45-0.49 axis range of
Marsili & Challet 2001 (PRE 64:056138, arXiv cond-mat/0102257, Fig. 2 bottom, p. 13), which also shows sigma^2/N
independent of Gamma above alpha_c.

**Units warning.** Our attendance n1 runs from 0 to N. Convert with A = 2 n1 - N, so the literature quantity is
sigma^2/N = <A^2>/N = 4 <(n1 - N/2)^2>/N. The catalog's dim1 minimum of 0.077 (detector units, var(n1)/N) is about
0.31 in literature units. The RS value at alpha_c is (1 - Q_c)/2 = 0.114.

### 2. Model match

| Item | Sources | Ours (`minority_game.py`) | Fix |
|---|---|---|---|
| Agents and strategies | N odd, S = 2, each strategy a uniform random table over P = 2^m histories | Same | None |
| Information | Real history of the last m minority signs (CM 1999; GP 2016). The CMZ theory assumes random mu (Cavagna 1999); MCZ (p. 29) finds De Bruijn effects not significant. | Real history (`history_bits`) | None. For T2 this is a known approximation; report it. |
| Score update | Sign payoff U_s += -a_s sign(A) (CM 1999 Eq. 4). Linear payoff U_s += -a_s A/P (CMZ Eq. 4). | +1 if the prediction equals the winning side, else -1. This is the sign payoff. | T1: none. T2: **code**, add `payoff='linear'` with `scores[:, s] += -(2*pred-1)*(2*attendance-N)`. The 1/P factor is irrelevant at Gamma = infinity. |
| Choice | Best strategy, Gamma = infinity (CM 1999; MCZ Fig. 2) | `argmax(scores)` | None |
| **Ties** | Coin toss (CM 1999 p. 2; GP 2016 p. 2; Savit et al. preprint adap-org/9712006 p. 5). Manuca et al. 2000 defaulted to "most recently played" (PDF p. 32). | `np.argmax` always returns strategy 0. **No source uses this rule.** | **Code:** add `tie_break='random'`, e.g. `np.where(s0 == s1, rng.integers(0, 2, N), argmax)`. Keep the default 'first' so catalog runs stay reproducible. |
| Initial scores | U = 0 (tabula rasa) | Zeros | None |
| Initial history | Arbitrary | Random m bits | None |
| N parity | Odd | Enforced | Grid uses odd N |
| alpha resolution | CMZ fixes P and varies N | Integer m only in the catalog runs | Parameter: fix P, vary odd N |

### 3. Protocol

**T1 (gating): sign payoff, `tie_break='random'`, S = 2.**

- **Sizes.**
  - P = 64, N in {427, 355, 305, 267, 247, 229, 213, 207, 201, 193, 189, 183, 177, 169, 161, 143, 121, 101, 81, 65}.
    This gives alpha from 0.150 to 0.985, with spacing about 0.01 near 0.34.
  - P = 128, N in {853, 711, 609, 533, 493, 457, 427, 413, 401, 387, 377, 365, 355, 337, 321, 285, 241, 203, 161, 129}.
- **Seeds.** 0-19 at every grid point.
- **Run length.** Burn-in 200P rounds, then 1000P measurement rounds.
- **Chunking.** Run in chunks of P rounds by setting `model.params.n_rounds = P` and calling `simulate()` repeatedly.
  State persists across calls. Keep only `attendance` and `winning_side`; discard the `choices` copies.
- **Estimator.**
  - Per seed: sigma^2/N = mean over the measurement rounds of (2 n1 - N)^2 / N.
  - Per alpha: the seed mean.
  - alpha* = argmin of the seed-mean curve on the grid.
  - 95% CI by bootstrap over seeds (2000 resamples, argmin recomputed each time).
- **Stationarity check (pre-registered).** At grid points with |alpha - 0.3374| < 0.03, compare the first and second
  halves of the measurement window. If they differ by more than 5%, double the burn-in at those points and rerun.
- **Sensitivity (non-gating).** Repeat the full grid with `tie_break='first'` (current catalog code), seeds 0-9.

**T2 (gating for the variant): `payoff='linear'`, `tie_break='random'`.**

- **Sizes.** (P, N) = (512, 511) and (512, 255) are gating. (256, 255) and (256, 127) are a finite-size check.
- **Seeds.** 0-9 per point.
- **Run length.** Burn-in 500P rounds, measurement 2000P rounds. Snapshot `Delta_i = scores[:, 0] - scores[:, 1]`
  every P rounds.
- **phi_hat.** The fraction of agents whose Delta_i has the same nonzero sign at all 2001 snapshots and satisfies
  |Delta_i(end)| > |Delta_i(start)|. Sensitivity: the same rule on each half-window.
- **H_hat/N.** (1/N) sum_mu f_mu <A | mu>^2. Here mu_t is rebuilt from the concatenated `winning_side` series,
  including the m rounds before the window, and f_mu is the empirical frequency of mu.
- **sigma_hat^2/N.** <A^2>/N (non-gating).
- **Sign-payoff companion (non-gating).** Run the same points with `payoff='sign'` and report phi_hat against the RS
  phi. GP 2016 says these should be close.

### 4. Pass rules (fixed)

**T1 passes** if alpha*(P = 128) lies in [0.29, 0.39], that is |alpha* - 0.3374| <= 0.05.

- **Source uncertainty.** The analytic value is exact. The sign-payoff numerics are quoted only to two significant
  figures (about 0.34).
- **Finite-size allowance.** At P = 64-128 (N about 190-380 near alpha_c), the cusp is rounded and the minimum can
  shift. +/-0.05 is about 15%, which is 3 grid spacings at P = 128.
- **Reported, not gated.** alpha*(P = 64) and the bootstrap CI. Also the minimum value against the RS 0.114, which
  is expected to be higher at finite N.
- **Flag.** If alpha* under 'first' differs from alpha* under 'random' by more than 0.05, the catalog's tie rule
  changes the physics. The P23 catalog entries must then be regenerated.

**T2 passes** if, at P = 512 and for both alpha ~ 1.002 and alpha ~ 2.008:

- |phi_hat - phi_RS| <= 0.03 (absolute), and
- |H_hat/N - H_RS| <= 0.10 H_RS.

Why these tolerances:

- The RS values are exact for N -> infinity with random information.
- The allowances cover three effects:
  - finite N of 255-511 (GP 2016 Figs. 2 and 8 show agreement with the replica phi for N = 101-2001);
  - endogenous histories;
  - snapshot-based frozen classification.
- Estimation bias in H from finite measurement length is about sigma^2/(2000N), which is negligible.
- These tolerances are author-set allowances informed by GP 2016. The source gives no finite-N error.

### 5. Compute and risks

**Compute.** Estimated from the code: about 25-40 µs per round, which is numpy overhead at N <= 1000.

| Run set | Runs | Cost |
|---|---|---|
| T1 | 20 alpha x 20 seeds x 2 P, 1200P rounds each | about 1 CPU-h |
| Tie-rule sensitivity | | about 0.4 CPU-h |
| T2 | 4 points x 10 seeds x 2 payoffs; P = 512 takes about 1.28M rounds, 45-60 s per run | about 1 CPU-h |
| **Total** | | **about 2.5 CPU-h** |

**Risks unrelated to model correctness.**

1. **Memory.** `simulate()` stores a dict and a `choices` copy every round. 1.3M rounds would be 0.5-1 GB.
   Mitigation: run in chunks, as specified above.
2. **Critical slowing near alpha_c.** Relaxation time diverges as alpha -> alpha_c from above. The stationarity
   check handles this.
3. **Frozen classification.** Agents near the frozen/fickle boundary have small drift, and the snapshot rule can
   misclassify them. The half-window sensitivity bounds this.
4. **Metric mismatch.** The P23 detector uses var(n1) about the sample mean, not <A^2> about N/2. Use the
   literature estimator defined above.
5. **Symmetric-phase dependence.** For alpha < alpha_c, sigma^2 depends on initial scores and on Gamma (CMZ p. 3;
   Marsili & Challet 2001). T1 uses tabula rasa, as the sources do. Do not mix initial conditions across runs.

---

## P27: Spatial reciprocity, Nowak & May spatial prisoner's dilemma

### 1. Published targets

**Background.** With payoffs R = 1, T = b, S = P = 0, Moore neighbourhood plus self-interaction and synchronous
best-neighbour imitation, the dynamics depend on b only through break points b = k2/k1.

- Szabó & Fáth 2007, p. 60: k1 = 2..z and k2 = 3..z+1, with z = 8 when self-interaction is included.
- Within each interval between break points the cellular automaton is identical.
- The spatially chaotic coexistence regime is 9/5 < b < 2 (Szabó & Fáth p. 61; Moskalenko et al. 2021,
  J. Phys. Conf. Ser. 1740:012057, p. 6).

**T1 (gating, figure reading).** Szabó G, Fáth G (2007). *Evolutionary games on graphs.* Phys. Rep. 446:97-216.
doi:10.1016/j.physrep.2007.04.004. arXiv cond-mat/0607344v3.

- **Fig. 26, p. 62.** Closed diamonds show Moore neighbourhood with self-interaction, random initial state with
  equal C and D, and about 10^6 sites. The cooperator frequency is flat across 9/5 < b < 2 at **rho_C = 0.316**.
- **Reading.** Pixel centroid of the filled-diamond row at 600 dpi = 0.3159.
  - 1 px = 0.00075.
  - Axis calibration uses the 0.0-1.0 ticks, 1337.5 px per unit.
  - Assigned reading uncertainty: **+/-0.005**.
- **Text, p. 61.** In this region the concentration does not depend on the starting state. Quote: "Within this
  region the average concentration of strategies becomes independent of the initial configuration".
- **Not stated** for Fig. 26: boundary conditions (Fig. 23 uses periodic) and run length.

**T1-analytic (not gating until the primary is verified): 12 ln 2 - 8 = 0.31777.**

- **Attribution.** Nowak MA, May RM (1992). *Evolutionary games and spatial chaos.* Nature 359:826-829.
  doi:10.1038/359826a0.
- **The Nature text could not be accessed** (paywall). A copy hosted on a third-party GitHub repo was not fetched,
  because it is an untrusted source.
- **Secondary attestation.**
  - Wu, Fu & Wang (2022), arXiv 2209.08267, p. 9: the value is an estimate made "on the basis of a very large
    symmetric pattern", meaning the single-defector kaleidoscope.
  - Huberman & Glance (1993), PNAS 90:7716, doi:10.1073/pnas.90.16.7716, arXiv chao-dyn/9307017 p. 2: Nowak & May
    found the cooperator frequency reaches "the same proportion" regardless of initial conditions. No number is
    given.
- **Important caveat.** Wu et al. ran the original 99x99 fixed-boundary single-defector case. It enters a period-4
  cycle with f_C = 96/9801 = 0.0098 after about 1.0 x 10^9 generations. So 0.318 is not the long-time value of a
  finite single-defector lattice.
- **Before citing 0.3178,** obtain the Nature PDF through library access. Record the exact sentence, the figure or
  page, the lattice, the boundary and the initial condition. T1 does not depend on this.

### 2. Model match

| Item | Sources | Ours (`nowak_may.py`) | Fix |
|---|---|---|---|
| Payoffs | R = 1, T = b, S = P = 0 ("weak" PD) | Same | None |
| Games | 8 Moore neighbours plus self (9 games) | `coop_neighbors` sums the 3x3 block including self. C earns 1 per C; D earns b per C. | None |
| Update | Synchronous; adopt the strategy of the highest scorer in the 3x3 block including self | Same | None |
| **b** | Strictly inside (9/5, 2) | Catalog and panel use **b = 1.8 = 9/5, a break point**. In float64, 5*1.8 == 9.0 exactly, so a D with 5 C neighbours ties a C with 8 C neighbours plus self. Ties resolve by strict '>' in scan order (self first, then (-1,-1), (-1,0), ...), which is anisotropic and matches neither adjacent interval. | **Parameter:** b = 1.9. No ties are possible: 1.9k is non-integer for k = 1..8. |
| Boundary | Fig. 23 periodic; Fig. 26 unstated | `boundary='periodic'` | Use periodic. Do not use `'fixed'`: it pads the border with cooperators that add payoff, which differs from Nowak & May's fixed boundary. |
| Initial state | Random, rho_C(0) = 0.5 | `init_mode='random'`, `init_coop_fraction=0.5` gives C with probability 0.5 | None |
| Size | About 10^6 sites | Any | Parameter: L = 1000 |
| Memory | | `run()` keeps every grid, which would be 1.5 GB at L = 1000 and 1500 steps | Use a `setup()`/`step()` loop and keep only `coop_fraction` |

### 3. Protocol

- **Gating runs.**
  - `NowakMayModel(rows=L, cols=L, b=1.9, init_mode='random', init_coop_fraction=0.5, boundary='periodic', seed=s)`.
  - L = 1000 (gating), plus L = 500 and L = 200 (finite-size check).
  - Seeds 0-9 for each L.
- **Run length.** 1500 generations. Burn-in is generations 1-500. Measure over 501-1500.
- **Estimator.**
  - Per seed: the time-mean of f_C.
  - Ensemble: the seed mean +/- SE.
  - Report the per-seed temporal SD.
- **Absorbing outcomes.** Seeds that reach f_C = 0 or 1 are counted and excluded from the mean. Fig. 26 plots
  "non-vanishing" frequency.
- **Break-point invariance (implementation check, pre-registered).** At L = 200, seed 0, 300 generations, runs at
  b = 1.85 and b = 1.95 must give grids bit-identical to b = 1.9 at every generation. This follows from the
  break-point structure and costs seconds.
- **Non-gating.**
  - b = 1.8 at L = 1000, seeds 0-9, to document the break-point behaviour. The catalog reported f_C = 0.408 at
    L = 100.
  - A comparison of the gating mean with 0.3178.

### 4. Pass rules (fixed)

**Pass** if all three hold:

1. |mean f_C(L = 1000) - 0.316| <= 0.010.
2. |mean f_C(L = 500) - mean f_C(L = 1000)| <= 0.005.
3. At most 1 of the 10 seeds at L = 1000 is absorbing.

**The invariance check must pass before the gate is evaluated.** A failure means an implementation bug, not a
reproduction failure.

**Why +/-0.010.**

- Reading uncertainty: +/-0.005.
- An equal allowance for unstated protocol details: Szabó & Fáth's run length, boundary and transient handling.
- The window includes 0.3178.
- Expected statistical error at L = 1000 is below 0.001.

### 5. Compute and risks

**Compute.** Vectorized numpy, roughly 0.1-0.2 s per generation at L = 1000, including Moran's I in
`_state_dict`.

- 1500 generations takes about 3-5 min per run.
- All runs total about 1 CPU-h.

**Risks.**

1. **Break points.** Any b in {1.5, 1.8, 2.0} sits on a break point, and the catalog's regime table uses all three.
   Only use values strictly inside intervals.
2. **Memory** if `run()` is used.
3. **Finite-lattice periodic cycles** (Wu et al.) could appear at small L. This is the reason for L >= 500 and the
   absorbing-seed count.
4. **Protocol gaps.** Szabó & Fáth's protocol is partly unstated.
5. **The 12 ln 2 - 8 value is unverified** against the primary text.

---

## P18: Consensus and coarsening, the 2D voter model

### 1. Published targets

**T1 (gating, published analytic): the asymptotic slope.**

- **Source.** Dornic I, Chaté H, Chave J, Hinrichsen H (2001). *Critical coarsening without surface tension: the
  universality class of the voter model.* Phys. Rev. Lett. 87:045701. doi:10.1103/PhysRevLett.87.045701.
  arXiv cond-mat/0101202v1.
- **Eq. (1), p. 1:** rho_m(t) = (1 - m^2) [2 pi D / ln t + O(1/ln^2 t)].
  - rho is the fraction of +- nearest-neighbour pairs.
  - The initial condition is random with magnetization m.
  - D is the diffusion constant of the underlying random walk. Quote: "D = 1/4 for the standard case of n.n.,
    square lattice walks". The same sentence specifies that each spin is updated once per unit time on average.
- **Consequence for m = 0:** d(1/rho)/d ln t -> 1/(2 pi D) = **2/pi = 0.6366**.

**T2 (gating, exact solution of the published equations).**

- **Source.** Frachebourg L, Krapivsky PL (1996). *Exact results for kinetics of catalytic reactions.* Phys. Rev. E
  53:R3009. doi:10.1103/PhysRevE.53.R3009. arXiv cond-mat/9508123v1.
- **Eqs. (9)-(11), p. 2.** The two-point function obeys a lattice diffusion equation with R_0 = 1 and an
  uncorrelated initial state.
- **Time unit.** FK set tau = 4/D(dimension) (p. 2), which in 2D is the standard rate-1 voter. This matches
  Dornic's unit: one Monte Carlo sweep.
- **Evaluation.** The equations were evaluated exactly via the equivalent coalescing-random-walk survival
  probability (`group_A_support/voter_exact_crw.py`).

Exact infinite-lattice values:

| t (sweeps) | rho exact, von Neumann | FK Eq. (23) pi/(2 ln t + ln 256) | Eq. (23) error | rho exact, Moore (non-gating) |
|---|---|---|---|---|
| 100 | 0.19153 | 0.21291 | +11.2% | 0.24385 |
| 300 | 0.16954 | 0.18531 | +9.3% | 0.21946 |
| 500 | 0.16089 | 0.17478 | +8.6% | 0.20967 |
| 1000 | 0.15044 | 0.16227 | +7.9% | 0.19768 |
| 2000 | 0.14124 | 0.15142 | +7.2% | 0.18696 |
| 3000 | 0.13635 | 0.14573 | +6.9% | 0.18120 |

**FK's published closed form, Eq. (23) on p. 3, is not usable as a tight gate.** It leaves out a 2 gamma term,
where gamma is the Euler constant. The corrected form pi/(2 ln t + ln 256 + 2 gamma) lies within 1.4-3.1% of the
exact values. Eq. (23) is reported only.

**Exact slopes** over t = 300-3000 (OLS of 1/rho on ln t at t = 300, 500, 1000, 2000, 3000):

- von Neumann: 0.6234, which is -2.1% from 2/pi.
- Moore: 0.4179, which is -1.5% from 4/(3 pi) = 0.4244.

**The Moore asymptote is a derivation.** D = 3/8: a lineage jumps to one of 8 Moore neighbours, so <dx^2> = 3/4
per jump and 2D = 3/4. Dornic states D only for the n.n. case. The Moore results are therefore non-gating.

**Rejected target.** Consensus time T ~ N ln N (Cox 1989). Prefactors converge logarithmically and the distribution
is broad. The catalog's own attempt found tau/(L^2 ln L) drifting from 0.40 to 0.24 over L = 16-32. No clean
published number with tolerance was found.

### 2. Model match

| Item | Sources | Ours (`voter.py`) | Fix |
|---|---|---|---|
| Neighbourhood | n.n. square lattice, D = 1/4 | Default `'moore'`; `'von_neumann'` available | **Parameter:** `neighborhood='von_neumann'` for the gate |
| Dynamics | A random site copies a random n.n. Unit time = one update per site on average. | One step = N random-site updates (with replacement), applied sequentially. This approximates the Poisson clock. | None. The discrepancy is O(1/N). |
| Initial state | Uncorrelated, m = 0 | `rng.random < 0.5`, so m ~ 0 +/- N^{-1/2} | None. The (1 - m^2) factor is 1 - O(1e-6). |
| Observable | Fraction of +- n.n. pairs | `_wall_density` over half-offsets, the fraction of disagreeing bonds. Moore averages n.n. and diagonal bonds. | None. Matches the exact calculation for both neighbourhoods. |
| Geometry | Infinite lattice | Periodic L x L | L >> sqrt(t). Torus images are negligible for L/sqrt(t) > 4. |
| Memory | | `run()` keeps every grid, about 0.8 GB at L = 512 and 3000 sweeps | Use a `setup()`/`step()` loop and keep `wall_density` at the recording times |

### 3. Protocol

- **Gating runs.**
  - `VoterModel(rows=512, cols=512, init_mode='random', neighborhood='von_neumann', boundary='periodic', seed=s)`.
  - Seeds 0-19.
  - 3000 sweeps.
  - Record `wall_density` at t in {100, 300, 500, 1000, 2000, 3000}. Also at 10 and 30, non-gating.
- **Finite-size check.** The same at L = 256, seeds 0-9.
- **Secondary (non-gating).** `neighborhood='moore'`, L = 256, seeds 0-9, compared with the Moore column above.
- **Estimators.**
  - rho_bar(t) is the seed mean, with SE = SD/sqrt(n).
  - The slope is the OLS of 1/rho_bar(t) on ln t over t = 300, 500, 1000, 2000, 3000. Its CI comes from a
    bootstrap over seeds (2000 resamples).

### 4. Pass rules (fixed)

**Precision requirement.** SE(t = 3000) <= 0.01 rho_exact. If this is not met, add seeds up to 40 before judging.

**T2 passes** if, for every t in {100, 300, 1000, 3000}, |rho_bar(t) - rho_exact(t)| <= 0.03 rho_exact(t), with
von Neumann at L = 512. Residual allowances:

- statistical, about 1% at 2 SE;
- Poisson vs. discrete sweeps, < 0.5% for t >= 100;
- finite L, < 0.1%.

**T1 passes** if the fitted slope lies in [0.586, 0.688], that is 2/pi +/- 8%. The allowance covers:

- the -2.1% finite-t offset of the exact solution;
- about +/-3% statistical error, estimated from 20 seeds at L = 512.

**Finite-size check passes** if rho_bar(3000) at L = 256 and L = 512 differ by <= 2%.

### 5. Compute and risks

**Compute.** `step()` is a pure-Python loop over N elementary updates. Estimated at 1.5-2.5 µs each, so 0.4-0.65 s
per sweep at L = 512.

| Run set | Cost |
|---|---|
| Gating set (20 runs x 3000 sweeps) | about 7-11 CPU-h |
| L = 256 von Neumann | about 1-1.5 CPU-h |
| Moore L = 256 | about 1-1.5 CPU-h |
| **Total** | **about 10-14 CPU-h.** Run in parallel on the VPS. |

**Risks.**

1. **Possible excess in an earlier catalog figure.** The Sprint 20 note reports wall density 0.22 +/- 0.01 for Moore
   at L = 128, averaged over t in [600, 800]. The exact Moore value at t ~ 700 is about 0.204, interpolated in ln t,
   so the catalog figure is about 8% higher. Before running, confirm whether that figure came from `_wall_density`
   or from the detector's own computation. If it is the model's, T2 will likely fail and the cause must be found.
2. **Memory** if `run()` is used.
3. **Long runtime** in pure Python. Do not "speed up" with the checkerboard variant the docstring warns against.
4. **Accidental Moore runs.** Using Moore, the catalog default, for the gate would compare against a derived D.

---

## P1: Similarity-driven aggregation, the Schelling model

### 1. Published target and its limits

**Search result.** No published quantitative result was found for the catalog's exact rule. Sources checked:

- Schelling 1971
- Vinković & Kirman 2006
- Gauvin et al. 2009
- Dall'Asta, Castellano & Marsili 2008 (J. Stat. Mech. L07002, arXiv 0707.1681)
- Singh, Vainchtein & Weiss 2009 (Demogr. Res. 21:341)
- Rogers & McKane 2011 (arXiv 1104.1971)

The catalog's rule is: only discontent agents move; each moves to a uniformly random vacancy regardless of
whether it would be satisfied there; all of a step's moves are applied as a batch. Every quantitative study found
either lets agents move only to satisfying or utility-improving sites, or starts from a perturbed checkerboard.

**Best available target (requires a variant).** Gauvin L, Vannimenus J, Nadal J-P (2009). *Phase diagram of a
Schelling segregation model.* Eur. Phys. J. B 70:293-304. doi:10.1140/epjb/e2009-00234-0. arXiv 0903.4694.

Model (pp. 2-3):

- L = 50 with free boundaries.
- Moore neighbourhood.
- An agent is satisfied if N_d <= T x (number of occupied neighbours).
- A randomly chosen agent, satisfied or not, moves to a vacancy where it would be satisfied. Otherwise it stays.
  Quote: "satisfied agents can also move, not just discontent ones".
- One MC step is one attempt per occupied site.
- Order parameter s = 2 sum_c n_c^2 / (L^2(1-rho))^2, with clusters defined by 4-connectivity (Eq. 2, p. 4).
- Averages over 30000 steps after equilibrium.

**Target A: segregated-to-mixed transition, Fig. 3 on p. 5 (reading at 400 dpi, +/-0.03).**

| rho (vacancies) | <s> at T = 5/7 | <s> at T = 3/4 |
|---|---|---|
| 2% | 1.00 | 0.03-0.09 |
| 6% | 1.00 | 0.03-0.09 |
| 12% | 0.96 | 0.03-0.09 |
| 18% | 0.90 | 0.03-0.09 |
| 26% | 0.75 | 0.03-0.09 |

At T = 3/4 all densities read 0.03-0.09, with reading uncertainty +/-0.02. The collapse therefore lies between the
adjacent meaningful tolerances 5/7 and 3/4. The conclusion (p. 12) says segregation extends "as high as 3/4"; the
figure shows 3/4 already mixed, and the figure is used.

**Target B: frozen-to-segregated line, Table 1 in Appendix A.4.6 on p. 20.**

- These are exact discrete values: the tolerance above which more than 95% of 100 runs segregate, with no figure
  reading.
- Unambiguous entries:
  - 1/2 at rho = 2% and 4%
  - 3/8 at rho = 10%, 12% and 14%
  - 1/3 at rho = 20%
  - 1/5 at rho = 30-40%
  - 1/8 at rho = 46%
- Other densities list a metastable range.

**Mapping.** Our threshold theta (minimum same-type fraction) corresponds to T = 1 - theta. The catalog's
theta = 3/8 is T = 5/8, and density 0.9 is rho = 10%. That point lies inside Gauvin's segregated phase.

### 2. Model match

| Item | Gauvin 2009 | Ours (`run_schelling`) | Fix |
|---|---|---|---|
| Who moves | Any randomly chosen agent | Only unhappy agents | **Code** |
| Destination | A uniformly chosen satisfying vacancy, else stay | Any random vacancy | **Code** |
| Update | Random sequential, one agent at a time | Batch: unhappy agents identified on the step-start grid, then all moved in one pass | **Code** |
| Satisfaction | N_d <= T N_occ (inclusive); no occupied neighbours counts as satisfied | Unhappy if same/occ < theta, which is satisfied iff N_d <= (1-theta) N_occ; isolated agents content | Parameter: T = 1 - theta |
| Neighbourhood, boundary | Moore, free | Moore; padding with 0 is excluded from counts, so free | None |
| Composition | Equal types | n_occ//2 and the remainder | None |
| Observable | s, 4-connected clusters | Moran's I and the detector's segregation index | **Analysis code:** `scipy.ndimage.label` with 4-connectivity, no wrap |
| Ambiguities, pre-registered | Not specified | | (i) The mover's current site counts as vacated when testing a destination. (ii) Choose uniformly among satisfying vacancies. Sensitivity: one random trial vacancy, accepted if satisfying. |

### 3. Protocol, if the variant is adopted

- **Code.** Add `run_schelling_gauvin(grid_size=50, vacancy=rho, T, n_steps, seed)` as a separate function. Do not
  change `run_schelling`. Keep incremental neighbour-count arrays so a satisfying-vacancy query costs O(V).
  numba is recommended.
- **Target A points.**
  - rho in {2%, 12%, 26%}.
  - T in {2/3, 5/7, 3/4, 4/5}. Only 5/7 and 3/4 are gating; 2/3 and 4/5 are reported.
  - Seeds 0-9.
- **Target B points.**
  - (rho, T) in {(2%, 3/7), (2%, 1/2), (12%, 1/3), (12%, 3/8), (30%, 1/6), (30%, 1/5)}.
  - Seeds 0-19.
- **Run length.** 2000 MC steps of equilibration, then 3000 MC steps measuring s every 10 steps. Sensitivity: 30000
  measurement steps at (12%, 5/7), (12%, 3/4) and (2%, 1/2).
- **Non-gating characterization.** Run the catalog rule at the same (rho, theta = 1 - T) points to document how it
  differs. For example, at high theta the catalog rule churns rather than freezing.

### 4. Pass rules (fixed)

**Target A passes** if both hold:

- The seed-mean <s>(5/7) is >= 0.90 for rho = 2% and 12%, and >= 0.60 for rho = 26%.
- The seed-mean <s>(3/4) is <= 0.15 for all three densities.

**Target B passes** if, at each density:

- at least 18 of 20 seeds end frozen (s < 0.5) at the lower T, and
- at least 18 of 20 end segregated (s >= 0.5) at the upper T.

If the true rate is 95%, P(>= 18/20) = 0.92.

**Overall P1 passes** if A and B both pass.

Why these thresholds: the gaps between adjacent discrete T values make the outcomes near-binary, so the thresholds
sit well inside the published values minus the reading uncertainty.

### 5. Compute and risks

**Compute.** Plain numpy costs about 30 µs per attempt and about 2400 attempts per MC step, so about 5 min per
5000-step run.

| Run set | Runs | Plain numpy | numba |
|---|---|---|---|
| Target A | 120 | | |
| Target B | 120 | | |
| **Total** | 240 | **about 20 CPU-h** | **under 1 h** |

The existing `run_schelling` characterization adds under 0.5 h.

**Risks.**

1. **This validates a variant, not the catalog model.** A pass would show that our lattice, neighbourhood and
   satisfaction code reproduce a published phase boundary. It would not support quantitative claims about the
   catalog rule.
2. **Unstated rule details** in Gauvin, pre-registered above.
3. **Metastability near T_f.** The paper itself reports this; only unambiguous Table 1 entries are used.
4. **Finite size** at L = 50 matches the paper's own L, which is intended.

**Recommendation.** Adopt the variant only if the paper needs a quantitative P1 anchor. Otherwise reclassify P1's
reproduction as qualitative and remove the "Quantitative (Moran's I)" label.

---

## Issues found in passing

These are outside this spec but affect the catalog.

1. **P27.** The canonical b = 1.8 and the regime table's b = 1.5 and 2.0 are all break points. Results there depend
   on scan-order tie-breaking. The panel's "stochastic cooperation collapse" at 50x50 (Sprint 39) was observed at
   b = 1.8.
2. **P23.** `np.argmax` tie-breaking toward strategy 0 is not the rule of any source.
3. **P18.** FK 1996 Eq. (23) should not be cited as matching simulations at t <= 10^4. Its truncation error is 7-11%.
4. **P1.** The catalog rule has no published quantitative counterpart.
