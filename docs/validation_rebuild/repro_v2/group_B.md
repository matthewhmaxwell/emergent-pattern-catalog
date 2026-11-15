# Group B literature-reproduction specs: P5, P6, P10, P19

Frozen 2026-11-15, before any simulation is run. No model or simulation code was executed to prepare this
document. Two kinds of computation were done, both on published material only:

- **Figure digitization.** Points were read from the vector content of the published PDFs (glyph and polyline
  coordinates), not by eye, except D'Orsogna Fig. 2(d), which is a scan.
- **Root-finding on published analytic equations.** D'Orsogna 2006 Eq. (4); Abrams & Strogatz 2004 Eqs. (12)-(13).

All scripts and outputs are in `group_B_support/`, so a reviewer can regenerate every frozen number from the
public PDFs. The PDFs themselves are not stored in the repo. Each section gives one verbatim quote. Quotes are kept
under 15 words because of the copyright policy this work runs under. Exact locations (equation, figure, page) are
given instead of longer excerpts. Equations, parameter values and numbers are transcribed exactly.

Standard applied: a reproduction counts only if a specific published number, or an analytic result stated in the
source, is compared with model output under a pass rule fixed here.

| Pattern | Primary frozen target | Verdict |
|---|---|---|
| P5 flocking (Vicsek 1995) | v_a(eta) at N=400, L=10 (rho=4), 9 noise values, digitized from Fig. 2a | Target found, feasible. Needs one small code option (Vicsek's position-update order). |
| P6 milling (D'Orsogna 2006) | Rotating-ring radius R(N) from Eq. (4) (analytic) at N=100-400. Secondary: energy exponent N^2.0 in the vortex regime (Fig. 5) | Target found, feasible. Parameter settings only. |
| P10 chimera (Abrams & Strogatz 2004) | Stable-chimera solution of Eqs. (12)-(13) at A=0.995, beta=0.18 (r, Delta, f_drift, R(x) range), plus saddle-node beta_SN=0.2635 (Fig. 2) | Target found, feasible, strongest of the four. Needs a dt reset (parameter) and the paper's IC (attribute override). |
| P19 leadership (Couzin 2005) | Accuracy vs p, Fig. **1a** (not 2a), digitized exactly | Weak. Our model is not Couzin's model, and no parameter setting closes the gap. Drop the quantitative claim, or implement the Couzin 2005 zonal model as a new model. |

---

## P5: Translational alignment / flocking (Vicsek et al. 1995)

### 1. Published target

**Source.** Vicsek T, Czirok A, Ben-Jacob E, Cohen I, Shochet O (1995). *Novel type of phase transition in a
system of self-driven particles.* Phys. Rev. Lett. 75(6):1226-1229. doi:10.1103/PhysRevLett.75.1226.
Read from the authors' arXiv posting cond-mat/0611743 (2006): text pp. 3-7, Fig. 2 on PDF page 20. The typeset
PRL figure was not compared (paywalled). The arXiv figure is a vector plot of the same data series.

**Where.** Fig. 2(a): v_a versus eta. Settings from the text and caption:

- Five sizes: N=40/L=3.1, N=100/L=5, N=400/L=10, N=4000/L=31.6, N=10000/L=50.
- Fixed density rho = N/L^2 = 4.0 (4.16 for N=40). The caption gives the sizes but not rho.
- v=0.03, r=1, Delta t=1.
- Noise uniform on [-eta/2, eta/2].
- Random initial positions and headings.

**Published uncertainty.** Text p. 6: "The uncertainity of the data points is within the range of the symbols" [sic].
The exception is N=4000 and N=10000 near the transition, where the paper reports about 5% statistical error from
5 runs. Symbol heights, measured from the glyph bitmaps, are 0.027 v_a units (x), 0.030 (square) and 0.038 (+).
That puts the implied published uncertainty at about +/-0.013 to +/-0.019.

**Digitization (exact, not by eye).** The PDF draws each marker as a positioned Type-3 bitmap glyph:

| Glyph | Series |
|---|---|
| square | N=40 |
| + | N=100 |
| x (cid:2) | N=400 |
| triangle | N=4000 |
| diamond | N=10000 |

Method:

- Axis ticks were located on a 300-dpi raster: eta=0 at 164.40 pt, eta=5 at 398.16 pt; v_a=1 at 96.96 pt, v_a=0 at 312.72 pt.
- Glyph bitmap centres come from the font CharProcs. The model was checked against label and legend glyphs (residual at most 0.16 pt, about 0.0007 v_a).
- Each series was anchored at v_a(eta=0.1) = 0.9995. The applied shifts were 0.004 or less.
- Points fall on the published eta = 0.1k grid to within +/-0.01.

Digitization uncertainty is +/-0.005 in v_a. This includes the glyph-centring ambiguity and the slope times eta error near the transition. Script and full output:
`group_B_support/p5_digitize_vicsek1995_fig2a.py`, `p5_vicsek1995_fig2a_digitized.json`.

**Frozen published values (Fig. 2a):**

| eta | N=400, L=10 (primary) | N=100, L=5 | N=40, L=3.1 | N=4000, L=31.6 (optional) |
|---|---|---|---|---|
| 0.5 | 0.981 | 0.983 | 0.985 | - |
| 1.0 | 0.920 | 0.932 | 0.936 | 0.884 (eta=1.1) |
| 1.5 | 0.830 | 0.848 | 0.863 | 0.787 |
| 2.0 | 0.709 | 0.741 | 0.767 | 0.645 |
| 2.5 | 0.566 | 0.621 | 0.654 | 0.463 |
| 3.0 | 0.380 | 0.461 | 0.515 | 0.142 |
| 3.5 | 0.133 | 0.278 | 0.380 | 0.045 |
| 4.0 | 0.081 | 0.157 | 0.260 | 0.028 (eta=3.9) |
| 4.5 | 0.062 | 0.121 | 0.199 | 0.021 |

**The exponent, considered and not used for pass/fail.**

- Vicsek report beta = 0.45 +/- 0.07 and delta = 0.35 +/- 0.06 (text p. 7; PRL p. 1228).
- They obtained beta by choosing eta_c(L) so that the log-log plot (Fig. 3a) was "straightest", which leaves a free parameter.
- The +/-0.07 band is wide enough to have little discriminating power.
- There is an internal inconsistency. Fig. 3a and the eta_c(infinity) = 2.9 +/- 0.05 statement are labelled rho = 0.4. Fig. 2a is at rho = 4, and its N=10000 curve already falls to about 0.10 at eta = 3.0. rho = 4 is the likelier reading, so eta_c(infinity) is not used either.

**Later values considered.**

- Czirok, Stanley & Vicsek 1997, J. Phys. A 30:1375: beta = 0.42 +/- 0.03 at rho = 2. That study used cell-list neighbourhoods and v0 = 0.1, a different interaction rule, so it is not a target for our model.
- Baglietto & Albano 2008, PRE 78:021125, as tabulated in Baglietto, Albano & Candia 2012 (arXiv:1206.2901, Table 1): beta = 0.45(3) and beta/nu = 0.275(5) from finite-size scaling.
  - That work uses a different noise convention (xi in [-pi, pi], eta in [0, 1]) and N up to 32768.
  - The order of the transition is contested (Gregoire & Chate 2004, PRL 92:025702; Chate et al. 2008, PRE 77:046113).
  - It is too expensive, and too contested, to serve as the headline test. It is listed only as a possible extension.

### 2. Model match (`epc/models/vicsek.py` vs Vicsek 1995)

| Item | Vicsek 1995 | Ours | Fix |
|---|---|---|---|
| Heading rule | arctg of mean sin / mean cos over neighbours within r, including self | arctan2 of the same, cKDTree periodic, includes self | Same |
| Noise | Uniform [-eta/2, eta/2] | Same | Same |
| **Position update** | Eq. (1): x(t+1) = x(t) + v(t) Delta t, using the **old** velocity | Uses the **new** heading theta(t+1) | **Code change.** Add an option such as `position_update="old_velocity"`: move with the current headings, then update the headings. About 3 lines; the default stays unchanged. Baglietto & Albano call the original rule the "backward update" and warn that update schemes can create artefacts. |
| Speed, radius, dt | v=0.03, r=1, Delta t=1 | Settable | Parameter |
| Box / density | N/L as captioned (rho=4) | Canonical panel run uses N=300, L=7 (rho=6.12, the Fig. 1 snapshots, which carry no numbers) | Parameter: run N/L exactly as in the Fig. 2 caption |
| Initial condition | Random positions and directions | `init_mode="random"` | Same |
| Order parameter | v_a = abs((1/Nv) sum of v_i) | `polarization()`, identical | Same. The estimator is the time mean of the magnitude, which matches the finite-N tails in Fig. 2a (for example N=40 at eta=4.9 gives 0.167). |
| Averaging and run length | Not stated | Our choice | Pre-registered below |

### 3. Protocol

- **Sizes.** Primary N=400, L=10. Secondary N=100, L=5 and N=40, L=3.1. Optional N=4000, L=31.6 at eta in {2.0, 2.5, 3.0}.
- **Noise grid.** eta in {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5}.
- **Update rule.** Primary runs use Vicsek's position-update order (`old_velocity`). The same grid at N=400 is also run with the current order and reported as a sensitivity result, not pass/fail.
- **Seeds.** 10 per (N, eta), seeds 1000-1009; 5 for optional N=4000.
- **Run length.**
  - Burn-in 2x10^4 steps. At N=400 this is 60 crossing times (L/v = 333 steps).
  - Measurement 8x10^4 steps, sampling phi every 10 steps.
  - Convergence check: if the first-half and second-half measurement means differ by more than 0.02 at any (N, eta, seed), extend that cell's measurement once to 2x10^5 steps and use the extended value.
- **Estimator.** Per seed, take the time mean of phi(t). Then v_sim(N, eta) is the mean over seeds. SE is the standard deviation of the seed means divided by sqrt(10). Report 95% t-intervals.
- **Finite-size handling.** Compare like with like. The published points are finite-N values at the same N and L, so no extrapolation is needed. The secondary finite-size check: the simulation must reproduce the published ordering v(40) > v(100) > v(400) at every eta >= 2.5.

### 4. Tolerance (fixed)

- **PASS (primary)** if |v_sim - v_pub| <= 0.04 at all 9 eta values for N=400.
- **Justification.** The components are the half-symbol published error (0.014), digitization (0.005), our SE (<= 0.01), and unknown averaging and run-length conventions near eta = 3-3.5, where Vicsek note slow convergence. Two times their root sum of squares is 0.036, rounded up to 0.04.
- **Discriminating power.** The canonical rho = 6.12 setting would miss at eta = 3.0 by +0.23. A 0.1 shift in eta_c moves v_a by about 0.05 near eta = 3.2.
- **Secondary.** The same rule applies at N=100 and N=40, plus the ordering check above.
- **Optional N=4000.** Tolerance +/-0.06 (the paper's own 5% error plus slow convergence). Reported, not pass/fail.
- **beta.** Computed descriptively by Vicsek's straightest-line procedure. It is not pass/fail.

### 5. Compute and non-model risks

- **Cost.** REPLICATION_NOTES gives 1.2 ms/step at N=300, so roughly 1.5 ms/step at N=400 and 0.5 ms/step at N <= 100.
  - Each cell is about 1.0x10^5 steps.
  - N=400 alone is 90 runs x ~150 s = about 4 h CPU. N=100 and N=40 are about 2.5 h CPU together.
  - Optional N=4000 is about 12 ms/step, 20 min/run, about 8 h CPU.
  - All runs are independent; run them on VPS cores with `OMP_NUM_THREADS=1`.
- **Memory.** `VicsekModel.run()` stores a full state dict every step: about 1.6 GB per 10^5-step run at N=400. Use a `step()` loop with a streaming phi estimator, or `record_interval >= 10`.
- **Boundary crash.** `positions %= L` can return exactly L in floating point, and `cKDTree(boxsize=L)` then raises ValueError. Guard by mapping values >= L to 0.
- **Figure provenance.** The arXiv posting is the authors' manuscript. Differences from the typeset PRL figure are possible but unverified.
- **Record correction.** `reproductions.json` cites "Fig. 3" for phi vs eta. The phi-vs-eta curve is Fig. 2a; Fig. 3 is the log-log exponent plot.

---

## P6: Milling / vortex formation (D'Orsogna et al. 2006)

### 1. Published target

**Source.** D'Orsogna MR, Chuang YL, Bertozzi AL, Chayes LS (2006). *Self-propelled particles with soft-core
interactions: patterns, stability, and collapse.* Phys. Rev. Lett. 96:104302. doi:10.1103/PhysRevLett.96.104302.
Read from the author-hosted scan of the published PRL (math.ucla.edu/~bertozzi/papers/PRL06.pdf), 4 pages.

**Primary: analytic ring radius, Eq. (4), p. 104302-2**

Transcribed exactly:

    m*alpha / (2 R beta) = sum_{n=1}^{N/2} sin(pi n / N) *
        [ (C_a/l_a) exp(-2R sin(pi n/N) / l_a) - (C_r/l_r) exp(-2R sin(pi n/N) / l_r) ]

- **Regime.** Region II of Fig. 1 (l = C, where C = C_r/C_a and l = l_r/l_a). The rings are shown in Fig. 2(c).
- **Parameters**, from the Fig. 2 caption: alpha = 1, beta = 0.5, C_a = l_a = 1, C_r = 0.5, l_r = 0.5, unit mass.
- **What the paper says.** The text calls Eq. (4) an approximate implicit expression from balancing the centrifugal and centripetal forces for equidistant particles. Quote, p. 104302-2: "Estimates of R as given by Eq. (4) match extremely well those obtained numerically". Fig. 2(d) plots numeric and analytic R against N. The caption states R ~ N^-0.5.

**Frozen values.** Root-finding by `group_B_support/p6_solve_dorsogna2006_eq4.py`; one root at each N:

| N | R from Eq. (4) as printed | R from the exact force balance (sum n=1..N-1) | Fig. 2(d) numeric (digitized from scan, +/-0.0015) |
|---|---|---|---|
| 100 | 0.17426 | 0.17630 | 0.1767 |
| 200 | 0.11503 | 0.11567 | 0.1159 |
| 300 | 0.09133 | 0.09166 | 0.0917 |
| 400 | 0.07783 | 0.07804 | 0.0778 |

- **Why the two analytic columns differ.** The printed sum to N/2 with the 2 on the left-hand side double-counts the diametrically opposite partner. The resulting gap is 1.2% at N=100 and 0.3% at N=400.
- **Cross-check of the reading.** Our Eq. (4) solutions match the paper's plotted analytic dots to within about 1-2% across N = 100-500, which is the scan distortion level.
- **Range of the numeric curve.** Fig. 2(d) numeric points cover N = 100-460 in steps of 10. The numeric curve ends with a vertical drop near N of about 460, which the text does not explain. It is not used.

**Secondary: energy scaling in the vortex (mill) regime, Fig. 5 caption, p. 104302-4**

- Catastrophic regime VII uses the Fig. 3 Morse parameters: C_a = 0.5, C_r = 1, l_a = 2, l_r = 0.5, beta = 0.5.
- The caption gives the scaling |U_total| ~ N^2.0 for alpha = 1 and for alpha = 0.
- The H-stable control (C_a = 0.5, C_r = 1, l_a = 2, l_r = 1.5) gives |U_total| ~ N^1.0 for alpha = 0 and 0.5.
- U_total = sum over i of U(x_i), with U as in Eq. (3) (pair sum, double counted).
- This regime is the same class as our canonical P6 positive: C l^2 = 0.056, region VII.

**Not used.** Fig. 4 vortex area and inner/outer radii. The text never defines those estimators, the figure is a
scan, and the reading uncertainty would dominate.

### 2. Model match (`epc/models/dorsogna_spp.py` vs PRL Eqs. 1-3)

| Item | PRL | Ours | Fix |
|---|---|---|---|
| Equations of motion | dx/dt = v; m dv/dt = (alpha - beta abs(v)^2) v - grad U | Same, with m = 1 implicit | Same. The force sign was checked: a = -sum_j U'(r) r_hat_ij with U' = (C_a/l_a)e^{-r/l_a} - (C_r/l_r)e^{-r/l_r}. |
| Potential | sum over j != i of [C_r e^{-r/l_r} - C_a e^{-r/l_a}] | Same, diagonal zeroed | Same |
| Integrator | 4th-order Adams-Bashforth | RK4 (fixed dt) | Same order; does not change attractors. dt is pre-registered with a dt/2 check. |
| Space | Unbounded | Open space | Same |
| Initial condition | "Localized particles of unitary mass and random velocity" (extent and distribution not stated) | `random`: uniform in [-R0, R0]^2, random direction at speed sqrt(alpha/beta) | Parameter: `init_radius = 1.0` (pre-registered); sensitivity at 0.5 and 2.0 |
| Parameters for the target | Region II ring set; Fig. 3/5 set | Canonical panel uses the Carrillo et al. 2009 set (C_a = 0.5, C_r = 1, l_a = 3, l_r = 0.5), not checked here | Parameter settings only. No code change. |

### 3. Protocol

**Primary, T6.1: ring radius**

- **Settings.** C_a = 1, l_a = 1, C_r = 0.5, l_r = 0.5, alpha = 1, beta = 0.5; N in {100, 200, 300, 400}.
- **Seeds.** 10 per N (seeds 2000-2009). The dynamics are deterministic; seeds vary only the initial condition.
- **Integration.** RK4 with dt = 0.01, T = 200, measuring over t in [150, 200], sampled every 0.1.
- **dt check.** Repeat one seed at N=100 and one at N=400 with dt = 0.005. R must change by less than 0.1%, or every N is rerun at dt = 0.005.
- **Estimator.** Per sample, R(t) = mean over i of abs(x_i - x_cm). R_seed is the time mean. R_sim(N) is the median over valid seeds; also report the inter-seed range.
- **Validity gates.** All are applied per seed, on time means over the window.
  - Radial CV = std over i of abs(x_i - x_cm), divided by R, must be below 0.05.
  - Normalized angular momentum abs(L) must exceed 0.95. Use `angular_momentum` with velocities in open space, without periodic unwrapping.
  - Mean speed must be within 2% of sqrt(2).
  - The largest abs(x_i - x_cm) must be below 2R. This catches escaped particles.
  - Convergence: the window's first-half and second-half R must differ by less than 0.5%, otherwise extend once to T = 400.
- **Handling failed seeds.** Seeds that fail a gate are reported along with the state they reached. They are not used for R_sim.

**Secondary, T6.3: energy exponent**

- **Settings.** Catastrophic set (C_a = 0.5, C_r = 1, l_a = 2, l_r = 0.5, beta = 0.5, alpha = 1) and H-stable set (l_r = 1.5, alpha = 0.5); N in {50, 100, 200, 400}; 5 seeds per point.
- **Run.** dt = 0.01, T = 100, `init_radius = 1.0`.
- **Estimator.** |U_total| is evaluated at the final state. The exponent is the least-squares slope of log of median |U_total| against log N.

### 4. Tolerance (fixed)

**T6.1 PASS** requires both conditions below.

- At least 8 of 10 seeds at every N pass all gates. That is, the ring state is reached and is a co-rotating ring.
- abs(R_sim/R_Eq4 - 1) <= 0.03 at every N in {100, 200, 300, 400}.

Justification for the 3% band:

- Eq. (4) is itself approximate, by up to 1.2% at N = 100.
- The paper's own numeric-vs-analytic agreement in Fig. 2(d) is within about 1-2%.
- Inter-seed spread should be tiny because the ring is a unique attractor.

**T6.2 (implementation check).** abs(R_sim/R_exact - 1) <= 0.01 at every N. If T6.2 fails while T6.1 passes, record
it as a model or integration discrepancy.

**T6.3 PASS.** The fitted exponent must lie in [1.9, 2.1] for the catastrophic alpha = 1 set and in [0.9, 1.1] for
the H-stable alpha = 0.5 set. The published "2.0" and "1.0" carry two significant figures, and fitting over N = 50-400
instead of the paper's 50 to about 1500 allows about +/-0.05 of curvature.

### 5. Compute and non-model risks

- **Cost.** REPLICATION_NOTES gives 8.5 ms/step at N = 100 with dt = 0.01 RK4; cost scales as N^2.
  - T = 200 is 2x10^4 steps.
  - Per run: about 3 min (N=100), 12 min (N=200), 26 min (N=300), 45 min (N=400).
  - T6.1 total for 10 seeds per N is about 14 h CPU. Runs are parallel across VPS cores.
  - T6.3 adds about 6 h CPU.
- **Memory.** `run()` stores every step by default: about 0.4 GB at N = 400 over 2x10^4 steps. Use `record_interval = 10`.
- **State selection.** Random initial conditions could reach clumped rings or translating clumps instead of a ring. The gates make that visible, and a shortfall counts as a fail rather than a silent exclusion.
- **Record corrections.** The canonical P6 parameters are attributed to Carrillo et al. 2009, Fig. 3.1. That attribution was not verified here; the targets above use the PRL's own parameter sets. In `reproductions.json`, "phi = 0.008" vs "0.018" is an existing conflict, unaffected by this spec.

---

## P10: Chimera states (Abrams & Strogatz 2004)

### 1. Published target

**Source.** Abrams DM, Strogatz SH (2004). *Chimera states for coupled oscillators.* Phys. Rev. Lett. 93:174102.
doi:10.1103/PhysRevLett.93.174102. Read from arXiv:nlin/0407045 (same text and figures; vector figures).

**Model, Eqs. (1)-(2).**

- dphi/dt = omega - integral over x' in [-pi, pi] of G(x - x') sin[phi(x) - phi(x') + alpha]
- G(x) = (1 + A cos x)/(2 pi)
- beta = pi/2 - alpha
- The Fig. 1 caption gives:
  - A = 0.995, beta = 0.18, N = 256.
  - Fixed-step Runge-Kutta with dt = 0.025 for 200,000 iterations, so T = 5000.
  - Initial condition phi(x) = a exp(-b x^2) r, with r uniform on [-1/2, 1/2], a = 6, b = 30, and x in [-pi, pi].

**Analytic results stated in the source.**

- Self-consistency Eqs. (12)-(13) for c (real), a = a_r + i a_i, and Delta = omega - Omega, with R(x)^2 = abs(c + a cos x)^2 (Eq. 10).
- Fig. 2 is the existence region in (beta, A), computed from Eqs. (12)-(13) by continuation. Quote (Fig. 2 caption): "A saddle-node bifurcation of a stable and unstable chimera occurs at the boundary."
- Perturbation theory (p. 3): (beta_1)_max of about 0.2205, giving a boundary slope 1/(beta_1)_max of about 4.535. The minimum drift fraction is about 44%. Chimera birth occurs at delta = 16/(9 pi^2).

**What we froze.**

- **(a) Stable branch of Eqs. (12)-(13).** Solved by `group_B_support/p10_solve_as2004_selfconsistency_branches.py` (200,000-point quadrature, principal square-root branch).
  - Two chimera solutions exist at A = 0.995, beta = 0.18. The stable one is identified by:
    - the paper's statement that the stable branch has the larger drift fraction (Fig. 3);
    - the R range in Fig. 1(b).
  - Fig. 1(b), digitized from vector rectangles, gives R_max of about 0.809 and R_min of about 0.474. The solution gives 0.8068 and 0.4732, a match within 0.003.
- **(b) Saddle-node location.** Our fold computation gives beta_SN(0.995) = 0.2635. The solid curve in Fig. 2, digitized from vector paths (`p10_digitize_as2004_fig2.py`), gives 0.2625.
  - Calibration check: the digitized dashed line has slope 4.535, which is the value printed in the text.
  - Folds at A = 0.95 and 0.90 are 0.2466 and 0.2293, against 0.2461 and 0.2291 digitized.

| A = 0.995 | c = global r | Delta (paper time units) | f_drift | R_min | R_max |
|---|---|---|---|---|---|
| beta = 0.18 (paper's Fig. 1) | 0.6287 | 0.7320 | 0.7008 | 0.4732 | 0.8068 |
| beta = 0.10 | 0.5933 | 0.7259 | 0.7884 | 0.4317 | 0.7623 |
| beta = 0.05 (our canonical positive) | 0.5760 | 0.7160 | 0.8526 | 0.4214 | 0.7326 |
| Fold beta_SN | 0.2635 (fold); digitized Fig. 2: 0.2625 +/- 0.001 | | | | |

Notes:

- c is the stationary global order parameter abs(<e^{i theta}>) in the N -> infinity limit.
- The locked cluster rotates at Omega = omega - Delta = -Delta, with omega = 0.
- The existing unregistered runs already gave r_global of 0.577-0.582 at beta = 0.05 and 0.629-0.633 at beta = 0.18 (N = 128). Those agree with c, but used the coarse dt described below.

### 2. Model match (`epc/models/kuramoto_nonlocal.py` vs AS2004)

| Item | AS2004 | Ours | Fix |
|---|---|---|---|
| Kernel, phase lag, omega | Cosine, alpha = pi/2 - beta, omega irrelevant | Same (omega = 0) | Same |
| **Quadrature weight / time units** | Riemann sum (2 pi/N) sum_j | Sum_j with no 2 pi/N. The right-hand side is N/(2 pi) times faster. | **Parameter.** Set dt_ours = 0.025 x 2 pi/N. At N = 256 that is 6.1359e-4, which equals the paper's dt = 0.025. Convert t_paper = t_ours x N/(2 pi) and Omega_paper = Omega_ours x 2 pi/N. **Our default dt = 0.025 equals 0.51 (N=128) to 1.02 (N=256) paper time units per RK4 step, 20-40x coarser than the paper's.** Existing results were all produced at that coarse step. |
| **Initial condition** | 6 exp(-30 x^2) times U[-1/2, 1/2], x in [-pi, pi] | `init_mode="paper"` uses 6 exp(-30 (x/pi - 1)^2) times N(0,1)/1.5 | **Not the paper's IC.** The envelope is pi times wider in x (b_eff = 30/pi^2 = 3.04) and the noise is Gaussian, not uniform. Fix without a code change: overwrite `m.theta` after construction with 6 exp(-30 (x_i - pi)^2) r_i, r_i ~ U[-1/2, 1/2] from `default_rng(seed)`. Alternatively add an `init_mode="as2004"` (about 4 lines). |
| N, A, beta, integrator | 256, 0.995, 0.18, fixed-step RK | Settable; RK4 | Parameter |
| Run length | 200,000 steps (T_paper = 5000) | `run(n_frames, ...)` | Parameter. T_ours = 5000 x 2 pi/N = 122.72 at N = 256. |
| Frequency measurement | Not given | `run()` stores theta mod 2 pi per frame; frames are aliased at record_dt = 1 | Use the raw unwrapped `m.theta` after successive `run()` calls (`run()` continues from `self.theta`). Set record_dt to an exact multiple of dt. |

### 3. Protocol

**Stage A: profile at the paper's parameters (T10.1)**

- **Run.** N = 256, A = 0.995, beta = 0.18, seeds 0-9, paper IC.
  - dt_paper = 0.025, i.e. dt_ours = 6.1359e-4.
  - 200,000 steps (T_paper = 5000).
  - Measurement window: T_paper in [2500, 5000].
- **dt check.** Two seeds are repeated at dt_paper = 0.0125. All Stage A estimators must change by less than half their tolerance.
- **Estimators** (paper units):
  - **r_bar.** The time mean of abs((1/N) sum_j e^{i theta_j}), sampled every paper time 1.0.
  - **Frequencies.** Omega_i = [theta_i(5000) - theta_i(2500)]/2500 from unwrapped phases.
  - **Locked set.** Bin the Omega_i at width 0.001. The locked set is the largest bin plus neighbours within 0.002.
  - **Delta_sim and drift fraction.** Delta_sim = -median(Omega_locked); f_drift = 1 - n_locked/N.
  - **R(x).** Every paper time 10, compute R(x) = abs(Z0 + A (Zc cos x + Zs sin x)), the exact first-harmonic form of Eq. (3), where Z0 = mean of e^{i theta_j}, Zc = mean of cos(x_j) e^{i theta_j}, and Zs = mean of sin(x_j) e^{i theta_j}. Evaluate on a 1024-point x grid. R_min and R_max are the snapshot minimum and maximum, averaged over snapshots. This handles the spatial drift of the chimera, logged as ADR 53.
- **Classification.** A seed is a chimera if, over the window, r_bar < 0.95 and 0.3 < f_drift < 0.95. Basin fraction is reported descriptively.
- **Fallback.** If fewer than 5 of 10 paper-IC seeds end in a chimera, add chimera realizations created as follows: the asymmetric-Gaussian IC at beta = 0.05, then an adiabatic ramp to beta = 0.18 in steps of 0.01, holding T_paper = 500 per step. Continue until 10 chimera realizations exist, and label which is which.

**Stage B: existence boundary (T10.2)**

- **Start.** Each Stage A chimera state at beta = 0.18.
- **Ramp.** Step beta to 0.20, 0.22 and 0.24, holding T_paper = 500 at each. Then step 0.2400 -> 0.2800 by 0.0025, holding T_paper = 1000 at each.
- **Collapse rule.** Collapse occurs when r > 0.95 for the last 200 paper time of a hold and f_drift < 0.02.
- **Estimator.** beta_collapse(seed) is the first beta at which collapse occurs; the summary is the median over seeds.
- **Secondary.** Repeat at N = 128 (10 seeds) and N = 512 (5 seeds, optional) to show the finite-N trend.
  - Wolfrum & Omel'chenko 2011 (PRE 84:015201(R)) show that finite-N chimeras are chaotic transients, with lifetimes growing exponentially in N. Their model used a step kernel, so lifetimes do not transfer.
  - Expect beta_collapse to approach 0.2635 from below as N grows.

**Stage C: canonical cross-check (T10.3).** Repeat Stage A at beta = 0.05, N = 256, against the beta = 0.05 row of the
table, with the same tolerances. This ties the analytic targets to the configuration the catalog uses.

### 4. Tolerance (fixed)

**T10.1 PASS.** For the median over chimera seeds, all five of the following hold.

- abs(r_bar - 0.6287) <= 0.02
- abs(Delta_sim - 0.7320) <= 0.010
- abs(f_drift - 0.7008) <= 0.03
- abs(R_min - 0.4732) <= 0.03
- abs(R_max - 0.8068) <= 0.03

Justification:

- Finite-N fluctuation of r(t) at N = 128-256 is about 0.025-0.035 (Sprint 18 logs), so the time-mean bias should be well inside 0.02.
- The f_drift quantum is 1/N = 0.004. The +/-0.03 band allows about 8 oscillators of misclassification at the locked-drifting edge.
- The paper's own single N = 256 run (Fig. 1(b)) matches the analytic R range within 0.003.
- Delta is the most precise observable, a long time mean of locked phases. The +/-0.01 band allows O(1/N) quadrature and finite-size shifts.

**T10.2 PASS.** Both conditions below hold at N = 256.

- The median beta_collapse lies in [0.2500, 0.2675].
- No seed still holds a chimera at beta >= 0.2700.

Justification: the continuum fold is 0.2635. The +0.004 margin is the sweep step plus discreteness. The -0.0135
margin allows early escape of finite-N transient chimeras near the fold, where the basin shrinks to zero. A collapse
at about 0.22, as earlier IC-basin scans suggested, would FAIL. That earlier result was a basin-volume effect, not the
existence boundary.

**T10.3 PASS.** Same rules as T10.1, using the beta = 0.05 values.

### 5. Compute and non-model risks

- **Cost.** Per step is O(N^2), with 4 kernel evaluations per RK4 step.
  - The documented N = 128 timings in REPLICATION_NOTES are mutually inconsistent: 0.45 ms/step (Sprint 27: 150 runs in 133 s), 10 ms/step (module docstring), and about 75 ms/step (the "3 s/frame" note).
  - A numpy operation-count estimate gives about 2 ms/step at N = 256. Budget at 10 ms/step to be safe.
  - Stage A is 2x10^5 steps: 7-33 min per seed.
  - Stage B is about 7.4x10^5 steps: 25 min to 2 h per seed.
  - N = 256 total for 10 seeds: 5-26 h CPU. Optional N = 512 is 4x per step.
  - Run a 1000-step timing pilot (not science) before launch.
- **Memory.** Keep n_frames modest. Use raw `m.theta` checkpoints for frequencies.
- **Record corrections**, to fix in `reproductions.json` and the paper's section 4.19:
  - The 2004 PRL reports **no lifetimes** and **no basin-vs-N claim**. The "~5 PDE time units at N=128 to ~50 at N=256" attributed to its Fig. 2 (REPLICATION_NOTES Sprint 27) is not in the paper. Fig. 2 is the existence region.
  - The claim "no numeric boundary map published" is wrong. Fig. 2's solid line is exactly that, a numerical solution of Eqs. (12)-(13).
  - Sprint 18 "Phase 1a, paper IC" did not use the paper's IC (section 2 above). Its "no chimera at the paper IC" conclusion is therefore untested.
  - The time-units note should also state that the coarse dt, not only the observation window, differs by N/(2 pi).

---

## P19: Emergent leadership / minority guidance (Couzin et al. 2005)

### 1. Published target

**Source.** Couzin ID, Krause J, Franks NR, Levin SA (2005). *Effective leadership and decision-making in animal
groups on the move.* Nature 433:513-516. doi:10.1038/nature03236. Read from the publisher PDF (vector figures).

**Where.** Accuracy versus proportion informed p for N = 10, 30, 50, 100, 200 is **Fig. 1a** (p. 514). Fig. 2a is
something else: accuracy and fragmentation versus the weighting omega at p = 0.02 (1 individual), N = 50.
`reproductions.json` and Sprint 69 cite the wrong panel.

**Parameters** (Fig. 1 caption and Methods):

- 400 replicates; omega = 0.5.
- Repulsion zone alpha = 1, interaction range rho = 6, information error gamma = 0.
- Delta t = 0.2 s; turning rate theta = 2 rad/s (0.4 rad/step); speed s = alpha per second.
- Wrapped-Gaussian angular noise sigma = 0.01 rad.
- t_f = 2500 steps.
- Group direction is the vector between the centroid at t_f - 50 and the centroid at t_f.

**Accuracy metric.**

- p. 513 defines accuracy as the "normalized angular deviation of group direction (see Methods) around the preferred direction g".
- Methods cite Batschelet (1981) pp. 34-36, the mean angular deviation s = sqrt(2(1 - r)), normalized so that random directions give 0 and perfect alignment gives 1.
- The formula is not printed. Our pre-registered reading: accuracy = 1 - s/sqrt(2) = 1 - sqrt(1 - r_bar), with r_bar = (1/400) sum_k cos(h_k - g).

**Digitization (exact).** Taken from the series polylines, whose vertices are the data points (`p19_digitize_couzin2005_fig1a.py`).

- Axes: p = 0 at 94.11 pt, p = 1 at 263.26 pt; accuracy 0 at 517.19 pt, 1 at 396.92 pt (tick lines).
- Digitization error is about +/-0.002. The p values sit 0.001 above k/N, a sub-point calibration offset; they are reported as nominal k/N.

Frozen values:

| N = 50: p | 0.02 | 0.04 | 0.06 | 0.08 | 0.10 | 0.20 | 0.30 | 0.50 | 1.00 |
|---|---|---|---|---|---|---|---|---|---|
| accuracy | 0.096 | 0.282 | 0.464 | 0.602 | 0.700 | 0.880 | 0.919 | 0.944 | 0.963 |

| N = 200: p | 0.005 | 0.010 | 0.015 | 0.020 | 0.025 | 0.030 | 0.05 | 0.10 | 0.20 | 0.50 | 1.00 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| accuracy | 0.057 | 0.110 | 0.218 | 0.492 | 0.737 | 0.826 | 0.912 | 0.931 | 0.946 | 0.968 | 0.985 |

Informed count at accuracy 0.5, n_50, by linear interpolation in n:

| N | 10 | 30 | 50 | 100 | 200 |
|---|---|---|---|---|---|
| n_50 | 1.6 | 3.2 | 3.3 | 3.8 | 4.0 |

The Fig. 1a inset is the same data re-plotted against the number of informed individuals; it adds no independent values.

### 2. Model match (`epc/models/informed_minority.py` vs Couzin 2005): fails

| Item | Couzin 2005 | Ours | Fix |
|---|---|---|---|
| Social rules | Repulsion within alpha (priority). Otherwise attraction (sum of unit vectors to neighbours) plus alignment (sum of unit headings) within rho = 6 | Alignment only, metric r = 1 | **New model required** |
| Space | Unbounded; cohesion from attraction | Periodic box L = 10 | **New model** |
| Turning-rate limit | theta Delta t = 0.4 rad/step | None | **New model** |
| Noise | Wrapped Gaussian, sigma = 0.01 rad | Uniform, eta = 0.1 (sd 0.029) | Code (distribution) |
| Informed weighting | d' = (d_hat + omega g)/abs(...), omega = 0.5 | (1 - w) s_hat + w g, normalized, w = 0.3 | **Parameter.** w = omega/(1 + omega) = 1/3 reproduces omega = 0.5 exactly. |
| Accuracy | 1 - sqrt(1 - r_bar) over 400 replicates, from centroid displacement over the last 50 steps | cos of instantaneous mean heading, 10 seeds | Analysis code |
| Initial condition | "Randomized" positions and orientations; area not stated | Uniform in the periodic box | Choice to pre-register |
| Dimension | Not stated in the Letter (Fig. 3 uses a planar angle) | 2D | Assumption to declare |

No parameter setting of our model reproduces Couzin's dynamics. Matching Fig. 1a would test a Couzin 2005
implementation, not the catalog's P19 model. That only supports the catalog claim if the P19 positive moves to that
model. No published number exists for our Vicsek-with-bias configuration (eta = 0.1, w = 0.3, L = 10).

### 3. Protocol (only if a `couzin2005` model is implemented)

**Settings and rules**

- 2D, unbounded space, N in {50, 200}, with p = k/N at the tabulated values.
- 400 replicates per (N, p), matching the paper. Seeds are replicate indices 0-399.
- t_f = 2500 steps, Delta t = 0.2.
- Alignment includes self: Eq. (2) sums j from 1, unlike attraction. This is declared as a reading of the source.

**Pre-registered unstated details**

- Initial positions are uniform in a square of side 2 sqrt(N) alpha, giving about 29 neighbours within rho, a cohesive start. Headings are uniform.
- Sensitivity check at sides sqrt(N) and 4 sqrt(N): 100 replicates at 3 p values. Report it; do not tune to it.
- All replicates count toward accuracy. The fragmentation fraction is reported separately, via connected components at range rho.

**Estimator.** Accuracy = 1 - sqrt(1 - r_bar). Bootstrap the 400 replicates (2000 resamples) for a 95% CI. n_50 is
interpolated in n.

### 4. Tolerance (fixed; applies only to a new implementation)

**PASS** requires all four conditions, for both N = 50 and N = 200.

- Plateau points (published p >= 0.2): max abs(Delta) <= 0.04.
- Steep region: abs(n_50,sim - n_50,pub) <= 1.0 informed individual. Published n_50 is 3.3 at N = 50 and 4.0 at N = 200.
- Accuracy rises monotonically in p.
- Accuracy at p = 1 lies within 0.03 of the published value: 0.963 at N = 50, 0.985 at N = 200.

Justification:

- The published SE at 400 replicates is about 0.002 on the plateau and about 0.02-0.03 in the steep region.
- In the steep region, one informed individual moves accuracy by about 0.25 at N = 200. A pointwise band there would be meaningless, so the test is in units of individuals.
- The plateau band absorbs the declared reading choices. If the p = 1 test fails, the accuracy-normalization reading is wrong; report it as such rather than tuning.

### 5. Compute and non-model risks

- **Cost** (new model). O(N^2) per step: about 0.3-0.5 ms at N = 200, so about 1 s per replicate.
  - 400 replicates x 11 p values (N = 200) is about 1.3 h.
  - N = 50 adds about 15 min; the sensitivity check about 20 min.
- **Risks.** Ambiguity in the accuracy normalization, initial area, dimension and fragmentation handling. With four unstated details, a pass or fail partly tests our readings.
- **Recommendation.** Do not count P19 as a quantitative reproduction in the current paper. Correct the record to Fig. 1a and "qualitative only". Treat a Couzin 2005 re-implementation as a separate, optional work item.

---

## Cross-cutting notes for the plan

- **Frozen-target provenance.** Every number above can be regenerated from `group_B_support/`. The scripts take the public PDF path as an argument; the URLs are in each section.
- **Parallel runs.** Run all simulation work on the VPS with `OMP_NUM_THREADS=1` per process, using one process per (pattern, N, parameter, seed) cell.
- **Results logging.** Record PASS or FAIL per target exactly as defined, including failures. Any protocol deviation, such as an extended run or a fallback chimera, must be logged as pre-registered or as a deviation.
