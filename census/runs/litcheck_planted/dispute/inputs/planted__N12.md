# planted:N12

## Target paper (ground truth; the gate never saw this)
- Title: Logarithmic Quasi-Long-Range Order and Novel Continuous Phase Transition in the Two-Dimensional $XY$ Model with $1/r^4$ Interaction
- arXiv: 2609.39267 (v1 2026-09-30); authors: Xianzhi Pan, Zhijie Fan, Youjin Deng
- Behaviour the authors report as new: The 2D XY model with 1/r^4 couplings has a low-temperature phase whose correlations decay as a power of log r, entered through a continuous transition with power-law divergent correlation length (neither BKT nor Landau type).
- Authors' novelty statement: "First, using Gaussian spin-wave theory, we show that the low-temperature phase at $\sigma=2$ is a novel logarithmic quasi-long-range ordered (log-QLRO) phase, with $C(r)\sim(\ln r)^{-\eta_{\ell}}$ and $\eta_{\ell}=1/(2\pi\kappa)\propto T$ , where $\kappa$ is the dimensionless LR coupling strength."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: an L x L square lattice with wrap-around edges (L up to 2048), each site holding a planar unit spin (an angle). Energy: minus the sum over all pairs of J(r) times the cosine of their angle difference, where every pair interacts with J(r) proportional to 1/r^4 (r = distance), normalised so that the total coupling per site is fixed. Rules: equilibrium Monte Carlo with cluster updates at inverse temperatures beta = 1, 2, 4, 8 and near the ordering point. Measure the squared magnetisation and the correlation length versus L.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Nearest-neighbour coupling 1.2 plus a weak 1/r^4 tail: the tail still changes the low-temperature scaling.

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - log_decay_exponent_beta1 = 0.34  (magnetisation decays as power of log L; exponent at beta=1)
  - log_decay_exponent_beta8 = 0.02  (same exponent at beta=8 (proportional to temperature))
  - corr_length_log_slope_beta1 = 0.48  (slope of squared correlation-length ratio versus log L, beta=1)
  - exponent_slope_product = 1.0  (2*pi times slope times exponent; 1 means spin-wave picture holds)
  - corr_length_divergence_exponent = 2.24  (power-law divergence exponent of correlation length at the transition)

This behaviour appears across 4 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (first pass): verdict KNOWN
- **standard_name**: 2D long-range XY model at sigma=2: log-corrected (ln L)^(-eta(T)) magnetization decay, spin-wave scaling with k^2 ln k dispersion
- **model_family**: statistical physics: 2D classical XY model with long-range algebraic couplings 1/r^(2+sigma), sigma=2 (marginal point), cluster Monte Carlo
- **notes**: Only abstracts/snippets were read. The exact log-exponent formula and the 2*pi*slope*exponent=1 identity were not verified in the full text; they are consistent with a simple spin-wave estimate (k^2 ln k dispersion giving <theta^2> ~ T ln ln L). The 'NN 1.2 + weak 1/r^4 tail' variant was not found in any source and is inferred from the same theory. The corr-length divergence exponent 2.24 was not compared with the published nu.
- **citations**: ["Xiao, Yao, Zhang, Fan, Deng 2025, Chin. Phys. Lett. 42, 070002 - 2D LR XY MC up to L=8192; ordered phase for sigma<2 and logarithmic scaling in the low-T phase at sigma=2", "Yao, Xiao, Zhang, Deng, Fan 2025, Phys. Rev. B (DOI 10.1103/yghs-p5mg) - comprehensive MC; ordered for sigma<=2, crossover to short-range at sigma=2, critical exponents", "Giachetti, Defenu, Ruffo, Trombettoni 2021, Phys. Rev. Lett. 127, 156801 - theory of BKT with long-range couplings r^(-2-sigma)", "Defenu, Trombettoni, Ruffo 2023, Rev. Mod. Phys. - review of long-range quantum systems (seen only as a citation in snippets, not fully verified)", "arXiv 2512.01956 (2025) - 2D long-range Heisenberg model, same sigma*=2 threshold, 1/ln L decay (close relative)"]

### N_blocked_r2 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: 2D long-range XY model at marginal sigma=2 (1/r^4): logarithmic low-temperature scaling from Goldstone spin waves with k^2 ln(1/k) dispersion
- **behaviour_documented**: True
- **notes**: Logarithmic scaling at sigma=2 is documented by the Deng group, who read it as true long-range order (constant + a/ln r). The observed reading is a power of ln L going to zero, with exponent proportional to T and 2*pi*slope*exponent = 1. That fits Gaussian spin-wave theory and Bruno's theorem (no order for alpha >= 2D), but I found it nowhere stated explicitly. L up to 2048 (ln L about 7.6 at most) cannot separate the two forms. The beta=8 exponent is noisy; the ratio of about 17 instead of 8 fits a stiffness reduced at beta=1. The 2.24 divergence exponent is unmatched and probably an effective power-law fit at a boundary case with log corrections. Not an artifact.
- **citations**: ["Xiao, Yao, Zhang, Fan, Deng, Chin. Phys. Lett. 42, 070002 (2025), arXiv:2404.08498", "Yao, Xiao, Zhang, Deng, Fan, Phys. Rev. B 112, 144429 (2025), arXiv:2411.01811", "Bruno, Phys. Rev. Lett. 87, 137203 (2001)", "Giachetti, Defenu, Ruffo, Trombettoni, EPL 133, 57004 (2021); PRL 127, 156801 (2021)", "Walther, Willsher, Knolle, Phys. Rev. Lett. 136, 227102 (2026)"]

### N_blocked_r3 (first pass): verdict KNOWN
- **standard_name**: 2D long-range XY model at the marginal sigma=2 point (spin-wave magnetization decay as power of log L, log-corrected Goldstone phase)
- **model_family**: statistical physics: 2D long-range XY / O(2) model with power-law couplings 1/r^(2+sigma), sigma=2 marginal case (BKT / long-range-to-short-range crossover)
- **notes**: Could not read the full PRB text, so the specific numerical values (e.g. nu-like divergence exponent 2.24, exact (ln L)^-x exponents) were not checked against it. The product relation 2*pi*slope*exponent=1 is my own spin-wave derivation, not a quoted result. The nearest-neighbour+weak-tail variant is expected to be the same universality at sigma=2 with non-universal prefactors changed.
- **citations**: ["Yao, Xiao, Zhang, Deng, Fan 2025, Phys. Rev. B (DOI 10.1103/yghs-p5mg; arXiv:2411.01811) - MC study of 2D LR XY: LRO for sigma<=2, correlation saturates as inverse log of distance at sigma=2, crossover to SR at sigma=2", "Giachetti, Defenu, Ruffo, Trombettoni 2021, Phys. Rev. Lett. 127, 156801 - BKT transition with r^(-2-sigma) couplings; quasi-ordered phase for 7/4<sigma<2", "Kosterlitz 1976, Phys. Rev. Lett. 37, 1577 - RG for long-range n-component chains near the r^-2 marginal range (1D analogue of the sigma=2 boundary)"]
