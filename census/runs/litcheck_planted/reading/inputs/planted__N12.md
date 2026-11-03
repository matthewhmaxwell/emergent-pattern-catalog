# planted:N12

**Reported behavior (what the authors report as new):** The 2D XY model with 1/r^4 couplings has a low-temperature phase whose correlations decay as a power of log r, entered through a continuous transition with power-law divergent correlation length (neither BKT nor Landau type).

**Rule on the card:** World: an L x L square lattice with wrap-around edges (L up to 2048), each site holding a planar unit spin (an angle). Energy: minus the sum over all pairs of J(r) times the cosine of their angle difference, where every pair interacts with J(r) proportional to 1/r^4 (r = distance), normalised so that the total coupling per site is fixed. Rules: equilibrium Monte Carlo with cluster updates at inverse temperatures beta = 1, 2, 4, 8 and near the ordering point. Measure the squared magnetisation and the correlation length versus L.

**Observed features on the card:**
  - log_decay_exponent_beta1 = 0.34 (magnetisation decays as power of log L; exponent at beta=1)
  - log_decay_exponent_beta8 = 0.02 (same exponent at beta=8 (proportional to temperature))
  - corr_length_log_slope_beta1 = 0.48 (slope of squared correlation-length ratio versus log L, beta=1)
  - exponent_slope_product = 1.0 (2*pi times slope times exponent; 1 means spin-wave picture holds)
  - corr_length_divergence_exponent = 2.24 (power-law divergence exponent of correlation length at the transition)

## Reply R06

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "2D long-range XY model at marginal sigma=2 (1/r^4): logarithmic low-temperature scaling from Goldstone spin waves with k^2 ln(1/k) dispersion", "notes": "Logarithmic scaling at sigma=2 is documented by the Deng group, who read it as true long-range order (constant + a/ln r). The observed reading is a power of ln L going to zero, with exponent proportional to T and 2*pi*slope*exponent = 1. That fits Gaussian spin-wave theory and Bruno's theorem (no order for alpha >= 2D), but I found it nowhere stated explicitly. L up to 2048 (ln L about 7.6 at most) cannot separate the two forms. The beta=8 exponent is noisy; the ratio of about 17 instead of 8 fits a stiffness reduced at beta=1. The 2.24 divergence exponent is unmatched and probably an effective power-law fit at a boundary case with log corrections. Not an artifact."}

**Reply text:**

I rate this **KNOWN-VARIANT**, which matches the first pass's label, though I've changed parts of its reasoning. The behaviour itself is documented: logarithmic finite-size scaling in the low-temperature phase of the 2D long-range XY model at exactly σ=2 (J ∝ 1/r⁴). What I couldn't find in print is the specific reading your simulation suggests. That reading is magnetisation falling to zero as a power of ln L, with an exponent proportional to temperature and the spin-wave "product = 1" check. The two groups that simulated this model read the same logarithmic scaling as true long-range order.

**What is documented**
- **Logarithmic scaling at σ=2.** Xiao, Yao, Zhang, Fan and Deng ran cluster Monte Carlo up to L = 8192 and find logarithmic scaling behaviors in the low-temperature phase at σ = 2 (Chin. Phys. Lett. 42, 070002, 2025).
- **Their interpretation.** The follow-up PRB (112, 144429, 2025) says that for σ≤2, the system undergoes a second-order phase transition into a ferromagnetic phase characterized by the emergence of long-range order. It also says the correlation function settles to a non-zero constant in the form of the inverse logarithm of distance for σ=2. In other words, it tends to a constant plus a/ln r, rather than to zero.
- **The rigorous bound favours your reading.** Bruno's theorem (PRL 87, 137203) states that for monotonically decreasing interactions, ferro- or antiferromagnetic long-range order at T>0 is ruled out if α≥2D. Here α = 4 and D = 2, so true order is excluded at σ=2. A review of the Deng group's related Heisenberg paper makes the same point: if the σ=2 claim were correct, it would overturn the extended Mermin–Wagner theorem as cited.
- **Theory context.** Giachetti, Defenu, Ruffo and Trombettoni found that for any σ>2 the BKT mechanism occurs, so σ=2 is the boundary case. Walther, Willsher and Knolle (PRL 2026) argue that the BKT transition persists for arbitrary long-range exponents.
- **Data can't easily separate the two readings.** In a related logarithmic-decay dispute, a reviewer noted that when ln L only spans about 3 to 6, a critical phase with G_y ~ (ln L)^{-p}, p small, would also be well fit by a + b L^{-q} with a positive intercept. Your L ≤ 2048 gives ln L of about 7.6 at most, so the same problem applies.

**The twist and whether it matters**
- **Not found in print:** an explicit statement that the magnetisation decays as (ln L)^(−T/2πJ_eff) at σ=2 with no true order. This follows from a standard Gaussian spin-wave calculation: the k² ln(1/k) dispersion gives an angle variance growing like (T/2πJ) ln ln L. Your 2π × slope × exponent ≈ 1 check also isn't in the papers I found.
- **Why it plausibly matters:** this phenomenology is the order-free alternative to the long-range-order claim, and it fits Bruno's theorem. But it is a short derivation, not a new effect, and the available sizes can't rule out "constant + a/ln L".
- **Temperature scaling:** the exponents at β=1 and β=8 differ by a factor of about 17, not 8. That is roughly what you'd expect if the effective stiffness is reduced at higher temperature, so it doesn't break the spin-wave picture, but the β=8 value is noisy.
- **The 2.24 divergence exponent:** I found no literature match. Since σ=2 is a boundary case with logarithmic corrections, a power-law fit there is likely an effective exponent and shouldn't be over-interpreted.

**Genuine or artifact?** Genuine. The order parameter decays slowly and smoothly, the exponent scales with temperature, and the spin-wave consistency check holds. None of that looks like a frozen, numerically diverged or inert state.

## Reply R20

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "2D long-range XY model at the marginal point sigma=2 (J~1/r^4): logarithmic finite-size behaviour from k^2 ln k Goldstone modes, ln ln L phase fluctuations", "notes": "Behaviour is documented. Twist: the observed m^2~(ln L)^(-eta(T)) decay to zero matches Bruno's theorem and spin-wave theory, but contradicts Yao/Xiao et al.'s LRO reading (c1+c2/ln L) of the same model; no paper found that explicitly fits (ln L)^(-eta(T)) with eta~T. This is an interpretive dispute, not new physics; the two forms are hard to separate at L<=2048. eta(beta=1)/eta(beta=8)=17, not 8, so eta~T holds only at low T. nu=2.24 not compared with the published sigma=2 value and may be an effective fit if the transition is out of a log-quasi-ordered phase. Only abstracts/summaries were read."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium confidence.** The first pass was right about what kind of system this is. But it left out an open dispute in the literature over whether this system orders at all, and the observed behaviour lands directly on that question.

**What is documented**
- **Same model, numerically.** Yao, Xiao, Zhang, Deng and Fan (PRB 112, 144429, October 2025) simulated this exact model: the 2D XY model with ~1/r^(2+σ) couplings, including σ=2. They report that in the low-temperature phase, due to the presence of the Goldstone mode, the correlation function saturates to a nonzero constant in the form of a power law for σ<2, with decaying exponent 2−σ, and in the form of the inverse logarithm of distance for σ=2. They also determined critical points and exponents for various σ. So slow, logarithmic finite-size behaviour at σ=2, plus a transition with a measured ν, is published for this system.
- **The theory.** Giachetti, Defenu, Ruffo and Trombettoni (PRL 127, 156801, 2021) analysed the renormalisation-group picture. They find the usual short-range behaviour, with a BKT transition, if σ>2, and an upper bound σ*=2 for the threshold above which the traditional BKT transition persists. That makes σ=2 the marginal point.
- **Rigorous result.** Bruno (PRL 87, 137203, 2001): for interactions that decrease monotonically, ferro- or antiferromagnetic long-range order at T>0 is ruled out if α≥2D. Here α=4=2D, so true long-range order is forbidden.
- **The same behaviour in an equivalent effective theory.** In 2D generalised nonlinear sigma models with a stiff phase, the variance of the transverse spin fluctuations scales with L as ln ln L. Equivalently, for N=2 in 2D, the equal-time spin-spin correlations decay as powers of the inverse logarithm of the separation, with model-dependent exponents. This is exactly the observed m² ~ (ln L)^(−η(T)). It is also what you get from the k² ln k spin-wave dispersion of a 1/r⁴ coupling. The relation 2π·slope·η = 1 is just the Gaussian spin-wave identity.

**The twist, and whether it matters**
Yao/Xiao et al. read their σ=2 data as true long-range order, <⟨M²⟩ = c₁ + c₂/ln(L/L₀)>, and the same group says the same for the Heisenberg case. A public critique states that the claim of finite-temperature long-range order at σ=2 in the 2D LR Heisenberg model contradicts Bruno's extension of the Mermin–Wagner theorem.

Your observation, a magnetisation that decays to zero as a power of ln L with an exponent that shrinks with temperature, matches Bruno and spin-wave theory, not the long-range-order reading. So the phenomenology isn't new. What it does is take a side in a live interpretive dispute. Over L ≤ 2048, c₁ + c₂/ln L and A(ln L)^(−η) are hard to tell apart numerically. I didn't find any paper that explicitly fits (ln L)^(−η(T)) with η ∝ T for this model.

**Points to check**
- **Exponent vs temperature.** η falls from 0.34 to 0.02 when β goes from 1 to 8. That's a factor of about 17, not the factor of 8 that "η ∝ T" implies. β=1 is probably close to the transition, where non-Gaussian corrections kick in. So "proportional to temperature" only holds at low T.
- **ν = 2.24.** I didn't compare this with Yao et al.'s published σ=2 value. If there really is no long-range order, the transition is out of a log-quasi-ordered phase. A finite power-law ν there may only be an effective fit to something closer to BKT-style (exponential) divergence.
- **Sources read.** I only saw abstracts and summaries, not full texts.

**Genuine or artifact:** genuine. Pure spin-wave theory with no free parameters reproduces the measured exponent-slope product of 1.0, and the trends are smooth and consistent with temperature. Nothing here looks like a frozen state, an inert coupled field or a numerical divergence.
