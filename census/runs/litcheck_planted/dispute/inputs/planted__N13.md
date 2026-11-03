# planted:N13

## Target paper (ground truth; the gate never saw this)
- Title: Opinion dynamics under large-scale exogenous events
- arXiv: 2609.39858 (v1 2026-09-30); authors: Gerard Argany Herrera, Lucila G. Alvarez-Zuzek, Oriol Artime
- Behaviour the authors report as new: A noisy voter model hit by random population-wide opinion shocks shows trimodal and asymmetric bimodal stationary opinion distributions in addition to the usual unimodal and bimodal ones.
- Authors' novelty statement: "For the case of the all-to-all noisy voter model under exogenous events, we unveil a rich phase diagram characterized by novel trimodal and asymmetric bimodal states."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: 500 agents, each holding opinion +1 or -1, every agent able to interact with every other. Rules each step: pick an agent; with probability eta (around 0.002, the inverse population; swept) it flips its opinion at random; otherwise it copies the opinion of a randomly chosen other agent. Separately, at rate r (swept), a population-wide shock occurs: opinions are instantly reassigned so that the average opinion jumps to 0.3 (65% holding +1). Measure the long-run distribution of the average opinion.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  (none distinct)

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - distribution_peak_count = 3  (peaks in long-run distribution of the average opinion)
  - middle_peak_position = 0.3  (average opinion at the extra, shock-induced peak)
  - distinct_regimes = 6  (distinct long-run regimes across the noise-versus-shock-rate plane)
  - noise_threshold_without_shocks = 0.002  (noise level separating one-peak and two-peak behaviour without shocks)

This behaviour appears across 1 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Mean-field noisy voter (Kirman) model finite-size bimodal-unimodal transition (a_c=2/(N+2)) + stochastic resetting of the magnetization (renewal-formula NESS, cf. Ising model with stochastic resetting)
- **behaviour_documented**: True
- **notes**: No-shock threshold matches a_c=2/(N+2) (~0.004 per flip-to-opposite, ~0.002 for random-flip convention at N=500). Third peak at the reset value 0.3 follows directly from the renewal formula when reset rate exceeds the O(1/N) voter drift rate, with consensus peaks at +/-1 surviving while noise < ~1/N; analogous rich (T,r) phase diagram documented for the Ising model with resetting. No paper found specifically on noisy voter/Kirman + global resetting; the 6-regime count is unverified and likely depends on how regimes were classified. Middle peak is an externally imposed imprint of the reset target, not self-organised: genuine but shallow.
- **citations**: ["Peralta/Carro/San Miguel/Toral, Herding and idiosyncratic choices... noisy voter model, C. R. Physique 2019 (S1631070519300313)", "Carro, Toral, San Miguel, The noisy voter model on complex networks, Sci. Rep. 2016", "Khalil/San Miguel/Toral, Zealots in the mean-field noisy voter model, PRE 97, 012310 (2018)", "Magoni, Majumdar, Schehr, Ising model with stochastic resetting, Phys. Rev. Research 2, 033182 (2020)", "Evans, Majumdar, Schehr, Stochastic resetting and applications, J. Phys. A 53, 193001 (2020)"]

### N_blocked_r2 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Mean-field noisy voter (Kirman) model finite-size bimodal-unimodal transition (a_c ~ h/N) under stochastic resetting of the magnetization (resetting-induced NESS)
- **behaviour_documented**: True
- **notes**: The noise threshold of 0.002 = 1/N simply reproduces the known Kirman/noisy-voter transition. The middle peak at 0.3 is the generic resetting-induced peak at the reset point, put in directly by the shock rule rather than emergent from the interactions; it shows up when the shock rate is comparable to the voter relaxation rate (~1/N per sweep). The closest documented equivalent is the Ising model with magnetization resetting (Magoni et al. 2020), which has a rich (T,r) phase diagram. No paper found on the exact noisy-voter-plus-resetting combination. The only twist is that the voter model's bistability is a finite-size effect with no deterministic pull, which may change regime boundaries but probably not the qualitative picture. The count of 6 regimes depends on how regimes are defined. Sources were checked through search snippets only, not full texts. Not a numerical artifact, but close to trivial.
- **citations**: ["Nature Sci. Rep. srep24775 - The noisy voter model on complex networks (finite-size bimodal/unimodal transition, a > h/N unimodal)", "Magoni, Majumdar, Schehr, Ising model with stochastic resetting, Phys. Rev. Research 2, 033182 (2020)", "Evans, Majumdar, Schehr, Stochastic resetting and applications, J. Phys. A 53, 193001 (2020)", "Noisy voter models in switching environments, Phys. Rev. E 108, 044301 (2023) (multipeak stationary distributions from external influencers)", "Herding and idiosyncratic choices: nonlinearity and aging-induced transitions in the noisy voter model, C. R. Physique (2019)"]

### N_blocked_r3 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Mean-field noisy voter (Kirman) model finite-size bimodal-unimodal transition (noise ~ 1/N) under global stochastic resetting of magnetization (resetting-induced peak at reset point)
- **behaviour_documented**: True
- **notes**: Noise threshold ~1/N reproduces the known Kirman criterion (1/N vs 2/N depends on whether 'random flip' is forced or a random reassignment). The third peak at the reset value 0.3 is the generic resetting signature (renewal/Laplace-transform form of the NESS). The specific 6-regime count across the noise-shock plane is unverified and likely reflects combinations of uni/bimodal background, presence of the reset peak, and which peak dominates (shock rate vs relaxation rate ~ max(eta, 1/N)); it is not a qualitatively new phenomenon. Not an artifact: the system stays ergodic and the shock peak is real, but it is expected by construction.
- **citations**: ["Granovsky & Madras 1995, The noisy voter model, Stoch. Proc. Appl. 55:23", "Kirman 1993 ant/herding model; Carro, Toral, San Miguel, Sci. Rep. 6:24775 (2016) noisy voter on networks", "Magoni, Majumdar, Schehr, Ising model with stochastic resetting, Phys. Rev. Research 2, 033182 (2020)", "Evans, Majumdar, Schehr, Stochastic resetting and applications, J. Phys. A 53, 193001 (2020)", "Gupta & Jayannavar, Stochastic resetting in interacting particle systems: a review, J. Phys. A (2023)", "Llabrés, San Miguel, Toral, Partisan voter model: noise-induced transitions (trimodal), PRE 2023"]
