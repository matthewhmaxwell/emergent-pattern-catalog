# planted:N02

## Target paper (ground truth; the gate never saw this)
- Title: Chemotaxis-induced linear instabilities and pattern formation in a reaction-diffusion model
- arXiv: 2609.01159 (v1 2026-09-01); authors: Mintu Karmakar, Abhik Basu
- Behaviour the authors report as new: In a two-species reaction-diffusion (glycolysis-type) model with mutual chemotaxis, steady Turing pattern amplitude varies non-monotonically with chemotactic strength (re-entrant uniform/patterned transitions) and spots turn into stripes as chemotaxis grows.
- Authors' novelty statement: "This has several novel features, which lead to some new physical phenomena."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: two continuous concentration fields u and v on a 128x128 square grid with wrap-around edges, started at their uniform steady state (u=b, v=b/(a+b^2)) plus small random noise. Rules, integrated forward in time: u changes at rate -u + a*v + u^2*v + D1*(Laplacian of u) + c1*div(u grad v); v changes at rate b - a*v - u^2*v + D2*(Laplacian of v) + c2*div(v grad u). The last terms make each species drift along the other's gradient, with strength and sign set by c1 and c2. Parameters: a=0.015, b=10, D1=0.01, D2=1000*D1. Main sweep: c1=c and c2=-c (opposite signs, so one species is drawn toward the other while the other is pushed away), c from 0.01 to about 0.5.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Same rules with D1=0.1 and D2=100.
  - Same rules with b=5, D1=0.03 and c1=-c2 negative: spots turn into stripes, then fade, as |c| grows.
  - Same rules with c1=c2=c (equal signs).

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - coupling_at_sharpest_pattern = 0.17  (drift coupling c where the steady spot pattern is most intense)
  - coupling_pattern_gone = 0.49  (drift coupling c at which the pattern has vanished (uniform))
  - amplitude_maxima_vs_coupling = 1  (maxima of pattern amplitude as coupling increases (rise then fall))
  - spot_to_stripe_transitions = 1  (spot-to-stripe morphology changes as coupling magnitude increases)
  - seeds_with_same_outcome = 10  (independent random seeds reproducing the morphology sequence)

This behaviour appears across 4 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r2 (strict pass): verdict KNOWN
- **standard_name**: Cross-diffusion-modulated Turing patterns (spot→stripe conversion and annihilation) in the Sel'kov (Sel'kov–Schnakenberg) reaction-diffusion model
- **behaviour_documented**: True
- **notes**: Upgrades the first pass from KNOWN-VARIANT to KNOWN and corrects the model name: the reaction terms are Sel'kov kinetics, not Schnakenberg. Cross-diffusion in Sel'kov–Schnakenberg is documented as controlling pattern type, converting spots to stripes and annihilating patterns. A rough linearisation of my own predicts a rise-then-fall in pattern strength: opposite-signed coupling adds a destabilising term proportional to -c k^2 and a stabilising term proportional to +c^2 k^4. The same analysis does not clearly predict a return to the uniform state at c=0.49; the fastest-growing pattern stays at high wavenumber and may not be resolved on the grid. That threshold could be an artifact of the grid or time step, or of lost positivity, NaN or clamped values, and should be checked by refining the grid and time step. Spots, the spot-to-stripe change and the amplitude peak look genuine. Sources were read as snippets and abstracts only, not full texts.
- **citations**: ["https://link.springer.com/article/10.1007/s10910-023-01452-0", "https://www.sciencedirect.com/science/article/abs/pii/S0016003225003692", "https://pubs.aip.org/aip/cha/article-abstract/35/7/073142/3355622", "https://www.sciencedirect.com/science/article/pii/S0307904X17306960", "https://www.sciencedirect.com/science/article/abs/pii/S0960077924009275", "https://www.sciencedirect.com/science/article/abs/pii/S0960077916302259"]
