# planted:N02

**Reported behavior (what the authors report as new):** In a two-species reaction-diffusion (glycolysis-type) model with mutual chemotaxis, steady Turing pattern amplitude varies non-monotonically with chemotactic strength (re-entrant uniform/patterned transitions) and spots turn into stripes as chemotaxis grows.

**Rule on the card:** World: two continuous concentration fields u and v on a 128x128 square grid with wrap-around edges, started at their uniform steady state (u=b, v=b/(a+b^2)) plus small random noise. Rules, integrated forward in time: u changes at rate -u + a*v + u^2*v + D1*(Laplacian of u) + c1*div(u grad v); v changes at rate b - a*v - u^2*v + D2*(Laplacian of v) + c2*div(v grad u). The last terms make each species drift along the other's gradient, with strength and sign set by c1 and c2. Parameters: a=0.015, b=10, D1=0.01, D2=1000*D1. Main sweep: c1=c and c2=-c (opposite signs, so one species is drawn toward the other while the other is pushed away), c from 0.01 to about 0.5.

**Observed features on the card:**
  - coupling_at_sharpest_pattern = 0.17 (drift coupling c where the steady spot pattern is most intense)
  - coupling_pattern_gone = 0.49 (drift coupling c at which the pattern has vanished (uniform))
  - amplitude_maxima_vs_coupling = 1 (maxima of pattern amplitude as coupling increases (rise then fall))
  - spot_to_stripe_transitions = 1 (spot-to-stripe morphology changes as coupling magnitude increases)
  - seeds_with_same_outcome = 10 (independent random seeds reproducing the morphology sequence)

## Reply R28

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "cross-diffusion-induced/modified Turing patterns in the Schnakenberg model (pattern morphology and amplitude vs cross-diffusion strength)", "notes": "No source found for the exact antisymmetric nonlinear-taxis sweep on Schnakenberg with a=0.015, b=10, D2/D1=1000. Qualitative phenomena (cross-diffusion alters or erases Turing patterns, spot/stripe changes, non-monotonic effect of coupling) are documented in close relatives. Most Schnakenberg cross-diffusion papers use linear constant-coefficient cross-diffusion, which could matter for large-amplitude behaviour. Some cited papers were seen only as abstracts or snippets."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium-low confidence.** This is the Schnakenberg Turing system with added nonlinear cross-taxis terms. Its spot and stripe patterns, and the fact that cross-diffusion can enhance or suppress them, are well documented. I found no source that reports your exact setup: an antisymmetric ±c taxis pair, a rise-then-fall pattern amplitude, and a spot-to-stripe-to-uniform sequence.

**1. Model family**
- It is a reaction-diffusion (Turing) system with Schnakenberg kinetics, an activator–substrate-depletion scheme.
- The terms c1·div(u∇v) and c2·div(v∇u) are nonlinear cross-diffusion, or cross-taxis. The same form appears in chemotaxis and in diffusiophoresis models.
- With c2 = −c1, the coupling is nonreciprocal: one species is attracted to the other while the other is repelled.

**2. What the literature shows**
- **Schnakenberg with linear cross-diffusion (closest earliest source).** Madzvamuse, Ndakwo and Barreira (2014/15, *J. Math. Biol.* 70:709) derive cross-diffusion-driven instability conditions for Gierer–Meinhardt, Prigogine–Lefever and Schnakenberg kinetics. Unlike classical systems, one species no longer has to diffuse much faster than the other. The cross-diffusion-induced parameter spaces are substantially larger and different from those without cross-diffusion.
- **Weakly nonlinear analysis of the same system.** Gambino/Lombardo-type follow-up work by Ricerche di Matematica (2016) found that the diffusion coefficients influence whether the Turing bifurcation is supercritical or subcritical. The patterns it classifies are squares, rhombi, rectangles and hexagons.
- **Spots, stripes and mixed patterns as the cross-diffusion coefficient varies.** Schnakenberg-type papers report this directly. Examples are the Nonlinear Dynamics 2022 paper "Cross-diffusion induced spatiotemporal patterns in Schnakenberg reaction–diffusion model" (hexagon-like, stripe-like and coexistence patterns) and a fractional Sel'kov–Schnakenberg study in which the cross-diffusion coefficient controls the type of Turing pattern.
- **Cross-diffusion can erase patterns.** In a ratio-dependent prey-predator model with nonlinear cross-diffusion, the cold-spot pattern may turn into a homogeneous steady state, a stripe/spot mixture or a labyrinth. I could not read the authors, journal or year from the search result, so treat this citation as incomplete.
- **Closest physical analogue, with the same div(u∇v) structure.** Menou, Luo and Zwicker (2023, *J. R. Soc. Interface* 20:20230244) add interactions between activator and inhibitor to Turing systems. Weak repulsive interactions stabilize patterns by inducing cross-diffusion, while strong interactions lead to phase separation. Repulsion promotes pattern formation and attraction suppresses it. This is a non-monotonic dependence on coupling strength. It differs from your system because the interaction is reciprocal and derived from a free energy, and the kinetics are Hill–Langmuir rather than Schnakenberg.
- **Reviews and general frameworks.** I found no dedicated review of nonreciprocal Schnakenberg taxis. The nearest are Madzvamuse, Barreira and Gerisch (2017), "Cross-Diffusion in Reaction-Diffusion Models: Analysis, Numerics, and Applications" (ECMI 2016 proceedings), and Villar-Sepúlveda and Champneys's bilayer reaction-cross-diffusion paper (*Bull. Math. Biol.* 2023, bioRxiv 2023.05.30.542795). I only saw the Madzvamuse reference as a citation in other papers.
- **Chemotaxis coupled to Turing systems.** "Novel Aspects in Pattern Formation Arise from Coupling Turing Reaction–Diffusion and Chemotaxis" (*Bull. Math. Biol.* 2023) is a relative, not an equivalent.

**3. What differs and whether it matters**
- Most Schnakenberg cross-diffusion papers use linear, constant-coefficient cross-diffusion. Yours is nonlinear, with fluxes proportional to u∇v and v∇u. Linear stability at the uniform state depends only on c·u* and c·v*. Late-stage behaviour such as amplitude saturation and the spot-to-stripe change can differ, because of the nonlinear flux and the depletion of the uniform state.
- I found no source describing the opposite-sign (nonreciprocal) pairing with a peak in pattern amplitude near c ≈ 0.17 and uniformity by c ≈ 0.49. A rise-then-fall of amplitude, and a spot-to-stripe-to-uniform sequence with increasing coupling, are qualitatively consistent with the cross-diffusion literature above. The specific numbers (0.17, 0.49) are parameter-dependent and not expected to appear in any paper.
- The parameters (a=0.015, b=10, D2/D1=1000) are unusual, and I did not find them in the literature.

**4. Confidence**
- I am confident that the model class and the qualitative phenomena are known.
- I am less sure about the exact antisymmetric nonlinear-taxis sweep. I ran only 10 searches and did not read the full text of any paper that sweeps this coupling in Schnakenberg. A relevant paper on nonreciprocal or chemotactic Turing systems could have been missed.
-
