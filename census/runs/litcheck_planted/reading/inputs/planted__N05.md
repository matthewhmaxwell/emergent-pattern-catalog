# planted:N05

**Reported behavior (what the authors report as new):** In a 1D two-species driven lattice gas on a fluctuating landscape, early ordered domains disintegrate and new ones form, giving two coarsening exponents (0.5 then 0.34) and a non-monotonic landscape length scale.

**Rule on the card:** World: a ring of 16,416 sites; between neighbouring sites is a bond tilted up or down (half of each, random), defining a height profile. Each site holds an A particle, a B particle or nothing (one third each, random). Rules, random sequential updates: (1) a particle swaps with a neighbouring site of different content; a move taking A downhill or B uphill succeeds with probability 0.9, the opposite move with 0.1; (2) a local peak (up-bond then down-bond) and a dip swap: a peak holding A becomes a dip with probability 0.5+b and back with 0.5-b, b=-0.4 (so A pushes ground up); a dip holding B becomes a peak with probability 0.5+b' and back with 0.5-b', b'=0.3; empty sites flip either way with 0.5. Track domain size over time.

**Observed features on the card:**
  - early_growth_exponent = 0.5 (power-law exponent of domain size growth before t~1000)
  - late_growth_exponent = 0.34 (power-law exponent of domain size growth at late times)
  - crossover_time = 1000 (time at which early domains dissolve and growth slows)
  - late_macroscopic_hills = 1 (macroscopic hills in the height profile in the late regime)

## Reply R03

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Coupled two-species particle-landscape model (light-heavy model) of Chakraborty-Chatterjee-Barma", "notes": "Only abstracts and reference lists were accessible, not full texts. Model identity and the macroscopic-hill phase are well supported. The specific early exponent 0.5, crossover near t~1000 and late exponent 0.34 (0.41 at b'=0.18) were not confirmed in any source I could read. The late exponent may be an effective or crossover value near a phase boundary."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The model and its ordered phase are described in the literature. I could not verify the specific two-stage coarsening numbers, namely the early exponent of 0.5, the crossover near t≈1000 and the late exponent of 0.34.

**1. Model family.** This is a driven lattice gas of two particle species plus vacancies, coupled two-way to a stochastic interface or landscape. It belongs to nonequilibrium statistical mechanics: "sliding particles on a fluctuating surface" and fluctuation-dominated phase ordering.

The two-way coupling is the key feature. Particles slide downhill on the height profile, and each species also biases the local peak/dip dynamics. In your rules that bias is b=-0.4 for A and b'=0.3 for B.

**2. Closest matches in the literature**
- **Earliest clear source for the exact model.** Chakraborty, Chatterjee and Barma, "Ordered phases in coupled nonequilibrium systems: Static properties", Phys. Rev. E 96, 022127 (2017).
  - It studies a coupled driven system in which two species of particles are advected by a fluctuating potential energy landscape.
  - In all the ordered phases the two species phase separate completely. The landscape either orders fully, with macroscopic regions of distinct average slope, or shows coexisting ordered and disordered regions. Which one occurs depends on the differential effect of the two species on the landscape.
  - This matches your "macroscopic hills" feature.
- **Companion paper.** Chakraborty, Chatterjee and Barma, "Ordered phases in coupled nonequilibrium systems: Dynamic properties", Phys. Rev. E 96, 022128 (2017). It shows that macroscopic particle clusters move over an ergodic time scale growing exponentially with system size. The ordered landscape moves faster, on a time scale growing as a power of system size.
- **Earlier letter by the same group.** The reference list of a later paper cites Chakraborty, Pal, Chatterjee and Barma, "Large compact clusters and fast dynamics in coupled nonequilibrium systems", Phys. Rev. E 93, 050102(R) (2016). I saw it only in that reference list, not as an abstract.
- **Most recent study of the coarsening.** Mahapatra, Ramola and Barma, "Light and heavy particles on a fluctuating surface: Bunchwise balance, irreducible sequences, and local density-height correlations", Phys. Rev. Research 2, 043279 (2020).
  - It studies the early-time and coarsening dynamics in the light-heavy model, with two species of particles coupled to a fluctuating surface described by tilt fields.
  - It predicts an intermediate diffusive (power-law) stretch that precedes the onset of coarsening, and it finds this in simulations of the ordered regime.
  - This is qualitatively like your early t^0.5 stage followed by slower coarsening.
- **Root of the field (single species, passive, one-way coupling).** Das and Barma, "Particles sliding on a fluctuating surface: phase separation and power laws", Phys. Rev. Lett. 85, 1602 (2000). Initially random particle density coarsens and obeys scaling with a growing length scale ∼t^{1/z}. A related paper is Das, Barma and Majumdar, Phys. Rev. E 64, 046126 (2001), whose title I confirmed only from a reference list.
- **Review.** I found no formal review. Barma's group page summarises the field, noting that when particles sliding on a fluctuating surface affect and are affected by the surface dynamics, there are phase transitions to a set of ordered phases.

**3. What differs and whether it matters.** I only had abstracts and reference lists, not full texts. So I could not confirm that these papers report the following for your parameters (b=-0.4, b'=0.3, one third each of A, B and vacancies):
- an early exponent of about 0.5;
- a dissolution of early domains near t≈1000;
- a late exponent of about 0.34, or a slower value of about 0.41 at b'=0.18.

The "dynamics" and "light-heavy" papers probably cover this, but that is unconfirmed. Your parameter values could also sit near a phase boundary between the fully ordered and mixed ordered/disordered landscape phases. If so, the 0.34 versus 0.41 shift may be a crossover effect rather than a distinct universal exponent. That is a plausible caveat, not a finding.

**4. Confidence: medium.** I am fairly sure the model is the same one and that the macroscopic-hill phase is documented. Confidence is low on the quantitative exponents, because the full texts were not accessible. Check the Mahapatra et al. 2020 paper and the 2017 dynamics paper for the growth-exponent data.

**

## Reply R16

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Coupled particle-landscape (light-heavy) model of Chakraborty-Chatterjee-Barma: strong phase separation / compact-cluster ordered phases with algebraic coarsening (not FDPO, which only occurs on the separatrix)", "notes": "Agrees with first pass on the model family (KNOWN-VARIANT), but corrects the label: generic couplings (b=-0.4, b'=0.3) put the system in an ordered phase-separated phase, not FDPO, which lives only on the order-disorder separatrix. Unverified twist: the explicit early t^0.5 -> late t^~0.34 crossover near t~1000 with early-domain dissolution; full texts were not accessible. The late exponent shifts with b' (0.34 vs 0.41), suggesting effective/crossover exponents rather than a new class. Mapping of these parameters to the compact-cluster phase is inferred, not verified."}

**Reply text:**

**Verdict: KNOWN-VARIANT, at medium confidence.** The behaviour (two particle species plus vacancies on a landscape they also reshape, phase separating into macroscopic clusters and macroscopic landscape order) is documented. The exact 0.5 → 0.34 growth crossover near t≈1000 is not confirmed in anything I could open. I partly overturn the first pass on the name: this is not best described as fluctuation-dominated phase ordering (FDPO).

**What this system is.** It is the coupled particle–landscape model of Chakraborty, Pal, Chatterjee and Barma, which the Mahapatra–Ramola–Barma line calls the light–heavy (LH) model. Your A particles are the heavy species (they slide downhill) and your B particles are the light species (they slide uphill). The documented phenomenology matches yours in kind:
- **Particles and landscape:** two species of particles are advected by a fluctuating landscape and modify the landscape in the process. The phase diagram is valid for all particle densities, so vacancies are covered.
- **Ordered states:** the two particle species phase separate completely from each other, but the underlying landscape may either show complete ordering, with macroscopic regions with distinct average slopes, or may show coexistence of ordered and disordered regions. Your late-time "macroscopic hills" fit this.
- **Coarsening:** in one phase, ordered landscape coexists with disorder, giving coarsening and steady state dynamics on time scales which grow algebraically with size. Algebraic coarsening is consistent with a late exponent around 1/3, i.e. dynamic exponent z≈3, but I only read the abstract and could not confirm that value.
- **Early times and a crossover:** Mahapatra et al. (2020) study the early-time and coarsening dynamics in the light-heavy model. They find an intermediate diffusive (power-law) stretch before coarsening starts, and identify nonlocal "irreducible sequences" that provide a coarsening length scale in the ordered phases at late times. A diffusive early stretch fits your exponent of 0.5.

**Correction to the first pass.** FDPO is not this model's generic ordered state. the separatrix between the disordered and ordered phases is a surface in parameter space on which the system exhibits fluctuation-dominated phase ordering. Your couplings (b = −0.4, b′ = 0.3) are well off that boundary. You are almost certainly in one of the strongly phase-separated ordered phases. Because A builds up the ground it slides down while B builds up the ground it climbs, the most likely phase is the "compact cluster" one, where ordered and disordered landscape coexist. That last step is my inference from the abstracts; I could not map your parameters onto the papers' phase diagram.

**The unverified twist.** The specific 1/2 → ~1/3 crossover at t≈1000, with early domains dissolving, is not stated in any abstract I could access. I could not open the full texts of Chakraborty et al. 2016/2017 or Mahapatra et al. 2020. This probably doesn't matter much:
- The late exponent moves with b′ (0.34 here, 0.41 at b′ = 0.18), so these look like effective, finite-time exponents rather than a new universality class.
- A diffusive early regime followed by slower coarsening is what the 2020 paper describes.

**Is it real or an artifact?** It looks like a genuine emergent pattern: clear growth regimes, a crossover, and macroscopic landscape order. Nothing suggests random colouring, an inert coupled field, numerical divergence or a frozen absorbing state.

**

## Reply R27

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Light-heavy (LH) model on a fluctuating surface: coarsening to macroscopic phase separation (FDPO)", "notes": "Only abstracts and snippets were read. The exact rate parametrization (0.9/0.1, b=-0.4, b'=0.3) and the specific exponents 0.5 -> 0.34 (0.41 for b'=0.18) with crossover near t=1000 were not confirmed in any source. 0.5 likely reflects z=2 surface diffusion and ~1/3 the conserved-order-parameter value, but this is my inference. Check the full text of Mahapatra et al. 2020 to close the gap."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The model family is known, and the coarsening mechanism in it is documented. I could not read any source that gives the specific exponents you observed (0.5 early, about 0.34 late, crossover near t≈1000), so those numbers are unconfirmed.

**1. Model family.** This is a driven or nonequilibrium lattice gas on a fluctuating surface. It has two species of hard-core particles (A and B) plus vacancies, coupled to a tilt or height field. The particles feel the local slope, and the particles in turn bias how the slope changes. The literature calls this family "fluctuation-dominated phase ordering" (FDPO). The closest named match is the **light–heavy (LH) model**. The surface-coupled lattice models in this family include "sedimenting lattices", sliding particles and the LH model.

**2. Specific behaviour and sources.**
- **Earliest clear source for the family.** Das & Barma, "Particles Sliding on a Fluctuating Surface: Phase Separation and Power Laws", Phys. Rev. Lett. 85, 1602 (2000). Hard-core particles slide locally downwards on a fluctuating 1D surface, and an initially random density coarsens with a growing length scale ~t^{1/z}. The result for the one-species model is that the exponent is set by the surface's dynamical exponent z. For a surface with z=2 this gives 1/2. Das, Barma & Majumdar (Phys. Rev. E 64, 046126, 2001) extend this. They report surface-fluctuation-driven coarsening with length scale growing as t^{1/z}.
- **Two coupled species with a tilt field.** Lahiri, Barma & Ramaswamy, "Strong phase separation in a model of sedimenting lattices", Phys. Rev. E 61, 1648 (2000). I identified the authors from memory and from the PDF's location on S. Ramaswamy's page. The search result did not show them. It maps the problem to a 1D lattice model with two coupled species of spins, and each species shows macroscopic phase separation in the steady state. This matches your "macroscopic hills" in the late regime.
- **Closest equivalent system.** Mahapatra, Ramola & Barma, "Light and heavy particles on a fluctuating surface: Bunchwise balance, irreducible sequences, and local density-height correlations", Phys. Rev. Research 2, 043279 (2020). It studies early-time and coarsening dynamics of two species coupled to a fluctuating surface described by tilt fields. It predicts an intermediate diffusive (power-law) stretch before coarsening begins, and finds it in simulations of the ordered regime. That is qualitatively like your early t^0.5 regime followed by a crossover. It also says irreducible sequences give a coarsening length scale in the ordered phases at late times. The steady-state phases, ordered or disordered, depend on the microscopic rates. That fits your two parameter sets giving different late exponents (0.41 and 0.34).
- **Review.** M. Barma, "Fluctuation-dominated Phase Ordering", in *50 Years of the Renormalization Group* (World Scientific, 2024). It appears in his publication list, but I did not read it. Barma, Majumdar & Mukamel, J. Phys. A 52, 254001 (2019) is a related overview. I only saw the title in the same list.
- A later paper, Prakash, Barma & Ramola, J. Stat. Mech. 103208 (2025), concerns the *scaled* LH model, which has a different phase diagram. I did not use it.

**3. What differs, and whether it matters.**
- I could not confirm that the published LH model uses exactly your rules: the 0.9/0.1 slope-biased hops, the peak/dip flips biased by b = −0.4 and b′ = 0.3, and vacancies at one third. Matching these details would need the Mahapatra et al. paper itself, which I did not read.
- The abstracts I saw do not give a late-time growth exponent near 1/3 or 0.41, or a crossover time near 10³. The late exponent 0.34 may simply be the conserved-order-parameter value 1/3. That is plausible but unverified. The standard textbook value for a conserved order parameter is t^{1/3} (Redner's coarsening chapter says conserved-order-parameter coarsening gives t^{1/3} growth). The drift to 0.34 and 0.41 may be finite-time effects.
- The early t^0.5 regime is consistent with z=2 diffusive surface dynamics, because t^{1/z} with z=2 gives 1/2. This is my inference, not something I found stated.
- It is unclear whether "early domains dissolve" at t≈1000 is the same event the LH paper calls the "onset of coarsening". The hills are presumably the macroscopic phase-separated state of the ordered phase.

**4. Confidence: medium-low.** The model family and the qualitative phenomenon (a surface-driven early regime, then coarsening, then macroscopic phase separation) are well documented. I only read abstracts and snippets, not the full texts. I found no quantitative statement of your exponents or of the crossover. That is why I chose KNOWN-VARIANT instead of KNOWN.

**

## Reply R41

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Light-heavy (LH) model / particles advected by and modifying a fluctuating landscape: early diffusive stretch followed by late-time coarsening toward ordered phases with macroscopic landscape order (or ordered/disordered coexistence)", "notes": "Documented qualitatively: an intermediate diffusive power-law stretch before coarsening (consistent with exponent 0.5), a correlation function with a broad maximum followed by late-time power-law decay, and macroscopic landscape ordering or coexistence. Not verified from full text: a late domain-growth exponent of about 1/3 and a crossover near t~1000. The exponent 0.41 reported for b'=0.18 suggests effective, parameter-dependent exponents. Both species push the surface up here, which may place the system in the documented coexistence phase or near the separatrix rather than a fully ordered phase, so the first pass's FDPO label may be imprecise. Not trivial: there is real coarsening, not a frozen absorbing state or a diverged field."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I'm keeping the first pass's verdict, with medium confidence. Each part of the behaviour is described in the light-heavy (LH) model papers. The one thing I couldn't confirm is the exact numbers: a late exponent of about 0.34 (roughly 1/3) and a crossover near t≈1000. I couldn't open the full text of any paper, including Mahapatra et al. 2020, so the gap the first pass flagged is still open.

**It's the same system.** Mahapatra, Ramola and Barma (2020) study the LH model, "a system consisting of two species of particles (light and heavy) coupled to a fluctuating surface (described by tilt fields)". They report that particle and tilt updates are coupled locally and give different ordered and disordered steady states depending on the rates. Your ring of up/down bonds, A/B/empty sites, particle hops biased by the slope, and peak↔dip flips that depend on which particle sits there is this model family. It was first studied by Chakraborty, Pal, Chatterjee and Barma (PRE 93, 050102(R), 2016).

**How each observed feature matches:**
- **Early exponent ≈ 0.5.** The 2020 paper predicts that between a lattice-cutoff timescale and the start of coarsening there is an intermediate diffusive (power-law) stretch, and says they see it in simulations of the ordered regime. A diffusive stretch fits a t^½ length.
- **Early domains dissolving, then slow coarsening.** In the same paper, a local correlation function rises linearly from a disordered start and has a broad maximum. It then decays as a power law at late times, which encodes the coarsening toward the ordered phases. That is a short-time build-up followed by a different late-time power law, which looks like your crossover.
- **Macroscopic hill.** The 2017 static-properties paper says that in all ordered phases the two species separate completely, while the landscape either orders fully into macroscopic regions with distinct slopes, or shows ordered regions coexisting with disordered ones, depending on how each species reshapes the surface. The 2016 paper links this coexistence to coarsening and steady-state dynamics on time scales that grow algebraically with size.

**What I couldn't confirm.** None of the abstracts I could reach gives a late-time domain-growth exponent of 1/3, or any crossover time. A 1/3 is plausible for diffusion-limited coarsening with conserved order, but that is my own inference, not something I verified. The exponent of 0.41 the first pass mentions for b'=0.18 also suggests the late exponent is an effective, parameter-dependent value rather than a clean 1/3. If so, that is ordinary crossover or finite-time behaviour, not new physics. It matters only if someone wants to claim a new universality class, and nothing here supports that.

**One possible mislabel.** Here both A and B push the surface up (b=−0.4 at peaks, b'=+0.3 at dips). Under the published sign conventions, that pairing may put the system in the coexistence phase or near the boundary where fluctuation-dominated ordering occurs, not in a fully ordered phase. "Coarsening to macroscopic phase separation (FDPO)" may therefore be the wrong phase name. The behaviour is still within the documented phase diagram.

**Not trivial or an artifact.** The growth changes regime, a large hill forms, and the exponents are finite and sensible. That is real coarsening, not a frozen state, a random pattern or a diverged run.

**Searches I ran:** LH-model coarsening crossover with macroscopic hills; the Chakraborty, Pal, Chatterjee and Barma 2016 paper; the Mahapatra, Ramola and Barma 2020 paper (and its arXiv entry, 2007.02971) with the exponents; the 2017 dynamic-properties paper; and a 1/3 exponent in the LH landscape.
