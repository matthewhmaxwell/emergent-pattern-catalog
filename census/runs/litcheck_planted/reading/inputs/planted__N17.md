# planted:N17

**Reported behavior (what the authors report as new):** In the coexistence phase of the 2D Vicsek model, side-by-side lanes of band trains travelling in opposite directions (and fat/inverse bands) form at particle level; long lanes buckle and the ordered band states are metastable.

**Rule on the card:** World: self-propelled point particles in a periodic 2D box (square of side 1200, or rectangular), density 1 per unit area. Rules, discrete time steps: (1) each particle sets its heading to the average heading of all particles within distance 1 (itself included), plus a uniformly random angle in [-eta*pi, eta*pi], with noise eta between 0.28 and 0.35; (2) it moves 0.5 along its new heading. Start from random positions and headings, or from alternating parallel strips of opposite headings. Run up to 10^8 steps and watch the dense travelling stripes that form.

**Observed features on the card:**
  - parallel_lane_count = 7 (side-by-side lanes of stripe trains travelling in opposite directions)
  - random_start_success_share = 0.042 (share of random starts ending in opposite-travelling lanes (1 of 24))
  - noise_window_low = 0.28 (lower noise of window where opposite-travelling lanes form)
  - noise_window_high = 0.35 (upper noise of window where opposite-travelling lanes form)
  - lane_persistence_steps = 100000000.0 (time steps lanes persisted without disappearing)

## Reply R05

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Vicsek-model travelling bands (microphase separation) arranged as counter-propagating lanes; single-species analogue of the antiparallel-flocking (APF) band state", "notes": "Bands are textbook. Counter-propagating band lanes are documented only in models with two species, two preferred directions or anti-alignment. The active Potts and bidirectional-model papers say lanes are absent in the standard Vicsek model. The twist: lanes of transverse-band trains travelling in opposite directions, persisting for 1e8 steps in the single-species metric Vicsek model. Most likely a geometry-protected metastable state: 1 of 24 random starts, strip starts build it in, and periodic boundaries make lane reversal rare. Checks needed: box size and aspect ratio, starts not aligned with the box. Full texts of Chaté 2008, Solon 2015 and Caussin 2014 not read."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** I agree with the first pass. The parts of this behaviour that are documented are known. One specific twist in what the stripes do, which goes beyond the rule choice, is something I couldn't find described anywhere.

**What is documented**
- **Bands and microphase separation.** Kürsten & Ihle (PRL 2020) list the standard Vicsek phases, including a phase of polar ordered regularly arranged high density bands with surrounding low density regions. They also show the standard Vicsek model has a fourth phase for large system sizes: a polar ordered cross sea phase. Your parameters (density 1, speed 0.5, noise about 0.3 in the ±ηπ convention) sit in the band region.
- **Head-on collisions destroy opposite bands.** In constant-speed noisy systems, travelling bands are found to be unstable in a head-on collision. That fits with stripe trains only surviving if they slide past each other in separate lanes.
- **Opposite-travelling band states and lanes exist, but only in related models.**
  - In the two-species Vicsek model there is an antiparallel flocking state, in which the bands of species A and species B move in opposite directions. It needs two species that actively anti-align with each other.
  - In the active Potts model, a switch from transversal band motion to longitudinal lane formation is found, which is absent in the Vicsek model and the active Ising model.
  - In the discretised Vicsek model, the authors note the counter-propagating dense bands find it difficult to reverse its flow.
  - Other lane papers get two opposite streams only by adding something the standard model lacks: anti-alignment (arXiv 2412.10577), or particles assigned one of two preferred directions (each particle is assigned one of two opposite desired directions). That second paper says its lane phase is absent in the standard Vicsek model.

**The twist I couldn't find described**
In the standard single-species Vicsek model, the box splits into 7 side-by-side lanes. Each lane carries a train of ordinary transverse bands, and neighbouring lanes travel in opposite directions. Nothing in the rules favours two directions; it is a spontaneously split version of the two-species antiparallel state. It persisted for 10⁸ steps. This is not the active Potts lane state, where the dense stripe itself lies along the direction of motion, and not the cross-sea pattern.

Does it matter? Possibly, but it is probably a long-lived metastable state rather than a new phase:
- Only 1 of 24 random starts (4.2%) ended up in it.
- Starting from alternating opposite strips essentially builds the answer in.
- Periodic boundaries can lock in whole lanes: once a lane spans the box, flipping it needs a collective fluctuation, and the two-species papers call that kind of whole-band reversal a rare event. So box geometry is a likely explanation.

**Genuine or artefact?** It doesn't look like the obvious artefacts: not a random colouring, not a numerical blow-up, and not frozen, since the stripes keep travelling. But the fact that random starts rarely reach it and the strip start forces it makes geometry-induced lock-in plausible. To settle it you would need checks such as different box sizes and aspect ratios, starts that don't align lanes with the box, and the shear interface between lanes.

**Limits of the search:** I read abstracts and snippets only, not the full texts of Chaté 2008, Solon–Chaté–Tailleur 2015 or Caussin 2014. A remark on long-lived antiparallel band domains could be buried in one of them, so confidence is medium.

## Reply R13

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Cluster-lane pattern: parallel lanes of polar band/cluster trains moving in opposite directions (overall nematic order) in the band-forming regime of polar-alignment (Vicsek-class) active matter", "notes": "Same verdict class as the first pass, but it corrects the first pass's claim that the phenomenology was not found: opposite-travelling cluster-lane trains are documented as a stable solution in the kinetic theory of Vicsek-like polar alignment (Thüroff et al. 2014). Twist: seen in particle simulations of the standard metric Vicsek model as a rare (1/24 random starts), long-lived (1e8 steps) state, mostly reached from strip starts; no particle-simulation report found. Possible box-size locking of the 7 lanes and the parameter regime should be checked. Full texts of Chaté 2008, Solon 2015 and Kürsten-Ihle not read."}

**Reply text:**

**Verdict: KNOWN-VARIANT, and closer to KNOWN than the first pass thought.** The behaviour itself, side-by-side lanes of travelling polar clusters moving in opposite directions, is documented for the kinetic (Boltzmann-equation) version of Vicsek-style polar alignment. What I could not find is a report of it in particle simulations of the standard Vicsek model itself.

**The key source (the first pass missed it)**
Thüroff, Weber & Frey, *Phys. Rev. X* 4, 041030 (2014), studied self-propelled particles with only polar (head-to-tail aware) pairwise alignment, the BDG model. That model is inspired by Vicsek's agent-based model. Within the density-segregated parameter regime, they find a novel stable limit-cycle solution of the Boltzmann equation, which consists of parallel lanes of polar clusters moving in opposite directions, so as to render the overall symmetry of the system's ordered state nematic. The full text adds that each lane hosts a train of equally spaced clusters. That matches your observation closely: dense travelling stripes, arranged in parallel lanes going opposite ways, inside the regime where dense bands form.

**Related work that doesn't quite match**
- **Active Potts model:** there is a band-to-lane reorientation, but the authors say it is absent in the Vicsek model and the active Ising model. Their lanes are single dense lanes running along the direction of motion, not trains of clusters.
- **Discrete-direction flocks:** small counterpropagating droplets disrupt the flocking phase by inducing a persistent sandwich state, where the droplet-induced opposite-polarity lane remains embedded within the original flock. This is a similar long-lived opposite lane, but in a different model.
- **Cross-sea phase:** Kürsten & Ihle showed the standard Vicsek model has a fourth phase for large system sizes: a polar ordered cross sea phase. That pattern is band trains crossing at an angle, not opposite-moving parallel lanes. Still, it shows that big-box Vicsek runs can produce states beyond simple bands.
- **Lane models with built-in opposite motion:** recent lane models get lanes by assigning particles opposite target directions or by anti-alignment. That's a different mechanism, so they're not the source.

**What isn't documented, and whether it matters**
- A report of the cluster-lane state in particle simulations of the standard Vicsek model (metric neighbourhood, angular noise of the form used here), as opposed to in the kinetic equation.
- Some statistics: it appeared in only 1 of 24 random starts, sits in a noise window of 0.28–0.35, and lasted 10⁸ steps. In the kinetic theory it is a stable solution. In your system it looks like a rare, long-lived state that is mostly reached from the strip starting condition.

This plausibly matters a bit. Confirming in particle simulations, with real fluctuations, that a prediction from kinetic theory survives is a worthwhile check. But it isn't a new phenomenon.

**Genuine or artifact?** It looks like a genuine emergent pattern, not a frozen state or a numerical blow-up. There are three caveats:
- **Starting condition:** the strip start effectively plants the lanes.
- **Box size:** seven lanes may be locked in place by the periodic box. Changing the box width would test this.
- **Noise range:** check that 0.28–0.35 at density 1 and speed 0.5 is really inside the band regime rather than near the cross-sea or ordered-liquid boundary.

Limits of my search: I only saw abstracts and snippets. I didn't read the full texts of Chaté et al. 2008, Solon et al. 2015 or Kürsten & Ihle, so an earlier particle-simulation report could still exist.

## Reply R43

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Vicsek model band (liquid-gas coexistence / microphase-separated) regime; counter-propagating lane ('sandwich') state of metastable flocks", "notes": "Bands, cross-sea and the metastability of flocks to counter-propagating droplets are documented. Persistent opposite-polarity lanes are documented in discrete-symmetry (active Potts 'sandwich state') and two-species/bidirectional variants, but not found for the single-species continuous-angle standard VM. The twist is a multi-lane (7) antiparallel band-train state lasting 1e8 steps that is reached rarely (1 of 24) from random starts. Plausibly matters as a long-lived metastable state of the coexistence regime. Must rule out pinning by periodic box axes: test tilted or incommensurate boxes and system-size scaling. Full texts and supplementary videos of the large-system VM papers were not read."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium confidence.** I agree with the first pass's verdict but not fully with its reasoning. This is the standard Vicsek model in its band (coexistence) regime. The bands themselves are textbook. Side-by-side lanes moving in opposite directions are documented in closely related models, but I found no report of them in the standard single-species, continuous-angle Vicsek model.

**What is documented**
- **Standard model:** the Vicsek model has a coexistence region where finite-size liquid bands emerge in a gas environment and move along their transverse direction. Large systems also have a polar ordered cross sea phase (Kürsten & Ihle, PRL 2020). Your noise window of 0.28–0.35 fits this region: in a closely related discretized model, counter-propagating dense bands appear, and there is a cross-sea-to-band transition near η = 0.35.
- **Flocks are metastable against reversed regions:** a perturbation can break order in the Vicsek model by triggering a counter-propagating front that grows and eventually reverses the global polarity. In the published work such regions grow and take over; they don't settle into lanes.
- **Closest match, discrete-direction model:** in the four-state active Potts model (Chatterjee et al., PRE 2025), small counterpropagating droplets disrupt the flocking phase by inducing a persistent sandwich state, where the droplet-induced opposite-polarity lane remains embedded within the original flock. That is essentially your phenomenon, but in a model with only four directions, and the authors link it to spatial anisotropy.
- **Two-species model:** there is an APF (antiparallel flocking) state in which the bands of species A and species B move in opposite directions. That relies on built-in anti-alignment between species, which your system lacks.
- **Other lane papers:** lanes also appear when particles are assigned opposite preferred directions. The authors say this lane phase is absent in the standard Vicsek model. Lanes also appear with deliberate anti-alignment (arXiv 2412.10577). Neither matches your rules.

**The twist I could not find**
In the standard model, with no discrete directions, no species and no anti-alignment, you see:
- a multi-lane state (7 lanes) of band trains moving in opposite directions;
- persistence for 10⁸ steps;
- formation from random starts (1 in 24), not just from prepared strips;
- a noise window of 0.28–0.35.

Whether this matters: the literature expects continuous-angle flocks to reverse or rotate rather than lock into opposed lanes. If these lanes are truly stable rather than slow-decaying, that would be a real addition to the standard model's state diagram. Even then, it is a new metastable state of a known phase, not a new kind of behaviour.

**Real pattern or artifact?** It looks like a real emergent pattern: dense bands that keep travelling, not a frozen or diverged state. The main risk is the box. Periodic boundaries are known to lock band directions to the box axes; in one study the band direction is in principle random but aligns with the edges of our finite simulation box. Lanes running exactly along the axis could be held in place by the box itself. That would also fit the low 1-in-24 success rate from random starts and the use of prepared-strip starts. Three checks would settle it:
- tilted or non-commensurate box sizes;
- scaling with system size;
- whether lanes coarsen (merge over time) in larger boxes.

Of the searches I ran, only one returned nothing ("antiparallel lanes … stationary state … stripes"). Full texts and supplementary videos of the large-system Vicsek papers weren't read, so this phenomenon may already appear in a figure somewhere.

## Reply R46

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Vicsek model, liquid-gas coexistence with microphase-separated travelling bands (cross-sea / antiparallel-band relatives)", "notes": "Model and bands are firmly known; a stable multi-lane state of opposite-travelling band trains in single-species Vicsek was not found (closest: cross-sea phase, two-species Vicsek APF state, bidirectional Vicsek lane preprint). Only abstracts and citing snippets were read, so full-text mentions may exist. The 1/24 random-start rate and 1e8-step persistence are not addressed in the sources found. Check that rho=1, v0=0.5 puts the chosen noise window inside the band/coexistence regime."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium-low confidence.** The model and its dense travelling bands are well documented. I found no source describing the specific state you observed: several side-by-side lanes of opposite-travelling band trains, stable in a single-species Vicsek model.

**1. Model family**
- This is the standard Vicsek model of dry aligning active matter.
- The rule is average heading within radius r, plus uniform angular noise, then a step of fixed length.
- Each particle takes the average direction of particles within a circle of radius r around it, itself included, plus noise, then moves at constant speed.
- Your parameters are r = 1, v0 = 0.5, density 1 and angular noise η. They match the large-speed, large-system regime of the Chaté et al. 2008 simulations. I did not confirm those parameter values in the search results.

**2. What is known**
- **Bands in the coexistence region.** The flocking transition is a liquid-gas transition with micro-phase separation. The phase diagram shows a first-order transition with microphase separation.
- **Earliest comprehensive study.** Chaté, Ginelli, Grégoire and Raynaud, "Collective motion of self-propelled particles interacting without cohesion", Phys. Rev. E 77, 046113 (2008). It is a comprehensive study of Vicsek-style models in two and three dimensions. I did not read the band details in that abstract. I know of its band results from memory and from the citing literature.
- **Solon, Chaté and Tailleur 2015.** "From phase to microphase separation in flocking models: the essential role of nonequilibrium fluctuations", Phys. Rev. Lett. 114, 068101 (2015). The citing reference list confirms the title, authors and journal details. I only saw this paper through citing lists, not its abstract.
- **Review.** Chaté's lecture notes, Eur. Phys. J. Special Topics 225, 2099 (2016). Bands extend transversally to the direction of motion. On very long timescales they settle into a regular smectic arrangement of travelling bands. I inferred the author from the notes' content and did not verify it.
- **Cross-sea phase.** Kürsten and Ihle, "Dry Active Matter Exhibits a Self-Organized Cross Sea Phase", Phys. Rev. Lett. 125, 188003 (2020). The standard Vicsek model has a fourth phase for large systems, a polar-ordered cross sea phase.
  - It looks like a superposition of bands tilted in two directions.
  - This is the closest known large-system, multi-direction band pattern. It is polar-ordered, though, and its bands are tilted rather than antiparallel.
- **Cross-sea states in related Vicsek variants.** In some delay-Vicsek runs, cross-sea states formed in which bands move in different directions. That paper also noted uncertainty about whether cross-sea configurations or parallel-band states are truly stationary.
- **Stability of the polarised state.** At small η the ordered phase is very stable, so counter-propagating dense bands find it hard to reverse the flow. This is from the discretized Vicsek model.

**Close relatives, and what differs**
- **Two-species Vicsek model.** Chatterjee, Mangeat, Noh and Rieger, Phys. Rev. E 107, 024607 (2023). I took the year from the search result's age field, not from the abstract. It has a parallel-flocking state and an antiparallel-flocking state in which the bands of species A and B move in opposite directions.
  - Here the particles have species labels with anti-alignment between species.
  - Your system has one species with pure ferromagnetic alignment, so opposite-travelling lanes would be self-organised rather than imposed. That difference probably matters mechanistically.
- **Bidirectional Vicsek-type lane model.** A preprint on Zenodo reports a lane phase. In that phase polar order is small, nematic order is high, and the phase is absent from the standard Vicsek model. The preprint assigns each particle one of two desired directions, so it also differs from your system. Its claim that lanes are absent in standard Vicsek is mildly relevant to your finding. I could not verify its authorship or peer-review status.
- **Counter-propagating macrobands.** Taming polar active matter with moving substrates (New J. Phys., 2018). Bands comoving with a substrate can later reverse, started by one counter-propagating band that consumes the others. This needs an external moving potential and does not produce stable lanes.
- **Band collisions.** Constant-velocity noisy bands were found to be unstable in head-on collisions. That is a different model from yours. The result suggests collisions between opposite bands in plain Vicsek are normally destructive, which is consistent with your very low success rate of 1 in 24. A persistent lane state would need the bands to be spatially separated.

**3.

## Reply R49

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": false, "undocumented_phenomenology": "", "standard_name": "Vicsek model liquid-gas coexistence (travelling bands / microphase separation) with counter-propagating band metastability", "notes": "Bands, band arrays, cross-sea and counter-propagating band transients/reversals are documented for single-species Vicsek models. Long-lived side-by-side antiparallel lanes of band trains (7 lanes, ~1e8 steps, 1/24 random starts) were not found; stable lane/antiparallel band states in the literature need extra ingredients (two species with anti-alignment, preferred directions, anti-alignment), and one 2026 preprint states lanes are absent in the standard Vicsek model. Twist plausibly matters as a rare long-lived non-polar state. Needs checks: strip starts impose the state; box-size commensurability; verify 1e8-step runs at N~1.4M; test larger and other aspect-ratio boxes. Could become LITERATURE-SILENT if full-text and supplementary-video checks of the large-system Vicsek papers come up empty."}

**Reply text:**

**Verdict: KNOWN-VARIANT, low-to-medium confidence.** I agree with the first pass, but only just. The building blocks are well documented. The specific state is not: side-by-side lanes of band trains moving in opposite directions, lasting for a long time, in the ordinary single-species Vicsek model.

**What is documented (the "KNOWN" part)**
- **Bands, band arrays and the cross-sea phase.** The standard model is known in three different phases: a polar ordered homogeneous phase, also called the Toner-Tu phase, a phase of polar ordered regularly arranged high density bands with surrounding low density regions without polar order, and a homogeneous phase without polar order. Large systems add a fourth phase for large system sizes: a polar ordered cross sea phase (Kürsten & Ihle, PRL 125, 188003, 2020).
- **Opposite-moving dense structures in single-species flocks.** Ordered flocks can be destabilised: liquid polar flocks are metastable to the presence of a small obstacle [31] or to the nucleation of an opposite-phase droplet [32] and triggers counter-propagating dense clusters leading to the reversal of the liquid flow. In that work, for small η, the ordered phase is very stable, and thus, the counter-propagating dense bands find it difficult to reverse its flow. (This is from a paper on a discretised Vicsek model.)
- **Stable antiparallel band states exist, but only with an extra ingredient.**
  - In the two-species Vicsek model, an APF (antiparallel flocking) state in which the bands of species A and species B move in opposite directions is a known steady state. It depends on anti-alignment between the species.
  - Lanes also appear in a 2026 preprint where each particle has a fixed preferred direction. Its authors state that the phase diagram shows a lane phase and a critical density that are absent in the standard Vicsek model.
  - Particles that align *opposite* to their neighbours spontaneously split into two groups moving in opposite directions along parallel lanes.
  - None of these is the plain single-species model.

**What I could not find (the twist)**
I found no description of a long-lived state in the plain single-species, metric-alignment Vicsek model where several (7) lanes of band trains run side by side in opposite directions. That includes it arising from random starts (1 of 24) and persisting about 10⁸ steps. The 2026 preprint's claim that lanes are absent in the standard model is weak evidence that the community does not recognise such a state.

It plausibly matters. In the standard model, opposite-moving bands are normally treated as temporary: they collide, merge or reverse the flow. A lane arrangement in which the boundary between opposite flows somehow survives would be a new kind of long-lived, near-zero-net-polarisation state. Because of the 1/24 rate, it would be a rare long-lived state, not a phase.

I did not read full texts or supplementary videos. Large-system Vicsek papers often only remark that band configurations are multistable, so this state could already appear in a figure somewhere. That is why I don't call it LITERATURE-SILENT.

**Is it real or an artefact?**
The pattern looks genuine, but the claim needs checking:
- **Strip starts don't count as evidence.** Starting from alternating strips imposes the lanes. Only the 1/24 random-start result shows them forming on their own.
- **Check the 10⁸-step claim.** A 1200×1200 box at density 1 is about 1.44 million particles, so 10⁸ steps is a very large run. Confirm the run length, the random-number generator, and floating-point drift in the periodic wrapping.
- **Check the box size.** Lanes might be held in place by the box: if lane width times 7 divides the box side exactly, or in rectangular boxes, the arrangement may be fixed by size. Repeat at a few other sizes.
- **Check the number of lanes.** Make sure the 7 lanes don't shift towards a cross-sea or a single-direction band state at larger sizes.

## Reply R55

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Vicsek model liquid-gas coexistence: travelling dense bands, microphase separation, cross-sea phase", "notes": "Bands, microphase separation, band arrays and cross-sea states are established. Stable side-by-side antiparallel lanes of band trains (including from alternating-strip starts, 10^8-step persistence, 1/24 random-start rate) were not found as an explicit prior result, but are likely a special case of the known band-reversal and cross-sea phenomenology. Full texts were not read; a figure or SI video in the large-system Vicsek literature could already show it."}

**Reply text:**

**Verdict: KNOWN-VARIANT, with medium-low confidence.** The model and the dense travelling stripes are well known. I found no source that clearly describes side-by-side lanes of stripe trains moving in opposite directions as a stable state.

**1. Model family.** This is the standard Vicsek model of dry aligning active matter. It uses metric-neighbour alignment, angular noise η, speed v0 = 0.5, and a periodic box. The continuous-time variant with sine alignment and Gaussian noise is the Langevin version of the same class. The "dense travelling stripes" are the liquid-gas coexistence bands (microphase separation) of the Vicsek model.

**2. What the literature covers**
- **Bands and their pattern.** Chaté, Ginelli, Grégoire and Raynaud (2008) is the earliest clear source I found. It is a comprehensive study of Vicsek-style self-propelled particle models in two and three space dimensions (Phys. Rev. E 77, 046113; DOI 10.1103/PhysRevE.77.046113). I only saw the abstract and bibliographic data, so I did not verify the band details in the full text. Chaté's review, "Dry Aligning Dilute Active Matter" (Annu. Rev. Condens. Matter Phys. 11, 189–212, 2020), is the review source. Its snippet says traveling bands eventually form a regular smectic A pattern after possibly very long transients. It also notes that bands interact strongly when they become numerous. I saw only fragments of the review text.
- **Bands in different directions.** Kürsten and Ihle (2020), "Dry Active Matter Exhibits a Self-Organized Cross Sea Phase" (Phys. Rev. L 125, 188003), is the closest published relative. The Vicsek model has regularly arranged high density bands at intermediate noise, and for large systems there is a fourth phase. This looks like a superposition of bands, tilted in two directions. The pattern is self-organised and separated from the band phase by discontinuous transitions.
- **Parallel and cross-sea bands at similar parameters.** A study of a discretised Vicsek model (New J. Phys. 2024, "Consequence of anisotropy on flocking: the discretized Vicsek model") uses v0 = 0.5 and noise around 0.3–0.4. It reports multiple bands moving parallelly or in a cross-sea pattern through the gas phase, and it also discusses counter-propagating bands and band reversal.
- **Large-system Vicsek simulations.** The delay-Vicsek paper (New J. Phys. 2025) says that in some simulation runs, this led to the formation of so-called cross-see states, where bands move in different directions.

**3. What differs, and whether it matters.**
- I found no source that reports a lane structure of antiparallel band trains in the plain Vicsek model. The reported cases are cross-sea patterns with tilted bands, or parallel bands moving in one direction.
- The "alternating opposite-heading strips" initial condition and the 10^8-step persistence do not appear in anything I found.
- The 1-in-24 random-start success rate and the 0.28–0.35 window are specific numbers I could not match to a source.
- The lane picture is plausibly a special case of the known band-reversal and cross-sea phenomenology in elongated or large periodic boxes. It may also be a long-lived metastable state, in which case the 10^8 steps would simply show long lifetimes. That is my inference, not a reported result.

Two sources I excluded as different systems. A 2025–26 Zenodo preprint on a "bidirectional Vicsek-type model" has particles with built-in opposite desired directions. A paper on antialigning Vicsek-type rules has antiparallel alignment. Neither is the standard Vicsek rule.

**4. Confidence.** Medium-low. The family and the bands are certain. The lane claim rests on abstracts and snippets, since I could not read full texts. It also depends on a negative search, so a specific figure or supplementary video could already show antiparallel lanes.

**
