# Stage 3 (b) review inputs: D_agents_fields_oscillators (16 phenomena)

## census:g1:agents-1

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. Rules, applied in order each step: (1) each agent turns toward the centre of any type agents within distance 1.0.

OBSERVED (features that passed the simulation checks):
  a_density_cv, a_step, a_nn_ratio, a_local_align

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_density_cv: median 2.2094 with interaction, 0.8089 with interaction knocked out (difference 1.4005, robust spread 0.0695); difference at 4x area 0.7113, at 3x run length 1.4494
  - a_step: median 0.2751 with interaction, 6.6006 with interaction knocked out (difference -6.3256, robust spread 0.1375); difference at 4x area -19.4756, at 3x run length -7.0327
  - a_nn_ratio (clumping (below 1 = clustered)): median 0.1752 with interaction, 0.99 with interaction knocked out (difference -0.8148, robust spread 0.0233); difference at 4x area -0.8191, at 3x run length -0.8327
  - a_local_align (local alignment): median -0.2056 with interaction, -0.0393 with interaction knocked out (difference -0.1663, robust spread 0.007); difference at 4x area -0.1679, at 3x run length -0.1586

### Stage 3 reply: verdict KNOWN, standard name: Local Attraction Model (Strömbom 2011) — attraction-only self-propelled particles, swarm/aggregation regime

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Rule is the Local Attraction Model (steer toward local centre of mass within R, constant speed, noise). Strong clustering, near-stationary aggregates (low net displacement) and low or negative local alignment match the documented 'swarm' regime at strong attraction. Negative nearest-neighbour alignment is not stated explicitly in the sources read, but follows geometrically from agents crossing through the local centre. Overturns 'variant' to 'known' because there is no meaningful twist in the phenomenology. Suggested check: measure rotation order to rule out mills.

Citations:
- Strömbom D. (2011) Collective motion from local attraction. J Theor Biol 283:145-151, doi:10.1016/j.jtbi.2011.05.019
- Strömbom D., Hassan T., Greis W.H., Antia A. (2019) Asynchrony induces polarization in attraction-based models of collective motion. R Soc Open Sci 6:190381
- Strömbom D., Antia A. Anticipation induces polarized collective motion in attraction based models, arXiv:1710.05692
- Strömbom D., Siljestam M., Park J., Sumpter D.J.T. (2015) The shape and dynamics of local attraction. Eur Phys J ST 224:3311-3323

Full reply text:

**Verdict: KNOWN.** This is Strömbom's Local Attraction Model (LAM) running in its "swarm" regime. I'm overturning the first pass. It's not really a variant: the rule is essentially the same as LAM, and every measured feature matches a documented behaviour of it.

**Why it's the same system.** In LAM, N particles move at constant speed in two dimensions and interact via local attraction only. On every time step, each particle calculates the position of the local centre of mass (LCM) of all particles within a distance of R from it (its neighbours). The new heading mixes the direction to the LCM, the current heading and a noise term (a linear combination of the normalized direction toward the LCM, its normalized current heading and a normalized uniform noise vector). Your rule is that: one type, constant speed 1.0, heading noise, and a turn toward the local centre within radius 1.0. Strömbom introduced it in J. Theor. Biol. 283:145–151 (2011). Later work calls it the simplest model known to produce all three of these groups (polarized groups, mills and swarms).

**How the observed features map onto known behaviour.**
- **Which regime you're in depends on turn strength.** Depending on the relative strength of local attraction (c) groups of different type will be produced by the model. In the synchronous version, if c < 0.2 no group will form. If 0.2 < c < 1.7 mills will form, and if c > 1.7 swarms will form.
- **Strong clumping and high density variation (a_nn_ratio ≈ 0.18, a_density_cv ≈ 2.2).** These match aggregated groups rather than the no-group state.
- **Very low net movement (a_step ≈ 0.28 vs 6.6 in the knock-out).** This matches the strong-attraction case, where for c ∈ [1.7, 2] relatively stationary swarms form (σ small and α small). A rule that turns straight toward the local centre is effectively at high c, so the swarm regime is what you'd expect. The spec doesn't give an exact turn rate, so I'm inferring this.
- **Negative local alignment (≈ −0.21).** Low swarm polarization is documented. A *negative* nearest-neighbour heading correlation isn't stated in what I read, but it follows directly from the swarm shape: agents overshoot the local centre and pass through it, so close neighbours tend to face opposite ways. I don't think this is a phenomenological twist that matters. Mills would show rotational order instead; you didn't measure rotation, which would be a useful check.

**Genuine or artifact?** Genuine. It's a real aggregation transition, not random colouring or a diverged field. The gaps from the knock-out are many robust spreads wide and hold up at 4× area and 3× run length. The one honest caveat is that these swarms are close to a stationary, kinetically trapped state, so "low displacement" is part of the known phenomenon rather than something new.

**Searches run:** Strömbom's local attraction model (2011); LAM regimes (swarm, mill, no group) versus attraction strength c. I checked sources including the Royal Society Open Science 2019 paper and the arXiv anticipation paper.

---

## census:g1:agents-2

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. Rules, applied in order each step: (1) each agent turns away from the centre of any type agents within distance 3.0.

OBSERVED (features that passed the simulation checks):
  a_step

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_step: median 0.784 with interaction, 6.5048 with interaction knocked out (difference -5.7208, robust spread 0.1458); difference at 4x area -19.127, at 3x run length -6.4383

### Stage 3 reply: verdict KNOWN, standard name: Repulsion-only (zone-of-repulsion) self-propelled particle model, i.e. a Couzin/Vicsek-type model with only the repulsion zone; disordered homogeneous active gas with reduced persistence

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: The rule is the Couzin zone-of-repulsion rule on its own. Pure repulsion is documented to give a disordered homogeneous state. Flocking from turning away requires positional repulsion, which is absent here. a_step is not precisely defined; it is probably a persistence or displacement measure. Its sharp drop compared with the knock-out matches neighbour-induced turning wiping out persistence at about 28 neighbours per interaction disc. That is a real but trivial direct effect of the interaction, not an emergent pattern. 'Trivial' here means not an emergent pattern; it is not a numerical artifact.

Citations:
- https://arxiv.org/pdf/1301.6294
- https://www.researchgate.net/publication/255989028_Self-propelled_particles_with_selective_attraction-repulsion_interaction_From_microscopic_dynamics_to_coarse-grained_theories
- https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.031008

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass's "variant" label. Nothing new happens in the behaviour here, and the measured effect is a trivial one rather than an emergent pattern.

**The system is a standard one.** Agents that only steer away from the centre of neighbours within a radius are using the "zone of repulsion" rule from the Couzin and Vicsek-type swarm models. In that rule, if other particles are present within the zone of repulsion, the individual attempts to keep its personal space and to avoid collision by moving away from the neighbours. This rule has been proposed by Couzin et al. Those models even use the same constant propulsion speed v0=1 and time step Δt=1. Your system is that model with the alignment and attraction zones removed.

**What pure repulsion is known to do.** In the closely related selective attraction–repulsion model, pure repulsion produces a disordered, homogeneous state, not bands or flocks. One mechanism by which turning away can produce flocking is ruled out here. That flocking happens when particles also push each other apart physically: the emergence of polar order from turn-away interactions requires particle repulsion, because repulsion pushes them apart before they can turn away too much, thus producing effective alignment. Your agents only turn; nothing displaces them. So the expected outcome is a disordered gas, and nothing in the evidence suggests otherwise.

**Why the a_step feature is trivial.** a_step drops from about 6.5 with the interaction knocked out to about 0.78 with it, and the gap grows at 4x area. I don't know exactly how the program defines a_step. I'm assuming it measures how far agents travel persistently, or how far their position or heading carries over across steps. The scaling fits that reading: free agents with 0.3 rad of noise keep their direction for a long time, while interacting agents keep getting turned. The density makes this almost certain. Alone, there are 400 agents on 20×20, so about 1 per unit area. A radius-3 neighbourhood then holds roughly 28 neighbours, so each agent is re-steered every step and mostly caged. Losing persistence and diffusing less because of neighbour-induced turning is the expected result of any repulsive or collision-avoiding interaction. It's a real effect, not a numerical artifact or a frozen state, but no pattern has formed. Without knowing what a_step measures, I can't completely rule out some other effect, though nothing points to one.

I ran a single targeted search, because the system matched the standard zone-of-repulsion model right away.

---

## census:g1:agents-3

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 3.14159 rad per step. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 1.0.

OBSERVED (features that passed the simulation checks):
  a_density_cv, a_local_align

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_density_cv: median 1.0545 with interaction, 0.7899 with interaction knocked out (difference 0.2646, robust spread 0.0369); difference at 4x area 0.297, at 3x run length 0.2398
  - a_local_align (local alignment): median 0.1717 with interaction, -0.0462 with interaction knocked out (difference 0.218, robust spread 0.0138); difference at 4x area 0.1949, at 3x run length 0.2155

### Stage 3 reply: verdict KNOWN, standard name: Vicsek model (1995): short-range orientational order and alignment-induced clustering/density fluctuations near the order-disorder transition

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: The rule set is the standard Vicsek model (unit radius, constant speed, angular noise, periodic box), not a variant, so 'VARIANT' is a small overturn to 'KNOWN'. A uniform noise width of pi (eta=pi) is consistent with the nonzero local alignment; full noise over [-pi,pi] would erase alignment, so check the noise convention in the code. With density about 1 and eta=pi the system is near the transition: weak local alignment (0.17) plus higher density CV matches documented short-range correlations and clustering/bands. The effect is robust against the knock-out, so not an artifact, but it is weak order, not strong flocking.

Citations:
- Vicsek, Czirok, Ben-Jacob, Cohen, Shochet, 'Novel Type of Phase Transition in a System of Self-Driven Particles', PRL 75, 1226 (1995)
- Chatterjee, 'The Flocking Transition: A Review of the Vicsek Model' (UCSD PHYS 563 essay, 2017)
- Active-learning mapping of the Vicsek model phase diagram, arXiv:2604.28167

Full reply text:

**Verdict: KNOWN.** This is the original Vicsek model, not a variant of it. Calling it a "variant" undersells how closely it matches. Moving it from VARIANT to KNOWN is a small correction.

**Same system.** The rules match the 1995 model point for point. Vicsek et al. used particles that are driven with a constant velocity and at each time step assume the average direction of motion of the particles in their neighborhood with some random perturbation added. Standard descriptions say each particle reorients itself along the average direction of propulsion of particles in a unit radius around it, with some external noise preventing perfect alignment, within a square domain with periodic boundary conditions. Your setup has the same ingredients: unit radius, constant speed, angular noise, wrap-around box. Only the parameters differ: density 1 (if the 20×20 box applies), speed 1.0, noise π.

**How to read the noise.** In the standard convention the noise is drawn uniformly from the interval [-η/2, η/2], with η going up to 2π. If your "3.14159 rad" meant a uniform draw from [−π, π], each new heading would be completely random and you couldn't see any alignment. Since you did see positive local alignment, the noise is probably the standard η = π (draws from [−π/2, π/2]). I'd confirm this in the code. For scale, the original paper reports a critical noise of about η_c(∞) = 2.9 ± 0.05 for a density of ρ = 0.4, and increasing density shifts the crossovers to larger noise because each particle has more neighbors within the interaction radius. So at density about 1 and η = π, the system is near the transition or in the coexistence region, where weak order is expected.

**Both features are documented:**
- **Local alignment (0.17 vs −0.05):** below the global transition, Vicsek systems keep short-range order. The pair orientation correlation exhibits exponential decay in the disordered phase, indicative of short-range order.
- **Density CV going up (1.05 vs 0.79):** clustering caused by alignment is a central Vicsek result. The original paper set out to study clustering, transport and phase transition in non-equilibrium systems. Near the transition, high density ordered bands traveling in a disordered background appear, and recent work describes dense, locally ordered bands coexisting with a dilute, weakly ordered background.

**Genuine or artifact?** Genuine, though weak. The no-interaction baseline (density CV 0.79, alignment about 0) is plain random clumping of independent walkers. The interacting run beats it by many spreads and holds up at 4× area and 3× run length. That rules out a frozen state, a random colouring or a numerical blow-up. But the size of the effect (alignment 0.17) points to disordered or near-transition short-range correlations, not strong flocking. There's no new behaviour here: the known Vicsek behaviour, at a parameter point inside the explored phase diagram.

---

## census:g1:agents-8

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0.

OBSERVED (features that passed the simulation checks):
  a_step

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_step: median 0.7902 with interaction, 19.7547 with interaction knocked out (difference -18.9645, robust spread 0.4109); difference at 4x area -19.0006, at 3x run length -22.2537

### Stage 3 reply: verdict KNOWN, standard name: Chemotactic gradient ascent / trapping of active (persistent random-walk) particles at local maxima of a static, passively diffusing chemoattractant landscape

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: No emission and no agent-agent sensing, so the field is inert and the agents are independent. The knock-out removes the only rule, leaving persistent random walkers (persistence ~20 steps, which fits net displacement ~19.75). With gradient-following, agents localize at residual maxima of the diffusing initial field (~0.79). The suppressed displacement is a trivial effect of following a gradient, not an emergent collective pattern. It probably persists late in the run because the steering responds to gradient direction regardless of how weak the gradient is, giving a near-frozen trapped state. Calling it a 'variant' of a collective chemotaxis model overstates it.

Citations:
- https://arxiv.org/html/2511.18745
- https://arxiv.org/pdf/1202.2442

Full reply text:

**Verdict: KNOWN. This is a trivial, artifact-like outcome. I'm overturning the first pass, which called it a variant of a known model.**

**What the system actually is.** There is one rule: each agent turns up the gradient of field 0. No rule makes agents emit into the field, and agents don't sense each other. So the field never reacts to the agents. It is a passive landscape that starts with some initial pattern and smooths out by diffusion. That rules out Keller–Segel aggregation, chemotactic collapse and self-trapping, which all need agents to emit. Here the agents are just doing noisy gradient-climbing on a fixed landscape that slowly flattens.

**What the knock-out measures.** Removing the "interaction" also removes the only rule, so the knock-out is plain persistent random walkers. With 0.3 rad of heading noise per step, a walker keeps its direction for roughly 20 steps. Over a short time window it moves about 20 units in net distance, which fits the knock-out value of 19.75. With gradient-following on, agents settle near local maxima of the field and wobble in place, which fits the value of 0.79. I'm assuming a_step is net displacement over a window, since a per-step speed of 1.0 can't give 19.75. The large gap between the two values is just "following a gradient vs. not following one", not something that emerges from interaction.

**Literature.** It is well documented that chemotactic agents get trapped at local maxima of a fixed chemoattractant landscape. A 2025 paper on chemotaxis and stochastic gradient ascent describes how local maxima or competing gradients may transiently retain, redirect, or effectively "trap" chemotactic cells away from a desired region. It also notes that additional local maxima present a more challenging environment, leading to more pronounced trapping of trajectories in secondary peaks. The active Brownian particle literature makes the same point: agents stay pinned at field maxima unless something weakens their response to the gradient, in which case the particles are able to leave the maxima of the field distribution.

**Why I call it trivial.**
- The field is inert, with no feedback from the agents.
- Nothing in the rules reads other agents, so the "interaction" label is misleading.
- The trapped state persists late in the run probably because the steering responds to the gradient's direction, not its size. Even a nearly flat field still pins agents to its leftover peaks. That looks like a frozen, nearly absorbing state, not a self-sustaining pattern.

**Twists I could not attribute.** None. Staying the same at 4× area and 3× run length is what you'd expect from agents acting independently on a slowly decaying landscape. One caveat: I only ran a single search, so the verdict rests on that plus first-principles reasoning about the rules.

---

## census:g1:network-phase-2

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) each node is a phase oscillator pulled toward its neighbours' phases with strength 1.0; (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (every step).

OBSERVED (features that passed the simulation checks):
  p_r_std, p_freq_spread

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - p_r_std (phase order fluctuates): median 0.1157 with interaction, 0.0289 with interaction knocked out (difference 0.0868, robust spread 0.0142); difference at 4x area 0.0864, at 3x run length 0.1611
  - p_freq_spread (spread of oscillator natural frequencies): median 0.0552 with interaction, 0.4987 with interaction knocked out (difference -0.4435, robust spread 0.0231); difference at 4x area -0.4565, at 3x run length -1.2561

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: p_r_std 0.0058; p_freq_spread 0.0643
  - rules 2 only: p_r_std None; p_freq_spread None
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Kuramoto model on a homophily-fragmented (modular/disconnected) network: intra-community frequency locking with beating/oscillating global order parameter (cluster/modular synchronization)

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Static types + always-on cross-type rewiring drive the graph to a frozen absorbing topology of ~3 disconnected same-type components; rule 2 then becomes inert. Kuramoto locks each component at its own mean frequency (low effective-frequency spread), and the global order parameter beats between the incommensurate cluster phasors (high r std). No ongoing feedback between the rules; growth of the effect with run length points to a frozen-state/measurement artifact rather than a novel pattern. Overturns 'variant' toward plain KNOWN with trivial mechanism.

Citations:
- Community Structure and Multi-Modal Oscillations in Complex Networks, PLoS ONE 2013 (10.1371/journal.pone.0075569)
- Rodrigues et al., The Kuramoto model in complex networks, arXiv:1511.07139
- Cluster Synchronization of Kuramoto Oscillators and the Method of Averaging, arXiv:2204.02627
- Papadopoulos et al., Chaos 27, 073115 (2017), arXiv:1706.08070 (adaptive Kuramoto rewiring -> modular organization, local not global synchrony)

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass's "variant" label. The behaviour is a standard, well-documented one, and what is happening underneath is largely trivial.

**What the system actually does.** The node types never change, and rule 2 replaces every link between different types, every step. So the network quickly ends up in a frozen state: roughly three separate groups, one per type, with no links between them. Once that happens, the second rule stops doing anything, and what's left is the ordinary Kuramoto oscillator model running on three unconnected random graphs. Each group syncs internally at its own average frequency. There's nothing to pull the three groups together, so their combined phases drift past each other ("beat"). That makes the overall order parameter rise and fall, which is p_r_std ≈ 0.116 versus about 0.006 for a single fully synced graph (rule 1 only).

**The two features:**
- **p_r_std** (phase order fluctuates): an oscillating global order parameter on modular or clustered networks is well documented. One study shows global order parameter time series with "multi-modal dynamics" on clustered random networks. It found that "the oscillating behaviour begins to break down as the hierarchical structure is decreased". A review notes "the sparser the connection between nodes in different communities, the harder for these nodes to lock in a common phase". Your case is the extreme version, with zero links between groups. The cluster synchronization literature also tracks separate per-cluster and network-wide order parameters for three clusters at different frequencies (arXiv 2204.02627).
- **p_freq_spread**: despite the name, this falls from about 0.50 to about 0.055 when coupling is on. That means it measures the spread of *effective* frequencies, not natural ones, since natural frequencies don't change. Its value is just ordinary frequency locking within each group, which is textbook Kuramoto. The small remaining spread (0.055) is close to the rule-1-only value (0.064), and is about what the gap between three group-average frequencies would give.
- **Adaptive-network context:** Kuramoto models with adaptive rewiring are known to produce modular structure. Reinforcing links between dynamically similar nodes has been shown to lead to the emergence of modular organization in Kuramoto networks, giving enhanced local rather than global synchrony. Here the modules come from fixed types rather than phase similarity. That changes the rule, not the behaviour you observe.

**Trivial or artifact?** Partly trivial. The dynamics are real, not a numerical failure, but the result follows directly from (a) homophilic rewiring with fixed types locking the network into separate same-type groups, plus (b) standard Kuramoto locking within each group. The two rules don't keep feeding back on each other: rule 2 stops mattering once the groups have separated. Two further signs of trivial scaling: the fluctuation difference roughly doubles at 3× run length, and the frequency-spread difference grows to −1.26. Both point to measurement-window effects in an already frozen state, not a stable emergent pattern. One caveat: I didn't find a paper using exactly this setup (Kuramoto plus fixed-type homophilic rewiring). By your own criteria, that missing rule combination doesn't count as a discovery.

---

## census:novelty:agents-16

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 1 type (random start). Agents have 4 types (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. 2 field(s) with diffusion [0.02, 0.02]. Rules, applied in order each step: (1) each agent turns up the gradient of field 1; (2) each agent turns left on an even-type cell and right on an odd-type cell.

OBSERVED (features that passed the simulation checks):
  a_density_cv

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_density_cv: median 1.8433 with interaction, 0.8 with interaction knocked out (difference 1.0433, robust spread 0.0529); difference at 4x area 1.0654, at 3x run length 1.0769

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: a_density_cv 6.4938
  - rules 2 only: a_density_cv 0.799
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Chiral active particles (circle swimmers) doing gradient ascent / trapping in a quenched random field landscape; chirality reduces accumulation

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Cells have only one type, so rule 2 is a constant left turn (chirality); alone it gives the same CV as the knock-out (0.799 vs 0.8). No emission rule is listed, so field 1 is probably decaying initial noise that acts as a fixed landscape; agents pile up at its peaks (CV 6.49 with rule 1 alone). The full program's CV of 1.84 is chirality partly breaking that trapping, a documented effect. The 'interaction' is agent-reads-environment, not agent-reads-agent. The clustering follows the field's initial noise rather than self-organisation, so it is close to a trivial outcome. Overturns 'variant' toward KNOWN/trivial. If agents actually emit field 1, reclassify as chiral Keller-Segel aggregation (still known).

Citations:
- https://www.researchgate.net/publication/254498325_Active_Brownian_particles_in_two-dimensional_traps
- https://www.researchgate.net/publication/331170019_Active_chiral_particles_under_confinement_Surface_currents_and_bulk_accumulation_phenomena
- https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/unsteady-chiral-swimmer-and-its-response-to-a-chemical-gradient/39FD23E37ED399BFD4324811DD5A1785

Full reply text:

**Verdict: KNOWN.** I agree that this is a known behaviour, but not that it is a "variant". I'd go further than the first pass and call the measured feature mostly trivial, not a new emergent pattern.

**What the rules actually do**
- **Rule 2 is just a constant left turn.** Cells have only one type, so every cell is type 0, which counts as "even". Every agent therefore turns left every step, so the agents are basically **chiral active particles (circle swimmers)** with heading noise. With only rule 2 the density CV is 0.799, the same as the knock-out (0.8). On its own it does nothing to how agents spread out.
- **Rule 1 means agents steer up a fixed background field.** No rule makes agents emit field 1 or consume it. So the field is probably just its random starting values, slowly smoothing out (diffusion 0.02). Agents steering up its gradient pile up at its local peaks. That explains the very high CV with rule 1 alone (6.49). It's single-agent trapping in a random landscape, not agents interacting with each other. The knock-out "interaction" here is agent-reads-environment, not agent-reads-agent.
- **Full program (CV 1.84):** the constant turning makes agents circle and orbit rather than settle at the peaks. So they cluster less than with rule 1 alone, but more than random. That 1.84 sits between the two sub-programs is what you'd expect from two separate effects competing. It isn't a new emergent mechanism.

**Literature**
- Non-interacting active particles in an external potential or landscape are known to pile up in a non-random way. Simulations of microswimmers exploring a heterogeneous, random landscape suggest a rich phenomenology due to the interference of the landscape with the persistence properties of the trajectories.
- Chiral turning is known to suppress clustering. In circle swimmers, active rotation can induce localized motion and, thus, generally oppose motility-induced clustering and phase separation. Also, by modulating rotational dynamics, it is possible to control the number fluctuation in clusters and retard or promote particles' aggregation.
- Chiral swimmers in chemical gradients have been studied directly. One study reports a complete state diagram representing the successful, unsuccessful and orbital states of a chiral swimmer responding to a chemical gradient.

**Is it an artifact?** Partly. The clustering comes from the field's initial random noise, which is close to "an inert field acting as a fixed landscape." Agents don't build or keep up the pattern themselves. The effect does hold up across seeds, at 4x area and at 3x run length. But it doesn't come from agents interacting with each other. The cell lattice adds nothing except a constant turn. The four agent types do nothing in these rules.

One check would settle the artifact question: confirm the field has no emission or decay. If agents do in fact emit field 1, this becomes Keller–Segel-type chemotactic aggregation with chirality, which is also documented (chiral swimmers in chemical gradients, as above).

---

## census:novelty:agents-17

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns away from the centre of any type agents within distance 2.0.

OBSERVED (features that passed the simulation checks):
  a_step, a_density_cv, a_nn_ratio, a_local_align

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_step: median 2.342 with interaction, 17.4443 with interaction knocked out (difference -15.1023, robust spread 0.3957); difference at 4x area -15.0047, at 3x run length -23.2641
  - a_density_cv: median 1.1785 with interaction, 0.792 with interaction knocked out (difference 0.3865, robust spread 0.0394); difference at 4x area 0.6708, at 3x run length 0.4286
  - a_nn_ratio (clumping (below 1 = clustered)): median 0.8208 with interaction, 0.9915 with interaction knocked out (difference -0.1707, robust spread 0.0362); difference at 4x area -0.203, at 3x run length -0.2331
  - a_local_align (local alignment): median -0.1074 with interaction, -0.0265 with interaction knocked out (difference -0.0809, robust spread 0.0197); difference at 4x area -0.0703, at 3x run length -0.0788

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: a_step 0.8494; a_density_cv 11.1953; a_nn_ratio 0.0416; a_local_align 0.0511
  - rules 2 only: a_step 4.6543; a_density_cv 0.5879; a_nn_ratio 1.2384; a_local_align -0.0621
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Autochemotactic (Keller-Segel) aggregation of self-propelled particles with short-range repulsion limiting chemotactic collapse (finite/dynamic clustering)

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Emission into field 0 is not stated in the rule but is implied by the rules-1-only near-total collapse (density_cv ~11, nn_ratio ~0.04), i.e. the documented chemotactic collapse. Adding short-range turn-away repulsion stops this at moderate clustering, as documented for chemotactic walkers with steric repulsion. Reduced step length follows from agents being trapped in clusters. Weak negative local alignment is a trivial side effect of turn-away repulsion (rule 2 alone gives -0.06) and is not new order. Robust across seeds, 4x area and 3x run length, so not an artifact. Only abstracts/snippets were checked, not full texts; confidence medium because the emission term is inferred.

Citations:
- Taktikos, Zaburdaev & Stark, Collective dynamics of model microorganisms with chemotactic signaling, Phys. Rev. E 85, 051901 (2012) https://pubmed.ncbi.nlm.nih.gov/23004782/
- Pohl & Stark, Dynamic clustering and chemotactic collapse of self-phoretic active particles, arXiv:1403.4063 https://arxiv.org/html/1403.4063
- Stark, Artificial Chemotaxis of Self-Phoretic Active Colloids: Collective Behavior, Acc. Chem. Res. 51, 2681 (2018) https://pubmed.ncbi.nlm.nih.gov/30346724/
- Liebchen et al., Synthetic Chemotaxis and Collective Behavior in Active Matter, arXiv:1810.01117
- Phase separation of chemokinetic active particles, Commun. Phys. (2025) https://www.nature.com/articles/s42005-025-02211-0

Full reply text:

**Verdict: KNOWN.** I agree with the first pass that this is the known model. I'd go slightly further: the behaviour itself is already documented, so it isn't really a "variant". This is self-propelled agents that follow a chemical they produce themselves (autochemotaxis, the Keller–Segel aggregation instability), with clumping held back by short-range repulsion.

**What the rules actually are.** The plain-words rule never says agents emit field 0, but they must. With only rule 1, the agents collapse almost completely (density CV about 11.2, nearest-neighbour ratio about 0.04). Without a source, gradient-following on a diffusing field couldn't do that. So the knock-out keeps emission and only removes the gradient-following and the repulsion.

**How each feature maps to the literature:**
- **Rule 1 alone gives total collapse.** That is the documented attract-only case. In a particle model of chemically interacting swimmers, "in case of mere attraction a chemotactic collaps occurs directly from the gas-like state before pronounced clusters are able to form", and this was supported by mapping the Langevin dynamics onto the Keller-Segel model for bacterial chemotaxis. A review states it directly: when particles release the chemical themselves, "it is well known that chemoattraction leads to complete phase separation".
- **Adding repulsion gives moderate, limited clustering.** In the full program, density CV is 1.18 and nearest-neighbour ratio is 0.82. That is collapse kept in check by repulsion, which is also documented. Taktikos, Zaburdaev & Stark (PRE 85, 051901, 2012) studied walkers whose velocity is modified when they overlap, via a linear force-velocity relation and a harmonic repulsion force. They described cluster formation in many-walker systems. Later work summarises them as showing that steric interactions can limit aggregation, as indeed they must. The Stark-group state diagrams also cover dynamic clustering between a gas and a collapsed state (dynamic clustering states 1 and 2 ... to different types of collapsed states).
- **Lower step length (a_step).** Agents are trapped inside clusters with their headings scrambled, so they travel less far. This follows directly from the aggregation.
- **Slightly negative local alignment (about −0.11).** This is a predictable side effect of the turn-away rule: neighbours that just repelled each other point in opposite directions. Rule 2 on its own already gives −0.06. It isn't a new kind of order. It's the only feature that is close to its single-rule value, so I don't count it as a twist.

**Genuine pattern or artifact?** It's a genuine emergent pattern, not an artifact. The full program is clearly separated from both single-rule versions and from the knock-out. The clustering signal stays when the area is 4x larger (the density-CV difference actually grows) and when the run is 3x longer. It isn't the frozen absorbing state, because that is the rule-1-only collapse, and repulsion keeps the full system out of it. Nothing suggests numerical divergence. But it's a textbook outcome, not a new one.

**Limits of this check.** I relied on abstracts and search snippets, not the full papers. I didn't check that the 2012 paper's cluster sizes or numbers match this run. The "emission is present" conclusion is inferred from the sub-program results, not stated in the rule.

---

## census:novelty:agents-8

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns to the average heading of agents within distance 2.0.

OBSERVED (features that passed the simulation checks):
  a_density_cv, a_step, a_local_align

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_density_cv: median 10.9707 with interaction, 0.8149 with interaction knocked out (difference 10.1559, robust spread 0.2616); difference at 4x area 6.3385, at 3x run length 11.6614
  - a_step: median 0.7435 with interaction, 17.4773 with interaction knocked out (difference -16.7339, robust spread 0.3644); difference at 4x area -16.7341, at 3x run length -26.341
  - a_local_align (local alignment): median 0.8728 with interaction, -0.0362 with interaction knocked out (difference 0.9091, robust spread 0.0182); difference at 4x area 0.92, at 3x run length 0.929

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: a_density_cv 11.2833; a_step 0.8847; a_local_align 0.0604
  - rules 2 only: a_density_cv 1.6327; a_step 18.3601; a_local_align 0.3412
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Autochemotactic (Keller-Segel-type) collapse of self-propelled particles with Vicsek alignment

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Rule 1 alone reproduces the clumping and low step; the knock-out implies agents emit field 0 (autochemotaxis), although no emission rule is listed. High local alignment inside the collapsed clumps follows almost automatically from Vicsek averaging once agents sit within radius 2; low step reflects agents trapped by their own chemical peak. Genuine emergent collapse, but close to an absorbing state; the alignment score adds little beyond that geometry. Overturns the first pass mildly (VARIANT -> KNOWN): no phenomenological twist beyond documented collapse plus alignment. Verified from search excerpts only.

Citations:
- Pohl & Stark, Dynamic clustering and chemotactic collapse of self-phoretic active particles, arXiv:1403.4063 / PRL 2014
- Saha, Golestanian & Ramaswamy, Clusters, asters and collective oscillations in chemotactic colloids, arXiv:1309.4947
- Boltzmann-Ginzburg-Landau theory for autochemotaxis of active rod-like particles, arXiv:2506.18263
- Stark, Artificial Chemotaxis of Self-Phoretic Active Colloids: Collective Behavior, PubMed 30346724

Full reply text:

**Verdict: KNOWN.** This is chemotactic collapse of self-propelled particles (Keller–Segel type) with Vicsek alignment layered on top. Both pieces, and the combination, are in the literature. I'm moving the first pass from "variant" to "known" because I don't see a new twist in the behaviour.

**What the sub-programs show**
- **Clumping comes from following the chemical alone.** With only rule 1, clumping (a_density_cv ≈ 11.3) and how far agents move per step (a_step ≈ 0.88) already match the full program (10.97 and 0.74). The knock-out drops rule 1, so reading the field counts as an interaction. That only makes sense if agents release the chemical themselves (self-chemotaxis). Your rule text doesn't actually list an emission rule; I'm inferring it from this.
- **Collapse into clumps is documented.** In self-propelled particles that attract each other chemically, Pohl & Stark found that with pure attraction a chemotactic collaps occurs directly from the gas-like state before pronounced clusters are able to form. They also predict [this] by mapping our Langevin dynamics on the Keller-Segel model for bacterial chemotaxis. Saha, Golestanian and Ramaswamy describe "clumping" and gravitational collapse in chemotactic colloids.
- **The combination with alignment has been studied directly.** A 2025 paper investigate[s] the interplay between chemotaxis and alignment interactions, starting from a discrete model of self-propelled rods with chemotactic responses. Its alignment is the rod (nematic) version of the Vicsek model, not purely polar like yours, so it's the closest match rather than an exact one. It also notes that earlier work focused on either chemotaxis or alignment effects in isolation, so this is a recognised regime, not unexplored.
- **Tight, highly aligned clumps are familiar from alignment models.** Related models report clusters that behave like a single, highly persistent, self-propelled particle placed at its center of mass.

**Possible twist (I don't think it matters):** the clumps are strongly aligned (0.87) yet the agents barely move. This follows from how the rules are built. Once nearly everyone sits within distance 2 of each other, the alignment rule hands them all roughly the same heading, so high local alignment is almost automatic. The self-made chemical peak then pulls the clump back, which keeps net movement small. Pohl & Stark also describe a collapsed cluster that starts to fluctuate at the rim, then oscillates, but under a modified interaction, so that's only loosely related.

**Genuine or artifact?** The collapse is real emergent behaviour: it survives at 4× area and 3× run length and differs sharply from the knock-out. However, it is close to an absorbing state (a near-point clump, with CV ~11). The alignment score is largely a by-product of that geometry rather than flocking in its own right; alignment alone only reaches 0.34. I've marked it not trivial, with that caveat.

**Searches run:** "self-propelled particles Vicsek alignment autochemotaxis clustering collapse" and "active particles chemotaxis and polar alignment collapse into dense clusters simulation". I worked from search-result excerpts and didn't open the full papers.

---

## census:novelty:field-1

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 1 field(s) with diffusion [0.02]. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 8.0; (2) each agent deposits 0.1 of field 0 where it stands.

OBSERVED (features that passed the simulation checks):
  f_std, f_n_blobs, f_peak_sharp

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - f_std: median 0.6339 with interaction, 0.261 with interaction knocked out (difference 0.3729, robust spread 0.0313); difference at 4x area 1.0491, at 3x run length 0.3167
  - f_n_blobs: median 2.8332 with interaction, 3.9978 with interaction knocked out (difference -1.1646, robust spread 0.241); difference at 4x area -2.0352, at 3x run length -0.8343
  - f_peak_sharp (sharp field pattern): median 0.8488 with interaction, 0.395 with interaction knocked out (difference 0.4538, robust spread 0.1099); difference at 4x area 0.4287, at 3x run length 0.3433

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: f_std 0.0124; f_n_blobs 1.3863; f_peak_sharp 0.6416
  - rules 2 only: f_std 0.261; f_n_blobs 3.9978; f_peak_sharp 0.395
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Vicsek model flocking / clustering and band formation, read out through a passive (non-feedback) deposited tracer field

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Nothing in the program reads field 0, so it has no effect on the agents and only records a diffusion-blurred, time-integrated map of agent density. Higher f_std, fewer blobs and sharper peaks all come from known Vicsek clustering/banding (density about 0.1 with about 20 neighbours inside radius 8 puts the system in the ordered/clustered regime). The knock-out keeps deposition and removes alignment, so the measured differences just quantify Vicsek density inhomogeneity. The rule-1-only field features are meaningless because no deposition happens. The underlying flocking is genuine and robust, but the field features are an inert-tracer readout. Overturns the first pass from VARIANT to KNOWN, with no twist in the behaviour.

Citations:
- https://guava.physics.ucsd.edu/~nigel/Courses/Web%20page%20563/Essays_2017/PDF/Chatterjee.pdf
- https://arxiv.org/pdf/1511.01451
- https://arxiv.org/pdf/2604.28167
- https://arxiv.org/html/2610.26569

Full reply text:

**Verdict: KNOWN.** This is the standard Vicsek flocking model. The chemical field doesn't change the agents' behaviour at all; it just records where they have been. I'm overturning the first pass's VARIANT call: there is no new twist in the behaviour.

**Why it's equivalent to a documented system**
- Rule 1 is exactly the Vicsek alignment rule. Rule 2 adds 0.1 to field 0 where each agent stands. Nothing in the program ever reads that field, so it has no effect on how agents move. It is a passive tracer: a running map of agent density, blurred by diffusion.
- The settings are well inside the Vicsek ordered or clustered regime. 400 agents on 64×64 is a density of about 0.1, and an alignment radius of 8 covers roughly 20 neighbours per agent. Heading noise of 1 rad is moderate for that many neighbours.
- In this model, alignment drives density inhomogeneity. Reviews describe a first-order transition with an intermediate coexistence regime of high density ordered bands traveling in a low density disordered background. They also report that at low densities and low noise (ρ = 0.48, η = 0.1) the particles get grouped into clusters which move in random directions. Bands are characterized by a well-defined width, so that it is possible to accommodate several in the same system. Even the fully ordered phase can therefore exhibit substantial density inhomogeneity, as expected from density-orientation coupling in polar active matter.
- The three field measurements simply follow from this. Agents moving as dense flocks or bands lay down concentrated trails. That raises the field's variation (f_std), concentrates it into fewer blobs (f_n_blobs), and gives sharper peaks (f_peak_sharp). The knock-out keeps deposition but removes alignment, so agents scatter evenly. The gap between the two runs measures Vicsek clustering and nothing else.

**Checking the sub-program results**
- The "rule 1 only" values (f_std ≈ 0.01) are meaningless for these features: without deposition the field is essentially empty.
- "Rule 2 only" is the same program as the knock-out, which is why its numbers match exactly.
- So "no sub-program reproduces the full program" just means you need both alignment (to cluster) and deposition (to make the field visible). It doesn't show any synergy.

**Genuine pattern or artifact?** The clustering itself is a real, robust emergent pattern: Vicsek ordering, consistent at 4× area and 3× run length. But the features measured belong to an inert coupled field that only echoes that known pattern. Taken as field phenomena, they fall into the trivial or artifact category. If you wanted a real check of the agents, you'd measure their order and density directly.

I didn't search for "Vicsek plus deposited trail" papers specifically, because the field has no effect on the agents and so can't add any new behaviour.

---

## census:novelty:grid-8

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 2 field(s) with diffusion [0.02, 0.1]. Rules, applied in order each step: (1) each agent advances the type of the cell under it by one (cyclically); (2) each agent turns up the gradient of field 0.

OBSERVED (features that passed the simulation checks):
  g_activity_trend, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_activity_trend (activity rising): median -3.2875 with interaction, 0.0 with interaction knocked out (difference -3.2875, robust spread 0.0657); difference at 4x area -3.5441, at 3x run length -3.4345
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_activity_trend 0.0017; g_reversals 1.0
  - rules 2 only: g_activity_trend 0.0; g_reversals 0.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Random-walk bit-flipping agents (turmite without colour feedback) plus chemotactic gathering at maxima of a passive, emitter-free field

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: g_reversals equals the rule-1-only value (1.0) and comes straight from 2-state cyclic increment = flip on revisits. g_activity_trend (negative) comes from agents gathering at maxima of field 0, which nobody emits into: fewer distinct cells get flipped, and same-step double flips cancel. No feedback from cells to agents or from agents to fields, so nothing self-organises; field 1 does nothing. Overturns the first pass only in the less-novel direction (KNOWN and trivial, not a meaningful variant).

Citations:
- https://arxiv.org/pdf/1702.05547
- https://en.wikipedia.org/wiki/Langton%27s_ant
- https://www.nature.com/articles/s41598-021-83963-x
- https://arxiv.org/pdf/1202.2442

Full reply text:

**Verdict: KNOWN.** Both features are mechanical results of the rules, not emergent patterns. This overturns the first pass in the "less novel" direction: it's not a variant worth flagging, and the measured features are trivial.

**1. g_reversals ("cells flip back and forth") comes straight from rule 1.** With only 2 cell types, "advance the type by one, cyclically" just means "flip the cell." So any cell an agent visits twice flips back. The sub-program data shows this: rule 1 alone gives the full 1.0, and rule 2 adds nothing. This is the cell-writing half of Langton's ant / turmites, where each agent flips the cell under it, but here the agents ignore the cell's colour. Turmites with more colours "augment the symbol by 1", which is exactly rule 1. Without the feedback from cell colour to turning, the complex Langton behaviour can't happen. What's left is random walkers flipping bits, so reversals are guaranteed. Calling this "interaction" only reflects the knock-out removing the agent→cell writing itself.

**2. g_activity_trend (negative, so activity is falling) is agents piling up on a field nobody feeds.** No rule adds anything to either field, so field 0 must be some initial pattern slowly smoothing out by diffusion. If it were flat, rule 2 would do nothing and the full program would match rule 1 alone, which it doesn't. Agents that steer up a gradient collect at local maxima. That's standard: "chemotactic cells do not disperse upon reaching a chemoattractant source... chemotactic bacteria seek out and accumulate at nutrient maxima". Once crowded, fewer distinct cells get flipped. Also, two agents on the same 2-state cell in one step flip it twice, leaving it unchanged. Both effects push net flip activity down over time. That's a geometric consequence of agents gathering at fixed hotspots, not self-organisation. There's no Keller–Segel-style feedback, because agents never emit anything. Field 1 does nothing at all.

**Is it trivial or an artifact?** Yes:
- Reversals follow from the 2-state flip rule.
- The activity decline is passive aggregation on a pre-set (initial) landscape, with field 1 doing nothing.
- The effects being stable across seeds and sizes just shows they're deterministic consequences, not that they're interesting.

**One caveat:** I don't have the field's initial conditions. If something not described (an initial random field, or hidden emission) is setting up the landscape, that would explain rule 2's effect. It wouldn't change the verdict.

**Searches:** chemotactic active particles collecting at maxima of an external field; turmite / Langton's ant agents flipping cells. Both behaviours are covered by those standard frameworks.

---

## qd:027

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.02]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns to the average heading of agents within distance 1.0.

OBSERVED (features that passed the simulation checks):
  a_density_cv, a_local_align, a_step, a_polar_std

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_density_cv: median 8.4486 with interaction, 0.796 with interaction knocked out (difference 7.6526, robust spread 0.1957); difference at 4x area 4.2028, at 3x run length 14.5057
  - a_local_align (local alignment): median 0.9106 with interaction, -0.0351 with interaction knocked out (difference 0.9457, robust spread 0.0204); difference at 4x area 0.942, at 3x run length 0.9817
  - a_step: median 0.098 with interaction, 5.9253 with interaction knocked out (difference -5.8273, robust spread 0.1205); difference at 4x area -5.6355, at 3x run length -17.1298
  - a_polar_std: median 0.1204 with interaction, 0.0133 with interaction knocked out (difference 0.1071, robust spread 0.0048); difference at 4x area 0.1067, at 3x run length 0.1249

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: a_density_cv 8.2145; a_local_align 0.0556; a_step 0.2416; a_polar_std 0.0393
  - rules 2 only: a_density_cv 4.1452; a_local_align 0.6942; a_step 5.9581; a_polar_std 0.0445
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Vicsek alignment + autochemotaxis (Keller-Segel-type chemotactic collapse into aligned, trapped aggregates / 'active droplets')

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: The emitter is not stated in the rules; agent emission is inferred because rule 1 alone produces strong clumping. The features describe chemotactic collapse into dense, internally aligned, nearly stationary clumps; rule 2 alone gives Vicsek bands. Density CV keeps rising at 3x run length, so this is ongoing collapse/coarsening rather than a steady state. High local alignment is partly automatic once agents sit within one alignment radius. Real emergent clustering, not an artifact, but close to a degenerate collapsed state. Overturns 'variant' only in that no extra phenomenological twist was found.

Citations:
- https://www.nature.com/articles/s41467-022-34484-2
- https://arxiv.org/html/2506.18263
- https://www.sciencedirect.com/science/article/pii/S0006349518312268
- https://arxiv.org/html/2011.04580v1

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass. Nothing about what you observed is a new behaviour that would justify calling it a variant.

**What the system is.** The rules leave out who produces field 0. The rule-1-only sub-program still clumps very strongly (density CV 8.2 against 0.8 for the knock-out), so the field must be produced by the agents themselves. That emission is a one-way rule, so it stays in the knock-out. This makes the system a Vicsek flock (self-propelled agents copying their neighbours' heading, with noise) plus autochemotaxis: agents lay down a chemical and steer toward it. Each part, and the two together, are well studied:
- **Chemotaxis plus alignment in one model:** a 2025 paper adds chemotaxis to a self-propelled rod model, which is itself a variant of the Vicsek model. It then derives coarse-grained equations for density, polar order, nematic order and the attractant field. "We investigate the interplay between chemotaxis and alignment interactions in active rod-like particles... Starting from a discrete model of self-propelled rods with chemotactic responses"
- **Agents that align with a signal field they make collectively:** Ziepke et al. (Nat. Commun. 2022) found that order appears at much lower densities than the Vicsek transition, due to the alignment of the polar particles with the collectively established signaling field. When self-propulsion is weak the steady-state is governed by active droplets, and for vanishing signal sensitivity, the model reduces to a Vicsek-type model, and only polar bands are found. That last case matches your rules-2-only result.
- **Vicsek-type agents with chemotaxis forming stable aggregates:** the *Myxococcus* agent-based model is based on the well-known Vicsek model. In it, a chemotaxis model with adaptation can reproduce the observed experimental results leading to the formation of stable aggregates and reproduces the experimentally observed patterns of cell alignment around aggregates.
- **Related field-theory result:** a chemotaxis model notes that the higher the density or the strength of chemotactic interactions, the greater the chance of forming large stationary aggregates, due to the cooperative aligning interactions among closely spaced neighbors.

**How your measurements fit.** Together the features describe agents that collapse into one or a few dense, nearly stationary, internally aligned clumps:
- density CV is very high;
- local alignment is about 0.91;
- step size is tiny (0.098), so the clump is trapped by its own chemical peak;
- global polar order is weak but present (0.12).

Each rule alone gives part of this: rule 1 gives the clumping without alignment, and rule 2 gives alignment and moving bands without trapping. The combination is the "active droplet" / aligned-aggregate regime described above. I found no extra phenomenon that needs explaining.

**Genuine or artifact.** It is real emergent clustering, not random noise or a numerical blow-up, but it is close to a degenerate state:
- **Collapse, not a steady state:** density CV keeps growing with run time (difference 14.5 at 3× run length against 7.7 at the original length). That looks like ongoing collapse or merging, the familiar Keller–Segel-style aggregation, rather than a settled pattern.
- **Alignment may be close to automatic:** if most agents end up inside one alignment radius (1.0), averaging headings will line them up almost regardless. Local alignment near 1 inside a collapsed clump says little on its own.
- **Missing evidence:** without snapshots I can't tell whether the clump mills (rotates), translates slowly, or sits frozen. A rotating clump would match the documented vortex/droplet states but would still not be new.

---

## qd:033

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 1 type (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns left on an even-type cell and right on an odd-type cell.

OBSERVED (features that passed the simulation checks):
  a_density_cv, a_nn_ratio, a_local_align

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_density_cv: median 1.1586 with interaction, 0.8089 with interaction knocked out (difference 0.3496, robust spread 0.0558); difference at 4x area 0.4652, at 3x run length 0.5886
  - a_nn_ratio (clumping (below 1 = clustered)): median 0.826 with interaction, 0.9926 with interaction knocked out (difference -0.1666, robust spread 0.0364); difference at 4x area -0.2046, at 3x run length -0.2699
  - a_local_align (local alignment): median 0.2566 with interaction, -0.0405 with interaction knocked out (difference 0.2971, robust spread 0.0259); difference at 4x area 0.3044, at 3x run length 0.3409

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: a_density_cv 11.2227; a_nn_ratio 0.0088; a_local_align -0.0074
  - rules 2 only: a_density_cv 0.8009; a_nn_ratio 0.9973; a_local_align -0.0321
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Chemotactic (autochemotactic / Keller-Segel-type) aggregation of active Brownian particles, softened by a fixed random torque field acting as extra direction noise, with emergent alignment from streaming toward aggregates

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Rule 2 alone matches the knock-out (inert fixed random torque = extra direction noise). Rule 1 alone gives chemotactic collapse (nn_ratio 0.009). The full program sits in between: noise-limited aggregation, with local alignment from agents streaming toward diffuse aggregates, which is documented in agent-based chemotaxis models. The field's source is not stated; if agents emit, this is standard autochemotaxis, and if the field is sourceless, it is gradient-climbing on a landscape. Either way the behaviour is known. The overturn is mild (VARIANT -> KNOWN) because no phenomenological twist was found.

Citations:
- https://arxiv.org/html/2409.04262v1
- https://arxiv.org/html/1403.4063v1
- https://www.sciencedirect.com/science/article/pii/S0010465526002675
- https://www.sciencedirect.com/science/article/pii/S0006349518312268
- https://arxiv.org/pdf/2303.07746

Full reply text:

**Verdict: KNOWN.** I'm moving the first pass from "variant" to "known". The behaviour is chemotactic aggregation of self-propelled agents, made weaker by the cell lattice. The alignment comes from agents streaming toward the aggregates, which is also documented. Nothing in the phenomenology looks new.

**What the sub-programs show**
- **Rule 2 does nothing on its own.** With only the cell-steering rule, the numbers (density CV 0.80, clumping ratio 0.997, alignment −0.03) match the knock-out (0.81, 0.993, −0.04). The cells are a fixed random checkerboard of left/right turns, so on their own they just add direction noise. That is "quenched disorder" in active-matter terms. Here it acts only as a weak pre-set turning force, not as obstacles or traps.
- **Rule 1 alone makes agents collapse into tight clusters.** Density CV is 11.2 and the clumping ratio is 0.009, which means near point-like clusters. The rule text never says what produces field 0. The most likely reading is that agents emit it, which makes this the standard self-attracting-chemical (autochemotactic) particle model. It is literally active Brownian particles emitting a diffusive chemical cue that serves as a chemotactic signal for all particles. It reproduces the known result that autochemotactic cells typically conglomerate in large aggregates or even exhibit a chemotactic collaps. If the field instead has no source and simply relaxes, agents are just climbing a fixed landscape, which is even more trivially known.
- **The full program sits in between.** The extra turning from rule 2 blocks the full collapse and leaves moderate clusters (CV 1.16, clumping ratio 0.83). This fits the general picture of chemotactically induced clustering, dynamic bands, motility-dependent phase separation, and collapse-like instabilities, where noise competes with aggregation.
- **The local alignment (about 0.26) is the only feature neither rule shows alone, and it is also documented.** Agents heading toward the same aggregate line up with each other. In tight clusters (rule 1 alone) headings scramble, so alignment drops to zero. In the full program, agents keep travelling toward diffuse aggregates and stay aligned. Agent-based chemotaxis models report this: one study says its model reproduces the experimentally observed patterns of cell alignment around aggregates, with cells lining up radially toward the nearest aggregate (M. xanthus). Alignment arising without any explicit alignment rule has also been shown for attracting active particles (even non-aligning attractive interactions can lead to a flocking state).

**Genuine pattern or artifact?** The effect is real, not an artifact. The gaps from the knock-out are 6 to 11 times the spread between seeds, and they get bigger with 4× area and 3× run length. It is not a frozen absorbing state, a diverged field or a random colouring. But the full program's interest is only that it lands partway between collapse and the noise baseline, because rule 2 is acting as extra noise. Two caveats:
- I couldn't confirm whether agents emit the field.
- The alignment measure may simply be picking up streams heading into clusters. I didn't find a paper using this exact checkerboard-turning setup, but per the brief, a new rule combination doesn't make a discovery.

**Searches run:** autochemotactic active Brownian particle clustering; active particles in quenched disorder or with environmental torques; alignment emerging from chemotactic streaming toward aggregates.

---

## qd:036

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.02]. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 3.0; (2) each agent deposits 0.1 of field 0 where it stands.

OBSERVED (features that passed the simulation checks):
  f_std, f_n_blobs

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - f_std: median 4.2511 with interaction, 0.4394 with interaction knocked out (difference 3.8117, robust spread 0.1164); difference at 4x area 5.0181, at 3x run length 3.4702
  - f_n_blobs: median 1.0986 with interaction, 4.0247 with interaction knocked out (difference -2.9261, robust spread 0.5274); difference at 4x area -3.3949, at 3x run length -2.627

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: f_std 0.0132; f_n_blobs 1.6094
  - rules 2 only: f_std 0.4394; f_n_blobs 4.0247
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Vicsek model flocking: density clustering / traveling bands / giant number fluctuations, recorded in a passive (inert, never-read) accumulated deposit field

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: No agent reads field 0 and it does not decay, so it is a time-integrated record of agent density. Higher f_std and fewer blobs only show that Vicsek alignment makes agents cluster into a few coherent flocks, which is standard Vicsek phenomenology. The flocking itself is real; the field features are a readout of an inert coupled field, so they are trivial. Because the field never decays, f_std also depends partly on run length. The first pass called this a VARIANT; it is plain KNOWN.

Citations:
- https://arxiv.org/html/2312.15525v2
- https://guava.physics.ucsd.edu/~nigel/Courses/Web%20page%20563/Essays_2017/PDF/Chatterjee.pdf
- https://arxiv.org/html/2604.28167v2
- https://arxiv.org/html/1709.01499
- https://www.researchgate.net/publication/283531761_The_Physics_of_the_Vicsek_model

Full reply text:

**Verdict: KNOWN. This overturns the first pass's "variant" label.** The behaviour is ordinary Vicsek flocking: agents bunch into dense groups that move together. The field doesn't add a new mechanism. It only records where those groups have been.

**1. The field is a passive tracer.** Rule 1 aligns each agent with its neighbours. Rule 2 has agents drop chemical into the field. Nothing ever reads the field, and it has no decay. So the field is just the agents' density, added up over time and blurred slightly by diffusion (0.02). Neither field measurement (f_std, f_n_blobs) can show anything the agents aren't already doing. The sub-programs confirm this:
- **Rule 1 only (alignment, no depositing):** f_std is about 0.01, because almost nothing gets deposited.
- **Rule 2 only (depositing, no alignment):** independent noisy walkers leave a fairly even trace, with many small blobs (f_n_blobs ≈ 4.0; if this is a log count, that's roughly 55 blobs).
- **Both rules:** agents cluster into a few moving flocks. They leave strong, concentrated tracks: f_std ≈ 4.25, and f_n_blobs ≈ 1.1, roughly 3 blobs on the same reading.

**2. That clustering is textbook Vicsek physics.** The density here is low (400 agents on 64×64), heading noise is low (0.3 rad), and the alignment radius is 3. In that regime the Vicsek model is known to split into dense ordered groups. The literature describes a coexistence phase that shows up as either multiple bands of particles moving collectively, a phenomenon known as microphase separation, or a polar-ordered cross sea phase, primarily driven by giant number fluctuations (GNF). The transition is actually 1st order, with an intermediate coexistence regime of high density ordered bands traveling in a low density disordered background. At very low noise, a single big cluster has also been reported ("Emergence of a single cluster in Vicsek's model at very low noise", arXiv 1709.01499). Even the ordered phase can have large density fluctuations while remaining locally aligned throughout the system. A trace field with high variance and a few blobs is exactly what you'd expect from any of these states.

**3. Genuine or artifact?** The underlying flocking is genuine, but it is already well documented. The two features that were actually measured are artifacts of an inert coupled field. Because the field never decays, f_std also depends on how long the run is, though the effect held up at 3× run length (difference 3.47 vs 3.81). I can't name a paper that looks at this exact trace-field picture, but that's only a difference in rules. The trace adds no new behaviour, so I don't see a twist that would justify KNOWN-VARIANT.

---

## qd:041

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent advances the type of the cell under it by one (cyclically).

OBSERVED (features that passed the simulation checks):
  g_activity_trend, g_activity_cv, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_activity_trend (activity rising): median -4.5106 with interaction, 0.0 with interaction knocked out (difference -4.5106, robust spread 0.0902); difference at 4x area -4.4088, at 3x run length -3.5965
  - g_activity_cv (bursty activity): median 1.5212 with interaction, 0.0 with interaction knocked out (difference 1.5212, robust spread 0.055); difference at 4x area 0.8032, at 3x run length 3.9367
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_activity_trend 0.0; g_activity_cv 0.0; g_reversals 0.0
  - rules 2 only: g_activity_trend -0.158; g_activity_cv 0.0781; g_reversals 1.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Gradient-ascent (chemotactic) aggregation of agents in a passive, decaying, non-emitted field, plus parity toggling of a 2-state substrate (Langton's-ant-like marking)

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Field 0 has no source, so it only diffuses from its random start; rule 1 is passive hill-climbing to the initial maxima. Falling activity = agents crowding onto the same cells, with toggles cancelling under the 2-state parity. High CV = small-mean relaxation fluctuations (effect halves at 4x area). g_reversals=1 is guaranteed for 2-state cyclic toggling (rule-2-only also gives 1.0). The knock-out deletes all cell updates, so its zero baseline is uninformative; the fair baseline is rule-2-only. 'Overturns' only because the first pass's VARIANT label overstates novelty; the family identification stands. Reasoning based on the rule description; not rerun.

Citations:
- Chatterjee & Goldenfeld, Field-theoretic model for chemotaxis in run and tumble particles, PRE 2021 (aggregation in externally imposed attractant field)
- Schweitzer, Brownian Agents and Active Particles (2007) / Helbing et al. PRE 56, 2527 (1997) - active walker framework, here lacking the feedback loop

Full reply text:

**Verdict: KNOWN. The interesting-looking features are mostly an artifact.** I agree with the first pass that this is a known model. But "variant" gives it too much credit: nothing in the behaviour is new, and most of the signal comes from how the features and the knock-out are set up. Below I reason from the rule description; I did not rerun the simulation.

**What is actually happening**
1. **Nothing feeds the field.** The rules list no emission, deposit or coupling into field 0. It only diffuses (D = 0.05) from its random starting state, so it is a passive, flattening landscape. Rule 1 is then just agents climbing toward the high points of a fixed, smoothing field. Agents aggregating in an imposed attractant field is standard: one paper says its model "exhibits particle aggregation in an externally imposed constant attractant field, as is observed for phototactic or thermotactic agents" . There is no feedback loop, so this is not even the "active walker" class, where the walkers modify their environment which, in turn, influences them at a later time. Here the cells are changed by the agents but are never read by them or by the field.
2. **Why activity falls (g_activity_trend about −4.5).** As the agents bunch up near the field's peaks, they keep hitting the same few cells. With two cell types, "advance cyclically" just means toggle, so two agents on one cell in the same step cancel out. Net cell changes per step therefore drop as crowding increases. That is crowding plus parity cancellation, not a new dynamic.
3. **Why activity looks bursty (g_activity_cv about 1.5).** Once activity is low, a small number of agents wandering near the flattening peaks makes the count fluctuate a lot relative to its mean. The CV is high because the mean is small. The values support this: the effect shrinks at 4x area (0.80) and grows at 3x run length (3.94), which looks like a relaxation process rather than a steady bursty state.
4. **Why cells "flip back and forth" (g_reversals = 1.0).** This is guaranteed. With only two types, any cell changed twice has reversed. The rule-2-only program also scores exactly 1.0, so this feature tells us nothing.
5. **The knock-out comparison is misleading.** Removing rule 2 removes all cell changes, so every activity measure is zero by construction. The fair baseline is the rule-2-only program, where agents wander randomly: steady, roughly flat toggling (trend −0.16, CV 0.08). The only real difference from adding rule 1 is that agents gather onto the passive field's peaks. That is ordinary gradient-ascent aggregation, which also reduces how many distinct cells get toggled.

**Phenomenology:** gradient-ascent aggregation of agents in a passive, decaying random field, plus Langton's-ant-style toggling of the cells beneath them. Nothing is left unexplained, and I found no new twist.

---

## qd:065

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 4 types (random start). Agents move at speed 2.0 with heading noise 0.0 rad per step. 2 field(s) with diffusion [0.02, 0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) field 1 grows autocatalytically by consuming field 0 (rate 1.0 x f x g^2).

OBSERVED (features that passed the simulation checks):
  a_density_cv, a_local_align, a_nn_ratio

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - a_density_cv: median 2.0015 with interaction, 0.8207 with interaction knocked out (difference 1.1808, robust spread 0.065); difference at 4x area 1.1839, at 3x run length 1.1505
  - a_local_align (local alignment): median 0.4133 with interaction, -0.0494 with interaction knocked out (difference 0.4627, robust spread 0.019); difference at 4x area 0.6058, at 3x run length 0.4377
  - a_nn_ratio (clumping (below 1 = clustered)): median 0.4609 with interaction, 0.9979 with interaction knocked out (difference -0.537, robust spread 0.0292); difference at 4x area -0.7118, at 3x run length -0.5709

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: a_density_cv 6.3703; a_local_align 0.7919; a_nn_ratio 0.0164
  - rules 2 only: a_density_cv 0.8207; a_local_align -0.0494; a_nn_ratio 0.9979
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Chemotactic gradient-ascent aggregation of non-interacting self-propelled particles in a one-way-coupled (externally evolving) attractant field; field = unfed Gray-Scott (U+2V->3V) depletion transient

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Agents only read field 0; they neither emit nor consume, and do not sense each other, so the knock-out removes single-agent chemotaxis, not agent-agent interaction. Clustering and alignment come from agents climbing the same field and gathering at its maxima. Rule 2 is Gray-Scott without feed or decay, so field 0 is used up and its gradients keep moving and weakening; this explains why the full program clusters less than rule 1 alone. With zero noise, late-time steering may follow near-round-off gradients; check the size of field 0's gradient over time.

Citations:
- https://arxiv.org/html/2011.04580v1
- https://itp.uni-frankfurt.de/~gros/StudentProjects/Projects_2020/projekt_schulz_kaefer/
- https://www.mrob.com/pub/comp/xmorphia/index.html

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass, which called this a variant. The behaviour is what happens when non-interacting particles climb the gradient of a field. It isn't a new collective pattern, and it has the marks of a trivial outcome.

**What the rules actually do**
- **The coupling only goes one way.** Agents read field 0 but nothing in the rules lets them emit, consume or otherwise change any field. They don't sense each other either. So there is no feedback loop of the Keller–Segel kind, where particles make the attractant they gather on. Each agent just climbs a field it has no effect on. The "interaction knock-out" removes that reading of the field, so the large differences it shows only mean "chemotaxis on vs. off". They don't show that agents interact with each other.
- **The field part is the Gray–Scott model without its feed and decay terms.** The reaction term f·g² is the Gray–Scott reaction U + 2V → 3V. In the full model, the amount of both substances has to be controlled to keep the reaction going, which is done by adding U at the "feed rate" F and removing V at the "kill rate" k. Here there is no F and no k, so field 0 is simply used up. That gives a transient that relaxes towards a flattening field 0, not the steady Gray–Scott spots or waves. Field 0's gradients therefore keep moving and weakening, and that alone explains why the full program clusters less than the chemotaxis-only program (density CV 2.0 vs 6.37, nearest-neighbour ratio 0.46 vs 0.016).

**Matching the three features to known behaviour**
- **Clustering (`a_density_cv`, `a_nn_ratio`):** Particles in a field they don't alter build up at its maxima. A run-and-tumble chemotaxis model exhibits particle aggregation in an externally imposed constant attractant field, as is observed for phototactic or thermotactic agents. Its particles start out moving randomly but eventually aggregate in the center of the arena where the attractant concentration is maximum.
- **Local alignment (`a_local_align`):** Agents near each other sense nearly the same gradient, so their headings line up. That alignment comes from the shared field, not from agents influencing each other. With no heading noise and speed 2, the rule-1-only run gives the strongest version (alignment 0.79).
- **Reaction–diffusion part:** this is the standard U + 2V → 3V autocatalysis that converts U into V; without the feed term it just runs down.

**Is anything new in the phenomenology?** No. The full program only weakens the rule-1-only clustering because the attractant is being depleted and moved around. No feature is new in kind.

**Trivial or artifact?** I'd call it trivial rather than a numerical artifact. It amounts to passive gradient-ascent particles following a field that is running down. There's also a possible numerical artifact: with zero noise and field 0 heading towards depletion, late-time steering may be following tiny, near-round-off gradients. This should be checked by tracking the size of field 0's gradient over the run.

---

## review:009

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=0.1  (compact grammar notation)

OBSERVED (features that passed the simulation checks):
  p_r, p_freq_spread

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - p_r (global phase order (full sync)): median 0.4459 with interaction, 0.0665 with interaction knocked out (difference 0.3794, robust spread 0.021); difference at 4x area 0.2922, at 3x run length 0.7323
  - p_freq_spread (spread of oscillator natural frequencies): median 0.0251 with interaction, 0.5146 with interaction knocked out (difference -0.4895, robust spread 0.0108); difference at 4x area -0.4672, at 3x run length -1.3797

### Stage 3 reply: verdict KNOWN, standard name: Locally (nearest-neighbour) coupled Kuramoto lattice with stochastic on-off (blinking) coupling: frequency entrainment with partial, size-dependent, slowly coarsening phase order

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Overturned from KNOWN-VARIANT to KNOWN: random gating with p=0.1 is just a blinking network, and in the fast-switching limit it behaves like a lattice Kuramoto model with effective coupling of about pK. No new phenomenology. The p_freq_spread metric must be measuring effective frequencies, not natural ones (natural frequencies can't change under coupling); check its label. Partial r that falls with area and rises with run length matches documented finite-size effects and vortex coarsening in the 2D lattice Kuramoto model.

Citations:
- Hong, Park & Choi, PRE 72, 036217 (2005) - collective synchronization in spatially extended oscillator lattices
- Noise-induced synchronization in the Kuramoto model on finite 2D lattice, arXiv:2004.00294
- L.-H. Tang, finite-size scaling in the Kuramoto model, arXiv:1011.5823
- Frequency Spirals in the 2D Kuramoto Lattice, arXiv:1512.02335
- Jeter & Belykh, Synchronization in On-Off Stochastic Networks, IEEE TCAS-I 62, 1260 (2015)

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass, which called this a variant. Once you average out the random switching, this is the nearest-neighbour (lattice) Kuramoto model. Both observed features are textbook results for that model: partial global phase order that depends on system size, and frequency entrainment.

**What the rule is.** Each oscillator, with probability 0.1 per step, nudges its phase toward one neighbour's with strength K=0.2. On average this is local Kuramoto coupling with effective strength of roughly pK ≈ 0.02 per step. Random on/off coupling like this has its own literature ("blinking" networks). The standard result is that when switching is fast, the network synchronizes as long as synchronization in the averaged network, obtained by replacing the random variables by their mean, becomes stable. So the stochastic gating changes the rule, not the kind of behaviour. That is the "different rules, same behaviour" situation you warned about.

**Why each feature is documented:**
- **p_freq_spread drops (≈0.51 → 0.025).** This is frequency entrainment. For local coupling, there exists a possibility that the frequencies of a large proportion of the oscillators are synchronized, but the phases remain spread out across this synchronized cluster. Hong, Park and Choi (PRE 72, 036217, 2005) found that phase synchronization occurs above a lower critical dimension of d=4, while frequency entrainment has a lower critical dimension of d=2. One caution: in a standard Kuramoto model the natural frequencies are fixed. If this metric really measured natural frequencies, the coupling couldn't shrink it. It is almost certainly measuring effective (observed) frequencies, so the label is probably wrong.
- **p_r is partial (≈0.45), falls at 4x area and rises at 3x run length.** This matches the known picture for finite lattices. Phase locking has been shown to be impossible on infinite lattices, yet simulations tended to show a strong degree of average synchronization and entrainment in finite grids. In 2D, the entrainment crossover takes place in a two-dimensional lattice via the stability of vortex-antivortex pairs formed in the phase-field. That fits slow coarsening: order keeps growing with run time and is weaker in bigger systems. Locally coupled oscillators exhibit two types of frequency entrainment, including clumping into compact synchronized domains in low dimensions (Tang).

**Trivial or artifact?** Neither. The effect beats the knock-out by more than 15 robust spreads, survives 4x area and 3x run length, and the frozen frequencies with partial phase order are the classic signature of a locally coupled oscillator lattice. It is not an inert field, a divergence or an absorbing state. The fall in r with system size is expected finite-size behaviour, not a new effect.

**What I couldn't check.** The grid type (square vs. other), the frequency distribution, the time step, and whether the 2π-jump correction is handled correctly. None of these would change the verdict.

---
