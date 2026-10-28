# Stage 3 (b): author review packet

Model labels are provisional (VERIFICATION_PLAN.md, Stage 3). For each phenomenon: check that the cited sources exist and describe the behaviour as stated, then record a final label (documented / closely related / undocumented candidate). Order: closely related first, then documented and judged non-trivial, then documented and judged trivial or artifactual.

84 phenomena: 80 documented, 4 closely related.

## review:046: CLOSELY RELATED (KNOWN-VARIANT)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (every step); (2) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step).
- Standard name given: Synchronous majority-vote coarsening interrupted by cyclic (rock-paper-scissors) invasion: arrested small domains with period-2-like boundary flipping (majority CA + cyclic dominance; cf. Potts-energy cyclic voter model)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 4
- Twist (if any): 
- Notes: Overturns the first pass's reasoning, not its 'known' conclusion. The rule-2-only numbers match the knock-out (zero reversals, no clustering), so the CCA spiral/demon explanation does not apply: rule 2 alone appears to advance the whole random colouring together, which is trivial. The full system's features come from majority smoothing (period-2 flipping, clustering) with cyclic invasion stopping coarsening at about 3 cells. Twist: deterministic synchronous arrest at very small scale instead of spirals or large domains. Probably a short repeating state on a finite board, not a new mechanism.
- Citations:
  - Nakar & Ron, arXiv:2205.08972 (synchronous majority CA: fixed points or period-2 cycles; Goles-Olivos 1981, Poljak-Sura 1983)
  - Wolfram Demonstration, Cellular Automata with Majority Rule (patchy coarsening)
  - Szabo & Szolnoki, arXiv:cond-mat/0110255 (three-state cyclic voter model with Potts energy)
  - Szabo & Borsos, arXiv:1508.03147 (synchronous logit RPS: stabilised oscillating domains)
  - Fisch, Gravner, Griffeath, Cyclic Cellular Automata in Two Dimensions (debris/droplets/defects/demons)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:050: CLOSELY RELATED (KNOWN-VARIANT)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Cyclic cellular automaton (Griffeath), 4 colours, Moore range-1, threshold 1 (incoherent period-4 local oscillation) coupled with Schelling-type relocation; nearest analogue: cyclic-dominance lattice models with mobility
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- Notes: Each part alone is textbook: rule 1 = CCA T=1 period-4 incoherent oscillation (activity 1.0, Moran 0.07); rule 2 = Schelling freezing (activity 0, Moran 0.23). Coupled system gives a lasting partly active (~61%), more clustered (Moran 0.32) state, probably because empty is a state in the cycle, so vacancies keep changing and Schelling cannot freeze. No paper found on CCA plus Schelling coupling specifically; mobility-in-cyclic-dominance literature covers the general effect. Overturns the first pass's claim that segregation 'adds little'; agrees with its identification of the rule-1 phase. Caveat: Moran's I on a cyclic, ordinal type label that includes the empty state is partly arbitrary; ~6 domains is a weak signal.
- Citations:
  - Shalizi et al., Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067
  - Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39 (1991)
  - Hiltemann, Multi-coloured Cellular Automata (LIACS thesis, 2008)
  - MCell rules lexicon, Cyclic CA (psoup.math.wisc.edu)
  - Reichenbach, Mobilia & Frey, Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, arXiv 0709.0217
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:051: CLOSELY RELATED (KNOWN-VARIANT)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity, g_interface
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell; (2) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step).
- Standard name given: Schelling segregation coupled to a 4-colour threshold-1 Moore cyclic cellular automaton (cyclic dominance plus an ordering drive, like the Potts + rock-paper-scissors model)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- Notes: The first pass explained everything as the pure CCA T=1 quenched period-4 state. That fits the rule-2-only run (activity 1.0, near-random boundary measure 0.72) but not the full program (activity 0.67, boundary measure 0.38, about the Schelling-only value, fewer domains), each off by many robust spreads. The full behaviour is a lasting clustered state: Schelling patches eaten from their edges by cyclic invasion and rebuilt by moves, with idle patch interiors. The general mechanism (ordering drive plus cyclic dominance gives a lasting changing pattern) is documented in Potts+RPS and in switching-agent Schelling models. The twist, vacancies that are themselves a cycling state, giving partial activity of about 2/3 and Schelling-scale domains, was not found but is likely minor. The pattern is genuine (not frozen, random or diverged) but the features are generic.
- Citations:
  - https://arxiv.org/pdf/nlin/0507067
  - https://arxiv.org/html/nlin/0507067
  - https://journals.aps.org/pre/abstract/10.1103/PhysRevE.71.027102
  - https://www.researchgate.net/publication/233947119_A_Schelling_model_with_switching_agents_Decreasing_segregation_via_random_allocation_and_social_mobility
  - https://en.wikipedia.org/wiki/Cyclic_cellular_automaton
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:064: CLOSELY RELATED (KNOWN-VARIANT)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_moran, g_activity
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 2 neighbours already are t+1 (every step); (2) neighbouring cells swap contents (every step).
- Standard name given: Griffeath cyclic cellular automaton (4 colours, Moore range-1, threshold 2) with superimposed neighbour-exchange (stirring) mobility; cf. mobility-enlarged spiral/domain patterns in cyclic (rock-paper-scissors) lattice models
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 4
- Twist (if any): 
- Notes: Downgraded from KNOWN to KNOWN-VARIANT. The first pass said swaps break up the CCA spirals, but the data show the opposite: swaps raise correlation length (2 to 4.5) and Moran's I (0.47 to 0.80) and halve activity (0.99 to 0.48). That direction fits the documented mobility-enlarges-patterns effect in cyclic competition models, but I found no paper adding moderate random swapping to the threshold CCA. The full-interaction knock-out is just the frozen random start, so the single-rule programs are the meaningful baselines; the full program beats both. Not trivial: no absorbing state, strong clustering, robust to size and run length.
- Citations:
  - Shalizi & Shalizi, Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067 (citing Fisch, Gravner, Griffeath 1991)
  - Reichenbach, Mobilia, Frey, Nature 448, 1046 (2007): spiral wavelength grows ~sqrt(mobility)
  - Reichenbach, Mobilia, Frey, PRL 99, 238105 (2007)
  - Durrett & Neuhauser, Ann. Probab. 22, 289 (1994); Durrett, Ten Lectures on Particle Systems (stirring = neighbour exchange)
  - Roman, Konrad, Pleimling, J. Stat. Mech. P07014 (2012): four-species cyclic competition with exchange
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:agents-1: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: a_density_cv, a_step, a_nn_ratio, a_local_align
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. Rules, applied in order each step: (1) each agent turns toward the centre of any type agents within distance 1.0.
- Standard name given: Attraction-only (cohesion-only) self-propelled particle aggregation / swarm phase (local attraction model; minimal cognitive flocking model with full vision cone)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: All four features (high density CV, NN ratio ~0.18, suppressed per-step displacement from cluster trapping, mildly negative local alignment from agents overshooting the centre and returning inside aggregates) are the standard disordered-aggregate/swarm phase of attraction-only SPPs. The negative local alignment value is not explicitly quantified in the sources I checked but follows directly from centre-seeking without alignment; it is not a new phenomenon. Robust to size and run length, so a genuine emergent pattern, not an artifact. The density-CV effect roughly halves at 4x area, consistent with a finite number of clusters in a larger box.
- Citations:
  - Strömbom D. (2011) Collective motion from local attraction, J. Theor. Biol. 283:145-151
  - Barberis L. & Peruani F. (2016) Large-scale patterns in a minimal cognitive flocking model: incidental leaders, nematic patterns, and aggregates, Phys. Rev. Lett. 117:248001
  - Negi R.S., Winkler R.G., Gompper G. (2022) Emergent collective behavior of active Brownian particles with visual perception, Soft Matter 18:6167
  - Emergence of Active Foams in a Minimal Cognitive Flocking Model, arXiv:2610.15835
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:agents-3: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: a_density_cv, a_local_align
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 3.14159 rad per step. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 1.0.
- Standard name given: Vicsek model (standard angular-noise version), disordered gas phase or near the order-disorder transition, with weak local alignment and alignment-induced density clustering
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: The rule is the textbook Vicsek model. At rho=1 and eta~pi the system is in or near the high-noise disordered gas: weak local alignment (~0.17) and mildly above-random density variation, both of which disappear when the interaction is knocked out. The effect is real but weak and well documented. Caveat: if the noise is uniform on the full [-pi,pi] and added after alignment, any measured alignment must come from measuring before the noise, which would make it partly true by construction. Not a discovery.
- Citations:
  - Vicsek, Czirok, Ben-Jacob, Cohen, Shochet, Phys. Rev. Lett. 75, 1226 (1995)
  - Active-learning mapping of the Vicsek model phase diagram, arXiv:2604.28167
  - Nature of the order-disorder transition in the Vicsek model (clustering and ordering vs noise), ResearchGate 43020002
  - Wikipedia: Vicsek model
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:network-phase-2: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: p_r_std, p_freq_spread
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) each node is a phase oscillator pulled toward its neighbours' phases with strength 1.0; (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (every step).
- Standard name given: Kuramoto locking within groups and beating between groups on a network that homophilic (type-sorting) rewiring has split into disconnected type-pure groups; an adaptive network that ends in a frozen state
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- Notes: Overturns the first pass's KNOWN-VARIANT hedge (its label already said KNOWN). Node types are frozen, so rewiring just sorts the links and soon stops, leaving about 3 type-pure groups with no links between them. Each group locks internally, which collapses the effective-frequency spread; groups drift at their own average frequencies, so global r beats, which raises r_std. The effects are real against the knock-out but close to trivial: after a short sorting phase the graph is frozen, and the r fluctuation is just beating between disconnected locked groups, not a new collective state. 'Natural frequency spread' must really mean effective frequency. The run-length dependence of the knock-out difference for that feature looks like a measurement-window effect worth checking.
- Citations:
  - Arenas, Diaz-Guilera, Perez-Vicente, PRL 96, 114102 (2006), Synchronization reveals topological scales in complex networks
  - Papadopoulos, Kim, Kurths, Bassett, Chaos 27, 073115 (2017), and its review of adaptive Kuramoto rewiring producing modular synchronized clusters
  - Holme & Newman, PRE 74, 056108 (2006), fragmentation by homophilic rewiring
  - Gross & Blasius, J. R. Soc. Interface 5, 259 (2008), adaptive networks review
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:agents-16: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: a_density_cv
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 1 type (random start). Agents have 4 types (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. 2 field(s) with diffusion [0.02, 0.02]. Rules, applied in order each step: (1) each agent turns up the gradient of field 1; (2) each agent turns left on an even-type cell and right on an odd-type cell.
- Standard name given: Chemotactic (Keller-Segel-type) aggregation of chiral active particles / chemically signalling active rotors; rotation partially arrests chemotactic collapse
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: First pass's verdict (KNOWN) stands but its mechanism is wrong: with 1 cell type, rule 2 is a uniform constant torque (circle swimmers), not quenched random-torque disorder, so Peruani-Aranson does not apply. Rule 2 alone = random baseline (CV 0.8); rule 1 alone = strong chemotactic collapse (CV 6.5, implies agents emit field 1, which is not stated in the rule text); full program = intermediate CV 1.84, i.e. rotation partly arrests collapse, as in Liebchen et al. 2016. Cannot tell arrested-cluster from travelling-wave regime with density CV alone. Agent types have no effect.
- Citations:
  - Liebchen, Cates & Marenduzzo, Pattern formation in chemically interacting active rotors with self-propulsion, Soft Matter 12, 7259 (2016), arXiv:1605.05982
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:agents-17: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: a_step, a_density_cv, a_nn_ratio, a_local_align
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns away from the centre of any type agents within distance 2.0.
- Standard name given: Autochemotactic active (Brownian) particles: chemotactic collapse arrested by short-range repulsion into finite/dynamic clusters (volume-filling Keller-Segel regime)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Rule-1-only gives full collapse, rule-2-only gives an anti-clustered gas, and the full program gives intermediate finite clusters with agents trapped inside them (a_step drops sharply). This matches the collapse / dynamic-clustering / gas regimes of Pohl-Stark and the blow-up prevention of volume-filling Keller-Segel. The weak negative local alignment fits agents pointing into clusters and repelling each other head-on (aster-like inward polarization) and is not a meaningful new feature. The verdict assumes agents deposit chemical; if the field is only seeded, the behaviour is still known (chemotaxis to a fixed attractant with exclusion).
- Citations:
  - Pohl & Stark, Dynamic clustering and chemotactic collapse of self-phoretic active particles, PRL 112, 238303 (2014)
  - Pohl & Stark, Eur. Phys. J. E 38, 93 (2015)
  - Stark, Artificial chemotaxis of self-phoretic active colloids: collective behavior, Acc. Chem. Res. 51, 2681 (2018)
  - Painter & Hillen, Volume-filling and quorum-sensing in models for chemosensitive movement, Can. Appl. Math. Q. 10, 501 (2002)
  - Saha, Golestanian & Ramaswamy, Clusters, asters, and collective oscillations in chemotactic colloids, PRE 89, 062316 (2014)
  - Schweitzer & Schimansky-Geier, Physica A 206, 359 (1994) (not directly verified)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:agents-8: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: a_density_cv, a_step, a_local_align
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns to the average heading of agents within distance 2.0.
- Standard name given: Autochemotactic (Keller-Segel) collapse of self-propelled particles with Vicsek alignment; locally aligned chemotactic clusters
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 4
- Twist (if any): 
- Notes: Sub-programs show clustering (density_cv ~11) and arrest (low a_step) come entirely from rule 1 (autochemotactic collapse; field emission must be implicit since rule-1-only clusters). Rule 2 only raises local alignment (0.34 -> 0.87), which is expected for Vicsek alignment at the high density inside collapsed clusters. Genuine but generic, near-absorbing collapsed state. The first pass's extra claims (rotating clusters, global order 0.39) are not among the verified features and should be disregarded.
- Citations:
  - Marsden, Valeriani, Sullivan, Cates, Marenduzzo, Soft Matter 10, 157 (2014) - CFV model (chemotactic force + Vicsek alignment)
  - Pohl & Stark, PRL 112, 238303 (2014) - chemotactic collapse of self-propelled particles mapped to Keller-Segel
  - Chavanis, stochastic Keller-Segel model, Commun. Nonlinear Sci. Numer. Simul. (2010) - chemotactic collapse
  - Agent-Based Modeling Reveals Possible Mechanisms for Observed Aggregation Cell Behaviors, Biophys. J. (2018) - Vicsek + chemotaxis, alignment around aggregates
  - arXiv:2506.18263 - BGL theory for autochemotaxis with alignment
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:field-1: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: f_std, f_n_blobs, f_peak_sharp
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 1 field(s) with diffusion [0.02]. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 8.0; (2) each agent deposits 0.1 of field 0 where it stands.
- Standard name given: Vicsek model flocking/clustering (polar clusters / travelling bands), read out through a passive one-way deposited diffusing field
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: The field features (more uneven field, fewer blobs, sharper peaks) are a smoothed, time-integrated picture of agent density. The underlying effect is documented Vicsek clustering, banding or flocking. Agents never respond to the field, so it cannot create new behaviour. The clustering is real (robust to seeds, 4x area and 3x run length), but the field itself is only a readout. Ignore the rule-1-only sub-program values: with no depositing the field is nearly empty, so those blob and sharpness numbers are noise. The contrast shrinks somewhat over longer runs, consistent with the field piling up a moving flock's trail. The first pass's claim about which Vicsek regime these parameters sit in was not verified, but every likely regime predicts the same field signature.
- Citations:
  - https://arxiv.org/pdf/1511.01451
  - https://arxiv.org/html/1709.01499
  - https://arxiv.org/pdf/2604.28167
  - https://guava.physics.ucsd.edu/~nigel/Courses/Web%20page%20563/Essays_2017/PDF/Chatterjee.pdf
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:field-3: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: f_std, f_n_blobs
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 2 field(s) with diffusion [0.02, 0.1]. Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.1 per step); (2) type-1 cells release 1.0 of field 0 per step.
- Standard name given: Zero-temperature majority-vote / Glauber-Ising coarsening (curvature-driven, non-conserved), read out by a passive, non-decaying diffusing field
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Only rule 1 involves interaction. The field is a one-way, non-decaying sum over time of the type-1 domain map, so high f_std and few blobs simply mean the domains have coarsened (or frozen into stripes or one phase) compared with a random colouring. f_std grows linearly with run length and never reaches a steady state. Agents and field 1 have no effect. The pattern is genuine coarsening, but the field features add no new phenomenology; the late state is often frozen or near-absorbing.
- Citations:
  - Spirin, Krapivsky, Redner, Fate of zero-temperature Ising ferromagnets, Phys. Rev. E 63, 036118 (2001), arXiv:cond-mat/0008019
  - Coarsening and Bifurcations in Wide-Range Two-Dimensional Totalistic Cellular Automata, arXiv:2604.10082
  - Bray, Theory of phase-ordering kinetics, Adv. Phys. 43, 357 (1994)
  - Superdiffusion-like behavior in zero-temperature coarsening of the d=3 Ising model (states t^1/2 non-conserved growth law)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:network-2: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: n_activity, n_type_modularity
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 2 types (random start). Rules, applied in order each step: (1) nodes play a two-strategy game with their neighbours (T=1.9, S=0.0) and copy the best-scoring neighbour; (2) each node adopts its neighbours' majority type (with probability 0.1 per step).
- Standard name given: Spatial/network Prisoner's Dilemma (Nowak-May best-neighbour imitation, b=1.9) with mixed payoff-driven + conformity (majority-rule) strategy updating; network reciprocity with cooperator clusters and fluctuating boundary strategists
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Rules-1-only gives type modularity 0.158, so 'type' is the game strategy. Modularity = cooperator clustering (network reciprocity). Activity = boundary nodes fluctuating, made larger by payoff and conformity rules disagreeing. Rules-2-only is a frozen state (0/0). The first pass's network-structure features are irrelevant to the passing features, but its KNOWN verdict holds. The difference from the literature (conformity as a per-step chance rather than a fixed fraction of conformist players) is a rule change, not new behaviour.
- Citations:
  - Gomez-Gardenes, Campillo, Floria, Moreno, PRL 98, 108103 (2007)
  - Gracia-Lazaro, Floria, Gomez-Gardenes, Moreno, Chaos Solitons Fractals 56, 188 (2013)
  - Szolnoki & Perc, J. R. Soc. Interface 12, 20141299 (2015)
  - Pena, Pestelacci, Tomassini, Volken, Conformity and Network Effects in the Prisoner's Dilemma; Pena et al. PRE 80, 016110 (2009)
  - Cui & Wu, Physica A 392, 1500 (2013)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:005: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_activity_trend
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.1 per step).
- Standard name given: Spatial snowdrift (hawk-dove) game with imitate-the-best updating plus voter-model noise; coarsening toward near-consensus
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: S>P makes this a snowdrift game with r=0.5 and additive payoffs (T-R=S-P). The game alone already gives one dominant type (0.815); the voter step strengthens coarsening (0.94). Activity falls, so the system is drifting toward a near-absorbing consensus. The knock-out removes both rules, leaving a frozen random colouring, so the big difference from it is unsurprising. Result is robust, not a numerical artifact, but nothing new.
- Citations:
  - Hauert & Doebeli 2004, Nature 428:643 (doi:10.1038/nature02360)
  - Szabo & Fath 2007, Evolutionary games on graphs, Phys. Rep. (arXiv cond-mat/0607344)
  - Cox, Durrett & Perkins voter model perturbations; Durrett, Spatial evolutionary games with small selection coefficients (arXiv 1406.5876)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:013: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_type_entropy
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) each cell adopts its neighbourhood's majority type (with probability 0.1 per step).
- Standard name given: Spatial Prisoner's Dilemma (Nowak-May imitate-the-best, Moore + self) with conformity/majority-rule updating (Szolnoki & Perc 2015; mixed payoff-conformity imitation)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Rule 1 alone: defection takes over (entropy 0.05). Rule 2 alone: majority coarsening. Together: conformity-stabilised cooperator/defector coexistence in compact clusters with persistent flipping, which is the documented 'conformity enhances network reciprocity' effect. The knock-out baseline is a frozen random lattice, so the comparison with the sub-programs is the meaningful evidence; that comparison supports a genuine emergent effect. Payoff point (S=-0.5), synchronous updating and per-step probabilistic conformity are rule differences, not new phenomenology. Sources were checked through search snippets only, not full texts.
- Citations:
  - Szolnoki A, Perc M. Conformity enhances network reciprocity in evolutionary social dilemmas. J R Soc Interface 12:20141299 (2015)
  - Szolnoki A, Perc M. Leaders should not be conformists in evolutionary social dilemmas. Sci Rep (2016)
  - Javarone MA, Antonioni A, Caravelli F. Conformity-driven agents support ordered phases in the spatial public goods game. EPL 114:38001 (2016)
  - Szolnoki A, Chen X. Competition and partnership between conformity and payoff-based imitations in social dilemmas
  - Nowak MA, May RM. Evolutionary games and spatial chaos. Nature (1992)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:015: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: n_oscillation, n_edge_turnover
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 2 types (random start). Rules, applied in order each step: (1) nodes play a two-strategy game with their neighbours (T=1.9, S=0.0) and copy the best-scoring neighbour; (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (with probability 0.003 per step).
- Standard name given: Nowak-May spatial Prisoner's Dilemma chaos (b=1.9, imitate-the-best) on a random graph with independent fixed-type homophilic rewiring (adaptive-network PD, Zimmermann-Eguiluz-San Miguel family)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: n_oscillation is Nowak-May chaotic C/D coexistence at b in (1.8,2); the game alone already gives 0.13. Rewiring raises it to 0.33 by repeatedly disturbing frozen or short-period patterns, which is an expected effect of network churn and a difference in amount only (possibly a long transient until cross-type links run out). n_edge_turnover is essentially a direct readout of rule 2 (0.38 with rule 2 alone vs 0.34 full), so that feature is near-trivial; the knock-out baseline of 0 only reflects that the rewiring rule reads neighbour type.
- Citations:
  - Nowak & May 1992, Nature 359:826 'Evolutionary games and spatial chaos'
  - Nowak, Bonhoeffer & May 1994, IJBC 'More spatial games' (random spatial arrays)
  - Zimmermann, Eguiluz & San Miguel, PRE 72, 056118 (2005); PRE 69, 065102(R) (2004)
  - Szabo & Fath 2007, Phys Rep 'Evolutionary games on graphs' (arXiv cond-mat/0607344)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:027: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: a_density_cv, a_local_align, a_step, a_polar_std
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.02]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns to the average heading of agents within distance 1.0.
- Standard name given: Chemotactic (Keller-Segel-type) aggregation/collapse combined with Vicsek polar alignment -- localized aligned clusters / active droplets in communicating active matter
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Gradient-climbing alone gives the clumping and the arrest (CV 8.2, step 0.24); alignment alone gives a moving Vicsek flock; together they give highly aligned, nearly arrested clumps. That is alignment on top of chemotactic collapse, which is documented. The rule text does not say what emits the field: if agents are not the source, this is trapping at the peaks of the initial field, which is even more trivially known. Low a_step cannot tell vortices from jammed aligned clumps, but both are documented. Nonzero polar_std and robust size/time scaling argue against a frozen artifact.
- Citations:
  - Ziepke, Maryshev, Aranson, Frey, Nat. Commun. 13, 6727 (2022) https://www.nature.com/articles/s41467-022-34484-2
  - Pohl & Stark, PRL 112, 238303 (2014) Dynamic clustering and chemotactic collapse of self-phoretic active particles
  - Collective chemotactic search strategies, arXiv:2409.04262
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:033: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: a_density_cv, a_nn_ratio, a_local_align
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 1 type (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns left on an even-type cell and right on an odd-type cell.
- Standard name given: Autochemotactic chiral active particles (chemically interacting active rotors): chemotactic clustering arrested/controlled by active rotation, with chemotaxis-induced local locking of propulsion direction
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Since there is only one cell type, rule 2 is a constant turn in one direction, so every agent is a noisy circle swimmer. Rule 1 alone collapses almost completely, which implies the agents emit field 0 (autochemotaxis), but no emission rule is listed, so check this. The verdict KNOWN is upheld. The first pass's 'mostly trivial' framing is overturned: local alignment appears only with both rules, an emergent effect matching the Liebchen et al. locking mechanism. I did not determine which regime (coarsening clusters vs travelling wave) these parameters fall into.
- Citations:
  - Liebchen, Cates & Marenduzzo, Pattern formation in chemically interacting active rotors with self-propulsion, Soft Matter 12, 7259 (2016), arXiv:1605.05982
  - van Teeffelen & Löwen, Dynamics of a Brownian circle swimmer, PRE 78, 020101 (2008)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:036: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: f_std, f_n_blobs
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.02]. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 3.0; (2) each agent deposits 0.1 of field 0 where it stands.
- Standard name given: Vicsek model flocking/clustering (ordered or coexistence phase, low-noise macrocluster), shown through a passive deposited field with no decay
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: The field is passive: it never acts back on the agents and does not decay. The f_std and f_n_blobs gains over the knock-out just show that alignment-driven clustering concentrates the deposits into one patch, versus near-uniform deposits from independent walkers. The only twist is geometric: a straight-moving flock on the torus keeps depositing on the same track. That follows from known behaviour and is not new. The features are genuine signal, but they come from known Vicsek clustering; f_std also partly grows simply because deposits accumulate with no decay.
- Citations:
  - Vicsek et al., PRL 75, 1226 (1995)
  - Grégoire & Chaté, Onset of collective and cohesive motion, PRL 92, 025702 (2004)
  - Emergence of a single cluster in Vicsek's model at very low noise, arXiv:1709.01499
  - Phase Diagram of the Two-Dimensional Vicsek Model, arXiv:2610.26569
  - Collective dynamics of trail-interacting particles, arXiv:2512.03950 (contrast: feedback trail models)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:040: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour; (2) a type-1 cell becomes type 0 if at least 2 neighbours are type 0 (with probability 0.03 per step).
- Standard name given: Nowak–May spatial Prisoner's Dilemma (imitate-the-best, Moore + self-interaction), low-b (b=1.3) cooperator-dominated quasi-static regime with blinking defector fragments
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Rule 1 alone reproduces every feature; rule 2 is a one-way D->C bias that only raises cooperator dominance (alone it drives the lattice to all-C with no structure). Flipping cells are the documented blinkers along defector fragments. The pattern is genuine game-driven clustering, but the state is nearly frozen; the knock-out baseline is an uncorrelated random lattice, which makes the contrasts look bigger than they are.
- Citations:
  - Nowak & May 1992, Nature 359:826 'Evolutionary games and spatial chaos'
  - Nowak & May 1993, Int. J. Bifurc. Chaos 3:35
  - Szabó & Fáth, 'Evolutionary games on graphs', arXiv:cond-mat/0607344 (Sec. 6.5)
  - Wu, Fu & Wang, 'Evolutionary games and spatial periodicity', arXiv:2209.08267
  - Mukherji et al. 1996 (noise extension of Nowak–May)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:044: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell copies the type of a random neighbour (with probability 0.3 per step); (2) each cell adopts its neighbourhood's majority type (with probability 0.2 per step); (3) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.
- Standard name given: Spatial snowdrift (hawk-dove) game with imitate-the-best on a Moore lattice, plus voter/majority noise: one dominant strategy with many small persistent minority clusters
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Verdict matches the first pass, but its reasoning is corrected: g_n_domains looks log-scaled (ln2 and ln4 appear among the sub-program values), so the full program has about 45 domains, not 3.4. The state is a stationary dominant phase (~90%) with small scattered minority clusters (corrlen ~2), set mainly by the snowdrift payoffs (T>R>S>P). The rule-3-only sub-program already shows it, and voter/majority noise only adjusts it. This is not coarsening toward consensus. The rules-2,3 sub-program reaching a uniform absorbing state is the only artifact-like outcome.
- Citations:
  - Hauert & Doebeli 2004, Nature 428:643 (https://www.nature.com/articles/nature02360)
  - Fu, Nowak, Hauert 2010 JTB, invasion and expansion of cooperators in lattice populations (https://personal.math.ubc.ca/~hauert/publications/reprints/fu_jtb10.pdf)
  - Szabo & Fath, Evolutionary games on graphs (arXiv cond-mat/0607344)
  - Wang/Perc et al., arXiv 1102.4570 (summary of snowdrift filament clusters)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:045: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_n_domains, g_corrlen, g_n_types, g_reversals, g_type_entropy, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) each cell copies the type of a random neighbour (with probability 0.1 per step); (2) a type-0 cell becomes type 2 where field 0 exceeds 0.5.
- Standard name given: Multi-state (Potts-like) voter model coarsening without surface tension on a 2D lattice; field rule acts as a one-shot elimination of type 0 into type 2
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Spatial features (moran, corrlen, n_domains, reversals) are the same in the rule-1-only and full programs, so they come from voter dynamics alone. g_n_types=3 and g_type_entropy~0.74 come from rule 2 converting nearly all type-0 cells to type 2 near the start (expected normalised entropy for proportions 1/2,1/4,1/4 over 4 types is 0.75). The coarsening is genuine but standard. The type-count and entropy changes are artifacts, not emergent, and there is no real synergy between the two rules. Agrees with the first pass, but the field rule is not 'probably irrelevant': it acts once, removes a type, and then does nothing.
- Citations:
  - Dornic, Chate, Chave, Hinrichsen, Phys. Rev. Lett. 87, 045701 (2001), arXiv:cond-mat/0101202
  - S. Redner, Reality-inspired voter models: A mini-review, C. R. Physique 20, 275 (2019)
  - arXiv:1412.6297 (Coexistence in neutral theories), voter-model criticality and logarithmic coarsening
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:046: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: f_std, f_n_blobs, f_peak_sharp
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) each cell copies the type of a random neighbour (with probability 0.2 per step); (2) type-0 cells release 1.0 of field 1 per step.
- Standard name given: Multi-state 2D voter model coarsening (tensionless domain growth) with a one-way, passively driven, non-decaying diffusive field that records type-0 domains
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: The cell layer is a pure voter model because the field never feeds back on the cells. f_n_blobs comes from voter domains (rule 1 alone gives about 1.1). f_std and f_peak_sharp come from emission piling up under persistent type-0 domains. With no decay the field grows without limit (f_std difference rises from 62.8 to 141 at 3x run length), so the field statistics depend on run time and are a passive record, not a pattern the field forms itself. The emergent pattern is genuine but textbook; the field is a slaved readout, not a numerical artifact.
- Citations:
  - Dornic, Chaté, Chave & Hinrichsen, Critical coarsening without surface tension: the universality class of the voter model, PRL 87, 045701 (2001), arXiv:cond-mat/0101202
  - Holley & Liggett, Ann. Probab. 3, 643 (1975) (seen in reference lists, not read)
  - Phase transition and power-law coarsening in Ising-doped voter model, arXiv:1706.03575
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:057: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 3 types (random start). Agents have 1 type (random start). Agents move at speed 0.5 with heading noise 0.0 rad per step. Rules, applied in order each step: (1) each agent turns away from the centre of any type agents within distance 1.0; (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.4 moves to a random empty cell; (3) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Noise-assisted coarsening in the Schelling model with vacancies (escape from frozen finite-cluster states; Schelling–Ising at finite temperature)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: The verdict label stands but the first pass's mechanism is reversed. The data show the agent forcing promotes order (largest domain 0.03 -> 0.47 when rule 3 is added to rule 2) by unfreezing zero-noise Schelling; it does not push the system into the mixed phase. g_reversals is largely trivial direct forcing (rule 3 alone gives 0.31). Agent repulsion (rule 1) only spreads the noise more evenly. Unverified twist: the cyclic 0->1->2->0 forcing is asymmetric and may bias toward one type or toward empty cells; check which type forms the largest domain.
- Citations:
  - Dall'Asta, Castellano & Marsili, J. Stat. Mech. L07002 (2008), arXiv:0707.1681
  - Stauffer & Solomon, Eur. Phys. J. B 57, 473 (2007), arXiv:physics/0701051
  - Müller, Schulze & Stauffer, arXiv:0706.2592
  - Gauvin, Vannimenus & Nadal, Phase diagram of a Schelling segregation model, arXiv:0903.4694
  - Ódor, arXiv:0710.1496
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:060: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity_trend, g_corrlen, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 3.14159 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent advances the type of the cell under it by one (cyclically); (3) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.
- Standard name given: Nowak–May spatial game (imitate-the-best, Moore neighbourhood + self) with snowdrift payoffs, plus agent-driven local strategy flips (noise)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Rule 3 alone already gives the observed signature (largest domain ~0.82, corrlen 2, reversals, ~4.8 domains); the full program differs only modestly. Correction to the first pass: the field-following rule is not inert. It changes the activity trend (rules 1+2 vs rule 2: -3.53 vs 0) and stops agent flipping from disrupting the lattice (largest domain 0.86 vs 0.34 without rule 1). Likely cause (inferred, not verified in code): agents gather at peaks of an unreplenished starting field, confining flips to a few spots. This is a mechanism detail, not new phenomenology. The knock-out baseline is a random colouring, which inflates the apparent interaction effect; the activity feature is a decline, not a rise.
- Citations:
  - Nowak & May 1992, Nature 359:826 'Evolutionary games and spatial chaos'
  - Hauert & Doebeli 2004, Nature 428:643 'Spatial structure often inhibits the evolution of cooperation in the snowdrift game'
  - Sandholm et al., Agent-Based Evolutionary Game Dynamics, ch. III-1 (Nowak–May replication)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:071: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_activity_trend, g_interface
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) field 1 grows autocatalytically by consuming field 0 (rate 0.5 x f x g^2); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.2 moves to a random empty cell; (3) type-2 cells release 1.0 of field 1 per step.
- Standard name given: Schelling segregation model (multi-group, vacancies, random relocation, low tolerance threshold) with a passive one-way-driven reaction-diffusion field
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Cells never read either field, so the cell dynamics are exactly Schelling relocation. The g_interface drop comes entirely from rule 2 (0.603 vs 0.615, within ~0.5 robust spread). The extra activity trend in the full program (-1.50 vs -1.08 for rule 2 alone) most likely comes from field changes driven by type-2 emitters moving (rules 1+3 without movement give 0). That is a measurement contribution from an inert slave field, not new behaviour. Segregation is weak because a 0.2 threshold is below the random same-type share (~1/3 with 3 occupied types). The segregation itself is genuine; the field coupling is trivial.
- Citations:
  - https://arxiv.org/pdf/0903.4694
  - https://www.nature.com/articles/s41598-022-23224-7
  - https://arxiv.org/pdf/2509.14462
  - https://www.nature.com/articles/s41598-025-26995-x
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:008: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_corrlen, g_n_domains, g_moran, g_cyclic_changes
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step).
- Standard name given: Cyclic cellular automaton (Griffeath CCA), rule '313' (R1/T3/C3/NM), with alpha-asynchronous / stochastic updating (discrete-time approximation of the threshold cyclic particle system): self-organised cycling wave/spiral domains
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- Notes: This exact rule is Griffeath's named 313 CCA. The four measured features (correlation length ~5, fewer domains than the knock-out, Moran's I ~0.79, cycling well above chance) are what the documented self-organising, locally periodic CCA regime looks like. Cycling rules out a frozen state, and the knock-out baseline (Moran ~0, correlation length 1) shows the pattern comes from the interaction. Updating each cell with probability 0.1 is the standard alpha-asynchronous scheme and roughly approximates the cyclic particle systems that the CCA was introduced as a discrete-time analogue of. Caveat: these summary statistics on a 64x64 grid can't tell crisp spirals from patchy travelling-wave domains, and any effect of the random updating on spiral shape isn't measured here, so no twist in the behaviour is shown.
- Citations:
  - http://psoup.math.wisc.edu/mcell/rullex_cycl.html
  - https://github.com/setanarut/ccae
  - https://arxiv.org/pdf/nlin/0507067
  - https://math.ucdavis.edu/~gravner/papers/cca.pdf
  - https://link.springer.com/chapter/10.1007/978-1-4612-0451-0_8
  - https://arxiv.org/pdf/1610.07320
  - https://projecteuclid.org/journals/experimental-mathematics/volume-2/issue-3/Asymptotic-behavior-of-excitable-cellular-automata/em/1062620830.full
  - https://arxiv.org/pdf/1911.03609
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:009: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: p_r, p_freq_spread
- Rule: N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=0.1  (compact grammar notation)
- Standard name given: Partial (frequency-entrained) Kuramoto synchronization on a stochastically blinking / fast-switching sparse coupling network
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- Notes: r~0.45 vs 0.07 null (~1/sqrt(N)); effective-frequency spread collapses 0.51->0.025 = entrainment. r gain shrinks at 4x area (finite-size/local-coupling effect) and grows at 3x run length (slow ordering, not yet stationary, not frozen). Feature label 'natural frequencies' is a misnomer; the metric must be effective frequencies. First-pass r~0.65 differs from the re-measured median 0.45; this does not affect the verdict. If kN=1 means lattice nearest-neighbour coupling, the behaviour is still documented (lattice Kuramoto, Sakaguchi-Shinomoto-Kuramoto lineage).
- Citations:
  - Acebrón et al., Rev. Mod. Phys. 77, 137 (2005)
  - Petit, Carletti, Asllani, Fanelli et al. / 'Synchronization in time-varying random networks with vanishing connectivity', Sci. Rep. (2019), arXiv:1811.09591
  - Barabash & Belykh, 'Synchronization Thresholds in an Ensemble of Kuramoto Phase Oscillators with Randomly Blinking Couplings', Radiophys. Quantum Electron. (2018)
  - Ghosh et al., 'The synchronized dynamics of time-varying networks', arXiv:2109.07618
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:012: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_interface, g_activity
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (with probability 0.3 per step).
- Standard name given: Stochastic (asynchronous / probabilistic-update) cyclic cellular automaton = cyclic particle system, kappa=4, threshold 1, Moore range 1; fluctuating/clustering 'forest-fire' stochastic wave-front regime
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Probabilistic firing is the discrete-time version of the asynchronous cyclic particle system that the CCA was designed to imitate, so it is not a twist. kappa=4 is on the fluctuating side (fixation needs >=5). The knock-out is a frozen random colouring (unlike-neighbour fraction about 0.75 = 1-1/4), so the interface effect is modest real coarsening, not an artifact. g_activity=0.417 is above the 0.3 per-step cap, which suggests the metric is scaled or summed differently; the difference from the knock-out is still valid.
- Citations:
  - Fisch, Gravner & Griffeath (1991), Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39
  - Fisch, Gravner & Griffeath (1991), Cyclic Cellular Automata in Two Dimensions, Spatial Stochastic Processes, Birkhauser
  - Bramson & Griffeath (1989), Flux and fixation in cyclic particle systems, Ann. Probab. 17:26-45
  - Fisch (1990), The one-dimensional cyclic cellular automaton, J. Theor. Probab. 3:311-338
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:016: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). 1 field(s) with diffusion [0.2]. Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.9, S=0.0) and copy the best-scoring neighbour.
- Standard name given: Nowak-May spatial Prisoner's Dilemma (synchronous imitate-the-best, Moore neighbourhood + self-interaction), chaotic regime 1.8<b<2
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: b=T=1.9 is inside the chaotic window. Reversals and fragmenting domains are the documented spatial chaos. Largest domain of about 0.68 matches the defector background implied by the known cooperator fraction of about 0.318 (12 ln2 - 8). The diffusing field is not coupled to anything, so it is inert. The knock-out removes the only rule, freezing the random initial state, so the comparison is trivially passed. Correlation length 2 means short-range correlation, not long-range. Small tori may end in periodic cycles, which is also documented.
- Citations:
  - Nowak & May, Nature 359:826-829 (1992), doi:10.1038/359826a0
  - Nowak & May, Int. J. Bifurcation Chaos 3:35-78 (1993)
  - Szabo & Fath, Evolutionary games on graphs, Phys. Rep. 446 (2007), arXiv:cond-mat/0607344 sec. 6.5
  - Evolutionary games and spatial periodicity, ScienceDirect S294985542300014X
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:017: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_reversals
- Rule: C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.5, 0.5)) p=0.1  (compact grammar notation)
- Standard name given: Nowak–May spatial evolutionary game (unconditional imitation / imitate-the-best) with partially asynchronous (p=0.1) updating; boundary strategy oscillation/flipping in C–D coexistence
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- Notes: Changed from KNOWN-VARIANT to KNOWN because there is no twist in the behaviour. The interaction knock-out is trivially frozen, since imitation is the only rule, so a reversal difference of 1.0 vs 0.0 only shows the dynamics are not frozen. A median of 1.0 with spread 0.02 suggests the measure is saturated or binary. The flipping is real but low-information and generic to this model class; it is not an artifact.
- Citations:
  - Nowak & May 1992, Nature 359:826 (Evolutionary games and spatial chaos)
  - Szabó & Fáth, Evolutionary games on graphs, arXiv cond-mat/0607344
  - Roca, Cuesta & Sánchez, arXiv 0911.1720 (unconditional imitation for PD and Snowdrift)
  - Scale-invariant behavior in a spatial game of prisoners' dilemma (periodic patterns)
  - arXiv 1408.1793 (asynchronous updating, strategy flippings)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:018: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: n_edge_turnover, n_type_modularity, n_giant
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) each node copies the type of a random neighbour (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (every step).
- Standard name given: Fragmentation transition in the multistate (G=3) adaptive/coevolving voter model with rewire-to-same (Holme-Newman model, rewiring-dominated regime)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Observed end state is the known absorbing fragmented state: about 3 like-type components (modularity 0.644 vs ideal 0.667 for 3 equal modules; giant about 0.40). Edge turnover looks like a running total of rewired original edges (rewiring-only run gives 0.6675, about 2/3 = starting share of mixed-type edges for 3 random types), not ongoing activity. Genuine interaction-driven pattern but near-trivial: rule 2 alone (homophilic sorting) reproduces nearly all of it; voter copying only slightly changes the final values. Not a discovery.
- Citations:
  - Holme & Newman, Phys. Rev. E 74, 056108 (2006), arXiv:physics/0603023
  - Böhme & Gross, Fragmentation transitions in multistate voter models, Phys. Rev. E 85, 066117 (2012)
  - Durrett et al. / Shi et al., Multi-opinion coevolving voter model with infinitely many phase transitions, arXiv:1303.7434 (2013)
  - Opinion dynamics: Statistical physics and beyond, arXiv:2507.11521
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:019: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: n_edge_turnover, n_consensus_gain
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) each node copies the type of a random neighbour (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (with probability 0.1 per step).
- Standard name given: Adaptive / coevolving (multistate) voter model with rewire-to-same, consensus (non-fragmented) phase
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Low rewiring (0.1) vs copying (1.0) puts the system in the consensus phase, below fragmentation. Consensus gain is essentially all from voter copying (rules-1-only 0.625 vs full 0.635, a gap of about one robust spread). Edge turnover is the direct output of rule 2, raised from 0.66 to 0.97 because copying keeps creating discordant links. Real model behaviour, not a numerical artifact, but the coupling adds almost nothing beyond each rule's direct effect. Synchronous updating and G=3 are minor differences with no qualitative effect here.
- Citations:
  - Holme & Newman, Phys. Rev. E 74, 056108 (2006)
  - Vazquez, Eguiluz & San Miguel, Phys. Rev. Lett. 100, 108702 (2008)
  - Böhme & Gross, Phys. Rev. E 83, 035101 (2011)
  - Böhme & Gross, Phys. Rev. E 85, 066117 (2012) - multistate incl. 3 opinions
  - Shi et al. (2013) multistate coevolving voter model
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:027: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) neighbouring cells swap contents (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.
- Standard name given: Schelling segregation model with vacancies and thermal/random-exchange noise (noisy Schelling / Schelling-Ising model)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Type 0 = vacancy, so this is the standard two-species Schelling model. The random swaps act as infinite-temperature exchange noise. Moran's I is higher with noise (0.36) than with the movement rule alone (0.19); this matches the documented finding that small nonzero noise lets the system escape frozen finite clusters and segregate further (Stauffer-Solomon; Gauvin et al.). Clustering is genuine but moderate, and the domain-count effect is small. g_reversals is largely a direct, trivial result of blind swapping (0.59 with the swap rule alone), not an emergent oscillation.
- Citations:
  - Stauffer & Solomon, Urban and Scientific Segregation: The Schelling-Ising Model, arXiv:0710.5237
  - Gauvin, Vannimenus & Nadal, Phase diagram of a Schelling segregation model, arXiv:0903.4694
  - Dall'Asta, Castellano & Marsili, Statistical physics of the Schelling model of segregation, J. Stat. Mech. (2008) L07002
  - Aging effects in Schelling segregation model, Sci. Rep. 2022 (s41598-022-23224-7)
  - Gauvin, Nadal & Vannimenus et al., A Schelling model with switching agents, arXiv:1212.4929
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:028: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_spec_peak
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell; (2) each cell copies the type of a random neighbour (every step).
- Standard name given: Synchronous multi-state (3-state, vacancy-as-state) voter model coarsening with Schelling-move noise; cf. Schelling-Voter model
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- Notes: Rule 2 (unconditional synchronous copy of a random neighbour, including empty cells) dominates. The voter-only sub-program gives most of the signal (domains 6.45, reversals 0.65, spec_peak 0.27); the move rule adds disorder (more reversals, fewer domains). The knock-out freezes the lattice, so the reversal and spectral-peak differences are guaranteed for any stochastic rule and are not evidence of a pattern. Domain persistence = slow 2D voter coarsening (logarithmic for two states, possibly algebraic for multi-state); reversals = noisy interfaces without surface tension. The spectral peak is unverified as a real oscillation and is most likely a synchronous-update or interface-noise effect; peak location not checked. Coarsening is genuine but known; two of three features are near-artifacts, so trivial_or_artifact is borderline.
- Citations:
  - https://arxiv.org/pdf/1506.05321
  - https://arxiv.org/html/1903.09198
  - https://www.researchgate.net/publication/11879762_Critical_Coarsening_without_Surface_Tension_The_Universality_Class_of_the_Voter_Model
  - https://www.sciencedirect.com/science/article/abs/pii/S0378437117306118
  - https://www.sciencedirect.com/science/article/abs/pii/S0960077913001720
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:034: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_n_domains, g_corrlen, g_activity_trend, g_moran
- Rule: C kC=4 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0  (compact grammar notation)
- Standard name given: 4-colour cyclic cellular automaton / 4-species cyclic dominance clustering with neutral-pair (difference-2) blockades; noisy-copy (voter-like) variant
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Same verdict as the first pass, but its 'ongoing coarsening' description is not supported: the correlation length stays near 3 and its gain over the knock-out halves (2.0 to 1.0) at 3x run length, while activity declines. The data fit short-range clustering that is settling or freezing behind non-interacting 0|2 and 1|3 borders, which is also documented. The pattern is genuine but weak and short-range; a partly frozen absorbing state should be checked (look at the final activity level and whether any adjacent next-colour pairs remain).
- Citations:
  - https://www.researchgate.net/publication/2895807_Quantifying_Self-Organization_in_Cyclic_Cellular_Automata
  - https://arxiv.org/pdf/1711.04741
  - https://arxiv.org/html/cond-mat/9606069
  - https://arxiv.org/abs/q-bio/0310017
  - https://arxiv.org/html/nlin/0507067
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:037: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_type_entropy
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.1 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour.
- Standard name given: Spatial Prisoner's Dilemma (Nowak-May imitate-the-best, Moore+self) with conformity / majority-rule mixing (conformity-enhanced network reciprocity)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- Notes: Behaviour = clusters of both types coexisting in a mixed state sustained by competition between conformity (smooths borders, coarsens) and payoff-driven imitation (roughens borders, nearly takes over the grid on its own at T=1.5,S=-0.5). Rule differences from the literature (synchronous update, every cell conforms with p=0.1 instead of a fixed conformist subset) are not differences in behaviour. Caveats: the knock-out baseline is a static random field, so the four features mostly show the pattern is not random noise; g_reversals already appears with the majority rule alone (its own randomness), so it is not an interaction effect. The real interaction effect is that the full system sits outside both single-rule versions. Unchecked: which type dominates, and whether conformity here preserves coexistence rather than boosting cooperation as Szolnoki-Perc report; a possible small twist, not a discovery.
- Citations:
  - Nowak & May 1992, Nature 359:826 'Evolutionary games and spatial chaos'
  - Szolnoki & Perc 2015, J. R. Soc. Interface 12:20141299 'Conformity enhances network reciprocity in evolutionary social dilemmas'
  - Szolnoki & Perc 2016, Sci. Rep. 6:23633 'Leaders should not be conformists in evolutionary social dilemmas'
  - Sandholm et al. textbook chapter III-1 'Spatial chaos in the Prisoner's Dilemma' (wisc.pb.unizin.org)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:041: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran
- Rule: C kC=2 | C[any] CNT_GE(b=0,th=3) -> IMITATE_BEST(g=(1.5, 0.5)) p=1.0  (compact grammar notation)
- Standard name given: Spatial 2x2 evolutionary game with imitate-the-best (Nowak-May cellular automaton), non-chaotic regime: one strategy invades and the lattice freezes with small leftover islands
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Agrees with the first pass. Fingerprint: about 89% one state, about 4 domains, correlation length 2, a single reversal per cell, Moran's I about 0.49. That is ordering followed by freezing, not chaos. The neighbour-count gate (at least 3 neighbours in state 0) changes the rules, not the behaviour; it only adds more frozen configurations. If g=(T,S)=(1.5,0.5) with R=1, P=0, the game is Snowdrift rather than PD; both are covered by the imitate-the-best payoff-plane maps. The pattern is a real interaction effect (stable against knock-out and scaling) but ends in a nearly absorbing frozen state, so it is uninteresting rather than an artifact. Confidence is medium because only one search was run and the meaning of g is unconfirmed.
- Citations:
  - Nowak & May 1992, Nature 359:826-829 (Evolutionary games and spatial chaos)
  - https://wisc.pb.unizin.org/agent-based-evolutionary-game-dynamics/chapter/iii-1/
  - https://github.com/giakoumoglou/game_theory_spatial_PD
  - https://arxiv.org/html/2607.24638
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:045: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_n_domains, g_moran, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) each cell copies the type of a random neighbour (every step); (2) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (every step).
- Standard name given: Cyclic cellular automaton (Griffeath '313', R1/T3/C3/NM) with voter-type neighbour copying, i.e. a cyclically dominated voter / rock-paper-scissors excitable medium with mixing
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- Notes: Mild overturn: first pass said KNOWN-VARIANT leaning KNOWN; no phenomenological twist found, so KNOWN. Rule 2 alone (313 CCA) gives most of the features (Moran 0.75, ~3.35 domains, corrlen 3, reversals 0.45). Adding voter copying modestly raises corrlen (3->5) and Moran (0.75->0.82), matching the documented effect that mixing/smoothing enlarges cyclic-dominance patterns. The all-rules-removed knock-out is a frozen random grid, so it is a weak comparison; sub-programs are the useful baseline. 'Reversals' likely measures forward cyclic state changes, not true back-and-forth flips. Activity persists at 3x run length, so not an absorbing state.
- Citations:
  - http://psoup.math.wisc.edu/mcell/rullex_cycl.html
  - https://arxiv.org/pdf/nlin/0507067
  - https://arxiv.org/pdf/cond-mat/0110255
  - https://epub.ub.uni-muenchen.de/19038/1/Frey_19038.pdf
  - https://arxiv.org/html/1608.08932
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:047: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_largest_domain, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.
- Standard name given: 3-colour threshold-1 cyclic cellular automaton (Griffeath CCA / rock-paper-scissors with empty as a species) with long-range random mixing: mobility-driven loss of diversity / global takeover
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Type 0 (empty) is inside the cycle, so this is a 3-species cyclic automaton, not Schelling segregation. Rule 1 alone gives the known threshold-1 CCA cycling texture (no reversals, tiny domains). Random relocation acts as strong mixing and produces the documented mobility-driven takeover (largest domain ~0.81). g_reversals is mostly bookkeeping from relocation (rule 2 alone gives 0.41; cyclic updates cannot produce A->B->A flips), so that feature is weak or partly an artifact. The dominant domain may be a phase snapshot of a near-synchronised global oscillation or a transient towards the uniform dead-end state, not a stable new pattern.
- Citations:
  - Fisch, Gravner, Griffeath (1991) Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39
  - Fisch, Gravner, Griffeath (1991) Cyclic cellular automata in two dimensions, Birkhauser
  - Shalizi et al., Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067
  - Reichenbach, Mobilia, Frey (2007) Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, Nature 448:1046
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:049: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Cyclic cellular automaton '313' (R1/T3/C3/NM, Griffeath), randomly updated, with unconstrained Schelling relocation acting as long-range mixing noise
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: The rule-1-only sub-program has MORE order (Moran 0.77, correlation length 4) than the full program (0.62, 3), so the clustering comes from the cyclic automaton and Schelling relocation only blurs it. The knock-out baseline is the random starting grid, so beating it is a weak test. A correlation length of 3 on a 64x64 grid is short-range patchiness, not long-range order. The pattern is real (not an artifact), but the measurements cannot tell spirals from turbulent or fixating patches.
- Citations:
  - http://psoup.math.wisc.edu/mcell/rullex_cycl.html
  - https://arxiv.org/html/nlin/0507067
  - http://psoup.math.wisc.edu/java/jcca.html
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9653388/
  - https://arxiv.org/pdf/0707.1681
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:052: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: n_consensus_gain, n_type_entropy, n_edge_turnover
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 4 types (random start). Rules, applied in order each step: (1) each node adopts its neighbours' majority type (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (with probability 0.03 per step).
- Standard name given: Majority-rule (majority-preference) opinion dynamics with homophilic rewiring on an adaptive random network (coevolving voter/majority model); slow-rewiring regime leading to consensus
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Sub-program results show consensus gain and low entropy come from majority rule alone (rule 1: 0.685/0.131 vs full 0.6725/0.1175). Edge turnover is just rule 2's rewiring, damped once types settle. None of the three features needs the two rules to interact. Not an artifact (real ordering vs knock-out) but standard majority-vote ordering. Final state looks near-settled with a few possibly stuck minority nodes. The 4 types, all-at-once updates and random-node rewiring don't change the behaviour qualitatively. The first-pass 'giant component ~90%' claim isn't among the measured features, so it's unchecked. Fu & Wang full text not read, abstract only.
- Citations:
  - Fu & Wang, Phys. Rev. E 78, 016104 (2008)
  - Benczik, Benczik, Schmittmann & Zia, Phys. Rev. E 79, 046104 (2009)
  - Vazquez et al. / 'Opinion Dynamics on Coevolving Networks' review (relative timescales decide consensus vs fragmentation)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:054: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity_trend, g_corrlen, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.003 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour.
- Standard name given: Nowak–May spatial Prisoner's Dilemma (deterministic imitate-the-best, Moore 8+1 neighbourhood, weak PD) in the low-b static/periodic regime, with weak majority-rule noise
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: At b=1.3 the game sits below the chaotic regime: defector webs and cooperator clusters freeze, with small period-2 blinkers (the 'reversals'). Majority flips at p=0.003 are a sparse perturbation that briefly restarts local activity. The falling activity trend means relaxation, not growth. A correlation length of 2 cells is short, not long. Fewer domains than the knock-out just means imitation clusters cells compared with a frozen random field. The pattern is real but weak and close to a frozen/periodic state; no new phenomenology.
- Citations:
  - Nowak & May, Nature 359, 826 (1992)
  - Nowak & May, Int. J. Bifurcation Chaos 3, 35 (1993)
  - Szabó & Fáth, Evolutionary games on graphs, Phys. Rep. 446, 97 (2007), arXiv:cond-mat/0607344, Sec. 6.5
  - Mean-field interactions in evolutionary spatial games, arXiv:2107.11088 (static D webs below b≈1.53)
  - Dynamic phase transition in PD on a lattice with stochastic modifications, arXiv:0910.0955 (noise variants, Mukherji et al.)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:056: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_type_entropy
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) neighbouring cells swap contents (with probability 0.05 per step).
- Standard name given: Spatial snowdrift (Hawk-Dove) game with imitate-best updating on a Moore lattice, plus random mobility/swaps (Hauert & Doebeli 2004; Vainstein et al. 2007; Sicardi et al. 2009)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Payoffs give snowdrift with cost-to-benefit ratio 0.5. Observed: one type takes ~87% of the lattice, the other survives as a few small patches, correlation length 2, some cells keep flipping; the game rule alone produces nearly all of this. Random swaps only shift the numbers (majority 0.81 to 0.87), and random mobility in spatial snowdrift games is well documented, though which way it shifts cooperation depends on the setup. The reversals feature is weak (binary, also present with swaps alone). 'Chaos' framing from the PD regime does not apply.
- Citations:
  - Hauert & Doebeli, Nature 428, 643-646 (2004), doi:10.1038/nature02360
  - Vainstein, Silva & Arenzon, J. Theor. Biol. 244, 722-728 (2007)
  - Sicardi, Fort, Vainstein & Arenzon, J. Theor. Biol. 256, 240-246 (2009), arXiv:0809.1479
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:060: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_corrlen, g_n_domains, g_activity_trend, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 4 types (random start). Agents have 3 types (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (every step).
- Standard name given: 2-D cyclic cellular automaton, Griffeath's 'Perfect' rule (R1/T3/C4/Moore): metastable turbulent soup with a large length scale that coarsens, with spiral cores possibly taking over
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Agents are inert (no cell-agent coupling), so they do not contribute. The knock-out baseline is a frozen random colouring, so the large effects mainly show that the CCA self-organises at all, which is expected. Negative activity trend = debris dying out during coarsening / turbulence relaxing on a finite 64x64 torus. This is consistent with Griffeath's description of the Perfect rule and with the metastability of T=3 turbulence on finite lattices noted by Shalizi et al.
- Citations:
  - http://psoup.math.wisc.edu/mcell/rullex_cycl.html
  - https://arxiv.org/abs/nlin/0507067
  - Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata, Statistics and Computing 1(1):23-39 (1991)
  - Fisch, Gravner & Griffeath, Cyclic cellular automata in two dimensions, Spatial Stochastic Processes (1991)
  - https://github.com/setanarut/ccae
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:062: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_n_domains, g_moran, g_cyclic_changes
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.3 per step); (2) neighbouring cells swap contents (every step).
- Standard name given: Cyclic cellular automaton (Griffeath 313: R1/T3/C3/Moore, asynchronous) with pair-exchange mobility, i.e. spatial rock-paper-scissors with mobility (Reichenbach-Mobilia-Frey regime below the critical mobility)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: All four features come from rule 1 (CCA 313 alone gives Moran 0.77, cycling 0.48, corrlen 4). Swapping alone gives nothing. Adding the swap lengthens the correlation length (4 to 6) and enlarges domains while slightly lowering Moran and cycling, which matches mobility-driven wavelength growth below the critical mobility (Reichenbach-Mobilia-Frey). Not an artifact: it keeps cycling, does not freeze, and the no-interaction baseline is random noise. Spiral forms were not seen directly; this is inferred from the measures.
- Citations:
  - Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39 (1991)
  - Fisch, Gravner & Griffeath, Cyclic Cellular Automata in Two Dimensions, Spatial Stochastic Processes (1991)
  - MCell Cyclic CA rule lexicon, rule 313 (psoup.math.wisc.edu/mcell/rullex_cycl.html)
  - Reichenbach, Mobilia & Frey, Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, Nature 448:1046-1049 (2007)
  - Szczesny, Mobilia & Rucklidge, The Influence of Mobility Rate on Spiral Waves in Spatial Rock-Paper-Scissors Games, Games 7:24 (2016)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:063: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_type_entropy, g_spec_peak
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.005 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.7 moves to a random empty cell.
- Standard name given: Schelling segregation with vacancies combined with majority/voter-type state change (Schelling-Voter class): coarsening domains and drift toward extinction of types
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Moving cells alone (rule 2) cannot change the type mix, so the entropy drop comes from the majority rule and is strongly sped up by movement. This matches the documented extinction/coarsening outcomes of Schelling-Voter models. The domain-count and entropy changes are genuine but still growing at 3x run length (not a steady state). g_spec_peak is likely a low-frequency artifact of slow drift (the knock-out is static, so it scores 0); detrend before treating it as an oscillation.
- Citations:
  - https://ui.adsabs.harvard.edu/abs/2017PhyA..487..125C/abstract
  - https://www.sciencedirect.com/science/article/abs/pii/S0960077913001720
  - https://arxiv.org/pdf/2204.01417
  - https://arxiv.org/pdf/0903.4694
  - https://arxiv.org/pdf/1311.5934
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:066: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_activity_trend, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 7 neighbours already are t+1 (with probability 0.03 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.
- Standard name given: Spatial snowdrift/Hawk-Dove game with Nowak-May synchronous imitate-the-best updating (Moore + self-interaction), perturbed by weak stochastic majority-rule noise
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- Notes: Agree with first-pass KNOWN, but reject its argument that anti-coordination gives many small interleaved patches: largest_domain ~0.92 shows one strategy (almost certainly C, from payoff counting with self-interaction) filling the lattice with small flickering islands of the other, matching Killingback-Doebeli/Nowak-May rather than Hauert-Doebeli stochastic results. Correlation length of 2 is short, not long. Rule 1 is a 2-type majority-vote noise rule, inert alone. The only synergistic feature (activity rising) fits noise preventing the deterministic CA from freezing into period-2 states; noisy Nowak-May games are studied (Mukherjee-Nowak-May 1996, cited from memory, not retrieved here). That feature halves at 3x run length, suggesting a transient or metric quirk. Not a discovery.
- Citations:
  - Nowak & May 1992, Nature 359:826-829 (Evolutionary games and spatial chaos)
  - Killingback & Doebeli 1996, Proc R Soc Lond B 263:1135-1144 (Spatial evolutionary game theory: hawks and doves revisited)
  - Hauert & Doebeli 2004, Nature 428:643-646 (Spatial structure often inhibits the evolution of cooperation in the snowdrift game)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:067: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_n_types, g_reversals, g_activity_trend
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) a type-1 cell spontaneously becomes type 0 (with probability 0.01 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour.
- Standard name given: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore neighbourhood + self), low-b coexistence regime: mostly static C/D mosaic with period-2 blinkers, plus a small amount of one-way D->C noise
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: The game rule alone reproduces every flagged feature; the one-way D->C noise only lowers the domain count and slows the decline in activity. A correlation length of 2 cells is short, so 'long correlation' overstates it. The state is a genuine but mostly frozen/periodic coexistence pattern (blinkers), not rich ongoing dynamics. The noise-only knock-out goes to the trivial all-cooperator state. I found no paper on noise that only turns defectors into cooperators, but Nowak et al. reported their results robust to the perturbations Mukherji et al. tested (random introduction of C and D; spontaneous C->D).
- Citations:
  - Nowak & May 1992, Nature 359:826 (https://www.nature.com/articles/359826a0)
  - Stanford Encyclopedia of Philosophy, Evolutionary Game Theory (blinkers / stable 2-cycle) https://plato.stanford.edu/Archives/spr2006/entries/game-evolutionary/
  - Mukherji, Rajan & Slagle 1996 Nature 379:125 and Nowak et al. reply, as summarised in arXiv:0910.0955
  - arXiv:0910.0961 (coexistence steady state for 1<b<2)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:069: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_types, g_type_entropy, g_interface, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (with probability 0.1 per step); (2) a type-0 cell spontaneously becomes type 1 (with probability 0.1 per step).
- Standard name given: Stochastic (asynchronous) 3-colour cyclic particle system / lattice rock-paper-scissors with unidirectional mutation: perpetual cycling ('fluctuation') with three-species coexistence
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: The observed features (3 types survive, mixing, rough interfaces, cells cycling) are generic, proven/documented properties of the 3-colour cyclic particle system; the 0->1 mutation only lowers the share of type 0 (entropy 1.0 -> 0.887), as expected for unidirectional mutation. The knock-out is a frozen two-type state, so the contrast is easy to beat; the full program does not exceed the rule-1-only sub-program on mixing or interface. Real dynamics, but generic; no coherent spirals are implied by these metrics.
- Citations:
  - Bramson & Griffeath, Flux and fixation in cyclic particle systems, Ann. Probab. 17 (1989) 26-45
  - Fisch, Gravner & Griffeath, Cyclic cellular automata in two dimensions (1991)
  - Mobilia, Oscillatory dynamics in rock-paper-scissors games with mutations, J. Theor. Biol. 264 (2010) 1-10, arXiv:0912.5179
  - Toupo & Strogatz, Nonlinear dynamics of the rock-paper-scissors game with mutations, arXiv:1502.03370
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:agents-2: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: a_step
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. Rules, applied in order each step: (1) each agent turns away from the centre of any type agents within distance 3.0.
- Standard name given: Self-trapping / caging (localisation within the interaction radius) of self-propelled particles with orientation-only turn-away (repulsive-Vicsek / anti-cohesion) interactions
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: True; web searches: 6
- Twist (if any): 
- Notes: Only a_step passed the checks: a strong, size-robust suppression of a step/displacement statistic (6.5 -> 0.78) caused by the turn-away rule. That is localisation/self-trapping, which is documented for orientation-only repulsive SPP models (agents confined to random motion within ~R, clustering at higher density). The first pass's variant twist (clumping plus local alignment without global order) relies on numbers outside the verified evidence, so the downgrade to KNOWN is mild. The behaviour is a real interaction effect but phenomenologically near-trivial: likely back-and-forth jitter at about 28 neighbours per agent, possibly amplified by synchronous discrete-time updating. Not a novel emergent pattern. Exact definition of a_step is unknown, which limits confidence.
- Citations:
  - arXiv:2412.10577 Lanes and lattice structures in a repulsive model for self-propelled agents
  - arXiv:2604.24232 Fava, Ginelli, Mahault - Density protected states in active matter under virtual confinement (discrete-time repelling-torque Vicsek-type model)
  - Das et al., Flocking by Turning Away, Phys. Rev. X 14, 031008 (2024), arXiv:2401.17153
  - arXiv:1403.0697 Aggregation and Segregation of Confined Active Particles (self-trapping suppresses effective self-propulsion)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:agents-8: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: a_step
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0.
- Standard name given: Chemotactic trapping / gradient-ascent localization of persistent random walkers (auto-chemotactic aggregation, Keller-Segel, if agents secrete the field)
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Agrees with the first pass's KNOWN label but rejects its basis: the density CV it cited is not among the features that passed checks, and no emission rule is listed. The only verified feature, a_step (0.79 with gradient-following vs 19.75 without, robust across size and run length), fits agents climbing to field maxima and staying there, which reduces displacement. This follows from single-agent gradient ascent on any persistent field (secreted or a slowly relaxing initial field), so it is likely a direct consequence of the rule rather than an emergent collective pattern. Test: check for a field source term and compare against a frozen copy of the initial field.
- Citations:
  - Keller & Segel 1970; Chavanis, A stochastic Keller-Segel model of chemotaxis, arXiv:0804.4425
  - Meyer, Schimansky-Geier, Romanczuk, Active Brownian agents with concentration-dependent chemotactic sensitivity, PRE 89, 022711 (arXiv:1310.8272)
  - Channel formation enhances target consumption by chemotactic ABPs, arXiv:2507.11328
  - Unusual diffusion of chiral active Brownian particles, Commun. Phys. (s42005-026-02591-x)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:grid-19: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_activity_cv, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.003 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour.
- Standard name given: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore+self) with weak noise: near-absorbing takeover by one strategy, perturbed by rare flips
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: With the game alone (rule 2), the lattice freezes into one strategy almost everywhere (largest domain 0.99, activity variability 0), as documented for deterministic imitate-the-best. The 0.003 majority-adoption rule acts as rare noise: high activity variability is just intermittent flips over near-zero mean activity, the largest domain reflects the near-takeover state, and the correlation length (5.5) reflects a few leftover minority clusters. Effect sizes shrink at larger area and longer runs, consistent with a slow approach to an absorbing state. Majority-vote noise in Nowak-May was not found specifically, but no new phenomenology is evident. The first pass's label was KNOWN but its text said KNOWN-VARIANT; KNOWN is confirmed, with the added judgement that the features are largely trivial.
- Citations:
  - https://arxiv.org/pdf/cond-mat/0607344
  - https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0013749
  - https://arxiv.org/html/0910.0955v2
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:grid-25: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_activity_cv
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) neighbouring cells swap contents (every step); (2) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step).
- Standard name given: Stochastic threshold cyclic cellular automaton (Fisch-Gravner-Griffeath CCA, kappa=4, Moore range 1, theta=3, 'turbulent'/metastable regime) with pair-exchange mobility in the high-mobility, well-mixed regime of spatial rock-paper-scissors (Reichenbach-Mobilia-Frey)
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Agrees with the first pass's KNOWN label but rejects the KNOWN-VARIANT wording in its text: the only novelty claimed was the rule combination. The knock-out removes both rules (swapping also reads neighbours), so it is a frozen absorbing state and the 3-spread differences are uninformative. Swap-only already gives reversals 0.40 and activity_cv 0.52 (higher than the full program's 0.44); CCA-only gives 0 reversals. The full program's domain count (6.5) lies between CCA-only (4.4) and swap-only (7.7), consistent with diffusion blurring CCA domains. Swapping every step against a 0.1 conversion probability puts the system in the known high-mobility, well-mixed regime with near-uniform stochastic cycling (mean-field: about 32% of cells meet the threshold when types are evenly mixed). The features are mostly mixing artifacts plus slow cyclic advance, not a new emergent pattern.
- Citations:
  - Fisch, Gravner & Griffeath (1991) Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39
  - Shalizi et al., Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067
  - Reichenbach, Mobilia & Frey (2007) Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, Nature 448:1046-1049 (arXiv 0709.0217)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:grid-6: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_largest_domain, g_n_domains
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) a type-0 cell becomes type 1 if at least 1 neighbours are type 1 (with probability 0.1 per step).
- Standard name given: Spatial SI epidemic cellular automaton / multi-seed Eden (Richardson) growth = site first-passage percolation (geometric waiting times), ending in the absorbing all-infected state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Rule is k=1 irreversible contact infection; on a finite torus from a ~50% random start it surely reaches the all-1 absorbing state. g_largest_domain=1.0 is that absorbing state. g_n_domains=0.6931=ln 2, most likely log(1+1), i.e. a single domain, so the 'many domains' label contradicts the value and is an artifact. The knock-out comparison is meaningless because the non-interacting program leaves the random colouring frozen. The first pass's verdict is correct, but its discussion of k=3 bootstrap percolation does not apply to this rule.
- Citations:
  - https://arxiv.org/html/1909.12001v1
  - https://link.springer.com/content/pdf/10.1007/s00454-022-00474-w.pdf
  - https://arxiv.org/html/1005.0649
  - https://link.springer.com/article/10.1007/s42979-021-00619-3
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:grid-4: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_types, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.03 per step).
- Standard name given: Nowak-May spatial Prisoner's Dilemma (unconditional imitation, Moore+self) with voter-model imitation noise; fixation to an absorbing monomorphic (likely all-C) state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Both features (n_types=1, largest_domain=1) are the consensus/fixation endpoint. At b=1.3 the game freezes into cooperator-dominated domains, and neighbour-copy noise (which cannot reintroduce lost types) erodes the residual defector strips until absorption. The knockout difference only reflects timescales: the voter step alone is too slow to fix on 64x64 within the run, and the game alone freezes. The 'transient' measurement label contradicts the reported values. No twist in the phenomenology beyond expected noise-assisted fixation.
- Citations:
  - Nowak & May 1992, Nature 359:826 (Evolutionary games and spatial chaos)
  - Schweitzer, Behera & Muhlenbein, arXiv cond-mat/0211605 (b-regimes: 1.0<b<1.33 large cooperator domains)
  - Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344
  - arXiv 0706.0648 (all-C/all-D absorbing under unconditional imitation without mutation)
  - Scale-invariant behavior in a spatial game of prisoners' dilemma (1<b<1.5 static, cooperator-dominated)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:grid-7: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_reversals, g_interface
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 2 field(s) with diffusion [0.02, 0.1]. Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.1 per step); (2) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Majority-vote model with noise (Oliveira 1992) / zero-T Glauber-Ising dynamics plus random spin flips, high-noise disordered phase
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: 2-state cyclic advance = spin flip; agents ignore the lattice, so they are just spatially uniform random noise at about 10% of cells per step, comparable to the 10% majority-update rate. Rough estimate puts effective noise far above q_c (about 0.075 for 4 neighbours), so the system is in the disordered phase: g_interface 0.43 vs 0.50 random means only weak short-range correlation, not arrested coarsening as the first pass described. g_reversals = 1.0 already in the rules-2-only sub-program, so it comes from the flip rule, not the interaction. Diffusing fields are inert. Agreeing with the KNOWN label, correcting the first pass's description.
- Citations:
  - M.J. de Oliveira, Isotropic majority-vote model on a square lattice, J. Stat. Phys. 66, 273 (1992)
  - Kwak, Yang, Sohn, Kim, Critical behavior of the majority voter model is independent of transition rates, PRE 75, 061110 (2007)
  - F.W.S. Lima, arXiv:1110.0251 / arXiv:1307.1776 (majority-vote model variants on square lattice)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:grid-8: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 2 field(s) with diffusion [0.02, 0.1]. Rules, applied in order each step: (1) each agent advances the type of the cell under it by one (cyclically); (2) each agent turns up the gradient of field 0.
- Standard name given: Multi-walker parity toggling (Langton's-ant/turmite-type site flipping) plus chemotactic accumulation of active particles at maxima of a passive diffusing field
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- Notes: Overturns the first pass's reasoning, not its overall verdict. The first pass claimed rule 2 is inert, but the sub-program data shows rule 2 causes the whole activity decline (rule 1 alone gives 0.0017; both rules give -3.29), so field 0 must start non-uniform. Reversals = 1 is automatic for parity flipping and goes to 0 when flipping is removed. The activity decline is most plausibly agents gathering at field peaks: fewer distinct cells flipped, and even numbers of co-located flips cancel. This was inferred, not checked with a cluster measurement. Because there is no feedback (cells don't affect agents, agents don't emit into the field), this is a relaxation transient, not self-organised lattice structure. The feature label 'activity rising' contradicts its negative value; check the sign convention.
- Citations:
  - https://en.wikipedia.org/wiki/Langton%27s_ant
  - https://tetrahedral.blogspot.com/2010/09/langtons-ant.html
  - https://www.nature.com/articles/s41598-021-83963-x
  - https://arxiv.org/html/2011.04580v1
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:004: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_activity_cv, g_corrlen, g_n_types, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell adopts its neighbourhood's majority type (with probability 0.1 per step).
- Standard name given: Fixation/consensus into an absorbing homogeneous state (full-C or full-D, bistable) in spatial evolutionary games mixing payoff-driven imitation and conformity (majority) updating
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Observed features describe collapse to a single uniform frozen type; corrlen and largest_domain follow trivially from uniformity. Diagnostics are internally inconsistent: transient-window values show the collapsed state, 'activity rising' is strongly negative, and corrlen median 0.0 contradicts the 'long correlation' label, suggesting metric degeneracy. Knock-out contrast is uninformative because the knock-out has no dynamics. Single-rule sub-programs not fixating fits known lattice snowdrift coexistence and majority-vote freezing; this is background knowledge, not separately searched.
- Citations:
  - Lin, Huang, Dai, Yang (2020) Chaos Solitons Fractals 140:110146 - Evolutionary game dynamics of combining the payoff-driven and conformity-driven update rules
  - Szolnoki & Perc (2018) arXiv:1809.05079 - Competition and partnership between conformity and payoff-based imitations in social dilemmas
  - Szolnoki & Perc (2015) J. R. Soc. Interface 12:20141299 - Conformity enhances network reciprocity in evolutionary social dilemmas
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:012: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_n_domains, g_reversals, g_activity
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) neighbouring cells swap contents (every step).
- Standard name given: Nowak-May spatial Snowdrift (Hawk-Dove) game with random mobility/diffusion (Vainstein-Silva-Arenzon / Sicardi et al. mobility model)
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- Notes: Payoffs T>R>S>P make this a Snowdrift game, not a PD. Most of the activity comes from swapping alone (0.33 of 0.45); the game alone nearly freezes (0.003). Reversals = 1.0 in both sub-programs, so that feature does not tell the game apart from swapping. The knock-out removes swapping, which inflates the contrast. The domain reduction is ordinary imitation-driven clustering competing with stirring. No twist in the behaviour beyond the documented mobility + spatial Snowdrift phenomenology; differences from the papers (Moore neighbourhood, swapping every step) are in the rules only.
- Citations:
  - Vainstein, Silva & Arenzon, Does mobility decrease cooperation?, J. Theor. Biol. 244, 722-728 (2006/2007)
  - Sicardi, Fort, Vainstein & Arenzon, Random mobility and spatial structure often enhance cooperation, J. Theor. Biol. 256, 240-246 (2009)
  - Vainstein & Arenzon, Spatial social dilemmas: dilution, mobility and grouping effects with imitation dynamics, Physica A 394 (2014)
  - Nowak & May, Evolutionary games and spatial chaos, Nature 359 (1992)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:016: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_largest_domain, g_n_types, g_reversals, g_n_domains, g_spec_peak
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell adopts its neighbourhood's majority type (with probability 0.01 per step).
- Standard name given: Fixation to a single absorbing state in a spatial snowdrift game with unconditional imitation plus conformity (majority) updating
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Same verdict as the first pass, but the game is a snowdrift game (T>R>S>P), not the S=0 Nowak-May prisoner's dilemma. With unconditional imitation and 8 neighbours at T=1.5, cooperation is high and mixed; the 1% majority rule removes the leftover minority pockets, giving full fixation, as in conformity-plus-payoff studies where full-C/full-D phases grow. n_domains=0.6931=ln2 means one domain. The end state is a frozen absorbing state: interaction-driven, but trivial as a pattern. The temporal peak is most likely a leftover transient (inferred, not checked).
- Citations:
  - Szolnoki & Perc, Competition and partnership between conformity and payoff-based imitations in social dilemmas, arXiv:1809.05079
  - Roca, Cuesta, Sanchez, Effect of spatial structure on the evolution of cooperation, arXiv:0806.1649
  - Roca, Cuesta, Sanchez, Evolutionary game theory: temporal and spatial effects beyond replicator dynamics, arXiv:0911.1720
  - Hauert & Doebeli, Nature 428, 643 (2004)
  - Nowak & May, Nature 359, 826 (1992)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:023: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_n_types, g_n_domains, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.002 per step).
- Standard name given: Spatial Prisoner's Dilemma (Nowak-May type, strict PD payoffs) - fixation/complete invasion of a single strategy (most likely all-defector absorbing state)
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Verdict label agrees with first pass, but first pass's characterization (chaotic coexistence, avalanches, persistent domains) is wrong: the measured features (1 surviving type, largest domain 1.0, domain score at floor) show total fixation to one strategy, a frozen absorbing state. Game rule alone gives ~99.6% takeover; rare neighbour-copy noise just mops up remnants, after which the uniform state is absorbing. 'Reversals' is likely transient churn before fixation. With S=-0.5, T=1.5 and random 50/50 start, defectors (payoff ~6 vs ~3 for cooperators) should invade; complete defector invasion is a documented regime for strict-PD payoffs (Schweitzer et al. regime E). Which type wins was not directly reported by the simulation, so identity of the winner is inferred.
- Citations:
  - https://arxiv.org/pdf/cond-mat/0211605
  - https://www.researchgate.net/publication/11497714_Scale-invariant_behavior_in_a_spatial_game_of_prisoners'_dilemma
  - https://arxiv.org/pdf/2411.13741
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:031: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_activity_cv, g_n_domains, g_n_types, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) a cell of type t becomes type t+1 (cyclically) if at least 6 neighbours already are t+1 (with probability 0.03 per step).
- Standard name given: Nowak-May imitate-the-best spatial prisoner's dilemma with conformity (majority-threshold) updating; coarsening to a single-strategy absorbing state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: With 2 types, the 'cyclic' rule 2 is just a random flip-to-the-majority (conformity) rule. Rule 1 alone already nearly reaches a single type (largest domain 0.989, n_types 1.5). Rule 2 removes the leftover minority clusters, so the full program reaches one type everywhere and freezes. activity_cv is most likely a transient burst followed by a frozen phase; reversals=1 fits a single takeover; corrlen=3.0 at every size and run length looks like a capped value on a uniform grid. Which type wins at T=1.5, S=-0.5 was not checked.
- Citations:
  - https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0013749
  - https://wisc.pb.unizin.org/agent-based-evolutionary-game-dynamics/chapter/iii-1/
  - https://core.ac.uk/outputs/42932284/
  - https://arxiv.org/html/1604.06438
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:037: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_type_entropy, g_reversals, g_interface
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 3 types (random start). Agents have 4 types (random start). Agents move at speed 0.0 with heading noise 1.0 rad per step. Rules, applied in order each step: (1) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.2 moves to a random empty cell; (2) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Schelling segregation model in the high-tolerance (non-segregating) regime, with fixed cyclic forcing sites (stationary agents acting as period-3 clocks)
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: With a 0.2 threshold almost no cell is ever unhappy (rule-1-only run matches the knock-out), so Schelling is close to inert. Speed-0 agents permanently force about 10% of cells around 0->1->2->0, which directly produces g_reversals. The drops in entropy, interface and domain count follow from the forcing changing the type mix (create/delete plus relocation to random empty sites); this mechanism is my inference and has not been checked against per-type fractions. The first pass's 'dynamical mixed phase' explanation is wrong: the churn comes from the agents, not from Schelling's own dynamics. Verdict unchanged.
- Citations:
  - Gauvin, Vannimenus & Nadal, Phase diagram of a Schelling segregation model, Eur. Phys. J. B 70, 293 (2009), arXiv:0903.4694
  - Dall'Asta, Castellano & Marsili, Statistical physics of the Schelling model of segregation, arXiv:0707.1681
  - Barthelemy, The statistical physics of cities, arXiv:1905.01953
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:039: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_n_domains, g_activity_cv, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) a type-0 cell becomes type 1 if at least 4 neighbours are type 1 (with probability 0.03 per step).
- Standard name given: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore+self) with cooperator-to-defector noise (Mukherji et al. 1996 type), collapsing to a near-absorbing all-defector state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Largest domain about 0.97 and n_domains about ln3, so the lattice is nearly uniform. Long corrlen is a by-product of uniformity; activity CV and 1 reversal reflect a transient plus residual stochastic flips. Rule 2 is one-way C->D pressure, and alone it fills the lattice (largest domain 1.0). Rule 1 alone already reaches 0.991. The knock-out baseline is a frozen random field, so the large differences are trivial. Only one search was run and the exact neighbour-threshold defection rule was not located in a paper; the verdict rests on the documented behaviour, not the rule tuple.
- Citations:
  - Nowak & May 1992, Nature 359:826 'Evolutionary games and spatial chaos'
  - https://wisc.pb.unizin.org/agent-based-evolutionary-game-dynamics/chapter/iii-1/
  - https://arxiv.org/pdf/0910.0955
  - https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0013749
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:041: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_activity_cv, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Noisy gradient climbing (chemotaxis) on a fixed or decaying random field, giving temporary trapping at local maxima; agents write a symbol-blind 2-state parity mark on the lattice (degenerate turmite)
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Same verdict as the first pass (KNOWN, trivial), but its mechanism is corrected: the gradient rule is not inert. The full program differs strongly from the rule-2-only sub-program, which fits agents being temporarily trapped at local maxima of the random starting field as it slowly smooths out (D=0.05). The falling activity is most likely agents crowding onto the same cells, where two hits cancel in a 2-state cycle. The high CV comes from dividing by a near-zero mean and halves at 4x area. Reversals=1 is automatic with 2 states. There is no feedback loop (agents never add to the field or read the cells). The knock-out test can't tell 'clustered' from 'random walkers', so passing it doesn't show emergence. Mechanism inferred from rules and summary statistics; not run.
- Citations:
  - https://arxiv.org/pdf/2511.18745
  - https://www.nature.com/articles/s41598-021-83963-x
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:050: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity_trend, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) a type-0 cell becomes type 3 where field 1 exceeds 0.05; (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Schelling segregation model (vacancy-limited / frozen-jammed regime)
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Inferred, not simulated, and field starting values were not given. Rule 1 seems to turn empty cells into type 3 over much of the grid, since field 1 has no source or sink. Each Schelling move frees a site that is then refilled, so empty cells run out and the system jams. That explains falling activity, partial coarsening that stalls (fewer domains than the knock-out) and the 'reversals' count (occupied -> empty -> type 3 turnover). The vacancy-refilling rule is unusual but adds no new behaviour beyond the documented low-vacancy frozen/jammed Schelling phase.
- Citations:
  - Gauvin, Vannimenus, Nadal 2009, Phase diagram of a Schelling segregation model, arXiv:0903.4694
  - Rogers & McKane 2012, Jamming and pattern formation in models of segregation, arXiv:1204.5400
  - Abella et al. 2022, Aging effects in Schelling segregation model, Sci Rep s41598-022-23224-7
  - Dall'Asta, Castellano, Marsili 2008, Statistical physics of the Schelling model, arXiv:0707.1681
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:052: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_largest_domain, g_n_types
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). 2 field(s) with diffusion [0.02, 0.01]. Rules, applied in order each step: (1) field 1 grows autocatalytically by consuming field 0 (rate 1.0 x f x g^2); (2) a type-0 cell becomes type 1 where field 1 exceeds 0.8.
- Standard name given: Front-driven fixation into an absorbing state: closed cubic autocatalysis (A+2B->3B, Gray-Scott without feed or decay) travelling fronts, recorded by an irreversible threshold switch of cell type
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: With no feed or decay, field 1 is never removed and simply converts field 0 until the substrate is used up; rule 2 only switches 0->1 and never back. The features (one type, largest domain = 1.0) therefore mean the fronts filled the whole torus and every cell was frozen as type 1, not that coarsening happened. Both rules are required only because one writes the field and the other reads it. Holding at larger size and longer runs reflects a guaranteed outcome, not self-organisation. If the initial field total were mostly below 0.8, the result would instead be a partly converted, frozen pattern.
- Citations:
  - Billingham & Needham, Proc. R. Soc. Lond. A 430 (1990) 'The development of travelling waves in a simple isothermal chemical system III. Cubic and mixed autocatalysis' https://royalsocietypublishing.org/rspa/article/430/1880/509/16458
  - Biscale Chaos in Propagating Fronts, arXiv chao-dyn/9601004 (closed A+2B->3B reaction-diffusion fronts in 2D)
  - Billingham & Needham, Dyn. Stab. Syst. 6(1) 1991, travelling-wave solutions in cubic autocatalysis
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:053: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_n_domains, g_n_types, g_activity_trend, g_type_entropy
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.001 per step); (2) a type-0 cell becomes type 1 where field 0 exceeds 0.5.
- Standard name given: Zero-temperature / low-noise q-state Potts (majority-rule) coarsening with pinned multi-domain states, plus a trivial irreversible field-threshold type conversion
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Rule-2-only entropy 0.7511 matches type 0 being fully folded into type 1 (shares 1/2,1/4,1/4 -> 1.5/2 bits = 0.75), so field 0 is above 0.5 at essentially every type-0 cell at some point and the passive field just converts them all once; field 1 is never read. The lost type and lower entropy are therefore built into the rules, not emergent. The fewer domains and falling activity are standard q>2 Potts coarsening with T-junction pinning. The knock-out removes both rules, leaving a frozen random lattice, which makes the interaction effects look larger than they are. Full-program entropy 0.65 is 3-type coarsening starting with type 1 at 50%.
- Citations:
  - arXiv:0912.3166 (Loureiro & Arenzon, curvature-driven coarsening in 2D Potts)
  - Olejarz, Krapivsky & Redner, Zero-temperature coarsening in the 2d Potts model (2013)
  - arXiv:2101.02425 (High-degeneracy Potts coarsening)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:056: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 3 types (random start). Agents have 3 types (random start). Agents move at speed 0.1 with heading noise 0.1 rad per step. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 1.0; (2) each cell adopts its neighbourhood's majority type (every step); (3) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Period-2 limit cycles of synchronous majority-rule cellular automata (Goles-Olivos / Poljak-Sura theorem), plus local agent forcing that majority voting undoes
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Same KNOWN verdict as the first pass, but for a different reason: the first pass explained domain structure via frozen T-junctions, which does not cover the flipping feature. Rule 2 alone gives more reversals (0.85) than the full program (0.58), so nothing new appears when rules are combined and the agents reduce flipping. The knock-out baseline of 0 is misleading because it removes the majority rule itself. Flocking (rule 1) has no measurable effect on this feature. Suggested check: rerun with cells updated one at a time; the reversals should largely vanish.
- Citations:
  - Goles & Olivos, Periodic behaviour of generalized threshold functions, Discrete Math. 30 (1980) 187-189
  - Goles & Olivos, Discrete Applied Math 3 (1981) 93-105; Poljak & Sura (multi-colour period <=2)
  - arXiv:1711.10920 (Biased) Majority Rule Cellular Automata
  - arXiv:2205.08972 (2-cycles exist only under synchronous update)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:064: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: f_mean, f_change, f_n_blobs, f_std
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 3 types (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) a type-1 cell becomes type 0 if at least 4 neighbours are type 0 (every step); (2) type-0 cells release 1.0 of field 0 per step.
- Standard name given: Threshold (r=4, Moore) bootstrap percolation driving an inert, source-only (no-decay) diffusing field
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Starting at 50% density with threshold 4/8 is far above criticality, so the lattice almost fully fills with type 0 (rules-1-only f_mean ~0.98, frozen absorbing state). The field has constant sources and no decay: f_change doubles (5->10) because source density doubles; f_mean grows linearly with run length (gap 500 -> 1499 at 3x), i.e. diverging; f_std ~0.035 and f_n_blobs ~ln2 mean a spatially uniform field. No emergent pattern; the knock-out difference is just random source map vs uniform source map. Agents are uncoupled.
- Citations:
  - Chalupa, Leath & Reich 1979, J. Phys. C 12 L31 (origin of bootstrap percolation, cited via Balogh-Peres-Pete arXiv math/0311125)
  - Aizenman & Lebowitz 1988, J. Phys. A 21 3801 (finite-box bootstrap percolation)
  - Balogh, Peres, Pete, 'Bootstrap percolation on infinite trees and non-amenable groups', arXiv:math/0311125
  - Balister et al., 'The time of bootstrap percolation in two dimensions' (memphis.edu/msci/people/pbalistr/timeboot.pdf)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:065: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: a_density_cv, a_local_align, a_nn_ratio
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 4 types (random start). Agents move at speed 2.0 with heading noise 0.0 rad per step. 2 field(s) with diffusion [0.02, 0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) field 1 grows autocatalytically by consuming field 0 (rate 1.0 x f x g^2).
- Standard name given: Chemotactic accumulation of non-interacting self-propelled tracers at maxima of an external (depleting, Gray-Scott-type substrate) field
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Agents are pure tracers: no agent-agent coupling, no deposition, so the knock-out only removes gradient-following. Rule 1 alone gives stronger clustering and alignment than the full program; rule 2 merely depletes and moves the field-0 peaks. Alignment comes from noise-free agents converging on, or orbiting, the same peak, not from any collective mechanism. As described, field 0 has no feed term, so the field is likely winding down. Verdict agrees with the first pass, but its Keller-Segel/autochemotaxis discussion does not apply to this rule set.
- Citations:
  - https://www.nature.com/articles/s41598-021-83963-x
  - https://arxiv.org/pdf/1202.2442
  - https://arxiv.org/pdf/2410.04735
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:011: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: n_components
- Rule: N kN=2 | N[0] ALWAYS() -> REWIRE_SAME() p=1.0  (compact grammar notation)
- Standard name given: Rewire-to-same fission / fragmented absorbing phase of the coevolving (adaptive) voter model, frozen-opinion (alpha=1, phi=1) limit
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Fixed binary types + unconditional rewire-to-same-type necessarily drive the graph to a state with no mixed-type edges, so at least one component per type (>=2). The observed ~2.9 vs ~1.5 (knock-out) is guaranteed by the rule plus leftover isolated or small pieces. Robust and not a numerical artifact, but a frozen absorbing state fixed in advance rather than a genuine emergent pattern.
- Citations:
  - Holme & Newman, Phys. Rev. E 74, 056108 (2006)
  - Durrett et al., Graph fission in an evolving voter model, PNAS 109:3682-3687 (2012)
  - Nonlinear coevolving voter model phase diagram incl. absorbing phase in a fragmented network (via ResearchGate summary)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:015: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: n_cyclic_changes
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) a node of type t becomes t+1 (cyclically) if at least 2 neighbours are t+1 (with probability 0.1 per step).
- Standard name given: Threshold cyclic cellular automaton (3-colour CCA, threshold 2) with asynchronous/random-timing updates on a random graph; fixation vs fluctuation dichotomy
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: The verdict agrees with the first pass, but I disagree with its claim that the system freezes: cycling persists at 3x run length (0.86 of original), which fits the documented fluctuating phase or slow fixation on a finite sparse graph. The measured feature 'n_cyclic_changes' just reads out the rule (the knock-out cannot change state at all), so it shows no emergent behaviour; measures like fixation time, cluster sizes or persistent loop waves would be needed.
- Citations:
  - Bello & Sivakoff, Cyclic Cellular Automata and Greenberg-Hastings Models on Regular Trees, arXiv:2108.06404
  - Gravner, Lyu & Sivakoff, Limiting behavior of 3-color excitable media on arbitrary graphs, Ann. Appl. Probab. 28(6) (arXiv:1610.07320)
  - Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata
  - Bramson & Griffeath, Flux and fixation in cyclic particle systems, Ann. Probab. 1989
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:021: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_activity
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell copies the type of a random neighbour (every step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.
- Standard name given: Nowak-May imitate-best spatial game with snowdrift (hawk-dove) payoffs plus voter-model noise (noisy Nowak-May / spatial snowdrift)
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Spatial features (corrlen 2, ~5 domains) are already present with the game rule alone (4.7 domains, corrlen 2.0), consistent with Hauert-Doebeli small filament clusters. Activity and reversals come from voter noise being repeatedly undone by deterministic imitate-best (known perturbation sensitivity); the frozen no-interaction baseline makes those differences guaranteed. Near-50/50 composition is more likely due to noise-driven mixing than the well-mixed snowdrift equilibrium, since the lattice snowdrift game is known to usually sit below it. The trivial/artifact flag is a judgement call (low-medium confidence): the clustering is genuine but known; the dynamics are mostly noise churn.
- Citations:
  - Nowak & May 1992, Nature 359:826, Evolutionary games and spatial chaos
  - Hauert & Doebeli 2004, Nature 428:643, doi:10.1038/nature02360
  - Mukherji, Rajan & Slagle 1996 (noise/update errors in Nowak-May), via wisc.pb.unizin.org chapter III-1
  - arXiv:2311.08119 (Hamming-distance sensitivity in deterministic spatial games)
  - arXiv:1102.4570 (filament-like cooperator clusters in spatial snowdrift)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:023: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_activity_trend
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.1 per step).
- Standard name given: Spatial snowdrift game with unconditional imitation mixed with voter-model social imitation (strategic + social imitation -> consensus)
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 4
- Twist (if any): 
- Notes: The verdict stays KNOWN, but the first pass's explanation (Hauert-Doebeli 'space suppresses cooperation') is probably wrong for this rule: with best-neighbour copying at T=1.5 on an 8-neighbour lattice the literature predicts cooperation is fostered. The full program beating both sub-programs on domain size is the documented result that voter-type copying acts as interface noise and drives the lattice to consensus (Vilone et al.). The observed features describe a near-absorbed consensus state on a finite torus, a relaxation towards an absorbing state rather than a lasting pattern; the frozen-random-start baseline is weak. Check which type dominates: if it is defectors, that would conflict with Roca et al.'s best-neighbour-copying result at T=1.5.
- Citations:
  - Hauert & Doebeli, Nature 428, 643 (2004)
  - Roca, Cuesta & Sanchez, arXiv:0911.1720 (Phys. Life Rev. 2009) - UI on k=8 lattice fosters cooperation in SD below T~1.7
  - Li, Ke, Lin & Hui, Phys. Rev. E 85, 021111 (2012)
  - Vilone, Ramasco, Sanchez & San Miguel, Phys. Rev. E 90, 022810 (2014), arXiv:1405.4997
  - Vilone et al., Sci. Rep. 2, 686 (2012)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:025: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_types
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Majority-vote (zero-temperature Ising-like) coarsening to consensus; noise/mobility-assisted escape from frozen two-domain stripe states
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- Notes: The only surviving feature is the number of surviving types. Full program vs knock-out (1 vs 3) is trivial because nothing can disappear without interaction. Majority vote alone gives 2, which matches known frozen two-stripe states of Ising-like majority coarsening. Schelling moves act as noise that lets the system escape those states and reach a single type, a standard noise-assisted consensus effect. The endpoint is a frozen absorbing state that may be all-empty (type 0 also takes part in the vote). The effect halves at 4x area, consistent with finite-size consensus. Ising results are for two-state random-sequential updates, so they are qualitative analogues for this 3-type synchronous rule.
- Citations:
  - Moore C. 1997, Majority-Vote Cellular Automata, Ising Dynamics, and P-Completeness, J. Stat. Phys. 88:795
  - Spirin, Krapivsky, Redner 2001/2002, Fate of / Freezing in zero-temperature Ising ferromagnets, PRE 65:016119 (arXiv cond-mat/0008019, cond-mat/0105037)
  - Barros, Krapivsky, Redner 2009, arXiv:0905.3521; Olejarz, Krapivsky, Redner 2012, arXiv:1208.2944
  - Dominguez Arroyo 2020, Artificial Life Simulation of the Majority Vote (RUG thesis)
  - arXiv:2505.03464, Barrier induced stalemate-consensus transition of self-propelled participants subject to majority rule
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:026: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_activity_cv
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) neighbouring cells swap contents (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.
- Standard name given: Noisy Schelling segregation with vacancies (Kawasaki/exchange noise), mixed/weakly clustered dynamical phase
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: The first-pass reasoning wrongly cites a COPY/Glauber rule that this system does not have; the KNOWN conclusion still holds because Kawasaki-type swap noise plus Schelling relocation is covered by Barret & Torri. g_activity_cv (0.756 swap-only vs 0.752 full) and g_reversals (0.50 swap-only) come from the swap noise; the knock-out's zeros are a frozen baseline created by removing the swap. Only g_n_domains (7.36 -> 5.9) reflects the interaction: modest extra clustering in a mixed phase dominated by noise. Marked trivial_or_artifact=true because 2 of 3 flagged features are noise artifacts; the weak clustering itself is real but documented.
- Citations:
  - Gauvin, Vannimenus, Nadal, Phase diagram of a Schelling segregation model, EPJB 70, 293 (2009), arXiv:0903.4694
  - Barret & Torri, Hydrodynamic limit of the Schelling model with spontaneous Glauber and Kawasaki dynamics, Electron. J. Probab. 29 (2024), arXiv:2302.09866
  - Abella et al., Aging effects in Schelling segregation model, Sci. Rep. (2022), arXiv:2204.01417
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:036: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: n_oscillation
- Rule: N kN=4 | N[any] ALWAYS() -> COPY_RAND() p=0.2  (compact grammar notation)
- Standard name given: Voter model (pull/imitation dynamics) on sparse random graph; mixed-state flipping before consensus
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Agrees with first-pass KNOWN, but the flagged 'n_oscillation' is not a genuine periodic/emergent oscillation: it is random opinion flipping in the long-lived mixed state of the voter model, which ends in a frozen consensus. The knock-out gives exactly 0 because without copying nothing changes. The signal halves at 3x run length, which fits slow relaxation towards consensus, not a sustained oscillation. The score holds at 4x area because consensus time grows with network size.
- Citations:
  - https://www.eurandom.tue.nl/wp-content/uploads/2023/01/VMRG.pdf
  - https://arxiv.org/html/1506.05321
  - https://arxiv.org/pdf/2507.11521
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:042: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_n_domains, g_corrlen, g_activity_trend, g_reversals
- Rule: C kC=2 | C[any] CNT_LE(b=0,th=3) -> IMITATE_BEST(g=(1.9, 0.0)) p=0.1  (compact grammar notation)
- Standard name given: Asynchronous / semi-synchronous Nowak-May spatial Prisoner's Dilemma (imitate-best, b~1.9): loss of spatial chaos, short-range coarsening into frozen C/D domains
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- Notes: Agrees with first-pass KNOWN. The gate (only cells with <=3 cooperator neighbours may update) is a rule change that only adds freezing; no new behaviour. Observed features (about 4 domains, correlation length 2, falling activity, about 1 reversal per cell) describe brief imitation-driven coarsening that then freezes, a near-absorbing static state. The clustering comes from the interaction but is weak and not self-sustaining, so it is marked trivial/frozen rather than a rich emergent pattern. The 'activity rising' label is the wrong way round: the measured trend is negative.
- Citations:
  - Nowak & May 1992 Nature 359:826 (Evolutionary games and spatial chaos)
  - Huberman & Glance 1993 PNAS 90:7716 doi:10.1073/pnas.90.16.7716
  - Nowak, Bonhoeffer & May 1994 PNAS 91:4877 doi:10.1073/pnas.91.11.4877
  - Prisoner's Dilemma with Semi-synchronous Updates: Evidence for a First Order Phase Transition, arXiv:0910.0961
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:053: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: n_edge_turnover, n_oscillation
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 4 types (random start). Rules, applied in order each step: (1) a node of type t becomes t+1 (cyclically) if at least 1 neighbours are t+1 (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (every step).
- Standard name given: Cyclic-dominance (cyclic CA / rock-paper-scissors-type) dynamics on a coevolving network with homophilic rewiring; active vs fragmented-frozen phases
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Agree with KNOWN, but the first pass's 'freezes in 1-2 steps' argument conflicts with the persistent late-run oscillation, so the implementation evidently leaves some mismatched links each step. Oscillation is inherited from rule 1 alone (0.82) and is weakened by rewiring, shrinking at larger size and longer runs: consistent with slow approach to the known fragmented absorbing phase. Edge turnover ~0.98 and unchanged with run length, so it probably counts all links ever rewired (mostly the initial burst). Knock-out comparison is weak because both features are zero by definition without interaction. No new phenomenology.
- Citations:
  - Demirel, Prizak, Reddy & Gross, Opinion formation and cyclic dominance in adaptive networks, EPJ B 84, 541 (2011), arXiv:1011.1124
  - Adaptive cyclically dominating game on co-evolving networks, EPJ B (2017), arXiv:1605.03370
  - Szolnoki et al., Cyclic dominance in evolutionary games: A review, arXiv:1408.6828
  - Holme & Newman coevolving voter model (rewire to same-opinion node)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:057: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a type-1 cell spontaneously becomes type 0 (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.
- Standard name given: Schelling segregation model (degenerate single-species case after decay of one type)
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- Notes: Type 1 decays irreversibly (gone in ~70-100 steps), leaving a single species at ~1/3 density. Persistent reversals (~0.7, still there at 3x run length) can only come from isolated agents (0/0 same-type share, apparently treated as unhappy) jumping at random forever: an artifact of the zero-division convention. The fall in n_domains vs the knock-out is isolated agents jumping until they land next to another agent (sticky relocation), which removes singletons. No genuine segregation or coarsening; correlation length ~1.
- Citations:
  - https://en.wikipedia.org/wiki/Schelling%27s_model_of_segregation
  - https://mesa.readthedocs.io/stable/examples/basic/schelling.html
  - https://jasss.soc.surrey.ac.uk/18/4/15.html
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:058: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a type-1 cell spontaneously becomes type 0 (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Schelling segregation model (vacancy-move dynamics) relaxing to an absorbing frozen state after one type is removed by one-way decay
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: Type 1 decays irreversibly (about 0.5% left after 50 steps). Then only type 2 and empty cells remain; every type-2 cell with an occupied neighbour has share 1 and never moves, and isolated movers (if share 0 is assigned to no-neighbour cells) only decrease, so the system goes to an absorbing state. g_activity_trend is negative (relaxation). g_reversals beats the knock-out only because vacancy moves necessarily make cells go occupied to empty and back, while the knock-out has no movement; it most likely stays high at long runs because it is a ratio that remains defined at very low activity. Recommend checking the absolute moves-per-step series; a non-vanishing plateau would point to an implementation bug, not new behaviour.
- Citations:
  - https://arxiv.org/pdf/2608.16557
  - https://arxiv.org/pdf/0707.1681
  - https://jasss.soc.surrey.ac.uk/12/1/6/appendixB/Schelling1969.html
  - https://arxiv.org/pdf/0903.4694
  - https://jasss.soc.surrey.ac.uk/15/1/6.html
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:068: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_types
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a type-0 cell becomes type 1 if at least 3 neighbours are type 1 (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Schelling segregation with vacancy depletion by (polluted) threshold-3 bootstrap percolation -> zero-vacancy frozen absorbing state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- Notes: n_types 3->2 is the empty state (type 0) running out: rule 1 only removes empty cells, rule 2 only moves them, and type 2 never changes. Rule 1 alone stalls because type-2 cells act as polluted/closed sites; Schelling moves keep breaking up the stalled empty pockets until all are filled, then the system freezes. The type count includes the empty state, so it is a bookkeeping artifact of reaching an absorbing state, not an emergent pattern. Robust across seeds and scales because it is guaranteed.
- Citations:
  - https://arxiv.org/pdf/1812.06008
  - https://arxiv.org/pdf/0903.4694
  - https://arxiv.org/pdf/0711.2212
  - https://www.researchgate.net/publication/395582600_Three-Neighbour_Bootstrap_Percolation_in_Thin_Three-Dimensional_Grids
  - https://etd.ohiolink.edu/acprod/odb_etd/ws/send_file/send?accession=osu166109554423144&disposition=inline
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

