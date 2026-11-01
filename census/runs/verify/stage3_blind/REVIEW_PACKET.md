# Stage 3 (b): author review packet

Model labels are provisional (VERIFICATION_PLAN.md, Stage 3). For each phenomenon: check that the cited sources exist and describe the behaviour as stated, then record a final label (documented / closely related / undocumented candidate). Order: closely related first, then documented and judged non-trivial, then documented and judged trivial or artifactual.

84 phenomena: 79 documented, 5 closely related.

## census:novelty:network-2: CLOSELY RELATED (KNOWN-VARIANT)

- Class: COMBINATION EFFECT; flagged: False; robust features: n_activity, n_type_modularity
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 2 types (random start). Rules, applied in order each step: (1) nodes play a two-strategy game with their neighbours (T=1.9, S=0.0) and copy the best-scoring neighbour; (2) each node adopts its neighbours' majority type (with probability 0.1 per step).
- Standard name given: Nowak-May spatial (weak) prisoner's dilemma with imitate-the-best on a random graph, plus conformity / majority-rule label dynamics (tag hitchhiking)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 4
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Type modularity is not an interaction effect: rule-1-only gives 0.158 vs 0.167 for the full program, within one robust spread, which implies best-neighbour copying carries the type along (tag hitchhiking). Only the activity (0.29 vs 0.02/0.0) needs both rules. It is plausibly persistent frustration between payoff-driven copying and majority-rule relaxation of the label, which is unremarkable. Twist not found in print: conformity acting on a label separate from the strategy rather than on the strategy. Probably not important. Check in the code whether imitation copies the type.
- Citations:
  - Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344 (Nowak-May weak PD, imitate best)
  - Pacheco et al., Conformity and Network Effects in the Prisoner's Dilemma
  - Evolutionary PD on complex networks with conformist mentality strategy, Sci China PMA 2012
  - Kim, A Tag-Based Evolutionary Prisoner's Dilemma Game on Networks with Different Topologies, JASSS 13(3)2
  - Jensen & Bornholdt, Imitating the winner leads to discrimination in spatial PD, Sci Rep 9:3776 (2019)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:037: CLOSELY RELATED (KNOWN-VARIANT)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_type_entropy, g_reversals, g_interface
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 3 types (random start). Agents have 4 types (random start). Agents move at speed 0.0 with heading noise 1.0 rad per step. Rules, applied in order each step: (1) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.2 moves to a random empty cell; (2) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Schelling segregation with spatially fixed type-switching sites (Hazan & Randon-Furling 2012 'switching agents'); noise-driven unfreezing of a low-threshold Schelling lattice
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Speed 0 makes the agents fixed sites that cycle cell type each step, which is equivalent to Hazan & Randon-Furling's fixed switching agents in a Schelling lattice. The domain and interface reductions are ordinary Schelling clustering, unlocked because the switching sites keep making cells unhappy (at threshold 0.2 with ~1/3 empties the knock-out barely moves). Twist: the coupling does not conserve counts. When an unhappy switching-site cell relocates, the agent refills the vacated site, steadily turning empty cells into occupied ones. This vacancy pumping most likely explains the entropy drop (a change in type counts, not spatial order). Non-zero reversals are built into the rules (fixed sites cycling plus relocations). Two of the four features are therefore partly artifacts of the coupling; the clustering features look genuine.
- Citations:
  - https://arxiv.org/pdf/1212.4929
  - https://arxiv.org/pdf/0903.4694
  - https://arxiv.org/pdf/2204.01417
  - https://arxiv.org/pdf/1504.03809
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:046: CLOSELY RELATED (KNOWN-VARIANT)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (every step); (2) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step).
- Standard name given: Synchronous majority-vote coarsening (Goles-Olivos period-2 blinking) perturbed by a threshold-1 cyclic cellular automaton (Griffeath CCA / rock-paper-scissors invasion at domain boundaries)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 4
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN-VARIANT
- Notes: Reversals (0.94) are mostly the known synchronous-majority period-2 artifact; majority alone gives 0.86. The knock-out removes both rules and leaves a frozen random start, so beating it is weak evidence. The real interaction effect is that cyclic invasion at domain edges roughly doubles clustering (Moran 0.36 -> 0.72) and raises the domain measure over majority alone. No spirals or travelling waves are reported; correlation length is only 3 cells. The exact majority + cyclic combination was not found in a paper, but the behaviour is a quantitative variant of known behaviours, not a new mechanism.
- Citations:
  - Goles & Olivos 1980, Periodic behaviour of generalized threshold functions, Discrete Math 30:187 (synchronous majority -> fixed points or 2-cycles)
  - Poljak & Sura (multi-colour majority period <=2), as cited in Biased Majority Rule Cellular Automata
  - Fisch, Gravner & Griffeath 1991, Cyclic Cellular Automata in Two Dimensions
  - Shalizi et al., Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067 (low threshold -> global periodic oscillation, near-random texture)
  - MCell Cyclic CA lexicon, psoup.math.wisc.edu (R1/T3/C3/NM is the only self-organizing 3-colour range-1 CCA)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:050: CLOSELY RELATED (KNOWN-VARIANT)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Cyclic cellular automaton (Griffeath; K=4, threshold 1, 8-neighbour) waves/spirals coupled to Schelling segregation moves; similar to cyclic dominance with mobility (Reichenbach-Mobilia-Frey)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN-VARIANT
- Notes: Rule 1 alone is the standard cyclic automaton: every cell changes every step, as in the known wave/demon-dominated state. Rule 2 alone is the standard Schelling model: it segregates and freezes. Combined, the moves partly break up the waves (activity about 0.61) and push up same-colour clustering (Moran about 0.32), the same kind of in-between regime as mobility breaking up spirals in rock-paper-scissors models. Twist not found in the literature: the empty state is part of the colour cycle, and moves are long-range jumps triggered by being unhappy with neighbours rather than local swaps. Likely changes the numbers, not the type of behaviour. The domain-count effect is small and the weakest of the three features.
- Citations:
  - Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata, Stat. Comput. 1 (1991) 23-39
  - Durrett & Griffeath, Asymptotic behavior of excitable cellular automata, Exp. Math. 2 (1993)
  - Reichenbach, Mobilia & Frey, Nature 448 (2007) 1046-1049
  - Szczesny, Mobilia & Rucklidge, The Influence of Mobility Rate on Spiral Waves in Spatial Rock-Paper-Scissors Games, Games 7 (2016) 24
  - Schelling (1971) segregation model with random relocation
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:051: CLOSELY RELATED (KNOWN-VARIANT)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity, g_interface
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell; (2) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step).
- Standard name given: Schelling segregation coupled to a Griffeath cyclic cellular automaton (κ=4, T=1, Moore range 1; locally periodic 'quench' regime) — cf. mobility in cyclic-dominance/RPS lattice models and Schelling variants with switching agents
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN-VARIANT
- Notes: Rule 2 alone reproduces the documented CCA T=1 Moore 4-colour locally periodic state (activity=1.0). The full program's domain/interface values are close to the Schelling-only sub-program; activity is the CCA cycling partly suppressed (~0.67), probably because relocations and the 3->empty step insert defects. The knock-out baseline is a frozen random state, so beating it is weak evidence; the sub-program comparison suggests superposition with mild interference. Twist not found in the literature: persistent cycling of Schelling-coarsened domains at intermediate activity. Probably not significant unless images or spiral/phase-defect statistics show new structure.
- Citations:
  - https://arxiv.org/html/nlin/0507067
  - http://psoup.math.wisc.edu/mcell/rullex_cycl.html
  - https://arxiv.org/pdf/1807.05419
  - https://arxiv.org/pdf/1212.4929
  - https://www.sciencedirect.com/science/article/abs/pii/S0378437117306118
  - https://arxiv.org/pdf/1608.08932
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:agents-1: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: a_density_cv, a_step, a_nn_ratio, a_local_align
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. Rules, applied in order each step: (1) each agent turns toward the centre of any type agents within distance 1.0.
- Standard name given: Local Attraction Model (Strömbom 2011) — attraction-only self-propelled particles, swarm/aggregation regime
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule is the Local Attraction Model (steer toward local centre of mass within R, constant speed, noise). Strong clustering, near-stationary aggregates (low net displacement) and low or negative local alignment match the documented 'swarm' regime at strong attraction. Negative nearest-neighbour alignment is not stated explicitly in the sources read, but follows geometrically from agents crossing through the local centre. Overturns 'variant' to 'known' because there is no meaningful twist in the phenomenology. Suggested check: measure rotation order to rule out mills.
- Citations:
  - Strömbom D. (2011) Collective motion from local attraction. J Theor Biol 283:145-151, doi:10.1016/j.jtbi.2011.05.019
  - Strömbom D., Hassan T., Greis W.H., Antia A. (2019) Asynchrony induces polarization in attraction-based models of collective motion. R Soc Open Sci 6:190381
  - Strömbom D., Antia A. Anticipation induces polarized collective motion in attraction based models, arXiv:1710.05692
  - Strömbom D., Siljestam M., Park J., Sumpter D.J.T. (2015) The shape and dynamics of local attraction. Eur Phys J ST 224:3311-3323
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:agents-3: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: a_density_cv, a_local_align
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 3.14159 rad per step. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 1.0.
- Standard name given: Vicsek model (1995): short-range orientational order and alignment-induced clustering/density fluctuations near the order-disorder transition
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The rule set is the standard Vicsek model (unit radius, constant speed, angular noise, periodic box), not a variant, so 'VARIANT' is a small overturn to 'KNOWN'. A uniform noise width of pi (eta=pi) is consistent with the nonzero local alignment; full noise over [-pi,pi] would erase alignment, so check the noise convention in the code. With density about 1 and eta=pi the system is near the transition: weak local alignment (0.17) plus higher density CV matches documented short-range correlations and clustering/bands. The effect is robust against the knock-out, so not an artifact, but it is weak order, not strong flocking.
- Citations:
  - Vicsek, Czirok, Ben-Jacob, Cohen, Shochet, 'Novel Type of Phase Transition in a System of Self-Driven Particles', PRL 75, 1226 (1995)
  - Chatterjee, 'The Flocking Transition: A Review of the Vicsek Model' (UCSD PHYS 563 essay, 2017)
  - Active-learning mapping of the Vicsek model phase diagram, arXiv:2604.28167
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:agents-17: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: a_step, a_density_cv, a_nn_ratio, a_local_align
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns away from the centre of any type agents within distance 2.0.
- Standard name given: Autochemotactic (Keller-Segel) aggregation of self-propelled particles with short-range repulsion limiting chemotactic collapse (finite/dynamic clustering)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Emission into field 0 is not stated in the rule but is implied by the rules-1-only near-total collapse (density_cv ~11, nn_ratio ~0.04), i.e. the documented chemotactic collapse. Adding short-range turn-away repulsion stops this at moderate clustering, as documented for chemotactic walkers with steric repulsion. Reduced step length follows from agents being trapped in clusters. Weak negative local alignment is a trivial side effect of turn-away repulsion (rule 2 alone gives -0.06) and is not new order. Robust across seeds, 4x area and 3x run length, so not an artifact. Only abstracts/snippets were checked, not full texts; confidence medium because the emission term is inferred.
- Citations:
  - Taktikos, Zaburdaev & Stark, Collective dynamics of model microorganisms with chemotactic signaling, Phys. Rev. E 85, 051901 (2012) https://pubmed.ncbi.nlm.nih.gov/23004782/
  - Pohl & Stark, Dynamic clustering and chemotactic collapse of self-phoretic active particles, arXiv:1403.4063 https://arxiv.org/html/1403.4063
  - Stark, Artificial Chemotaxis of Self-Phoretic Active Colloids: Collective Behavior, Acc. Chem. Res. 51, 2681 (2018) https://pubmed.ncbi.nlm.nih.gov/30346724/
  - Liebchen et al., Synthetic Chemotaxis and Collective Behavior in Active Matter, arXiv:1810.01117
  - Phase separation of chemokinetic active particles, Commun. Phys. (2025) https://www.nature.com/articles/s42005-025-02211-0
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:agents-8: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: a_density_cv, a_step, a_local_align
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns to the average heading of agents within distance 2.0.
- Standard name given: Autochemotactic (Keller-Segel-type) collapse of self-propelled particles with Vicsek alignment
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 1 alone reproduces the clumping and low step; the knock-out implies agents emit field 0 (autochemotaxis), although no emission rule is listed. High local alignment inside the collapsed clumps follows almost automatically from Vicsek averaging once agents sit within radius 2; low step reflects agents trapped by their own chemical peak. Genuine emergent collapse, but close to an absorbing state; the alignment score adds little beyond that geometry. Overturns the first pass mildly (VARIANT -> KNOWN): no phenomenological twist beyond documented collapse plus alignment. Verified from search excerpts only.
- Citations:
  - Pohl & Stark, Dynamic clustering and chemotactic collapse of self-phoretic active particles, arXiv:1403.4063 / PRL 2014
  - Saha, Golestanian & Ramaswamy, Clusters, asters and collective oscillations in chemotactic colloids, arXiv:1309.4947
  - Boltzmann-Ginzburg-Landau theory for autochemotaxis of active rod-like particles, arXiv:2506.18263
  - Stark, Artificial Chemotaxis of Self-Phoretic Active Colloids: Collective Behavior, PubMed 30346724
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:012: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_n_domains, g_reversals, g_activity
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) neighbouring cells swap contents (every step).
- Standard name given: Spatial snowdrift (hawk-dove) game with deterministic best-neighbour imitation (Nowak-May lattice dynamics) plus random mobility/strategy diffusion
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. The swap rule alone produces most of the activity (0.33 of 0.45), and g_reversals is saturated (1.0) in both sub-programs, so it carries no information. The genuine interaction effect is modest: mobility keeps the otherwise near-frozen imitation lattice churning, and imitation slightly coarsens domains (5.2 vs 6.4). Both are documented mobility-plus-imitation phenomenology. Not a numerical artifact, but largely driven by the non-game mixing rule. Upgraded from KNOWN-VARIANT because there is no phenomenological twist beyond the parameter choice.
- Citations:
  - Nowak & May, Evolutionary games and spatial chaos, Nature 359:826-829 (1992)
  - Sicardi, Fort, Vainstein & Arenzon, Random mobility and spatial structure often enhance cooperation, J. Theor. Biol. 256:240-246 (2009)
  - Vainstein, Silva & Arenzon, Does mobility decrease cooperation?, J. Theor. Biol. 244:722-728 (2007)
  - Vainstein & Arenzon, Spatial social dilemmas: dilution, mobility and grouping effects with imitation dynamics, Physica A 394 (2014)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:013: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_type_entropy
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) each cell adopts its neighbourhood's majority type (with probability 0.1 per step).
- Standard name given: Nowak-May spatial prisoner's dilemma (imitate-best, Moore neighbourhood + self) with conformity/majority-rule updating (Szolnoki & Perc conformity-enhanced network reciprocity)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule-level differences (all cells update at once, every cell conforms with probability 0.1 each step instead of a fixed conformist fraction, S=-0.5) do not change the behaviour, which is documented. Mild overturn: VARIANT -> KNOWN, because no new behaviour was identified. The effects are real but small (correlation length 2 vs 1; reversals looks like a binary flag). Full papers were not opened, so parameter regimes were not checked.
- Citations:
  - Nowak & May 1992, Nature 359:826 (spatial chaos, C/D coexistence in fluctuating clusters)
  - Szolnoki & Perc 2015, J. R. Soc. Interface 12:20141299, arXiv:1412.4113 (conformity in spatial social dilemmas; surface tension, smooth domain interfaces)
  - Szolnoki & Perc 2016, Sci. Rep. 6:23633 (conformists; strategy-neutral state at full conformity)
  - Chaos Solitons Fractals 2018 S0960077918305629 (time-alternating payoff-driven/conformity-driven updating)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:015: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: n_oscillation, n_edge_turnover
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 2 types (random start). Rules, applied in order each step: (1) nodes play a two-strategy game with their neighbours (T=1.9, S=0.0) and copy the best-scoring neighbour; (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (with probability 0.003 per step).
- Standard name given: Coevolutionary (adaptive-network) Prisoner's Dilemma with copy-the-best-neighbour updating and partner switching (Zimmermann-Eguiluz / Nowak-May lineage)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 4
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Overturns VARIANT to KNOWN: no twist in the behaviour itself was found. Oscillation is the known non-settling behaviour of synchronous copy-the-best weak PD at T=1.9; it is already present with the game alone (0.127), and rewiring raises it modestly (0.33). Edge turnover nearly reads out rule 2 directly (rule 2 alone gives 0.38, more than the full model's 0.34), so it is close to trivial rather than emergent. The knock-out test cannot show emergence because both features are interactions by definition. Unclear point: whether 'type' is the strategy (then it is the documented mismatched-link rewiring model) or a fixed label (then it is two documented processes running side by side). Oscillation is a real dynamical pattern, not a numerical artifact, but it is documented.
- Citations:
  - Zimmermann & Eguiluz, Phys. Rev. E 72, 056118 (2005)
  - Zimmermann, Eguiluz & San Miguel, Phys. Rev. E 69, 065102 (2004)
  - Abramson & Kuperman, Phys. Rev. E 63, 030901 (2001)
  - Perc & Szolnoki, BioSystems 99, 109 (2010)
  - Evolutionary prisoner's dilemma games coevolving on adaptive networks, J. Complex Networks 6, 1 (2018), arXiv:1702.05119
  - Szabo & Fath, Phys. Rep. 446, 97 (2007)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:027: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: a_density_cv, a_local_align, a_step, a_polar_std
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.02]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns to the average heading of agents within distance 1.0.
- Standard name given: Vicsek alignment + autochemotaxis (Keller-Segel-type chemotactic collapse into aligned, trapped aggregates / 'active droplets')
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The emitter is not stated in the rules; agent emission is inferred because rule 1 alone produces strong clumping. The features describe chemotactic collapse into dense, internally aligned, nearly stationary clumps; rule 2 alone gives Vicsek bands. Density CV keeps rising at 3x run length, so this is ongoing collapse/coarsening rather than a steady state. High local alignment is partly automatic once agents sit within one alignment radius. Real emergent clustering, not an artifact, but close to a degenerate collapsed state. Overturns 'variant' only in that no extra phenomenological twist was found.
- Citations:
  - https://www.nature.com/articles/s41467-022-34484-2
  - https://arxiv.org/html/2506.18263
  - https://www.sciencedirect.com/science/article/pii/S0006349518312268
  - https://arxiv.org/html/2011.04580v1
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:033: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: a_density_cv, a_nn_ratio, a_local_align
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 1 type (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent turns left on an even-type cell and right on an odd-type cell.
- Standard name given: Chemotactic (autochemotactic / Keller-Segel-type) aggregation of active Brownian particles, softened by a fixed random torque field acting as extra direction noise, with emergent alignment from streaming toward aggregates
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 2 alone matches the knock-out (inert fixed random torque = extra direction noise). Rule 1 alone gives chemotactic collapse (nn_ratio 0.009). The full program sits in between: noise-limited aggregation, with local alignment from agents streaming toward diffuse aggregates, which is documented in agent-based chemotaxis models. The field's source is not stated; if agents emit, this is standard autochemotaxis, and if the field is sourceless, it is gradient-climbing on a landscape. Either way the behaviour is known. The overturn is mild (VARIANT -> KNOWN) because no phenomenological twist was found.
- Citations:
  - https://arxiv.org/html/2409.04262v1
  - https://arxiv.org/html/1403.4063v1
  - https://www.sciencedirect.com/science/article/pii/S0010465526002675
  - https://www.sciencedirect.com/science/article/pii/S0006349518312268
  - https://arxiv.org/pdf/2303.07746
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:039: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_n_domains, g_activity_cv, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) a type-0 cell becomes type 1 if at least 4 neighbours are type 1 (with probability 0.03 per step).
- Standard name given: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore+self) with one-sided C->D noise (cf. SPD(P) stochastic variant): near-absorbing dominance with noise-seeded transient clusters
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Overturns 'variant' only mildly: the twist (C->D conversion gated by >=4 defector neighbours) changes where noise acts (interfaces only), not the kind of behaviour. Features are largely consequences of near-consensus (largest domain ~0.97); g_reversals appears binary. Knock-out removes both rules, leaving a frozen random field, so knock-out contrasts are uninformative. Strict-PD S=-0.5 assumed qualitatively like the weak PD per Nowak-May; the dominant type was not determined from the summaries.
- Citations:
  - Nowak & May 1992, Nature 359:826 (as summarised in Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344)
  - Dynamic phase transition in Prisoner's Dilemma on a lattice with stochastic modifications, arXiv 0910.0955 (SPD(P): cooperators defect with prob p under unconditional imitation)
  - Szabo, Vukov, Szolnoki 2005, Phase diagrams for PD on 2D lattices, arXiv cond-mat/0506433
  - Helbing & Yu, The outbreak of cooperation among success-driven individuals under noisy conditions, arXiv 0903.4054
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:040: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour; (2) a type-1 cell becomes type 0 if at least 2 neighbours are type 0 (with probability 0.03 per step).
- Standard name given: Nowak-May spatial Prisoner's Dilemma (weak PD, b=1.3, Moore neighbourhood + self-interaction, synchronous imitate-the-best): cooperator-dominated state with stable defector lines and blinkers
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 1 alone is the textbook Nowak-May model and reproduces every feature in the sub-program test. Rule 2 (D->C conversion with probability 0.03 given at least 2 C neighbours) is a one-way bias toward cooperators that only shifts the magnitudes slightly. The first pass called this a variant; I upgrade it to plain KNOWN, a mild overturn. The pattern is genuine but near-static: a cooperator sea with defector lines and documented blinkers that account for the reversals. It is not an artifact, since the knock-out is just the random initial colouring.
- Citations:
  - Nowak & May 1992, Evolutionary games and spatial chaos, Nature 359:826 (https://www.researchgate.net/publication/216634494_Evolutionary_Games_and_Spatial_Chaos)
  - Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344, Sec. 6.5
  - Newth & Cornforth, Asynchronous spatial evolutionary games: spatial patterns, diversity and chaos, IEEE CEC 2007 (b=1.3: stable defector lines + blinkers)
  - Sigmund, Univie notes on Nowak-May blinkers/gliders
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:044: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell copies the type of a random neighbour (with probability 0.3 per step); (2) each cell adopts its neighbourhood's majority type (with probability 0.2 per step); (3) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.
- Standard name given: Spatial snowdrift/hawk-dove game with imitate-best update (Nowak-May type) plus voter noise and conformity: near-dominance of one strategy with noise-maintained small minority clusters
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Downgraded from first-pass 'variant' to KNOWN: the only novelty is which rules are combined, not the behaviour. Rules 2+3 alone lock into a single type; voter noise (rule 1) keeps flipping cells and maintains small minority islands. Correlation length is short (2 cells), not long. The knock-out baseline (pure random noise) is weak, so beating it on clustering features is expected for any copy-a-neighbour or majority dynamics. Not an artifact, but a low-interest near-absorbing state kept alive by noise.
- Citations:
  - Nowak & May 1992, Nature 359:826 (Evolutionary games and spatial chaos)
  - Hauert & Doebeli 2004, Nature 428:643 (Spatial structure often inhibits the evolution of cooperation in the snowdrift game)
  - Szolnoki & Perc 2015, J. R. Soc. Interface 12:20141299 (Conformity enhances network reciprocity in evolutionary social dilemmas)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:045: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_n_domains, g_corrlen, g_n_types, g_reversals, g_type_entropy, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) each cell copies the type of a random neighbour (with probability 0.1 per step); (2) a type-0 cell becomes type 2 where field 0 exceeds 0.5.
- Standard name given: Multi-state voter model coarsening on a 2D lattice (logarithmically slow, pre-consensus transient), with one type removed by a one-time threshold conversion
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 1 is the voter model and alone reproduces the clustering, correlation length, reversals and domain count. Rule 2 permanently converts all type-0 cells to type 2 (rule-2-only entropy 0.753 matches frequencies (0,0.25,0.5,0.25)); since the voter rule cannot recreate an extinct type, this alone explains the 3 surviving types and the lower entropy. No emitters are described, so the fields presumably only diffuse: field 1 is inert and field 0 is a passive threshold. The knock-out removes the voter rule itself, so its large differences reflect the voter model's presence, not an interaction between the rules. Overturns the first pass in that this is not a meaningful variant: the spatial pattern is genuine but textbook voter coarsening, and the coupling adds only a trivial recolouring.
- Citations:
  - Dornic, Chaté, Chave, Hinrichsen, Critical coarsening without surface tension: the universality class of the voter model, PRL 87, 045701 (2001)
  - Percolation and coarsening in the bidimensional voter model, arXiv:1506.05321
  - Ben-Naim, Frachebourg, Krapivsky, Coarsening and persistence in the voter model, PRE (cnls.lanl.gov/People/ebn/pubs/voter/voter-pub.pdf)
  - Dall'Asta & Galla, Algebraic coarsening in voter models with intermediate states, arXiv:0806.0817
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:050: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity_trend, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) a type-0 cell becomes type 3 where field 1 exceeds 0.05; (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Schelling segregation model (multi-type, random relocation, Moore neighbourhood) with field-gated vacancy filling acting as a vacancy-density/composition shift
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Downgraded from KNOWN-VARIANT to KNOWN: no new behaviour. Rule 2 is the Schelling rule. Rule 1 is a one-way field threshold with no emission, so the field just smooths out and its only effect is filling vacancies with type-3 cells, i.e. a shift in vacancy density and type mix, both standard Schelling parameters. Fewer domains = segregation/coarsening. The 'activity rising' value is negative, so activity is falling, which is normal relaxation (feature mislabelled). Flip-backs = relocation churn of unsatisfied agents. The field part is nearly inert; the segregation domains are genuine.
- Citations:
  - https://arxiv.org/pdf/2608.16557
  - https://arxiv.org/html/0711.2212v1
  - https://arxiv.org/pdf/0903.4694
  - https://www.nature.com/articles/s41598-022-23224-7
  - https://arxiv.org/pdf/2509.14462
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:057: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 3 types (random start). Agents have 1 type (random start). Agents move at speed 0.5 with heading noise 0.0 rad per step. Rules, applied in order each step: (1) each agent turns away from the centre of any type agents within distance 1.0; (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.4 moves to a random empty cell; (3) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Noise-driven coarsening in the Schelling segregation model (noisy / Glauber-Kawasaki Schelling; the zero-noise model freezes, small noise gives large domains)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Agent repulsion (rule 1) is inert: rules 1+2 equal rule 2, and rules 2+3 alone reproduce or exceed the effect. The mechanism is Schelling relocation (rule 2) plus non-conservative type-cycling noise from walkers (rule 3), which removes the zero-noise freezing and lets domains coarsen; this is documented. The noise source (spatially correlated, one-directional cycling) differs in the rules but no new phenomenology was identified. g_reversals is largely a trivial consequence of walkers rewriting cells. Check whether the large domain is the empty type (type 0); if so it may be a vacancy artifact rather than segregation.
- Citations:
  - Dall'Asta, Castellano, Marsili, J. Stat. Mech. (2008) L07002, arXiv:0707.1681
  - Stauffer & Solomon / Schelling-Ising model, arXiv:0710.5237
  - Barmpalias et al., Unperturbed Schelling segregation in two or three dimensions, arXiv:1504.03809
  - Hydrodynamic limit of the Schelling model with spontaneous Glauber and Kawasaki dynamics, arXiv:2302.09866
  - Aging effects in Schelling segregation model (noisy constrained Schelling, power-law coarsening)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:060: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity_trend, g_corrlen, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 3.14159 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent advances the type of the cell under it by one (cyclically); (3) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.
- Standard name given: Nowak-May spatial game with imitate-the-best on an 8+1 Moore lattice (snowdrift/hawk-dove payoffs, T>R>S>P), plus agent-driven local strategy flips acting as noise
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 3 alone reproduces largest_domain (0.82 vs 0.86), n_domains, corrlen and reversals. Rule 1 is inert for the cells (sub-program 1,3 is identical to 3), likely because the field has no source. Rule 2 is localized random flipping, i.e. noise. The negative activity trend comes mainly from rules 1+2 together and is likely a transient from the decaying initial field. The knock-out baseline (frozen random lattice) inflates apparent emergence. Feature labels are misleading: activity is falling, not rising, and corrlen is about 2 cells. Partly artifactual: near-uniform end state with an inert field. No new phenomenology found.
- Citations:
  - Nowak & May 1992, Nature 359:826 (spatial games, imitate-the-best)
  - Hauert & Doebeli 2004, Nature 428:643, doi:10.1038/nature02360 (spatial snowdrift game)
  - Grujic et al. 2010 PLoS ONE 5:e13749 (imitate-the-best goes to cooperation near 0 or near 1 depending on initial condition)
  - Cuesta et al. arXiv:0805.2071 (Nowak-May 8+1 neighbourhood)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:071: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_activity_trend, g_interface
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) field 1 grows autocatalytically by consuming field 0 (rate 0.5 x f x g^2); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.2 moves to a random empty cell; (3) type-2 cells release 1.0 of field 1 per step.
- Standard name given: Schelling segregation model (multi-type, Moore neighbourhood, vacancies) with coarsening / decreasing interface length and relaxing move activity
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Both features come only from rule 2 (Schelling relocation with a 0.2 threshold); sub-programs without rule 2 equal the knock-out exactly. The fields are one-way: no rule reads them, so the reaction-diffusion part is inert with respect to the cells. The gap between the full program and rule 2 alone is probably caused by the activity measure counting field changes, or by RNG-stream effects, not by real coupling. Calling it a 'variant' gives the inert field too much credit; the behaviour is plain Schelling coarsening.
- Citations:
  - https://www.demographic-research.org/articles/volume/21/12
  - https://www.nature.com/articles/s41598-022-23224-7
  - https://arxiv.org/pdf/1505.00429
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:008: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_corrlen, g_n_domains, g_moran, g_cyclic_changes
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step).
- Standard name given: Cyclic cellular automaton (Fisch-Gravner-Griffeath), rule '313' (range 1, threshold 3, 3 colours, Moore), stochastic/asynchronous update
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule is exactly Griffeath's named CCA rule 313 (R1/T3/C3/NM), documented as the only 3-colour range-1 CCA rule with significant self-organization (spirals). Only difference is the 0.1 flip probability; random versions of CCA/GH are discussed by FGG and Griffeath (stochastic wave fronts). Clustering, domain coarsening, finite correlation length and above-chance cyclic changes are the standard CCA signatures. The knock-out is trivially frozen because there is no non-interacting rule, but the interacting system keeps cycling and the features are robust to size and run length. The overturn is mild: the first pass said 'variant', but the behaviour is not a variant in phenomenology. On 64x64 the data can't tell spirals from turbulent waves; both are documented.
- Citations:
  - http://psoup.math.wisc.edu/mcell/rullex_cycl.html
  - https://link.springer.com/article/10.1007/BF01890834
  - https://arxiv.org/pdf/nlin/0507067
  - https://www.stat.cmu.edu/~cshalizi/462/readings/Self-Organization%20of%20Random%20Cellular%20Automata-%20Four%20Snapshots%20(Griffeath).pdf
  - https://link.springer.com/chapter/10.1007/978-1-4612-0451-0_8
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:009: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: p_r, p_freq_spread
- Rule: N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=0.1  (compact grammar notation)
- Standard name given: Locally (nearest-neighbour) coupled Kuramoto lattice with stochastic on-off (blinking) coupling: frequency entrainment with partial, size-dependent, slowly coarsening phase order
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Overturned from KNOWN-VARIANT to KNOWN: random gating with p=0.1 is just a blinking network, and in the fast-switching limit it behaves like a lattice Kuramoto model with effective coupling of about pK. No new phenomenology. The p_freq_spread metric must be measuring effective frequencies, not natural ones (natural frequencies can't change under coupling); check its label. Partial r that falls with area and rises with run length matches documented finite-size effects and vortex coarsening in the 2D lattice Kuramoto model.
- Citations:
  - Hong, Park & Choi, PRE 72, 036217 (2005) - collective synchronization in spatially extended oscillator lattices
  - Noise-induced synchronization in the Kuramoto model on finite 2D lattice, arXiv:2004.00294
  - L.-H. Tang, finite-size scaling in the Kuramoto model, arXiv:1011.5823
  - Frequency Spirals in the 2D Kuramoto Lattice, arXiv:1512.02335
  - Jeter & Belykh, Synchronization in On-Off Stochastic Networks, IEEE TCAS-I 62, 1260 (2015)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:012: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_interface, g_activity
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (with probability 0.3 per step).
- Standard name given: Cyclic cellular automaton (Griffeath CCA, 4 colours, 8-cell neighbourhood, threshold 1) with random update, i.e. a discrete-time form of the Bramson-Griffeath cyclic particle system
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Shalizi & Shalizi describe exactly these settings (4 colours, 8-cell neighbourhood, threshold 1): every cell keeps cycling, and the field shows only a slight departure from randomness as twisting single-colour bands. This matches both measured features (interface 0.65 vs 0.75 random; sustained activity). Bramson-Griffeath / Fisch results say the 4-colour case keeps fluctuating rather than freezing. The p=0.3 random update is a rule detail, not new behaviour. The pattern is weak; the activity knock-out is trivially frozen; an activity of 0.417 above p=0.3 points to one-cell-at-a-time updating or a multi-step measure and should be checked.
- Citations:
  - https://arxiv.org/pdf/nlin/0507067
  - https://arxiv.org/html/1711.04741
  - https://en.wikipedia.org/wiki/Cyclic_cellular_automaton
  - https://www.stat.cmu.edu/~cshalizi/462/readings/Cyclic%20Cellular%20Automata%20in%20Two%20Dimensions%20(Fisch,%20Gravner%20and%20Griffeath).pdf
  - https://commons.wikimedia.org/wiki/Category:Cyclic_cellular_automata
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:015: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: n_cyclic_changes
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) a node of type t becomes t+1 (cyclically) if at least 2 neighbours are t+1 (with probability 0.1 per step).
- Standard name given: Threshold cyclic cellular automaton (Griffeath CCA / Greenberg-Hastings-type excitable medium), asynchronous update on a random graph
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Overturns 'variant' toward plain KNOWN: a random graph plus asynchronous updates adds no new phenomenology in the measured feature. The 'cyclic changes' metric is cyclic by construction (t->t+1 is the only allowed move), so it only measures persistent activity versus fixation. The knock-out value of -0.5 looks like a placeholder for zero activity, so the contrast is active versus completely frozen. Persistent cycling at theta=2, 3 states, mean degree 4 fits the known persistent regime: rotating cycles form by chance at low threshold. Not a numerical artifact, but the measured feature is close to trivial.
- Citations:
  - Fisch, Gravner, Griffeath (1991) Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39
  - Fisch, Gravner, Griffeath (1991) Cyclic cellular automata in two dimensions, Spatial Stochastic Processes
  - Limiting behavior of 3-color excitable media on arbitrary graphs, arXiv:1610.07320
  - Durrett & Griffeath, Asymptotic behavior of excitable cellular automata / Medley of spirals from cyclic cellular automata
  - https://en.wikipedia.org/wiki/Cyclic_cellular_automaton
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:016: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). 1 field(s) with diffusion [0.2]. Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.9, S=0.0) and copy the best-scoring neighbour.
- Standard name given: Nowak-May spatial Prisoner's Dilemma (weak PD, Moore neighbourhood + self, imitate-best), spatial chaos regime 1.8<b<2
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Exact Nowak-May parameters (R=1, S=P=0, T=1.9, 8 neighbours + self, deterministic imitate-best). Largest domain ~0.68 matches 1-0.318 (asymptotic defector fraction), i.e. the defector background, not one type taking over. Reversals and domain churn = documented spatiotemporal chaos. Correlation length 2 is short, not long-range. The diffusing field is inert (read by no rule). The knock-out (frozen random start) is a trivial baseline. Small 64x64 tori may enter periodic cycles (Wu et al.), which doesn't change the verdict. Upgraded from KNOWN-VARIANT to KNOWN because the system and parameters are identical to the original model.
- Citations:
  - Nowak MA & May RM (1992) Evolutionary games and spatial chaos, Nature 359:826-829, doi:10.1038/359826a0
  - Nowak MA & May RM (1993) The spatial dilemmas of evolution, Int J Bifurcation Chaos 3:35-78
  - Nowak & Sigmund (2000) Games on grids, in The Geometry of Ecological Interactions
  - Wu, Fu & Wang (2023) Evolutionary games and spatial periodicity (arXiv:2209.08267)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:017: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_reversals
- Rule: C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.5, 0.5)) p=0.1  (compact grammar notation)
- Standard name given: Nowak-May spatial Prisoner's Dilemma (imitate-the-best-neighbour), with partially asynchronous / stochastic updating; dynamic C/D coexistence ('spatial chaos')
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Downgrades the first pass from KNOWN-VARIANT to KNOWN: two-strategy imitate-best with update probability 0.1 is the Nowak-May spatial game with partial asynchrony, and persistent strategy flipping (dynamic coexistence) is its signature result. The g_reversals feature is weak evidence: the knock-out removes the only rule, so a score of 0 there is guaranteed, and a yes/no flip measure cannot tell real turnover apart from small blinking spots or flips caused by random tie-breaking. Interpreting g=(1.5,0.5) as PD payoffs (b=1.5) is an assumption; the verdict holds for any two-strategy imitate-best lattice game that has ongoing turnover.
- Citations:
  - Nowak & May 1992, Nature 359:826 'Evolutionary games and spatial chaos'
  - Huberman & Glance 1993, PNAS 90:7716 'Evolutionary games and computer simulations'
  - Szabo & Fath 2007, Phys. Rep. 446:97 'Evolutionary games on graphs' (arXiv cond-mat/0607344)
  - Grilo & Correia, 'Asynchronous stochastic dynamics and the spatial Prisoner's Dilemma game' / 'Effects of asynchronism on evolutionary games' JTB 2011
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:018: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: n_edge_turnover, n_type_modularity, n_giant
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) each node copies the type of a random neighbour (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (every step).
- Standard name given: Adaptive/coevolving voter model (multistate), fragmented phase / fragmentation transition
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Both rules fire every step, so rewiring dominates and the system lands deep in the known fragmented phase: like-type, internally uniform components, giant component ~1/3 for 3 types. Rewiring alone (rule 2) already gives nearly the same modularity and giant fraction, so most of the structure is homophilic rewiring on random types. Copying mainly raises turnover. The late state is probably a frozen absorbing fragmented configuration, and turnover may be a cumulative count. Genuine but standard behaviour, not an artifact. Tightens 'variant' to 'known'.
- Citations:
  - Holme & Newman, Phys. Rev. E 74, 056108 (2006), arXiv:physics/0603023
  - Böhme & Gross, Phys. Rev. E 85, 066117 (2012), Fragmentation transitions in multistate voter models
  - Böhme & Gross, Phys. Rev. E 83, 035101 (2011), arXiv:1012.1213
  - Vazquez, Eguíluz & San Miguel, PRL 100, 108702 (2008)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:019: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: n_edge_turnover, n_consensus_gain
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) each node copies the type of a random neighbour (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (with probability 0.1 per step).
- Standard name given: Coevolving / adaptive voter model (Holme-Newman rewire-to-same), 3-opinion version, consensus side of the fragmentation transition
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Low rewiring rate (0.1 vs copying every step) puts the system below the fragmentation threshold, so it reaches global consensus. The consensus gain comes almost entirely from rule 1 alone (0.625 vs 0.635), i.e. plain voter-model absorption; edge turnover follows directly from rule 2. No synergistic emergent feature. Synchronous updating and separate rewiring step are implementation details that don't change the known phenomenology. Not an artifact, but low novelty; the final state is the documented absorbing consensus state.
- Citations:
  - Holme & Newman, Phys. Rev. E 74, 056108 (2006)
  - Durrett et al., PNAS 2012, Graph fission in an evolving voter model
  - Shi, Mucha, Durrett, Multi-opinion coevolving voter model with infinitely many phase transitions (PMC5131864)
  - Böhme & Gross 2011 (three-opinion adaptive voter model, as cited in PMC5131864)
  - Vazquez, Eguíluz, San Miguel PRL 2008
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:021: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_activity
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell copies the type of a random neighbour (every step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.
- Standard name given: Spatial snowdrift (hawk-dove) game with synchronous imitate-the-best updating plus voter-model imitation noise (voter-model perturbation / noisy Nowak-May-type lattice game)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. The fine-grained small domains (corrlen ~2) fit known filament-like anti-coordinated snowdrift patterns. The voter step is standard imitation noise. Activity and reversals come from the voter step disrupting the anti-coordinated state and best-imitation restoring it, which is ordinary noisy steady-state behaviour. The 'long correlation length' is only 2 cells. The no-interaction knock-out is empty, since both rules read neighbours, so the comparison with each rule alone is the meaningful one. Overturned in the direction of KNOWN rather than VARIANT: no twist in the behaviour was identified. Genuine interplay of the two rules but low novelty.
- Citations:
  - Hauert & Doebeli, Nature 428, 643-646 (2004), Spatial structure often inhibits the evolution of cooperation in the snowdrift game
  - Szabó & Fáth, Physics Reports 446, 97-216 (2007), Evolutionary games on graphs
  - Cox, Durrett & Perkins, Voter Model Perturbations and Reaction Diffusion Equations, AMS (2013)
  - arXiv:1102.4570 (summary of filament-like snowdrift clusters / anti-coordination)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:023: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_activity_trend
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.1 per step).
- Standard name given: Spatial snowdrift (hawk-dove) game with imitate-the-best updating (Hauert-Doebeli parameterisation, r=0.5) plus voter-model noise / voter-model perturbation; coarsening toward single-strategy dominance
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The payoffs T=1.5,R=1,S=0.5,P=0 make this a snowdrift game, exactly on the Hauert-Doebeli line with r=0.5, not a prisoner's dilemma. Adding 10% neighbour copying is voter-model noise. The observed coarsening, ~94% single-type dominance, falling activity and noise-driven flips are what this documented system is expected to do. Moved from VARIANT to KNOWN because I found no new behaviour. The pattern is not a numerical artifact but is close to trivial: it is drift toward a near-absorbing single-type state, and the knock-out baseline is just a frozen random colouring. The feature labelled 'activity rising' is actually negative (activity falls). I did not check which strategy wins or how fast fixation happens.
- Citations:
  - Hauert & Doebeli 2004, Nature 428:643, doi:10.1038/nature02360
  - Cox, Durrett & Perkins, Voter model perturbations and reaction diffusion equations, Asterisque 349 (arXiv:1103.1676)
  - Durrett 2014, Spatial evolutionary games with small selection coefficients, EJP 19 (arXiv:1406.5876)
  - The effect of dispersal and neighbourhood in games of cooperation, J. Theor. Biol. (S0022519308001112)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:027: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) neighbouring cells swap contents (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.
- Standard name given: Schelling segregation model (two agent types + vacancies, tolerance 0.5, random relocation) with Kawasaki-type swap noise / noisy Schelling model
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Type 0 = vacancy, so this is literally Schelling's two-group lattice model at T=0.5 plus random neighbour swaps (temperature-like Kawasaki noise). Clustering and coarsening (g_moran, g_n_domains) are canonical. Noise-assisted extra clustering compared with rule 2 alone matches the documented role of noise in escaping large-but-finite cluster states. g_reversals is a near-trivial side effect of swap noise and unhappy-agent churn (each rule alone already gives about 0.6). Overturn is mild: from 'variant' to 'known', because no new phenomenology appears. Caveat: Moran's I on nominal labels 0/1/2 is a questionable metric, but the effect size makes real segregation clear.
- Citations:
  - https://arxiv.org/abs/0903.4694 (Gauvin, Vannimenus, Nadal - Phase diagram of a Schelling segregation model)
  - https://www.nature.com/articles/s41598-022-23224-7 (Aging effects in Schelling segregation model; noisy Schelling phases)
  - https://arxiv.org/pdf/1504.03809 (Barmpalias et al.; swap dynamic = Kawasaki)
  - Dall'Asta, Castellano, Marsili 2008, Statistical physics of the Schelling model of segregation
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:028: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_spec_peak
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell; (2) each cell copies the type of a random neighbour (every step).
- Standard name given: Schelling-Voter model (voter-model interfacial-noise coarsening with Schelling relocation)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The knock-out is a frozen grid (both rules read neighbours), so the large differences for reversals and spectral peak are inflated. Copying alone accounts for most of reversals (0.65 of 0.89) and the spectral peak (0.27 of 0.31). The spectral peak is probably a synchronous-update flicker artifact rather than a real oscillation. The domain count reflects unfinished, logarithmically slow voter coarsening on a 64x64 grid. Differences from the published Schelling-Voter model are in the rules only (unconditional synchronous copying, empty cells copied), not in the behaviour.
- Citations:
  - https://www.sciencedirect.com/science/article/abs/pii/S0378437117306118
  - https://www.sciencedirect.com/science/article/abs/pii/S0960077913001720
  - https://arxiv.org/html/1506.05321
  - https://arxiv.org/pdf/1506.05321
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:034: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_n_domains, g_corrlen, g_activity_trend, g_moran
- Rule: C kC=4 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0  (compact grammar notation)
- Standard name given: Multi-state (q=4) voter model / neutral copy-a-random-neighbour dynamics: voter-class coarsening toward consensus
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The rule is the voter model with 4 states. The NEXT_GE(th=1) gate only removes copies that would change nothing, so it at most rescales time. Observed features: fewer domains, longer correlation length, falling activity (the feature label says rising but the value is negative) and positive Moran's I. These are the textbook signs of voter coarsening in 2D (interface density falls as about 1/ln t, length scale grows with z=2). The pattern is real (Moran ~0 and corrlen 1 with interaction knocked out) but transient, heading for an absorbing consensus state. First pass called it a VARIANT; nothing about the behaviour is new, so the verdict is upgraded to KNOWN.
- Citations:
  - Starnini, Baronchelli, Pastor-Satorras, Ordering dynamics of the multi-state voter model, J. Stat. Mech. 2012 (arXiv:1207.5810)
  - Dornic, Chate, Chave, Hinrichsen, Critical coarsening without surface tension: the universality class of the voter model, PRL 87, 045701 (2001)
  - Tartaglia, Cugliandolo, Picco, Percolation and coarsening in the bidimensional voter model, arXiv:1506.05321
  - Dall'Asta & Galla, Algebraic coarsening in voter models with intermediate states, arXiv:0806.0817
  - Ben-Naim, Frachebourg, Krapivsky, Coarsening and persistence in the voter model, PRE 1996
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:037: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_type_entropy
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.1 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour.
- Standard name given: Spatial prisoner's dilemma (Nowak-May imitate-the-best) with conformity / majority-vote hybrid updating
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Both rules read neighbours, so the knock-out is just the frozen random initial state and is a weak control; the meaningful evidence is that the full program differs from each single-rule sub-program (majority-only coarsens, corrlen 10.5; game-only nearly absorbs, entropy 0.06). The combination gives the documented conformity-plus-reciprocity regime: biased but persistent coexistence (~27/73), small domains, boundary cells that keep flipping. Effects are modest (corrlen 1->2, integer-valued). Overturns first pass only in the sense of upgrading VARIANT to KNOWN: no phenomenological twist found.
- Citations:
  - Cui P.-B. & Wu Z.-X., Impact of conformity on the evolution of cooperation in the prisoner's dilemma game, Physica A 392, 1500-1509 (2013)
  - Szolnoki A. & Perc M., Conformity enhances network reciprocity in evolutionary social dilemmas, J. R. Soc. Interface 12, 20141299 (2015)
  - Szolnoki A. & Perc M., Leaders should not be conformists in evolutionary social dilemmas, Sci. Rep. (2016), arXiv:1604.06438
  - Nowak M.A. & May R.M., Evolutionary games and spatial chaos, Nature 359 (1992) [from background knowledge, not retrieved]
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:041: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran
- Rule: C kC=2 | C[any] CNT_GE(b=0,th=3) -> IMITATE_BEST(g=(1.5, 0.5)) p=1.0  (compact grammar notation)
- Standard name given: Nowak-May spatial Prisoner's Dilemma with deterministic imitate-the-best (unconditional imitation), b<1.8 regime where cooperators take over with leftover defector patches and blinkers
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The count gate (CNT_GE th=3) only restricts which cells update; no new phenomenology. Effects are robust vs knock-out, but the correlation length only rises from 1 to 2 (short range), and binary reversals likely reflect known period-2 blinkers. Overturns VARIANT to KNOWN because no phenomenological twist was identified. Not verified: that the (1.5, 0.5) payoff tuple maps exactly onto Nowak-May (b,R,P,S), and the 1992 paper was not opened directly.
- Citations:
  - Nowak & May 1992, Nature 359:826, Evolutionary games and spatial chaos
  - Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344, Sec 6.5
  - Roca et al., arXiv 0901.0355
  - arXiv 2607.24638 spectral cartography of spatiotemporal chaos in 2x2 games
  - github giakoumoglou/game_theory_spatial_PD (b-regime summary)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:042: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_n_domains, g_corrlen, g_activity_trend, g_reversals
- Rule: C kC=2 | C[any] CNT_LE(b=0,th=3) -> IMITATE_BEST(g=(1.9, 0.0)) p=0.1  (compact grammar notation)
- Standard name given: Nowak-May spatial (weak) Prisoner's Dilemma with imitate-the-best, semi-synchronous/probabilistic updating; C/D cluster coexistence with fluctuating boundaries (spatial chaos regime, b=1.9)
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Read the rule as imitate-the-best weak PD with b=1.9 and p=0.1 partial updating. CNT_LE gate (update only if <=3 neighbours in state 0) changes the rule but not the kind of behaviour. Features (domains merging, short correlation length of 2, persistent flipping) match documented C/D cluster dynamics. Activity trend is negative, i.e. settling or slow drift (consistent with the asynchronous-update literature), so the 'activity rising' label is wrong. Correlation length 2 means short-range clustering only. Overturn is mild: first pass said variant; no new phenomenology found. Did not confirm whether this exact gated version ends in lasting coexistence or slow defector takeover; either outcome is documented.
- Citations:
  - Nowak & May 1992, Nature 359:826 (doi:10.1038/359826a0)
  - Szabo & Fath 2007, Phys. Rep. 446:97 (arXiv cond-mat/0607344)
  - Huberman & Glance 1993, PNAS 90:7716
  - Prisoner's Dilemma with Semi-synchronous Updates, arXiv 0910.0961
  - Spectral cartography of spatiotemporal chaos in 2x2 games, arXiv 2607.24638
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:045: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_n_domains, g_moran, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) each cell copies the type of a random neighbour (every step); (2) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (every step).
- Standard name given: Noisy Griffeath cyclic cellular automaton '313' (R1/T3/C3/Moore) combined with voter-model copying, i.e. a three-state cyclic voter / spatial rock-paper-scissors (cyclic Lotka-Volterra) lattice model
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 2 alone is exactly Griffeath's named 313 CCA, documented as self-organising into spirals and domains. Rule 1 (synchronous random-neighbour copying) adds voter-model noise and coarsening, giving a standard three-state cyclic voter / RPS-type system. The measured features (clustering, fewer and larger domains, sustained cycling) are generic for this family. 'Reversals' in a 3-cycle largely counts cyclic advancement. The knock-out is a frozen random grid because both rules read neighbours, so the sub-program comparison is more informative. The first pass called it a VARIANT; the twist is only in the rule combination and the magnitude of the measures, not in the phenomenology, so KNOWN.
- Citations:
  - http://psoup.math.wisc.edu/mcell/rullex_cycl.html
  - https://github.com/setanarut/ccae
  - https://arxiv.org/pdf/nlin/0507067
  - Szabo & Szolnoki, Three-state cyclic voter model extended with Potts energy, PRE 65, 036115 (2002)
  - https://arxiv.org/html/1408.6828v1
  - https://doi.org/10.3390/g7030024
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:047: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_largest_domain, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.
- Standard name given: Mobility/long-range-mixing-induced global oscillation and loss of coexistence in cyclic-dominance (rock-paper-scissors / cyclic CA, Greenberg-Hastings-type) lattice dynamics
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: True; web searches: 4
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 1 alone is the 3-colour, threshold-1, 8-neighbour cyclic CA (locally periodic, quenched, fine-grained; zero reversals by construction). Rule 2 (Schelling relocation to random empty cells anywhere) acts as long-range mixing. Mixing is documented to synchronise cyclic dominance into global oscillations or near-fixation, which explains the ~81% largest single-type domain and the fewer domains. g_reversals is largely an artifact of measuring at fixed cells: vacating and refilling a cell reads as t->0->t, and rule 2 alone already gives 0.41. Exact rule pairing (empty state inside the cycle, Schelling-driven jumps) not found in the literature, but the behaviour is standard. Unresolved: whether the dominant type cycles over time (global oscillation) or stays fixed (near-absorbing); either way the behaviour is documented.
- Citations:
  - Reichenbach, Mobilia & Frey, Nature 448, 1046 (2007) https://www.nature.com/articles/nature06095
  - Biodiversity in models of cyclic dominance is preserved by heterogeneity in site-specific invasion rates (review of shortcut/mobility-induced global oscillations) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5137108/
  - Fisch, Gravner & Griffeath, Cyclic cellular automata in two dimensions (1991); phenomenology summarised in https://arxiv.org/html/nlin/0507067
  - Szabo, Szolnoki & Sznaider, PRE 76, 051921 (2007) https://arxiv.org/pdf/0710.3371
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:049: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Threshold cyclic cellular automaton / Greenberg-Hastings-type excitable CA (range 1, Moore, kappa=3, theta=3; turbulence/fixation domain formation), with Schelling-type random relocation acting as mixing that coarsens less
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Measured features (domain count, correlation length, Moran's I) are generic clustering signatures. Rule 1 alone (cyclic threshold CA, empty is part of the cycle) yields stronger clustering (Moran 0.77, corrlen 4); the Schelling relocation only dilutes it (full program sits between the two sub-programs on every feature), consistent with known mobility-induced breakdown of cyclic-competition patterns. No synergistic phenomenology. Mild overturn: VARIANT -> KNOWN, since the only difference is in the rules, not the behaviour. Pattern is genuine but weak (corrlen about 3 cells) and possibly near-fixated; whether it moves over time was not measured.
- Citations:
  - Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata, Statistics and Computing 1 (1991) 23-39
  - Fisch, Gravner & Griffeath, Cyclic Cellular Automata in Two Dimensions (1991), https://link.springer.com/chapter/10.1007/978-1-4612-0451-0_8
  - Shalizi et al., Quantifying Self-Organization in Cyclic Cellular Automata, https://arxiv.org/abs/nlin/0507067
  - Szczesny, Mobilia & Rucklidge, The influence of mobility rate on spiral waves in spatial rock-paper-scissors games, https://arxiv.org/abs/1608.08932
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:052: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: n_consensus_gain, n_type_entropy, n_edge_turnover
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 4 types (random start). Rules, applied in order each step: (1) each node adopts its neighbours' majority type (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (with probability 0.03 per step).
- Standard name given: Coevolving/adaptive majority-rule (adaptive voter-type) model with homophilic rewiring; Holme-Newman coevolution of networks and opinions; majority rule + rewiring on coevolving networks (Fernandez-Gracia et al.)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Consensus gain and entropy drop are essentially reproduced by majority rule alone (0.685/0.131 vs 0.6725/0.1175), so the coupling adds little to these measures. Edge turnover follows directly from the rewiring rule and drops in the full program because majority rule removes discordant links. Not a numerical artifact, but the features are standard majority-rule ordering plus rewiring that happens by construction. Incomplete consensus (~0.67) fits documented majority-rule traps on sparse networks. The measured features don't probe network fragmentation, which is the distinctive behaviour of the coupled system and is itself documented. The knock-out baseline (a frozen random colouring) is a weak control.
- Citations:
  - https://arxiv.org/pdf/1403.6106
  - https://arxiv.org/abs/physics/0603023
  - https://arxiv.org/html/2507.11521v2
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:054: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_activity_trend, g_corrlen, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.003 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour.
- Standard name given: Nowak-May spatial Prisoner's Dilemma (weak PD, imitate-the-best, Moore + self-interaction), low-b regime (b<1.75) with weak majority-vote noise: cooperator cluster growth/coarsening relaxing to static or period-2 (blinker) configurations
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Mild overturn: first pass said KNOWN-VARIANT; no phenomenological twist found, so KNOWN. At b=1.3 C clusters grow and D clusters shrink (documented thresholds b<1.8, b<2); deterministic imitate-best with self-interaction provably ends in fixed or periodic states, which explains falling activity (the 'activity rising' label has the wrong sign) and blinker reversals. Correlation length 2 is short, not long. Knock-out removes both rules (both read neighbours), leaving a frozen random field, so the contrast is against a trivial baseline. Rule 1 (p=0.003 majority vote) only adds slight smoothing and rare perturbations; each sub-program reproduces the same qualitative features. Genuine dynamics, but a relaxational, near-absorbing textbook outcome, not novel emergence. Sources verified via search snippets only; the Nature paper itself was not opened.
- Citations:
  - Nowak & May 1992, Nature 359:826 'Evolutionary games and spatial chaos' https://www.nature.com/articles/359826a0
  - Wu, Fu & Wang, 'Evolutionary games and spatial periodicity' arXiv:2209.08267
  - Nowak & May 1993, Int. J. Bifurc. Chaos 3:35 'The spatial dilemmas of evolution'
  - Szabo & Fath 2007, Physics Reports 446:97 'Evolutionary games on graphs'
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:056: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_type_entropy
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) neighbouring cells swap contents (with probability 0.05 per step).
- Standard name given: Spatial snowdrift (hawk-dove) game with Nowak-May best-takes-over imitation plus random mobility/diffusion
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Payoffs T=1.5>R=1>S=0.5>P=0 make this a snowdrift game, not a Prisoner's Dilemma. The observed outcome (one strategy at roughly 85% with the other persisting as a minority, short-range clusters) is the documented spatial snowdrift result. The swaps only modestly shift the strategy balance (largest domain 0.81 to 0.87), consistent with the mobility studies. g_reversals comes from the swaps alone and carries no information about the game. The knock-out baseline (an inert 50/50 mix) makes the significance tests easy to pass. The first pass's 'variant' label is too generous because no new behaviour was found; the variant looks closer to a reparameterisation of known work.
- Citations:
  - Hauert & Doebeli 2004, Nature 428:643 (doi:10.1038/nature02360)
  - Sicardi, Fort, Vainstein & Arenzon 2009, J. Theor. Biol. 256:240 (arXiv:0809.1479)
  - Vainstein, Silva & Arenzon 2007, J. Theor. Biol. 244:722
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:057: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a type-1 cell spontaneously becomes type 0 (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.
- Standard name given: Schelling segregation model (two-type, high vacancy density, random-destination relocation): diluted segregated phase with persistent churn
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 1 (type-1 cells become empty with p=0.1) is a one-way transient: type 1 is effectively gone after about 50 steps, so late-time behaviour is two-type Schelling with threshold 0.5 at about 67% vacancy, which is the documented diluted-segregation regime. Fewer domains than the knock-out = small like-type clusters; reversals = ongoing relocation churn, which is expected because moves go to random empty cells rather than satisfying ones. The knock-out removes all movement, so reversals=0 there by construction. 'Overturns first pass' only in the sense of upgrading VARIANT to KNOWN; no phenomenological twist found. Caveat: at high vacancy many agents have no occupied neighbours (0/0 share); how the code handles that probably sets much of the reversal rate, a partial implementation artifact.
- Citations:
  - Gauvin, Vannimenus & Nadal, Phase diagram of a Schelling segregation model, Eur. Phys. J. B 70, 293-304 (2009), arXiv:0903.4694
  - A Schelling model with switching agents: decreasing segregation via random allocation and social mobility, arXiv:1212.4929
  - Aging effects in Schelling segregation model, Sci. Rep. (2022), s41598-022-23224-7
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:060: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_corrlen, g_n_domains, g_activity_trend, g_moran
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 4 types (random start). Agents have 3 types (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (every step).
- Standard name given: Cyclic cellular automaton (Griffeath CCA), rule 1/3/4/Moore: metastable 'turbulent' regime with domain coarsening / possible spiral entrainment
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Only the cell rule is active; it is exactly the 4-colour, range-1 Moore, threshold-3 CCA, the case studied explicitly by Shalizi et al. citing Fisch-Gravner-Griffeath (T=3 = metastable turbulence that may become spirals; T>=4 fixation). The 400 agents are inert (no coupling) and do not affect classification. Measured clustering, ~10-cell correlation length, fewer domains and falling activity match documented self-organisation and relaxation on a finite lattice. Activity trend is negative despite the 'rising' label. Not a variant: no new phenomenology.
- Citations:
  - Fisch, Gravner & Griffeath, 'Cyclic cellular automata in two dimensions' (1991), https://math.ucdavis.edu/~gravner/papers/cca.pdf
  - Shalizi, Haslinger, Rouquier, Klinkner & Moore, 'Quantifying Self-Organization in Cyclic Cellular Automata', arXiv:nlin/0507067
  - MCell rules lexicon, Cyclic CA, http://psoup.math.wisc.edu/mcell/rullex_cycl.html
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:062: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_n_domains, g_moran, g_cyclic_changes
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.3 per step); (2) neighbouring cells swap contents (every step).
- Standard name given: Stochastic threshold cyclic cellular automaton (Fisch-Gravner-Griffeath CCA / Greenberg-Hastings family) with neighbour-exchange mobility, as in spatial rock-paper-scissors with mobility (Reichenbach-Mobilia-Frey)
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Changes the first pass from VARIANT to KNOWN. Rule 1 alone already gives cyclic wave domains (documented CCA self-organisation). Swapping makes the domains larger (higher corrlen) and noisier (slightly lower Moran and cycling, more domains), which matches the documented finding that mobility lengthens spiral wavelength below the critical mobility. Not an artifact: no freezing, not random colouring, holds up at larger size and longer runs. Unchecked: whether these are true spirals or drifting bands, and where threshold 3 / 3 colours / 8 neighbours sits in the CCA phase diagram; neither would make the behaviour new.
- Citations:
  - Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39 (1991)
  - Fisch, Gravner & Griffeath, Cyclic cellular automata in two dimensions, Spatial Stochastic Processes (1991)
  - Reichenbach, Mobilia & Frey, Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, Nature 448:1046 (2007)
  - Szolnoki et al., Cyclic dominance in evolutionary games: a review, J R Soc Interface 11:20140735 (2014)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:064: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_moran, g_activity
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 2 neighbours already are t+1 (every step); (2) neighbouring cells swap contents (every step).
- Standard name given: Cyclic cellular automaton (Griffeath CCA, 1/2/4/M; spiral/travelling-wave regime) with stirring/mobility, cf. mobility-enlarged spiral waves in cyclic-dominance lattice models
- Trivial or artifact: False; confidence: medium; overturns earlier verdict: False; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN-VARIANT
- Notes: Rule 1 alone is the documented CCA 1/2/4/M (persistent locally periodic waves; matches activity ~0.99). Adding swaps raises correlation length and Moran's I and halves activity, consistent with the documented growth of spiral/wave size with mobility. Only the exact pairing (deterministic synchronous CCA plus swaps every step) is not in the literature; the behaviour is. Patterns are sustained, clustered and stable across sizes, so not an artifact. Caveat: correlation length ~4.5 on 64x64 is modest; not checked by eye whether these are clean spirals or wave fragments. Sources checked via search excerpts, not full papers.
- Citations:
  - Fisch, Gravner, Griffeath (1991) Cyclic cellular automata in two dimensions / Threshold-range scaling of excitable cellular automata, Stat. Comput. 1:23-39
  - Shalizi et al., Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067 (range-1 Moore: T=2 spiral waves)
  - Reichenbach, Mobilia, Frey (2007) Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, Nature 448:1046
  - Reichenbach, Mobilia, Frey (2008) Self-organization of mobile populations in cyclic competition, J. Theor. Biol. 254:368
  - Roman, Konrad, Pleimling (2012) Cyclic competition of four species: domains and interfaces, J. Stat. Mech. P07014
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:066: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_activity_trend, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 7 neighbours already are t+1 (with probability 0.03 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.
- Standard name given: Spatial snowdrift (hawk-dove) game with Nowak-May synchronous imitate-the-best dynamics on a Moore lattice with self-interaction, plus weak threshold flip noise
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. Rule 2 alone already gives near-dominance (0.82), clustering and flip-back oscillations. Rule 1 (with 2 types it just means switch to the other type) is a 3%/step noise that removes isolated minority cells and keeps activity going; this explains the small extra gain in dominance and the rising activity. That activity signal halves at 3x run length, so it looks transient. Correlation length 1 -> 2 sites is short-range. 'Overturn' only means moving from KNOWN-VARIANT to KNOWN: no twist in the behaviour was found that would make it a variant.
- Citations:
  - Nowak & May 1992, Evolutionary games and spatial chaos, Nature 359:826
  - Hauert & Doebeli 2004, Spatial structure often inhibits the evolution of cooperation in the snowdrift game, Nature 428:643
  - Szabo & Fath 2007, Evolutionary games on graphs, arXiv cond-mat/0607344 (sec. 6.5)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:067: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_n_types, g_reversals, g_activity_trend
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) a type-1 cell spontaneously becomes type 0 (with probability 0.01 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour.
- Standard name given: Nowak-May spatial Prisoner's Dilemma (weak PD, imitate-the-best, Moore neighbourhood + self-interaction), low-b regime: cooperator-dominated coexistence with static defector clusters and blinkers; here with one-way D->C mutation noise
- Trivial or artifact: False; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Agrees with first pass on substance; relabels VARIANT->KNOWN because the only novelty is the D->C noise rule, which per the sub-program data only changes amounts (fewer domains: 3.18 vs 4.20) and leaves the corrlen/types/reversals phenomenology unchanged. Pattern is genuine but low-complexity (b=1.3 is well below the chaotic window ~1.8-2 for this neighbourhood).
- Citations:
  - Nowak & May 1992, Nature 359:826 (https://www.nature.com/articles/359826a0)
  - Szabo & Fath 2007, Evolutionary games on graphs, Phys. Rep. (arXiv cond-mat/0607344)
  - Stanford Encyclopedia of Philosophy, Evolutionary Game Theory (blinkers / two-cycle coexistence)
  - Nowak, Bonhoeffer & May 1994, More Spatial Games, IJBC 4:33 (errors/noise variants; not directly verified here)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:069: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_types, g_type_entropy, g_interface, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (with probability 0.1 per step); (2) a type-0 cell spontaneously becomes type 1 (with probability 0.1 per step).
- Standard name given: 3-colour cyclic cellular automaton / Bramson-Griffeath cyclic particle system (threshold 1, Moore neighbourhood), i.e. spatial rock-paper-scissors with one asymmetrically enhanced invasion rate
- Trivial or artifact: False; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Upgraded from KNOWN-VARIANT to KNOWN: rule 1 is the canonical stochastic 3-colour threshold-1 cyclic CA, proven to fluctuate (perpetual cycling, coexistence of all three colours) in 2D. Rule 2 is equivalent to raising the 1-invades-0 rate, an asymmetric RPS that gives uneven but persistent coexistence, which is documented. The knock-out baseline is a trivial frozen state (all 0s turn into 1s, then nothing changes), so the large feature differences partly reflect a degenerate baseline.
- Citations:
  - Fisch, Gravner & Griffeath (1991), Cyclic Cellular Automata in Two Dimensions, Spatial Stochastic Processes, Progress in Probability 19
  - Bramson & Griffeath (1989), Flux and fixation in cyclic particle systems, Ann. Probab. 17
  - arXiv:2108.06404 (CCA/GHM review: fluctuation on Z^2 for theta=1, all kappa>=3)
  - Frean & Abraham (2001), Rock-scissors-paper and the survival of the weakest, Proc. R. Soc. B 268:1323
  - arXiv:1907.11280, Predominance of the weakest species in LV and May-Leonard RPS
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:agents-2: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: a_step
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. Rules, applied in order each step: (1) each agent turns away from the centre of any type agents within distance 3.0.
- Standard name given: Repulsion-only (zone-of-repulsion) self-propelled particle model, i.e. a Couzin/Vicsek-type model with only the repulsion zone; disordered homogeneous active gas with reduced persistence
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The rule is the Couzin zone-of-repulsion rule on its own. Pure repulsion is documented to give a disordered homogeneous state. Flocking from turning away requires positional repulsion, which is absent here. a_step is not precisely defined; it is probably a persistence or displacement measure. Its sharp drop compared with the knock-out matches neighbour-induced turning wiping out persistence at about 28 neighbours per interaction disc. That is a real but trivial direct effect of the interaction, not an emergent pattern. 'Trivial' here means not an emergent pattern; it is not a numerical artifact.
- Citations:
  - https://arxiv.org/pdf/1301.6294
  - https://www.researchgate.net/publication/255989028_Self-propelled_particles_with_selective_attraction-repulsion_interaction_From_microscopic_dynamics_to_coarse-grained_theories
  - https://journals.aps.org/prx/abstract/10.1103/PhysRevX.14.031008
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:agents-8: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: a_step
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0.
- Standard name given: Chemotactic gradient ascent / trapping of active (persistent random-walk) particles at local maxima of a static, passively diffusing chemoattractant landscape
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: No emission and no agent-agent sensing, so the field is inert and the agents are independent. The knock-out removes the only rule, leaving persistent random walkers (persistence ~20 steps, which fits net displacement ~19.75). With gradient-following, agents localize at residual maxima of the diffusing initial field (~0.79). The suppressed displacement is a trivial effect of following a gradient, not an emergent collective pattern. It probably persists late in the run because the steering responds to gradient direction regardless of how weak the gradient is, giving a near-frozen trapped state. Calling it a 'variant' of a collective chemotaxis model overstates it.
- Citations:
  - https://arxiv.org/html/2511.18745
  - https://arxiv.org/pdf/1202.2442
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:grid-19: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_activity_cv, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.003 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour.
- Standard name given: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore neighbourhood plus self-play) with conformity-driven (majority) noise; near-absorbing single-strategy domination with perturbation-triggered avalanches
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: False; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Downgrades the first pass's VARIANT to KNOWN. Rule 2 alone gives 0.99 domination (near-absorbing frozen state). Rule 1 adds rare local perturbations that trigger short imitation cascades, which shows up as bursty activity whose burstiness score roughly halves at 4x area, consistent with independent local events rather than criticality. The correlation length mostly reflects the giant domain. The knock-out (static random colouring) is a weak baseline. The exact rule pair (global conformity at p=0.003 plus synchronous imitate-the-best with S=-0.5) was not found, but this is a parameter difference, not a new behaviour.
- Citations:
  - Nowak & May 1992, Nature, Evolutionary games and spatial chaos
  - Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344
  - Lim et al., Scale-invariant behavior in a spatial game of prisoners' dilemma, PRE 2002
  - Szolnoki & Perc 2015, Conformity enhances network reciprocity, J R Soc Interface rsif.2014.1299
  - Szolnoki & Chen 2018, Competition and partnership between conformity and payoff-based imitations, NJP, arXiv 1809.05079
  - Asynchronous spatial evolutionary games, BioSystems (frozen static configurations)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:grid-25: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_activity_cv
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) neighbouring cells swap contents (every step); (2) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step).
- Standard name given: Stochastic cyclic cellular automaton (Griffeath CCA, range-1 Moore, threshold 3, 4 colours, 'turbulent' regime) with stirring/mobility (cyclic dominance with exchange; Durrett-Neuhauser fast stirring; Reichenbach-Mobilia-Frey mobility vs pattern formation)
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Upgrades the first pass from VARIANT to plain KNOWN: no phenomenological twist found. The knock-out also removes swapping, so it is a frozen lattice and the reversal/activity-CV contrasts are inflated. Swap-only already gives reversals 0.40 and activity CV 0.52 (higher than the full system's 0.44), so these features are mostly stirring noise. The only real interaction effect is that stirring partly suppresses CCA domain growth (domain count 6.5, between 4.4 and 7.7), which is the documented mobility-vs-pattern trade-off. None of the features directly tests for spirals or waves.
- Citations:
  - Fisch, Gravner, Griffeath (1991) Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39
  - Shalizi, Shalizi, Haslinger (2004) Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067 (T=3, 4-colour, Moore: metastable turbulence)
  - Durrett (1999) Stochastic spatial models, SIAM Review (rapid stirring limits)
  - Reichenbach, Mobilia, Frey (2007) Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, Nature 448:1046
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:grid-6: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: g_largest_domain, g_n_domains
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) a type-0 cell becomes type 1 if at least 1 neighbours are type 1 (with probability 0.1 per step).
- Standard name given: SI epidemic probabilistic cellular automaton / Eden growth (irreversible contact spreading to all-type-1 absorbing state)
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Single irreversible rule 0->1 with probability 0.1 if any of 8 neighbours is type 1, from a ~50% random start: fills the lattice within tens of steps. g_largest_domain=1.0 is the all-type-1 absorbing state. g_n_domains=0.6931=ln2 probably means log(1+1), i.e. one domain, so the 'many domains' label looks wrong (inference about metric scaling). The knock-out removes the only rule, leaving a frozen random field, so the comparison is uninformative. No phenomenological twist found; first-pass 'variant' downgraded to known and trivial.
- Citations:
  - https://link.springer.com/chapter/10.1007/978-3-319-65558-1_12
  - https://www.researchgate.net/publication/323313442_Epidemic_Automaton_and_the_Eden_Model_Various_Aspects_of_Robustness
  - https://arxiv.org/pdf/2010.09677v1
  - https://en.wikipedia.org/wiki/Epidemic_models_on_lattices
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:g1:network-phase-2: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: p_r_std, p_freq_spread
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) each node is a phase oscillator pulled toward its neighbours' phases with strength 1.0; (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (every step).
- Standard name given: Kuramoto model on a homophily-fragmented (modular/disconnected) network: intra-community frequency locking with beating/oscillating global order parameter (cluster/modular synchronization)
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Static types + always-on cross-type rewiring drive the graph to a frozen absorbing topology of ~3 disconnected same-type components; rule 2 then becomes inert. Kuramoto locks each component at its own mean frequency (low effective-frequency spread), and the global order parameter beats between the incommensurate cluster phasors (high r std). No ongoing feedback between the rules; growth of the effect with run length points to a frozen-state/measurement artifact rather than a novel pattern. Overturns 'variant' toward plain KNOWN with trivial mechanism.
- Citations:
  - Community Structure and Multi-Modal Oscillations in Complex Networks, PLoS ONE 2013 (10.1371/journal.pone.0075569)
  - Rodrigues et al., The Kuramoto model in complex networks, arXiv:1511.07139
  - Cluster Synchronization of Kuramoto Oscillators and the Method of Averaging, arXiv:2204.02627
  - Papadopoulos et al., Chaos 27, 073115 (2017), arXiv:1706.08070 (adaptive Kuramoto rewiring -> modular organization, local not global synchrony)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:agents-16: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: a_density_cv
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 1 type (random start). Agents have 4 types (random start). Agents move at speed 1.0 with heading noise 0.3 rad per step. 2 field(s) with diffusion [0.02, 0.02]. Rules, applied in order each step: (1) each agent turns up the gradient of field 1; (2) each agent turns left on an even-type cell and right on an odd-type cell.
- Standard name given: Chiral active particles (circle swimmers) doing gradient ascent / trapping in a quenched random field landscape; chirality reduces accumulation
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Cells have only one type, so rule 2 is a constant left turn (chirality); alone it gives the same CV as the knock-out (0.799 vs 0.8). No emission rule is listed, so field 1 is probably decaying initial noise that acts as a fixed landscape; agents pile up at its peaks (CV 6.49 with rule 1 alone). The full program's CV of 1.84 is chirality partly breaking that trapping, a documented effect. The 'interaction' is agent-reads-environment, not agent-reads-agent. The clustering follows the field's initial noise rather than self-organisation, so it is close to a trivial outcome. Overturns 'variant' toward KNOWN/trivial. If agents actually emit field 1, reclassify as chiral Keller-Segel aggregation (still known).
- Citations:
  - https://www.researchgate.net/publication/254498325_Active_Brownian_particles_in_two-dimensional_traps
  - https://www.researchgate.net/publication/331170019_Active_chiral_particles_under_confinement_Surface_currents_and_bulk_accumulation_phenomena
  - https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/unsteady-chiral-swimmer-and-its-response-to-a-chemical-gradient/39FD23E37ED399BFD4324811DD5A1785
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:field-1: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: f_std, f_n_blobs, f_peak_sharp
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 1 field(s) with diffusion [0.02]. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 8.0; (2) each agent deposits 0.1 of field 0 where it stands.
- Standard name given: Vicsek model flocking / clustering and band formation, read out through a passive (non-feedback) deposited tracer field
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Nothing in the program reads field 0, so it has no effect on the agents and only records a diffusion-blurred, time-integrated map of agent density. Higher f_std, fewer blobs and sharper peaks all come from known Vicsek clustering/banding (density about 0.1 with about 20 neighbours inside radius 8 puts the system in the ordered/clustered regime). The knock-out keeps deposition and removes alignment, so the measured differences just quantify Vicsek density inhomogeneity. The rule-1-only field features are meaningless because no deposition happens. The underlying flocking is genuine and robust, but the field features are an inert-tracer readout. Overturns the first pass from VARIANT to KNOWN, with no twist in the behaviour.
- Citations:
  - https://guava.physics.ucsd.edu/~nigel/Courses/Web%20page%20563/Essays_2017/PDF/Chatterjee.pdf
  - https://arxiv.org/pdf/1511.01451
  - https://arxiv.org/pdf/2604.28167
  - https://arxiv.org/html/2610.26569
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:field-3: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: f_std, f_n_blobs
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 2 field(s) with diffusion [0.02, 0.1]. Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.1 per step); (2) type-1 cells release 1.0 of field 0 per step.
- Standard name given: Majority-vote cellular automaton / zero-temperature Glauber (Ising) coarsening, read out by a passive, non-decaying emitted diffusive field
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Only interaction is majority-rule cell updating, i.e. known Ising-like coarsening into few large domains (lower f_n_blobs). Field 0 is emitted one-way, never read back, and never decays, so f_std grows without limit as it builds up over large domains (difference roughly triples at 3x run length): an inert field plus a build-up artifact. Field 1 and the agents are not used by any rule. End state tends to a frozen absorbing configuration. First pass's 'variant' label is overturned to KNOWN, with the features largely trivial.
- Citations:
  - https://arxiv.org/pdf/2604.10082
  - https://www.sciencedirect.com/science/article/abs/pii/S0304397521004370
  - https://link.springer.com/article/10.1007/s00220-018-3180-2
  - https://arxiv.org/pdf/cond-mat/9701118
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:grid-4: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_types, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.03 per step).
- Standard name given: Nowak-May spatial Prisoner's Dilemma (unconditional imitation, Moore + self) mixed with voter-model imitation; finite-size consensus/fixation to an absorbing state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The two features (one surviving type, largest domain = 1.0) just describe finite-size collapse to all-C or all-D consensus. Vilone et al. 2014 study the same combination (unconditional imitation in the PD plus voter imitation with probability q) and report that voter dynamics strengthens consensus, i.e. uniform asymptotic states. At b=1.3 with self-interaction, Nowak-May is outside its chaotic regime (9/5<b<2) and tends to freeze; voter noise then drives coarsening to consensus. The knock-out comparison is uninformative because the interaction rules are what produce consensus. This is a frozen absorbing state, not an emergent pattern; the first pass's 'variant' label is not justified.
- Citations:
  - Vilone, Ramasco, Sanchez, San Miguel, 'Social imitation versus strategic choice, or consensus versus cooperation, in the networked Prisoner's Dilemma', PRE 90, 022810 (2014), arXiv:1405.4997
  - Nowak & May, 'Evolutionary games and spatial chaos', Nature 359, 826 (1992)
  - Szabo & Fath, 'Evolutionary games on graphs', Phys. Rep. 446, 97 (2007), arXiv:cond-mat/0607344
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:grid-7: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_reversals, g_interface
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 2 field(s) with diffusion [0.02, 0.1]. Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.1 per step); (2) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Majority-vote model with noise (Oliveira 1992) / kinetic Ising coarsening interrupted by noise, disordered high-noise phase
- Trivial or artifact: True; confidence: high; overturns earlier verdict: False; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Agents flipping cells act as random-walking noise sources. They flip ~10% of cells per step while the majority rule updates ~10% per step, so effective noise is ~0.5, far above q_c≈0.075: a disordered phase. g_reversals is trivial (rules-2-only already gives 1.0). g_interface (0.43) sits just below random colouring (0.50), i.e. weak short-range clustering competing with noise. The chemical fields are inert. Agreeing with 'known' but sharpening it: this is not a meaningful variant, and the features are mostly trivial.
- Citations:
  - de Oliveira, M.J., Isotropic majority-vote model on a square lattice, J. Stat. Phys. 66, 273-281 (1992), doi:10.1007/BF01060069
  - Kwak et al., Critical behavior of the majority voter model is independent of transition rates, PRE 75, 061110 (2007)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## census:novelty:grid-8: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 2 field(s) with diffusion [0.02, 0.1]. Rules, applied in order each step: (1) each agent advances the type of the cell under it by one (cyclically); (2) each agent turns up the gradient of field 0.
- Standard name given: Random-walk bit-flipping agents (turmite without colour feedback) plus chemotactic gathering at maxima of a passive, emitter-free field
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: g_reversals equals the rule-1-only value (1.0) and comes straight from 2-state cyclic increment = flip on revisits. g_activity_trend (negative) comes from agents gathering at maxima of field 0, which nobody emits into: fewer distinct cells get flipped, and same-step double flips cancel. No feedback from cells to agents or from agents to fields, so nothing self-organises; field 1 does nothing. Overturns the first pass only in the less-novel direction (KNOWN and trivial, not a meaningful variant).
- Citations:
  - https://arxiv.org/pdf/1702.05547
  - https://en.wikipedia.org/wiki/Langton%27s_ant
  - https://www.nature.com/articles/s41598-021-83963-x
  - https://arxiv.org/pdf/1202.2442
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:004: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_activity_cv, g_corrlen, g_n_types, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell adopts its neighbourhood's majority type (with probability 0.1 per step).
- Standard name given: Spatial snowdrift/hawk-dove game with imitate-the-best (Nowak-May-type) plus majority-vote coarsening; fixation to a single-type absorbing state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. The lattice collapses to one type (absorbing state). n_types=1, largest_domain=1 and corrlen=0 just describe that uniform lattice, and the negative activity trend is the decay into it. The knock-out and sub-program comparisons are against 'None' (no dynamics), so they don't show emergence. The first pass's 'variant' label overstates novelty: there's no twist in the phenomenology.
- Citations:
  - Nowak & May 1992, Nature 359:826 (Evolutionary games and spatial chaos)
  - Hauert & Doebeli 2004, Nature 428:643 (Spatial structure often inhibits the evolution of cooperation in the snowdrift game)
  - Killingback & Doebeli 1996, Proc R Soc B 263:1135 (Spatial evolutionary game theory: Hawks and Doves revisited)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:005: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_activity_trend
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.1 per step).
- Standard name given: Spatial snowdrift (hawk-dove) game with unconditional imitation mixed with voter-model social imitation; voter-driven coarsening / finite-size consensus (Vilone et al. hybrid UI+VM dynamics)
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Payoffs T>R>S>P make this the snowdrift game (Hauert-Doebeli parametrisation, r=0.5), not PD. Mixing best-neighbour imitation with voter copying is the Vilone et al. hybrid; the voter part acts as interface noise and drives finite-size symmetry breaking toward one dominant type. Observed few domains, ~94% one type and falling activity match that picture. Both rules preserve absorbing states and there is no mutation, so this is slow drift toward fixation in a 64x64 world. The knock-out is a frozen random colouring, so the knock-out contrast is trivially satisfied. Mislabelled features: 'activity rising' is negative (falling), 'long corrlen' is 2 lattice units, 'reversals' is a binary flag. The snowdrift-vs-PD/coordination difference is a rule difference, not a new phenomenology; each rule alone already shows most of the ordering.
- Citations:
  - Hauert & Doebeli, Nature 428, 643 (2004) doi:10.1038/nature02360
  - Vilone, Ramasco, Sanchez, San Miguel, Phys. Rev. E 90, 022810 (2014) arXiv:1405.4997
  - Vilone, Ramasco, Sanchez, San Miguel, Sci. Rep. 2, 686 (2012) 'Social and strategic imitation: the way to consensus'
  - Roca, Cuesta, Sanchez, arXiv:0911.1720 (review, UI snowdrift on k=8 lattice)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:016: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_largest_domain, g_n_types, g_reversals, g_n_domains, g_spec_peak
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell adopts its neighbourhood's majority type (with probability 0.01 per step).
- Standard name given: Spatial snowdrift/hawk-dove game with imitate-the-best (Nowak-May) plus conformity/majority-rule updating; coarsening to a frozen all-C or all-D absorbing state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. The end state is a single-type, frozen absorbing monoculture: n_types=1, largest_domain=1, activity falling. The reversals and spectral-peak scores are probably short-lived effects of the coarsening phase. The knock-out baseline is a static random colouring, which makes the robustness differences look bigger than they mean. Overturned only in degree: from variant to plain known/trivial.
- Citations:
  - Nowak & May 1992, Nature 359:826 (Evolutionary games and spatial chaos)
  - Hauert & Doebeli 2004, Nature 428:643 (spatial snowdrift game)
  - Evolutionary game dynamics of combining the payoff-driven and conformity-driven update rules, Chaos Solitons Fractals 2020 (S0960077920305427) - always absorbed into full-C or full-D frozen state
  - Szolnoki & Perc, Conformity enhances network reciprocity in evolutionary social dilemmas, arXiv:1412.4113
  - Szolnoki & Perc 2018, Competition and partnership between conformity and payoff-based imitations, New J. Phys.
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:023: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_n_types, g_n_domains, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.002 per step).
- Standard name given: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore neighbourhood with self-interaction), donation-game payoffs, plus weak voter-model noise; outcome = fixation to a single-strategy absorbing state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 1 alone already fixes to one type (largest domain 0.996). Rule 2 (copy a random neighbour) cannot reintroduce an extinct type, so the end state is a frozen uniform absorbing state. g_reversals=1 most likely reflects flips during the transient, not a sustained dynamic. The knock-out comparison only contrasts fixation with slow voter-model coarsening. The first pass's 'variant' label overstates novelty: the different payoffs (S=-0.5, T=1.5) just place the system in the parameter region where one strategy takes over, a region already mapped in the follow-up literature.
- Citations:
  - Nowak & May 1992, Nature 359:826 (https://www.nature.com/articles/359826a0)
  - Roca, Cuesta & Sanchez 2009, Phys Rev E 80:046106 (T-S plane outcomes under unconditional imitation on lattices)
  - https://wisc.pb.unizin.org/agent-based-evolutionary-game-dynamics/chapter/iii-1/
  - https://arxiv.org/pdf/0901.0355
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:031: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_corrlen, g_activity_cv, g_n_domains, g_n_types, g_reversals, g_largest_domain
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) a cell of type t becomes type t+1 (cyclically) if at least 6 neighbours already are t+1 (with probability 0.03 per step).
- Standard name given: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore+self) with conformity/majority-rule noise: fixation/near-consensus of one strategy
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 3
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: With 2 types the 'cyclic' rule is just a stochastic 6-of-8 majority (conformity) rule. Rule 1 alone already yields largest domain 0.989 and correlation length 3, so the knock-out contrast overstates emergence. The coupling only removes residual minorities (1.5 -> 1 surviving types). High activity CV and reversals come from sporadic flips in a near-frozen, nearly uniform state, which is close to an absorbing state, so the outcome is trivial. First pass said VARIANT; it is plainly a known behaviour.
- Citations:
  - Nowak & May 1992, Nature 359:826 (doi:10.1038/359826a0)
  - Szabo & Fath 2007, Evolutionary games on graphs, arXiv:cond-mat/0607344
  - Grujic et al. 2010, PLoS ONE e13749 (imitate-the-best Moore lattice gives cooperation ~0 or ~1 depending on initial condition)
  - Szolnoki & Perc 2015, J. R. Soc. Interface 12:20141299 (conformity enhances network reciprocity)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:036: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: f_std, f_n_blobs
- Rule: World: 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.02]. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 3.0; (2) each agent deposits 0.1 of field 0 where it stands.
- Standard name given: Vicsek model flocking: density clustering / traveling bands / giant number fluctuations, recorded in a passive (inert, never-read) accumulated deposit field
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: No agent reads field 0 and it does not decay, so it is a time-integrated record of agent density. Higher f_std and fewer blobs only show that Vicsek alignment makes agents cluster into a few coherent flocks, which is standard Vicsek phenomenology. The flocking itself is real; the field features are a readout of an inert coupled field, so they are trivial. Because the field never decays, f_std also depends partly on run length. The first pass called this a VARIANT; it is plain KNOWN.
- Citations:
  - https://arxiv.org/html/2312.15525v2
  - https://guava.physics.ucsd.edu/~nigel/Courses/Web%20page%20563/Essays_2017/PDF/Chatterjee.pdf
  - https://arxiv.org/html/2604.28167v2
  - https://arxiv.org/html/1709.01499
  - https://www.researchgate.net/publication/283531761_The_Physics_of_the_Vicsek_model
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:041: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_activity_cv, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Gradient-ascent (chemotactic) aggregation of agents in a passive, decaying, non-emitted field, plus parity toggling of a 2-state substrate (Langton's-ant-like marking)
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Field 0 has no source, so it only diffuses from its random start; rule 1 is passive hill-climbing to the initial maxima. Falling activity = agents crowding onto the same cells, with toggles cancelling under the 2-state parity. High CV = small-mean relaxation fluctuations (effect halves at 4x area). g_reversals=1 is guaranteed for 2-state cyclic toggling (rule-2-only also gives 1.0). The knock-out deletes all cell updates, so its zero baseline is uninformative; the fair baseline is rule-2-only. 'Overturns' only because the first pass's VARIANT label overstates novelty; the family identification stands. Reasoning based on the rule description; not rerun.
- Citations:
  - Chatterjee & Goldenfeld, Field-theoretic model for chemotaxis in run and tumble particles, PRE 2021 (aggregation in externally imposed attractant field)
  - Schweitzer, Brownian Agents and Active Particles (2007) / Helbing et al. PRE 56, 2527 (1997) - active walker framework, here lacking the feedback loop
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:046: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: f_std, f_n_blobs, f_peak_sharp
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) each cell copies the type of a random neighbour (with probability 0.2 per step); (2) type-0 cells release 1.0 of field 1 per step.
- Standard name given: 2D (multi-state) voter model coarsening with a passive diffusing tracer of one opinion
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 1 is the neighbour-copying voter model (4 states, Moore neighbourhood), which coarsens slowly (active interfaces ~1/ln t) towards absorbing consensus. The field never feeds back. Field 1 is a non-decaying integral of the type-0 indicator, so f_std grows without bound with run length (62.8 -> 141.4 at 3x run length). Fewer blobs and sharper peaks just reflect voter domain coarsening. Field 0 has no source and is inert. The domains are genuine voter coarsening, but the field features are a passive accumulating readout heading for a consensus absorbing state, not a new pattern.
- Citations:
  - Dornic, Chate, Chave, Hinrichsen, PRL 87, 045701 (2001) - Critical coarsening without surface tension: the universality class of the voter model
  - Tartaglia et al., Percolation and coarsening in the bidimensional voter model, arXiv:1506.05321
  - Roy & Sen, Interplay of interfacial noise and curvature driven dynamics in two dimensions, arXiv:1701.09117
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:052: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_largest_domain, g_n_types
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). 2 field(s) with diffusion [0.02, 0.01]. Rules, applied in order each step: (1) field 1 grows autocatalytically by consuming field 0 (rate 1.0 x f x g^2); (2) a type-0 cell becomes type 1 where field 1 exceeds 0.8.
- Standard name given: Cubic autocatalysis (A+2B->3B) reaction-diffusion front driving an irreversible threshold conversion to an absorbing single-type state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Field rule is closed cubic autocatalysis (Gray-Scott without feed/kill): autocatalyst invades and saturates everywhere via documented travelling fronts. Cell rule is a one-way type-0 to type-1 flip above threshold, with no reverse process and no feedback onto fields, so the lattice converges to an all-type-1 absorbing state. largest_domain=1 and n_types=1 are the trivial endpoint; knock-out differences only show that the threshold is reached when the field can grow. Not a variant: no new phenomenology. The first pass overrated novelty.
- Citations:
  - https://royalsocietypublishing.org/rspa/article/430/1880/509/16458
  - https://arxiv.org/pdf/chao-dyn/9601004
  - https://link.springer.com/article/10.1007/s10910-026-01776-7
  - https://royalsocietypublishing.org/doi/10.1098/rsta.1991.0098
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:053: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: g_n_domains, g_n_types, g_activity_trend, g_type_entropy
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.001 per step); (2) a type-0 cell becomes type 1 where field 0 exceeds 0.5.
- Standard name given: Zero-temperature q-state Potts / majority-vote coarsening (q=4) from a biased initial composition; type 0 removed once by a passive diffusing field that does not respond to the cells
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 1 alone gives the coarsening, the falling activity and the domain counts. Rule 2 alone gives 3 surviving types and entropy of about 0.75 (simple arithmetic: type 0 gone, type 1 about 0.5, types 2 and 3 about 0.25 each). The field has no source and no feedback, so it acts only as a one-time bias on the starting mix. The full program's lower entropy and fewer domains are just coarsening enlarging the type that already has the most cells. The first pass's 'variant' label overstates this: there is no new twist in the behaviour. The knock-out counts the passive field as an interaction, which inflates the interaction effect. The coarsening is genuine but textbook; the drop in types and entropy is an artifact of the passive field.
- Citations:
  - http://physics.bu.edu/~redner/pubs/pdf/potts2d.pdf
  - https://arxiv.org/html/cond-mat/9503163
  - https://arxiv.org/html/cond-mat/9705155
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:056: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 3 types (random start). Agents have 3 types (random start). Agents move at speed 0.1 with heading noise 0.1 rad per step. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 1.0; (2) each cell adopts its neighbourhood's majority type (every step); (3) each agent advances the type of the cell under it by one (cyclically).
- Standard name given: Period-2 oscillation of the synchronous majority-rule cellular automaton (Goles–Olivos / Poljak–Sůra theorem), plus agents perturbing cells and the majority rule restoring them
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The majority rule alone gives a flip score of 0.85, higher than the full program's 0.58, so the feature comes from the known period-2 cycles of updating every cell at once. Agents bumping cell types (rule 3 alone: 0.36) add flips and break up the cycles; heading alignment (rule 1) has no effect. The knock-out removes the majority rule, so its zero score is guaranteed and the 0.58 gap says nothing new. The first pass's VARIANT verdict is downgraded to KNOWN because there is no twist in the behaviour itself, and the feature is largely an artifact of synchronous updating.
- Citations:
  - Goles & Olivos 1980, Periodic behaviour of generalized threshold functions, Discrete Math 30:187-189
  - Poljak & Sůra (period <=2 for majority dynamics)
  - arXiv:1711.10920 (Biased) Majority Rule Cellular Automata
  - arXiv:2610.04311 The Two-Dimensional Majority Rule is P-Complete
  - Frischknecht, Keller, Wattenhofer 2013, Convergence in (Social) Influence Networks
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:064: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: f_mean, f_change, f_n_blobs, f_std
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 3 types (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) a type-1 cell becomes type 0 if at least 4 neighbours are type 0 (every step); (2) type-0 cells release 1.0 of field 0 per step.
- Standard name given: Irreversible threshold growth (threshold-4 bootstrap percolation, 8-neighbour lattice) reaching the all-type-0 absorbing state, with a field that only accumulates (no decay) and grows linearly
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule 1 can only turn cells type 0, never back. From a random 50/50 start the 64x64 grid almost certainly fills (rule-1-only f_mean about 0.98). Every cell then emits 1.0 per step with no decay, so the field is uniform (f_std about 0) and grows linearly: f_mean about 1000, roughly double the knock-out's 500 because all cells emit instead of about half; the gap triples at 3x run length. f_n_blobs = ln 2 means a single grid-wide blob. The field never acts back on the cells and the agents take part in no rule. The features are bookkeeping of a frozen absorbing state, not an emergent pattern.
- Citations:
  - https://arxiv.org/pdf/2501.18976
  - https://link.springer.com/article/10.1007/s00493-026-00213-2
  - https://arxiv.org/pdf/2605.01151
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## qd:065: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: False; robust features: a_density_cv, a_local_align, a_nn_ratio
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 4 types (random start). Agents move at speed 2.0 with heading noise 0.0 rad per step. 2 field(s) with diffusion [0.02, 0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) field 1 grows autocatalytically by consuming field 0 (rate 1.0 x f x g^2).
- Standard name given: Chemotactic gradient-ascent aggregation of non-interacting self-propelled particles in a one-way-coupled (externally evolving) attractant field; field = unfed Gray-Scott (U+2V->3V) depletion transient
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Agents only read field 0; they neither emit nor consume, and do not sense each other, so the knock-out removes single-agent chemotaxis, not agent-agent interaction. Clustering and alignment come from agents climbing the same field and gathering at its maxima. Rule 2 is Gray-Scott without feed or decay, so field 0 is used up and its gradients keep moving and weakening; this explains why the full program clusters less than rule 1 alone. With zero noise, late-time steering may follow near-round-off gradients; check the size of field 0's gradient over time.
- Citations:
  - https://arxiv.org/html/2011.04580v1
  - https://itp.uni-frankfurt.de/~gros/StudentProjects/Projects_2020/projekt_schulz_kaefer/
  - https://www.mrob.com/pub/comp/xmorphia/index.html
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:011: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: n_components
- Rule: N kN=2 | N[0] ALWAYS() -> REWIRE_SAME() p=1.0  (compact grammar notation)
- Standard name given: Rewiring-only (phi=1) limit of the coevolving/adaptive voter model; homophilic network segregation / graph fission
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Two fixed node kinds plus always-rewire-to-same-kind is the opinion-frozen, rewiring-only limit of the Holme-Newman/adaptive voter model. Cross-kind links are removed and never re-created, so the graph must split into at least one component per kind (2) plus small leftover pieces, matching the measured ~2.9 vs 1.5. The split is real and robust, but it follows almost directly from the rule and can't be undone (a frozen segregated state even though links keep moving), so it is not a meaningful emergent pattern. First pass's VARIANT label is downgraded to KNOWN because no twist in the behaviour was found.
- Citations:
  - Holme & Newman, Phys. Rev. E 74, 056108 (2006)
  - Henry, Pralat & Zhang, PNAS 108, 8605 (2011)
  - Durrett et al., Graph fission in an evolving voter model, PNAS 109, 3682 (2012)
  - Coevolving nonlinear voter model phase diagram incl. absorbing fragmented phase
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:025: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_types
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Majority-rule / zero-temperature Potts coarsening to consensus (absorbing homogeneous state), with Schelling-type relocation acting as interface noise
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Only feature is g_n_types -> 1, i.e. fixation into a single-type absorbing state (all-empty or all-one-type, after which relocation is inert). Rule-1-only freezing at 2 types matches known majority/Potts freezing (q>=3); adding relocation noise removes the frozen interfaces and lets the system reach consensus, also documented. The 4x-area result of 2 types fits finite-size consensus timing. The knock-out comparison is trivially guaranteed. Frozen absorbing state, not an emergent pattern.
- Citations:
  - https://arxiv.org/pdf/1811.11888
  - https://arxiv.org/html/2604.10082
  - https://arxiv.org/pdf/2208.08770
  - https://arxiv.org/html/physics/0609079
  - https://arxiv.org/html/cond-mat/0407622
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:026: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_reversals, g_activity_cv
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) neighbouring cells swap contents (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.
- Standard name given: Schelling segregation model with vacancies and noise (random neighbour exchange ~ infinite-temperature Kawasaki noise); segregated/mixed competition
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Only g_n_domains reflects genuine (weak) Schelling segregation partly undone by swap noise; rule-2-only already gives 6.40 vs 5.89 full. g_activity_cv is entirely produced by the unconditional swap rule (rule-1-only 0.756 vs full 0.752), and g_reversals is mostly swap noise (rule-1-only 0.50) plus relocation/swap undoing. The knock-out scores of 0 for those two features are artifacts of the knock-out removing the swap rule. Not a new variant.
- Citations:
  - Gauvin, Vannimenus & Nadal, Phase diagram of a Schelling segregation model, EPJ B 70:293 (2009), arXiv:0903.4694
  - Dall'Asta, Castellano & Marsili, Statistical physics of the Schelling model of segregation, arXiv:0707.1681
  - Aging effects in Schelling segregation model, Sci. Rep. (2022), s41598-022-23224-7
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:036: DOCUMENTED (KNOWN)

- Class: SINGLE MECHANISM; flagged: True; robust features: n_oscillation
- Rule: N kN=4 | N[any] ALWAYS() -> COPY_RAND() p=0.2  (compact grammar notation)
- Standard name given: 2D voter model (nearest-neighbour, asynchronous/random-sequential); interfacial-noise-driven coarsening with 1/ln t decay of active interfaces
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Rule is exactly the voter model (copy random von Neumann neighbour with prob 0.2 per step = rescaled time), not a variant. 'n_oscillation' is aperiodic stochastic flipping at domain interfaces; its decline at 3x run length (0.72->0.37) matches logarithmic coarsening toward an absorbing consensus state, and insensitivity to 4x area matches domains smaller than system size. The knock-out removes the only dynamics, so the contrast is trivially guaranteed. Overturns first pass only in calling it the canonical model rather than a variant; no novel phenomenology. The 'oscillation' label is an artefact of the metric (no periodicity).
- Citations:
  - Dornic, Chaté, Chave, Hinrichsen, 'Critical Coarsening without Surface Tension: the Voter Universality Class', PRL 2001, arXiv:cond-mat/0101202
  - Tartaglia, Cugliandolo, Picco, 'Percolation and coarsening in the bidimensional voter model', arXiv:1506.05321
  - Ben-Naim, Frachebourg, Krapivsky, 'Coarsening and persistence in the voter model' (opinion-change statistics), PRE 1996
  - Interplay of interfacial noise and curvature driven dynamics in two dimensions, arXiv:1701.09117
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:053: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: n_edge_turnover, n_oscillation
- Rule: World: a random network of 200 nodes (average 4 links each). Nodes have 4 types (random start). Rules, applied in order each step: (1) a node of type t becomes t+1 (cyclically) if at least 1 neighbours are t+1 (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (every step).
- Standard name given: Cyclic-dominance (rock-paper-scissors / cyclic CA) dynamics on an adaptive network with Holme-Newman-type same-type rewiring (adaptive RPS / coevolving cyclic voter model)
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: The two measured features (edge turnover, state oscillation) are each produced by a single rule on its own (rule 2 alone gives turnover 0.74; rule 1 alone gives oscillation 0.82), so the knock-out and sub-program tests only show the two rules adding up, not interacting. Combining the rules lowers oscillation, and it falls further at larger size and longer runs, consistent with a slow drift towards the documented frozen, fragmented phase at high rewiring strength. The 4 types, synchronous threshold-1 update and every-step rewiring change the rules but add no new behaviour. The first pass's 'variant' verdict is downgraded to KNOWN because no twist in the behaviour was identified.
- Citations:
  - Demirel, Prizak, Reddy, Gross, 'Cyclic dominance in adaptive networks', EPJ B 84, 541 (2011), arXiv:1011.1124
  - 'Adaptive cyclically dominating game on co-evolving networks: numerical and analytic results', EPJ B (2017), arXiv:1605.03370
  - Holme & Newman, 'Nonequilibrium phase transition in the coevolution of networks and opinions', PRE 74, 056108 (2006)
  - Lutz, Risau-Gusman, Arenzon, 'Intransitivity and coexistence in four species cyclic games', J. Theor. Biol. 317, 286 (2013)
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:058: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_activity_trend, g_reversals
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a type-1 cell spontaneously becomes type 0 (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Schelling segregation model (two agent types + vacancies, random-vacancy relocation) with one-way removal of one population
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 1
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Type 0 = vacancy, so rule 2 is standard Schelling relocation at threshold 0.3 (Moore neighbourhood, random vacancy). Rule 1 removes type-1 agents at 10% per step, so they are essentially extinct within about 50 steps. After that every type-2 agent is fully satisfied and the lattice freezes into an absorbing state. g_reversals is just Schelling relocation seen at the cell level (empty/occupied toggling). It is nonzero only while decay keeps creating dissatisfied minorities. g_activity_trend is negative (activity falls, not rises), consistent with relaxation into the frozen single-type state. Passing the knock-out and the size/length checks reflects a size-independent per-cell mechanism, not genuine emergent structure. The first pass's 'variant' label is too generous: no new phenomenology was found.
- Citations:
  - https://www.nature.com/articles/s41598-022-23224-7
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC9653388
  - https://link.gale.com/apps/doc/A348647898/AONE
  - https://arxiv.org/pdf/0711.2212
  - https://arxiv.org/pdf/1002.3758
  - https://en.wikipedia.org/wiki/Schelling%27s_model_of_segregation
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:063: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_domains, g_type_entropy, g_spec_peak
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.005 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.7 moves to a random empty cell.
- Standard name given: Schelling-Voter model (Schelling segregation + majority/voter type-adoption), multi-type; coarsening toward type extinction/consensus
- Trivial or artifact: True; confidence: medium; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Downgraded from KNOWN-VARIANT: no phenomenological twist found. Rule-2-only (threshold 0.7, 3 occupied types) stays in the documented disordered/no-coarsening Schelling regime; the majority rule adds coarsening toward one type taking over (documented extinction states). Entropy and domain-count effects grow strongly with run length, so this is a transient approach to an absorbing state, not a stationary pattern. g_spec_peak is very likely low-frequency trend power (knock-out = 0 because it is frozen), not a real oscillation. Domain/entropy changes are genuine interaction effects; the spectral feature is likely an artifact.
- Citations:
  - https://www.sciencedirect.com/science/article/abs/pii/S0378437117306118
  - https://www.sciencedirect.com/science/article/abs/pii/S0960077913001720
  - https://arxiv.org/pdf/1905.01953
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

## review:068: DOCUMENTED (KNOWN)

- Class: COMBINATION EFFECT; flagged: True; robust features: g_n_types
- Rule: World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a type-0 cell becomes type 1 if at least 3 neighbours are type 1 (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.
- Standard name given: Vacancy exhaustion: polluted/threshold bootstrap percolation (Moore, r=3) filling the vacancies of a Schelling model, ending in a frozen jammed state
- Trivial or artifact: True; confidence: high; overturns earlier verdict: True; web searches: 2
- Twist (if any): 
- First run (cards also showing the earlier verdict): KNOWN
- Notes: Type-2 count is conserved (moves only) and type-1 count only grows, so the lost type must be type 0 (empty cells). Rule 1 alone is polluted bootstrap percolation: type 2 acts as blockers that shield some empty cells forever, so 3 types survive. Schelling moves relocate the blockers and expose those empty cells to threshold-3 growth until none remain. With no empty cells, Schelling dynamics cannot move, giving a frozen absorbing state. The interaction effect is real but trivial; the outcome follows directly from the rules. This overturns the first pass only in the sense of KNOWN instead of a novel variant, and adds the trivial/absorbing flag.
- Citations:
  - https://arxiv.org/pdf/1705.01652
  - https://www.math.ucdavis.edu/~gravner/papers/repl.pdf
  - https://www.nature.com/articles/s41598-022-23224-7
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC9653388
  - https://arxiv.org/pdf/0707.1681
- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________

