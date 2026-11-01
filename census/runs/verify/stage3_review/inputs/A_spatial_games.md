# Stage 3 (b) review inputs: A_spatial_games (25 phenomena)

## census:g1:grid-19

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.003 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_corrlen, g_activity_cv, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_corrlen (long spatial correlation length): median 5.5 with interaction, 1.0 with interaction knocked out (difference 4.5, robust spread 0.3706); difference at 4x area 5.0, at 3x run length 8.0
  - g_activity_cv (bursty activity): median 2.2848 with interaction, 0.0 with interaction knocked out (difference 2.2848, robust spread 0.0719); difference at 4x area 1.2061, at 3x run length 1.189
  - g_largest_domain (one cell-type fills most of the lattice): median 0.9585 with interaction, 0.0419 with interaction knocked out (difference 0.9166, robust spread 0.02); difference at 4x area 0.9383, at 3x run length 0.7355

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_corrlen 3.0; g_activity_cv 0.7329; g_largest_domain 0.3707
  - rules 2 only: g_corrlen 3.0; g_activity_cv 0.0; g_largest_domain 0.9926
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore neighbourhood plus self-play) with conformity-driven (majority) noise; near-absorbing single-strategy domination with perturbation-triggered avalanches

Trivial or artifact: True; confidence: medium

Twist (if any): 

Notes: Downgrades the first pass's VARIANT to KNOWN. Rule 2 alone gives 0.99 domination (near-absorbing frozen state). Rule 1 adds rare local perturbations that trigger short imitation cascades, which shows up as bursty activity whose burstiness score roughly halves at 4x area, consistent with independent local events rather than criticality. The correlation length mostly reflects the giant domain. The knock-out (static random colouring) is a weak baseline. The exact rule pair (global conformity at p=0.003 plus synchronous imitate-the-best with S=-0.5) was not found, but this is a parameter difference, not a new behaviour.

Citations:
- Nowak & May 1992, Nature, Evolutionary games and spatial chaos
- Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344
- Lim et al., Scale-invariant behavior in a spatial game of prisoners' dilemma, PRE 2002
- Szolnoki & Perc 2015, Conformity enhances network reciprocity, J R Soc Interface rsif.2014.1299
- Szolnoki & Chen 2018, Competition and partnership between conformity and payoff-based imitations, NJP, arXiv 1809.05079
- Asynchronous spatial evolutionary games, BioSystems (frozen static configurations)

Full reply text:

**Verdict: KNOWN.** I'd downgrade the first pass's "variant" call. Nothing in the behaviour needs new physics. It's a standard spatial Prisoner's Dilemma that freezes into one dominant type, plus a known conformity rule that adds rare noise and short bursts of activity.

**What the system is**
- **Rule 2 is the Nowak–May game.** Players sit on a square lattice and play one round against each of their 8 neighbours. Playing against oneself as well is one of the standard options in that model. Each cell holds one player, who plays one round of the Prisoner's Dilemma against each neighbour: either the four to the north, south, east and west, or the eight a chess king can reach. The update is "imitate the best": each player next plays whatever action earned the highest payoff in her neighbourhood, her own included. Your payoffs (S = −0.5 instead of 0) make this the "strong" version of the dilemma, a variant that has been well studied.
- **Rule 1 is conformity.** Conformists simply adopt whichever strategy is most common within their interaction range at any given time, regardless of the expected payoff. Payoff imitation combined with conformity has been modelled directly (Szolnoki & Perc 2015; Szolnoki & Chen 2018). One such model combines the payoff-driven and conformity-driven rules and always ends up in a frozen state.

**How each feature maps to known behaviour**
- **Largest domain (0.96):** With rule 2 alone, one type takes over 0.99 of the lattice, which is a near-absorbing state. Freezing under imitate-the-best is documented. For example, at b=1.1 every updating scheme settled into a static configuration. Rule 1 only slightly erodes that dominance.
- **Bursty activity:** Rare majority flips act as small local disturbances to a frozen game, and each one sets off a short cascade before the lattice freezes again. Lim et al. (2002) describe this "frozen state plus local perturbation gives avalanches" picture for the same model with randomness added: after a local perturbation, the critical states evolve through avalanches of many sizes.
- **Correlation length (5.5):** This mostly just measures the one giant domain. It is not a separate pattern.

**Trivial or artifact? Largely yes.**
- The main structure (one type filling the lattice) is close to a frozen absorbing state. The game alone produces it.
- The burstiness score is high because activity is almost zero, apart from occasional noise-triggered events. The score roughly halves at 4× area (2.28 → 1.21), which is what you'd expect from many independent local events averaging out. It doesn't look like scale-free avalanches. I haven't seen event-size distributions, so I can't rule out criticality.
- The interaction knock-out doesn't help settle this. It removes both rules, leaving a static random colouring with zero activity, so beating it is easy.

**Possible twist:** I didn't find this exact rule pair (conformity applied to everyone at probability 0.003, then imitate-the-best with synchronous updates and S = −0.5). It changes the parameters, not the behaviour, so I don't think it matters.

**Searches run:** Nowak–May imitate-the-best basics; noise and avalanches in the spatial dilemma; conformity combined with payoff imitation.

---

## census:novelty:grid-4

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.03 per step).

OBSERVED (features that passed the simulation checks):
  g_n_types, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
These grid features are measured on the transient: the stretch of the run before the lattice collapsed to a single type.
  - g_n_types (number of surviving types): median 1.0 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length
  - g_largest_domain (one cell-type fills most of the lattice): median 1.0 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_types None; g_largest_domain None
  - rules 2 only: g_n_types None; g_largest_domain None
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (unconditional imitation, Moore + self) mixed with voter-model imitation; finite-size consensus/fixation to an absorbing state

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: The two features (one surviving type, largest domain = 1.0) just describe finite-size collapse to all-C or all-D consensus. Vilone et al. 2014 study the same combination (unconditional imitation in the PD plus voter imitation with probability q) and report that voter dynamics strengthens consensus, i.e. uniform asymptotic states. At b=1.3 with self-interaction, Nowak-May is outside its chaotic regime (9/5<b<2) and tends to freeze; voter noise then drives coarsening to consensus. The knock-out comparison is uninformative because the interaction rules are what produce consensus. This is a frozen absorbing state, not an emergent pattern; the first pass's 'variant' label is not justified.

Citations:
- Vilone, Ramasco, Sanchez, San Miguel, 'Social imitation versus strategic choice, or consensus versus cooperation, in the networked Prisoner's Dilemma', PRE 90, 022810 (2014), arXiv:1405.4997
- Nowak & May, 'Evolutionary games and spatial chaos', Nature 359, 826 (1992)
- Szabo & Fath, 'Evolutionary games on graphs', Phys. Rep. 446, 97 (2007), arXiv:cond-mat/0607344

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass, which called it a variant. This behaviour is documented, and the reported features are just the system ending in a trivial absorbing state.

**What the system is.** Rule 1 is the Nowak–May spatial Prisoner's Dilemma. Each cell plays its 8 neighbours plus itself and copies the best scorer, with a weak dilemma: R=1, T=b=1.3, S=P=0. In that model, the strategy at the next step is the strategy within N(I) that receives the highest payoff. Rule 2 adds voter-model copying, at probability 0.03 per step.

**The same model has been published.** Vilone et al. (PRE 90, 022810, 2014; arXiv:1405.4997) study a Prisoner's Dilemma where decisions are made either based on the replication of the most successful neighbor's strategy (unconditional imitation) or by pure social imitation following an update rule inspired by the voter model. Each agent updates her strategy simply imitating a neighbor at random, that is, following the voter model (VM) dynamics with probability q. That is the same pair of rules as here. Their main findings:
- Voter dynamics produces an enhancement of the final consensus, i.e., asymptotic states are generally uniform.
- In finite systems, the fluctuations will eventually drive the system to consensus, which is an absorbing state.
- voter dynamics acts as an interface noise and has a similar effect to a pure random noise

Your two features (one surviving type; one type filling the lattice) are exactly this consensus outcome on a finite 64×64 grid.

**Why the rule-1-only run doesn't show it.** At b=1.3 the deterministic Nowak–May model is outside its chaotic regime. That regime is 9/5 < b < 2 with self-interactions included; outside it, the time-dependence of the strategy frequencies and the limit values depend on the initial configuration. Rule 1 alone therefore tends to freeze into a mixed pattern. The voter noise then pushes the boundaries around until one type takes over, which is the documented result.

**Is it a real pattern or a trivial outcome?** It's a trivial absorbing state. "One type" and "largest domain = 1.0" both just mean the grid has gone to all-cooperate or all-defect. The knock-out comparison doesn't show anything new. The interaction rules are what drive the grid to consensus, so of course the version without them never gets there. A pure voter model on a finite grid also reaches consensus eventually, so even the two-rule combination isn't needed for this outcome in principle. The game only affects how fast it happens and which type tends to win, and that is also covered by Vilone et al. I found no difference in the observed behaviour itself that would make this a variant.

**Searches run:** (1) Nowak–May 1992 regimes by b value; (2) spatial Prisoner's Dilemma with imitation plus voter noise, fixation and absorbing states.

---

## qd:004

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell adopts its neighbourhood's majority type (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  g_activity_trend, g_activity_cv, g_corrlen, g_n_types, g_reversals, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
These grid features are measured on the transient: the stretch of the run before the lattice collapsed to a single type.
  - g_activity_trend (activity rising): median -7.2083 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length
  - g_activity_cv (bursty activity): median 0.5469 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length
  - g_corrlen (long spatial correlation length): median 0.0 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length
  - g_n_types (number of surviving types): median 1.0 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length
  - g_reversals (cells flip back and forth): median 1.0 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length
  - g_largest_domain (one cell-type fills most of the lattice): median 1.0 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_activity_trend None; g_activity_cv None; g_corrlen None; g_n_types None; g_reversals None; g_largest_domain None
  - rules 2 only: g_activity_trend None; g_activity_cv None; g_corrlen None; g_n_types None; g_reversals None; g_largest_domain None
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial snowdrift/hawk-dove game with imitate-the-best (Nowak-May-type) plus majority-vote coarsening; fixation to a single-type absorbing state

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. The lattice collapses to one type (absorbing state). n_types=1, largest_domain=1 and corrlen=0 just describe that uniform lattice, and the negative activity trend is the decay into it. The knock-out and sub-program comparisons are against 'None' (no dynamics), so they don't show emergence. The first pass's 'variant' label overstates novelty: there's no twist in the phenomenology.

Citations:
- Nowak & May 1992, Nature 359:826 (Evolutionary games and spatial chaos)
- Hauert & Doebeli 2004, Nature 428:643 (Spatial structure often inhibits the evolution of cooperation in the snowdrift game)
- Killingback & Doebeli 1996, Proc R Soc B 263:1135 (Spatial evolutionary game theory: Hawks and Doves revisited)

Full reply text:

**Verdict: KNOWN.** What this system does is fixation into a single uniform type, an absorbing state. That's a trivial outcome, not a new emergent pattern, so I'm overturning the first pass's "variant" call, though only slightly.

**What the system actually is**
- **Rule 1 is the standard Nowak–May spatial game.** In their model, every agent plays the game with all his neighbors (once with each neighbor) and with himself (Moore neighborhood). Then everyone updates at the same time by "imitate the best neighbor": consider the set of all your neighbors plus yourself; then adopt the strategy of one of the agents in this set who has obtained the greatest payoff.
- **The payoffs make it a snowdrift game, not a Prisoner's Dilemma.** The ordering is T=1.5 > R=1 > S=0.5 > P=0, which is the snowdrift (hawk–dove) ordering. Hauert & Doebeli studied exactly this case: in unstructured snowdrift games, intermediate levels of cooperation persist. Unexpectedly, spatial structure reduces the proportion of cooperators for a wide range of parameters. On lattices, cooperators form only small filament-like clusters, and cooperation can die out (their paper shows lattice configurations near the extinction threshold of cooperators).
- **Rule 2 is a noisy majority-vote rule.** Majority rules smooth out domain boundaries and coarsen the lattice. On a finite 64×64 torus, that kind of dynamics is well known to end in consensus. Adding it to a dynamic that already pushes one strategy out will speed fixation up; it doesn't create anything new. I'm relying on general knowledge of majority-vote models here and didn't pull a specific citation.

**Why the observed features are trivial or artifacts**
- n_types = 1, largest_domain = 1.0 and corrlen = 0 are simply what a uniform lattice looks like. A correlation length of 0 on a lattice filled by one type means there's no domain structure left to measure, which is the opposite of "long spatial correlation."
- The features are said to be measured on the transient, yet n_types = 1 is the absorbing value. Either the measurement window includes the collapsed state, or the transient is very short. Both undercut the features.
- A negative activity trend (−7.2) is just the decay into a frozen state. The burstiness (CV ≈ 0.55) and a single reversal per cell fit a short, decaying transient.
- **The controls don't show what they appear to.** The knock-out "never produces this view" because without interaction the lattice never changes, so nothing is measured; it doesn't show emergence. The sub-programs return "None" too, so the 3-robust-spread test compared real numbers against missing values. That makes it meaningless.

**Searches run:** Hauert–Doebeli spatial snowdrift on a lattice; Nowak–May imitate-the-best spatial games. There's no remaining behaviour that needs a separate search: the outcome is fixation to one type, which is fully explained by the documented snowdrift-on-lattice results plus majority-rule coarsening.

---

## qd:005

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_activity_trend

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 3.6109 with interaction, 6.3716 with interaction knocked out (difference -2.7607, robust spread 0.1997); difference at 4x area -2.8673, at 3x run length -2.734
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.9397 with interaction, 0.0393 with interaction knocked out (difference 0.9004, robust spread 0.0196); difference at 4x area 0.929, at 3x run length 0.9053
  - g_activity_trend (activity rising): median -0.625 with interaction, 0.0 with interaction knocked out (difference -0.625, robust spread 0.0125); difference at 4x area -0.5871, at 3x run length -0.425

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 4.7406; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.8151; g_activity_trend -0.761
  - rules 2 only: g_n_domains 5.4702; g_corrlen 3.0; g_reversals 1.0; g_largest_domain 0.433; g_activity_trend -0.3092
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial snowdrift (hawk-dove) game with unconditional imitation mixed with voter-model social imitation; voter-driven coarsening / finite-size consensus (Vilone et al. hybrid UI+VM dynamics)

Trivial or artifact: True; confidence: medium

Twist (if any): 

Notes: Payoffs T>R>S>P make this the snowdrift game (Hauert-Doebeli parametrisation, r=0.5), not PD. Mixing best-neighbour imitation with voter copying is the Vilone et al. hybrid; the voter part acts as interface noise and drives finite-size symmetry breaking toward one dominant type. Observed few domains, ~94% one type and falling activity match that picture. Both rules preserve absorbing states and there is no mutation, so this is slow drift toward fixation in a 64x64 world. The knock-out is a frozen random colouring, so the knock-out contrast is trivially satisfied. Mislabelled features: 'activity rising' is negative (falling), 'long corrlen' is 2 lattice units, 'reversals' is a binary flag. The snowdrift-vs-PD/coordination difference is a rule difference, not a new phenomenology; each rule alone already shows most of the ordering.

Citations:
- Hauert & Doebeli, Nature 428, 643 (2004) doi:10.1038/nature02360
- Vilone, Ramasco, Sanchez, San Miguel, Phys. Rev. E 90, 022810 (2014) arXiv:1405.4997
- Vilone, Ramasco, Sanchez, San Miguel, Sci. Rep. 2, 686 (2012) 'Social and strategic imitation: the way to consensus'
- Roca, Cuesta, Sanchez, arXiv:0911.1720 (review, UI snowdrift on k=8 lattice)

Full reply text:

**Verdict: KNOWN.** I'm moving the first pass from "variant" to "known". The behaviour is a lattice slowly ordering toward one type taking over, which is a documented outcome. Mostly it's a finite-size drift toward an absorbing state (a fixed state the system can't leave), not a new pattern.

**What the system is.** With R=1, T=1.5, S=0.5, P=0 we get T > R > S > P. That is a snowdrift (hawk–dove) game, not a Prisoner's Dilemma. These numbers match the setup used by Hauert & Doebeli with cost-to-benefit ratio r = 0.5. Rule 1 is "copy the best-scoring neighbour" (unconditional imitation) with 8 neighbours. In unstructured snowdrift games, intermediate levels of cooperation persist. Unexpectedly, spatial structure reduces the proportion of cooperators for a wide range of parameters. The same game on an 8-neighbour lattice under unconditional imitation is covered in later reviews. One finds that with unconditional imitation the crossover between both behaviors takes place at T ≈ 1.7.

**The mix of the two rules is also published.** Vilone, Ramasco, Sánchez & San Miguel (PRE 90, 022810, 2014) study exactly this kind of hybrid. After every interaction, an individual updates her strategy simply imitating a neighbor at random, that is, following the voter model (VM) dynamics with probability q. Otherwise, with probability 1-q, an update rule based on the actual fitness attained by herself and her neighbors is implemented. They use unconditional imitation as the payoff-based rule, with a Prisoner's Dilemma. They report that voter dynamics acts as an interface noise and has a similar effect as a pure random noise. Their earlier paper (Sci. Rep. 2012) mixes the same two rules with a coordination game. On low-dimensional lattices it finds that competition between neighboring domains with opposite choices then takes place until eventually the fluctuations and finite system size lead to a symmetry breaking and to the selection of a single option as dominant. The introduction of the voter dynamics is the main responsible for the symmetry breaking, irrespective of how low q is in so far as it is non-zero. That matches what you observed: few domains, one type filling about 94% of the lattice, and activity falling.

**The only twist** is that the payoff game here is snowdrift, while those papers use PD or coordination games. I don't think it changes the phenomenology. Random-neighbour copying alone (rule 2 only) already coarsens, with a largest domain of 0.43. Imitation alone (rule 1 only) already reaches 0.82 dominance. Combining them just speeds up ordering toward one type, as Vilone et al. describe.

**Trivial or artifact?** Partly. It's not a numerical blow-up, but it does look like drift toward a frozen absorbing state:
- Neither rule can bring back a type once it's gone, so a largest domain of 0.94 is a system partway to full takeover in a 64×64 world.
- The knock-out comparison means little. With both rules reading neighbours, removing them leaves a frozen random colouring. Almost any dynamics would differ from that by "3 robust spreads".
- Three feature labels don't match the numbers:
  - "Activity rising" is negative (−0.625), so activity is actually falling, as expected while the system orders.
  - "Long correlation length" is 2 lattice units, which is short.
  - "Reversals" is a 0/1 flag. It shows flipping happens, not that a sustained pattern of flipping back and forth exists.

**Searches I ran:** Hauert–Doebeli spatial snowdrift; spatial snowdrift with voter-model or mixed imitation and coarsening; Vilone et al. on social versus strategic imitation. I didn't find a paper that runs snowdrift plus voter-style copying in particular.

---

## qd:012

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) neighbouring cells swap contents (every step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_reversals, g_activity

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 5.2013 with interaction, 6.3639 with interaction knocked out (difference -1.1626, robust spread 0.2313); difference at 4x area -1.2413, at 3x run length -1.0782
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_activity (many cells change per step): median 0.4514 with interaction, 0.0 with interaction knocked out (difference 0.4514, robust spread 0.009); difference at 4x area 0.4542, at 3x run length 0.4568

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 4.7621; g_reversals 1.0; g_activity 0.0032
  - rules 2 only: g_n_domains 6.3535; g_reversals 1.0; g_activity 0.3345
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial snowdrift (hawk-dove) game with deterministic best-neighbour imitation (Nowak-May lattice dynamics) plus random mobility/strategy diffusion

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. The swap rule alone produces most of the activity (0.33 of 0.45), and g_reversals is saturated (1.0) in both sub-programs, so it carries no information. The genuine interaction effect is modest: mobility keeps the otherwise near-frozen imitation lattice churning, and imitation slightly coarsens domains (5.2 vs 6.4). Both are documented mobility-plus-imitation phenomenology. Not a numerical artifact, but largely driven by the non-game mixing rule. Upgraded from KNOWN-VARIANT because there is no phenomenological twist beyond the parameter choice.

Citations:
- Nowak & May, Evolutionary games and spatial chaos, Nature 359:826-829 (1992)
- Sicardi, Fort, Vainstein & Arenzon, Random mobility and spatial structure often enhance cooperation, J. Theor. Biol. 256:240-246 (2009)
- Vainstein, Silva & Arenzon, Does mobility decrease cooperation?, J. Theor. Biol. 244:722-728 (2007)
- Vainstein & Arenzon, Spatial social dilemmas: dilution, mobility and grouping effects with imitation dynamics, Physica A 394 (2014)

Full reply text:

**Verdict: KNOWN.** This is a standard, well-studied model: a spatial game on a lattice where cells copy their best-scoring neighbour, plus random movement. I'm moving the first pass from "variant" to "known" because I couldn't find any behaviour here that isn't already described.

**What the rules are**
- **Rule 1 is the Nowak–May setup.** Cells are pure cooperators or defectors on a grid. Each plays its 3×3 block (8 neighbours plus itself) and copies the best scorer. That is exactly the spatial prisoner's dilemma territory of Nowak & May, "Evolutionary games and spatial chaos," Nature 359, 826–829 (1992). Nowak and May showed that cooperators and defectors can coexist indefinitely in this kind of setup. For some payoffs the result is a world constantly in flux, with each type invading the other's regions.
- **The payoffs make it a Snowdrift (Hawk–Dove) game, not a Prisoner's Dilemma.** With T=1.5 > R=1 > S=0.5 > P=0, defecting against a defector scores worst. That ordering defines the Snowdrift game, which is also well studied on lattices.
- **Rule 2 (neighbour swaps every step) is random mobility.** Mobility combined with imitation in the Prisoner's Dilemma, Snowdrift and Stag Hunt games is studied directly by Sicardi, Fort, Vainstein & Arenzon. They report that random dilution and mobility may suppress the inhibiting factors of the spatial structure in the Snowdrift game. Vainstein & Arenzon later mapped the Prisoner's Dilemma and Snowdrift games on a square lattice under a synchronous, noiseless imitation dynamics with mobility. They found intermediate values that may increase the amount of cooperators in the system and new phases appear. The earlier paper is Vainstein, Silva & Arenzon, "Does mobility decrease cooperation?" J. Theor. Biol. 244 (2007).

**Checking the three features**
- **Activity (many cells change per step), 0.45:** this is mostly caused by the swapping itself. Swaps alone already give 0.33. Imitation alone settles into an almost frozen state (0.003), as deterministic Nowak–May lattices often do. So the "emergent" part is only the extra ~0.12 of churn that appears when moving cells keep breaking up the frozen clusters. That is the same mobility-plus-imitation mechanism the papers above study.
- **Reversals (cells flip back and forth):** this measure tells us nothing. It reads 1.0 in both sub-programs and 0 only when nothing moves at all. It just separates "something changes" from "nothing changes."
- **Number of domains, about 5.2 vs 6.4:** a modest but robust drop. Imitation joins regions into larger patches, mixing breaks them up, and the full program sits between the two (rule 1 alone gives 4.8, rule 2 alone 6.4). This is the expected outcome of clustering competing with diffusion, not a new pattern.

**Trivial or artifact?** Not an artifact: the field doesn't blow up numerically, isn't a random colouring, and doesn't get stuck. But it's largely trivial: most of the activity and reversal signal would appear even if the cells didn't interact. The real interaction effect is a small, ordinary combination of clustering and mixing.

**Searches run:** random mobility and spatial structure in the Snowdrift game (Sicardi/Arenzon), and Nowak & May 1992 spatial chaos. I checked the search-result abstracts and summaries, not the full papers.

---

## qd:013

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) each cell adopts its neighbourhood's majority type (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_type_entropy

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 3.6362 with interaction, 6.3725 with interaction knocked out (difference -2.7363, robust spread 0.2002); difference at 4x area -2.8329, at 3x run length -2.6123
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_type_entropy (cell types evenly mixed): median 0.8205 with interaction, 0.9999 with interaction knocked out (difference -0.1794, robust spread 0.0364); difference at 4x area -0.1836, at 3x run length -0.164

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 1.3863; g_corrlen 3.0; g_reversals 0.0; g_type_entropy 0.05
  - rules 2 only: g_n_domains 1.7918; g_corrlen 11.0; g_reversals 1.0; g_type_entropy 0.9773
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial prisoner's dilemma (imitate-best, Moore neighbourhood + self) with conformity/majority-rule updating (Szolnoki & Perc conformity-enhanced network reciprocity)

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Rule-level differences (all cells update at once, every cell conforms with probability 0.1 each step instead of a fixed conformist fraction, S=-0.5) do not change the behaviour, which is documented. Mild overturn: VARIANT -> KNOWN, because no new behaviour was identified. The effects are real but small (correlation length 2 vs 1; reversals looks like a binary flag). Full papers were not opened, so parameter regimes were not checked.

Citations:
- Nowak & May 1992, Nature 359:826 (spatial chaos, C/D coexistence in fluctuating clusters)
- Szolnoki & Perc 2015, J. R. Soc. Interface 12:20141299, arXiv:1412.4113 (conformity in spatial social dilemmas; surface tension, smooth domain interfaces)
- Szolnoki & Perc 2016, Sci. Rep. 6:23633 (conformists; strategy-neutral state at full conformity)
- Chaos Solitons Fractals 2018 S0960077918305629 (time-alternating payoff-driven/conformity-driven updating)

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass's "variant" label, though only mildly. What you observed is cooperators and defectors coexisting in clusters that keep forming and breaking up, with more structure than random and cells switching back and forth. That is documented behaviour for the spatial prisoner's dilemma, both on its own and with a conformity rule added.

**What the system is.** Rule 1 is the Nowak–May spatial prisoner's dilemma. Each cell plays its 8 neighbours and itself, then copies the best scorer. The payoffs (T=1.5 > R=1 > P=0 > S=−0.5, with 2R > T+S) make it a proper prisoner's dilemma. Rule 2 is a noisy majority vote. A spatial prisoner's dilemma where players can follow the local majority is the "conformity" model of Szolnoki & Perc.

**Sources:**
- **Nowak–May:** Nowak and May observed that the deterministic "imitate-the-best" update rule generates chaotic dynamics in the spatial Prisoner's Dilemma game. Simple local interactions give rise to complex and persistent spatial patterns, allowing cooperators and defectors to coexist through the constant formation and breakup of clusters. In the standard replication, clusters of both types of players keep expanding, colliding and fragmenting. Your "reversals" and "many domains" features are this.
- **Conformity added:** Szolnoki & Perc (J. R. Soc. Interface 12, 20141299, 2015; arXiv 1412.4113) found that an appropriate fraction of conformists within the population introduces an effective surface tension around cooperative clusters and ensures smooth interfaces between different strategy domains. In their model, conformists simply choose the strategy that is most prevalent in their interaction range rather than the best-paying one. They also note that if conformity is too strong, a strategy neutral state is reached. That matches your rule-2-only result: large single-colour domains (correlation length 11) and evenly mixed types.
- **Follow-up work** varies the same mix: the alternation of the strategy-updating rule between payoff-driven and conformity-driven rules with time. That is very close to your "imitate, then conform with probability 0.1 each step" scheme.

**Differences, all in the rules rather than the behaviour:**
- Every cell updates at the same time, and imitation is deterministic.
- Every cell conforms with probability 0.1 each step, instead of a fixed share of players always conforming.
- The sucker's payoff is negative.

None of these produces behaviour that hasn't already been described. The full system sits between the two single-rule cases: imitation alone almost wipes out one type (entropy 0.05), and majority vote alone coarsens into big domains. That is the payoff-versus-conformity interplay the literature describes.

**Genuine or artifact?** It's a genuine interaction effect, not a random colouring or a frozen state. It differs from the no-interaction control and from both single-rule versions, and it holds up at 4× area and 3× run length. The effects are small, though: correlation length 2 versus 1, and "reversals" looks like a yes/no flag.

**Limits of the check:** I worked from search-result text and didn't open the full papers, so I can't confirm their exact parameter regimes match yours.

---

## qd:015

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a random network of 200 nodes (average 4 links each). Nodes have 2 types (random start). Rules, applied in order each step: (1) nodes play a two-strategy game with their neighbours (T=1.9, S=0.0) and copy the best-scoring neighbour; (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (with probability 0.003 per step).

OBSERVED (features that passed the simulation checks):
  n_oscillation, n_edge_turnover

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - n_oscillation (node states oscillate over time): median 0.3337 with interaction, 0.0 with interaction knocked out (difference 0.3337, robust spread 0.018); difference at 4x area 0.3341, at 3x run length 0.4021
  - n_edge_turnover (edges rewire): median 0.3355 with interaction, 0.0 with interaction knocked out (difference 0.3355, robust spread 0.0119); difference at 4x area 0.3645, at 3x run length 0.4734

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: n_oscillation 0.1274; n_edge_turnover 0.0
  - rules 2 only: n_oscillation 0.0; n_edge_turnover 0.3833
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Coevolutionary (adaptive-network) Prisoner's Dilemma with copy-the-best-neighbour updating and partner switching (Zimmermann-Eguiluz / Nowak-May lineage)

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Overturns VARIANT to KNOWN: no twist in the behaviour itself was found. Oscillation is the known non-settling behaviour of synchronous copy-the-best weak PD at T=1.9; it is already present with the game alone (0.127), and rewiring raises it modestly (0.33). Edge turnover nearly reads out rule 2 directly (rule 2 alone gives 0.38, more than the full model's 0.34), so it is close to trivial rather than emergent. The knock-out test cannot show emergence because both features are interactions by definition. Unclear point: whether 'type' is the strategy (then it is the documented mismatched-link rewiring model) or a fixed label (then it is two documented processes running side by side). Oscillation is a real dynamical pattern, not a numerical artifact, but it is documented.

Citations:
- Zimmermann & Eguiluz, Phys. Rev. E 72, 056118 (2005)
- Zimmermann, Eguiluz & San Miguel, Phys. Rev. E 69, 065102 (2004)
- Abramson & Kuperman, Phys. Rev. E 63, 030901 (2001)
- Perc & Szolnoki, BioSystems 99, 109 (2010)
- Evolutionary prisoner's dilemma games coevolving on adaptive networks, J. Complex Networks 6, 1 (2018), arXiv:1702.05119
- Szabo & Fath, Phys. Rep. 446, 97 (2007)

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass's "variant" verdict, partly. Neither observed behaviour is new. The 4x-area and 3x-run-length checks show the effects are robust, but robust doesn't mean undocumented.

**What the system is.** This is a coevolutionary Prisoner's Dilemma on an adaptive network, which has been studied for over 20 years. The weak-PD payoffs (T=1.9, S=0), synchronous play and copying the best-scoring neighbour are the Nowak–May setup. Moved onto networks, it's the standard Zimmermann–Eguíluz setup: each agent plays a Prisoner's Dilemma game with its local neighbors, collects an aggregate payoff, and imitates the strategy of its best neighbor, and the network is rewired at the same time. Variants that rewire based on mismatched neighbours are also documented. In one, the edge chosen is between a cooperator and a defector, and rewiring occurs with probability 1 − w with cooperating node X dismissing its defecting neighbor Y and rewiring to a random node Z. The topic has a well-cited review (Perc & Szolnoki 2010), which notes that coevolutionary rules may affect the interaction network, among other properties.

**Feature by feature:**
- **Node states oscillate.** With the game alone, states already oscillate (0.127). This is the expected result: copy-the-best PD at high temptation keeps shifting between cooperator and defector patterns rather than settling. On networks, different topologies, ranging from regular lattices to random graphs, produce a variety of emergent behaviors. Others report very high values resulting great oscillations in cooperation frequencies at temptations of 1.7–1.95. Adding rewiring raises oscillation to 0.33. That fits the known finding that rewiring keeps the network from freezing, so I don't see this as a new kind of behaviour.
- **Edges rewire.** This is close to trivial: it reads out rule 2 almost directly. Rule 2 alone already gives turnover of 0.38, higher than the full model's 0.34. The game somewhat reduces turnover; it doesn't create it.

**Unclear point about "types".** I can't tell from the description whether a node's type is its game strategy or a separate fixed label.
- If it's the strategy, this is essentially the documented model where mismatched links get rewired.
- If it's a fixed label, rewiring is a slow sorting process that doesn't depend on the game. The model is then two documented processes running side by side, and the only connection between them is that rewiring changes who plays whom. That only strengthens the "known" verdict.

**Weak evidence of interaction.** The knock-out test removes every rule that reads other nodes. Both features are interactions by definition, so of course they drop to 0; that test can't show anything emergent. The more useful comparison is with the single-rule runs. The only sign of the two rules affecting each other is that rewiring boosts oscillation by about 2.6x. That's a modest amount, in the direction you'd expect.

**Searches I ran:** Zimmermann–Eguíluz coevolutionary PD with copy-the-best; Perc–Szolnoki coevolutionary games review; Abramson–Kuperman PD on small-world networks with oscillations; PD with rewiring of mismatched links.

---

## qd:016

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell adopts its neighbourhood's majority type (with probability 0.01 per step).

OBSERVED (features that passed the simulation checks):
  g_activity_trend, g_largest_domain, g_n_types, g_reversals, g_n_domains, g_spec_peak

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_activity_trend (activity rising): median -4.2082 with interaction, 0.0 with interaction knocked out (difference -4.2082, robust spread 0.0842); difference at 4x area -4.2156, at 3x run length -2.7438
  - g_largest_domain (one cell-type fills most of the lattice): median 1.0 with interaction, 0.056 with interaction knocked out (difference 0.944, robust spread 0.0211); difference at 4x area 0.9802, at 3x run length 0.944
  - g_n_types (number of surviving types): median 1.0 with interaction, 2.0 with interaction knocked out (difference -1.0, robust spread 0.06); difference at 4x area -1.0, at 3x run length -1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_n_domains (many separate domains): median 0.6931 with interaction, 6.3552 with interaction knocked out (difference -5.6621, robust spread 0.141); difference at 4x area -7.0096, at 3x run length -5.6621
  - g_spec_peak (a dominant temporal frequency): median 0.0931 with interaction, 0.0 with interaction knocked out (difference 0.0931, robust spread 0.0039); difference at 4x area 0.1118, at 3x run length 0.0938

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_activity_trend -0.806; g_largest_domain 0.8267; g_n_types 2.0; g_reversals 1.0; g_n_domains 4.7318; g_spec_peak 0.0269
  - rules 2 only: g_activity_trend -1.8779; g_largest_domain 0.3944; g_n_types 2.0; g_reversals 1.0; g_n_domains 2.7081; g_spec_peak 0.3098
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial snowdrift/hawk-dove game with imitate-the-best (Nowak-May) plus conformity/majority-rule updating; coarsening to a frozen all-C or all-D absorbing state

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. The end state is a single-type, frozen absorbing monoculture: n_types=1, largest_domain=1, activity falling. The reversals and spectral-peak scores are probably short-lived effects of the coarsening phase. The knock-out baseline is a static random colouring, which makes the robustness differences look bigger than they mean. Overturned only in degree: from variant to plain known/trivial.

Citations:
- Nowak & May 1992, Nature 359:826 (Evolutionary games and spatial chaos)
- Hauert & Doebeli 2004, Nature 428:643 (spatial snowdrift game)
- Evolutionary game dynamics of combining the payoff-driven and conformity-driven update rules, Chaos Solitons Fractals 2020 (S0960077920305427) - always absorbed into full-C or full-D frozen state
- Szolnoki & Perc, Conformity enhances network reciprocity in evolutionary social dilemmas, arXiv:1412.4113
- Szolnoki & Perc 2018, Competition and partnership between conformity and payoff-based imitations, New J. Phys.

Full reply text:

**Verdict: KNOWN.** The end state is a frozen state where one strategy has taken over the whole lattice. That outcome is already documented for this kind of system, and it is a trivial absorbing state rather than a new pattern. I'm moving the first pass's "variant" call to plain "known": I found no twist in the observed behaviour itself, as opposed to the rules.

**What the system is**
- **The game:** with T=1.5 > R=1 > S=0.5 > P=0, this is a snowdrift (hawk–dove) game, not a prisoner's dilemma. The two games that get the most attention are the Prisoner's Dilemma and the snowdrift game (also known as hawk–dove or chicken). The spatial version is well studied: "spatial structure reduces the proportion of cooperators for a wide range of parameters."
- **Rule 1** is the Nowak–May setup exactly: each cell plays all its neighbours and itself, then every cell at once adopts the strategy of whoever in that 3×3 block scored highest.
- **Rule 2** is a majority-vote ("conformity-driven") update. Mixing it with payoff copying is an established line of work: one paper looks at the case where a player follows the local majority, with both learning methods present and competing.

**Is the behaviour documented?**
Yes. A spatial game combining payoff-driven and conformity-driven updates was found to be "always absorbed into the frozen state". More precisely, it always ends either all cooperators or all defectors. The same group of papers also describes conformity pulling players into compact clusters with smooth borders, which then pushes one strategy to near-complete dominance. Your sub-program results match this: payoff copying alone, or majority voting alone, each keeps both types alive. Combining them gives coarsening into one uniform domain, and with deterministic copying a uniform lattice can never change again.

**Trivial or artifact?**
It's an absorbing-state outcome:
- One type survives, it fills the lattice, there is about one domain left, and activity falls toward zero.
- The "flip back and forth" and "dominant frequency" scores probably come from the early, short-lived coarsening phase, not from lasting oscillations. A dominant frequency of about 0.09 in a run that freezes is not real periodic behaviour.
- The knock-out comparison is weak. Without interactions the lattice just stays as its random starting colouring, so almost any interacting rule would beat it on these measures.

**Possible twist (minor):** the papers I found mostly use the prisoner's dilemma, not the snowdrift game. A snowdrift game on its own tends to let both strategies coexist, and here the majority rule overrides that. But conformity pushing the system to a single strategy is the documented mechanism. I didn't verify which strategy wins here, and I wouldn't count that as new behaviour.

---

## qd:023

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.002 per step).

OBSERVED (features that passed the simulation checks):
  g_n_types, g_n_domains, g_reversals, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_types (number of surviving types): median 1.0 with interaction, 2.0 with interaction knocked out (difference -1.0, robust spread 0.06); difference at 4x area -1.0, at 3x run length -1.0
  - g_n_domains (many separate domains): median 0.6931 with interaction, 6.3792 with interaction knocked out (difference -5.6861, robust spread 0.1414); difference at 4x area -5.9258, at 3x run length -5.6861
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 1.0 with interaction, 0.0425 with interaction knocked out (difference 0.9575, robust spread 0.0208); difference at 4x area 0.9746, at 3x run length 0.9575

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_types 1.0; g_n_domains 1.2425; g_reversals 0.0; g_largest_domain 0.9957
  - rules 2 only: g_n_types 2.0; g_n_domains 5.9441; g_reversals 1.0; g_largest_domain 0.1122
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore neighbourhood with self-interaction), donation-game payoffs, plus weak voter-model noise; outcome = fixation to a single-strategy absorbing state

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Rule 1 alone already fixes to one type (largest domain 0.996). Rule 2 (copy a random neighbour) cannot reintroduce an extinct type, so the end state is a frozen uniform absorbing state. g_reversals=1 most likely reflects flips during the transient, not a sustained dynamic. The knock-out comparison only contrasts fixation with slow voter-model coarsening. The first pass's 'variant' label overstates novelty: the different payoffs (S=-0.5, T=1.5) just place the system in the parameter region where one strategy takes over, a region already mapped in the follow-up literature.

Citations:
- Nowak & May 1992, Nature 359:826 (https://www.nature.com/articles/359826a0)
- Roca, Cuesta & Sanchez 2009, Phys Rev E 80:046106 (T-S plane outcomes under unconditional imitation on lattices)
- https://wisc.pb.unizin.org/agent-based-evolutionary-game-dynamics/chapter/iii-1/
- https://arxiv.org/pdf/0901.0355

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass's "variant" label. This is the spatial Prisoner's Dilemma of Nowak and May (1992) on a square lattice, with a small amount of voter-model copying added. What it ends up doing is falling into a single-type absorbing state. That outcome is documented and close to trivial.

**Why it's the same system.** Rule 1 is the Nowak–May model in its standard form. Each player plays the game with all his neighbors (once with each neighbor) and with himself (Moore neighborhood), then everyone updates at once by imitate the best neighbor. The only change is the payoffs. Nowak and May used the "weak" Prisoner's Dilemma, with CD-payoff = DD-payoff = 0, CC-payoff = 1, DC-payoff = 1.85. Here the sucker's payoff is S = −0.5 and the temptation is T = 1.5. That is a strict Prisoner's Dilemma, specifically the donation game with benefit 1.5 and cost 0.5. Mapping how the outcome depends on T and S, for this update rule on lattices with 8 neighbours, is a standard follow-up (Roca, Cuesta & Sánchez 2009, PRE 80, 046106, Effect of spatial structure on the evolution of cooperation). The field as a whole grew out of Nowak and May's work, which opened the way to a large number of studies focused on different games, different evolutionary rules, and different lattices or networks. I didn't open Roca et al. to confirm which type wins at exactly these payoffs, but it doesn't change the verdict.

**What was actually observed:**
- **One surviving type, largest domain = 1.0, number-of-domains measure ≈ ln 2.** Rule 1 alone already gets there (1 type, largest domain 0.996). It's the well-known fixation into all-cooperators or all-defectors that this update rule produces once S < 0 shifts the balance. The famous "spatial chaos" coexistence happens in a different payoff regime and isn't what we see here.
- **Rule 2 can't sustain anything.** "Copy a random neighbour" can't create a new type. Once the lattice is uniform it does nothing at all, and in a uniform state there is no game-driven change either. So the end state is a frozen absorbing state, not a noise-maintained steady state.
- **The reversals feature (1.0) should be treated with suspicion.** The final state is a single type, and neither rule can reintroduce the missing one, so no flipping can continue. The flag most likely counts flips during the transient before fixation, or it's a measurement artifact. The game-only sub-program shows 0 reversals, and the rule-2-only sub-program shows 1, which points to transient voter-style churn rather than any lasting oscillation.
- **The knock-out comparison isn't informative.** Removing the interaction leaves rule 2's voter copying, which coarsens very slowly and still has many domains at this timescale. "Interaction causes a single domain" is just the difference between fixation and slow coarsening.

**Is it genuine emergence or an artifact?** The effect itself is real: interaction drives the lattice to consensus. But the end state is a frozen, uniform absorbing state, which is the kind of trivial outcome the checklist warns about. Nothing here (no coexistence, clusters, waves or chaos) goes beyond the standard spatial Prisoner's Dilemma.

---

## qd:031

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) a cell of type t becomes type t+1 (cyclically) if at least 6 neighbours already are t+1 (with probability 0.03 per step).

OBSERVED (features that passed the simulation checks):
  g_corrlen, g_activity_cv, g_n_domains, g_n_types, g_reversals, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_corrlen (long spatial correlation length): median 3.0 with interaction, 1.0 with interaction knocked out (difference 2.0, robust spread 0.08); difference at 4x area 2.0, at 3x run length 2.0
  - g_activity_cv (bursty activity): median 1.7393 with interaction, 0.0 with interaction knocked out (difference 1.7393, robust spread 0.0348); difference at 4x area 1.2595, at 3x run length 1.4044
  - g_n_domains (many separate domains): median 1.6094 with interaction, 6.3708 with interaction knocked out (difference -4.7613, robust spread 0.1596); difference at 4x area -5.5216, at 3x run length -4.7613
  - g_n_types (number of surviving types): median 1.0 with interaction, 2.0 with interaction knocked out (difference -1.0, robust spread 0.06); difference at 4x area -1.0, at 3x run length -1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.9904 with interaction, 0.0533 with interaction knocked out (difference 0.937, robust spread 0.0209); difference at 4x area 0.9743, at 3x run length 0.937

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_corrlen 3.0; g_activity_cv 0.0; g_n_domains 1.7006; g_n_types 1.5; g_reversals 0.0; g_largest_domain 0.9893
  - rules 2 only: g_corrlen 2.0; g_activity_cv 0.0; g_n_domains 4.2185; g_n_types 2.0; g_reversals 0.0; g_largest_domain 0.2542
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore+self) with conformity/majority-rule noise: fixation/near-consensus of one strategy

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: With 2 types the 'cyclic' rule is just a stochastic 6-of-8 majority (conformity) rule. Rule 1 alone already yields largest domain 0.989 and correlation length 3, so the knock-out contrast overstates emergence. The coupling only removes residual minorities (1.5 -> 1 surviving types). High activity CV and reversals come from sporadic flips in a near-frozen, nearly uniform state, which is close to an absorbing state, so the outcome is trivial. First pass said VARIANT; it is plainly a known behaviour.

Citations:
- Nowak & May 1992, Nature 359:826 (doi:10.1038/359826a0)
- Szabo & Fath 2007, Evolutionary games on graphs, arXiv:cond-mat/0607344
- Grujic et al. 2010, PLoS ONE e13749 (imitate-the-best Moore lattice gives cooperation ~0 or ~1 depending on initial condition)
- Szolnoki & Perc 2015, J. R. Soc. Interface 12:20141299 (conformity enhances network reciprocity)

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass's "variant" label. What happens here is one strategy taking over the whole lattice in a Nowak–May spatial Prisoner's Dilemma, nudged by a local "follow the majority" rule. Both pieces are already described in the literature. The measured features mostly describe a nearly frozen, almost uniform end state, not a new pattern.

**What the rules really are**
- **Rule 1 is the standard Nowak–May model.** Players on a square lattice play themselves and their 8 neighbours, add up their payoffs, and then adopt the strategy of one of the agents in this set who has obtained the greatest payoff. One review describes this exact setup on a lattice with self-, nearest-, and next-nearest-neighbor interactions) as suggested by Nowak and May (1992). The payoffs here (S = −0.5) are a harsher version of their usual weak-PD choice (S = 0), but the rule is the same.
- **Rule 2 is not really cyclic.** With only two types, "become type t+1" just means "switch to the other type". So the rule says: with probability 0.03, switch if at least 6 of your 8 neighbours are the other type. That is a noisy majority (conformity) rule. Calling it cyclic suggests rock-paper-scissors dynamics that cannot happen with two types.

**Each behaviour is documented**
- **One type fills almost the whole lattice (largest domain 0.99, one surviving type).** A study that ran this exact imitate-the-best Moore-neighbourhood model found the asymptotic level of cooperation is either 0 or a large value close to 1, depending on the initial condition. That fits your rule-1-only run: largest domain 0.989, with a median of 1.5 surviving types, meaning some seeds keep a small minority and some don't.
- **What adding rule 2 does.** Mixing conformity into payoff-based imitation on lattices has been studied: conformists simply adopt whichever strategy is most common within their interaction range at any given time, regardless of the expected payoff. Such players are known to introduce an effective surface tension around cooperative clusters and ensures smooth interfaces between different strategy domains. In your runs this "surface tension" wipes out the leftover minority specks, so every seed ends with one type. That is the only clear effect of the coupling (1.5 → 1 surviving types).
- **Bursty activity and cells flipping back and forth.** In a nearly uniform lattice almost nothing changes, and occasionally a few minority cells flip. Synchronous imitate-the-best is also known to produce small oscillating clusters. A low average with rare spikes automatically gives a high coefficient of variation (1.74). And "reversals = 1.0" probably just means "at least one cell flipped back", which a few residual cells or the random majority flips easily produce.

**Trivial or genuine? Mostly trivial (close to a frozen absorbing state)**
- Rule 1 alone already gives the long correlation length (3.0), the near-total takeover (0.989) and about 1.7 domains. The knock-out comparison makes these look like effects of coupling only because the knock-out has no interaction at all: it is a random mix. The fair baseline is rule 1 alone.
- A correlation length of 3 cells on a 64×64 lattice that is 99% one type is short. It mostly reflects the last specks, not long-range order.
- The features that do depend on the coupling (bursty activity, reversals, one surviving type) are leftover effects of noisy majority flips on a nearly consensus state.

**Searches:** Nowak–May spatial chaos with Moore neighbourhood and self-interaction; spatial PD with a negative sucker's payoff and imitate-the-best; conformity/majority rules in the lattice PD (Szolnoki–Perc).

**Two caveats:**
- I didn't find a paper using exactly S = −0.5, T = 1.5 with this probabilistic 6-of-8 threshold. But the behaviour (one strategy takes over, conformity smooths domains and removes leftover minorities) is the standard result.
- The rule description doesn't say which type wins. Under these payoffs it is most likely defection, but the verdict doesn't depend on that.

---

## qd:039

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour; (2) a type-0 cell becomes type 1 if at least 4 neighbours are type 1 (with probability 0.03 per step).

OBSERVED (features that passed the simulation checks):
  g_corrlen, g_n_domains, g_activity_cv, g_reversals, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_corrlen (long spatial correlation length): median 6.5 with interaction, 1.0 with interaction knocked out (difference 5.5, robust spread 0.3706); difference at 4x area 10.0, at 3x run length 7.5
  - g_n_domains (many separate domains): median 1.0986 with interaction, 6.4036 with interaction knocked out (difference -5.305, robust spread 0.15); difference at 4x area -5.5743, at 3x run length -5.1611
  - g_activity_cv (bursty activity): median 0.8042 with interaction, 0.0 with interaction knocked out (difference 0.8042, robust spread 0.0161); difference at 4x area 0.5086, at 3x run length 0.5431
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.9661 with interaction, 0.0433 with interaction knocked out (difference 0.9227, robust spread 0.0214); difference at 4x area 0.9293, at 3x run length 0.8499

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_corrlen 3.0; g_n_domains 1.7006; g_activity_cv 0.0; g_reversals 0.0; g_largest_domain 0.991
  - rules 2 only: g_corrlen 0.0; g_n_domains 0.6931; g_activity_cv 0.0; g_reversals 0.0; g_largest_domain 1.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (imitate-the-best, Moore+self) with one-sided C->D noise (cf. SPD(P) stochastic variant): near-absorbing dominance with noise-seeded transient clusters

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Overturns 'variant' only mildly: the twist (C->D conversion gated by >=4 defector neighbours) changes where noise acts (interfaces only), not the kind of behaviour. Features are largely consequences of near-consensus (largest domain ~0.97); g_reversals appears binary. Knock-out removes both rules, leaving a frozen random field, so knock-out contrasts are uninformative. Strict-PD S=-0.5 assumed qualitatively like the weak PD per Nowak-May; the dominant type was not determined from the summaries.

Citations:
- Nowak & May 1992, Nature 359:826 (as summarised in Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344)
- Dynamic phase transition in Prisoner's Dilemma on a lattice with stochastic modifications, arXiv 0910.0955 (SPD(P): cooperators defect with prob p under unconditional imitation)
- Szabo, Vukov, Szolnoki 2005, Phase diagrams for PD on 2D lattices, arXiv cond-mat/0506433
- Helbing & Yu, The outbreak of cooperation among success-driven individuals under noisy conditions, arXiv 0903.4054

Full reply text:

**Verdict: KNOWN.** This is the classic Nowak–May spatial Prisoner's Dilemma with a small, one-sided "cooperator turns defector" noise term added. What you observe is the standard result for that kind of model: one strategy holds almost the whole lattice (about 97%), with small clusters that keep appearing, flipping back and getting wiped out. I'm overturning the first pass only slightly, from "variant" to "known", because I couldn't find any feature of the behaviour (as opposed to the rules) that is new.

**Rule 1 is Nowak–May exactly.** In Nowak and May's model, players on a square lattice play a Prisoner's Dilemma with all their neighbors, and the payoffs are summed up. After this, each player imitates the strategy of those neighboring players (including herself) who has scored the highest payoff. In most cases the authors considered the Moore neighborhood (involving nearest and next-nearest neighbors) with or without self-interactions. Your payoffs (S = −0.5 rather than 0) are the "strict" version of the game. That doesn't change the qualitative picture: Nowak and May found the weak version has the same qualitative properties as the typical version with c<0, at least when |c|<<1. At |S| = 0.5 I'm relying on that rather than having it confirmed. Your rule-1-only run settles to one type covering 99% of the lattice, which fits the frozen regime at low temptation.

**The full system (rule 1 plus one-sided noise) is also documented.** The SPD(P) model adds exactly this kind of noise to unconditional imitation: after each Monte Carlo step, each cooperator may choose to change his state to defector state with probability p, while each agent imitate the best (unconditional imitation) strategy of neighboring agents. In that model only one absorbing state is possible, namely, in which all players choose to defect at all times. Coexistence of D and C is considered as active phase or fluctuating phase. Helbing and Yu, and the noise studies reviewed by Szabó and Fáth, describe the same mix of dominance, mutation-seeded invasions and recovery. A related absorbing-state result: This system exhibits a phase transition from a mixed state of cooperators and defectors to a homogeneous one where only the defectors remain alive.

**The only twist is in the rules, and it probably doesn't matter.** Your conversion only fires when a cooperator has at least 4 defector neighbours, so it acts only at boundaries and around leftover defectors, never inside a solid cooperator block. That changes where the noise acts, not what kind of behaviour you get. I found no new behaviour from it: the bursts, back-and-forth flips and near-total dominance are the usual signature of noise acting on a nearly absorbing state.

**Trivial or genuine?** It's a real emergent effect, since neither rule alone produces bursts or flip-backs. But it's a modest one:
- **The lattice is close to uniform.** "Long correlation length" and "few domains" (a median of 1.10, which looks like the log of about 3 domains) mostly just reflect near-consensus.
- **"Reversals" looks like a yes/no measure.** Its median of exactly 1.0 suggests it says flips happen, not how many.
- **The knock-out comparison is weak.** Removing every rule that reads neighbours deletes both rules, so the knock-out is just the frozen random starting pattern. Every feature will differ from that, so the large gaps add little.

One thing I couldn't settle from the summaries alone is which type ends up dominant. That changes how you'd describe it ("cooperator sea with defector bursts" vs. "defector sea") but not the verdict.

**Searches run:** Nowak–May imitate-best with Moore neighbourhood and self-interaction; spatial PD with imitation, noise and mutation; unconditional imitation with mutations, avalanches and bursts.

---

## qd:040

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour; (2) a type-1 cell becomes type 0 if at least 2 neighbours are type 0 (with probability 0.03 per step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 3.1355 with interaction, 6.374 with interaction knocked out (difference -3.2385, robust spread 0.1902); difference at 4x area -3.3646, at 3x run length -3.5426
  - g_corrlen (long spatial correlation length): median 3.0 with interaction, 1.0 with interaction knocked out (difference 2.0, robust spread 0.08); difference at 4x area 2.0, at 3x run length 1.5
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.9426 with interaction, 0.0424 with interaction knocked out (difference 0.9003, robust spread 0.0197); difference at 4x area 0.9197, at 3x run length 0.8978
  - g_moran (like-valued cells cluster spatially): median 0.5727 with interaction, 0.0001 with interaction knocked out (difference 0.5726, robust spread 0.0169); difference at 4x area 0.5692, at 3x run length 0.5639

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 4.151; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.89; g_moran 0.5599
  - rules 2 only: g_n_domains 0.6931; g_corrlen 0.0; g_reversals 0.0; g_largest_domain 1.0; g_moran 0.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (weak PD, b=1.3, Moore neighbourhood + self-interaction, synchronous imitate-the-best): cooperator-dominated state with stable defector lines and blinkers

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Rule 1 alone is the textbook Nowak-May model and reproduces every feature in the sub-program test. Rule 2 (D->C conversion with probability 0.03 given at least 2 C neighbours) is a one-way bias toward cooperators that only shifts the magnitudes slightly. The first pass called this a variant; I upgrade it to plain KNOWN, a mild overturn. The pattern is genuine but near-static: a cooperator sea with defector lines and documented blinkers that account for the reversals. It is not an artifact, since the knock-out is just the random initial colouring.

Citations:
- Nowak & May 1992, Evolutionary games and spatial chaos, Nature 359:826 (https://www.researchgate.net/publication/216634494_Evolutionary_Games_and_Spatial_Chaos)
- Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344, Sec. 6.5
- Newth & Cornforth, Asynchronous spatial evolutionary games: spatial patterns, diversity and chaos, IEEE CEC 2007 (b=1.3: stable defector lines + blinkers)
- Sigmund, Univie notes on Nowak-May blinkers/gliders

Full reply text:

**Verdict: KNOWN.** This is the Nowak–May spatial Prisoner's Dilemma (1992) in its weak-dilemma regime with b = 1.3. The extra conversion rule only nudges the result a bit further toward cooperators.

**Why it's an equivalent system**
- Rule 1 is the original Nowak–May setup almost exactly. Players are on a square lattice and play only "always cooperate" or "always defect". The payoffs are R=1, T=b, S=P=0. They used this simplified payoff matrix, and the classic version includes self-, nearest-, and next-nearest-neighbor interactions (the 8 neighbours plus self here). The update is the same too: players adopt the strategy that has received the highest payoff among its neighbors including themself.
- Your own sub-program test shows that rule 1 alone produces every observed feature: 2 out of 3 parts in the largest domain, 2–3 times the knock-out correlation length, ongoing reversals, and clustering. Rule 2 just shifts the numbers slightly (largest domain 0.89 → 0.94, fewer domains, correlation length 2 → 3). Run on its own, rule 2 ends in a frozen all-cooperator state.

**The behaviour is documented at this exact b value**
- A comparison of update schemes reports that when b=1.3 and b=1.5, the synchronous and self-synchronizing schemes created spatial patterns with stable lines of defectors, and small sub-populations that "flip" between cooperation and defection (blinkers). That matches what you measured: one cooperator domain covering about 94%, a few separate defector structures, and cells that keep flipping back and forth.
- More generally, Nowak and May's cellular automata can display most of the gadgetry of Conway's Game of Life – periodic blinkers included. The basic finding is that spatial structure enabled the cooperators to form clusters on the square lattice and so protect themselves against the exploitation by defectors.

**Is there a twist? No real one.** Rule 2 is a one-way, probabilistic bias from defector to cooperator. It makes the cooperator domain bigger and thins out the defector lines, but it doesn't create any new kind of behaviour. Adding a cooperation bias or mutation to Nowak–May is common in the follow-up literature. The first pass called this a variant; I'd upgrade that to plain KNOWN, which is a mild overturn.

**Genuine pattern or artifact? A genuine pattern, but a near-static one.**
- It isn't an artifact. The knock-out removes every interaction, so it's just the random initial colouring. The full system isn't absorbing, since reversals are still happening late in the run. It also isn't numerically diverged.
- What it actually is: an almost frozen configuration of a well-known model, with oscillators. The "reversals" feature most likely comes from the documented blinkers rather than anything chaotic.

---

## qd:044

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell copies the type of a random neighbour (with probability 0.3 per step); (2) each cell adopts its neighbourhood's majority type (with probability 0.2 per step); (3) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 3.7954 with interaction, 6.3869 with interaction knocked out (difference -2.5915, robust spread 0.2036); difference at 4x area -2.6903, at 3x run length -2.8315
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.9025 with interaction, 0.0444 with interaction knocked out (difference 0.858, robust spread 0.0189); difference at 4x area 0.9005, at 3x run length 0.8796
  - g_moran (like-valued cells cluster spatially): median 0.6071 with interaction, -0.002 with interaction knocked out (difference 0.609, robust spread 0.0122); difference at 4x area 0.6144, at 3x run length 0.6174

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 3 rules):
  - rules 1 only: g_n_domains 5.3952; g_corrlen 3.0; g_reversals 1.0; g_largest_domain 0.527; g_moran 0.501
  - rules 2 only: g_n_domains 1.3863; g_corrlen 25.0; g_reversals 1.0; g_largest_domain 0.5886; g_moran 0.9517
  - rules 3 only: g_n_domains 4.7999; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.8257; g_moran 0.5412
  - rules 1, 2 only: g_n_domains 2.8888; g_corrlen 14.0; g_reversals 1.0; g_largest_domain 0.7977; g_moran 0.9046
  - rules 1, 3 only: g_n_domains 4.682; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.6965; g_moran 0.5485
  - rules 2, 3 only: g_n_domains 0.6931; g_corrlen 0.0; g_reversals 1.0; g_largest_domain 1.0; g_moran 0.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial snowdrift/hawk-dove game with imitate-best update (Nowak-May type) plus voter noise and conformity: near-dominance of one strategy with noise-maintained small minority clusters

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Downgraded from first-pass 'variant' to KNOWN: the only novelty is which rules are combined, not the behaviour. Rules 2+3 alone lock into a single type; voter noise (rule 1) keeps flipping cells and maintains small minority islands. Correlation length is short (2 cells), not long. The knock-out baseline (pure random noise) is weak, so beating it on clustering features is expected for any copy-a-neighbour or majority dynamics. Not an artifact, but a low-interest near-absorbing state kept alive by noise.

Citations:
- Nowak & May 1992, Nature 359:826 (Evolutionary games and spatial chaos)
- Hauert & Doebeli 2004, Nature 428:643 (Spatial structure often inhibits the evolution of cooperation in the snowdrift game)
- Szolnoki & Perc 2015, J. R. Soc. Interface 12:20141299 (Conformity enhances network reciprocity in evolutionary social dilemmas)

Full reply text:

**Verdict: KNOWN.** The behaviour is documented. I'm moving it one step down from the first pass's "variant" call, because I found no twist in the behaviour itself, only in which rules are combined.

**What the system is.** Rule 3 is the Nowak–May spatial game: 8 neighbours plus playing against yourself, then copy the best-scoring neighbour. Every agent plays the game with all his neighbors and with himself (Moore neighborhood), and the total payoff is the sum of these encounters. Then each agent adopts the strategy of one of the agents in this set who has obtained the greatest payoff. The payoffs (T=1.5 > R=1 > S=0.5 > P=0) make this a snowdrift (hawk–dove) game, not a prisoner's dilemma. In the snowdrift game, spatial structure reduces the proportion of cooperators for a wide range of parameters. The result is that cooperators form only small filament-like clusters, and a transition to the pure D phase (extinction of cooperators) is a standard finding. Rule 1 is voter-model copying, which here acts as noise. Rule 2 is conformity / majority vote. Mixing payoff-based imitation with conformity on lattices is an established topic. Szolnoki and Perc found that conformists introduce an effective surface tension around cooperative clusters and ensure smooth interfaces between different strategy domains, and models where the imitation mechanism consists of two combined rules: payoff-based imitation and conformity have been studied.

**What was observed:**
- One type covers about 90% of the lattice.
- The rest is many small minority islands (about 44 domains, if g_n_domains is a log count).
- The correlation length is short (2 cells), so calling it "long" is misleading.
- Moran's I (a measure of like cells clustering together) is moderate, about 0.6.
- Cells keep flipping back and forth.

This is the textbook picture for a noisy imitation game on a lattice: one strategy takes over and noise keeps small minority clusters alive. The sub-programs show how it works. With rules 2 and 3 only, the system locks into a single type (largest domain 1.0, Moran 0). That is a frozen state where everyone is one type. Adding the voter noise (rule 1) keeps re-creating small minority patches, so the system stays near that state without freezing. That is the documented behaviour of near-pure-defector (or near-pure-cooperator) spatial games with mutation or imitation noise. It is not a new kind of pattern.

**Trivial or artifact?** It is not a numerical artifact and not frozen, because the flipping continues. But it is low-interest. The knock-out comparison is weak: removing every interaction leaves random noise, so any copy-a-neighbour or majority rule would beat it on these five features. The "no sub-program matches" result just reflects tuning between the frozen 2+3 state and the noisier states with rule 1. Nothing new appears.

**Searches I ran:** Hauert–Doebeli on spatial snowdrift; Nowak–May imitate-best on the lattice; conformity/majority rule combined with spatial games (Szolnoki–Perc).

**Caveat:** I couldn't tell from these features which type wins. Hauert–Doebeli would predict defectors, but either way the behaviour is the same documented one.

---

## qd:060

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 0.3 with heading noise 3.14159 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) each agent turns up the gradient of field 0; (2) each agent advances the type of the cell under it by one (cyclically); (3) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_activity_trend, g_corrlen, g_reversals, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.4659 with interaction, 6.3483 with interaction knocked out (difference -1.8824, robust spread 0.2163); difference at 4x area -1.9168, at 3x run length -1.8939
  - g_activity_trend (activity rising): median -2.666 with interaction, 0.0 with interaction knocked out (difference -2.666, robust spread 0.0847); difference at 4x area -2.9408, at 3x run length -2.7764
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.864 with interaction, 0.0569 with interaction knocked out (difference 0.8071, robust spread 0.0184); difference at 4x area 0.8517, at 3x run length 0.8192

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 3 rules):
  - rules 1 only: g_n_domains 6.3483; g_activity_trend 0.0; g_corrlen 1.0; g_reversals 0.0; g_largest_domain 0.0569
  - rules 2 only: g_n_domains 6.3953; g_activity_trend 0.0001; g_corrlen 1.0; g_reversals 1.0; g_largest_domain 0.0386
  - rules 3 only: g_n_domains 4.7664; g_activity_trend -0.5022; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.8235
  - rules 1, 2 only: g_n_domains 6.3578; g_activity_trend -3.5294; g_corrlen 1.0; g_reversals 1.0; g_largest_domain 0.0505
  - rules 1, 3 only: g_n_domains 4.7664; g_activity_trend -0.5022; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.8235
  - rules 2, 3 only: g_n_domains 4.8402; g_activity_trend -0.0038; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.3419
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial game with imitate-the-best on an 8+1 Moore lattice (snowdrift/hawk-dove payoffs, T>R>S>P), plus agent-driven local strategy flips acting as noise

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Rule 3 alone reproduces largest_domain (0.82 vs 0.86), n_domains, corrlen and reversals. Rule 1 is inert for the cells (sub-program 1,3 is identical to 3), likely because the field has no source. Rule 2 is localized random flipping, i.e. noise. The negative activity trend comes mainly from rules 1+2 together and is likely a transient from the decaying initial field. The knock-out baseline (frozen random lattice) inflates apparent emergence. Feature labels are misleading: activity is falling, not rising, and corrlen is about 2 cells. Partly artifactual: near-uniform end state with an inert field. No new phenomenology found.

Citations:
- Nowak & May 1992, Nature 359:826 (spatial games, imitate-the-best)
- Hauert & Doebeli 2004, Nature 428:643, doi:10.1038/nature02360 (spatial snowdrift game)
- Grujic et al. 2010 PLoS ONE 5:e13749 (imitate-the-best goes to cooperation near 0 or near 1 depending on initial condition)
- Cuesta et al. arXiv:0805.2071 (Nowak-May 8+1 neighbourhood)

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass. Almost everything observed comes from rule 3, and rule 3 is a textbook model. The agents add only a known kind of perturbation, and the chemical field does essentially nothing.

**1. What rule 3 actually is.** Each cell plays its 8 neighbours and itself, and then copies the best scorer. That is the Nowak–May "imitate the best" lattice game. In that model agents on a square lattice play their Moore neighbourhood and copy the neighbour with the largest payoff in the previous round (including themselves). The "8+1" detail matches too: Nowak and May's simulations do not correspond exactly to a Moore neighborhood in so far as every player plays also with herself, i.e., it is an 8+1 neighborhood.

The payoffs here (T=1.5 > R=1 > S=0.5 > P=0) are not a Prisoner's Dilemma, though. They define a snowdrift (hawk–dove) game. That case is covered by Hauert & Doebeli (2004): "spatial structure reduces the proportion of cooperators for a wide range of parameters". Also, in the spatial snowdrift game cooperators form only small filament-like clusters.

With deterministic imitate-the-best, the population usually ends up mostly one strategy: the asymptotic level of cooperation is either 0 or a large value close to 1, depending on the initial condition. That fits g_largest_domain ≈ 0.86, a few domains, and a short correlation length.

**2. Your sub-program data says the same thing.**
- **Rule 3 alone** already gives largest domain 0.82, 4.8 domains, correlation length 2 and reversals = 1. The full program gives 0.86, 4.5, 2 and 1. The only feature rule 3 doesn't reproduce is the activity trend (–0.50 vs –2.67).
- **Rule 1 does nothing to the cells.** "Rules 1, 3" is identical to "rule 3" on every feature. The field seems to have no source, so it presumably diffuses flat and its gradient carries no information.
- **Rule 2** (agents flipping the cell under them) is localized random flipping, i.e. noise added to Nowak–May dynamics. It explains the reversals and the small shifts in domain statistics.
- **The activity trend** shows up mainly when rules 1 and 2 are combined (–3.53). The likeliest explanation is that agent motion changes on the decaying initial field, which reduces flips over time. That is a transient, not a cell-pattern phenomenon.

**3. The interaction knock-out is a misleading baseline.** Removing every interaction rule leaves the lattice frozen at its random start. That is why every feature looks strongly "emergent" against it. Measured against "rule 3 alone", most of the effect disappears.

**4. The feature labels don't match the numbers.** "Activity rising" has a negative value, so activity is actually falling. "Long correlation length" is 2 lattice cells. "Many separate domains" is about 4–5 domains with one covering 86% of the lattice. What the run shows is coarsening towards a near-uniform state, with agents keeping some cells flipping.

**5. Genuine or artifact?** Not trivial: the domain coarsening is a real, well-known effect of imitation dynamics. But part of it is artifact: the field is inert for the cells, and the end state is close to a single strategy with ongoing local perturbation. I didn't find a twist in the behaviour worth calling a new variant.

I ran two searches, one on Nowak–May imitate-the-best and one on the spatial snowdrift game. I didn't search specifically for "noise or random walkers perturbing a spatial game", but given rule 3 alone reproduces nearly everything, I don't think it would change the verdict.

---

## review:016

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). 1 field(s) with diffusion [0.2]. Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.9, S=0.0) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.3498 with interaction, 6.3404 with interaction knocked out (difference -1.9906, robust spread 0.2138); difference at 4x area -2.1213, at 3x run length -2.0442
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.6847 with interaction, 0.0475 with interaction knocked out (difference 0.6372, robust spread 0.0214); difference at 4x area 0.6621, at 3x run length 0.6554

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (weak PD, Moore neighbourhood + self, imitate-best), spatial chaos regime 1.8<b<2

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Exact Nowak-May parameters (R=1, S=P=0, T=1.9, 8 neighbours + self, deterministic imitate-best). Largest domain ~0.68 matches 1-0.318 (asymptotic defector fraction), i.e. the defector background, not one type taking over. Reversals and domain churn = documented spatiotemporal chaos. Correlation length 2 is short, not long-range. The diffusing field is inert (read by no rule). The knock-out (frozen random start) is a trivial baseline. Small 64x64 tori may enter periodic cycles (Wu et al.), which doesn't change the verdict. Upgraded from KNOWN-VARIANT to KNOWN because the system and parameters are identical to the original model.

Citations:
- Nowak MA & May RM (1992) Evolutionary games and spatial chaos, Nature 359:826-829, doi:10.1038/359826a0
- Nowak MA & May RM (1993) The spatial dilemmas of evolution, Int J Bifurcation Chaos 3:35-78
- Nowak & Sigmund (2000) Games on grids, in The Geometry of Ecological Interactions
- Wu, Fu & Wang (2023) Evolutionary games and spatial periodicity (arXiv:2209.08267)

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass's "variant" label. This system is the spatial Prisoner's Dilemma of Nowak & May (1992), run with exactly their parameters. Nothing in the observed behaviour needs a "variant" twist to explain it.

**Same system, same settings.** In the 1992 paper, each player "plays the game" with its immediate neighbours, and then each site is kept by its owner or taken by whichever neighbour scored highest that round . Your rule does the same, with payoffs R=1, T=b=1.9 and S=P=0. Nowak & May studied the case *P* = *S* , with players who interact with their eight neighbors and with themselves. Your b=1.9 falls inside their chaotic range, 1.8 < b < 2. The reason: when b > 1.8, a 2x2 or larger cluster of D will continue to grow at the corners, and when b < 2, a 2x2 or larger cluster of C will continue to grow.

**How each measured feature matches the paper:**
- **Reversals and many domains:** in this range, cooperators and defectors coexist indefinitely in a chaotically shifting balance, and the equilibrium frequency of cooperators is (almost) completely independent of the initial conditions. Clusters of both types keep expanding, colliding and fragmenting. That accounts for g_reversals = 1 and the persistent multiple domains.
- **Largest domain ≈ 0.68:** For many initial conditions, the frequency of cooperators converges to 0.3178... for 1.8 < T < 2.0. So defectors make up about 68% of the lattice. Your largest domain (0.6847) almost exactly matches 1 − 0.318, which suggests the "dominant domain" is the connected background of defectors. Note that it is *not* one type taking over the lattice.
- **Correlation length = 2:** this is short and fits the paper's small-scale chaotic patterns. Calling it "long-range" would be an overstatement.

**One caveat from the literature.** On small finite lattices, the chaotic regime with 1.8 < b < 2 can end up in periodic cycles after long transients , . A 64×64 run that has settled into a cycle would still show reversals, so this doesn't change the verdict.

**Genuine or artifact?** The pattern itself is a real emergent effect: it isn't a frozen state, a random colouring or a numerical blow-up. Two parts of the setup are artifacts, though:
- The diffusing chemical field is never read by any rule, so it does nothing.
- The knock-out baseline is trivial. With imitation removed, the random starting colouring just freezes, which is why the knock-out shows a correlation length of 1, no reversals and a largest domain of about 0.05. The 3-spread differences therefore only show that the imitation rule does something. They don't point to anything new.

---

## review:017

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.5, 0.5)) p=0.1  (compact grammar notation)

OBSERVED (features that passed the simulation checks):
  g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (imitate-the-best-neighbour), with partially asynchronous / stochastic updating; dynamic C/D coexistence ('spatial chaos')

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Downgrades the first pass from KNOWN-VARIANT to KNOWN: two-strategy imitate-best with update probability 0.1 is the Nowak-May spatial game with partial asynchrony, and persistent strategy flipping (dynamic coexistence) is its signature result. The g_reversals feature is weak evidence: the knock-out removes the only rule, so a score of 0 there is guaranteed, and a yes/no flip measure cannot tell real turnover apart from small blinking spots or flips caused by random tie-breaking. Interpreting g=(1.5,0.5) as PD payoffs (b=1.5) is an assumption; the verdict holds for any two-strategy imitate-best lattice game that has ongoing turnover.

Citations:
- Nowak & May 1992, Nature 359:826 'Evolutionary games and spatial chaos'
- Huberman & Glance 1993, PNAS 90:7716 'Evolutionary games and computer simulations'
- Szabo & Fath 2007, Phys. Rep. 446:97 'Evolutionary games on graphs' (arXiv cond-mat/0607344)
- Grilo & Correia, 'Asynchronous stochastic dynamics and the spatial Prisoner's Dilemma game' / 'Effects of asynchronism on evolutionary games' JTB 2011

Full reply text:

**Verdict: KNOWN.** The first pass called this a variant. I'm overturning that only in strength: nothing in the observed behaviour goes beyond the documented model. It is the standard spatial game where each cell copies its best-scoring neighbour, with partly random update timing. The one measured feature, cells still flipping late in the run, is that model's best-known result.

**What the rule is.** Each cell has one of two strategies and always copies the best-scoring cell around it, with probability 0.1 per step. That is the Nowak–May setup. Players sit on a square lattice and use one of two fixed strategies: AllD (to be denoted simply as D) and AllC (or simply C), i.e., unconditional defection and unconditional cooperation. After this, each player imitates the strategy of those neighboring players (including herself) who has scored the highest payoff. I'm assuming g=(1.5, 0.5) is a game payoff, probably the temptation b=1.5 in their simplified one-parameter matrix (which only contains one free parameter b). I couldn't confirm the notation, but the verdict doesn't depend on it.

**The flipping is documented.** Nowak and May reported that this rule can generate chaotically changing spatial patterns, in which cooperators and defectors both persist indefinitely (in fluctuating proportions about predictable long-term averages). A recent paper summarises it as cooperators and defectors that coexist through the constant formation and breakup of clusters. Cells changing strategy back and forth is exactly that.

**The p=0.1 timing is documented too.** Cells updating at different times is the change Huberman and Glance (PNAS 1993) studied. It has since been explored with "asynchronous stochastic dynamics", which allows the exploration of all the space between synchronous and sequential updating. Results differ on whether this helps or hurts cooperators. Lindgren and Nordahl stress that asynchronous updating has little effect on the simulations (although it excludes spiral waves). Either way, ongoing turnover under partly random updating is part of the known picture. I found no twist in the behaviour, so KNOWN-VARIANT doesn't fit.

**Genuine or trivial?** It's not an artifact: there's no numerical blow-up, no frozen state, no random colouring. But the evidence is weak:
* The knock-out removes the only rule there is, so nothing can ever change. A score of 0 there is guaranteed, and 1.0 against 0 shows only that the system isn't frozen.
* The measure appears to be a yes/no flag (median exactly 1.0). It doesn't tell real, roughly chaotic turnover apart from small blinking spots, or from flips caused only by random tie-breaking between equal scores.

The flipping is real but uninformative. Nothing here suggests anything new.

**Searches I ran:** Nowak–May 1992 spatial chaos with imitate-the-best, and Huberman–Glance 1993 asynchronous updating.

---

## review:021

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell copies the type of a random neighbour (every step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_activity

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.9972 with interaction, 6.3673 with interaction knocked out (difference -1.3701, robust spread 0.2273); difference at 4x area -1.4316, at 3x run length -1.2861
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_activity (many cells change per step): median 0.4623 with interaction, 0.0 with interaction knocked out (difference 0.4623, robust spread 0.0092); difference at 4x area 0.4621, at 3x run length 0.4619

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 5.7421; g_corrlen 1.0; g_reversals 1.0; g_activity 0.263
  - rules 2 only: g_n_domains 4.727; g_corrlen 2.0; g_reversals 1.0; g_activity 0.0031
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial snowdrift (hawk-dove) game with synchronous imitate-the-best updating plus voter-model imitation noise (voter-model perturbation / noisy Nowak-May-type lattice game)

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. The fine-grained small domains (corrlen ~2) fit known filament-like anti-coordinated snowdrift patterns. The voter step is standard imitation noise. Activity and reversals come from the voter step disrupting the anti-coordinated state and best-imitation restoring it, which is ordinary noisy steady-state behaviour. The 'long correlation length' is only 2 cells. The no-interaction knock-out is empty, since both rules read neighbours, so the comparison with each rule alone is the meaningful one. Overturned in the direction of KNOWN rather than VARIANT: no twist in the behaviour was identified. Genuine interplay of the two rules but low novelty.

Citations:
- Hauert & Doebeli, Nature 428, 643-646 (2004), Spatial structure often inhibits the evolution of cooperation in the snowdrift game
- Szabó & Fáth, Physics Reports 446, 97-216 (2007), Evolutionary games on graphs
- Cox, Durrett & Perkins, Voter Model Perturbations and Reaction Diffusion Equations, AMS (2013)
- arXiv:1102.4570 (summary of filament-like snowdrift clusters / anti-coordination)

Full reply text:

**Verdict: KNOWN.** I'm partly overturning the first pass. This isn't a new behaviour that resembles a known model. It is a standard spatial snowdrift game with imitation noise, and that setup is well documented.

**What the system actually is**
- **The game is a snowdrift game, not a Prisoner's Dilemma.** The payoffs are T=1.5 > R=1 > S=0.5 > P=0. Because S > P, this is the snowdrift (hawk–dove) ordering. The "cooperate/defect" labels in the rule description are misleading here.
- **Rule 2 is the classic deterministic rule:** everyone copies the best-scoring neighbour at the same time, as in Nowak–May lattice games. Hauert & Doebeli studied the snowdrift game on lattices with this kind of updating. They found that while in the spatial prisoner's dilemma cooperators can survive by forming large, compact clusters, in the spatial snowdrift game cooperators form only small filament-like clusters. This happens because the payoffs make it advantageous to adopt strategies that are opposite to neighboring strategies. That fits what you measured: a fine-grained mix of many small domains, with a correlation length of only about 2 cells.
- **Rule 1 is the voter model.** Copying a random neighbour every step is just the voter model, used here as noise that ignores payoffs. Combining voter dynamics with game-based selection is a studied framework (J. T. Cox, R. Durrett, E. A. Perkins, Voter Model Perturbations and Reaction Diffusion Equations). Noisy versions of imitation dynamics on lattices are also covered in the Szabó & Fáth review. It notes that the new (target) strategy can be imitated from the neighborhood or chosen randomly.

**How the features arise**
- **Rule 2 alone almost freezes** (activity 0.003). Deterministic best-neighbour imitation is known to settle into static or near-static patterns.
- **Rule 1 alone flickers** (activity 0.26). Updating the voter model for every cell at once makes many cells change each step.
- **Together, activity rises to 0.46.** The voter step pushes cells away from the snowdrift's favoured "opposite to my neighbours" pattern, and the best-imitation step pushes them back. That accounts for the back-and-forth flips (g_reversals) and the extra activity. I found no twist in the behaviour beyond this.

**Is it real or a trivial outcome?** It's not a trivial outcome in the strict sense: it isn't frozen, diverged, or a random colouring. The two rules genuinely interact, since the full program beats each rule on its own. But it's an ordinary noise-driven steady state, and several labels overstate it:
- "Long correlation length" means 2 lattice cells versus 1.
- The no-interaction knock-out is empty, because both rules read neighbours. So that comparison is weak. The comparison with each rule alone is the meaningful one.

---

## review:023

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) each cell copies the type of a random neighbour (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_activity_trend

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 3.6109 with interaction, 6.3716 with interaction knocked out (difference -2.7607, robust spread 0.1997); difference at 4x area -2.8673, at 3x run length -2.734
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.9397 with interaction, 0.0393 with interaction knocked out (difference 0.9004, robust spread 0.0196); difference at 4x area 0.929, at 3x run length 0.9053
  - g_activity_trend (activity rising): median -0.625 with interaction, 0.0 with interaction knocked out (difference -0.625, robust spread 0.0125); difference at 4x area -0.5871, at 3x run length -0.425

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 4.7406; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.8151; g_activity_trend -0.761
  - rules 2 only: g_n_domains 5.4702; g_corrlen 3.0; g_reversals 1.0; g_largest_domain 0.433; g_activity_trend -0.3092
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial snowdrift (hawk-dove) game with imitate-the-best updating (Hauert-Doebeli parameterisation, r=0.5) plus voter-model noise / voter-model perturbation; coarsening toward single-strategy dominance

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: The payoffs T=1.5,R=1,S=0.5,P=0 make this a snowdrift game, exactly on the Hauert-Doebeli line with r=0.5, not a prisoner's dilemma. Adding 10% neighbour copying is voter-model noise. The observed coarsening, ~94% single-type dominance, falling activity and noise-driven flips are what this documented system is expected to do. Moved from VARIANT to KNOWN because I found no new behaviour. The pattern is not a numerical artifact but is close to trivial: it is drift toward a near-absorbing single-type state, and the knock-out baseline is just a frozen random colouring. The feature labelled 'activity rising' is actually negative (activity falls). I did not check which strategy wins or how fast fixation happens.

Citations:
- Hauert & Doebeli 2004, Nature 428:643, doi:10.1038/nature02360
- Cox, Durrett & Perkins, Voter model perturbations and reaction diffusion equations, Asterisque 349 (arXiv:1103.1676)
- Durrett 2014, Spatial evolutionary games with small selection coefficients, EJP 19 (arXiv:1406.5876)
- The effect of dispersal and neighbourhood in games of cooperation, J. Theor. Biol. (S0022519308001112)

Full reply text:

**Verdict: KNOWN.** This is a standard spatial snowdrift game with imitate-the-best updating, plus voter-model copying as noise. Its behaviour (domains merging and growing until one type takes over most of the lattice) is documented. That's a step beyond the first pass's "variant", but only a small one.

**Mapping to the literature**
- **The payoffs are a snowdrift game, not a prisoner's dilemma.** The order is T > R > S > P (1.5 > 1 > 0.5 > 0), which defines the snowdrift game. In fact the values sit exactly on the line Hauert and Doebeli studied, with r = 0.5. One analysis notes they did not consider the whole ST-plane but only the line given by T = 1 + r, S = 1 - r. Their Nature paper (2004, vol. 428, p. 643) found that spatial structure reduces the proportion of cooperators for a wide range of parameters. They also found that in the spatial snowdrift game cooperators form only small filament-like clusters. Rule 1 on its own is therefore a well-studied system.
- **Adding random neighbour copying is a voter-model perturbation, which is well studied.** In that framework, when w=0 either version of the dynamics reduces to the voter model, a system in which each site at rate 1 changes its state to that of a randomly chosen neighbor. Cox, Durrett and Perkins give conditions for the existence of non-trivial stationary distributions, and for extinction of one type, and one of their worked examples is an evolution-of-cooperation model (a model of the evolution of cooperation of Ohtsuki, Hauert, Lieberman, and Nowak). Separately, adding mixing to the Hauert–Doebeli grid has been studied directly: the authors repeat the grid models of the PD and SD games with the same update technique by which the main results of Hauert and Doebeli (2004)... were obtained but with an added mixing (dispersal) effect. That rigorous theory assumes weak selection, whereas here selection is strong. Even so, the plain voter model's coarsening toward consensus on a finite torus is textbook behaviour.

**What the features show**
- The pattern is a few large domains, with one type covering about 94% of the lattice. The game alone reaches 0.82 and voter copying alone 0.43, so the combination fixes faster.
- Activity falls over time. The feature is labelled "activity rising", but the value is negative.
- Cells flip back and forth, which is expected from the 10% copying noise.

Two coarsening mechanisms are being added together. The game picks a winning type and voter drift speeds up its takeover on a 64×64 lattice. The full program differs from both sub-programs, but I see nothing new in that difference.

**Real or artifact?** It's not a numerical artifact, but it is close to trivial. The system is heading toward a near-absorbing, almost single-type state; on a finite lattice the voter noise makes this inevitable. The knock-out baseline is also weak: with no interaction, nothing changes, so the comparison is against a frozen random colouring. Beating that baseline doesn't make the pattern interesting. The "correlation length" goes from 1 to 2, which is minimal structure.

**Possible twist:** I didn't check which type wins (cooperators or defectors) or how fast fixation happens at r = 0.5 with Moore neighbourhood plus self-interaction. That is a parameter detail, not new behaviour.

---

## review:037

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.1 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=-0.5) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_type_entropy

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 3.6762 with interaction, 6.3775 with interaction knocked out (difference -2.7013, robust spread 0.2011); difference at 4x area -2.9673, at 3x run length -2.7276
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_type_entropy (cell types evenly mixed): median 0.8378 with interaction, 0.9998 with interaction knocked out (difference -0.162, robust spread 0.0368); difference at 4x area -0.1681, at 3x run length -0.1696

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 1.7918; g_corrlen 10.5; g_reversals 1.0; g_type_entropy 0.9207
  - rules 2 only: g_n_domains 1.4979; g_corrlen 3.0; g_reversals 0.0; g_type_entropy 0.0581
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial prisoner's dilemma (Nowak-May imitate-the-best) with conformity / majority-vote hybrid updating

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Both rules read neighbours, so the knock-out is just the frozen random initial state and is a weak control; the meaningful evidence is that the full program differs from each single-rule sub-program (majority-only coarsens, corrlen 10.5; game-only nearly absorbs, entropy 0.06). The combination gives the documented conformity-plus-reciprocity regime: biased but persistent coexistence (~27/73), small domains, boundary cells that keep flipping. Effects are modest (corrlen 1->2, integer-valued). Overturns first pass only in the sense of upgrading VARIANT to KNOWN: no phenomenological twist found.

Citations:
- Cui P.-B. & Wu Z.-X., Impact of conformity on the evolution of cooperation in the prisoner's dilemma game, Physica A 392, 1500-1509 (2013)
- Szolnoki A. & Perc M., Conformity enhances network reciprocity in evolutionary social dilemmas, J. R. Soc. Interface 12, 20141299 (2015)
- Szolnoki A. & Perc M., Leaders should not be conformists in evolutionary social dilemmas, Sci. Rep. (2016), arXiv:1604.06438
- Nowak M.A. & May R.M., Evolutionary games and spatial chaos, Nature 359 (1992) [from background knowledge, not retrieved]

Full reply text:

**Verdict: KNOWN.** This overturns the first pass's "variant" call. What you're seeing is a spatial prisoner's dilemma where players copy their best-scoring neighbour, mixed with a probabilistic majority ("conformity") rule. That hybrid has been studied directly, and the behaviour you measured falls inside what those studies report. I found no new phenomenology.

**The equivalent systems in the literature**
- **Rule 2 alone** is the classic Nowak–May spatial game: 8 neighbours plus self, copy the best scorer. I'm citing that from memory; neither search pulled it up. Your sub-program result fits it. With S = −0.5, one strategy nearly takes over (type entropy 0.058), which is close to a frozen end state.
- **Rule 1 alone** is the noisy majority-vote / voter-type model. Its signature is slow coarsening into large domains, which matches your correlation length of 10.5.
- **The combination has been modelled before:**
  - Cui & Wu (2013) built an evolutionary model of the prisoner's dilemma game taking into account two cognitive mechanisms: (1) payoff-biased strategy..., with conformity as the second mechanism. They found when the temptation to defect is low and both conformism and local interactions are present, the system can reach high levels of cooperation or even a full cooperation state.
  - Szolnoki & Perc (2015, J. R. Soc. Interface 12:20141299) found that a proper proportion of conformists in the population can form an efficient surface around cooperative clusters and guarantee smooth interfaces between diverse strategy domains of payoff-driven individuals.
  - Their follow-up even uses a 0.1 noise parameter: K=0.1 introduces some uncertainty to the process, such that it is not completely impossible for a conformist to adopt a strategy that is in the minority among its neighbors.
  - Normative conformity has also been reported to produce within-cluster behavioral homogeneity and between-cluster behavioral diversity.

**How your four features map onto this**
- **Fewer, somewhat larger domains** than random, but much smaller than under pure majority vote.
- **Short correlation length** (2): domain size is set by the competition between the two rules.
- **Cells flipping back and forth:** the majority rule keeps knocking cells off the payoff-imitation fixed points, so domain borders keep moving instead of freezing.
- **Uneven but persistent mix** (entropy 0.84, roughly a 27/73 split): conformity stops the game from driving one strategy to near-extinction.

This is the documented "conformity + network reciprocity" regime. I didn't find a paper using this exact payoff set with best-neighbour copying on this lattice. As you noted, that's a rule-tuple gap, not a behaviour gap.

**Genuine or artifact?** It's a genuine pattern, but a modest one.
- The interaction knock-out is a weak control here. Both rules read neighbours, so removing them leaves the frozen random starting state. The fact that the full program differs from it says little.
- The meaningful evidence is that the full program differs from both single-rule versions. It does, robustly, so the coexistence comes from the two rules interacting rather than from one of them.
- The effects are small: correlation length goes from 1 to 2 (whole-number steps), and domain count drops by about 2.7.
- It isn't diverged or random colouring. It is a known mixed, actively fluctuating steady state.

I ran two searches: (1) conformity + network reciprocity in the spatial prisoner's dilemma (Szolnoki & Perc), and (2) Cui & Wu 2013 conformity / majority rule in the prisoner's dilemma. I verified these papers through abstracts, citing papers and excerpts, not the full texts.

---

## review:041

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
C kC=2 | C[any] CNT_GE(b=0,th=3) -> IMITATE_BEST(g=(1.5, 0.5)) p=1.0  (compact grammar notation)

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.4006 with interaction, 6.357 with interaction knocked out (difference -1.9564, robust spread 0.2152); difference at 4x area -1.9715, at 3x run length -1.9626
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.8939 with interaction, 0.0475 with interaction knocked out (difference 0.8464, robust spread 0.0188); difference at 4x area 0.872, at 3x run length 0.8457
  - g_moran (like-valued cells cluster spatially): median 0.4901 with interaction, 0.0062 with interaction knocked out (difference 0.4839, robust spread 0.0193); difference at 4x area 0.5004, at 3x run length 0.4872

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma with deterministic imitate-the-best (unconditional imitation), b<1.8 regime where cooperators take over with leftover defector patches and blinkers

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: The count gate (CNT_GE th=3) only restricts which cells update; no new phenomenology. Effects are robust vs knock-out, but the correlation length only rises from 1 to 2 (short range), and binary reversals likely reflect known period-2 blinkers. Overturns VARIANT to KNOWN because no phenomenological twist was identified. Not verified: that the (1.5, 0.5) payoff tuple maps exactly onto Nowak-May (b,R,P,S), and the 1992 paper was not opened directly.

Citations:
- Nowak & May 1992, Nature 359:826, Evolutionary games and spatial chaos
- Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344, Sec 6.5
- Roca et al., arXiv 0901.0355
- arXiv 2607.24638 spectral cartography of spatiotemporal chaos in 2x2 games
- github giakoumoglou/game_theory_spatial_PD (b-regime summary)

Full reply text:

**Verdict: KNOWN.** I'd move this up from the first pass's KNOWN-VARIANT. The behaviour itself is already documented. The only new thing is the exact rule combination, and that doesn't count as a discovery.

**What the rule is.** It has two cell types, payoff parameters (1.5, 0.5), and a deterministic "copy the highest-scoring neighbour" update. That is the spatial Prisoner's Dilemma of Nowak & May (1992). They introduced evolution in this setup by using the unconditional imitation rule (also known as "imitate-the-best"), where each player chooses the strategy of the neighbor with largest payoff, provided this payoff is greater than the player's. Their payoff table has one free parameter, b (the temptation to defect). It only contains one free parameter b, but is expected to preserve the essence of the Prisoner's Dilemma. Here b is about 1.5. The extra gate (a cell only updates if at least 3 neighbours are type 0) just limits where updates happen. It doesn't add a new kind of behaviour.

**Why the behaviour is the documented one.** What happens depends on b:
- For b < 1.8, a big D cluster will shrink, and for b < 2, a 2x2 or larger cluster of C will continue to grow.
- In that range, cooperators take over most of the lattice, leaving small leftover defector patches and lines plus some cells that flip back and forth with period 2 ("blinkers").

That matches every measured feature:
- One type fills about 89% of the lattice (g_largest_domain).
- There are few domains, and they are strongly clustered (Moran's I ≈ 0.49).
- Cells keep flipping (g_reversals = 1), which fits the known oscillating boundary sites.

The mechanism is also known: cooperators survived by self-organizing in clusters, where the interactions within the clusters yielded larger payoffs to cooperators than those obtained by defectors at the boundaries of cooperators' clusters. Nowak & May's famous chaotic coexistence only appears at higher b (about 1.8–2). There, clusters of both types of players expand, collide and fragment. This run is below that range, in the regime where one type dominates. Recent work is still mapping these regimes as phase boundaries across the payoff plane: analytical instability thresholds for local motifs including invaders, cooperative pairs, stripe interfaces, and cooperative cores.

**Is it real or an artifact?** It's a genuine pattern from the interaction, not an artifact. The knock-out stays at random-noise values (Moran ≈ 0.006, largest domain ≈ 0.05), and the effects hold at 4× area and 3× run length. Two caveats, though:
- The correlation length only goes from 1 to 2. That is a short range, so it shouldn't be described as "long".
- "Reversals = 1" looks like a yes/no flag. It probably reflects the known period-2 blinker sites rather than large-scale reversals.

So this is a real effect, but it's mostly one type taking over plus small leftover patches. It isn't rich pattern formation.

I ran one search (Nowak–May imitate-the-best regimes). I didn't open the original 1992 paper directly; the b-regime details come from the secondary sources cited above. I also couldn't confirm that the (1.5, 0.5) tuple maps exactly onto Nowak–May's (b, R, P, S) payoffs. That's the main remaining uncertainty, and it's why confidence is medium-high rather than high.

---

## review:042

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
C kC=2 | C[any] CNT_LE(b=0,th=3) -> IMITATE_BEST(g=(1.9, 0.0)) p=0.1  (compact grammar notation)

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_activity_trend, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.3819 with interaction, 6.3456 with interaction knocked out (difference -1.9637, robust spread 0.2146); difference at 4x area -2.2018, at 3x run length -1.9513
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_activity_trend (activity rising): median -1.0836 with interaction, 0.0 with interaction knocked out (difference -1.0836, robust spread 0.0217); difference at 4x area -1.128, at 3x run length -1.0771
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial (weak) Prisoner's Dilemma with imitate-the-best, semi-synchronous/probabilistic updating; C/D cluster coexistence with fluctuating boundaries (spatial chaos regime, b=1.9)

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Read the rule as imitate-the-best weak PD with b=1.9 and p=0.1 partial updating. CNT_LE gate (update only if <=3 neighbours in state 0) changes the rule but not the kind of behaviour. Features (domains merging, short correlation length of 2, persistent flipping) match documented C/D cluster dynamics. Activity trend is negative, i.e. settling or slow drift (consistent with the asynchronous-update literature), so the 'activity rising' label is wrong. Correlation length 2 means short-range clustering only. Overturn is mild: first pass said variant; no new phenomenology found. Did not confirm whether this exact gated version ends in lasting coexistence or slow defector takeover; either outcome is documented.

Citations:
- Nowak & May 1992, Nature 359:826 (doi:10.1038/359826a0)
- Szabo & Fath 2007, Phys. Rep. 446:97 (arXiv cond-mat/0607344)
- Huberman & Glance 1993, PNAS 90:7716
- Prisoner's Dilemma with Semi-synchronous Updates, arXiv 0910.0961
- Spectral cartography of spatiotemporal chaos in 2x2 games, arXiv 2607.24638

Full reply text:

**Verdict: KNOWN.** This is a version of the Nowak–May spatial Prisoner's Dilemma, with partial (semi-synchronous) updating. The behaviour it shows is classic spatial-PD behaviour. I'm overturning the first pass's "variant" label, but only slightly. I didn't open the full texts; I checked the abstracts and secondary sources the search returned.

**How I read the rule.** There are two strategies (kC=2), and players copy the best-scoring neighbour ("imitate-the-best"). The payoffs g=(1.9, 0.0) look like Nowak–May's one-parameter "weak" Prisoner's Dilemma with temptation b=1.9. In each step only 10% of cells update (p=0.1), which makes the updating partly random rather than all-at-once. The CNT_LE condition (update only if at most 3 neighbours are in state 0) restricts which cells can update. That changes the rule, not the kind of behaviour you get.

**What's documented:**
- Nowak and May put players on a square lattice and used this exact copy rule. Cooperators survived by self-organizing in clusters, where the interactions within the clusters yielded larger payoffs to cooperators than those obtained by defectors at the boundaries (Nature 359, 826, 1992).
- The flipping and many-domain features are the known "spatial chaos" regime: cooperators and defectors coexist through the constant formation and breakup of clusters. b=1.9 sits inside the chaotic window. In a common write-up of their analysis, for b > 1.8, defector clusters of 2x2 or larger keep growing at the corners, while cooperator clusters of 2x2 or larger keep growing for b < 2. There is also a 2026 paper that maps where this chaotic phase begins and ends (arXiv 2607.24638).
- The paper uses the payoffs (0, b; 0, 1), the same weak dilemma as yours. It notes the weak version has the same qualitative properties as the typical version with c<0, at least when |c|<<1 (Szabó & Fáth 2007).
- How the update scheme matters is itself well studied. Huberman & Glance (PNAS 1993) found that with fully one-at-a-time updates the system rapidly converges to steady state where all the players become defectors. Partial updating, as here, has been studied directly (arXiv 0910.0961). Results depend a lot on how the randomness is introduced: it is the way how asynchronous behaviour is introduced that determines the persistence or breakdown of cooperation.

**Genuine pattern or artifact?** Real, but weak. It isn't a random colouring or a frozen state: with no interaction there is no flipping, and turning interaction on gives persistent flipping, fewer domains (merging) and a longer correlation length. Three caveats:
- The correlation length is only 2 lattice units, against 1 without interaction. That is short-range clustering, not long-range order.
- The "activity rising" feature actually has a *negative* value (about −1.08). Activity is falling over time, so the label is misleading. The data fit gradual settling of domains or a slow drift toward one strategy, which is the documented outcome under random updating. They don't show increasing activity.
- I couldn't confirm whether this particular gated, 10%-update version ends in lasting coexistence or slowly heads toward all-defectors. Whichever it is, that outcome is already documented in the update-scheme literature, so it doesn't add anything new.

---

## review:054

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.003 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_activity_trend, g_corrlen, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 3.3842 with interaction, 6.3716 with interaction knocked out (difference -2.9873, robust spread 0.1951); difference at 4x area -3.0383, at 3x run length -4.7621
  - g_activity_trend (activity rising): median -2.7171 with interaction, 0.0 with interaction knocked out (difference -2.7171, robust spread 0.0543); difference at 4x area -2.8581, at 3x run length -2.3989
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 4.2554; g_activity_trend -1.3135; g_corrlen 3.0; g_reversals 1.0
  - rules 2 only: g_n_domains 4.1026; g_activity_trend -1.6204; g_corrlen 2.0; g_reversals 1.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (weak PD, imitate-the-best, Moore + self-interaction), low-b regime (b<1.75) with weak majority-vote noise: cooperator cluster growth/coarsening relaxing to static or period-2 (blinker) configurations

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Mild overturn: first pass said KNOWN-VARIANT; no phenomenological twist found, so KNOWN. At b=1.3 C clusters grow and D clusters shrink (documented thresholds b<1.8, b<2); deterministic imitate-best with self-interaction provably ends in fixed or periodic states, which explains falling activity (the 'activity rising' label has the wrong sign) and blinker reversals. Correlation length 2 is short, not long. Knock-out removes both rules (both read neighbours), leaving a frozen random field, so the contrast is against a trivial baseline. Rule 1 (p=0.003 majority vote) only adds slight smoothing and rare perturbations; each sub-program reproduces the same qualitative features. Genuine dynamics, but a relaxational, near-absorbing textbook outcome, not novel emergence. Sources verified via search snippets only; the Nature paper itself was not opened.

Citations:
- Nowak & May 1992, Nature 359:826 'Evolutionary games and spatial chaos' https://www.nature.com/articles/359826a0
- Wu, Fu & Wang, 'Evolutionary games and spatial periodicity' arXiv:2209.08267
- Nowak & May 1993, Int. J. Bifurc. Chaos 3:35 'The spatial dilemmas of evolution'
- Szabo & Fath 2007, Physics Reports 446:97 'Evolutionary games on graphs'

Full reply text:

**Verdict: KNOWN. This is the Nowak–May spatial Prisoner's Dilemma at low temptation, with a little majority-vote noise added. I'm moving it from the first pass's KNOWN-VARIANT to KNOWN.**

**Why it matches Nowak–May**
- Rule 2 is Nowak–May's model almost exactly. Each cell plays the game with its 8 neighbours and itself, and each lattice-site is occupied by the player with the highest score among the previous owner and the immediate neighbors. The payoffs (R=1, T=b=1.3, P=S=0) are their "weak dilemma" form, which has a single parameter.
- The behaviour at b=1.3 is documented. Nowak and May describe how (b < 1.8) big D cluster will shrink; (b < 2) 2x2 or larger cluster of C will continue to grow. So at b=1.3, groups of cooperators expand, defectors are pushed into small remnants, and the run settles.
- What it settles into is also proven. For deterministic "imitate the best" updating with self-interaction, the process after a transient phase of finite generations, either be stabilized at one state, or enter periodic cycle almost everywhere (Wu, Fu & Wang).

**How each observed feature maps onto that**
- **Activity falling (g_activity_trend about −2.7):** this is the transient dying out into a frozen or periodic state. Note that the feature is labelled "activity rising" but its sign means activity is *falling*.
- **Fewer domains and a longer correlation length:** clusters of cooperators grow and merge (coarsening). A correlation length of 2 lattice spacings is short, not "long".
- **Cells flipping back and forth (g_reversals):** these are the standard period-2 "blinker" structures where defector remnants sit next to cooperators. They come from the periodic part of the theorem above, plus occasional flips from rule 1.
- **Rule 1 (majority vote with probability 0.003):** this is a weak, slow smoothing push. It slightly reinforces coarsening and sometimes perturbs the near-frozen state. Each rule alone gives similar numbers to the full program, so the two together are not producing anything new.

**Is it genuine or an artifact?** Genuine, but ordinary. It isn't diverged or random. However, the "interaction knocked out" baseline is weak. Both rules read neighbours, so removing them leaves the random starting pattern frozen in place: activity 0, correlation length 1. Beating that baseline only shows that the dynamics do something. The late state is a near-absorbing, mostly cooperative pattern with blinkers, which is exactly what the textbook predicts for b below about 1.75–1.8.

**Is there any twist?** None that matters. Spatial games with added noise or extra update rules are heavily studied (see the Szabó & Fáth review). No feature measured here (domain count, falling activity, short correlations, blinkers) goes beyond the documented low-b Nowak–May behaviour.

I ran one search, which turned up the original Nature paper, the Wu–Fu–Wang periodicity theorem and an implementation summary. I verified these through search snippets only and did not open the Nature paper itself.

---

## review:056

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour; (2) neighbouring cells swap contents (with probability 0.05 per step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_type_entropy

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.4773 with interaction, 6.3776 with interaction knocked out (difference -1.9003, robust spread 0.2171); difference at 4x area -1.8924, at 3x run length -1.8953
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.8746 with interaction, 0.0454 with interaction knocked out (difference 0.8292, robust spread 0.0184); difference at 4x area 0.8449, at 3x run length 0.813
  - g_type_entropy (cell types evenly mixed): median 0.5675 with interaction, 0.9999 with interaction knocked out (difference -0.4324, robust spread 0.0313); difference at 4x area -0.4297, at 3x run length -0.4257

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 4.7999; g_corrlen 2.0; g_reversals 1.0; g_largest_domain 0.8131; g_type_entropy 0.6954
  - rules 2 only: g_n_domains 6.35; g_corrlen 1.0; g_reversals 1.0; g_largest_domain 0.0381; g_type_entropy 0.9999
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial snowdrift (hawk-dove) game with Nowak-May best-takes-over imitation plus random mobility/diffusion

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Payoffs T=1.5>R=1>S=0.5>P=0 make this a snowdrift game, not a Prisoner's Dilemma. The observed outcome (one strategy at roughly 85% with the other persisting as a minority, short-range clusters) is the documented spatial snowdrift result. The swaps only modestly shift the strategy balance (largest domain 0.81 to 0.87), consistent with the mobility studies. g_reversals comes from the swaps alone and carries no information about the game. The knock-out baseline (an inert 50/50 mix) makes the significance tests easy to pass. The first pass's 'variant' label is too generous because no new behaviour was found; the variant looks closer to a reparameterisation of known work.

Citations:
- Hauert & Doebeli 2004, Nature 428:643 (doi:10.1038/nature02360)
- Sicardi, Fort, Vainstein & Arenzon 2009, J. Theor. Biol. 256:240 (arXiv:0809.1479)
- Vainstein, Silva & Arenzon 2007, J. Theor. Biol. 244:722

Full reply text:

**Verdict: KNOWN.** This overturns the first pass's "variant" call. I ran two searches; the results show the setup and its behaviour have both been studied for about 20 years.

**What the system actually is.** The payoffs give T (1.5) > R (1) > S (0.5) > P (0). That ordering is the **snowdrift game** (also called hawk–dove or chicken), not the Prisoner's Dilemma. The update rule (play all 8 neighbours plus yourself, then copy the best-scoring neighbour, on a 64x64 wrap-around lattice) is the classic Nowak–May lattice setup, here applied to the snowdrift game.

- **The game without movement (rule 1 only):** Hauert & Doebeli (Nature 2004) studied the spatial snowdrift game and found that in unstructured snowdrift games, intermediate levels of cooperation persist. Unexpectedly, spatial structure reduces the proportion of cooperators for a wide range of parameters. Their supplementary material also provides a detailed account of the method of pair approximation applied to game theory in spatially structured populations together with results for alternative update rules.
- **Adding random swaps (rule 2):** Sicardi, Fort, Vainstein & Arenzon (J. Theor. Biol. 2009) studied the effects of an unconditional move rule in the spatial Prisoner's Dilemma, Snowdrift and Stag Hunt games. They showed that random dilution and mobility may suppress the inhibiting factors of the spatial structure in the Snowdrift game. This builds on Vainstein, Silva & Arenzon, "Does mobility decrease cooperation?" J Theor Biol 244: 722–728. One difference: they used diluted lattices with moves into empty sites, whereas this system swaps neighbours on a fully occupied lattice. That's a change in the rules, not the behaviour.

**Matching the observed features to known behaviour:**
- **Largest domain ≈ 0.87 and type entropy ≈ 0.57:** one strategy holds roughly 85% of the lattice while the other survives as a minority. Stable coexistence at a frequency other than 50/50 is the standard result for the spatial snowdrift game. The data don't say which type is the majority.
- **Movement matters only a little:** the largest domain is 0.81 with the game alone and 0.87 with swaps added. The swaps nudge the balance between the two strategies rather than changing the kind of pattern, which fits the papers above.
- **Correlation length 2 and about 4–5 domains:** short-range clusters, typical of lattice snowdrift games. Snowdrift rewards doing the opposite of your neighbours, so it forms filament-like clusters rather than large compact blocks.
- **Reversals = 1:** this tells us nothing. The swaps-only version scores 1.0 too, so it comes from the movement, not the game.

**Real pattern or artifact?** It's real but generic: lasting coexistence of two strategies, with one in the majority. It isn't a frozen state, a diverged simulation or a random colouring, because the scores hold up at 4x area and 3x run length. However, the knock-out comparison is easy to pass. With no imitation, random swaps simply keep the starting 50/50 mix, so almost any working game rule would differ from it. I found nothing in the behaviour that the snowdrift-with-mobility literature doesn't already cover.

---

## review:066

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 7 neighbours already are t+1 (with probability 0.03 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.5, S=0.5) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_activity_trend, g_reversals, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.0859 with interaction, 6.3456 with interaction knocked out (difference -2.2597, robust spread 0.2086); difference at 4x area -2.075, at 3x run length -2.119
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_activity_trend (activity rising): median 1.5217 with interaction, 0.0 with interaction knocked out (difference 1.5217, robust spread 0.0624); difference at 4x area 1.6803, at 3x run length 0.8214
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_largest_domain (one cell-type fills most of the lattice): median 0.9237 with interaction, 0.0446 with interaction knocked out (difference 0.8792, robust spread 0.0194); difference at 4x area 0.8918, at 3x run length 0.8721

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 6.1014; g_corrlen 1.0; g_activity_trend -1.9136; g_reversals 0.0; g_largest_domain 0.0476
  - rules 2 only: g_n_domains 4.7706; g_corrlen 2.0; g_activity_trend -1.009; g_reversals 1.0; g_largest_domain 0.8198
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Spatial snowdrift (hawk-dove) game with Nowak-May synchronous imitate-the-best dynamics on a Moore lattice with self-interaction, plus weak threshold flip noise

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Payoffs T>R>S>P make this a snowdrift game, not a PD. Rule 2 alone already gives near-dominance (0.82), clustering and flip-back oscillations. Rule 1 (with 2 types it just means switch to the other type) is a 3%/step noise that removes isolated minority cells and keeps activity going; this explains the small extra gain in dominance and the rising activity. That activity signal halves at 3x run length, so it looks transient. Correlation length 1 -> 2 sites is short-range. 'Overturn' only means moving from KNOWN-VARIANT to KNOWN: no twist in the behaviour was found that would make it a variant.

Citations:
- Nowak & May 1992, Evolutionary games and spatial chaos, Nature 359:826
- Hauert & Doebeli 2004, Spatial structure often inhibits the evolution of cooperation in the snowdrift game, Nature 428:643
- Szabo & Fath 2007, Evolutionary games on graphs, arXiv cond-mat/0607344 (sec. 6.5)

Full reply text:

**Verdict: KNOWN.** This is a documented model with a small noise term added. I'd move the first pass from "variant" to plain "known": none of the five features is a new kind of behaviour.

**1. What rule 2 actually is.** With R=1, T=1.5, S=0.5, P=0 the payoff order is T>R>S>P. That is the **snowdrift (hawk–dove / chicken) game**, not a prisoner's dilemma, even though the rule calls the types cooperate and defect. (In the usual snowdrift notation this is b=1.5, c=1, a cost-to-benefit ratio of 0.5.) The update itself is the classic Nowak–May cellular automaton. Each cell plays the game with all his neighbors (once with each neighbor) and with himself (Moore neighborhood). Then every cell updates at once: it adopts the strategy of whichever cell in its neighbourhood, itself included, has obtained the greatest payoff. Because the update is deterministic and simultaneous, this evolutionary rule defines a cellular automaton. The well-known runs of that automaton show clusters of both types of players expanding, colliding and fragmenting.

The snowdrift version on lattices is also well studied. Hauert & Doebeli (Nature 2004) found that spatial structure reduces the proportion of cooperators for a wide range of parameters, and later work extends this across networks and update rules. My searches turned up the Nowak–May and Hauert–Doebeli papers plus follow-ups; I didn't find a paper on this exact rule pair, but that doesn't matter here.

**2. What rule 1 actually is.** With only two types, "become t+1" just means "switch to the other type". So rule 1 makes a cell flip, 3% of the time, when 7 or more of its 8 neighbours are the other type. It's slow noise that erases isolated minority cells. On its own it leaves the random start nearly frozen: largest domain about 0.05, correlation length 1, activity falling.

**3. How the features map onto known behaviour.**
- **One type fills most of the lattice (0.92), and fewer domains:** rule 2 alone already gives 0.82 and 4.8. Rule 1 then removes stray single cells, so the gain to 0.92 and 4.1 is the expected result of adding that noise.
- **Cells flip back and forth:** the deterministic Nowak–May automaton is known for blinkers and oscillating boundaries. Rule 2 alone already gives the same value (1.0).
- **Activity rising:** the only thing that needs rule 1. Occasional random flips keep re-starting imitation fronts that would otherwise freeze. That noise stops lattice games from freezing is standard, though I didn't pull a citation for this exact point. The signal also halves at 3× run length (1.52 → 0.82), which looks like a slow transient rather than a lasting trend.
- **Long correlation length:** going from 1 to 2 lattice sites is not long-range order. It's what you get from clustering in a near-uniform state.

**4. Genuine or artifact?** It's genuine game dynamics: the field isn't inert, diverged or randomly coloured, and the knock-out really does change the outcome. But a strict check should note:
- The "full program beats both sub-programs" result comes from simple additive effects: noise cleans up and keeps re-starting a nearly settled snowdrift automaton.
- The labels oversell two features: "long correlation length" is only 2 sites, and "activity rising" is fading.
- At these payoffs the state is close to one-type dominance, so it's partly near-absorbing, not a rich self-organised pattern.

---

## review:067

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) a type-1 cell spontaneously becomes type 0 (with probability 0.01 per step); (2) cells play a two-strategy game with their 8 neighbours and themselves (type 0 cooperates, 1 defects; payoffs R=1, P=0, T=1.3, S=0.0) and copy the best-scoring neighbour.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_n_types, g_reversals, g_activity_trend

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 3.1781 with interaction, 0.6931 with interaction knocked out (difference 2.4849, robust spread 0.0774); difference at 4x area 3.7492, at 3x run length 2.2743
  - g_corrlen (long spatial correlation length): median 2.0 with interaction, 0.0 with interaction knocked out (difference 2.0, robust spread 0.04); difference at 4x area 2.5, at 3x run length 2.0
  - g_n_types (number of surviving types): median 2.0 with interaction, 1.0 with interaction knocked out (difference 1.0, robust spread 0.06); difference at 4x area 1.0, at 3x run length 1.0
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_activity_trend (activity rising): median -0.798 with interaction, -3.9248 with interaction knocked out (difference 3.1268, robust spread 0.1199); difference at 4x area 3.3463, at 3x run length 3.1889

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 0.6931; g_corrlen 0.0; g_n_types 1.0; g_reversals 0.0; g_activity_trend -3.9248
  - rules 2 only: g_n_domains 4.1972; g_corrlen 2.0; g_n_types 2.0; g_reversals 1.0; g_activity_trend -1.4476
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Nowak-May spatial Prisoner's Dilemma (weak PD, imitate-the-best, Moore neighbourhood + self-interaction), low-b regime: cooperator-dominated coexistence with static defector clusters and blinkers; here with one-way D->C mutation noise

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Agrees with first pass on substance; relabels VARIANT->KNOWN because the only novelty is the D->C noise rule, which per the sub-program data only changes amounts (fewer domains: 3.18 vs 4.20) and leaves the corrlen/types/reversals phenomenology unchanged. Pattern is genuine but low-complexity (b=1.3 is well below the chaotic window ~1.8-2 for this neighbourhood).

Citations:
- Nowak & May 1992, Nature 359:826 (https://www.nature.com/articles/359826a0)
- Szabo & Fath 2007, Evolutionary games on graphs, Phys. Rep. (arXiv cond-mat/0607344)
- Stanford Encyclopedia of Philosophy, Evolutionary Game Theory (blinkers / two-cycle coexistence)
- Nowak, Bonhoeffer & May 1994, More Spatial Games, IJBC 4:33 (errors/noise variants; not directly verified here)

Full reply text:

**Verdict: KNOWN.** I agree with the first pass that this behaviour is documented. Where I differ is the label: I'd call it simply known rather than a variant, because I found nothing new in the behaviour itself. Only the noise rule is unusual.

**What the system is.** Rule 2 on its own is the Nowak–May spatial Prisoner's Dilemma. In that model each individual plays its immediate neighbours each round, and each site then goes either to its owner or to whichever neighbour scored highest. Your payoffs (R=1, P=0, S=0, T=b=1.3) match the standard simplified matrix, which has one free parameter b and is meant to keep the essence of the Prisoner's Dilemma. Your neighbourhood (8 neighbours plus self, synchronous updates) is the one they used: self, nearest and next-nearest neighbours.

**How each observed feature maps to the literature:**
- **Two types survive in many separate domains:** This model is known for spatial structure in which cooperators and defectors coexist, with patterns ranging from static fractal states to chaotic evolving structures depending on b. In particular, large defector clusters shrink when b < 1.8, while 2x2 or larger cooperator clusters keep growing when b < 2. At b=1.3 that gives a cooperator background holding small, stable defector clusters. That fits a correlation length of about 2 and many domains.
- **Cells flipping back and forth:** These are the known "blinkers." For a nearby weak-dilemma payoff set, the dynamics settle into a stable cycle between two states, with some regions containing "blinkers" that oscillate between defectors and cooperators.
- **Falling activity:** This is the expected settling into a steady state. After a transient, the system reaches a steady state whose nature depends strongly on b.

**What rule 1 adds.** The defector-to-cooperator flips are one-way noise. Nowak–May models with errors or mutation were studied soon afterwards (Nowak, Bonhoeffer & May 1994, "More Spatial Games"). I did not confirm that paper's contents in this pass. Your own sub-program results show the noise only changes amounts, not the kind of pattern. With rule 2 alone you get more domains (4.20 vs 3.18) and the same correlation length, number of types and reversals. Rule 1 alone gives a frozen all-cooperator state. So the noise just thins out defector clusters, and the defectors re-invade lone cooperators. I could find no change in the phenomenology to count as a twist.

**Real or artefact?** It's a real (if well-known) emergent pattern, not an artefact. With the interaction knocked out, you get the trivial absorbing state: all cells end up as one type. The pattern itself is mostly near-static clusters plus blinkers, though, not rich spatial chaos (which in this neighbourhood needs b of about 1.8–2).

---
