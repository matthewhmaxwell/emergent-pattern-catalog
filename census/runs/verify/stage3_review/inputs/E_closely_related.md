# Stage 3 (b) review inputs: E_closely_related (5 phenomena)

## review:046

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (every step); (2) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_reversals, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 3.9703 with interaction, 7.3558 with interaction knocked out (difference -3.3855, robust spread 0.2265); difference at 4x area -3.735, at 3x run length -3.4642
  - g_corrlen (long spatial correlation length): median 3.0 with interaction, 1.0 with interaction knocked out (difference 2.0, robust spread 0.08); difference at 4x area 2.0, at 3x run length 2.0
  - g_reversals (cells flip back and forth): median 0.9407 with interaction, 0.0 with interaction knocked out (difference 0.9407, robust spread 0.0188); difference at 4x area 0.9382, at 3x run length 0.8389
  - g_moran (like-valued cells cluster spatially): median 0.7195 with interaction, -0.0023 with interaction knocked out (difference 0.7218, robust spread 0.0164); difference at 4x area 0.7208, at 3x run length 0.7241

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 1.0397; g_corrlen 1.5; g_reversals 0.8591; g_moran 0.3564
  - rules 2 only: g_n_domains 7.3398; g_corrlen 1.0; g_reversals 0.0; g_moran 0.0092
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN-VARIANT, standard name: Synchronous majority-vote coarsening (Goles-Olivos period-2 blinking) perturbed by a threshold-1 cyclic cellular automaton (Griffeath CCA / rock-paper-scissors invasion at domain boundaries)

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Reversals (0.94) are mostly the known synchronous-majority period-2 artifact; majority alone gives 0.86. The knock-out removes both rules and leaves a frozen random start, so beating it is weak evidence. The real interaction effect is that cyclic invasion at domain edges roughly doubles clustering (Moran 0.36 -> 0.72) and raises the domain measure over majority alone. No spirals or travelling waves are reported; correlation length is only 3 cells. The exact majority + cyclic combination was not found in a paper, but the behaviour is a quantitative variant of known behaviours, not a new mechanism.

Citations:
- Goles & Olivos 1980, Periodic behaviour of generalized threshold functions, Discrete Math 30:187 (synchronous majority -> fixed points or 2-cycles)
- Poljak & Sura (multi-colour majority period <=2), as cited in Biased Majority Rule Cellular Automata
- Fisch, Gravner & Griffeath 1991, Cyclic Cellular Automata in Two Dimensions
- Shalizi et al., Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067 (low threshold -> global periodic oscillation, near-random texture)
- MCell Cyclic CA lexicon, psoup.math.wisc.edu (R1/T3/C3/NM is the only self-organizing 3-colour range-1 CCA)

Full reply text:

**Verdict: KNOWN-VARIANT, confidence medium-low.** I agree with the first pass, but for a different reason. The behaviour isn't new in any meaningful way. It is two well-documented behaviours added together, and the knock-out tests overstate how special it is.

**What the two rules are**
- **Rule 2** is the standard cyclic cellular automaton (Griffeath's model) with 3 types, threshold 1 and the 8-cell neighbourhood. In that model a cell of state i changes to state (i+1) mod k if at least τ neighbouring cells are already in state i+1. Its behaviour at this setting is documented: with a threshold this low, the entire lattice eventually oscillates periodically. With the 8-cell neighbourhood and threshold 1, every cell steps through the colours, with only slight departure from uniform randomness. For the 3-type, range-1 rules, the documented self-organizing one is threshold 3, not threshold 1 ("the only three-color rule in the range-1 CCA rulespace that exhibits significant self-organization"). This fits your rule-2-only run: clustering near 0 (0.009) and no back-and-forth flips.
- **Rule 1** is the synchronous majority rule. Goles and Olivos showed that every trajectory on a finite torus eventually reaches a fixed point or a cycle of length two. The extension to multiple colours is documented: a large class of majority-based models always reach a cycle of period one or two.

**What explains each feature**
- **Reversals (0.94):** mostly the known period-2 blinking of synchronous majority voting. Majority alone already gives 0.86. This feature comes from updating every cell at once; it isn't a new mechanism.
- **Clustering and domains (Moran 0.72, about 4 domains, correlation length 3):** majority voting coarsens the lattice into domains, and the cyclic rule only acts at domain edges, where a neighbour of the next type exists. This is the "competition between smoothing and cyclic invasion" picture. A correlation length of 3 cells on a 64×64 grid is short, so these are small domains, not large structures. The run doesn't show spirals; with majority voting applied first each step, spirals are presumably suppressed.
- **The knock-out comparison is weak evidence.** Both rules read neighbours, so the knock-out removes everything and leaves a frozen random start. Beating that baseline by 3 spreads proves almost nothing. The meaningful comparison is with the single-rule runs. There the only real effect is that the combination roughly doubles clustering (0.36 to 0.72) and raises the domain measure (1.04 to 3.97). Domain count and correlation length are also partly inherited from rule 2 alone (7.34 and 1.0).

**The twist I couldn't find described:** majority-vote domains with period-2 blinking cells that are sharpened, rather than destroyed, by cyclic invasion at their edges. This gives clustered domains that flicker back and forth, not travelling waves or spirals. I didn't find a paper on exactly this combination of majority and cyclic rules. However, I think it only plausibly matters as a quantitative change, not a new kind of behaviour.

**Genuine pattern or artifact:** the clustering is a real but modest effect of the rules interacting. The reversal feature is largely an artifact of updating all cells at once.

**Searches run:**
- Cyclic CA with 3 types, threshold 1, 8-cell neighbourhood (Fisch, Gravner, Griffeath).
- Majority rule combined with cyclic dominance / rock-paper-scissors.
- Synchronous majority 2-cycles (Goles–Olivos).
- Cyclic CA with majority smoothing and flickering domain edges.

---

## review:050

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_activity, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 6.4474 with interaction, 7.665 with interaction knocked out (difference -1.2176, robust spread 0.2822); difference at 4x area -1.3034, at 3x run length -1.3588
  - g_activity (many cells change per step): median 0.6131 with interaction, 0.0 with interaction knocked out (difference 0.6131, robust spread 0.0123); difference at 4x area 0.6208, at 3x run length 0.6229
  - g_moran (like-valued cells cluster spatially): median 0.3197 with interaction, 0.0098 with interaction knocked out (difference 0.3099, robust spread 0.0213); difference at 4x area 0.3401, at 3x run length 0.3122

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 7.553; g_activity 1.0; g_moran 0.0723
  - rules 2 only: g_n_domains 6.6114; g_activity 0.0; g_moran 0.23
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN-VARIANT, standard name: Cyclic cellular automaton (Griffeath; K=4, threshold 1, 8-neighbour) waves/spirals coupled to Schelling segregation moves; similar to cyclic dominance with mobility (Reichenbach-Mobilia-Frey)

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Rule 1 alone is the standard cyclic automaton: every cell changes every step, as in the known wave/demon-dominated state. Rule 2 alone is the standard Schelling model: it segregates and freezes. Combined, the moves partly break up the waves (activity about 0.61) and push up same-colour clustering (Moran about 0.32), the same kind of in-between regime as mobility breaking up spirals in rock-paper-scissors models. Twist not found in the literature: the empty state is part of the colour cycle, and moves are long-range jumps triggered by being unhappy with neighbours rather than local swaps. Likely changes the numbers, not the type of behaviour. The domain-count effect is small and the weakest of the three features.

Citations:
- Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata, Stat. Comput. 1 (1991) 23-39
- Durrett & Griffeath, Asymptotic behavior of excitable cellular automata, Exp. Math. 2 (1993)
- Reichenbach, Mobilia & Frey, Nature 448 (2007) 1046-1049
- Szczesny, Mobilia & Rucklidge, The Influence of Mobility Rate on Spiral Waves in Spatial Rock-Paper-Scissors Games, Games 7 (2016) 24
- Schelling (1971) segregation model with random relocation

Full reply text:

**Verdict: KNOWN-VARIANT.** I agree with the first pass, and I mean that as a real check rather than a default. The behaviour is a cyclic cellular automaton (cells cycling through colours, producing waves) partly broken up by Schelling-style moving. Both parts are well described in the literature. The only new thing is how the two mix, and I don't think that twist matters much.

**Rule 1 alone is a standard cyclic cellular automaton.** Griffeath's cyclic automaton is defined so that a cell advances to (s+1) mod K if at least θ of its neighbours already hold that successor state. Your rule 1 is exactly the K=4 colours, threshold 1, 8-neighbour version. Threshold 1 is the most excitable, giving fewer but larger spirals, and starting from random noise these systems show nucleation of locally periodic spirals and a variety of equilibria that display large-scale stochastic wave fronts. Your rule-1-only run has every cell changing every step (activity 1.0) and almost no same-colour clustering (Moran's I = 0.07, where 0 means no spatial correlation). That fits the known end state where waves and self-sustaining wave sources ("demons") take over the whole grid: every cell always has a successor-coloured neighbour, and the wave bands are only a few cells wide. Main sources: Fisch, Gravner & Griffeath (1991); Durrett & Griffeath (1993).

**Rule 2 alone is the Schelling model with random relocation.** Agents stay put when this fraction is above a given tolerance threshold but seek a new location if the fraction is below the threshold, and the new location is selected randomly from the set of unoccupied locations with uniform probability. Your rule-2-only run segregates into clusters (Moran 0.23) and then freezes (activity 0). That is the textbook outcome.

**The combination is cyclic dominance plus mobility, which is also documented.** Reichenbach, Mobilia & Frey (Nature 2007) found that when mobility exceeds a certain value, biodiversity is jeopardized and lost. In contrast, below this critical threshold all subpopulations coexist and an entanglement of travelling spiral waves forms. Later work showed that mobility can make spiraling patterns become convectively unstable and that a far-field breakup occurs. Your full system sits in that kind of in-between regime:
- Activity drops from 1.0 to about 0.61, meaning the waves are partly broken up.
- Clustering rises from 0.07 to about 0.32, because moving builds same-colour patches that the cycling then erodes.
- The number of domains goes down slightly.

That is what you would expect when cycling waves and segregation-by-moving compete.

**The twist I couldn't find described:**
1. The "empty" state (type 0) is itself one of the four colours in the cycle. So the spaces Schelling agents move into are themselves cycling, and they get filled and emptied by the waves, not only by moves.
2. Movement depends on how happy a cell is with its neighbours and sends it anywhere on the grid. Standard rock-paper-scissors mobility models instead use random swaps between neighbours.

I didn't find this exact coupling. But it probably only shifts numbers (how much activity drops, how much clustering rises), not the kind of behaviour. My searches were: cyclic automata with spirals and threshold scaling; rock-paper-scissors with mobility and spiral breakup; and Schelling combined with cyclic dominance. The last one turned up nothing directly on point.

**It looks like a genuine emergent pattern, not an artifact.** The activity is high and steady across seeds, larger grids and longer runs, so the system isn't frozen. Clustering is well above both the knocked-out version and either rule alone, so it isn't a random colouring. One caution: the domain-count difference is small (−1.2 domains) and goes in the same direction as rule 2 alone, so it is the weakest of the three features.

---

## review:051

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell; (2) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_activity, g_interface

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 6.2442 with interaction, 7.6693 with interaction knocked out (difference -1.4251, robust spread 0.2783); difference at 4x area -1.4038, at 3x run length -1.316
  - g_activity (many cells change per step): median 0.6728 with interaction, 0.0 with interaction knocked out (difference 0.6728, robust spread 0.0135); difference at 4x area 0.6756, at 3x run length 0.6735
  - g_interface (long/rough domain boundaries): median 0.3755 with interaction, 0.7498 with interaction knocked out (difference -0.3743, robust spread 0.0225); difference at 4x area -0.3757, at 3x run length -0.3808

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 6.0958; g_activity 0.0; g_interface 0.3552
  - rules 2 only: g_n_domains 7.5535; g_activity 1.0; g_interface 0.7156
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN-VARIANT, standard name: Schelling segregation coupled to a Griffeath cyclic cellular automaton (κ=4, T=1, Moore range 1; locally periodic 'quench' regime) — cf. mobility in cyclic-dominance/RPS lattice models and Schelling variants with switching agents

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Rule 2 alone reproduces the documented CCA T=1 Moore 4-colour locally periodic state (activity=1.0). The full program's domain/interface values are close to the Schelling-only sub-program; activity is the CCA cycling partly suppressed (~0.67), probably because relocations and the 3->empty step insert defects. The knock-out baseline is a frozen random state, so beating it is weak evidence; the sub-program comparison suggests superposition with mild interference. Twist not found in the literature: persistent cycling of Schelling-coarsened domains at intermediate activity. Probably not significant unless images or spiral/phase-defect statistics show new structure.

Citations:
- https://arxiv.org/html/nlin/0507067
- http://psoup.math.wisc.edu/mcell/rullex_cycl.html
- https://arxiv.org/pdf/1807.05419
- https://arxiv.org/pdf/1212.4929
- https://www.sciencedirect.com/science/article/abs/pii/S0378437117306118
- https://arxiv.org/pdf/1608.08932

Full reply text:

**My verdict is KNOWN-VARIANT, with medium confidence.** The behaviour comes from two well-studied models run together, and the measurements look roughly like the two added up. I agree with the first pass, but for weaker reasons than it gave.

**What the system is made of**
- **Rule 2 is Griffeath's cyclic cellular automaton (CCA) at its most excitable setting:** 4 colours, threshold 1, 8-neighbour (Moore) range. In the CCA, each color can only advance to the next, the last cycling to 0, and a cell advances by 1 if at least Threshold cells of the next color are in its neighbourhood. The literature describes exactly this setting. With a Moore range-1 neighbourhood, T=1 and T=2 are both locally periodic, but T=2 produces spiral waves, while T=1 "quenches" incoherent local oscillations. The reference figures use the same setup as yours: 4 colors, a range-one box-shaped neighbourhood, and a uniform random start. Your "rule 2 only" result (activity exactly 1.0, every cell cycling every step) matches this documented locally periodic behaviour. That is a good sign the simulator is correct.
- **Rule 1 is Schelling segregation** with threshold 0.5 and moves to random empty cells. Schelling exhibits large-scale self-organizing neighborhoods, due to agents' desire to live close to their own kind.
- **Mixing the two kinds of process is also documented:**
  - Schelling with agents that change type: a fraction of agents able to switch type spontaneously, which can produce intermediate patterns (e.g. mixture of patches, fuzzy interfaces).
  - The Schelling–Voter model: an unhappy agent either migrates or changes her own type.
  - Cyclic dominance with mobility: spiral waves appear in lattice simulations only within a finite range of the mobility rate.

**Reading the numbers adversarially**
- **The domain and interface features belong to Schelling.** The full program (6.24 domains, interface 0.376) is close to Schelling alone (6.10, 0.355) and far from CCA alone (7.55, 0.716). It is not within 3 robust spreads on interface, but it is clearly Schelling-like.
- **The activity feature belongs to the CCA,** held down from 1.0 to about 0.67. Most likely the random relocations, and the cyclic step that turns type 3 into an empty cell, insert defects that break up the synchronised oscillation.
- **The knock-out test is weak here.** Both rules read neighbours, so removing interactions leaves a frozen random state. Beating that baseline by 3 spreads shows almost nothing. The sub-program comparison is the meaningful test, and it suggests superposition with partial interference rather than a new pattern.

**The twist I couldn't find described:** a persistent state where Schelling-coarsened domains keep changing colour by excitable cycling, with about two-thirds of cells changing per step. This is in between the documented CCA "quench" state (activity 1) and Schelling's frozen state (activity 0). I didn't find this exact Schelling + CCA mix described. It probably doesn't matter much: it's what you would expect from mobility acting as noise on a locally periodic CCA. The measured features don't show any new structure, such as spirals, travelling fronts or a change in the scaling of domain size. Turning it into a finding would need images, autocorrelation or spiral-tip counts.

**Genuine or artifact?** Not trivial in the narrow sense: the state isn't frozen, absorbing or numerically diverged, and it stays stable at 4× area and 3× run length. But it is only weakly emergent, mostly the sum of the two sub-programs.

**Searches run:** CCA threshold-1 Moore spirals (Fisch, Gravner, Griffeath); Schelling with cyclic dominance / rock-paper-scissors mobility spirals; Schelling combined with excitable dynamics or type-changing agents.

---

## census:novelty:network-2

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a random network of 200 nodes (average 4 links each). Nodes have 2 types (random start). Rules, applied in order each step: (1) nodes play a two-strategy game with their neighbours (T=1.9, S=0.0) and copy the best-scoring neighbour; (2) each node adopts its neighbours' majority type (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  n_activity, n_type_modularity

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - n_activity: median 0.2895 with interaction, 0.0 with interaction knocked out (difference 0.2895, robust spread 0.0058); difference at 4x area 0.3294, at 3x run length 0.3777
  - n_type_modularity (network splits into like-type modules): median 0.167 with interaction, -0.0031 with interaction knocked out (difference 0.17, robust spread 0.0137); difference at 4x area 0.1657, at 3x run length 0.1415

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: n_activity 0.0233; n_type_modularity 0.1579
  - rules 2 only: n_activity 0.0; n_type_modularity 0.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN-VARIANT, standard name: Nowak-May spatial (weak) prisoner's dilemma with imitate-the-best on a random graph, plus conformity / majority-rule label dynamics (tag hitchhiking)

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Type modularity is not an interaction effect: rule-1-only gives 0.158 vs 0.167 for the full program, within one robust spread, which implies best-neighbour copying carries the type along (tag hitchhiking). Only the activity (0.29 vs 0.02/0.0) needs both rules. It is plausibly persistent frustration between payoff-driven copying and majority-rule relaxation of the label, which is unremarkable. Twist not found in print: conformity acting on a label separate from the strategy rather than on the strategy. Probably not important. Check in the code whether imitation copies the type.

Citations:
- Szabo & Fath, Evolutionary games on graphs, arXiv cond-mat/0607344 (Nowak-May weak PD, imitate best)
- Pacheco et al., Conformity and Network Effects in the Prisoner's Dilemma
- Evolutionary PD on complex networks with conformist mentality strategy, Sci China PMA 2012
- Kim, A Tag-Based Evolutionary Prisoner's Dilemma Game on Networks with Different Topologies, JASSS 13(3)2
- Jensen & Bornholdt, Imitating the winner leads to discrimination in spatial PD, Sci Rep 9:3776 (2019)

Full reply text:

**Verdict: KNOWN-VARIANT, medium confidence. I agree with the first pass, but for different reasons.** This is the Nowak–May "imitate the best neighbour" prisoner's dilemma on a random graph, with a conformity (majority-rule) process added to a second label. Both observed features follow from well-documented mechanisms. The one thing I couldn't find described is churn produced specifically by the label clash (see the twist below). I ran 4 searches.

**What each piece already matches**
- **The game rule.** It is the Nowak–May model exactly. They used the "weak" payoffs (temptation b, sucker payoff 0), and each player imitates the strategy of those neighboring players (including herself) who has scored the highest payoff. The case c=0 is sometimes called a "Weak" Prisoner's Dilemma, and it behaves like the standard version. The rule has been studied on random graphs; for example, one paper ran binomial random-graphs with 1000 nodes and average connectivity 4. On such graphs it usually settles into a near-frozen state. That fits your rule-1-only activity of 0.023.
- **Payoff imitation combined with majority-following.** This is an established topic. Pacheco and colleagues confirm previous results on the improved levels of cooperation when both network reciprocity and conformity are present. Other models update strategies using not only the payoff differences between players and the neighbors, but also the majority of neighbors.
- **Labels copied from successful neighbours.** Tag models on networks cover this: agents imitate tags and learn tolerance from more successful neighbors. They also report endless cycles between cooperative societies and betrayal societies in small populations. Jensen and Bornholdt (Sci Rep 2019) used a random binary label with best-neighbour imitation on the same kind of random graph, though the labels described in this model never change.

**The two features, checked against your evidence**
- **Type modularity is not caused by the two rules interacting.** Running rule 1 alone gives 0.158, against 0.167 for the full program. That gap (0.009) is under one robust spread (0.0137).
  - Rule 1 alone can only cluster the types if copying the best-scoring neighbour copies the whole node, type included. In other words, types ride along with successful strategies, as in the tag models above. If types were truly separate from the game, rule 1 alone should give about 0.
  - **You should check this in the code.** If the copy includes type, the "like-type modules" are simply payoff-imitation clusters, which is known behaviour.
- **Activity is the only effect that needs both rules.** It rises from 0.02 or 0.00 with either rule alone to 0.29 with both, and it holds at 4x area and 3x run length. The obvious mechanism is frustration: payoff imitation keeps rewriting types, and the 10% majority step keeps pulling them back. Under the conformity framing, ongoing churn is expected. The paper on conformity plus Q-learning, for instance, notes that in some settings the overall population becomes highly unstable.

**The twist**
In the papers I found, conformity acts on the strategy itself. Here it acts on a separate label that may only be linked to the game by being copied along with it. I couldn't find a paper documenting lasting high activity that comes specifically from a label caught between payoff copying and majority rule. It probably doesn't matter much: it is an unsurprising result of two update rules pulling against each other, not a new pattern.

**Genuine or artifact?**
It is a genuine interaction effect, not a frozen or diverged artifact: it is robust, it scales, and the knock-out gives zero. But it is mundane. The modularity feature on its own is close to trivial, because it is almost entirely explained by rule 1. Rule 2 alone gives zero for both features, which most likely means the majority rule reaches consensus or freezes on 200 nodes.

---

## qd:037

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 3 types (random start). Agents have 4 types (random start). Agents move at speed 0.0 with heading noise 1.0 rad per step. Rules, applied in order each step: (1) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.2 moves to a random empty cell; (2) each agent advances the type of the cell under it by one (cyclically).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_type_entropy, g_reversals, g_interface

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 6.0811 with interaction, 7.3664 with interaction knocked out (difference -1.2853, robust spread 0.2689); difference at 4x area -1.2924, at 3x run length -1.2876
  - g_type_entropy (cell types evenly mixed): median 0.7364 with interaction, 0.9998 with interaction knocked out (difference -0.2634, robust spread 0.0347); difference at 4x area -0.263, at 3x run length -0.2624
  - g_reversals (cells flip back and forth): median 0.5455 with interaction, 0.0 with interaction knocked out (difference 0.5455, robust spread 0.0109); difference at 4x area 0.5525, at 3x run length 0.5758
  - g_interface (long/rough domain boundaries): median 0.3748 with interaction, 0.6701 with interaction knocked out (difference -0.2953, robust spread 0.0209); difference at 4x area -0.2916, at 3x run length -0.2982

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 7.1778; g_type_entropy 0.9998; g_reversals 0.0; g_interface 0.6271
  - rules 2 only: g_n_domains 7.3572; g_type_entropy 0.9999; g_reversals 0.0; g_interface 0.6681
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN-VARIANT, standard name: Schelling segregation with spatially fixed type-switching sites (Hazan & Randon-Furling 2012 'switching agents'); noise-driven unfreezing of a low-threshold Schelling lattice

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Speed 0 makes the agents fixed sites that cycle cell type each step, which is equivalent to Hazan & Randon-Furling's fixed switching agents in a Schelling lattice. The domain and interface reductions are ordinary Schelling clustering, unlocked because the switching sites keep making cells unhappy (at threshold 0.2 with ~1/3 empties the knock-out barely moves). Twist: the coupling does not conserve counts. When an unhappy switching-site cell relocates, the agent refills the vacated site, steadily turning empty cells into occupied ones. This vacancy pumping most likely explains the entropy drop (a change in type counts, not spatial order). Non-zero reversals are built into the rules (fixed sites cycling plus relocations). Two of the four features are therefore partly artifacts of the coupling; the clustering features look genuine.

Citations:
- https://arxiv.org/pdf/1212.4929
- https://arxiv.org/pdf/0903.4694
- https://arxiv.org/pdf/2204.01417
- https://arxiv.org/pdf/1504.03809

Full reply text:

**Verdict: KNOWN-VARIANT. This agrees with the first pass, for a more specific reason.**

**1. What the system actually is.** The agents have speed 0, so heading noise does nothing. The 400 agents are effectively fixed lattice sites, each pushing the cell under it through 0 → 1 → 2 → 0 every step. That makes the system Schelling relocation dynamics plus a fixed scattered set of sites that keep changing type. That exact setup has been published. Hazan & Randon-Furling (2012) add a given fraction f of fixed agents able to switch, that is: to change spontaneously from one of the two types of agents to the other. They describe this as imposing a noise or a perturbation in the form of a spatially-fixed random background within the Schelling dynamics, and they find it can lead to the emergence of intermediate patterns (e.g. mixture of patches, fuzzy interfaces).

**2. Why the knock-out and the "Schelling only" run look so mixed.** The baseline is known. With a random start, a 0.2 threshold (8 neighbours) and about one-third empty cells, almost every cell is already satisfied, so Schelling alone hardly moves (interface 0.627 vs 0.670). In the full system, the switching sites keep producing unhappy cells, so relocation and clustering continue. Fewer domains and shorter boundaries (interface 0.37) are ordinary Schelling cluster formation. It is known that the vacancy density and tolerance set whether a Schelling lattice segregates or stays mixed: For each density of vacancies there exist two critical values of the tolerance, and For low tolerance values, the system freezes in a disordered state. That last case is the inverse convention of the 0.2 threshold here. Adding noise to Schelling to change its outcome is also an established idea (Young... made another highly influential innovation by introducing noise to the model).

**3. The twist, and whether it matters.** Unlike Hazan's switchers, which conserve counts, this coupling creates occupied cells. When a switching site is type 1 or 2 and unhappy, it relocates to a random empty cell. That leaves the site empty, and the agent then turns it into type 1. Each such event turns an empty cell into an occupied one, so the system pumps away its vacancies. This most likely explains the drop in type entropy (0.9998 → 0.736): it is a shift in the type counts (fewer empties, uneven types), not spatial organisation. The non-zero reversal count is likewise built in: switching sites plus relocations away from and onto them. Sub-programs can't show it because their reversals are zero by construction. I didn't find this "vacancy pumping" effect described anywhere. It matters only as a caveat: two of the four features are mostly bookkeeping effects of the coupling, not emergent patterns.

**4. Real pattern or artifact?** It's mixed:
- **Likely real:** the domain count and interface reduction (Schelling clustering unfrozen by the switching sites).
- **Mostly built in by the rules:** entropy and reversals.
- **Not seen:** any numerical divergence or fully frozen state. Strictly, I can't rule out that the late-run state still has some frozen character.

---
