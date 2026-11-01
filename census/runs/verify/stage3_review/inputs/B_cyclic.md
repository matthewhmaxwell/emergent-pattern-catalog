# Stage 3 (b) review inputs: B_cyclic (12 phenomena)

## census:g1:grid-25

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) neighbouring cells swap contents (every step); (2) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_reversals, g_activity_cv

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 6.5043 with interaction, 7.6688 with interaction knocked out (difference -1.1645, robust spread 0.2835); difference at 4x area -1.1769, at 3x run length -1.3032
  - g_reversals (cells flip back and forth): median 0.7457 with interaction, 0.0 with interaction knocked out (difference 0.7457, robust spread 0.0149); difference at 4x area 0.7084, at 3x run length 0.8615
  - g_activity_cv (bursty activity): median 0.4376 with interaction, 0.0 with interaction knocked out (difference 0.4376, robust spread 0.0088); difference at 4x area 0.5215, at 3x run length 0.4455

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 7.6641; g_reversals 0.3962; g_activity_cv 0.519
  - rules 2 only: g_n_domains 4.3757; g_reversals 0.0; g_activity_cv 0.1396
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Stochastic cyclic cellular automaton (Griffeath CCA, range-1 Moore, threshold 3, 4 colours, 'turbulent' regime) with stirring/mobility (cyclic dominance with exchange; Durrett-Neuhauser fast stirring; Reichenbach-Mobilia-Frey mobility vs pattern formation)

Trivial or artifact: True; confidence: medium

Twist (if any): 

Notes: Upgrades the first pass from VARIANT to plain KNOWN: no phenomenological twist found. The knock-out also removes swapping, so it is a frozen lattice and the reversal/activity-CV contrasts are inflated. Swap-only already gives reversals 0.40 and activity CV 0.52 (higher than the full system's 0.44), so these features are mostly stirring noise. The only real interaction effect is that stirring partly suppresses CCA domain growth (domain count 6.5, between 4.4 and 7.7), which is the documented mobility-vs-pattern trade-off. None of the features directly tests for spirals or waves.

Citations:
- Fisch, Gravner, Griffeath (1991) Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39
- Shalizi, Shalizi, Haslinger (2004) Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067 (T=3, 4-colour, Moore: metastable turbulence)
- Durrett (1999) Stochastic spatial models, SIAM Review (rapid stirring limits)
- Reichenbach, Mobilia, Frey (2007) Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, Nature 448:1046

Full reply text:

**Verdict: KNOWN.** The behaviour itself is documented. It also looks mostly trivial: the "reversals" and "bursty activity" signals come largely from the swapping rule acting alone. I'm moving the first pass from "variant" to plain KNOWN because I found no new twist in the behaviour.

**What the system is.** Rule 2 alone is Griffeath's cyclic cellular automaton (CCA) with 4 colours, an 8-neighbour range-1 box neighbourhood and threshold 3, updated at random with probability 0.1. Rule 1 adds neighbour swaps, called "stirring" or mobility. Both parts are standard:

- **The cyclic automaton.** In a CCA, colors advance (k to k+1 mod κ) either automatically or by contact with at least a threshold number of successor colors in a prescribed local neighborhood. Your exact setting (4 colours, box neighbourhood, threshold 3) is a benchmark case. For range-one box neighbourhoods, at T=3, one encounters "turbulence," which is actually meta-stable — spiral waves can form and entrain the entire CA, but this does not always happen on finite lattices, and in any case the turbulent phase can persist for very long times. Fixation occurs with T≥4. That paper's figures use exactly 4 colors, a range-one box-shaped neighborhood, and the initial conditions are a uniform random mixture of colors.
- **Adding swaps.** Durrett and Neuhauser formalised swaps as stirring events and showed that with fast stirring the particle system on the rescaled lattice becomes deterministic and looks like solutions of the PDE.
- **Swaps plus cyclic dominance.** Reichenbach, Mobilia and Frey (Nature 2007) found a critical threshold of motility is identified, above which biodiversity is lost, but below which an entanglement of travelling spiral waves emerges. They also state this phenomenon is robust, it does not depend on the details of cyclic competition or spatial environment.

**Reading the evidence, feature by feature:**

- **g_n_domains.** Swaps alone give 7.66 (a random-like mixture) and the cyclic rule alone gives 4.38 (coarser domains). The full system sits in between at 6.50. This is the expected result when mixing works against domain growth, as in the mobility papers above. Nothing new here.
- **g_reversals.** The knock-out score of 0 is misleading. The knock-out removes swapping too, so it is a frozen lattice by construction. Swaps alone already give 0.40, and stirring by definition makes a site's colour go back and forth. The full system's 0.75 adds the cyclic flips, plus some swap-driven back-and-forth flipping near domain boundaries. That is expected, not new.
- **g_activity_cv.** Again the knock-out is frozen. Swaps alone give 0.52, which is *higher* than the full system's 0.44. So the cyclic rule makes activity *less* bursty, not more.

**Genuine pattern or artifact?** The three features passed the checks mostly because they are compared against a frozen baseline. Two of them are dominated by trivial stirring noise. The one real interaction effect is that swapping partly suppresses domain growth, and that is textbook. None of the three features directly measures spirals or travelling waves. It would take snapshots or a spatial correlation measure to tell whether turbulent wave fragments survive at this swap rate.

**Searches run:** cyclic CA with threshold and stirring (Durrett–Neuhauser); Fisch–Gravner–Griffeath threshold-range work; Reichenbach–Mobilia–Frey mobility in rock-paper-scissors.

---

## review:008

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  g_corrlen, g_n_domains, g_moran, g_cyclic_changes

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_corrlen (long spatial correlation length): median 5.0 with interaction, 1.0 with interaction knocked out (difference 4.0, robust spread 0.12); difference at 4x area 3.5, at 3x run length 3.0
  - g_n_domains (many separate domains): median 4.1023 with interaction, 7.3547 with interaction knocked out (difference -3.2523, robust spread 0.2291); difference at 4x area -3.3784, at 3x run length -3.228
  - g_moran (like-valued cells cluster spatially): median 0.7923 with interaction, 0.0024 with interaction knocked out (difference 0.7899, robust spread 0.0159); difference at 4x area 0.7744, at 3x run length 0.7787
  - g_cyclic_changes (cells cycle through states above chance): median 0.4877 with interaction, -0.5 with interaction knocked out (difference 0.9877, robust spread 0.0198); difference at 4x area 0.9878, at 3x run length 0.9801

### Stage 3 reply: verdict KNOWN, standard name: Cyclic cellular automaton (Fisch-Gravner-Griffeath), rule '313' (range 1, threshold 3, 3 colours, Moore), stochastic/asynchronous update

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Rule is exactly Griffeath's named CCA rule 313 (R1/T3/C3/NM), documented as the only 3-colour range-1 CCA rule with significant self-organization (spirals). Only difference is the 0.1 flip probability; random versions of CCA/GH are discussed by FGG and Griffeath (stochastic wave fronts). Clustering, domain coarsening, finite correlation length and above-chance cyclic changes are the standard CCA signatures. The knock-out is trivially frozen because there is no non-interacting rule, but the interacting system keeps cycling and the features are robust to size and run length. The overturn is mild: the first pass said 'variant', but the behaviour is not a variant in phenomenology. On 64x64 the data can't tell spirals from turbulent waves; both are documented.

Citations:
- http://psoup.math.wisc.edu/mcell/rullex_cycl.html
- https://link.springer.com/article/10.1007/BF01890834
- https://arxiv.org/pdf/nlin/0507067
- https://www.stat.cmu.edu/~cshalizi/462/readings/Self-Organization%20of%20Random%20Cellular%20Automata-%20Four%20Snapshots%20(Griffeath).pdf
- https://link.springer.com/chapter/10.1007/978-1-4612-0451-0_8

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass, which called this a variant. The rule is not a modification of a known model. It is a specific, named member of the Cyclic Cellular Automaton (CCA) family, with only the update made random.

**The rule is a documented rule.** In the standard CCA, a cell of color k will change its color to k+1 mod κ if there are already at least T cells of that color in its neighborhood. Your settings are 8 neighbours at distance 1, threshold 3 and 3 types. That is exactly Griffeath's rule "313". MCell's rule list gives it as 313 · R1/T3/C3/NM · Cyclic · "This is the only three-color rule in the range-1 CCA rulespace that exhibits significant self-organization. The same entry adds that it sports a number of distinct spiral types. The general theory is in Fisch, Gravner & Griffeath, "Cyclic Cellular Automata in Two Dimensions" (1991) and "Threshold-range scaling of excitable cellular automata."

**The randomness is not new either.** Your one difference from the classic rule is that a qualifying cell flips only with probability 0.1 per step, instead of every cell updating together. Fisch, Gravner & Griffeath already discuss random versions: started from random configurations the same rules exhibit complex self-organization, typically characterized by nucleation of locally periodic 'ram's horns' or spirals. Corresponding random processes give rise to a variety of 'forest fire' equilibria that display large-scale stochastic wave fronts. Griffeath also has a survey on randomly updated versions ("Self-organization of random cellular automata: four snapshots").

**The measurements match the known behaviour.** Your four features are what you'd expect from this documented self-organisation:
- like-coloured cells cluster (Moran's I about 0.79),
- there are fewer, larger domains,
- the correlation length is finite and about 5 cells,
- cells keep cycling through the three types more than chance would give.

These are the usual signs of travelling cyclic waves and spirals. I found nothing in the measurements that falls outside the documented behaviour.

**Real pattern, not an artifact.** The knock-out baseline is trivial: this system has no rule that doesn't read the neighbours, so removing interaction freezes every cell. That explains the −0.5 cycling score and the Moran's I of about 0. The interacting system itself is not frozen, though. It keeps cycling, and the features hold up at 4× area and 3× run length.

**Caveats:**
- I did not find a source that specifically describes 313 with flip probability 0.1. Random updating could make spirals look rounder or noisier, or break them up.
- A 64×64 grid with correlation length about 5 is small. Your measurements can't tell apart fully formed spirals, a turbulent wave state, or patchy oscillation. Even so, all three are documented CCA regimes.

---

## review:012

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (with probability 0.3 per step).

OBSERVED (features that passed the simulation checks):
  g_interface, g_activity

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_interface (long/rough domain boundaries): median 0.6521 with interaction, 0.7493 with interaction knocked out (difference -0.0972, robust spread 0.028); difference at 4x area -0.0982, at 3x run length -0.0961
  - g_activity (many cells change per step): median 0.4174 with interaction, 0.0 with interaction knocked out (difference 0.4174, robust spread 0.0083); difference at 4x area 0.4177, at 3x run length 0.4177

### Stage 3 reply: verdict KNOWN, standard name: Cyclic cellular automaton (Griffeath CCA, 4 colours, 8-cell neighbourhood, threshold 1) with random update, i.e. a discrete-time form of the Bramson-Griffeath cyclic particle system

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Shalizi & Shalizi describe exactly these settings (4 colours, 8-cell neighbourhood, threshold 1): every cell keeps cycling, and the field shows only a slight departure from randomness as twisting single-colour bands. This matches both measured features (interface 0.65 vs 0.75 random; sustained activity). Bramson-Griffeath / Fisch results say the 4-colour case keeps fluctuating rather than freezing. The p=0.3 random update is a rule detail, not new behaviour. The pattern is weak; the activity knock-out is trivially frozen; an activity of 0.417 above p=0.3 points to one-cell-at-a-time updating or a multi-step measure and should be checked.

Citations:
- https://arxiv.org/pdf/nlin/0507067
- https://arxiv.org/html/1711.04741
- https://en.wikipedia.org/wiki/Cyclic_cellular_automaton
- https://www.stat.cmu.edu/~cshalizi/462/readings/Cyclic%20Cellular%20Automata%20in%20Two%20Dimensions%20(Fisch,%20Gravner%20and%20Griffeath).pdf
- https://commons.wikimedia.org/wiki/Category:Cyclic_cellular_automata

Full reply text:

**Verdict: KNOWN.** I'd move the first pass's "variant" call to plain "known". This is Griffeath's cyclic cellular automaton (CCA), also studied as the cyclic particle system, run with 4 colours, an 8-cell neighbourhood and threshold 1. The only change is that cells update at random (probability 0.3 per step) instead of all at once. Both behaviours you observed have been described for this exact setting.

**What the rule is.** The standard CCA rule is that each step a cell advances by 1 if there are at least Threshold cells of the next color within its neighbour set. Your rule is that rule with threshold 1, 4 states and the 8-cell neighbourhood, plus the random update. Shalizi & Shalizi (arXiv nlin/0507067) used exactly these settings: the number of colors is fixed to 4, the rule uses a range-one box-shaped neighborhood, and the initial conditions are a uniform random mixture of colors.

**Your two features match what that paper describes:**
- **Long, rough boundaries (g_interface):** they report a slight departure from uniform randomness, in the shape of connected twisting single-color bands. That fits your result: 0.65 of neighbour pairs differ, against 0.75 for a random 4-colour pattern, so there is weak clustering but no large domains.
- **High activity (g_activity):** they note that the CA as a whole oscillates with period 4, each cell cycling through the colors.

**The random-update version is also covered.** Bramson & Griffeath (1989) introduced the stochastic counterpart. In one dimension they showed that the process fluctuates if κ≤4 and fixates otherwise. For the deterministic automaton, Fisch showed it fluctuates for κ∈{3,4}, clusters for κ=3, and fixates for κ≥5. With 4 colours you are on the side where cells keep changing for ever, as observed. More generally, for this family in 2D, patches of color spread and spirals form. Slowing the updates with p=0.3 only changes the timing and blurs perfect synchrony. I found no paper on this exact p=0.3 asynchronous version, but that is a rule-tuple gap, not missing behaviour.

**Genuine or artifact?** It is a real dynamical effect, but a weak pattern. The activity is sustained cycling, not divergence or a frozen state. Two cautions:
- The knock-out removes the only rule, so it can never change and activity is trivially 0. The activity "difference" really just says the dynamics exist at all.
- An activity of 0.417 is more than the 0.3 maximum you'd expect if each cell flips at most once per step. Either the rules are applied one cell at a time within a step, or "activity" is measured over more than one step. That's worth checking in the code, though it doesn't change the verdict.

---

## review:015

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) a node of type t becomes t+1 (cyclically) if at least 2 neighbours are t+1 (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  n_cyclic_changes

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - n_cyclic_changes (node states cycle): median 0.443 with interaction, -0.5 with interaction knocked out (difference 0.943, robust spread 0.0189); difference at 4x area 0.9396, at 3x run length 0.8619

### Stage 3 reply: verdict KNOWN, standard name: Threshold cyclic cellular automaton (Griffeath CCA / Greenberg-Hastings-type excitable medium), asynchronous update on a random graph

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Overturns 'variant' toward plain KNOWN: a random graph plus asynchronous updates adds no new phenomenology in the measured feature. The 'cyclic changes' metric is cyclic by construction (t->t+1 is the only allowed move), so it only measures persistent activity versus fixation. The knock-out value of -0.5 looks like a placeholder for zero activity, so the contrast is active versus completely frozen. Persistent cycling at theta=2, 3 states, mean degree 4 fits the known persistent regime: rotating cycles form by chance at low threshold. Not a numerical artifact, but the measured feature is close to trivial.

Citations:
- Fisch, Gravner, Griffeath (1991) Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39
- Fisch, Gravner, Griffeath (1991) Cyclic cellular automata in two dimensions, Spatial Stochastic Processes
- Limiting behavior of 3-color excitable media on arbitrary graphs, arXiv:1610.07320
- Durrett & Griffeath, Asymptotic behavior of excitable cellular automata / Medley of spirals from cyclic cellular automata
- https://en.wikipedia.org/wiki/Cyclic_cellular_automaton

Full reply text:

**Verdict: KNOWN.** The first pass called this a variant. I'm overturning that slightly: it isn't even a variant. It is the standard threshold cyclic cellular automaton, here with 3 states, threshold 2, random asynchronous updates, and a random network instead of a grid.

**Why this is the same system.** In the standard model, each cell moves forward one step through a fixed cycle of colours if there are at least Threshold cells of the next color within its neighbour set. Your rule matches this except for the setting. The threshold version has been part of the model from the start: more general forms of this type of rule also include a threshold parameter. Durrett and Griffeath studied the model family parametrized by the range ρ of interaction, l p shape of its neighbor set, threshold θ for contact updating, and number κ of possible states per site. Fisch, Gravner and Griffeath mapped this parameter space in "Threshold-range scaling of excitable cellular automata" (1991). There is also a paper on exactly your colour count on general networks: "Limiting behavior of 3-color excitable media on arbitrary graphs" (arXiv:1610.07320). It builds on Clustering in the one-dimensional three-color cyclic cellular automaton. That paper studies the threshold-1 version, but the question it asks, whether a 3-colour cyclic network keeps changing or freezes, is the same as yours. On grids, the known results already separate dying, temporary and persistent behaviour: Gravner and Griffeath's 1993 analysis splits the (K, θ) parameter space into a dying phase, a chaotic transient phase, and a persistent phase. The random updating (probability 0.1 per step) is just the standard asynchronous version. Moving from a grid to a random network doesn't create any new behaviour that you measured.

**Trivial or artifact?** The feature is close to trivial.
1. **"Node states cycle" is built into the rule.** The only change the rule allows is t → t+1 (wrapping round), so every change is a cyclic change. The measurement can only tell you whether the system is still active or has frozen.
2. **The knock-out value of −0.5 looks like an "undefined" placeholder.** With interaction removed, no rule is left, so nothing ever changes. The big difference you found (0.943) is mostly "active system versus completely frozen system".
3. **The genuine finding is narrow.** Activity does persist at threshold 2 on a network with average 4 links per node, and it holds up at 4× size and 3× run length. That fits the known persistent regime: with low thresholds and few states, cycles that keep rotating can form by chance from the random start. On a random network, small loops of nodes in all three states can keep rotating. Nodes with fewer than 2 links can never change at all, so some freezing is expected too.

So the behaviour is a real, ongoing excitable-medium process, not a numerical failure. But the measured feature doesn't show anything beyond the textbook behaviour.

---

## review:045

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) each cell copies the type of a random neighbour (every step); (2) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (every step).

OBSERVED (features that passed the simulation checks):
  g_corrlen, g_n_domains, g_moran, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_corrlen (long spatial correlation length): median 5.0 with interaction, 1.0 with interaction knocked out (difference 4.0, robust spread 0.12); difference at 4x area 4.0, at 3x run length 4.0
  - g_n_domains (many separate domains): median 3.5259 with interaction, 7.3524 with interaction knocked out (difference -3.8265, robust spread 0.2176); difference at 4x area -4.1255, at 3x run length -3.841
  - g_moran (like-valued cells cluster spatially): median 0.8219 with interaction, -0.0025 with interaction knocked out (difference 0.8244, robust spread 0.0165); difference at 4x area 0.8243, at 3x run length 0.8143
  - g_reversals (cells flip back and forth): median 0.5017 with interaction, 0.0 with interaction knocked out (difference 0.5017, robust spread 0.01); difference at 4x area 0.5008, at 3x run length 0.4277

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_corrlen 1.0; g_n_domains 6.367; g_moran 0.3305; g_reversals 0.6472
  - rules 2 only: g_corrlen 3.0; g_n_domains 3.3498; g_moran 0.7469; g_reversals 0.4484
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Noisy Griffeath cyclic cellular automaton '313' (R1/T3/C3/Moore) combined with voter-model copying, i.e. a three-state cyclic voter / spatial rock-paper-scissors (cyclic Lotka-Volterra) lattice model

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Rule 2 alone is exactly Griffeath's named 313 CCA, documented as self-organising into spirals and domains. Rule 1 (synchronous random-neighbour copying) adds voter-model noise and coarsening, giving a standard three-state cyclic voter / RPS-type system. The measured features (clustering, fewer and larger domains, sustained cycling) are generic for this family. 'Reversals' in a 3-cycle largely counts cyclic advancement. The knock-out is a frozen random grid because both rules read neighbours, so the sub-program comparison is more informative. The first pass called it a VARIANT; the twist is only in the rule combination and the magnitude of the measures, not in the phenomenology, so KNOWN.

Citations:
- http://psoup.math.wisc.edu/mcell/rullex_cycl.html
- https://github.com/setanarut/ccae
- https://arxiv.org/pdf/nlin/0507067
- Szabo & Szolnoki, Three-state cyclic voter model extended with Potts energy, PRE 65, 036115 (2002)
- https://arxiv.org/html/1408.6828v1
- https://doi.org/10.3390/g7030024

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass slightly. The measured behaviour (cyclically advancing, spatially clustered domains with cells flipping repeatedly) is documented. Only the exact pairing of rules is new, and you said that doesn't count.

**What the rules are:**
- **Rule 2 alone is a textbook model.** It is Griffeath's cyclic cellular automaton with 3 types, threshold 3 and an 8-cell neighbourhood, the named rule "313" (R1/T3/C3/NM). It is described as "the only three-color rule in the range-1 CCA rulespace that exhibits significant self-organization", with a number of distinct spiral types. The general behaviour is documented: from a random start, waves of excitation self-organize into large-scale spiral pairs that stabilize in a locally periodic state. Lists of named rules confirm 313 | David Griffeath | 1/3/3/M. Your "rule 2 only" run (Moran's I, a clustering score, 0.75; about 0.45 reversals; correlation length 3) fits this.
- **Rule 1 is the voter model** (copy a random neighbour), run with all cells updating at once. Adding it makes the 313 rule noisy, so the whole system is a stochastic three-type cyclic-dominance (rock-paper-scissors) model on a lattice. That family is heavily studied. It includes Szabó & Szolnoki, "Three-state cyclic voter model extended with Potts energy", Tainaka's three-candidate voter model and the Reichenbach–Mobilia–Frey cyclic Lotka-Volterra model. The 2014 review covers pattern formation, the impact of mobility, and the spontaneous emergence of cyclic dominance. It also notes that the characteristic length scale of spiral waves increases with increasing mobility.

**Your numbers against the literature:**
- The full program has a slightly larger correlation length (5 vs 3) and stronger clustering (0.82 vs 0.75) than rule 2 alone. Coarse domains keep cycling, with about half of cells reversing.
- That is what you'd expect when the cyclic rule's waves sit on top of voter-model coarsening. In this family, mixing and noise are known to change the domain or wave size. That isn't a new kind of behaviour.
- None of the four measured features (correlation length, domain count, clustering, reversals) captures anything beyond "clustered domains that keep cycling", which is already documented. Getting a spread of values across sub-programs is a change of degree, not new phenomenology.

**Genuine or artifact:** genuine, but standard. Three caveats:
1. The knock-out removes both rules, because both read neighbours. So it is just the frozen random starting grid (clustering about 0, reversals 0), a weak baseline. The sub-program comparison is the more useful one.
2. In a 3-type cycle, "reversals" partly just counts cells advancing round the cycle t→t+1→t+2. It isn't necessarily flipping back and forth.
3. A correlation length of 5 on a 64×64 grid is modest. Still, the effects are robust across 4× area and 3× run length, and nothing diverged or froze.

---

## review:047

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_largest_domain, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 5.3701 with interaction, 7.3658 with interaction knocked out (difference -1.9957, robust spread 0.2547); difference at 4x area -2.3367, at 3x run length -2.2084
  - g_largest_domain (one cell-type fills most of the lattice): median 0.814 with interaction, 0.0076 with interaction knocked out (difference 0.8064, robust spread 0.0164); difference at 4x area 0.8742, at 3x run length 0.8568
  - g_reversals (cells flip back and forth): median 0.5601 with interaction, 0.0 with interaction knocked out (difference 0.5601, robust spread 0.0112); difference at 4x area 0.5577, at 3x run length 0.6379

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 7.3457; g_largest_domain 0.0085; g_reversals 0.0
  - rules 2 only: g_n_domains 6.3994; g_largest_domain 0.0574; g_reversals 0.4143
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Mobility/long-range-mixing-induced global oscillation and loss of coexistence in cyclic-dominance (rock-paper-scissors / cyclic CA, Greenberg-Hastings-type) lattice dynamics

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Rule 1 alone is the 3-colour, threshold-1, 8-neighbour cyclic CA (locally periodic, quenched, fine-grained; zero reversals by construction). Rule 2 (Schelling relocation to random empty cells anywhere) acts as long-range mixing. Mixing is documented to synchronise cyclic dominance into global oscillations or near-fixation, which explains the ~81% largest single-type domain and the fewer domains. g_reversals is largely an artifact of measuring at fixed cells: vacating and refilling a cell reads as t->0->t, and rule 2 alone already gives 0.41. Exact rule pairing (empty state inside the cycle, Schelling-driven jumps) not found in the literature, but the behaviour is standard. Unresolved: whether the dominant type cycles over time (global oscillation) or stays fixed (near-absorbing); either way the behaviour is documented.

Citations:
- Reichenbach, Mobilia & Frey, Nature 448, 1046 (2007) https://www.nature.com/articles/nature06095
- Biodiversity in models of cyclic dominance is preserved by heterogeneity in site-specific invasion rates (review of shortcut/mobility-induced global oscillations) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5137108/
- Fisch, Gravner & Griffeath, Cyclic cellular automata in two dimensions (1991); phenomenology summarised in https://arxiv.org/html/nlin/0507067
- Szabo, Szolnoki & Sznaider, PRE 76, 051921 (2007) https://arxiv.org/pdf/0710.3371

Full reply text:

**Verdict: KNOWN.** This overturns the first pass in the stricter direction. What the system does is a documented behaviour: mobility, or long-range mixing, stops cyclic dominance from forming local spatial patterns, and the lattice ends up dominated by one type at a time. One of the three features (reversals) is largely a measurement artifact.

**What the rules are**
- **Rule 1 alone** is the classic Greenberg–Hastings/Griffeath cyclic cellular automaton with 3 states, threshold 1 and 8 neighbours. Here the empty state is simply one of the three steps in the cycle.
- That case is well characterised. With an 8-neighbour neighbourhood, T=1 and T=2 are both locally periodic, but T=2 produces spiral waves, while T=1 "quenches" incoherent local oscillations.
- This matches your rule-1-only numbers: many small domains (n_domains 7.3), tiny largest domain (0.0085), zero reversals. Zero reversals is guaranteed, because cells only ever move forward 0→1→2→0, so a cell can't flip straight back. Nothing pathological there.
- **Rule 2** adds Schelling-style relocation to *random* empty cells anywhere on the lattice. That is long-range mixing, not movement to a neighbouring cell.

**Why the main feature (largest domain 0.81) is known**
- Reichenbach, Mobilia and Frey (Nature 2007) report that in cyclic dominance when mobility exceeds a certain value, biodiversity is jeopardized and lost. In contrast, below this critical threshold all subpopulations coexist and an entanglement of travelling spiral waves forms.
- More specifically, on lattices, if shortcuts or long-range interactions are introduced to the lattice, or if the original lattice is simply replaced by a small-world network, then initially locally occurring oscillations do synchronize, leading to global oscillations and to the accidental extinction of one species in the loop. Also, the introduction of mobility of players, also promotes the emergence of global oscillations that jeopardize biodiversity.
- Random jumps to empty cells are exactly such long-range links. One type taking about 81% of the lattice, with fewer domains, is what a global oscillation looks like in a snapshot: one type dominates, then gives way to the type that beats it. It is not a new ordered phase.
- It also matches the mean-field picture: for finite population the system evolves into one of the homogeneous state via a Moran process.

**Real pattern or artifact?**
- **Largest domain and domain count:** a genuine emergent effect of the two rules interacting. Neither rule alone produces it, and it holds at 4× area and 3× run length. But it is the known mechanism, not something new.
- **Reversals (0.56):** mostly an artifact of measuring at fixed cells. When an occupant leaves, its cell becomes empty (type 0), and when another arrives it changes again. So a cell can read t→0→t from moving agents, with no actual oscillation. Rule 2 alone already gives 0.41. Only the extra ~0.15 could come from the interaction.

**What I couldn't settle**
- I didn't find a paper with this exact pairing: Schelling relocation driven by unhappiness, with empty cells as a member of the cycle. As you said, that alone doesn't make it new.
- One open question is whether the dominant type cycles over time (global oscillation) or stays put (approaching a frozen absorbing state). Either way the behaviour is documented, so it wouldn't change the verdict. A time series of which type is largest would tell them apart.

Searches I ran: mobility in rock-paper-scissors (Reichenbach et al.), cyclic automaton threshold/range phase diagrams (Fisch–Gravner–Griffeath), Schelling combined with cyclic dominance, and cyclic automata with mixing or small-world links leading to global synchronisation.

---

## review:049

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.7348 with interaction, 7.3575 with interaction knocked out (difference -2.6227, robust spread 0.2418); difference at 4x area -2.7592, at 3x run length -2.7131
  - g_corrlen (long spatial correlation length): median 3.0 with interaction, 1.0 with interaction knocked out (difference 2.0, robust spread 0.08); difference at 4x area 2.0, at 3x run length 2.0
  - g_moran (like-valued cells cluster spatially): median 0.6221 with interaction, -0.0056 with interaction knocked out (difference 0.6277, robust spread 0.0169); difference at 4x area 0.6371, at 3x run length 0.6277

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 4.1584; g_corrlen 4.0; g_moran 0.7744
  - rules 2 only: g_n_domains 6.8216; g_corrlen 1.0; g_moran 0.0835
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Threshold cyclic cellular automaton / Greenberg-Hastings-type excitable CA (range 1, Moore, kappa=3, theta=3; turbulence/fixation domain formation), with Schelling-type random relocation acting as mixing that coarsens less

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Measured features (domain count, correlation length, Moran's I) are generic clustering signatures. Rule 1 alone (cyclic threshold CA, empty is part of the cycle) yields stronger clustering (Moran 0.77, corrlen 4); the Schelling relocation only dilutes it (full program sits between the two sub-programs on every feature), consistent with known mobility-induced breakdown of cyclic-competition patterns. No synergistic phenomenology. Mild overturn: VARIANT -> KNOWN, since the only difference is in the rules, not the behaviour. Pattern is genuine but weak (corrlen about 3 cells) and possibly near-fixated; whether it moves over time was not measured.

Citations:
- Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata, Statistics and Computing 1 (1991) 23-39
- Fisch, Gravner & Griffeath, Cyclic Cellular Automata in Two Dimensions (1991), https://link.springer.com/chapter/10.1007/978-1-4612-0451-0_8
- Shalizi et al., Quantifying Self-Organization in Cyclic Cellular Automata, https://arxiv.org/abs/nlin/0507067
- Szczesny, Mobilia & Rucklidge, The influence of mobility rate on spiral waves in spatial rock-paper-scissors games, https://arxiv.org/abs/1608.08932

Full reply text:

**Verdict: KNOWN.** This partly overturns the first pass, which called it a variant. The three measured features are just generic spatial clustering. That clustering comes from a well-studied rule, and the second rule only weakens it. Nothing new appears when the two rules run together.

**What the system actually is**
- **Rule 1 is a cyclic cellular automaton with a threshold.** A cell is "consumed" by its successor type when enough neighbours have that type. Because the cycle includes type 0 (empty → 1 → 2 → empty), "empty" behaves as just a third colour here. This is the cyclic cellular automaton / Greenberg–Hastings family studied by Fisch, Gravner and Griffeath, where each family … [is] parametrized by the range ρ of interaction, l p shape of its neighbor set, threshold θ for contact updating, and number κ of possible states per site. Our case is range 1 with 8 neighbours, threshold 3, 3 states, and random 10%-per-step updating instead of all cells updating together.
- **What that family is known to do.** It has three generic forms of long-term behavior, depending on the size of the threshold relative to the range. At high thresholds it forms homogeneous blocks of solid colors, which are completely static — so-called fixation behavior. In the 4-colour, range-1 version, at T=3, one encounters "turbulence," which is actually meta-stable. Our system sits near that boundary between turbulence and fixation. Either way you get single-colour patches: a few domains, a correlation length of a few cells, and a high Moran's I. That matches the rule-1-only run (Moran 0.77, correlation length 4).
- **Rule 2 is Schelling relocation with tolerance 0.3.** Unhappy cells jump to a random empty site. On top of the cyclic rule, this acts as long-range mixing. Mixing is known to break up domains in cyclic-competition lattice models: below a certain mobility threshold species coexist by forming intriguing spiraling patterns … whereas there is loss of biodiversity above that threshold. Pair-exchange mobility is also known to cause far-field breakup of spiral waves.

**Why there is no new phenomenology**
- Every measure for the full program falls between the two single-rule runs:
  - Moran's I: rule 1 only 0.77, full 0.62, rule 2 only 0.08.
  - Correlation length: 4, then 3, then 1.
  - Domain count: 4.2, then 4.7, then 6.8.
- So the coupling just weakens rule 1's clustering. Nothing appears that neither rule produces alone.
- Relocation reducing order in a cyclic lattice model is documented qualitatively, as cited above. The "fails 3 robust spreads against every sub-program" check only shows the numbers move between the two; it doesn't show a new behaviour.
- The knock-out comparison only shows that clustering needs the neighbour-reading rule, which is true by construction.

**Genuine or artifact?** It's a genuine but generic pattern: patches formed by cyclic consumption. It isn't random colouring, a field that has blown up numerically, or an inert coupled field. It may be close to a frozen fixation-like state (blocks of colour that stop changing), and a correlation length of 3 cells on a 64×64 grid is weak structure. I didn't measure whether the pattern moves over time (turbulent) or freezes.

**Searches run:** cyclic CA threshold-range phase diagrams (Fisch–Gravner–Griffeath, plus a later paper quantifying their self-organisation), and Schelling or mobility combined with rock-paper-scissors cyclic dominance. I found no paper with this exact Schelling + cyclic-CA combination. Per the brief, that alone doesn't make it a discovery.

---

## review:053

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a random network of 200 nodes (average 4 links each). Nodes have 4 types (random start). Rules, applied in order each step: (1) a node of type t becomes t+1 (cyclically) if at least 1 neighbours are t+1 (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (every step).

OBSERVED (features that passed the simulation checks):
  n_edge_turnover, n_oscillation

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - n_edge_turnover (edges rewire): median 0.982 with interaction, 0.0 with interaction knocked out (difference 0.982, robust spread 0.0196); difference at 4x area 0.9941, at 3x run length 0.982
  - n_oscillation (node states oscillate over time): median 0.5207 with interaction, 0.0 with interaction knocked out (difference 0.5207, robust spread 0.1544); difference at 4x area 0.3567, at 3x run length 0.2706

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: n_edge_turnover 0.0; n_oscillation 0.8184
  - rules 2 only: n_edge_turnover 0.7411; n_oscillation 0.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Cyclic-dominance (rock-paper-scissors / cyclic CA) dynamics on an adaptive network with Holme-Newman-type same-type rewiring (adaptive RPS / coevolving cyclic voter model)

Trivial or artifact: True; confidence: medium

Twist (if any): 

Notes: The two measured features (edge turnover, state oscillation) are each produced by a single rule on its own (rule 2 alone gives turnover 0.74; rule 1 alone gives oscillation 0.82), so the knock-out and sub-program tests only show the two rules adding up, not interacting. Combining the rules lowers oscillation, and it falls further at larger size and longer runs, consistent with a slow drift towards the documented frozen, fragmented phase at high rewiring strength. The 4 types, synchronous threshold-1 update and every-step rewiring change the rules but add no new behaviour. The first pass's 'variant' verdict is downgraded to KNOWN because no twist in the behaviour was identified.

Citations:
- Demirel, Prizak, Reddy, Gross, 'Cyclic dominance in adaptive networks', EPJ B 84, 541 (2011), arXiv:1011.1124
- 'Adaptive cyclically dominating game on co-evolving networks: numerical and analytic results', EPJ B (2017), arXiv:1605.03370
- Holme & Newman, 'Nonequilibrium phase transition in the coevolution of networks and opinions', PRE 74, 056108 (2006)
- Lutz, Risau-Gusman, Arenzon, 'Intransitivity and coexistence in four species cyclic games', J. Theor. Biol. 317, 286 (2013)

Full reply text:

**Verdict: KNOWN.** This overturns the first pass's "variant" verdict, but only mildly. Both measured behaviours are documented properties of models where node states and network links change together and states beat each other in a cycle (rock-paper-scissors style). Neither measured feature shows anything that only appears when the two rules are combined.

**What the system is**
- **Rule 1** is a cyclic cellular automaton: a node moves to the next type if at least one neighbour already has it (the Fisch–Gravner–Griffeath / Greenberg–Hastings family). I'm citing this from memory; I didn't search for it. With 4 types it's a four-species cyclic game, which has also been studied: Lütz, A., Risau-Gusman, S. & Arenzon, J. Intransitivity and coexistence in four species cyclic games.
- **Rule 2** is the Holme–Newman rewiring step: a vertex i is selected at random and one of its edges (in this case the edge (i,j)) is rewired to a new vertex j' holding the same opinion as i.
- **The combination has been published.** Demirel, Prizak, Reddy & Gross ran rock-paper-scissors on a network where the loser either adopts the opinion of the winner or rewires the link. They found that nonequilibrium phase transitions occur as a function of the rewiring strength. The transitions separate four distinct phases which differ in the observed dynamics of opinions and topology. A summary of that work lists the outcomes as stationary coexistence, oscillatory cycling, heteroclinic consensus, and fragmentation—depending on adoption versus segregation probability.
- **A closely related model** (ARPS, EPJ B 2017) shows two phases. In the active one, one connected network of agents using different strategies who are continually interacting and taking adaptive actions. A frozen phase for p > p_c has three separate clusters of agents using only R, P, and S, respectively with terminated adaptive actions.

Your observations, constant edge turnover plus states cycling over time, match the documented active, cycling phase. The differences are 4 types instead of 3, a threshold-1 update applied to all nodes at once, and rewiring on every step. These are changes to the rules, not new behaviour.

**Trivial or emergent?** The features look trivial: each is just the sum of what the two rules do separately.
- Rule 2 alone already gives edge turnover (0.74), and rule 1 alone already gives oscillation (0.82). The knock-out test removes both, so beating it says nothing about the rules interacting.
- The "no sub-program matched on every feature" check passes automatically, because each sub-program scores exactly 0 on one of the two features.
- Combining the rules *lowers* oscillation (0.82 → 0.52). It falls further with size (0.36 at 4× area) and with run length (0.27 at 3× length, just above the half-size cutoff). This looks like a slow drift towards the documented frozen, split-up state. Same-type rewiring on every step corresponds to a very high rewiring strength, which is where the literature puts fragmentation.
- If the 0.98 turnover counts edges that have *ever* rewired, it would reach that level even if the network froze later.

What would be worth checking: does oscillation eventually reach zero in long runs, does the network end up split into single-type clusters, and does the sustained or frozen outcome depend on the number of types or the threshold?

**Searches run:** "cyclic dominance adaptive network rewiring rock-paper-scissors coevolution fragmentation" and "Holme Newman nonequilibrium phase transition coevolution networks opinions rewiring".

---

## review:060

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 4 types (random start). Agents have 3 types (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (every step).

OBSERVED (features that passed the simulation checks):
  g_corrlen, g_n_domains, g_activity_trend, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_corrlen (long spatial correlation length): median 10.0 with interaction, 1.0 with interaction knocked out (difference 9.0, robust spread 0.22); difference at 4x area 10.0, at 3x run length 7.0
  - g_n_domains (many separate domains): median 3.5773 with interaction, 7.669 with interaction knocked out (difference -4.0917, robust spread 0.2249); difference at 4x area -4.3613, at 3x run length -4.2077
  - g_activity_trend (activity rising): median -0.5756 with interaction, 0.0 with interaction knocked out (difference -0.5756, robust spread 0.0312); difference at 4x area -0.6441, at 3x run length -0.5527
  - g_moran (like-valued cells cluster spatially): median 0.8744 with interaction, -0.0026 with interaction knocked out (difference 0.877, robust spread 0.0175); difference at 4x area 0.8691, at 3x run length 0.8555

### Stage 3 reply: verdict KNOWN, standard name: Cyclic cellular automaton (Griffeath CCA), rule 1/3/4/Moore: metastable 'turbulent' regime with domain coarsening / possible spiral entrainment

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Only the cell rule is active; it is exactly the 4-colour, range-1 Moore, threshold-3 CCA, the case studied explicitly by Shalizi et al. citing Fisch-Gravner-Griffeath (T=3 = metastable turbulence that may become spirals; T>=4 fixation). The 400 agents are inert (no coupling) and do not affect classification. Measured clustering, ~10-cell correlation length, fewer domains and falling activity match documented self-organisation and relaxation on a finite lattice. Activity trend is negative despite the 'rising' label. Not a variant: no new phenomenology.

Citations:
- Fisch, Gravner & Griffeath, 'Cyclic cellular automata in two dimensions' (1991), https://math.ucdavis.edu/~gravner/papers/cca.pdf
- Shalizi, Haslinger, Rouquier, Klinkner & Moore, 'Quantifying Self-Organization in Cyclic Cellular Automata', arXiv:nlin/0507067
- MCell rules lexicon, Cyclic CA, http://psoup.math.wisc.edu/mcell/rullex_cycl.html

Full reply text:

**Verdict: KNOWN.** I'm changing the first pass's call from "variant" to "known".

**What the system is.** Only one rule actually runs: a cell moves on to the next of 4 types when at least 3 of its 8 neighbours already have that type. That is Griffeath's standard cyclic cellular automaton (CCA). In the usual range/threshold/colours/neighbourhood notation it is rule 1/3/4/Moore. In the CCA, each update a cell's color advances by 1 if there are at least Threshold cells of the next color within its neighbour set. The 400 agents do nothing to the cells: no rule reads them or is read by them. So they don't change the classification. The real question is whether the agent rules matter at all, and here they don't.

**The exact settings are already published.** Shalizi et al. ("Quantifying Self-Organization in Cyclic Cellular Automata") studied this parameter set as one of their main examples. In this and the subsequent figures, the number of colors is fixed to 4, the rule uses a range-one box-shaped neighborhood, and the initial conditions are a uniform random mixture of colors. Citing Fisch, Gravner and Griffeath, they describe threshold 3 this way: At T=3, one encounters "turbulence," which is actually meta-stable — spiral waves can form and entrain the entire CA, but this does not always happen on finite lattices, and in any case the turbulent phase can persist for very long times. Threshold 3 sits between T=2 produces spiral waves and Fixation occurs with T≥4.

**How the measurements fit.**
- **Correlation length ≈10, Moran's I ≈0.87, fewer domains (≈3.6 vs ≈7.7):** this is the expected self-organisation from random noise into large coloured patches and waves, rather than random colouring.
- **Activity falling (trend ≈ −0.58):** this fits a metastable turbulent phase settling down, or patches of spiral waves taking over regions of a finite 64×64 lattice. The feature label says "activity rising", but the measured sign shows activity is going down.
- **Knock-out comparison:** removing the neighbour rule freezes the random start, which gives correlation length 1 and Moran's I ≈ 0. So the order comes entirely from the cell rule.

**Is it real or an artifact?** The pattern is a genuine emergent one. It is not a random colouring, a numerical blow-up or a trivially frozen state, and it holds up at 4× area and 3× run length. The only "artifact" is the agent layer, which is inert decoration.

**No new behaviour found.** The domains, clustering, correlation length and relaxing activity are all part of the documented T=3 behaviour. I didn't find any twist in the observed behaviour that the literature doesn't already describe.

---

## review:062

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.3 per step); (2) neighbouring cells swap contents (every step).

OBSERVED (features that passed the simulation checks):
  g_corrlen, g_n_domains, g_moran, g_cyclic_changes

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_corrlen (long spatial correlation length): median 6.0 with interaction, 1.0 with interaction knocked out (difference 5.0, robust spread 0.14); difference at 4x area 5.0, at 3x run length 4.0
  - g_n_domains (many separate domains): median 5.4048 with interaction, 7.345 with interaction knocked out (difference -1.9402, robust spread 0.255); difference at 4x area -1.7152, at 3x run length -1.5566
  - g_moran (like-valued cells cluster spatially): median 0.6708 with interaction, 0.0083 with interaction knocked out (difference 0.6625, robust spread 0.0162); difference at 4x area 0.6181, at 3x run length 0.5373
  - g_cyclic_changes (cells cycle through states above chance): median 0.3924 with interaction, -0.5 with interaction knocked out (difference 0.8924, robust spread 0.0178); difference at 4x area 0.8981, at 3x run length 0.519

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_corrlen 4.0; g_n_domains 3.9597; g_moran 0.774; g_cyclic_changes 0.4831
  - rules 2 only: g_corrlen 1.0; g_n_domains 7.3553; g_moran -0.0038; g_cyclic_changes 0.0004
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Stochastic threshold cyclic cellular automaton (Fisch-Gravner-Griffeath CCA / Greenberg-Hastings family) with neighbour-exchange mobility, as in spatial rock-paper-scissors with mobility (Reichenbach-Mobilia-Frey)

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Changes the first pass from VARIANT to KNOWN. Rule 1 alone already gives cyclic wave domains (documented CCA self-organisation). Swapping makes the domains larger (higher corrlen) and noisier (slightly lower Moran and cycling, more domains), which matches the documented finding that mobility lengthens spiral wavelength below the critical mobility. Not an artifact: no freezing, not random colouring, holds up at larger size and longer runs. Unchecked: whether these are true spirals or drifting bands, and where threshold 3 / 3 colours / 8 neighbours sits in the CCA phase diagram; neither would make the behaviour new.

Citations:
- Fisch, Gravner & Griffeath, Threshold-range scaling of excitable cellular automata, Statistics and Computing 1:23-39 (1991)
- Fisch, Gravner & Griffeath, Cyclic cellular automata in two dimensions, Spatial Stochastic Processes (1991)
- Reichenbach, Mobilia & Frey, Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, Nature 448:1046 (2007)
- Szolnoki et al., Cyclic dominance in evolutionary games: a review, J R Soc Interface 11:20140735 (2014)

Full reply text:

**Verdict: KNOWN.** I'm changing the first pass's "variant" label to "known". Both halves of this system are standard models. The only new thing is the exact rule combination, and by the task's own criteria that doesn't count.

**What the system is.** Rule 1 alone is a cyclic cellular automaton with a threshold, made random by firing with probability 0.3. In Fisch, Gravner & Griffeath's framework, each cell of a two-dimensional lattice is painted one of κ colors, arranged in a 'color wheel'. The colors advance (k to k+1 mod κ)… by contact with at least a threshold number of successor colors in a prescribed local neighborhood. A system in which all colors update by contact is called a cyclic cellular automation (CCA). They report that started from random configurations the same rules exhibit complex self-organization, typically characterized by nucleation of locally periodic 'ram's horns' or spirals. Corresponding random processes give rise to a variety of 'forest fire' equilibria that display large-scale stochastic wave fronts. Your setup (3 colours, threshold 3, 8 neighbours, random firing) is one point in that parameter space.

**What the swapping adds.** Rule 2 is the mobility (neighbour-exchange) process from spatial rock-paper-scissors models. Reichenbach, Mobilia & Frey found that a critical threshold of motility is identified, above which biodiversity is lost, but below which an entanglement of travelling spiral waves emerges. They add that this phenomenon is robust, it does not depend on the details of cyclic competition or spatial environment, and that the spirals' wavelength λ rises with the individuals' mobility.

**How the four features match the literature:**
- **Spatial clustering (g_moran) and cycling (g_cyclic_changes):** Rule 1 alone already produces both (0.77 and 0.48). This is the documented CCA self-organisation into rotating wave domains.
- **Correlation length (g_corrlen):** Adding the swaps raises it from 4 to 6. That matches the documented result that mobility makes spiral wavelengths longer.
- **Domain count (g_n_domains):** It goes slightly up (3.96 to 5.40), while clustering and cycling go slightly down. That fits mobility acting as noise that blurs wavefronts: below the critical mobility, waves persist but are messier.
- **Swaps alone:** These give the expected mixed, uncorrelated state (clustering about 0, correlation length 1).

Nothing in these four measurements goes beyond "a cyclic excitable lattice, with mixing, forms wave domains whose size changes with the mixing rate".

**Trivial or artifact?** No. The pattern isn't frozen: cells cycle well above chance, and the knock-out gives −0.5. It isn't random colouring: clustering is 0.67 against 0.008 for the knock-out. And the differences hold up at 4× area and 3× run length.

**Caveats:**
- I haven't seen snapshots, so I can't confirm these are true spirals rather than drifting cyclic bands. Both are documented CCA behaviours.
- The sources I found don't check whether threshold 3 with 3 colours on 8 neighbours lands in the spiral phase of the deterministic CCA. The random 0.3 firing could move it between phases.
- Neither issue would make the behaviour new.

---

## review:064

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 2 neighbours already are t+1 (every step); (2) neighbouring cells swap contents (every step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_moran, g_activity

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 4.2194 with interaction, 7.6723 with interaction knocked out (difference -3.4529, robust spread 0.2378); difference at 4x area -2.9926, at 3x run length -3.9028
  - g_corrlen (long spatial correlation length): median 4.5 with interaction, 1.0 with interaction knocked out (difference 3.5, robust spread 0.3706); difference at 4x area 3.5, at 3x run length 4.5
  - g_moran (like-valued cells cluster spatially): median 0.7982 with interaction, 0.002 with interaction knocked out (difference 0.7962, robust spread 0.016); difference at 4x area 0.7116, at 3x run length 0.8267
  - g_activity (many cells change per step): median 0.4775 with interaction, 0.0 with interaction knocked out (difference 0.4775, robust spread 0.0095); difference at 4x area 0.4759, at 3x run length 0.4701

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 4.4427; g_corrlen 2.0; g_moran 0.4701; g_activity 0.9897
  - rules 2 only: g_n_domains 7.6725; g_corrlen 1.0; g_moran 0.0057; g_activity 0.6376
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Cyclic cellular automaton (Griffeath CCA, 1/2/4/M; spiral/travelling-wave regime) with stirring/mobility, cf. mobility-enlarged spiral waves in cyclic-dominance lattice models

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Rule 1 alone is the documented CCA 1/2/4/M (persistent locally periodic waves; matches activity ~0.99). Adding swaps raises correlation length and Moran's I and halves activity, consistent with the documented growth of spiral/wave size with mobility. Only the exact pairing (deterministic synchronous CCA plus swaps every step) is not in the literature; the behaviour is. Patterns are sustained, clustered and stable across sizes, so not an artifact. Caveat: correlation length ~4.5 on 64x64 is modest; not checked by eye whether these are clean spirals or wave fragments. Sources checked via search excerpts, not full papers.

Citations:
- Fisch, Gravner, Griffeath (1991) Cyclic cellular automata in two dimensions / Threshold-range scaling of excitable cellular automata, Stat. Comput. 1:23-39
- Shalizi et al., Quantifying Self-Organization in Cyclic Cellular Automata, arXiv nlin/0507067 (range-1 Moore: T=2 spiral waves)
- Reichenbach, Mobilia, Frey (2007) Mobility promotes and jeopardizes biodiversity in rock-paper-scissors games, Nature 448:1046
- Reichenbach, Mobilia, Frey (2008) Self-organization of mobile populations in cyclic competition, J. Theor. Biol. 254:368
- Roman, Konrad, Pleimling (2012) Cyclic competition of four species: domains and interfaces, J. Stat. Mech. P07014

Full reply text:

**Verdict: KNOWN.** The rule is a known model with mixing added, and its behaviour is the documented spiral-wave / travelling-wave state. I'm keeping the first pass's "variant" reading, but I'd tighten it to KNOWN, because I found no new behaviour that the literature leaves unexplained.

**What the system is.** Rule (1) is exactly Griffeath's cyclic cellular automaton (CCA). In standard notation (range / threshold / number of colours / neighbourhood) it is 1/2/4/M. Each update a cell's color advances by 1 if there are at least Threshold cells of the next color within its neighbour set. For range-1 Moore neighbourhoods the phenomenology is documented: T=1 and T=2 are both locally periodic, but T=2 produces spiral waves, while T=3 is "turbulent" and Fixation occurs with T≥4. Threshold 2 is therefore the persistent, locally periodic wave regime. That matches what you measured for rule 1 alone: activity of about 0.99 (nearly every cell advances every step) with moderate clustering (Moran ≈ 0.47).

Rule (2), neighbour swapping, is mobility (stirring) added to a cyclic-dominance system. That setting is well studied. Reichenbach, Mobilia and Frey found that below a critical mobility all subpopulations coexist and an entanglement of travelling spiral waves forms. They also report that this phenomenon is robust, it does not depend on the details of cyclic competition or spatial environment. As mobility increases, the typical size of the spiral waves rises, up to a critical value; beyond that value one species takes over. Four-species versions of cyclic competition (Szabó & Sznaider 2004; Roman, Konrad & Pleimling 2012) are also established (Roman, A., Konrad, D., Pleimling, M., 2012. Cyclic competition of four species: domains and interfaces.).

**How your data fit.** Adding swapping to the CCA:
- raises the correlation length (2 → 4.5) and Moran's I (0.47 → 0.80);
- lowers the domain count slightly;
- halves the activity (0.99 → 0.48), but it stays far from zero.

That is the documented effect: mobility makes the cyclic wave patterns coarser and larger. The halved activity is what you'd expect from a longer wavelength, since fewer cells sit on a wavefront at any moment. Most of the knock-out gap is just the CCA alone. Swapping only enlarges the patterns, which is the known scaling of spiral size with mobility.

**Trivial or artifact? No.**
- It is not frozen: activity of about 0.48 is sustained, and the CCA at threshold 2 is known not to fixate.
- It is not random colouring: Moran ≈ 0.80, against ≈ 0 for the knock-out and for swapping alone.
- It holds up at 4× area and 3× run length.
- One caveat: a correlation length of about 4.5 cells on a 64×64 grid is a modest scale. I couldn't check by eye that these are clean spirals rather than shorter, ragged wave fragments. Either way it falls within the documented CCA-plus-mobility behaviour.

**What I didn't find.** No paper studies this exact pairing of a deterministic synchronous CCA with swaps every step. Mobility has been studied in stochastic rock-paper-scissors / May–Leonard models, not in Griffeath's CCA. That gap is about the rule combination, not the behaviour, so it doesn't justify KNOWN-VARIANT or LITERATURE-SILENT.

**Searches:** cyclic CA threshold/range/Moore phenomenology (Fisch–Gravner–Griffeath), and cyclic dominance with mobility/exchange and spiral waves, including four species. I relied on search-result excerpts and did not read the full papers.

---

## review:069

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 1 neighbours already are t+1 (with probability 0.1 per step); (2) a type-0 cell spontaneously becomes type 1 (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  g_n_types, g_type_entropy, g_interface, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_types (number of surviving types): median 3.0 with interaction, 2.0 with interaction knocked out (difference 1.0, robust spread 0.1); difference at 4x area 1.0, at 3x run length 1.0
  - g_type_entropy (cell types evenly mixed): median 0.8868 with interaction, 0.5824 with interaction knocked out (difference 0.3045, robust spread 0.0294); difference at 4x area 0.3063, at 3x run length 0.3057
  - g_interface (long/rough domain boundaries): median 0.5269 with interaction, 0.4447 with interaction knocked out (difference 0.0822, robust spread 0.0194); difference at 4x area 0.0817, at 3x run length 0.0837
  - g_reversals (cells flip back and forth): median 0.4959 with interaction, 0.0 with interaction knocked out (difference 0.4959, robust spread 0.0099); difference at 4x area 0.4951, at 3x run length 0.5953

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_types 3.0; g_type_entropy 1.0; g_interface 0.6402; g_reversals 0.4285
  - rules 2 only: g_n_types 2.0; g_type_entropy 0.5824; g_interface 0.4447; g_reversals 0.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: 3-colour cyclic cellular automaton / Bramson-Griffeath cyclic particle system (threshold 1, Moore neighbourhood), i.e. spatial rock-paper-scissors with one asymmetrically enhanced invasion rate

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Upgraded from KNOWN-VARIANT to KNOWN: rule 1 is the canonical stochastic 3-colour threshold-1 cyclic CA, proven to fluctuate (perpetual cycling, coexistence of all three colours) in 2D. Rule 2 is equivalent to raising the 1-invades-0 rate, an asymmetric RPS that gives uneven but persistent coexistence, which is documented. The knock-out baseline is a trivial frozen state (all 0s turn into 1s, then nothing changes), so the large feature differences partly reflect a degenerate baseline.

Citations:
- Fisch, Gravner & Griffeath (1991), Cyclic Cellular Automata in Two Dimensions, Spatial Stochastic Processes, Progress in Probability 19
- Bramson & Griffeath (1989), Flux and fixation in cyclic particle systems, Ann. Probab. 17
- arXiv:2108.06404 (CCA/GHM review: fluctuation on Z^2 for theta=1, all kappa>=3)
- Frean & Abraham (2001), Rock-scissors-paper and the survival of the weakest, Proc. R. Soc. B 268:1323
- arXiv:1907.11280, Predominance of the weakest species in LV and May-Leonard RPS

Full reply text:

**Verdict: KNOWN.** The first pass called this a variant. I'd go a step further: this is the standard 3-colour cyclic cellular automaton with threshold 1, run in stochastic (partly asynchronous) form, plus one rule that speeds up one of the three cyclic steps. All four measured features are textbook behaviour for that family.

**What the system is**
- **Rule 1 is the cyclic CA / cyclic particle system.** In that model a cell changes to the next colour (k+1 mod κ) once it has at least T neighbours of that colour, and otherwise stays the same. A cell of color k will change its color to k+1 mod κ if there are already at least T cells of that color in its neighborhood... Otherwise, the cell retains its current color. Here κ=3, T=1, with the 8-cell (Moore) neighbourhood. Because each cell only updates with probability 0.1 per step, it sits between the synchronous CCA and the asynchronous cyclic particle system of Bramson & Griffeath. The two versions are known to behave alike: the CCA exhibits the same phase transition as was proved earlier by Bramson and Griffeath for the cyclic particle systems on Z, wherein vertices update their colors asynchronously.
- **The 2D case is settled.** On Z² with θ = 1, Fisch, Gravner and Griffeath proved that the CCA and GHM fluctuate for all κ ≥ 3. "Fluctuate" means every cell keeps cycling through the colours forever and all three colours coexist. That covers g_reversals ≈ 0.5, g_n_types = 3, high type entropy and long, rough domain boundaries (g_interface). Their result is for synchronous updating on an infinite lattice, but the agreement between the synchronous and asynchronous versions cited above suggests it carries over here. In ecology the same dynamics is the spatial rock-paper-scissors model (type y can eat type x at a neighboring site, provided that y = x + 1 mod N).
- **Rule 2 just biases one step.** A type-0 cell spontaneously becoming type 1 is the same as making the "1 invades 0" step faster. That turns it into an asymmetric-rate rock-paper-scissors game, which is a well-studied topic. Frean & Abraham showed that the species that is least competitive is expected to have the largest population. Lattice versions with one species' predation rate changed have also been studied in Lotka-Volterra and May-Leonard forms (one of the species has a reduced predation probability). The full program's entropy (0.887, versus 1.0 with rule 1 alone) shows exactly this kind of uneven but stable coexistence. I didn't check which species ends up most abundant here, but this effect is documented.

**Genuine or artifact?**
- **The emergent behaviour is genuine.** It is the well-known cyclic-dominance dynamics: domains keep invading one another and cells keep cycling.
- **The comparison baseline is trivial.** With the interaction removed, rule 2 turns every 0 into a 1 and then everything freezes. Only two types are left, nothing flips back, and the low interface value just comes from the random starting pattern of 1s and 2s. So the large gaps between the full program and the knock-out partly show that the knock-out is trivial, not that something new has appeared.
- **The sub-program results fit the standard picture.** Rule 1 alone is the symmetric cyclic model (entropy 1.0). Adding rule 2 only shifts the balance between colours.
- **One untested point.** On a 64×64 lattice at this update rate, I'd expect patchy, noisy fronts more than clean spiral waves, but nothing in the evidence tests that.

**Twist assessment:** The only extra feature is the uneven type balance caused by the one-sided spontaneous flip. That is covered by the asymmetric-RPS and survival-of-the-weakest literature, so it doesn't justify KNOWN-VARIANT.

---
