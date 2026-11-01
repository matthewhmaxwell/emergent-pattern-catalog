# Stage 3 (b) review inputs: C_schelling_voter_majority (26 phenomena)

## census:g1:grid-6

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 2 types (random start). Rules, applied in order each step: (1) a type-0 cell becomes type 1 if at least 1 neighbours are type 1 (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  g_largest_domain, g_n_domains

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
These grid features are measured on the transient: the stretch of the run before the lattice collapsed to a single type.
  - g_largest_domain (one cell-type fills most of the lattice): median 1.0 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length
  - g_n_domains (many separate domains): median 0.6931 with interaction; with interaction knocked out the program never produces this view, at the original size, at 4x area or at 3x run length

### Stage 3 reply: verdict KNOWN, standard name: SI epidemic probabilistic cellular automaton / Eden growth (irreversible contact spreading to all-type-1 absorbing state)

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Single irreversible rule 0->1 with probability 0.1 if any of 8 neighbours is type 1, from a ~50% random start: fills the lattice within tens of steps. g_largest_domain=1.0 is the all-type-1 absorbing state. g_n_domains=0.6931=ln2 probably means log(1+1), i.e. one domain, so the 'many domains' label looks wrong (inference about metric scaling). The knock-out removes the only rule, leaving a frozen random field, so the comparison is uninformative. No phenomenological twist found; first-pass 'variant' downgraded to known and trivial.

Citations:
- https://link.springer.com/chapter/10.1007/978-3-319-65558-1_12
- https://www.researchgate.net/publication/323313442_Epidemic_Automaton_and_the_Eden_Model_Various_Aspects_of_Robustness
- https://arxiv.org/pdf/2010.09677v1
- https://en.wikipedia.org/wiki/Epidemic_models_on_lattices

Full reply text:

**Verdict: KNOWN.** This is the textbook SI ("Epidemic") probabilistic cellular automaton, the finite-grid version of Eden growth. What the checks measured is mostly the system settling into the all-type-1 state it can never leave. I'm overturning the first pass's "variant" call: I found no twist in the behaviour itself.

**What the system is**
- The rule is one-way and never reverses: a type-0 cell turns type 1 with probability 0.1 per step if any of its 8 neighbours is type 1. That is susceptible→infected spread with no recovery.
- This model is studied directly in the literature: "The two-dimensional probabilistic cellular automaton Epidemic models the spread of an epidemic without recovering on graph." The same survey covers properties of Epidemic on a finite grid and its analogous on the infinite square lattice: the Eden model.
- Lattice SI versions in cellular-automaton form are also standard (e.g. cellular automata implementation of the three basic epidemic models SI, SIS and SIR on a two-dimensional regular lattice).
- Using 8 neighbours instead of 4, synchronous updates and 64×64 wrap-around are implementation details. They don't change what the system does.

**Why the features look trivial**
- The start is random, so about half the cells are already type 1, spread everywhere. With no rule turning cells back to type 0, the type-0 pockets just fill in. In expectation each type-0 cell next to a type-1 cell converts within about 10 steps, so the lattice should be all type 1 within a few tens of steps. The simulation evidence confirms the run "collapsed to a single type."
- **g_largest_domain = 1.0** is that final all-type-1 state. It's the expected outcome of an irreversible SI process on a finite lattice, not a self-organised domain pattern.
- **g_n_domains median 0.6931** equals ln 2. That's what you'd get from log(1 + 1), i.e. one domain on a log(1+n) scale. That reading is my inference about how the metric is scaled, but if right, the "many separate domains" label is wrong: it's recording one domain.
- **The knock-out comparison means little.** Removing the only rule makes the system do nothing, so the frozen random start stays as it is. Any statistic will differ between "spread happened" and "nothing happens."

**Searches run:** one combined search on SI epidemic cellular automata, irreversible infection, absorbing states and Eden growth. It returned the Epidemic/Eden survey and SI lattice-model papers directly. I didn't do more searching because the model matches an existing one exactly, not loosely.

---

## census:novelty:field-3

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 2 field(s) with diffusion [0.02, 0.1]. Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.1 per step); (2) type-1 cells release 1.0 of field 0 per step.

OBSERVED (features that passed the simulation checks):
  f_std, f_n_blobs

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - f_std: median 289.899 with interaction, 33.8377 with interaction knocked out (difference 256.0613, robust spread 6.4747); difference at 4x area 273.6838, at 3x run length 693.0557
  - f_n_blobs: median 1.3863 with interaction, 3.1355 with interaction knocked out (difference -1.7492, robust spread 0.2325); difference at 4x area -2.2914, at 3x run length -1.4263

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: f_std 0.0124; f_n_blobs 1.2425
  - rules 2 only: f_std 33.8377; f_n_blobs 3.1355
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Majority-vote cellular automaton / zero-temperature Glauber (Ising) coarsening, read out by a passive, non-decaying emitted diffusive field

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Only interaction is majority-rule cell updating, i.e. known Ising-like coarsening into few large domains (lower f_n_blobs). Field 0 is emitted one-way, never read back, and never decays, so f_std grows without limit as it builds up over large domains (difference roughly triples at 3x run length): an inert field plus a build-up artifact. Field 1 and the agents are not used by any rule. End state tends to a frozen absorbing configuration. First pass's 'variant' label is overturned to KNOWN, with the features largely trivial.

Citations:
- https://arxiv.org/pdf/2604.10082
- https://www.sciencedirect.com/science/article/abs/pii/S0304397521004370
- https://link.springer.com/article/10.1007/s00220-018-3180-2
- https://arxiv.org/pdf/cond-mat/9701118

Full reply text:

**Verdict: KNOWN.** The behaviour is majority-vote (zero-temperature Ising) coarsening, read out through a field that only goes one way and never decays. I'd also class it as mostly trivial. I'd overturn the first pass's "variant" label: nothing about the behaviour itself is new.

**What the system actually does**
- **Rule 1** is the classic majority-vote cellular automaton, here updated at random rather than all at once. It can also be considered as the zero-temperature limit of the Ising Model. In two dimensions, the dynamics of the nearest-neighbors majority model, like that of the Ising model, is determined by the Cahn-Allen theory of coarsening due to the curvature of the front and a kind of surface tensions. On a torus with an 8-neighbour (Moore) neighbourhood, the majority model is known to end in final complete occupancy by white, stable coexistence of both colors, and complete occupancy by black, depending on the starting density. The random-update form is the standard coarsening model (zero-temperature Ising Glauber dynamics).
- **Rule 2:** type-1 cells add field 0 every step. The field diffuses, but nothing in the rules makes it decay, and nothing reads it back. So the field just copies the domain map onto a continuous grid, and the amount in it keeps growing.
- **Never used:** field 1 and the 400 agents. No rule reads them, so they play no part in the observed features.

**Why the two features are expected**
- **f_n_blobs (fewer blobs, about 1.4 vs 3.1):** this is just coarsening. Once random noise has merged into one or two large domains, the field (which mirrors them) shows one or two blobs. In the knock-out, the random static colouring keeps the field blotchy.
- **f_std (290 vs 34):**
  - **Partly real coarsening:** large uniform domains stop diffusion from evening out the field. With random cells, nearby sources blur together.
  - **Partly the missing decay:** the difference roughly triples at 3× run length (256 to 693). That fits a field that keeps building up over large domains without limit. It does not fit a steady pattern.
- **Rule-1-only result:** in that sub-program field 0 is never fed, so its standard deviation is about 0. Its blob count of about 1.24 is probably measuring an essentially empty field, so it means little.

**Trivial or artifact?**
- The only genuine emergent process is majority-vote coarsening, and it is textbook.
- The measured features come from an inert, one-way field. Their size is inflated by build-up with no decay, which is not a stable pattern.
- The end state heads towards a frozen, absorbing configuration: stuck stripes or domains, or one colour everywhere.

There is no feedback between field and cells and no motility coupling, so I found no twist in the phenomenology to support "variant". I ran one search, for majority-rule coarsening and zero-temperature Glauber dynamics on a Moore lattice. That was enough because the behaviour maps directly onto a standard model.

---

## census:novelty:grid-7

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 1 type (random start). Agents move at speed 1.0 with heading noise 1.0 rad per step. 2 field(s) with diffusion [0.02, 0.1]. Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.1 per step); (2) each agent advances the type of the cell under it by one (cyclically).

OBSERVED (features that passed the simulation checks):
  g_reversals, g_interface

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_reversals (cells flip back and forth): median 1.0 with interaction, 0.0 with interaction knocked out (difference 1.0, robust spread 0.02); difference at 4x area 1.0, at 3x run length 1.0
  - g_interface (long/rough domain boundaries): median 0.4315 with interaction, 0.498 with interaction knocked out (difference -0.0665, robust spread 0.0186); difference at 4x area -0.0682, at 3x run length -0.0676

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_reversals 1.0; g_interface 0.0272
  - rules 2 only: g_reversals 1.0; g_interface 0.4999
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Majority-vote model with noise (Oliveira 1992) / kinetic Ising coarsening interrupted by noise, disordered high-noise phase

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Agents flipping cells act as random-walking noise sources. They flip ~10% of cells per step while the majority rule updates ~10% per step, so effective noise is ~0.5, far above q_c≈0.075: a disordered phase. g_reversals is trivial (rules-2-only already gives 1.0). g_interface (0.43) sits just below random colouring (0.50), i.e. weak short-range clustering competing with noise. The chemical fields are inert. Agreeing with 'known' but sharpening it: this is not a meaningful variant, and the features are mostly trivial.

Citations:
- de Oliveira, M.J., Isotropic majority-vote model on a square lattice, J. Stat. Phys. 66, 273-281 (1992), doi:10.1007/BF01060069
- Kwak et al., Critical behavior of the majority voter model is independent of transition rates, PRE 75, 061110 (2007)

Full reply text:

**Verdict: KNOWN.** I agree with the first pass that this is a known model, but I'd go further: it isn't a meaningful variant. It is a majority-vote (Ising-type) lattice pushed so hard by noise that it stays disordered. The two features that passed the checks look mostly trivial rather than like a real emergent pattern.

**What the system actually is**
- Rule 1 (adopt the neighbourhood majority with probability 0.1) is a zero-noise majority-vote rule. With only two types, rule 2 ("advance the cell type by one, cyclically") just flips the cell. The agents turn at random by up to 1 rad per step, so they act as moving, randomly placed noise sources.
- Together that is the standard majority-vote model with noise. The two-state majority-vote model with noise, proposed by Oliveira in 1992, is one of the simplest nonequilibrium generalizations of the Ising model. It presents up-down symmetry and a continuous order-disorder phase transition at a critical value of noise. In that model each agent assumes the state of the majority of its neighbors with probability 1 −q, and assumes the opposite state with probability q.
- On the square lattice the critical exponents ν, γ, and β are the same as those of the Ising model and the critical noise parameter is q_c=0.075±0.001. Moving or spatially correlated noise sources don't change this picture: the system is still up-down symmetric and in the same universality class.

**Rough noise estimate**
- 400 agents on 4,096 cells flip about 10% of cells per step.
- Rule 1 updates about 10% of cells per step.
- So the effective noise is roughly 0.5, far above the critical value of 0.075. The model should sit deep in its disordered phase.
- That matches the numbers. The full program's g_interface (0.43) is close to random colouring (≈0.50, the agents-only value) and far from coarsened domains (0.027, the majority-rule-only value). That is the weak short-range clustering of a disordered Ising-like phase, not a new kind of boundary structure.

**Why I call the features trivial or artifacts**
- **g_reversals:** the "rules 2 only" sub-program already gives 1.0. Cells flip back and forth simply because agents keep flipping them. Knocking out the interaction removes the agent-to-cell writes, so the knock-out value of 0 tells you nothing about emergence.
- **g_interface:** the shift of -0.067 is real and survives at larger size and longer runs. But it is just the expected in-between value when coarsening competes with noise; nothing new appears.
- **Chemical fields:** no rule reads or writes them, so they are inert.
- **Agents:** they don't interact with each other and nothing feeds back into how they move. There is no collective motion.

**What I couldn't confirm:** I didn't find a paper on this exact setup, where the noise comes from wandering agents. Under your criteria that doesn't matter, because the behaviour (a disordered phase with short-range domains above the critical noise) is standard. My search was a single query on Oliveira's majority-vote model.

---

## qd:045

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) each cell copies the type of a random neighbour (with probability 0.1 per step); (2) a type-0 cell becomes type 2 where field 0 exceeds 0.5.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_n_types, g_reversals, g_type_entropy, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 5.9242 with interaction, 7.6811 with interaction knocked out (difference -1.7569, robust spread 0.2721); difference at 4x area -1.8094, at 3x run length -1.8836
  - g_corrlen (long spatial correlation length): median 3.0 with interaction, 1.0 with interaction knocked out (difference 2.0, robust spread 0.08); difference at 4x area 2.0, at 3x run length 2.0
  - g_n_types (number of surviving types): median 3.0 with interaction, 4.0 with interaction knocked out (difference -1.0, robust spread 0.14); difference at 4x area -1.0, at 3x run length -1.0
  - g_reversals (cells flip back and forth): median 0.7641 with interaction, 0.0 with interaction knocked out (difference 0.7641, robust spread 0.0153); difference at 4x area 0.7688, at 3x run length 0.7499
  - g_type_entropy (cell types evenly mixed): median 0.7429 with interaction, 0.9998 with interaction knocked out (difference -0.2568, robust spread 0.0349); difference at 4x area -0.26, at 3x run length -0.258
  - g_moran (like-valued cells cluster spatially): median 0.5183 with interaction, -0.0033 with interaction knocked out (difference 0.5217, robust spread 0.0144); difference at 4x area 0.5062, at 3x run length 0.5397

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 6.143; g_corrlen 3.0; g_n_types 4.0; g_reversals 0.6558; g_type_entropy 0.976; g_moran 0.5167
  - rules 2 only: g_n_domains 7.2284; g_corrlen 1.0; g_n_types 3.0; g_reversals 0.0; g_type_entropy 0.7531; g_moran -0.0029
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Multi-state voter model coarsening on a 2D lattice (logarithmically slow, pre-consensus transient), with one type removed by a one-time threshold conversion

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Rule 1 is the voter model and alone reproduces the clustering, correlation length, reversals and domain count. Rule 2 permanently converts all type-0 cells to type 2 (rule-2-only entropy 0.753 matches frequencies (0,0.25,0.5,0.25)); since the voter rule cannot recreate an extinct type, this alone explains the 3 surviving types and the lower entropy. No emitters are described, so the fields presumably only diffuse: field 1 is inert and field 0 is a passive threshold. The knock-out removes the voter rule itself, so its large differences reflect the voter model's presence, not an interaction between the rules. Overturns the first pass in that this is not a meaningful variant: the spatial pattern is genuine but textbook voter coarsening, and the coupling adds only a trivial recolouring.

Citations:
- Dornic, Chaté, Chave, Hinrichsen, Critical coarsening without surface tension: the universality class of the voter model, PRL 87, 045701 (2001)
- Percolation and coarsening in the bidimensional voter model, arXiv:1506.05321
- Ben-Naim, Frachebourg, Krapivsky, Coarsening and persistence in the voter model, PRE (cnls.lanl.gov/People/ebn/pubs/voter/voter-pub.pdf)
- Dall'Asta & Galla, Algebraic coarsening in voter models with intermediate states, arXiv:0806.0817

Full reply text:

**Verdict: KNOWN.** I overturn the first pass's "variant" label. What you observed is ordinary voter-model coarsening. The only effect of the chemical rule is to remove type 0 and fold its cells into type 2. It adds no new behaviour.

**What's driving each feature**
- **Rule 1 is the voter model.** Copying a random neighbour's type is the standard definition: spin-like variables align with the state of a randomly chosen nearest neighbor. On its own it already produces the clustering, correlation length, reversals and domain count you measured: Moran 0.517 against 0.518 for the full program, corrlen 3 against 3, reversals 0.66 against 0.76, domains 6.1 against 5.9.
- **The clustering is documented.** The voter model's evolution shows spatial clustering of similar opinions. It approaches one of the two absorbing states with complete consensus via a coarsening process in d≤2. In 2D that process is slow: the density of active interfaces vanishes only as 1/ln(t/τ). So at 64×64 you see a long-lived stage with several domains and short correlation lengths, which is what you measured.
- **Rule 2 just deletes type 0.** Running rule 2 alone gives 3 types and an entropy of 0.753. That matches type frequencies of about (0, 0.25, 0.5, 0.25), which have normalised entropy 1.5/2 = 0.75. In other words, every type-0 cell became type 2. Once type 0 is gone, rule 1 can never bring it back, because a cell can only copy a type that still exists. That one-time removal is the only reason the full program differs from rule 1 alone on type count (3 vs 4) and entropy (0.74 vs 0.98). The rule says nothing about anything producing either field, so they presumably just diffuse. If so, field 1 does nothing and field 0 matters only as a threshold that flips type 0 to type 2.
- **"Reversals" are standard voter-model behaviour.** Voters flipping back and forth while domains coarsen has been studied directly: Ben-Naim et al. track the fraction of voters who changed their opinion exactly n times during the time interval (0,t).

**The knock-out test is misleading here.** The knock-out removes rule 1 and keeps rule 2, so the large gaps you see on clustering, correlation length and reversals simply measure "voter model vs no voter model". They are not evidence that the two rules interact. Your own sub-program check makes the same point: rule 1 alone matches the full program on the spatial features. The "no sub-program came within 3 robust spreads" result happens only because rule 2 changes the type count, not because the combination produces anything new.

**Genuine or artifact?** The spatial clustering is a genuine, well-understood emergent pattern: voter coarsening on a finite lattice before it reaches consensus. The field coupling is effectively a trivial artifact: it acts as a one-off recolouring rather than a dynamical feedback. The result is a 3-state voter model with uneven starting fractions (0.25, 0.5, 0.25), and I found no new behaviour to account for.

---

## qd:046

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) each cell copies the type of a random neighbour (with probability 0.2 per step); (2) type-0 cells release 1.0 of field 1 per step.

OBSERVED (features that passed the simulation checks):
  f_std, f_n_blobs, f_peak_sharp

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - f_std: median 93.5807 with interaction, 30.8284 with interaction knocked out (difference 62.7523, robust spread 5.2733); difference at 4x area 101.275, at 3x run length 141.4264
  - f_n_blobs: median 1.6094 with interaction, 3.6506 with interaction knocked out (difference -2.0411, robust spread 0.2509); difference at 4x area -2.3348, at 3x run length -1.4966
  - f_peak_sharp (sharp field pattern): median 0.7451 with interaction, 0.3226 with interaction knocked out (difference 0.4224, robust spread 0.0606); difference at 4x area 0.2276, at 3x run length 0.3556

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: f_std 0.0065; f_n_blobs 1.0986; f_peak_sharp 0.6863
  - rules 2 only: f_std 30.8284; f_n_blobs 3.6506; f_peak_sharp 0.3226
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: 2D (multi-state) voter model coarsening with a passive diffusing tracer of one opinion

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Rule 1 is the neighbour-copying voter model (4 states, Moore neighbourhood), which coarsens slowly (active interfaces ~1/ln t) towards absorbing consensus. The field never feeds back. Field 1 is a non-decaying integral of the type-0 indicator, so f_std grows without bound with run length (62.8 -> 141.4 at 3x run length). Fewer blobs and sharper peaks just reflect voter domain coarsening. Field 0 has no source and is inert. The domains are genuine voter coarsening, but the field features are a passive accumulating readout heading for a consensus absorbing state, not a new pattern.

Citations:
- Dornic, Chate, Chave, Hinrichsen, PRL 87, 045701 (2001) - Critical coarsening without surface tension: the universality class of the voter model
- Tartaglia et al., Percolation and coarsening in the bidimensional voter model, arXiv:1506.05321
- Roy & Sen, Interplay of interfacial noise and curvature driven dynamics in two dimensions, arXiv:1701.09117

Full reply text:

**Verdict: KNOWN.** This overturns the first pass's "variant" call. The behaviour is the textbook 2D voter model, plus a diffusing field that just records where the voter domains are. Nothing in the field behaviour falls outside the documented voter phenomenology.

**What the rules amount to.** Rule 1 says each cell copies a random neighbour. That is the voter model: an agent follows the state of a randomly chosen neighbour at each time step. Using 4 types instead of 2 and 8 neighbours instead of 4 is a standard multi-state version and doesn't change the qualitative behaviour. The field never affects the cells, so there is no feedback in either direction. Rule 2 only lets type-0 cells "write" into field 1.

**Why the features are explained:**
- **Fewer blobs and sharper peaks.** These come from type-0 domains coarsening. The voter model's evolution shows spatial clustering of similar opinions. It approaches one of the two absorbing states with complete consensus via a coarsening process in d≤2. The coarsening is slow: the density of active interfaces vanishes only as 1/ln(t/τ). That fits the blob-count difference shrinking at 3x run length (−2.04 → −1.50), as domains merge towards consensus on a finite 64×64 torus. A higher peak-sharpness median with rule 1 alone (0.69) than with rule 2 alone (0.32) also points to the domain structure as the source.
- **Domain shape.** Voter domains are very jagged and preserve their fractal geometry even at the late stages of evolution, because phase ordering takes place without surface tension (Dornic et al., PRL 2001).
- **Large field spread (f_std).** Emission at a constant rate with no decay means field 1 piles up wherever type-0 domains last longest. Its spread therefore grows with how unevenly type 0 is spread out and with run time. That is why the difference climbs from 62.8 to 141.4 at 3x run length rather than settling. It is accumulation, not a stable emergent pattern. With rule 1 alone the spread is about 0.0065, essentially zero, because nothing emits. Field 0 has no source at all and is inert.

**Is it trivial or an artifact?** Partly. Voter coarsening is a genuine emergent effect: removing the copying rule changes the features, as the knock-out shows. But the "pattern" being measured is a passive, ever-growing smoothed picture of where the type-0 domains are. It has no feedback, no decay and no steady state, and the system is heading for an absorbing consensus state. So "interaction matters" is true only in the trivial sense that the copying rule creates the domains the field is drawn over.

---

## qd:050

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) a type-0 cell becomes type 3 where field 1 exceeds 0.05; (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_activity_trend, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 5.5946 with interaction, 7.6667 with interaction knocked out (difference -2.072, robust spread 0.2652); difference at 4x area -1.1906, at 3x run length -2.0738
  - g_activity_trend (activity rising): median -2.3317 with interaction, 0.0 with interaction knocked out (difference -2.3317, robust spread 0.2155); difference at 4x area -3.6791, at 3x run length -1.3904
  - g_reversals (cells flip back and forth): median 0.767 with interaction, 0.0 with interaction knocked out (difference 0.767, robust spread 0.0158); difference at 4x area 0.4686, at 3x run length 0.6482

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 7.636; g_activity_trend -1.1874; g_reversals 0.0
  - rules 2 only: g_n_domains 6.6201; g_activity_trend -3.1844; g_reversals 0.3099
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Schelling segregation model (multi-type, random relocation, Moore neighbourhood) with field-gated vacancy filling acting as a vacancy-density/composition shift

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Downgraded from KNOWN-VARIANT to KNOWN: no new behaviour. Rule 2 is the Schelling rule. Rule 1 is a one-way field threshold with no emission, so the field just smooths out and its only effect is filling vacancies with type-3 cells, i.e. a shift in vacancy density and type mix, both standard Schelling parameters. Fewer domains = segregation/coarsening. The 'activity rising' value is negative, so activity is falling, which is normal relaxation (feature mislabelled). Flip-backs = relocation churn of unsatisfied agents. The field part is nearly inert; the segregation domains are genuine.

Citations:
- https://arxiv.org/pdf/2608.16557
- https://arxiv.org/html/0711.2212v1
- https://arxiv.org/pdf/0903.4694
- https://www.nature.com/articles/s41598-022-23224-7
- https://arxiv.org/pdf/2509.14462

Full reply text:

**Verdict: KNOWN.** I agree with the first pass that this is Schelling's segregation model, but I'd go further: none of the observed behaviour goes beyond documented Schelling dynamics, so I'm putting it one step lower than "variant".

**What the system is.** Rule 2 is Schelling's rule almost exactly. An occupied cell checks the same-type share among its occupied 8 neighbours. If that share is below a threshold (here 0.3), the cell moves to a random empty site. The literature uses the same setup: an L×L square lattice with periodic boundary conditions, where all unsatisfied agents are identified, shuffled into a random order, and relocated sequentially to uniformly chosen empty sites. Here there are three agent types instead of two, which is a standard generalisation.

Rule 1 only reads a field, never other cells. No rule emits into field 1, so it just spreads out from its random start toward a uniform value. In effect it turns some empty cells into type-3 cells for a limited time, which lowers the vacancy density and changes the mix of types. Vacancy density and type mix are the main control parameters studied in the Schelling literature, so this is a parameter shift, not a new mechanism.

**Feature by feature:**
- **Many separate domains (fewer than in the knock-out):** This is the classic segregation result. Segregation is routinely measured by cluster count: in maximally segregated states, clusters grow larger, thus reducing the total number of clusters. Cluster count is also known to depend on vacancy density: the number of clusters in the final states, N_C, decreases as v decreases. That is exactly the quantity rule 1 shifts.
- **"Activity rising":** The measured value is negative (−2.33), so activity is actually *falling*. The feature label is misleading. Falling activity is ordinary Schelling relaxation. In one study, after 10 time steps, the proportion of satisfied agents is already very close to 1, after which the clusters keep coarsening, in a process strongly reminiscent of coarsening effects in alloys.
- **Cells flip back and forth:** At the site level, this is just unhappy agents moving in and out of cells, the ongoing churn that random relocation produces. Known Schelling regimes include frozen, segregated and mixed phases. Persistent turnover with a low threshold, three types and random relocation is expected.

**Sub-programs.** That the full program differs from "rule 2 only" is expected, because rule 1 changes how many empty cells and type-3 cells there are. I don't count that as a new interaction between the rules.

**Genuine or artifact?** The domains and their coarsening are a genuine emergent pattern: textbook Schelling segregation. The field part is close to inert. It doesn't respond to the cells and only acts as a one-off, decaying trigger for filling vacancies.

---

## qd:052

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). 2 field(s) with diffusion [0.02, 0.01]. Rules, applied in order each step: (1) field 1 grows autocatalytically by consuming field 0 (rate 1.0 x f x g^2); (2) a type-0 cell becomes type 1 where field 1 exceeds 0.8.

OBSERVED (features that passed the simulation checks):
  g_largest_domain, g_n_types

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_largest_domain (one cell-type fills most of the lattice): median 1.0 with interaction, 0.0405 with interaction knocked out (difference 0.9595, robust spread 0.0208); difference at 4x area 0.9829, at 3x run length 0.9595
  - g_n_types (number of surviving types): median 1.0 with interaction, 2.0 with interaction knocked out (difference -1.0, robust spread 0.06); difference at 4x area -1.0, at 3x run length -1.0

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_largest_domain 0.0405; g_n_types 2.0
  - rules 2 only: g_largest_domain 0.0405; g_n_types 2.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Cubic autocatalysis (A+2B->3B) reaction-diffusion front driving an irreversible threshold conversion to an absorbing single-type state

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Field rule is closed cubic autocatalysis (Gray-Scott without feed/kill): autocatalyst invades and saturates everywhere via documented travelling fronts. Cell rule is a one-way type-0 to type-1 flip above threshold, with no reverse process and no feedback onto fields, so the lattice converges to an all-type-1 absorbing state. largest_domain=1 and n_types=1 are the trivial endpoint; knock-out differences only show that the threshold is reached when the field can grow. Not a variant: no new phenomenology. The first pass overrated novelty.

Citations:
- https://royalsocietypublishing.org/rspa/article/430/1880/509/16458
- https://arxiv.org/pdf/chao-dyn/9601004
- https://link.springer.com/article/10.1007/s10910-026-01776-7
- https://royalsocietypublishing.org/doi/10.1098/rsta.1991.0098

Full reply text:

**Verdict: KNOWN, and the result is trivial. It's an absorbing state, not a genuine emergent pattern. This overturns the first pass, which rated it more novel than it is.**

**What the system actually is**

Rule 1 is the closed cubic autocatalysis reaction A + 2B → 3B (rate f·g²), with no feed and no decay. Without those terms, the Gray-Scott model loses everything that makes its patterns. The textbook form is ∂α/∂t = −αβ² + D_A Δα, ∂β/∂t = αβ² + D_B Δβ. In this setting a localized introduction of the autocatalyst B into a uniform background of substrate A can trigger a travelling reaction front that converts the A-rich state into a B-rich one. These fronts were worked out decades ago. For pure cubic autocatalysis the wave has, at large times, a constant asymptotic speed v0 (where v0 = 1/√2 in the appropriate dimensionless units). Unequal diffusion rates, as here (0.02 vs 0.01), are also covered (Billingham and Needham, Phil. Trans. R. Soc. A, 1991). The final state is uniform: field 1 ends up holding nearly all the mass, about f+g everywhere. If that total is above 0.8, the threshold is crossed everywhere.

Rule 2 is a one-way, irreversible switch from type 0 to type 1 wherever field 1 is above 0.8. Nothing ever turns a cell back into type 0, and the cells have no effect on the fields.

**Why the measured features are trivial**

- **Largest domain = 1.0 and surviving types = 1.** This is exactly what you'd expect when an irreversible conversion is driven by a field that rises above its threshold everywhere. The lattice ends in a frozen absorbing state (all type 1). There's no competition, no coexistence and no domain shapes to speak of. The only thing that varies is how fast the front sweeps across, and that's the documented front speed.
- **The knock-out and sub-program tests don't count as evidence of emergence.**
  - Rule 2 alone never fires because field 1 never reaches 0.8 without growing.
  - Rule 1 alone grows the field but never touches the cells.
  - So the "interaction" is just one field switching on a one-way flip. That passes the "needs both rules" test without any collective behaviour.
- **Size and run-length robustness is expected.** An absorbing state that covers the whole lattice will look the same at 4× the area and 3× the run length.
- **The one interesting thing in this setup is ruled out by its own parameters.** In 2D, a planar cubic-autocatalysis front is unstable to transverse perturbations when D_A is sufficiently larger than D_B. Here D_A/D_B = 2, which is probably too small for that instability, and the cell features couldn't detect it anyway.

**Searches run:** cubic autocatalysis travelling fronts (one search; the result was conclusive).

---

## qd:053

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.001 per step); (2) a type-0 cell becomes type 1 where field 0 exceeds 0.5.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_n_types, g_activity_trend, g_type_entropy

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 6.271 with interaction, 7.6716 with interaction knocked out (difference -1.4006, robust spread 0.2789); difference at 4x area -1.4264, at 3x run length -2.7229
  - g_n_types (number of surviving types): median 3.0 with interaction, 4.0 with interaction knocked out (difference -1.0, robust spread 0.14); difference at 4x area -1.0, at 3x run length -1.0
  - g_activity_trend (activity rising): median -0.6227 with interaction, 0.0 with interaction knocked out (difference -0.6227, robust spread 0.0125); difference at 4x area -0.5981, at 3x run length -1.2725
  - g_type_entropy (cell types evenly mixed): median 0.6487 with interaction, 0.9998 with interaction knocked out (difference -0.3511, robust spread 0.033); difference at 4x area -0.3465, at 3x run length -0.5449

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 6.6438; g_n_types 4.0; g_activity_trend -0.6347; g_type_entropy 0.999
  - rules 2 only: g_n_domains 7.2101; g_n_types 3.0; g_activity_trend 0.0; g_type_entropy 0.7511
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Zero-temperature q-state Potts / majority-vote coarsening (q=4) from a biased initial composition; type 0 removed once by a passive diffusing field that does not respond to the cells

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Rule 1 alone gives the coarsening, the falling activity and the domain counts. Rule 2 alone gives 3 surviving types and entropy of about 0.75 (simple arithmetic: type 0 gone, type 1 about 0.5, types 2 and 3 about 0.25 each). The field has no source and no feedback, so it acts only as a one-time bias on the starting mix. The full program's lower entropy and fewer domains are just coarsening enlarging the type that already has the most cells. The first pass's 'variant' label overstates this: there is no new twist in the behaviour. The knock-out counts the passive field as an interaction, which inflates the interaction effect. The coarsening is genuine but textbook; the drop in types and entropy is an artifact of the passive field.

Citations:
- http://physics.bu.edu/~redner/pubs/pdf/potts2d.pdf
- https://arxiv.org/html/cond-mat/9503163
- https://arxiv.org/html/cond-mat/9705155

Full reply text:

I think the first pass was too generous, though not completely wrong. My verdict is **KNOWN**: there is no "variant" here. The behaviour is two well-known things that don't interact in any real way. One is ordinary majority-rule coarsening. The other is a one-time removal of type 0, driven by a field that never changes in response to the cells, which then tilts the type mix. I ran one search (zero-temperature Potts coarsening). It confirmed the coarsening references, and I didn't need more for this verdict.

**What each part does**
- **Rule 1 on its own is the standard q-state Potts model at zero temperature (majority voting), here with q = 4.** In 2D, domains grow from a random start, driven by surface tension at their boundaries, and an intricate coarsening domain mosaic emerges shortly after the quench. This mosaic has a characteristic length scale that grows in time as t^1/2. Studies of this model use uncorrelated initial conditions where each of the q phases is present with equal density 1/q, which is the same random start used here. For q ≥ 3, frozen and "blinker" final states are also documented: a subset of spins can flip ad infinitum. This rule alone explains the falling activity (−0.63 alone vs −0.62 in the full program) and the domain counts.
- **Rule 2 does not depend on the cells at all.** The field has no source and nothing feeds back into it. It is just a random starting field smoothing out by diffusion. Wherever it starts above 0.5, which is roughly half the grid, type-0 cells switch to type 1. The arithmetic matches: type 0 gone, type 1 at about 0.5, types 2 and 3 at about 0.25 each gives a normalised entropy of 1.5/2 = 0.75. Rule 2 alone measured 0.751. It alone also produces the surviving-type count of 3.
- **The combination:** coarsening that starts from a lopsided mix enlarges the type that is already in the majority. That explains why the full program's entropy is lower (0.65) and it has fewer domains than either rule alone. In the Potts literature an uneven starting mix has a known equivalent: a kinetic Ising model evolving with a fixed magnetization, meaning one state starts out with more cells than the others. Uneven starting conditions are standard there. So "no sub-program reproduces the full program" is just the two effects adding up. It is not a new behaviour.

**Real pattern or artifact?** Partly an artifact. The coarsening is real, but it is textbook behaviour. The drop from 4 types to 3 and the lower entropy come from a field that doesn't respond to the cells, acting as a one-time bias. The "knock-out" comparison also counts rule 2 as an interaction, even though the field takes nothing from the cells. That inflates how much the interaction seems to matter. The true two-way interaction between the field and the cells is zero.

---

## qd:056

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 3 types (random start). Agents have 3 types (random start). Agents move at speed 0.1 with heading noise 0.1 rad per step. Rules, applied in order each step: (1) each agent turns to the average heading of agents within distance 1.0; (2) each cell adopts its neighbourhood's majority type (every step); (3) each agent advances the type of the cell under it by one (cyclically).

OBSERVED (features that passed the simulation checks):
  g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_reversals (cells flip back and forth): median 0.5807 with interaction, 0.0 with interaction knocked out (difference 0.5807, robust spread 0.0116); difference at 4x area 0.5789, at 3x run length 0.71

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 3 rules):
  - rules 1 only: g_reversals 0.0
  - rules 2 only: g_reversals 0.8523
  - rules 3 only: g_reversals 0.3575
  - rules 1, 2 only: g_reversals 0.8523
  - rules 1, 3 only: g_reversals 0.3274
  - rules 2, 3 only: g_reversals 0.8042
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Period-2 oscillation of the synchronous majority-rule cellular automaton (Goles–Olivos / Poljak–Sůra theorem), plus agents perturbing cells and the majority rule restoring them

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: The majority rule alone gives a flip score of 0.85, higher than the full program's 0.58, so the feature comes from the known period-2 cycles of updating every cell at once. Agents bumping cell types (rule 3 alone: 0.36) add flips and break up the cycles; heading alignment (rule 1) has no effect. The knock-out removes the majority rule, so its zero score is guaranteed and the 0.58 gap says nothing new. The first pass's VARIANT verdict is downgraded to KNOWN because there is no twist in the behaviour itself, and the feature is largely an artifact of synchronous updating.

Citations:
- Goles & Olivos 1980, Periodic behaviour of generalized threshold functions, Discrete Math 30:187-189
- Poljak & Sůra (period <=2 for majority dynamics)
- arXiv:1711.10920 (Biased) Majority Rule Cellular Automata
- arXiv:2610.04311 The Two-Dimensional Majority Rule is P-Complete
- Frischknecht, Keller, Wattenhofer 2013, Convergence in (Social) Influence Networks

Full reply text:

**Verdict: KNOWN.** The flipping is well-documented behaviour, and most of it is a side effect of every cell updating at the same moment rather than a new emergent pattern. This mildly overturns the first pass: there's no real twist here that would justify calling it a variant.

**Where the flipping comes from**
- **The majority rule alone already does it.** With only rule 2 running, the flip score is 0.85 — higher than the full program's 0.58. Rule 2 updates every cell at once, and that is known to produce cells that flip every step. Goles and Olivos (1980), and independently Poljak and Sůra, proved that a large class of majority-based models, including (biased) majority model, always reach a cycle of period one or two. Their result covers multi-colour versions too: they consider settings where every v∈V has an initial color from set {0,1,··· , l}, so it applies to 3 types. In practice, after a certain amount of synchronous rounds, each participant has either a fixed opinion or changes its mind in every round. A 2D torus is covered explicitly: every trajectory on a finite torus eventually reaches a fixed point or a cycle of length two.
- **The agents add flips of their own, and also break up the cycles.** With only rule 3 (agents bumping the cell under them), the score is 0.36. That's simply agents writing to cells, and some cells switching back as agents move on. The full program's 0.58 lies between the two single-rule values. That fits two known effects mixing: agents knock cells out of their 2-cycles, and the majority rule repairs those cells. Rule 1 (agents aligning headings) changes almost nothing: 0.33 with rules 1 and 3, versus 0.36 with rule 3 alone.
- **The knock-out comparison overstates the result.** The knock-out removes the majority rule and agent-to-cell writing, so its score of 0 is guaranteed. The 0.58 gap is therefore not evidence of new collective behaviour. The full program also doesn't match any sub-program, but that's expected when two known mechanisms are combined.

**Trivial or artifact?** Yes, largely. The core is the period-2 artifact of updating all cells simultaneously; updating cells one at a time would remove most of it. The rest is agents perturbing cells and the majority rule pulling them back.

**Search:** I ran one targeted search, on synchronous majority-rule cellular automata and period-2 cycles. It turned up several independent sources confirming the theorem. I didn't search further for agents on a lattice combined with majority voting, because the sub-program results already show where the flipping comes from.

---

## qd:057

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields). Cells have 3 types (random start). Agents have 1 type (random start). Agents move at speed 0.5 with heading noise 0.0 rad per step. Rules, applied in order each step: (1) each agent turns away from the centre of any type agents within distance 1.0; (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.4 moves to a random empty cell; (3) each agent advances the type of the cell under it by one (cyclically).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_reversals, g_largest_domain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 6.3044 with interaction, 7.3569 with interaction knocked out (difference -1.0525, robust spread 0.2732); difference at 4x area -1.0676, at 3x run length -1.0644
  - g_reversals (cells flip back and forth): median 0.7029 with interaction, 0.0 with interaction knocked out (difference 0.7029, robust spread 0.0141); difference at 4x area 0.7068, at 3x run length 0.7418
  - g_largest_domain (one cell-type fills most of the lattice): median 0.3966 with interaction, 0.0077 with interaction knocked out (difference 0.3889, robust spread 0.0081); difference at 4x area 0.3625, at 3x run length 0.4049

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 3 rules):
  - rules 1 only: g_n_domains 7.3569; g_reversals 0.0; g_largest_domain 0.0077
  - rules 2 only: g_n_domains 6.6005; g_reversals 0.3333; g_largest_domain 0.0326
  - rules 3 only: g_n_domains 7.3418; g_reversals 0.3081; g_largest_domain 0.0071
  - rules 1, 2 only: g_n_domains 6.6005; g_reversals 0.3333; g_largest_domain 0.0326
  - rules 1, 3 only: g_n_domains 7.3576; g_reversals 0.47; g_largest_domain 0.0073
  - rules 2, 3 only: g_n_domains 5.7823; g_reversals 0.7992; g_largest_domain 0.4655
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Noise-driven coarsening in the Schelling segregation model (noisy / Glauber-Kawasaki Schelling; the zero-noise model freezes, small noise gives large domains)

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: Agent repulsion (rule 1) is inert: rules 1+2 equal rule 2, and rules 2+3 alone reproduce or exceed the effect. The mechanism is Schelling relocation (rule 2) plus non-conservative type-cycling noise from walkers (rule 3), which removes the zero-noise freezing and lets domains coarsen; this is documented. The noise source (spatially correlated, one-directional cycling) differs in the rules but no new phenomenology was identified. g_reversals is largely a trivial consequence of walkers rewriting cells. Check whether the large domain is the empty type (type 0); if so it may be a vacancy artifact rather than segregation.

Citations:
- Dall'Asta, Castellano, Marsili, J. Stat. Mech. (2008) L07002, arXiv:0707.1681
- Stauffer & Solomon / Schelling-Ising model, arXiv:0710.5237
- Barmpalias et al., Unperturbed Schelling segregation in two or three dimensions, arXiv:1504.03809
- Hydrodynamic limit of the Schelling model with spontaneous Glauber and Kawasaki dynamics, arXiv:2302.09866
- Aging effects in Schelling segregation model (noisy constrained Schelling, power-law coarsening)

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass's "variant" call in the direction of "known": the behaviour is Schelling segregation made to coarsen by noise. Only my last point, about which cell type forms the big domain, could change that.

**What the sub-programs show**
- Rule 1, agents steering away from each other, does almost nothing. Rules 1+2 give exactly the same numbers as rule 2 alone. Rules 2+3 without rule 1 already reproduce the effect, and slightly more strongly (largest domain 0.47, versus 0.40 for the full program). So agent–agent interaction is not the mechanism.
- Rule 2 alone is the plain Schelling model with no noise. It freezes into small domains: largest domain 0.03, reversals pinned at 0.333.
- Rule 3 adds what is effectively random noise. The 400 wandering agents step cells to the next type, which converts them, creates occupants and removes them. With that noise, the largest domain grows to about 40–47% of the lattice.

**Why this is documented**
- Dall'Asta, Castellano and Marsili treat the Schelling model with statistical physics and show it has nontrivial coarsening like in driven-particle systems and percolation-related phenomena.
- Most directly, the Schelling–Ising paper finds that for zero noise the Schelling and Ising models do not lead to large domains... Thus only small but nonzero noise can produce segregation into ghettos. That is the same contrast as rule 2 alone versus rules 2+3.
- Noise that changes a cell's type, rather than moving the occupant, is also studied. In Glauber-type Schelling models a single randomly selected unhappy agent switches type. Combined Glauber-plus-Kawasaki Schelling dynamics with random perturbations have been analysed too (arXiv 2302.09866).
- A noisy Schelling variant with aging likewise shows segregated states with a slow power-law coarsening process.

The noise here is spatially correlated, carried by moving agents, and pushes types in one direction around the cycle. That is a difference in the rules. I found nothing showing it changes the behaviour: the domains look like ordinary noise-driven segregation.

**Is it a real pattern or an artifact?**
- **Reversals:** largely trivial. Every agent rewrites the cell under it, so cells will flip back and forth whatever else happens. Rule 3 alone already gives 0.31, and rule 2 alone sits at about 1/3.
- **Number of domains:** the drop from 7.36 to 6.30 is small. It is just the other side of coarsening.
- **Largest domain:** this looks like genuine emergent coarsening. It is robust at 4× the area and 3× the run length.
- **Caveat:** please check which cell type forms the large domain. If it is type 0 (empty), the "domain" could be a region of vacancies the agents keep clearing, not segregation of occupants. That would make it closer to an artifact.

---

## qd:064

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus 400 self-propelled agents moving in a wrap-around square (20x20 when agents are alone, 64x64 when they share the lattice or fields); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 2 types (random start). Agents have 3 types (random start). Agents move at speed 0.3 with heading noise 0.3 rad per step. 1 field(s) with diffusion [0.05]. Rules, applied in order each step: (1) a type-1 cell becomes type 0 if at least 4 neighbours are type 0 (every step); (2) type-0 cells release 1.0 of field 0 per step.

OBSERVED (features that passed the simulation checks):
  f_mean, f_change, f_n_blobs, f_std

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - f_mean: median 1000.7294 with interaction, 501.1056 with interaction knocked out (difference 499.6237, robust spread 30.0367); difference at 4x area 500.5152, at 3x run length 1499.3796
  - f_change: median 10.0 with interaction, 5.0012 with interaction knocked out (difference 4.9988, robust spread 0.3); difference at 4x area 5.0076, at 3x run length 14.9963
  - f_n_blobs: median 0.6931 with interaction, 2.3979 with interaction knocked out (difference -1.7047, robust spread 0.1352); difference at 4x area -2.0571, at 3x run length -1.0075
  - f_std: median 0.035 with interaction, 21.7242 with interaction knocked out (difference -21.6892, robust spread 0.4352); difference at 4x area -24.2093, at 3x run length -32.38

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: f_mean 0.9833; f_change 0.0001; f_n_blobs 1.4979; f_std 0.0083
  - rules 2 only: f_mean 501.1056; f_change 5.0012; f_n_blobs 2.3979; f_std 21.7242
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Irreversible threshold growth (threshold-4 bootstrap percolation, 8-neighbour lattice) reaching the all-type-0 absorbing state, with a field that only accumulates (no decay) and grows linearly

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Rule 1 can only turn cells type 0, never back. From a random 50/50 start the 64x64 grid almost certainly fills (rule-1-only f_mean about 0.98). Every cell then emits 1.0 per step with no decay, so the field is uniform (f_std about 0) and grows linearly: f_mean about 1000, roughly double the knock-out's 500 because all cells emit instead of about half; the gap triples at 3x run length. f_n_blobs = ln 2 means a single grid-wide blob. The field never acts back on the cells and the agents take part in no rule. The features are bookkeeping of a frozen absorbing state, not an emergent pattern.

Citations:
- https://arxiv.org/pdf/2501.18976
- https://link.springer.com/article/10.1007/s00493-026-00213-2
- https://arxiv.org/pdf/2605.01151

Full reply text:

**Verdict: KNOWN.** This overturns the first pass's "variant" call. The cell rule is irreversible threshold growth, a standard bootstrap percolation model, and it fills the lattice. After that, the field just piles up evenly forever. The four features are artifacts of that, not an emergent pattern.

**What is actually happening**

- **Rule 1 is bootstrap percolation.** Type 0 is the "infected" state and it never reverts. A type-1 cell turns type 0 once at least 4 of its 8 neighbours are type 0. That is threshold-4 bootstrap percolation on the 8-neighbour square lattice, a family studied in its own right. Recent work proves results for a wide family of threshold bootstrap percolation rules on the 2D square lattice, including neighbourhoods given by large ℓ^p balls (p up to ∞). The ∞ case is the 8-neighbour square. That paper also discusses the 4-neighbour model on the square-diagonal lattice (the ⊠ model). The standard rule is: an inactive site becomes active if at least k of its neighbors are active, while active sites remain active. Starting from a random half-and-half mix on a 64×64 wrap-around grid, the type-0 share is far above what is needed to take over. So the grid almost certainly locks into the all-type-0 absorbing state within a few steps.
- **Rule 2 only adds up emissions.** The field has no decay, so once every cell is type 0, each cell adds 1.0 per step. The field becomes uniform and rises in a straight line.

**How the numbers follow from this**

- **f_mean:** about 1000 is roughly 4096 cells × 1.0 per step × the number of steps, spread over the grid. The knock-out gives about 500 because only about half the cells are type 0 there. It also scales exactly with run length: the gap is 499.6 at the normal length and 1499.4 at 3× length.
- **f_change:** the field rises by a fixed 10 versus 5. Those are just the emission rates for all cells versus half of them.
- **f_std:** about 0.035 means the field is flat, because every cell emits the same amount. The knock-out's std of 21.7 comes from its frozen random emitter pattern, so it is the knock-out that has spatial structure, not the full system.
- **f_n_blobs:** 0.693 = ln 2, which presumably means "one blob" (the whole grid). The knock-out's value of 2.398 = ln 11 likely reflects scattered patches.
- **Sub-program check:** with rule 1 alone, f_mean is 0.98, which looks like the fraction of cells that end up type 0. That confirms the grid fills almost completely.

**Is it genuine or an artifact?** Artifact: a frozen absorbing state plus an unbounded field that only accumulates. The difference from the knock-out is real and robust. But it comes entirely from rule 1 turning half the grid into emitters, and rule 2 just records that. There is no feedback from the field to the cells, and the agents play no part in any rule. The "passes at 4× area and 3× run length" result only reflects linear growth in time, not a stable pattern.

---

## qd:071

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each); plus continuous chemical field(s) diffusing on a 64x64 wrap-around grid. Cells have 4 types (random start). 2 field(s) with diffusion [0.05, 0.02]. Rules, applied in order each step: (1) field 1 grows autocatalytically by consuming field 0 (rate 0.5 x f x g^2); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.2 moves to a random empty cell; (3) type-2 cells release 1.0 of field 1 per step.

OBSERVED (features that passed the simulation checks):
  g_activity_trend, g_interface

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_activity_trend (activity rising): median -1.5028 with interaction, 0.0 with interaction knocked out (difference -1.5028, robust spread 0.0301); difference at 4x area -1.5348, at 3x run length -1.5028
  - g_interface (long/rough domain boundaries): median 0.6146 with interaction, 0.7487 with interaction knocked out (difference -0.1341, robust spread 0.0273); difference at 4x area -0.1302, at 3x run length -0.1341

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 3 rules):
  - rules 1 only: g_activity_trend 0.0; g_interface 0.7487
  - rules 2 only: g_activity_trend -1.0829; g_interface 0.6031
  - rules 3 only: g_activity_trend 0.0; g_interface 0.7487
  - rules 1, 2 only: g_activity_trend -1.0829; g_interface 0.6031
  - rules 1, 3 only: g_activity_trend 0.0; g_interface 0.7487
  - rules 2, 3 only: g_activity_trend -1.0829; g_interface 0.6031
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Schelling segregation model (multi-type, Moore neighbourhood, vacancies) with coarsening / decreasing interface length and relaxing move activity

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Both features come only from rule 2 (Schelling relocation with a 0.2 threshold); sub-programs without rule 2 equal the knock-out exactly. The fields are one-way: no rule reads them, so the reaction-diffusion part is inert with respect to the cells. The gap between the full program and rule 2 alone is probably caused by the activity measure counting field changes, or by RNG-stream effects, not by real coupling. Calling it a 'variant' gives the inert field too much credit; the behaviour is plain Schelling coarsening.

Citations:
- https://www.demographic-research.org/articles/volume/21/12
- https://www.nature.com/articles/s41598-022-23224-7
- https://arxiv.org/pdf/1505.00429

Full reply text:

**Verdict: KNOWN.** This overturns the first pass. Both observed features come from Schelling segregation (rule 2), and the chemical part of the system is inert.

**What drives the behaviour.** Rule 2 is the Schelling model: agents of 3 types plus empty sites on a wrap-around grid with 8 neighbours each. A cell moves to a random empty site when fewer than 20% of its occupied neighbours share its type. No rule makes a cell read either field. Rules 1 and 3 change the fields but never feed back into what the cells do, so the fields are passive passengers.

The sub-program results agree:
- Every program that includes rule 2 shows falling activity (-1.08) and shorter boundaries (0.603).
- Every program without rule 2 matches the knock-out exactly (0.0 and 0.7487).

**How the observed features match the literature.**
- **Boundaries shortening (g_interface below the knock-out):** This is the classic Schelling result. It has been proved that as the Schelling model evolves, the total perimeter of the interface between the different agents decreases. Physics studies describe it as coarsening: the interface density ⟨ρ(t)⟩ decreases as ⟨ρ(t)⟩∼t^-α so the domains should grow in our model following a power-law with time. They use the same kind of grid, a L×L Moore lattice with 8 neighbors per site and periodic boundary conditions.
- **Activity falling (negative g_activity_trend):** Once the domains form, unhappy agents settle and movement winds down. In the coarsening regime, the size of the interface monotonically decreases in time as t^(-1/z)... This is the typical signature of a coarsening phenomenon, that has already been observed in the classical Schelling model.

**The gap between the full program and rule 2 alone.** The full program differs from rule 2 on its own (-1.50 vs -1.08 for activity, 0.615 vs 0.603 for boundaries). This does not look like a real coupling, because no rule lets the fields affect the cells. The likely explanations are:
- the activity measure also counts changes in the fields (field 1 is fed by type-2 cells and grows autocatalytically), or
- rules 1 and 3 use random numbers and shift the random sequence that rule 2 draws from.

So the "interaction" is just Schelling segregation with an inert field attached. A useful check would be to re-measure activity using only cell moves: I expect the full program to then match rule 2 alone.

**Is it genuine or an artifact?** The segregation and coarsening are a real emergent pattern, but a textbook one. The reaction-diffusion part (field 1 growing by consuming field 0, plus emission) is a trivial add-on: it has no effect on the cells, so it is not a coupled pattern.

---

## review:011

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
N kN=2 | N[0] ALWAYS() -> REWIRE_SAME() p=1.0  (compact grammar notation)

OBSERVED (features that passed the simulation checks):
  n_components

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - n_components (number of connected components): median 2.9174 with interaction, 1.4979 with interaction knocked out (difference 1.4195, robust spread 0.3233); difference at 4x area 1.5063, at 3x run length 1.4195

### Stage 3 reply: verdict KNOWN, standard name: Rewiring-only (phi=1) limit of the coevolving/adaptive voter model; homophilic network segregation / graph fission

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Two fixed node kinds plus always-rewire-to-same-kind is the opinion-frozen, rewiring-only limit of the Holme-Newman/adaptive voter model. Cross-kind links are removed and never re-created, so the graph must split into at least one component per kind (2) plus small leftover pieces, matching the measured ~2.9 vs 1.5. The split is real and robust, but it follows almost directly from the rule and can't be undone (a frozen segregated state even though links keep moving), so it is not a meaningful emergent pattern. First pass's VARIANT label is downgraded to KNOWN because no twist in the behaviour was found.

Citations:
- Holme & Newman, Phys. Rev. E 74, 056108 (2006)
- Henry, Pralat & Zhang, PNAS 108, 8605 (2011)
- Durrett et al., Graph fission in an evolving voter model, PNAS 109, 3682 (2012)
- Coevolving nonlinear voter model phase diagram incl. absorbing fragmented phase

Full reply text:

**Verdict: KNOWN.** This is the "rewiring-only" limit of the coevolving or adaptive voter model, also known as homophilic network segregation. The result itself is real, but it follows almost directly from the rule, so it isn't a meaningful emergent pattern. I'm overturning the first pass's "variant" label.

**How I read the rule.** There are N nodes of two kinds. Every node, at every step, cuts one of its links and reconnects to a node of its own kind. Node kinds never change. Without the rule, the graph stays random and static.

**Where the behaviour is documented**
- **Holme & Newman (2006).** Their model has one parameter that sets the balance between nodes copying a neighbour's opinion and nodes rewiring links to people who share their opinion. As that parameter varies, the model undergoes a continuous phase transition from a regime in which opinions are arbitrarily diverse to one in which most individuals hold the same opinion. At the rewiring-only end, opinions never change and the network breaks into separate components, each holding a single opinion. That is this system with two kinds. Note: the search results only show the abstract. The rewiring-only limit is from my knowledge of the paper, not something I checked in these sources.
- **Later work on fragmentation.** Other adaptive-network models show the same split: a phase diagram with three phases: dynamically active coexistence in a single-component network, absorbing consensus in a single-component network, and an absorbing phase in a fragmented network. Durrett et al., "Graph fission in an evolving voter model" (PNAS 2012), studies the same splitting; I only saw it cited, not the paper itself.
- **Henry, Prałat & Zhang (PNAS 2011).** Their actors have a small aversion from being connected to others who are dissimilar to themselves, and yet no actor strictly prefers a segregated network. They show that segregated networks always emerge regardless of the level of aversion. The rule here is the extreme case of that: every rewiring step is homophilic.

**Why it's close to trivial.** Rewiring keeps adding same-kind links and never adds a cross-kind link. Once the original cross-kind links have been rewired away, the graph has to split into at least one component per kind (2), plus whatever small pieces or isolated nodes are left. That matches what you measured: about 2.9 components versus about 1.5 for the static random graph. The links keep moving, but the split by kind can't be undone. No twist in the behaviour stands out: the count holds steady when you scale the area or run length, which is what a split by kind should do.

**Searches I ran:**
1. Holme–Newman coevolution and fragmentation.
2. Henry–Prałat–Zhang segregation in evolving networks.

---

## review:018

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) each node copies the type of a random neighbour (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (every step).

OBSERVED (features that passed the simulation checks):
  n_edge_turnover, n_type_modularity, n_giant

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - n_edge_turnover (edges rewire): median 0.921 with interaction, 0.0 with interaction knocked out (difference 0.921, robust spread 0.0184); difference at 4x area 0.9198, at 3x run length 0.921
  - n_type_modularity (network splits into like-type modules): median 0.6443 with interaction, -0.0062 with interaction knocked out (difference 0.6504, robust spread 0.0176); difference at 4x area 0.658, at 3x run length 0.6504
  - n_giant (fraction in the giant component): median 0.4025 with interaction, 0.985 with interaction knocked out (difference -0.5825, robust spread 0.0278); difference at 4x area -0.5819, at 3x run length -0.5825

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: n_edge_turnover 0.0; n_type_modularity 0.0; n_giant 0.985
  - rules 2 only: n_edge_turnover 0.6675; n_type_modularity 0.6625; n_giant 0.3775
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Adaptive/coevolving voter model (multistate), fragmented phase / fragmentation transition

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Both rules fire every step, so rewiring dominates and the system lands deep in the known fragmented phase: like-type, internally uniform components, giant component ~1/3 for 3 types. Rewiring alone (rule 2) already gives nearly the same modularity and giant fraction, so most of the structure is homophilic rewiring on random types. Copying mainly raises turnover. The late state is probably a frozen absorbing fragmented configuration, and turnover may be a cumulative count. Genuine but standard behaviour, not an artifact. Tightens 'variant' to 'known'.

Citations:
- Holme & Newman, Phys. Rev. E 74, 056108 (2006), arXiv:physics/0603023
- Böhme & Gross, Phys. Rev. E 85, 066117 (2012), Fragmentation transitions in multistate voter models
- Böhme & Gross, Phys. Rev. E 83, 035101 (2011), arXiv:1012.1213
- Vazquez, Eguíluz & San Miguel, PRL 100, 108702 (2008)

Full reply text:

**Verdict: KNOWN.** This is the adaptive (coevolving) voter model with 3 opinion states, running in its well-documented fragmented phase. The first pass found the right family. I'd tighten "variant" to "known", because nothing in the observed behaviour goes beyond what's already published.

**Why it's the same system**
- **Same two rules.** Holme & Newman's model combines exactly these two processes: copying a neighbour's opinion, and cutting links to disagreeing neighbours in favour of like-minded nodes. One parameter sets the balance between them. The model undergoes a continuous phase transition as this parameter is varied, from a regime in which opinions are arbitrarily diverse to one in which most individuals hold the same opinion. Their version also allows more than two opinions: a network with N nodes ... and G different opinions.
- **Three states specifically.** Böhme & Gross (2012) studied the multistate version, in which the underlying social network breaks into disconnected components. In particular, they analytically studied the three-opinion case.
- **The observed outcome is the documented one.** In the fragmented phase, both opinions survive in disconnected network components that are internally in consensus.

**Matching each feature**
- **Type modularity ≈ 0.64:** like-type, internally uniform components. That is the fragmented phase.
- **Giant component ≈ 0.40:** close to 1/3 plus fluctuations, which is what you'd expect if the network splits into roughly one block per type. That's the fragmentation transition.
- **Edge turnover:** rewiring is one of the model's two rules, so it's built in. The interaction adds turnover (0.92 vs 0.67 with rewiring alone), most likely because copying keeps creating new mismatched links before things settle. That's expected and not new.
- **Your run's settings:** both rules fire every step, so rewiring is strong. That puts the model on the fragmentation side of the transition, which is the regime you observed.

**Is it genuine or trivial?** It's real emergent behaviour, but a weak case:
- **Rewiring alone gets you most of the way.** Rule 2 by itself gives modularity 0.66 and giant 0.38, about the same as the full program. The type-segregated network mostly comes from like-with-like rewiring on fixed random types. Copying mainly changes turnover and drives which component structure the network freezes into.
- **The end state is probably frozen.** Once components are uniform, nothing more can happen, so the late-run averages likely describe a frozen absorbing state, not ongoing dynamics.
- **Turnover may be cumulative.** If the 0.92 counts edges ever rewired, not current rewiring, it says little about dynamics.

None of this makes it an artifact. It's the standard frozen fragmented state, not a numerical failure or a random colouring.

**Twist check:** I found no phenomenology beyond the documented fragmentation and like-type segregation. Synchronous updates and always-on rewiring only choose where on the known phase diagram the run sits.

---

## review:019

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a random network of 200 nodes (average 4 links each). Nodes have 3 types (random start). Rules, applied in order each step: (1) each node copies the type of a random neighbour (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (with probability 0.1 per step).

OBSERVED (features that passed the simulation checks):
  n_edge_turnover, n_consensus_gain

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - n_edge_turnover (edges rewire): median 0.9695 with interaction, 0.0 with interaction knocked out (difference 0.9695, robust spread 0.0194); difference at 4x area 0.9938, at 3x run length 0.9695
  - n_consensus_gain (nodes converge to agreement): median 0.635 with interaction, 0.0 with interaction knocked out (difference 0.635, robust spread 0.0127); difference at 4x area 0.6444, at 3x run length 0.635

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: n_edge_turnover 0.0; n_consensus_gain 0.625
  - rules 2 only: n_edge_turnover 0.6596; n_consensus_gain 0.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Coevolving / adaptive voter model (Holme-Newman rewire-to-same), 3-opinion version, consensus side of the fragmentation transition

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Low rewiring rate (0.1 vs copying every step) puts the system below the fragmentation threshold, so it reaches global consensus. The consensus gain comes almost entirely from rule 1 alone (0.625 vs 0.635), i.e. plain voter-model absorption; edge turnover follows directly from rule 2. No synergistic emergent feature. Synchronous updating and separate rewiring step are implementation details that don't change the known phenomenology. Not an artifact, but low novelty; the final state is the documented absorbing consensus state.

Citations:
- Holme & Newman, Phys. Rev. E 74, 056108 (2006)
- Durrett et al., PNAS 2012, Graph fission in an evolving voter model
- Shi, Mucha, Durrett, Multi-opinion coevolving voter model with infinitely many phase transitions (PMC5131864)
- Böhme & Gross 2011 (three-opinion adaptive voter model, as cited in PMC5131864)
- Vazquez, Eguíluz, San Miguel PRL 2008

Full reply text:

**Verdict: KNOWN.** This is the three-opinion coevolving (adaptive) voter model with rewire-to-same rewiring. The first pass called it a variant; I'm overturning that, because both observed behaviours are documented and nothing new shows up in the behaviour itself.

**Identification.** The coevolving voter model consists of two dynamical processes: copying and rewiring (Holme and Newman, 2006; Vazquez et al., 2008; Durrett et al., 2012). A node can adopt the opinion of one of its neighbors chosen randomly. Concurrently, the connections between voters can be rewired, depending on the agreement or disagreement of states between two connected voters. In Holme and Newman's original version, whenever there is an interaction between two vertices (agents) with different opinion, with probability α ∈ (0,1), one of the vertices breaks the link and connects to a different vertex of the same opinion. That is exactly your rule 2. The three-type case has its own literature: Böhme and Gross have studied the three-opinion version of the coevolving voter model with rewire-to-same dynamics.

**Why the observed outcome is the expected one.** The model has a known split between two final states. Below a critical rewiring rate everything converges to one opinion; above it, the network breaks apart. The frozen phase corresponds to a fragmented network where each component is in a consensus state. The fragmentation transition between active and fragmented phases is a peculiar characteristic of the coevolving voter model. Your settings put rewiring well below copying: rewiring happens with probability 0.1, while copying happens for every node at every step. So you'd expect the consensus side, and that's what you see. A consensus gain of about 0.635 fits: the starting majority of about 0.35–0.37 rises to full agreement (1.0). The edge turnover of about 0.97 is expected too, because discordant links keep getting rewired while the opinions are still coarsening.

**Differences that don't change the behaviour:**
- **Implementation details.** Updating is synchronous, copying is unconditional, and rewiring is a separate step. Durrett et al. note their own implementation differences and still find the same kind of transition: Despite the differences in implementation, this rewire-to-same model has a phase transition similar to that of Holme and Newman.
- **No synergy between the rules.** The sub-program results show the two rules don't combine into anything new. Copying alone gives a consensus gain of 0.625 against 0.635 for the full program; that's the ordinary voter model on a finite graph. Edge turnover needs rule 2 by construction. The full program differs from each sub-program, but only because each sub-program is missing one rule's direct effect.

**Trivial or artifact?** Neither, strictly. Nothing diverged numerically, there's no inert field, and the coloring isn't random. But both features are close to automatic. Edge turnover is close to a direct readout of rule 2, and consensus gain is plain voter-model absorption into a frozen consensus state. That frozen state is the documented absorbing state of the model, not a simulation bug. The "knock-out = 0" comparison is a weak test, because with every interaction removed nothing can change at all.

The search was a single query, but the identification is direct: the rule matches the Holme–Newman rewire-to-same model and the three-opinion case is studied by name.

---

## review:025

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_types

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_types (number of surviving types): median 1.0 with interaction, 3.0 with interaction knocked out (difference -2.0, robust spread 0.08); difference at 4x area -1.0, at 3x run length -2.0

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_types 2.0
  - rules 2 only: g_n_types 3.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Majority-rule / zero-temperature Potts coarsening to consensus (absorbing homogeneous state), with Schelling-type relocation acting as interface noise

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Only feature is g_n_types -> 1, i.e. fixation into a single-type absorbing state (all-empty or all-one-type, after which relocation is inert). Rule-1-only freezing at 2 types matches known majority/Potts freezing (q>=3); adding relocation noise removes the frozen interfaces and lets the system reach consensus, also documented. The 4x-area result of 2 types fits finite-size consensus timing. The knock-out comparison is trivially guaranteed. Frozen absorbing state, not an emergent pattern.

Citations:
- https://arxiv.org/pdf/1811.11888
- https://arxiv.org/html/2604.10082
- https://arxiv.org/pdf/2208.08770
- https://arxiv.org/html/physics/0609079
- https://arxiv.org/html/cond-mat/0407622

Full reply text:

**Verdict: KNOWN.** The only feature observed is that the number of surviving types drops to 1. That is the well-documented outcome of majority-rule dynamics: domains coarsen until one type takes over the whole lattice. It is a frozen absorbing state, so it is a trivial result rather than an emergent pattern.

**Why it's documented**
- Rule 1 is the standard majority-vote cellular automaton. It is known as the zero-temperature limit of the Ising model, a cell takes on the state held by the majority of its neighbors… analogous to the voter model… also… the zero-temperature limit of the Ising Model. With 3 types it is the zero-temperature 3-state Potts model.
- The end state is textbook. This dynamics leads to a coarsening mosaic of ↑ and ↓ domains that usually, but not always (see e.g., [54, 55]), ends in consensus. Also, in any majority model, there are at least two absorbing states, constituted by a homogeneous configuration of sites.
- The sub-program results fit known mechanisms:
  - **Rule 1 alone (2 types survive):** majority dynamics can get stuck before consensus. For the 3-state Potts model on a square lattice, freezing is observed for any value of q ≥ 3 for systems which are large enough.
  - **Adding rule 2 (1 type survives):** the relocation rule acts as noise that breaks up these frozen interfaces. This is also documented: straight interfaces, which are stable in the kinetic Ising model and prevent the ground state from being reached, are unstable in majority rule. Thus consensus is always achieved.
- Rule 2 is the Schelling relocation rule. It is well known to map onto Ising-type models: Schelling's model [1] of residential segregation is mathematically equivalent to the zero-temperature spin-exchange Kinetic Ising model with vacancies. Combining a majority-rule state change with a second rule that rewires who is next to whom has also been studied (the majority rule for the change in the automata state and the rule for the change in topological relations between automata).

**Why it's trivial**
- The surviving type is either type 0 (the whole lattice is empty) or a single occupied type (no empty cells are left, so rule 2 can no longer do anything). Either way the system is frozen and homogeneous.
- The 4x-area result (2 types survive, difference −1.0) fits ordinary finite-size consensus timing: bigger lattices take longer to reach consensus. It is not a new effect.
- The comparison against the knocked-out version is uninformative. Without interaction, nothing changes type, so all 3 types survive automatically.

**Overturning the first pass:** "variant" is too generous. The only measured feature, ending on one type, has no twist in its behaviour beyond the textbook result.

---

## review:026

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) neighbouring cells swap contents (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_reversals, g_activity_cv

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 5.893 with interaction, 7.3553 with interaction knocked out (difference -1.4623, robust spread 0.265); difference at 4x area -1.6757, at 3x run length -1.498
  - g_reversals (cells flip back and forth): median 0.9417 with interaction, 0.0 with interaction knocked out (difference 0.9417, robust spread 0.0188); difference at 4x area 0.9538, at 3x run length 0.9591
  - g_activity_cv (bursty activity): median 0.7518 with interaction, 0.0 with interaction knocked out (difference 0.7518, robust spread 0.015); difference at 4x area 0.5731, at 3x run length 0.7554

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 7.3582; g_reversals 0.5002; g_activity_cv 0.7561
  - rules 2 only: g_n_domains 6.4019; g_reversals 0.5357; g_activity_cv 0.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Schelling segregation model with vacancies and noise (random neighbour exchange ~ infinite-temperature Kawasaki noise); segregated/mixed competition

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Only g_n_domains reflects genuine (weak) Schelling segregation partly undone by swap noise; rule-2-only already gives 6.40 vs 5.89 full. g_activity_cv is entirely produced by the unconditional swap rule (rule-1-only 0.756 vs full 0.752), and g_reversals is mostly swap noise (rule-1-only 0.50) plus relocation/swap undoing. The knock-out scores of 0 for those two features are artifacts of the knock-out removing the swap rule. Not a new variant.

Citations:
- Gauvin, Vannimenus & Nadal, Phase diagram of a Schelling segregation model, EPJ B 70:293 (2009), arXiv:0903.4694
- Dall'Asta, Castellano & Marsili, Statistical physics of the Schelling model of segregation, arXiv:0707.1681
- Aging effects in Schelling segregation model, Sci. Rep. (2022), s41598-022-23224-7

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass. This isn't a meaningful variant. It's the Schelling segregation model with random noise added, and that model is well studied.

**What the system is**
- Rule 2 is the classic Schelling rule on a Moore (8-neighbour) lattice with vacancies, tolerance 0.5, and unhappy agents jumping to any empty cell. Gauvin, Vannimenus and Nadal studied exactly this kind of variant. In their version agents move one by one to *any* vacancy which has a satisfying neighborhood - this is equivalent to long-range diffusion in physical terms. They also let satisfied agents move, and note that satisfied agents can also move, not just discontent ones. That rule introduces some noise.
- Rule 1 swaps neighbours every step whatever their state. In physics terms that's spin-exchange (Kawasaki-type) dynamics at infinite temperature, i.e. plain mixing noise. The Kawasaki idea is that nearest-neighbor spins are exchanged at a rate $e^{-\beta\Delta E}$ … the order parameter remains conserved throughout the dynamics.
- The combination, one force that segregates plus one that mixes, is the standard "noisy Schelling" setup. Later work describes it as one where the movement of unsatisfied agents behaves as a driver for the system dynamics, while the motion of satisfied agents plays the role of noise. The model is known to have 3 different phases: frozen, segregated and mixed. Its segregation is understood as coarsening: static phase transitions typical of kinetic constrained models, nontrivial coarsening like in driven-particle systems and percolation-related phenomena.

**Going through each feature**
- **g_n_domains (fewer domains):** This is real but weak. Domains drop from 7.36 to 5.89, and rule 2 alone gives 6.40. That fits Schelling segregation being partly undone by the swap noise, which is the documented segregated-versus-mixed competition. It's genuine, but it's the textbook behaviour.
- **g_reversals (cells flip back and forth):** This is mostly an artifact. Rule 1 alone already gives 0.50. The knock-out scores 0.0 because removing every rule that reads other cells also removed the swap rule. The extra amount in the full program (0.94 vs about 0.5) is just swaps and relocations undoing each other, which is expected flickering in a noisy steady state.
- **g_activity_cv (bursty activity):** This is an artifact. Rule 1 alone gives 0.756, the same as the full program's 0.752. Rule 2 alone gives 0.0. So the "burstiness" comes entirely from the swap rule, not from any interaction. The knock-out comparison hides this for the same reason as above.

**Twist check:** I found nothing in how the system behaves that isn't covered by noisy Schelling or segregation-as-coarsening. The only new thing is the particular combination of rules, which doesn't count. I ran two searches and drew on four sources: the papers by Gauvin et al. and Dall'Asta et al., a 2022 study of the noisy Schelling model, and a review of Kawasaki dynamics.

---

## review:027

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) neighbouring cells swap contents (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_reversals, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 6.0856 with interaction, 7.3473 with interaction knocked out (difference -1.2617, robust spread 0.2687); difference at 4x area -1.3146, at 3x run length -1.272
  - g_reversals (cells flip back and forth): median 0.9531 with interaction, 0.0 with interaction knocked out (difference 0.9531, robust spread 0.0191); difference at 4x area 0.9569, at 3x run length 0.9706
  - g_moran (like-valued cells cluster spatially): median 0.3568 with interaction, 0.0006 with interaction knocked out (difference 0.3562, robust spread 0.0168); difference at 4x area 0.3767, at 3x run length 0.3786

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 7.3495; g_reversals 0.5946; g_moran -0.0079
  - rules 2 only: g_n_domains 6.4094; g_reversals 0.6; g_moran 0.1902
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Schelling segregation model (two agent types + vacancies, tolerance 0.5, random relocation) with Kawasaki-type swap noise / noisy Schelling model

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Type 0 = vacancy, so this is literally Schelling's two-group lattice model at T=0.5 plus random neighbour swaps (temperature-like Kawasaki noise). Clustering and coarsening (g_moran, g_n_domains) are canonical. Noise-assisted extra clustering compared with rule 2 alone matches the documented role of noise in escaping large-but-finite cluster states. g_reversals is a near-trivial side effect of swap noise and unhappy-agent churn (each rule alone already gives about 0.6). Overturn is mild: from 'variant' to 'known', because no new phenomenology appears. Caveat: Moran's I on nominal labels 0/1/2 is a questionable metric, but the effect size makes real segregation clear.

Citations:
- https://arxiv.org/abs/0903.4694 (Gauvin, Vannimenus, Nadal - Phase diagram of a Schelling segregation model)
- https://www.nature.com/articles/s41598-022-23224-7 (Aging effects in Schelling segregation model; noisy Schelling phases)
- https://arxiv.org/pdf/1504.03809 (Barmpalias et al.; swap dynamic = Kawasaki)
- Dall'Asta, Castellano, Marsili 2008, Statistical physics of the Schelling model of segregation

Full reply text:

**Verdict: KNOWN.** This is the standard Schelling segregation model with added noise. The first pass called it a variant, but nothing in the behaviour itself is new, so I'm overturning that, though only slightly.

**What the system actually is.** Type 0 means empty, so the "3 types" are really two kinds of agent plus vacancies. That is Schelling's own setup. In Schelling's model, individuals and vacancies are randomly arranged on a 2D lattice, and each individual has a tolerance threshold of 0.50: it is unhappy when more than half its neighbours are unlike it. Rule 2 is exactly this, using 8 neighbours and moves to random empty cells. Rule 1, the random neighbour swaps, is a known type of noise. In the statistical-physics treatments, Young's pairwise-swapping dynamic is known to statistical physicists as the *Kawasaki dynamic*. Also, including the motion of satisfied agents leads to a noisy effect playing the role of temperature in a statistical physics approach.

**Checking each observed feature against the literature:**
- **g_moran (clustering, 0.36 vs 0.0006 with interaction removed):** This is Schelling segregation itself. At tolerance 0.5 with few vacancies, the system evolves from a random configuration – where the vacancies and the two types of agents are intimately mixed – to a completely segregated configuration. One detail looks like it might be new but isn't. Rules 1+2 cluster more than rule 2 alone (0.36 vs 0.19), and that is also documented. Gauvin et al. describe noise as useful to avoid a quirk of the original model, namely that the system may end up in states where the clusters are large but finite, so that strictly speaking no large-scale segregation occurs.
- **g_n_domains (somewhat fewer domains):** This is ordinary coarsening. Clusters become more and more compact and well separated spatially, their surface gets less corrugated, in a process strongly reminiscent of coarsening effects in alloys. The effect here is small (6.1 vs 7.3), which fits a 64×64 grid in a noisy stationary state, where small clusters are continuously created and annihilated.
- **g_reversals (cells flip back and forth):** This is mostly a side effect of the noise and the churn of unhappy agents, not a separate pattern. Each rule alone already gives about 0.6, because swaps and relocations constantly change what occupies a cell. Having both rules just pushes it to about 0.95. The knock-out score of 0 only reflects that every rule there reads neighbours, so removing them leaves nothing that moves.

**Genuine pattern or artifact?** The clustering is a real emergent pattern: there is no frozen state, divergence or random colouring. The reversals feature is close to trivial. One caveat: if g_moran treats the type labels 0/1/2 as numbers, it is measuring categories as if they were values. The large effect would almost certainly survive a proper same-type measure, though.

**Searches I ran:** "Schelling segregation model noise random swaps multiple groups statistical physics" and "Schelling model three groups vacancies tolerance 0.5 cluster formation noise".

---

## review:028

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell; (2) each cell copies the type of a random neighbour (every step).

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_reversals, g_spec_peak

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 5.8607 with interaction, 7.3366 with interaction knocked out (difference -1.4759, robust spread 0.2639); difference at 4x area -1.6347, at 3x run length -1.5138
  - g_reversals (cells flip back and forth): median 0.8921 with interaction, 0.0 with interaction knocked out (difference 0.8921, robust spread 0.0178); difference at 4x area 0.9122, at 3x run length 0.9071
  - g_spec_peak (a dominant temporal frequency): median 0.3103 with interaction, 0.0 with interaction knocked out (difference 0.3103, robust spread 0.0088); difference at 4x area 0.3828, at 3x run length 0.3663

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 6.4465; g_reversals 0.4365; g_spec_peak 0.02
  - rules 2 only: g_n_domains 6.4481; g_reversals 0.6464; g_spec_peak 0.2728
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Schelling-Voter model (voter-model interfacial-noise coarsening with Schelling relocation)

Trivial or artifact: False; confidence: medium

Twist (if any): 

Notes: The knock-out is a frozen grid (both rules read neighbours), so the large differences for reversals and spectral peak are inflated. Copying alone accounts for most of reversals (0.65 of 0.89) and the spectral peak (0.27 of 0.31). The spectral peak is probably a synchronous-update flicker artifact rather than a real oscillation. The domain count reflects unfinished, logarithmically slow voter coarsening on a 64x64 grid. Differences from the published Schelling-Voter model are in the rules only (unconditional synchronous copying, empty cells copied), not in the behaviour.

Citations:
- https://www.sciencedirect.com/science/article/abs/pii/S0378437117306118
- https://www.sciencedirect.com/science/article/abs/pii/S0960077913001720
- https://arxiv.org/html/1506.05321
- https://arxiv.org/pdf/1506.05321

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass (which said VARIANT), though only slightly. What you're seeing is ordinary voter-model coarsening plus the noisy flipping at domain borders that the voter model is known for, with Schelling-style moving added. A 3-state voter model is one where each cell copies a random neighbour's type, and here "empty" counts as one of the three types. The emergent feature is real, but the main headline numbers are inflated by a knock-out that does nothing.

**The system is already in the literature.** It's the Schelling–Voter hybrid. That model is a mix of the classical Schelling segregation model with the Voter model, in which an unhappy agent can change her place or her state to be happy. It has been studied on 2D lattices, where the outcomes include the extinction of one of the states, the coexistence of both, segregated patterns with conglomerated clusters of agents of the same state, and a diluted region. The earlier paper on it also used a square lattice of side L. Your version differs only in the rules: every cell copies a neighbour every step, all at once rather than only when unhappy, and empty cells get copied too. That's a different rule combination, not a different behaviour.

**What each feature is:**
- **g_n_domains (several domains):** This is just unfinished coarsening on a 64×64 grid. In the 2D voter model, there is a coarsening process which brings the system towards the single-domain state, but it is very slow, since the density of active interfaces vanishes only as 1/ln(t/τ). The domains stay very jagged and preserve their fractal geometry even at the late stages of evolution. Having fewer domains than the knock-out simply means some ordering happened.
- **g_reversals (cells flipping back and forth):** This is the voter model's border noise. There, the dynamics are driven by interfacial noise. Because every cell copies a neighbour every step, border cells flicker constantly. Copying alone already gives 0.65 of the full 0.89.
- **g_spec_peak (a dominant frequency):** This is most likely an artifact of every cell updating at the same moment, which tends to produce period-2 flicker or a spectrum shaped by the update scheme. It's probably not a real oscillation: copying alone already gives 0.27 of the 0.31, and moving alone gives nearly nothing. I haven't checked the spectrum directly, so this is an inference.

**Why the comparisons are weak:**
- **The knock-out is empty.** Both rules read neighbours, so removing them leaves a grid that never changes. That explains why the knock-out values for reversals and spectral peak are exactly 0. Beating 3 spreads against a frozen grid tells you almost nothing. The fair baselines are the single-rule runs.
- **Against those, the extra effect of combining the rules is modest.** Moving cells around mixes the grid and adds border activity, which is the documented interplay in the Schelling–Voter papers. The effect is real but small.

**Genuine or trivial?** There is real but well-known emergent structure: coarsening domains with noisy borders. The frequency peak looks like an update-scheme artifact, and the size of the knock-out differences is misleading. No part of the behaviour appears to be undocumented.

I ran two searches: one for Schelling combined with the voter model on a lattice, and one for synchronous voter-model coarsening and border noise.

---

## review:034

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
C kC=4 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0  (compact grammar notation)

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_corrlen, g_activity_trend, g_moran

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 5.4183 with interaction, 7.6785 with interaction knocked out (difference -2.2602, robust spread 0.2619); difference at 4x area -2.3227, at 3x run length -2.3303
  - g_corrlen (long spatial correlation length): median 3.0 with interaction, 1.0 with interaction knocked out (difference 2.0, robust spread 0.08); difference at 4x area 2.0, at 3x run length 1.0
  - g_activity_trend (activity rising): median -0.7415 with interaction, 0.0 with interaction knocked out (difference -0.7415, robust spread 0.0155); difference at 4x area -0.9015, at 3x run length -0.8283
  - g_moran (like-valued cells cluster spatially): median 0.5784 with interaction, -0.0017 with interaction knocked out (difference 0.58, robust spread 0.016); difference at 4x area 0.6014, at 3x run length 0.5275

### Stage 3 reply: verdict KNOWN, standard name: Multi-state (q=4) voter model / neutral copy-a-random-neighbour dynamics: voter-class coarsening toward consensus

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: The rule is the voter model with 4 states. The NEXT_GE(th=1) gate only removes copies that would change nothing, so it at most rescales time. Observed features: fewer domains, longer correlation length, falling activity (the feature label says rising but the value is negative) and positive Moran's I. These are the textbook signs of voter coarsening in 2D (interface density falls as about 1/ln t, length scale grows with z=2). The pattern is real (Moran ~0 and corrlen 1 with interaction knocked out) but transient, heading for an absorbing consensus state. First pass called it a VARIANT; nothing about the behaviour is new, so the verdict is upgraded to KNOWN.

Citations:
- Starnini, Baronchelli, Pastor-Satorras, Ordering dynamics of the multi-state voter model, J. Stat. Mech. 2012 (arXiv:1207.5810)
- Dornic, Chate, Chave, Hinrichsen, Critical coarsening without surface tension: the universality class of the voter model, PRL 87, 045701 (2001)
- Tartaglia, Cugliandolo, Picco, Percolation and coarsening in the bidimensional voter model, arXiv:1506.05321
- Dall'Asta & Galla, Algebraic coarsening in voter models with intermediate states, arXiv:0806.0817
- Ben-Naim, Frachebourg, Krapivsky, Coarsening and persistence in the voter model, PRE 1996

Full reply text:

**Verdict: KNOWN.** This is the standard multi-state (q = 4) voter model on a lattice, and what you observed is its well-documented coarsening. I'm overturning the first pass's "variant" label: nothing about the behaviour itself counts as a twist.

**Identifying the system.** Each cell has one of 4 states. If at least one neighbour qualifies under the threshold, the cell copies the state of a random neighbour, with probability 1. Copying a neighbour that already matches changes nothing. So whether the threshold means "at least one neighbour" or "at least one *different* neighbour", the dynamics are the same as the voter model, just possibly running at a different speed. The voter model is defined as at each step a randomly chosen spin copies the state of one of its nearest neighbor, chosen at random. The version with more than two states is studied directly in "Ordering dynamics of the multi-state voter model" (Starnini et al., arXiv:1207.5810). That paper perform[s] a numerical study of the model in low dimensional lattices, comparing the case of multiple states with the usual binary voter model. It also notes the model has been treated in terms of mappings of coarsening of the Potts model.

**How each observed feature matches the known behaviour:**
- **Fewer domains (g_n_domains lower than the knock-out) and longer correlation length (g_corrlen 3 vs 1):** domains merge over time. The voter model's evolution shows spatial clustering of similar opinions. It approaches one of the two absorbing states with complete consensus via a coarsening process in d≤2. The typical domain size grows: the dynamical exponent, relating to the correlation length and properties of the dynamical structure factor, is still found to be z=2.
- **Falling activity (g_activity_trend −0.74):** the feature is labelled "activity rising", but the value is negative, so activity is actually falling. That fits the voter model's slow loss of active boundaries between domains. In 2D there is a coarsening process which brings the system towards the single-domain state, but it is very slow, since the density of active interfaces vanishes only as 1/ln(t/τ).
- **Like-valued cells clustering (g_moran ≈ 0.58):** this follows directly from the domain growth above.
- **Having 4 states instead of 2:** this doesn't change the kind of behaviour. Starnini et al. report how many states survive in 2D: s(t) ∼ t −1 log t for d = 2. Algebraic (faster) coarsening has been reported for voter models that add intermediate states (the logarithmic coarsening of the conventional voter model in two dimensions is eliminated in favour of an algebraic decay). Your rule has no intermediate states, so that doesn't apply here.

**Is it genuine or an artifact?** It's a genuine interaction-driven pattern: Moran's I is about 0 and the correlation length is 1 when the interaction is knocked out. It is not a random colouring, a diverged field or a frozen state. But it is a slow transient on the way to an absorbing state: given enough time a finite lattice ends up all one state. The measured correlation length is short (3), and the gap shrinks at 3x run length (2.0 → 1.0). That probably reflects the domain count and the late-time averaging rather than any non-standard physics, but it is at least consistent with ongoing coarsening. None of it is new.

**Searches run:** voter-model coarsening and logarithmic interface decay in 2D; multi-state / q-state voter model coarsening on lattices.

---

## review:036

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
N kN=4 | N[any] ALWAYS() -> COPY_RAND() p=0.2  (compact grammar notation)

OBSERVED (features that passed the simulation checks):
  n_oscillation

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - n_oscillation (node states oscillate over time): median 0.7199 with interaction, 0.0 with interaction knocked out (difference 0.7199, robust spread 0.0144); difference at 4x area 0.6997, at 3x run length 0.3663

### Stage 3 reply: verdict KNOWN, standard name: 2D voter model (nearest-neighbour, asynchronous/random-sequential); interfacial-noise-driven coarsening with 1/ln t decay of active interfaces

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Rule is exactly the voter model (copy random von Neumann neighbour with prob 0.2 per step = rescaled time), not a variant. 'n_oscillation' is aperiodic stochastic flipping at domain interfaces; its decline at 3x run length (0.72->0.37) matches logarithmic coarsening toward an absorbing consensus state, and insensitivity to 4x area matches domains smaller than system size. The knock-out removes the only dynamics, so the contrast is trivially guaranteed. Overturns first pass only in calling it the canonical model rather than a variant; no novel phenomenology. The 'oscillation' label is an artefact of the metric (no periodicity).

Citations:
- Dornic, Chaté, Chave, Hinrichsen, 'Critical Coarsening without Surface Tension: the Voter Universality Class', PRL 2001, arXiv:cond-mat/0101202
- Tartaglia, Cugliandolo, Picco, 'Percolation and coarsening in the bidimensional voter model', arXiv:1506.05321
- Ben-Naim, Frachebourg, Krapivsky, 'Coarsening and persistence in the voter model' (opinion-change statistics), PRE 1996
- Interplay of interfacial noise and curvature driven dynamics in two dimensions, arXiv:1701.09117

Full reply text:

**Verdict: KNOWN.** I agree that this is a known system, but I'd go further than the first pass. It isn't a variant. It is the ordinary voter model, and the "oscillation" it reports is just individual sites changing state while domains slowly merge.

**What the rule is.** Each node on a square lattice with 4 neighbours, with probability 0.2 per step, copies the state of one randomly chosen neighbour. That is the definition of the 2D voter model: a caricatural process in which sites on a square lattice adopt the opinion of a randomly-chosen neighbor. The p=0.2 update probability only rescales time (asynchronous updating). It doesn't change the model's class.

**Why the observed behaviour is documented and expected.**
- **The behaviour is known.** Its coarsening process, which gives rise to patterns with clusters of all sizes between 1 and √t, is characterized by a slow, logarithmic decay of the density of interfaces ρ∼1/ln t. Sites can only flip at domain boundaries, so the flip activity follows the interface density. Changes of opinion per site have been studied exactly: exact results are obtained for the average number of opinion changes (from a preprint listed on the ResearchGate page for Dornic et al.).
- **The decay fits.** In 2D, there is a coarsening process which brings the system towards the single-domain state, but it is very slow, since the density of active interfaces vanishes only as 1/ln(t/τ). Your measurements show this. The signal is about the same at 4x area (0.70 vs 0.72), as expected while domains are still much smaller than the lattice. It roughly halves at 3x run length (0.37). That is the system moving towards consensus, which is the all up and all down absorbing states, reached in a time that typically behaves as L^2 log L.
- **The knock-out comparison is trivially true.** Removing the copy rule leaves no dynamics at all, so 0 activity is guaranteed. The test only shows that the one rule does something.

**Genuine pattern or artefact?** The coarsening is real, well-understood voter-model behaviour. But labelling it "oscillation" is an artefact. There is no periodicity, no limit cycle and no sustained oscillation. It is random flipping at domain walls that dies away as the system heads to a frozen consensus state. As an "oscillation" feature it should count as trivial. Running longer, or measuring a power spectrum or autocorrelation, should show no characteristic frequency and activity falling roughly as 1/ln t.

---

## review:052

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a random network of 200 nodes (average 4 links each). Nodes have 4 types (random start). Rules, applied in order each step: (1) each node adopts its neighbours' majority type (every step); (2) a node linked to a different-type neighbour cuts that link and links to a random same-type node (with probability 0.03 per step).

OBSERVED (features that passed the simulation checks):
  n_consensus_gain, n_type_entropy, n_edge_turnover

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - n_consensus_gain (nodes converge to agreement): median 0.6725 with interaction, 0.0 with interaction knocked out (difference 0.6725, robust spread 0.0134); difference at 4x area 0.6931, at 3x run length 0.6725
  - n_type_entropy (node types evenly mixed): median 0.1175 with interaction, 1.3777 with interaction knocked out (difference -1.2602, robust spread 0.0299); difference at 4x area -1.2363, at 3x run length -1.2602
  - n_edge_turnover (edges rewire): median 0.1948 with interaction, 0.0 with interaction knocked out (difference 0.1948, robust spread 0.0285); difference at 4x area 0.2416, at 3x run length 0.1948

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: n_consensus_gain 0.685; n_type_entropy 0.1311; n_edge_turnover 0.0
  - rules 2 only: n_consensus_gain 0.0; n_type_entropy 1.3777; n_edge_turnover 0.7629
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Coevolving/adaptive majority-rule (adaptive voter-type) model with homophilic rewiring; Holme-Newman coevolution of networks and opinions; majority rule + rewiring on coevolving networks (Fernandez-Gracia et al.)

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Consensus gain and entropy drop are essentially reproduced by majority rule alone (0.685/0.131 vs 0.6725/0.1175), so the coupling adds little to these measures. Edge turnover follows directly from the rewiring rule and drops in the full program because majority rule removes discordant links. Not a numerical artifact, but the features are standard majority-rule ordering plus rewiring that happens by construction. Incomplete consensus (~0.67) fits documented majority-rule traps on sparse networks. The measured features don't probe network fragmentation, which is the distinctive behaviour of the coupled system and is itself documented. The knock-out baseline (a frozen random colouring) is a weak control.

Citations:
- https://arxiv.org/pdf/1403.6106
- https://arxiv.org/abs/physics/0603023
- https://arxiv.org/html/2507.11521v2

Full reply text:

**Verdict: KNOWN.** This overturns the first pass's VARIANT call. I can't find any twist in the observed behaviour that the literature doesn't already describe.

**What the system is.** Majority-rule opinion updating combined with homophilic link rewiring is a standard coevolving (adaptive) network model. The closest documented match is the link-state study by Fernández-Gracia and co-authors (arXiv:1403.6106), which has the same two ingredients: majority rule plus rewiring. They report that while large systems evolving under the majority rule alone always fall into disordered topological traps composed by frustrated links, any amount of rewiring is able to drive the network to complete order. Depending on relative rates, it reaches either a one-component network with all links in the same state or a network fragmented in two components with opposite states. The opinion-update rule differs: here it's majority rule rather than copying one neighbour. Otherwise the multi-opinion version is the Holme–Newman model, which undergoes a continuous phase transition as this parameter is varied, from a regime in which opinions are arbitrarily diverse to one in which most individuals hold the same opinion. In the general adaptive voter model, an agent either adopts the opinion of a neighbor or drops this connection to create a new one with an individual having the same opinion. Rewiring until consensus is reached is also documented.

**Your sub-program data says the coupling adds almost nothing new:**
- **Consensus and entropy come from rule 1 alone.** Majority rule by itself gives consensus 0.685 and entropy 0.131, against 0.6725 and 0.1175 for the full program. That is majority-rule ordering on a random network, which is textbook. Consensus stops at about 0.67 rather than 1, which fits the documented "topological traps" and the stuck or flip-flopping states that simultaneous majority updates produce on sparse networks.
- **Edge turnover is built in.** Rule 2 rewires disagreeing links by definition. In the full program, turnover falls from 0.76 to 0.19 simply because majority rule removes most disagreeing links. That outcome is expected, not emergent.
- **The control is weak.** The knock-out removes both rules, which leaves a frozen random colouring. Beating that baseline only shows that the rules do something.

**Is it real or an artifact?** Not a numerical artifact, but not a new emergent pattern either. The ordering is genuine majority-rule coarsening, and the turnover is a direct readout of the rewiring rule. The features measured don't capture the one thing that would distinguish the coupled system: whether the network splits into same-type components (the fragmentation transition). Even that is already documented.

**Searches run:** "majority rule opinion dynamics adaptive network rewiring fragmentation consensus" and "Holme Newman nonequilibrium phase transition coevolution networks opinions".

---

## review:057

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a type-1 cell spontaneously becomes type 0 (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.5 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 5.4681 with interaction, 6.3612 with interaction knocked out (difference -0.8932, robust spread 0.2366); difference at 4x area -0.9016, at 3x run length -0.8932
  - g_reversals (cells flip back and forth): median 0.7064 with interaction, 0.0 with interaction knocked out (difference 0.7064, robust spread 0.0141); difference at 4x area 0.7069, at 3x run length 0.6711

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 6.3612; g_reversals 0.0
  - rules 2 only: g_n_domains 6.3986; g_reversals 0.5
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Schelling segregation model (two-type, high vacancy density, random-destination relocation): diluted segregated phase with persistent churn

Trivial or artifact: False; confidence: high

Twist (if any): 

Notes: Rule 1 (type-1 cells become empty with p=0.1) is a one-way transient: type 1 is effectively gone after about 50 steps, so late-time behaviour is two-type Schelling with threshold 0.5 at about 67% vacancy, which is the documented diluted-segregation regime. Fewer domains than the knock-out = small like-type clusters; reversals = ongoing relocation churn, which is expected because moves go to random empty cells rather than satisfying ones. The knock-out removes all movement, so reversals=0 there by construction. 'Overturns first pass' only in the sense of upgrading VARIANT to KNOWN; no phenomenological twist found. Caveat: at high vacancy many agents have no occupied neighbours (0/0 share); how the code handles that probably sets much of the reversal rate, a partial implementation artifact.

Citations:
- Gauvin, Vannimenus & Nadal, Phase diagram of a Schelling segregation model, Eur. Phys. J. B 70, 293-304 (2009), arXiv:0903.4694
- A Schelling model with switching agents: decreasing segregation via random allocation and social mobility, arXiv:1212.4929
- Aging effects in Schelling segregation model, Sci. Rep. (2022), s41598-022-23224-7

Full reply text:

**Verdict: KNOWN.** The first pass called this a variant of a known model. I'm going a step further: there's no new behaviour here. What's left once the decay rule has played out is a standard two-type Schelling model with a high vacancy rate, run in the diluted, constantly moving regime. That regime is documented.

**What the system turns into.** Rule 1 is one-way and has an end. Each type-1 cell survives a step with probability 0.9, so after about 50 steps under 1% are left. The "late part of the run" is therefore a 2-type Schelling model: threshold 0.5, 8 neighbours, unhappy agents jump to a random empty cell. If the three types start at roughly a third each, about two thirds of the grid ends up empty. The decay rule is a startup transient that removes one population and raises the vacancy rate. It doesn't add any ongoing interaction.

**How the two features match the literature.**
- **Many separate domains (g_n_domains).** Gauvin, Vannimenus and Nadal (2009) mapped the Schelling model against tolerance and vacancy density. At high vacancy they describe a diluted segregated phase where the two types of agents are not mixed but there are domains with many small clusters of a same color in a sea of vacancies. Later summaries place it for high values of vacancies (V > 0.65) and low values of tolerance (S > 0.5). At about 67% empty with threshold 0.5, this system sits right at that boundary. Having fewer domains than the knock-out is what you'd expect once like-type agents cluster.
- **Cells flipping back and forth (g_reversals).** In this rule set, unhappy agents move to any random empty cell, not only cells where they'd be happy. So the dynamics don't freeze. A study of a Schelling model with agents who switch type says the same: unsatisfied agents move because moving does not depend on their finding a more "welcoming" site. In related versions the dynamics settle into a state where small clusters are continuously created and annihilated. Sites switching between empty and occupied is simply this ongoing churn.

**Why the knock-out comparison overstates the result.**
- Rule 2 is the only rule that reads neighbours. Removing it removes all movement, so g_reversals = 0 in the knock-out is guaranteed by construction.
- The "rule 2 only" run still shows churn (0.5). The full program's 0.71 most likely comes from the higher vacancy rate after type-1 cells decay: more isolated or poorly surrounded agents keep moving. Changing the vacancy density is a parameter shift within the known phase diagram, not new behaviour.

**Possible artifact.** With about 67% vacancies, many agents have no occupied neighbours, so "same-type share among occupied neighbours" is 0/0. The code has to choose whether those agents count as happy or unhappy. If they count as unhappy, much of g_reversals may be isolated agents wandering endlessly rather than anything collective. The domain structure is a genuine (if weak) pattern. The reversal rate is real ongoing movement, but its size depends on that implementation choice.

---

## review:058

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a type-1 cell spontaneously becomes type 0 (with probability 0.1 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_activity_trend, g_reversals

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_activity_trend (activity rising): median -4.0035 with interaction, -3.3698 with interaction knocked out (difference -0.6337, robust spread 0.1475); difference at 4x area -0.5777, at 3x run length -0.5706
  - g_reversals (cells flip back and forth): median 0.8364 with interaction, 0.0 with interaction knocked out (difference 0.8364, robust spread 0.0167); difference at 4x area 0.8356, at 3x run length 0.6714

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_activity_trend -3.3698; g_reversals 0.0
  - rules 2 only: g_activity_trend -0.8338; g_reversals 0.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Schelling segregation model (two agent types + vacancies, random-vacancy relocation) with one-way removal of one population

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Type 0 = vacancy, so rule 2 is standard Schelling relocation at threshold 0.3 (Moore neighbourhood, random vacancy). Rule 1 removes type-1 agents at 10% per step, so they are essentially extinct within about 50 steps. After that every type-2 agent is fully satisfied and the lattice freezes into an absorbing state. g_reversals is just Schelling relocation seen at the cell level (empty/occupied toggling). It is nonzero only while decay keeps creating dissatisfied minorities. g_activity_trend is negative (activity falls, not rises), consistent with relaxation into the frozen single-type state. Passing the knock-out and the size/length checks reflects a size-independent per-cell mechanism, not genuine emergent structure. The first pass's 'variant' label is too generous: no new phenomenology was found.

Citations:
- https://www.nature.com/articles/s41598-022-23224-7
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9653388
- https://link.gale.com/apps/doc/A348647898/AONE
- https://arxiv.org/pdf/0711.2212
- https://arxiv.org/pdf/1002.3758
- https://en.wikipedia.org/wiki/Schelling%27s_model_of_segregation

Full reply text:

**Verdict: KNOWN.** I'm overturning the first pass's "variant" label in the direction of "known". Both observed features also look like a trivial relaxation transient, not a genuine emergent pattern.

**What the system is.** Type 0 is the empty cell, so there are two agent types (1 and 2) plus vacancies. Rule 2 is Schelling's segregation model in its "blind move" form. Each site has three states, σi=±1 for the two kinds of agents and σi=0 for vacancies. A non-zero vacancy density is needed for agents to change their location. In the blind version, an agent who moves goes to a site picked uniformly at random among vacant sites. The neighbourhood is the standard one, where agents check whether the fraction of neighbours matching their group—ignoring empty spaces—meets the threshold. Rule 1 just deletes type-1 agents at random, at a rate of 10% per step. That is a one-way removal of one population, not a new interaction. Schelling's own setup also removes agents, which makes some of the remaining agents unhappy and drives the evolution. Variants with a variable number of vacancies are already studied, for example the generalization to an open domain, with a parameter acting as a chemical potential for the number of vacancies.

**Why the features are trivial:**
- **g_reversals:** A cell "flips back and forth" because a dissatisfied agent leaves it and another agent later lands in it. That is just Schelling relocation seen at the cell level. It is zero with the rules that read neighbours removed, because nothing moves then. It is also zero with rule 2 alone, presumably because a two-type Schelling system at a 0.3 threshold quickly reaches a state where everyone is satisfied. The decay keeps creating type-1 minorities that become dissatisfied, so moves continue for a while.
- **g_activity_trend:** This is negative, so activity is falling, not "rising" as the label says. With a 10% per-step loss, about 0.9^50 ≈ 0.5% of type-1 agents remain after 50 steps. Once they are gone, every remaining type-2 agent has a same-type share of 1, so movement stops completely. Activity therefore falls towards a frozen state with only one type left. The interaction only changes how fast it falls, by about −0.6.
- **Robustness checks:** The size and run-length checks pass because the mechanism is per-cell and independent of system size. Passing them doesn't make the behaviour non-trivial.

I found no distinct phenomenology to call a twist. There are no oscillations, travelling structures or steady state that is sustained without being driven. The interaction knock-out passes only because Schelling relocation stops when rule 2's neighbour check is removed. I ran one search, and textbook Schelling dynamics are enough to explain everything.

---

## review:063

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 4 types (random start). Rules, applied in order each step: (1) each cell adopts its neighbourhood's majority type (with probability 0.005 per step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.7 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_domains, g_type_entropy, g_spec_peak

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_domains (many separate domains): median 6.3526 with interaction, 7.6662 with interaction knocked out (difference -1.3136, robust spread 0.2804); difference at 4x area -1.107, at 3x run length -5.1813
  - g_type_entropy (cell types evenly mixed): median 0.8484 with interaction, 0.9998 with interaction knocked out (difference -0.1513, robust spread 0.037); difference at 4x area -0.0853, at 3x run length -0.5177
  - g_spec_peak (a dominant temporal frequency): median 0.6579 with interaction, 0.0 with interaction knocked out (difference 0.6579, robust spread 0.0132); difference at 4x area 0.6041, at 3x run length 0.6592

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_domains 4.8559; g_type_entropy 0.9922; g_spec_peak 0.3811
  - rules 2 only: g_n_domains 7.6041; g_type_entropy 0.9998; g_spec_peak 0.0935
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Schelling-Voter model (Schelling segregation + majority/voter type-adoption), multi-type; coarsening toward type extinction/consensus

Trivial or artifact: True; confidence: medium

Twist (if any): 

Notes: Downgraded from KNOWN-VARIANT: no phenomenological twist found. Rule-2-only (threshold 0.7, 3 occupied types) stays in the documented disordered/no-coarsening Schelling regime; the majority rule adds coarsening toward one type taking over (documented extinction states). Entropy and domain-count effects grow strongly with run length, so this is a transient approach to an absorbing state, not a stationary pattern. g_spec_peak is very likely low-frequency trend power (knock-out = 0 because it is frozen), not a real oscillation. Domain/entropy changes are genuine interaction effects; the spectral feature is likely an artifact.

Citations:
- https://www.sciencedirect.com/science/article/abs/pii/S0378437117306118
- https://www.sciencedirect.com/science/article/abs/pii/S0960077913001720
- https://arxiv.org/pdf/1905.01953

Full reply text:

**Verdict: KNOWN.** I'm not keeping the first pass's "variant" label. As far as I can tell, this is a multi-type version of the documented Schelling–voter model. The features look more like a slow drift toward a state where one type takes over and nothing changes than like a stable new pattern.

**Why it's the known model**
- This is a Schelling relocation rule combined with a rule that copies the neighbourhood majority. That combination has been studied: a mix of the classical Schelling segregation model with the Voter model, in which an unhappy agent can change her place or her state to be happy. The stationary states reached entail the extinction of one of the states, the coexistence of both, segregated patterns with conglomerated clusters of agents of the same state, and a diluted region. The language-competition version is another source for the same pairing (a mix between Schelling's segregation model and the Voter model, where a probability parameter p governs which dynamics is used).
- Rule 2 alone does nothing measurable: entropy is 0.9998 and the domain count matches the knock-out. That fits known results. In the heavily studied two-type case, demanding thresholds stop clustering: no coarsening takes place for extremely intolerant individuals… the system remains trapped in a disordered state. Here a 0.7 own-type requirement with 3 occupied types and about 25% empty cells is very demanding, so cells keep moving without ever forming domains.
- The full program's changes match what a voter-type rule is known to do. Fewer domains and lower type entropy are signs of clusters merging and growing toward one type. These changes grow sharply with run length (entropy difference −0.15 at the original length, −0.52 at 3× length). That points to an unfinished slide toward one type dying out or taking over, which is the extinction outcome in the citation above, not a steady state.

**Trivial or artifact?** Partly.
- *Domain count and entropy:* real coupled behaviour, since neither sub-program reproduces them. But they are transient: what's being measured is a snapshot of an ongoing drift toward a frozen end state, not a stationary pattern.
- *Dominant frequency (g_spec_peak):* probably an artifact. The knock-out scores exactly 0.0, and the full program's high value barely changes with run length. Both fit a spectrum dominated by the lowest frequencies, meaning the slow downward trend, not a real oscillation. Unless the spectrum was detrended and shows a peak at a non-zero frequency, I wouldn't count this as oscillatory behaviour.

**Searches I ran:** "Schelling segregation model combined with voter model opinion dynamics lattice" and "Schelling model high tolerance threshold no segregation unstable dynamic state multiple groups". The results I cite only cover two types, so the extension to 3 types plus empty cells, with the majority rule applied to all cells, is my inference.

**Caveat:** I didn't find a paper with exactly 3 occupied types, a 0.7 threshold and this update order. That is a difference in the rules, not in what the system does.

---

## review:068

### Card given to the Stage 3 model (rule + simulation evidence)

RULE (plain words):
World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a type-0 cell becomes type 1 if at least 3 neighbours are type 1 (every step); (2) an occupied cell (type 0 = empty) whose same-type share among occupied neighbours is below 0.3 moves to a random empty cell.

OBSERVED (features that passed the simulation checks):
  g_n_types

SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, averaged over the late part of the run. The interaction knock-out is the same program with every rule that reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the original size, and keeps its sign and at least half its size at 4x area and at 3x run length.
  - g_n_types (number of surviving types): median 2.0 with interaction, 3.0 with interaction knocked out (difference -1.0, robust spread 0.1); difference at 4x area -1.0, at 3x run length -1.0

SUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature above; the full program has 2 rules):
  - rules 1 only: g_n_types 3.0
  - rules 2 only: g_n_types 3.0
No sub-program came within 3 robust spreads of the full program on every feature above.

### Stage 3 reply: verdict KNOWN, standard name: Vacancy exhaustion: polluted/threshold bootstrap percolation (Moore, r=3) filling the vacancies of a Schelling model, ending in a frozen jammed state

Trivial or artifact: True; confidence: high

Twist (if any): 

Notes: Type-2 count is conserved (moves only) and type-1 count only grows, so the lost type must be type 0 (empty cells). Rule 1 alone is polluted bootstrap percolation: type 2 acts as blockers that shield some empty cells forever, so 3 types survive. Schelling moves relocate the blockers and expose those empty cells to threshold-3 growth until none remain. With no empty cells, Schelling dynamics cannot move, giving a frozen absorbing state. The interaction effect is real but trivial; the outcome follows directly from the rules. This overturns the first pass only in the sense of KNOWN instead of a novel variant, and adds the trivial/absorbing flag.

Citations:
- https://arxiv.org/pdf/1705.01652
- https://www.math.ucdavis.edu/~gravner/papers/repl.pdf
- https://www.nature.com/articles/s41598-022-23224-7
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9653388
- https://arxiv.org/pdf/0707.1681

Full reply text:

**Verdict: KNOWN.** The result is a trivial absorbing state, not an emergent pattern. I agree with the first pass that this is a known system, but not with calling it a variant: the behaviour is just two documented models running together, with no new twist.

**Which type disappears.** The rules alone tell us. Rule 2 only moves type-2 cells around, so the number of type-2 cells never changes. Rule 1 only turns empty cells into type 1, so the number of type-1 cells can only go up. So the only type that can drop to zero is type 0, the empty cells. "2 surviving types" means every empty cell has been filled.

**What the two rules are on their own:**
- **Rule 1 is polluted bootstrap percolation.** In that model, sites start occupied, empty, or blocked. At subsequent steps, a vertex becomes occupied if it is not closed and it has at least r occupied neighbors. Here type 1 is "occupied", empty cells are "open", and type 2 acts as the blocker. The model's well-known result is that some open sites stay empty forever because blockers shield them. That's why rule 1 alone leaves all 3 types. The 3-state lattice version of this process is still being studied (a 3-state version is the subject of the present paper. We call our rule *two-stage bootstrap percolation*).
- **Rule 2 is the standard Schelling model with empty cells.** It only works if some cells are empty: We also need a non-zero vacancy density, ρv>0, for agents to change their location.

**Why combining them kills the empty type.** Unhappy agents moving at random keep relocating cells, so any empty cell left next to type-1 clusters can be filled. One-way growth then uses up every empty cell. With none left, no agent can move and the system stops. The Schelling literature already documents frozen, blocked states where unhappy agents remain but there are no vacancies which could accommodate them.

**Trivial or genuine?** Trivial. It's a frozen absorbing state, and the end result follows directly from the rules. The only measured feature is the type count, and it moves from 3 to 2 because the empty type is used up. It isn't evidence of a new spatial pattern. The fact that neither rule alone produces the effect is real, but it's simple to explain: movement keeps exposing empty cells to growth.

**Searches run:** the Schelling model with zero empty cells / frozen states, and bootstrap percolation (Moore neighbourhood, threshold 3, polluted variants). I didn't find a paper that runs polluted bootstrap growth with mobile Schelling blockers. Since the outcome follows directly from the rules, that's a missing rule combination, not missing behaviour.

---
