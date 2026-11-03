# planted:N19

**Reported behavior (what the authors report as new):** A Keller-Segel model with chemotaxis toward a self-made attractant and a consumed nutrient produces a dense aggregate coexisting with travelling bands of cells.

**Rule on the card:** World: three continuous fields on a 50x50 periodic square: cell density b, an attractant a made by the cells, and a nutrient n eaten by the cells; all start uniform at 1 plus small noise. Rules (dimensionless): (1) b diffuses and drifts up the gradient of a with strength 5*b/(1+b^2) and up the gradient of n with strength 14*b/(1+b^2); cell number is conserved; (2) a diffuses 2.5 times faster than cells, is produced at rate b and decays at rate 1; (3) n diffuses 0.75 times as fast as cells, is consumed at rate 0.02*b*n and resupplied at rate 0.02 everywhere.

**Observed features on the card:**
  - coexisting_structure_types = 2 (kinds of structure coexisting: one dense blob and travelling bands)
  - dense_aggregate_count = 1 (dense aggregates sitting in or carried by the bands)
  - attractant_to_nutrient_strength_ratio = 0.36 (attractant/nutrient drift-strength ratio where this appears (5/14))

## Reply R08

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Keller-Segel attractant-secretion + nutrient-taxis patterning: depletion-driven travelling swarm bands with chemotactic-collapse aggregates (Brenner-Levitov-Budrene 1998 mechanism; no-growth conserved variant as in Saragosti et al. 2011)", "notes": "Twist not found described: bands appear on their own from a uniform start on a periodic domain with nutrient topped up everywhere (literature bands start from an inoculum or imposed gradient), and a single bounded blob (bounded because b/(1+b^2) sensitivity weakens at high density) persists alongside the bands instead of the bands breaking into an aggregate array. Both look like expected consequences of known mechanisms, not new mechanisms; the 5/14 ratio is a parameter value, not a behaviour. Artifact risk: one blob on a 50x50 periodic grid is the usual Keller-Segel coarsening endpoint and the bands may be a long transient; check long-time persistence, grid and time-step refinement, a larger domain, and conservation of total b. Primary sources were read via abstracts and review excerpts only."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium confidence. This agrees with the first pass's family, but I'd describe the mechanism more precisely and add one caveat.**

**What's documented**
- **Same ingredients, same structures.** Brenner, Levitov & Budrene (1998) model bacteria that make their own attractant while using up a substrate. They report that through chemoattractant signaling, the bacteria organize into swarm rings and aggregates. The ring motion is caused by the depletion of a substrate that is necessary to produce attractant. Aggregates form when the mass density of the swarm ring exceeds a threshold, the ring collapses cylindrically and then destabilizes into aggregates. A later review sums up the split: the ring moves because nutrient runs out locally and holds together because attractant is high there, while the aggregates formed in the ring results from fluctuations near the unstable uniform cell density. That is the same division of labour as your 14 (nutrient) vs 5 (attractant) drift strengths.
- **No cell growth, cell number conserved.** This setup is also studied. The Saragosti/Calvez line of work models bacteria that produce the same chemoattractant which gathers the population and at the same time they consume the common nutrient which triggers the motion, and finds robust and conserved travelling waves propagating at a constant velocity.
- **Several aggregates collapsing to one.** Keller–Segel models are known to produce aggregates that merge, and sometimes new ones that appear (Painter & Hillen 2011). In the absence of the production terms, the spikes in the original KS model either merge or drift. In a conserved system with no growth, ending up with a single blob is the expected outcome.

**The twist I couldn't find described (phenomenology, not rules)**
1. The bands appear on their own from a uniform start on a periodic domain where nutrient is topped up everywhere. In the literature they are launched from an inoculum or an imposed gradient.
2. One dense blob persists alongside the travelling bands, rather than the bands breaking up into an array of aggregates.

Does this matter? Probably only a little. Nutrient being topped up everywhere is what lets depletion waves keep re-forming. The b/(1+b²) sensitivity weakens at high density, so the blob can't collapse to a point, which keeps it finite. Both are expected consequences of known mechanisms. Neither looks like a new mechanism. The 5/14 ratio is a parameter value, not a behaviour.

**Genuine pattern or artifact?** It doesn't look like noise, a field that does nothing, or a frozen state: there are two kinds of moving structure and the nutrient coupling is clearly active. However, on a 50×50 periodic grid, "one blob" is the standard Keller–Segel coarsening endpoint, and the bands may be a long transient. Before treating the coexistence as a lasting regime, check:
- that it persists over much longer runs;
- that it survives a finer grid and smaller time step;
- that it holds on a larger domain;
- that the total amount of b stays constant.

**Limits of this check:** I read summaries, abstracts and review excerpts, not the full Brenner et al. or Keller–Segel 1971 papers. I found no paper that specifically shows one aggregate coexisting with spontaneously repeating depletion bands under uniform nutrient top-up.

**

## Reply R15

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Budrene-Berg-type chemotactic pattern formation: depletion-driven travelling swarm bands plus self-attractant aggregates (Keller-Segel with excreted attractant and consumed nutrient)", "notes": "Aggregates, depletion-driven moving bands and their link are documented. Twist not found: travelling bands arising spontaneously from a noisy uniform state with no growth and constant nutrient resupply, coexisting with one long-lived dense blob on a periodic 2D domain, rather than an inoculum ring breaking into spots. Plausibly a delayed effective self-repulsion from nutrient depletion competing with fast self-attraction (my reasoning, not from a source). The 0.36 ratio is an input parameter, not an emergent measurement. Risk: the single blob may be a long transient of vanishing-sensitivity coarsening on a small 50x50 grid; extend the run and check on finer and larger grids."}

**Reply text:**

**Verdict: KNOWN-VARIANT. This agrees with the first pass; I found no reason to overturn it.** The phenomena involved are well documented. What I could not find is a paper showing this exact combination arising spontaneously from a uniform start.

**What is documented**
- **Same model family:** Brenner, Levitov and Budrene couple cell movement to equations for both an attractant and a nutrient, and their model leads to analytic solutions corresponding to both rings and aggregates.
- **Why the bands move:** the movement of the swarm rings is due to the depletion of succinate around the band of bacteria: chemotaxis holds the band together, but the net motion is not caused directly by chemotactic fluxes. This matches your rule set, where cells consume the nutrient and climb its gradient.
- **Prior simulations:** they consistently show (one-dimensional) traveling rings destabilizing via transverse instability into spots of bacteria. Tyson, Lubkin and Murray give a mathematical model which captures all three of the pattern forming processes seen in E. coli and Salmonella.
- **No growth needed:** other work notes that more complex nutrient terms can give rise to travelling waves, even when growth is absent.
- **Moving patterns with two signals:** Painter (2009) studied cells responding to two chemical signals on a periodic domain. In 2D, these temporal patterns take the form of a travelling aggregation.
- **Your sensitivity function:** b/(1+b²) shrinks at high density. This puts it close to models with a vanishing chemotaxis coefficient at high population densities, which show spatio-temporal patterns which allow for ultra-long transients and merging or coarsening.

**The twist I could not find described**
- In the classic work, rings start from a small inoculum, travel outward and break up into spots. Here there is no inoculum: the noisy uniform state breaks up by itself into travelling bands, and these coexist with a single long-lived dense blob on a periodic domain.
- A plausible mechanism, from my own reasoning rather than a source: when the nutrient is depleted and slowly resupplied, climbing its gradient pushes cells away from crowded areas, but with a delay. That delayed push competes with the fast self-attraction and could drive the moving patterns. It fits known attraction–repulsion and two-signal behaviour.
- This twist matters only modestly. It is a different setting for familiar structures, not a new kind of structure.

**Real pattern or artifact?** Probably real, not trivial: the drift strengths are large, and the density-limited sensitivity prevents the blow-up that would otherwise make the blob a numerical artifact. Three caveats:
- On a 50×50 grid with vanishing high-density sensitivity, "one blob plus bands" may be a long transient on the way to a fully merged state, which the vanishing-sensitivity results above predict.
- Two structure types and one aggregate are low numbers and could be a single snapshot.
- The 0.36 attractant/nutrient ratio is just the input parameter 5/14, not something the simulation produced.

The run should be extended and checked on a finer grid and a larger domain before reading anything into it.

**Search limits:** I read abstracts and fragments of the papers above, not full texts. I ran six searches: coexistence of aggregates and bands with attractant plus nutrient; Brenner–Levitov–Budrene; two-signal periodic-domain patterns; the Tyson–Lubkin–Murray models; attraction–repulsion travelling waves; and Painter–Hillen spatio-temporal patterns.

## Reply R24

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Chemotactic aggregate formation plus nutrient-depletion-driven travelling bands/swarm rings (Budrene-Berg / Brenner-Levitov-Budrene / Tyson-Lubkin-Murray patterns)", "notes": "No source found that reproduces this exact setup (periodic 50x50, no growth, constant nutrient resupply, b/(1+b^2) sensitivity for both signals, uniform noisy start). The aggregate-plus-band coexistence in 2D is expected from the literature but not directly confirmed there. I read only abstracts and fragments of the Woodward and Tyson papers. The Hillen & Painter review was seen only as abstract and references."}

**Reply text:**

**Model family.** This is a continuum chemotaxis system of the Keller–Segel type, from bacterial pattern-formation modelling. It has three fields: cells with a receptor-saturating sensitivity b/(1+b²), a self-produced attractant, and a consumed nutrient. Cells are attracted to both the attractant and the nutrient, which fits the usual picture of swarm rings and aggregates in E. coli and Salmonella.

**Is the behaviour already described? Yes, in the same model class.** The observed state is a dense aggregate together with travelling bands, in a system with both attractant taxis and nutrient taxis. Four sources cover most of this:

- **Woodward et al. 1995 (Biophys J 68:2181).** The Salmonella model couples cells, a self-excreted attractant and a nutrient. It produces swarm rings and aggregates, with a saturating chemotactic sensitivity. I could only read fragments of the paper in the search results.
- **Brenner, Levitov & Budrene 1998 (Biophys J 74:1677–1693).** Through chemoattractant signalling, the bacteria organize into swarm rings and aggregates. Ring motion comes from depletion of a nutrient around the band; chemotaxis holds the band together. Rings can collapse into aggregates when their mass exceeds a threshold. Their stated condition is that aggregates can form as long as the bacteria produce attractant.
- **Tyson, Lubkin & Murray 1999 (Proc R Soc B 266:299–304).** Their model captures all three of the pattern-forming processes seen experimentally in E. coli and S. typhimurium. The search result gave only a fragment of the abstract.
- **Closest match in the search results: "Chemotactic Sensing towards Ambient and Secreted Attractant Drives Collective Behaviour of E. coli" (PLoS ONE 2013, DOI 10.1371/journal.pone.0074878).** Dynamics and steady-state patterns depend strongly on the ratio between the efficiencies of chemotaxis to nutrient and to attractant gradients. For slightly stronger nutrient sensing, spots still form but are mobile and glide up the nutrient gradient. That is the same kind of result as an attractant/nutrient strength ratio selecting a mixed aggregate-plus-band state. It is a different, agar-plate, ring-geometry model. I did not see its author list, so I leave it uncredited.

**Coexistence of a stationary cluster and a travelling pulse** is also reported in a related chemotaxis model. Saragosti, Calvez, Bournaveas, Buguin, Silberzan & Perthame 2010 (PLoS Comput Biol 6(8):e1000890) say a fraction remains trapped at the boundary as a stationary profile and a fraction travels across the channel at constant speed. That model is kinetic and one-dimensional, so it is only a loose analogue.

**Reviews.** The best review I found is Hillen & Painter 2009, "A user's guide to PDE models for chemotaxis" (J Math Biol 57). It reviews the formulation of these models, contrasts their patterning properties, and classifies their solution form. I only saw its abstract and reference list. I did not find it describing this exact coexistence.

**Verdict: KNOWN-VARIANT, leaning KNOWN.** The mechanism is established. Nutrient depletion drives and carries bands, self-attractant production collapses density into aggregates, and the attractant-to-nutrient taxis ratio selects the pattern.

What I did not find is a paper that simulates this exact system on a periodic 2D square, from a uniform start with noise. It uses a nutrient with constant resupply, no cell growth and b/(1+b²) sensitivity for both signals. In 2D, "a blob plus travelling bands" is a mixed pattern that the Brenner/Woodward/Tyson models would be expected to produce for some parameters. I did not verify it as an explicit 2D result. The differences I see (no growth, a periodic box, conserved cell number) most likely change where the pattern appears in parameter space, not whether it exists.

**Confidence: medium.** The mechanism and the ratio dependence are well supported. I could read only abstracts and fragments of the original papers, and I found no source reproducing the exact parameters.

## Reply R45

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Attractant-nutrient chemotaxis model (Tyson-Lubkin-Murray / Woodward et al. type); aggregates plus travelling nutrient-taxis bands", "notes": "Ingredients (nutrient-taxis bands, autochemotaxis aggregates, two-field attractant+nutrient models) are all documented. No source found reporting the specific 2D coexistence of a dense blob and travelling bands in a growth-free, nutrient-resupplied system. Did not read the Tyson equations, so the b/(1+b^2) sensitivity and parameter-regime match are unverified."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium-low confidence.** The model is a well-studied family, but I found no source that reports this exact coexistence.

**1. Model family**
This is a Keller–Segel-type chemotaxis system (a continuum reaction–diffusion–advection PDE). It has three fields:
- Cells (b) are conserved and drift up two gradients with a saturating, density-dependent sensitivity.
- An attractant (a) is secreted by the cells, so the cells attract themselves (autochemotaxis).
- A nutrient (n) is consumed by the cells, resupplied everywhere, and also attracts them.

The closest named relatives are the "attractant + nutrient" bacterial pattern models of Tyson, Lubkin and Murray. These were built to describe *E. coli* and *Salmonella* patterns. The sources I read describe such models in words and do not print the equations. I did not confirm that the b/(1+b²) sensitivity form is the one those papers use.

**2. What the literature says about the specific behaviour**
- **Travelling bands come from nutrient taxis.** Cells that consume a nutrient and climb its gradient form travelling bands or rings. In Adler's capillary assay, E. coli cells move up a nutrient gradient and the population forms moving bands or rings. The original theory is Keller and Segel's. They recovered a traveling-wave solution from reaction–diffusion–convection equations for bacterial density and the attractant.
- **Aggregates come from autochemotaxis.** Chemotactic *S. typhimurium* aggregate in response to gradients of attractants they excrete themselves. This is the Woodward et al. 1995 paper (Biophys J 68:2181–2189), which pairs experiments with a cell–chemotaxis model.
- **The Tyson models use two coupled fields.** Two review sources describe them as coupling the PKS chemotaxis equation with reaction-diffusion equations for both the attractant and nutrient, proposing a minimal mechanism for swarm rings and aggregates. The citations are Tyson, Lubkin and Murray 1999 (Proc R Soc B 266:299–304; J Math Biol 38:359–375) and Tyson, Stern and LeVeque 2000 (J Math Biol 41:455–475). I verified the titles, authors, years and venues in bibliographic listings, but I did not read the full texts.
- **Both taxis types appear together in other models.** A PLoS ONE 2013 paper, "Chemotactic Sensing towards Ambient and Secreted Attractant Drives Collective Behaviour of E. coli", combines both. It reports swarming rings when cells respond to the nutrient gradient only, and droplet formation then coarsening to a single stable spot from the secreted attractant. In that paper the two kinds of structure appear in different strains or conditions, not coexisting.
- **Two-species travelling pulses are known in 1D.** A PLoS Comput Biol two-species model has bacteria reacting to both a chemoattractant and a nutrient. This is the paper "Traveling Pulses for a Two-Species Chemotaxis Model", and it is 1D.
- **Coarsening of aggregates is general.** Mass-conserving Keller–Segel models coarsen. A PRE paper, "Coarsening dynamics of chemotactic aggregates", states that coarsening is driven by self-amplifying mass transport and aggregate coalescence. I only saw this paper's abstract and reference list, so I could not check its author list.
- **A recent review.** Hillen and Painter 2009, "A user's guide to PDE models for chemotaxis" (J Math Biol 58:183–217), is the standard review of this model class. I verified its title and journal only. For the pattern-formation side, Painter and Hillen 2011 (Physica D 240:363–375) covers merging and emerging aggregates, but only for 1D with logistic growth.

**3. What differs from the observed behaviour**
I found no paper that reports a 2D, cell-conserving, nutrient-resupplied system in which a dense blob coexists with travelling bands that carry aggregates. The differences from the published models are:
- **No cell growth.** The classic band and ring models rely on growth or on a singular sensitivity. This system conserves cell number.
- **Constant nutrient resupply.** Resupply is uniform, not a depleting reservoir.
- **Strong nutrient drift.** The nutrient drift is stronger than the attractant drift (ratio 0.36).
- **2D with a saturating sensitivity.** Most analyses of travelling waves are 1D, and I could not confirm that the published models use this sensitivity form.

These differences probably matter for the exact parameter regime. I expect the mechanism itself to be standard: nutrient-driven travelling fronts plus self-attraction collapse, followed by coarsening.

**4. Confidence and caveats**
- **Why medium-low.** I could not read the Tyson papers' equations, so I cannot say whether their sensitivity and parameter ranges match this system. My search budget was limited, and I did not find a source describing the coexistence directly.
- **Why not LITERATURE-SILENT.** The ingredients and the model family are clearly documented.

**

## Reply R51

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Chemotactic swarm rings/travelling bands with aggregates (Budrene-Berg patterns; motile spots guided by nutrient-driven rings)", "notes": "Twist: no growth, uniform nutrient resupply, periodic domain, uniform noisy start, so the bands appear spontaneously rather than spreading from an inoculum. Same behaviour otherwise. The nutrient/attractant ratio dependence of moving spots carried by rings is documented by Marenduzzo et al. Possible transient on the way to merging; long run times not checked. Only snippets of the sources were seen."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium confidence. I agree with the first pass.**

**What is documented**

- **Aggregates and travelling rings from the same ingredients.** Brenner, Levitov & Budrene (Biophys. J. 1998) say the model leads to analytic solutions corresponding to both rings and aggregates. They report that aggregate formation and ring instability depend only on coarse properties, such as chemotaxis and production of attractant by bacteria. In their mechanism, swarm ring migration is caused by the depletion of a chemical used in the process of attractant production, and a secondary instability causes the breakdown into aggregates.
- **Tyson, Lubkin & Murray** studied the same three-field setup (cells, nutrient, attractant) as a nonlinear partial differential equation system.
- **Closest match: a dense blob coexisting with bands, controlled by the attractant-to-nutrient ratio.** Marenduzzo and co-authors (PLoS ONE 2013) report that both the dynamics and the steady-state patterns are therefore much affected by the ratio between the efficiency of chemotaxis to nutrient or chemoattractant gradients. When nutrient sensing is strong enough, the bacterial spots become motile and the eventual pattern is a ring made up or merging moving spots rather than an amorphous or regular array of static ones. They also note that the positions of the spots is guided by the swarming ring. That is essentially the observed behaviour: aggregates sitting in or carried by nutrient-driven bands, at a high nutrient-to-attractant ratio (here 14/5, or 0.36 the other way round).
- **Nutrient-depletion bands in general.** Self-generated travelling bands go back to Keller and Segel, who found a traveling-wave solution using a pair of reaction–diffusion–convection equations.

**The twist, and whether it matters**

In the papers above, cells grow and the rings spread outward from a central inoculum. Here, cell number is conserved, nutrient is resupplied everywhere, and the domain is periodic with a uniform noisy start. So the bands appear spontaneously from a Turing-like instability rather than from an expanding front. That is a real difference in setup. But the behaviour (nutrient-taxis bands carrying or coexisting with an attractant-driven aggregate) is the same. It probably matters little unless the bands turn out to keep going indefinitely rather than being transient.

**Trivial or artifact?**

It is probably a genuine pattern, not an artifact. The b/(1+b²) sensitivity weakens at high density, so the aggregate stays finite instead of blowing up numerically, and the nutrient field is clearly active. Three cautions:
- The metrics are thin: one aggregate and two structure types.
- On a 50×50 grid this could be a passing stage on the way to merging into a single aggregate.
- It needs a check over long run times.

**What I couldn't check:** I only saw snippets of these papers, not full texts. None of them describes this exact setup (no growth, uniform resupply, periodic box, uniform start).
