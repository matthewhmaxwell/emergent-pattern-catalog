# planted:N19

## Target paper (ground truth; the gate never saw this)
- Title: Multi-scale patterns arising from competing chemotactic interactions
- arXiv: 2610.32591 (v1 2026-10-23); authors: Agniva Datta, Robert Großmann, Carsten Beta
- Behaviour the authors report as new: A Keller-Segel model with chemotaxis toward a self-made attractant and a consumed nutrient produces a dense aggregate coexisting with travelling bands of cells.
- Authors' novelty statement: "We show that this can give rise to diverse patterns including long wavelength uninterrupted coarsening, Turing-like arrested coarsening, traveling waves as well as a novel class of mixed multi-scale oscillatory coarsening patterns, where a dense aggregate coexist with a dynamically spreading part of the population."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: three continuous fields on a 50x50 periodic square: cell density b, an attractant a made by the cells, and a nutrient n eaten by the cells; all start uniform at 1 plus small noise. Rules (dimensionless): (1) b diffuses and drifts up the gradient of a with strength 5*b/(1+b^2) and up the gradient of n with strength 14*b/(1+b^2); cell number is conserved; (2) a diffuses 2.5 times faster than cells, is produced at rate b and decays at rate 1; (3) n diffuses 0.75 times as fast as cells, is consumed at rate 0.02*b*n and resupplied at rate 0.02 everywhere.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  (none distinct)

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - coexisting_structure_types = 2  (kinds of structure coexisting: one dense blob and travelling bands)
  - dense_aggregate_count = 1  (dense aggregates sitting in or carried by the bands)
  - attractant_to_nutrient_strength_ratio = 0.36  (attractant/nutrient drift-strength ratio where this appears (5/14))

This behaviour appears across 1 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Chemotactic swarm rings/travelling bands with aggregates (Budrene-Berg patterns; motile spots guided by nutrient-driven rings)
- **behaviour_documented**: True
- **notes**: Twist: no growth, uniform nutrient resupply, periodic domain, uniform noisy start, so the bands appear spontaneously rather than spreading from an inoculum. Same behaviour otherwise. The nutrient/attractant ratio dependence of moving spots carried by rings is documented by Marenduzzo et al. Possible transient on the way to merging; long run times not checked. Only snippets of the sources were seen.
- **citations**: ["Brenner, Levitov & Budrene, Biophys J 1998 (S0006-3495(98)77880-4)", "Tyson, Lubkin & Murray, J Math Biol 38:359 (1999)", "Tyson, Lubkin & Murray, Biophys J 68:2181 (1995)", "Marenduzzo et al., PLoS ONE 8:e74878 (2013)"]

### N_blocked_r2 (strict pass): verdict KNOWN
- **standard_name**: Coexistence of stationary/motile chemotactic aggregate(s) and nutrient-driven travelling pulses/bands in Keller-Segel models with autochemotaxis plus nutrient taxis (Budrene-Berg type)
- **behaviour_documented**: True
- **notes**: Mild overturn: first pass said VARIANT, but the phenomenology (one aggregate coexisting with travelling nutrient-depletion bands, appearing when nutrient taxis dominates) is explicitly documented in 1D (PLoS Comput Biol 2010) and in 2D (Curk 2013, including coarsening to a single cluster). Remaining differences (uniform nutrient resupply, no growth, b/(1+b^2) sensitivity, periodic domain) are rule-level; resupply mainly lets bands recur instead of being transient. Saturating sensitivity prevents blow-up, so the pattern is likely genuine, but the metrics are sparse: check that bands keep moving rather than the run relaxing to a frozen single-blob state.
- **citations**: ["Mathematical Description of Bacterial Traveling Pulses, PLoS Comput Biol 2010, doi:10.1371/journal.pcbi.1000890 (authors from memory: Saragosti, Calvez, Bournaveas, Buguin, Silberzan, Perthame) - 'Coexistence of a stationary cluster and a traveling pulse'", "Curk et al., Chemotactic Sensing towards Ambient and Secreted Attractant Drives Collective Behaviour of E. coli, PLoS ONE 2013, doi:10.1371/journal.pone.0074878 - 2D, motile vs fixed spots set by nutrient/attractant sensitivity ratio, coarsening to single cluster", "Brenner, Levitov & Budrene, Physical Mechanisms for Chemotactic Pattern Formation by Bacteria, Biophys J 1998 - swarm rings plus aggregates"]

### N_blocked_r3 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Keller-Segel attractant-secretion + nutrient-taxis patterning: depletion-driven travelling swarm bands with chemotactic-collapse aggregates (Brenner-Levitov-Budrene 1998 mechanism; no-growth conserved variant as in Saragosti et al. 2011)
- **behaviour_documented**: True
- **notes**: Twist not found described: bands appear on their own from a uniform start on a periodic domain with nutrient topped up everywhere (literature bands start from an inoculum or imposed gradient), and a single bounded blob (bounded because b/(1+b^2) sensitivity weakens at high density) persists alongside the bands instead of the bands breaking into an aggregate array. Both look like expected consequences of known mechanisms, not new mechanisms; the 5/14 ratio is a parameter value, not a behaviour. Artifact risk: one blob on a 50x50 periodic grid is the usual Keller-Segel coarsening endpoint and the bands may be a long transient; check long-time persistence, grid and time-step refinement, a larger domain, and conservation of total b. Primary sources were read via abstracts and review excerpts only.
- **citations**: ["Brenner, Levitov, Budrene, Biophys J 74:1677 (1998) doi:10.1016/S0006-3495(98)77880-4", "Xue & Othmer, Multiscale models of taxis-driven patterning in bacterial populations (review), PMC2752049", "Saragosti et al., PNAS 2011, doi:10.1073/pnas.1101996108", "Emako et al., Traveling pulses for a two-species chemotaxis model, PLoS Comput Biol 2016", "Painter & Hillen, Spatio-temporal chaos in a chemotaxis model, Physica D 240:363 (2011)"]
