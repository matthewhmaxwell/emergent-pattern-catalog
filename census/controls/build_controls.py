"""Positive and negative controls for the literature gate (POST-REGISTRATION, EXPLORATORY; 2026-10-26, after review).

Writes controls.json (40 cards in the exact census card schema) and truth.json (kept out of every prompt).

KNOWN controls (K01-K20): real, documented phenomena that are NOT in the census reference library, written in the
census card style (plain-language rule system, parameter variants, standout fingerprint features). Truth: KNOWN.

FABRICATED controls (F01-F20): well-documented rule systems paired with observed features that CONTRADICT the
system's documented behaviour. The reported behaviour is therefore undocumented for that system. Truth: not KNOWN.
These test whether the gate checks the observed behaviour or only recognises the model family. A gate that answers
KNOWN by citing the system's standard behaviour would also label a genuinely new behaviour in a familiar-looking
system as known. Scoring rules are fixed in SCORING.md before any card is run.
"""
import json, os

L64 = "World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). "
AG = "World: 400 self-propelled agents moving in a wrap-around 20x20 square. "
NET = "World: a random network of 400 nodes (mean degree 6). "

K = [
 ("K01", "grid", "Langton's ant: highway after a chaotic transient (Langton 1986)",
  "World: a 100x100 lattice of cells (wrap-around edges), all cells white at the start, plus one mobile walker on a "
  "cell with a heading. Rules each step: (1) on a white cell the walker turns 90 degrees right, on a black cell 90 degrees "
  "left; (2) it flips the colour of the cell it is on; (3) it moves forward one cell.",
  ["Same rules on a 200x200 lattice.", "Same rules, walker starts on a random black-and-white lattice."],
  [("g_spec_peak", 6.1), ("g_moran_trend", 2.4), ("g_corrlen", 3.0), ("g_activity_cv", 2.2)]),
 ("K02", "grid", "Conway's Game of Life: decay to still lifes and oscillators (ash), with gliders",
  L64 + "Cells have 2 types (dead/alive; 30% alive at random start). Rules each step, all cells at once: a dead cell "
  "with exactly 3 alive neighbours becomes alive; an alive cell with 2 or 3 alive neighbours stays alive; every other "
  "cell becomes dead.",
  ["Same rules, 10% alive at start.", "Same rules on a 128x128 lattice."],
  [("g_activity_trend", -3.5), ("g_n_domains", 6.2), ("g_activity_cv", 2.8), ("g_reversals", 4.1), ("g_spec_peak", 3.0)]),
 ("K03", "grid", "Brian's Brain: abundant gliders/spaceships (Silverman)",
  L64 + "Cells have 3 types (ready, firing, refractory; random start). Rules each step, all cells at once: a ready "
  "cell fires if exactly 2 of its 8 neighbours are firing; a firing cell becomes refractory; a refractory cell becomes ready.",
  ["Same rules, 5% firing at start."],
  [("g_cyclic_changes", 7.0), ("g_activity", 3.9), ("g_activity_cv", -1.2), ("g_n_domains", 5.6), ("g_moran", 0.3)]),
 ("K04", "grid", "Drossel-Schwabl forest-fire model: self-organized critical fire sizes, spiral fire fronts",
  L64 + "Cells have 3 types (empty, tree, burning; random start). Rules each step: (1) a burning cell becomes empty; "
  "(2) a tree with a burning neighbour starts burning; (3) a tree catches fire spontaneously with probability 0.0001; "
  "(4) an empty cell grows a tree with probability 0.01.",
  ["Same rules, spontaneous fire probability 0.00005.", "Same rules, tree growth probability 0.02."],
  [("g_activity_cv", 6.0), ("g_cyclic_changes", 5.1), ("g_corrlen", 2.7), ("g_n_domains", 3.8)]),
 ("K05", "grid", "Bak-Sneppen model: punctuated equilibrium, self-organized fitness threshold near 0.667",
  "World: a ring of 200 cells (each with 2 neighbours). Each cell holds a number between 0 and 1 (uniform random start). "
  "Rules each step: the cell with the lowest number and its two neighbours each get a new uniform random number.",
  ["Same rules on a ring of 500 cells."],
  [("g_activity_cv", 4.5), ("g_corrlen", 2.0), ("g_moran_trend", 1.6)]),
 ("K06", "grid", "Axelrod model of cultural dissemination: frozen cultural regions",
  "World: a 32x32 lattice of cells (no wrap-around, 4 neighbours each). Each cell holds 5 features, each with one of 10 "
  "traits (random start). Rules each step: a random cell picks a random neighbour; with probability equal to the "
  "fraction of features they share, it copies one feature on which they differ.",
  ["Same rules, 15 traits per feature.", "Same rules, 3 features."],
  [("g_n_domains", 5.0), ("g_largest_domain", 2.1), ("g_activity_trend", -4.0), ("g_moran", 4.4)]),
 ("K07", "network", "Deffuant bounded-confidence model: opinions settle into a few separated clusters",
  "World: 500 nodes, every pair can meet (well mixed). Each node holds an opinion between 0 and 1 (uniform random start). "
  "Rules each step: a random pair meets; if their opinions differ by less than 0.2, each moves halfway toward the other.",
  ["Same rules, threshold 0.3.", "Same rules, each moves a quarter of the way."],
  [("n_consensus_gain", 3.1), ("n_type_entropy", -2.2), ("n_oscillation", -1.0)]),
 ("K08", "grid", "Nagel-Schreckenberg traffic model: spontaneous jams and backward-moving jam waves",
  "World: a ring road of 1000 cells with 300 cars, each with a speed from 0 to 5. Rules each step, all cars at once: "
  "(1) speed up by 1 if below 5; (2) slow to the number of empty cells ahead if that is smaller; (3) with probability "
  "0.3 slow down by 1; (4) move forward by the speed.",
  ["Same rules, 200 cars.", "Same rules, braking probability 0.1."],
  [("g_activity_cv", 3.0), ("g_n_domains", 4.2), ("g_moran", 3.6), ("g_spec_peak", 2.5)]),
 ("K09", "grid", "Daisyworld: planetary temperature regulation by daisy albedo (Watson & Lovelock 1983)",
  L64 + "Cells have 3 types (bare ground, dark daisy, light daisy; random start). Each cell's local temperature is set "
  "by incoming sunlight and the albedo of its own type and its neighbours (dark absorbs more). Rules each step: a daisy "
  "dies with probability 0.3; a bare cell next to a daisy becomes that daisy type with probability given by a growth "
  "curve that peaks at 22.5 C; sunlight increases slowly over the run.",
  ["Same rules, faster increase of sunlight."],
  [("g_n_types", 3.0), ("g_moran", 3.0), ("g_largest_domain", 1.5), ("g_type_entropy", 1.4)]),
 ("K10", "grid", "Contact process near its critical point: absorbing-state transition, directed-percolation class",
  "World: a 128x128 lattice of cells (wrap-around edges, 4 neighbours each). Cells have 2 types (healthy, infected; "
  "half infected at start). Rules, in continuous time: an infected cell recovers at rate 1 and infects each healthy "
  "neighbour at rate 0.41 (total infection rate 1.64).",
  ["Same rules, infection rate 0.40 per neighbour.", "Same rules on a 256x256 lattice."],
  [("g_activity_cv", 4.0), ("g_corrlen", 5.5), ("g_activity_trend", -1.8), ("g_n_domains", 3.0)]),
 ("K11", "network-phase", "Explosive synchronization: degree-frequency correlation on a scale-free network (Gomez-Gardenes et al. 2011)",
  "World: a scale-free network of 1000 nodes (preferential attachment, 3 links per new node). Each node is a phase "
  "oscillator whose natural frequency equals its number of links. Rules each step: each oscillator's phase advances by "
  "its natural frequency plus coupling 1.5 times the mean sine of the phase differences to its neighbours.",
  ["Same rules, coupling 1.2.", "Same rules, coupling 2.0."],
  [("p_r", 4.5), ("p_r_std", 3.8), ("p_freq_spread", 5.0)]),
 ("K12", "grid", "Rock-paper-scissors with mobility: spirals, and loss of coexistence above a critical mobility (Reichenbach, Mobilia & Frey 2007)",
  "World: a 128x128 lattice of cells (wrap-around edges, 4 neighbours each). Cells have 4 types (empty and species A, B, "
  "C; random start). Rules, in continuous time with random neighbour pairs: A replaces a neighbouring B with an empty "
  "cell, B does the same to C, C to A (rate 1 each); any species reproduces into a neighbouring empty cell (rate 1); "
  "any two neighbours swap places at rate 50.",
  ["Same rules, swap rate 5 (low mobility)."],
  [("g_cyclic_changes", 4.0), ("g_n_types", 1.0), ("g_largest_domain", 6.0), ("g_activity_trend", -3.0)]),
 ("K13", "network", "Cooperation promoted on scale-free networks (Santos & Pacheco 2005)",
  "World: a scale-free network of 1000 nodes (preferential attachment, 2 links per new node). Node types: cooperator or "
  "defector (half each at random). Rules each step: every node plays a prisoner's dilemma with each neighbour "
  "(temptation 1.6); then each node copies a random neighbour's type with probability proportional to how much higher "
  "that neighbour's total payoff is.",
  ["Same rules, temptation 1.8.", "Same rules on a random network of the same mean degree."],
  [("n_consensus_gain", 2.5), ("n_degree_cv", 4.0), ("n_type_entropy", -2.5)]),
 ("K14", "agents", "Swarmalators: agents that both swarm and synchronize (O'Keeffe, Ha & Strogatz 2017)",
  "World: 500 agents in the plane, each with a position and an internal phase. Rules each step: each agent moves toward "
  "others with an attraction that is stronger when their phases are similar (strength 1) and away from very close "
  "agents; each phase is pulled toward neighbours' phases with coupling -0.5, weighted by inverse distance.",
  ["Same rules, phase coupling -0.75."],
  [("a_milling", 3.2), ("p_local_r", 3.5), ("p_r", 0.1), ("a_nn_ratio", 0.4)]),
 ("K15", "grid", "Ant brood sorting / clustering (Deneubourg et al. 1991)",
  "World: a 64x64 lattice of cells (wrap-around edges) with 400 scattered items of 2 kinds and 20 mobile walkers. Rules "
  "each step: each walker steps to a random neighbouring cell; an unladen walker picks up an item with probability "
  "(0.1/(0.1+f))^2, and a laden walker drops its item with probability (f/(0.3+f))^2, where f is the fraction of items "
  "of the same kind it met in its last 10 steps.",
  ["Same rules, 1 kind of item.", "Same rules, 40 walkers."],
  [("g_moran", 5.5), ("g_n_domains", -2.0), ("g_largest_domain", 2.0), ("g_moran_trend", 3.1)]),
 ("K16", "agents", "Couzin zonal model: milling (torus) state (Couzin et al. 2002)",
  AG + "Rules each step: each agent turns away from agents within distance 1; otherwise it aligns with agents between 1 "
  "and 2 and moves toward agents between 2 and 15, ignoring agents in a 90-degree blind zone behind it; speed 3, "
  "turning rate limited to 40 degrees per step, heading noise 0.05 rad.",
  ["Same rules, alignment zone up to 2.5."],
  [("a_milling", 6.5), ("a_polar", 0.1), ("a_nn_ratio", 0.3)]),
 ("K17", "lattice-phase", "Chimera state: coexisting coherent and incoherent domains in identical oscillators (Kuramoto & Battogtokh 2002)",
  "World: a ring of 256 identical phase oscillators. Rules each step: each phase advances by 1 plus a coupling to every "
  "other oscillator, weighted by exp(-4 x distance along the ring) and shifted by a phase lag of 1.46 rad.",
  ["Same rules, phase lag 1.40.", "Same rules, 512 oscillators."],
  [("p_local_r", 4.0), ("p_r_std", 2.0), ("p_r", 1.2)]),
 ("K18", "field", "Lenia: self-organizing moving localized patterns (gliders such as Orbium; Chan 2019)",
  "World: a 128x128 field (wrap-around edges) holding a value between 0 and 1 per cell (random patch at start). Rules "
  "each step: each cell's value changes by 0.1 times a growth function of the ring-weighted average of values within "
  "radius 13 (growth peaks when that average is 0.15, width 0.015), clipped to [0, 1].",
  ["Same rules, growth peak 0.14."],
  [("f_moran", 5.0), ("f_peak_sharp", 4.5), ("g_activity", 1.8)]),
 ("K19", "grid", "Elementary cellular automaton rule 110: localized gliders on a periodic background (Cook 2004)",
  "World: a ring of 400 cells (2 neighbours each), 2 types (random start). Rules each step, all cells at once: a cell's "
  "new type is given by Wolfram's elementary rule 110 applied to itself and its two neighbours.",
  ["Same rule on a ring of 1000 cells."],
  [("g_n_domains", 4.0), ("g_activity_cv", 2.5), ("g_reversals", 3.0), ("g_spec_peak", 2.0)]),
 ("K20", "grid", "Sznajd model: consensus through outflow of agreeing pairs (Sznajd-Weron & Sznajd 2000)",
  L64 + "Cells have 2 types (random start). Rules each step: pick a random pair of neighbouring cells; if they have the "
  "same type, all of the pair's other neighbours adopt that type.",
  ["Same rules, 3 types."],
  [("g_largest_domain", 5.0), ("g_activity_trend", -3.0), ("g_moran", 4.0)]),
]

F = [
 ("F01", "grid", "majority rule (asynchronous, 2 types): reported sustained cycling at a dominant frequency, no clustering",
  L64 + "Cells have 2 types (random start). Rules each step: each cell, with probability 0.1, adopts the type held by "
  "the majority of its 8 neighbours.",
  ["Same rules, probability 0.2."],
  [("g_cyclic_changes", 6.0), ("g_spec_peak", 5.5), ("g_moran", -0.5), ("g_n_domains", 7.0), ("g_activity_trend", 2.0)]),
 ("F02", "grid", "voter model (2 types): reported anti-clustering (checkerboard-like), no domain growth",
  L64 + "Cells have 2 types (random start). Rules each step: each cell copies the type of one random neighbour.",
  ["Same rules, 4 neighbours."],
  [("g_moran", -2.0), ("g_largest_domain", -3.0), ("g_n_domains", 6.5), ("g_activity_trend", 0.0)]),
 ("F03", "grid", "Schelling (threshold 0.3, 10% vacancies): reported persistent mixing with rising activity",
  L64 + "Cells have 3 types (type 0 = empty, 10% of cells; random start). Rules each step: an occupied cell whose "
  "same-type share among occupied neighbours is below 0.3 moves to a random empty cell.",
  ["Same rules, threshold 0.4."],
  [("g_moran", -1.5), ("g_type_entropy", 1.0), ("g_activity", 4.0), ("g_activity_trend", 3.0)]),
 ("F04", "network-phase", "Kuramoto, strong coupling on a random network: reported global order below chance, strong fluctuation",
  NET + "Each node is a phase oscillator with natural frequency drawn from a normal distribution (mean 0, sd 1). Rules "
  "each step: each phase advances by its natural frequency plus coupling 5 times the mean sine of the phase "
  "differences to its neighbours.",
  ["Same rules, coupling 4."],
  [("p_r", -2.5), ("p_local_r", -1.0), ("p_r_std", 3.5)]),
 ("F05", "lattice-phase", "uncoupled oscillators: reported global and local phase order and phase spirals",
  "World: a 64x64 lattice of phase oscillators. Natural frequencies drawn from a normal distribution (mean 0, sd 1). "
  "Rules each step: each phase advances by its own natural frequency; there is no coupling between oscillators.",
  ["Same rules on a 32x32 lattice."],
  [("p_r", 5.0), ("p_local_r", 5.5), ("p_vortex_density", 3.0)]),
 ("F06", "agents", "Vicsek, maximal noise: reported strong polar order",
  AG + "Rules each step: each agent takes the average heading of agents within distance 1, adds uniform heading noise "
  "spanning the full circle (2 pi rad), and moves at speed 0.3.",
  ["Same rules, 800 agents."],
  [("a_polar", 6.0), ("a_local_align", 5.0)]),
 ("F07", "agents", "Vicsek, low noise, dense: reported no polar order but strong rotational milling and clumping",
  AG + "Rules each step: each agent takes the average heading of agents within distance 1, adds heading noise of up to "
  "0.1 rad, and moves at speed 0.3.",
  ["Same rules, 800 agents."],
  [("a_polar", -2.0), ("a_milling", 6.0), ("a_nn_ratio", 0.2)]),
 ("F08", "grid", "cyclic CA (4 types, threshold 1): reported activity dying out with one type taking over",
  L64 + "Cells have 4 types (random start), ordered in a cycle. Rules each step, all cells at once: a cell adopts the "
  "next type in the cycle if at least 1 of its 8 neighbours already has that next type.",
  ["Same rules, 5 types."],
  [("g_cyclic_changes", -1.0), ("g_largest_domain", 6.5), ("g_activity_trend", -4.0)]),
 ("F09", "grid", "spatial prisoner's dilemma, low temptation 1.2: reported cooperators driven extinct, frozen lattice",
  L64 + "Cells have 2 types (cooperator, defector; half each at random). Rules each step: every cell plays a prisoner's "
  "dilemma with its 8 neighbours and itself (temptation 1.2, sucker payoff 0); then every cell copies the type of the "
  "highest-scoring cell in its neighbourhood.",
  ["Same rules, temptation 1.3."],
  [("g_n_types", 1.0), ("g_largest_domain", 7.0), ("g_activity_trend", -5.0)]),
 ("F10", "network", "co-evolving voter model, rewiring probability 0.95: reported global consensus in one connected component",
  NET + "Node types: 2 opinions (random start). Rules each step: a random node picks a random neighbour of a different "
  "opinion; with probability 0.95 it cuts that link and links to a random node of its own opinion, otherwise it adopts "
  "the neighbour's opinion.",
  ["Same rules, rewiring probability 0.9."],
  [("n_consensus_gain", 4.0), ("n_giant", 1.0), ("n_components", -2.0)]),
 ("F11", "network", "network voter model: reported periodic oscillation of opinions, no convergence",
  NET + "Node types: 2 opinions (random start). Rules each step: a random node copies the opinion of a random neighbour.",
  ["Same rules, mean degree 10."],
  [("n_oscillation", 5.0), ("n_cyclic_changes", 4.5), ("n_consensus_gain", -2.0)]),
 ("F12", "agents", "type-blind attraction: reported strong sorting of agents by type",
  AG + "Agents have 2 types (random start). Rules each step: each agent turns toward the centre of mass of all agents "
  "within distance 3 regardless of type, with heading noise 0.3 rad, and moves at speed 0.3.",
  ["Same rules, 3 types."],
  [("a_type_segregation", 6.0), ("a_nn_ratio", 0.3)]),
 ("F13", "agents", "pure repulsion: reported strong clumping and polar order",
  AG + "Rules each step: each agent turns directly away from the nearest agent within distance 1.5 (otherwise keeps its "
  "heading), with heading noise 0.3 rad, and moves at speed 0.3.",
  ["Same rules, repulsion distance 2."],
  [("a_nn_ratio", 0.1), ("a_polar", 4.0), ("a_local_align", 3.5)]),
 ("F14", "field", "pure diffusion, no reaction: reported sharp, growing spatial pattern",
  "World: a 64x64 grid of 2 chemical fields (wrap-around edges), random start. Rules each step: each field diffuses with "
  "coefficient 0.1; there are no reactions, sources or sinks.",
  ["Same rules, coefficient 0.05."],
  [("f_moran", 5.0), ("f_peak_sharp", 6.0)]),
 ("F15", "grid", "majority rule (3 types): reported cells rotating through all types at a dominant period",
  L64 + "Cells have 3 types (random start). Rules each step: each cell adopts the type held by the plurality of its 8 "
  "neighbours (ties keep the current type).",
  ["Same rules, 4 types."],
  [("g_cyclic_changes", 5.5), ("g_spec_peak", 6.0), ("g_type_entropy", 1.0)]),
 ("F16", "grid", "BTW sandpile, slow drive: reported perfectly regular, non-bursty periodic toppling",
  L64 + "Rules each step: one grain drops on a random cell; a cell reaching 4 grains topples one grain to each of its 4 "
  "neighbours (grains fall off the edge), repeated until no cell has 4 or more.",
  ["Same rules on a 128x128 lattice."],
  [("g_activity_cv", -2.0), ("g_spec_peak", 6.0)]),
 ("F17", "network", "random rewiring: reported all node degrees converging to exactly the mean (degree spread vanishes)",
  NET + "Rules each step: a random link is cut and one of its ends is reconnected to a node chosen uniformly at random.",
  ["Same rules, mean degree 10."],
  [("n_degree_cv", -3.0), ("n_edge_turnover", 4.0)]),
 ("F18", "lattice-phase", "weakly coupled lattice oscillators with frequency spread: reported full global synchrony",
  "World: a 64x64 lattice of phase oscillators. Natural frequencies drawn from a normal distribution (mean 0, sd 1). "
  "Rules each step: each phase advances by its natural frequency plus coupling 0.1 times the mean sine of the phase "
  "differences to its 8 neighbours.",
  ["Same rules, coupling 0.2."],
  [("p_r", 4.5), ("p_local_r", 4.0), ("p_r_std", -1.5)]),
 ("F19", "grid", "Schelling, very intolerant (threshold 0.9), 5% vacancies: reported fast freezing into two large domains",
  L64 + "Cells have 3 types (type 0 = empty, 5% of cells; random start). Rules each step: an occupied cell whose "
  "same-type share among occupied neighbours is below 0.9 moves to a random empty cell.",
  ["Same rules, threshold 0.85."],
  [("g_moran", 5.0), ("g_largest_domain", 4.0), ("g_activity_trend", -5.0)]),
 ("F20", "grid", "elementary CA rule 30: reported freezing into a uniform state",
  "World: a ring of 400 cells (2 neighbours each), 2 types (random start). Rules each step, all cells at once: a cell's "
  "new type is given by Wolfram's elementary rule 30 applied to itself and its two neighbours.",
  ["Same rule on a ring of 1000 cells."],
  [("g_activity_trend", -4.0), ("g_largest_domain", 6.0), ("g_spec_peak", 0.0)]),
]

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    cards, truth = [], {}
    for kind, rows in (("known", K), ("fabricated", F)):
        for cid, fam, desc, plain, variants, obs in rows:
            cards.append({"class_id": f"ctrl:{cid}", "stage": "control", "family": fam, "size": 5, "shortest_len": 0,
                          "program_plain": plain, "also_in_class_plain": variants,
                          "observed": [{"feature": f, "value": v} for f, v in obs]})
            truth[f"ctrl:{cid}"] = {"kind": kind, "intended": desc}
    json.dump(cards, open(os.path.join(here, "controls.json"), "w"), indent=1)
    json.dump(truth, open(os.path.join(here, "truth.json"), "w"), indent=1)
    print(len(cards), "cards written")
