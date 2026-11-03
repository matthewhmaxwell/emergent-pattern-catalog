# planted:N15

## Target paper (ground truth; the gate never saw this)
- Title: Emergence of Active Foams in a Minimal Cognitive Flocking Model
- arXiv: 2610.15835 (v1 2026-10-11); authors: Aryan Dwivedi, Diganta Bhaskar, Anil Kumar Dasanna
- Behaviour the authors report as new: Self-propelled particles that only steer toward neighbours inside a narrow finite-range vision cone form a system-spanning, continuously restructuring fractal-like network of particle streams (dimension ~1.5) that never coarsens.
- Authors' novelty statement: "We find that this regime exhibits a rich state diagram that combines characteristics of both limiting cases while giving rise to qualitatively new collective states."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: 10,000 point particles in a 100x100 periodic box (density 1), with no collisions. Rules: (1) each moves forward at constant speed 0.1 along its heading (small inertia); (2) it sees only particles within distance R (5 to 10) inside a vision cone of half-angle pi/8 around its heading; (3) it turns toward the positions of the visible particles, turning rate 3 times the average sine of the angle between its heading and the bearing to each (no copying of their headings); (4) its heading gets rotational noise of strength 0.08. Seeing is one-sided: being seen does not imply seeing back. Run 10 million steps of size 0.005.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Same rules with vision half-angle pi/12.
  - Same rules with turning strength 1 (the state extends to wider vision angles).

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - pattern_fractal_dimension = 1.5  (fractal dimension of particle pattern (1 = lines, 2 = uniform))
  - void_area_scaling_exponent = 3  (areas of empty voids scale as vision range to this power)
  - commonest_junction_degree = 3  (number of particle streams meeting at the most common junction)
  - max_common_junction_degree = 4  (largest junction degree seen; five or more essentially absent)

This behaviour appears across 3 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Minimal cognitive flocking / vision-cone position-based (non-reciprocal) attraction model (Barberis & Peruani 2016); worm/polar-file and interconnected nematic-band regimes; iABP family (Negi, Winkler, Gompper)
- **behaviour_documented**: True
- **undocumented_phenomenology**: No source reports fractal (D~1.5) branching stream networks, void-area ~ R^3 scaling, or degree-3-dominated junction statistics for vision-cone attraction models; the R^3 exponent implies a second length scale (turning radius or interparticle spacing) or unfinished coarsening.
- **notes**: Same model class as Barberis-Peruani. Files led by incidental leaders, and small nematic bands that interconnect during a transient, are documented; 3-way merging junctions are expected from follower convergence. The twist is the quantitative network statistics; the R^3 void scaling plausibly matters but needs finite-size and run-length checks to exclude a slow coarsening transient. Full texts of Barberis-Peruani and Brigatti-Peruani 2026 were not read.
- **citations**: ["Barberis & Peruani, PRL 117, 248001 (2016); arXiv:1912.07929", "Negi, Winkler, Gompper, Soft Matter 18, 6167 (2022)", "Saavedra & Peruani, PRE 110, 064602 (2024)", "Binary mixtures of iABPs with visual perception, New J. Phys. (2025) doi:10.1088/1367-2630/ae140c", "Directed swarm assembly due to mixed misaligned perception-dependent motility, Commun. Phys. (2025)"]

### N_blocked_r2 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Position-based vision-cone ('minimal cognitive flocking') active particles forming files/nematic streams; network morphology analogous to Jones-type Physarum agent transport networks
- **behaviour_documented**: True
- **notes**: The behaviour (vision-cone, non-reciprocal, position-based attraction producing files/streams that form a network) is documented. The twist is the quantified network shape: fractal dimension ~1.5, void area ~R^3, mostly degree-3 junctions with degree 5+ absent. Mostly degree-3 junctions is generic for coarsening 2D networks. R^3 void scaling is the only possibly meaningful item, but it rests on a factor-2 range of R in a small box and may reflect unfinished coarsening. The fractal dimension covers only about one decade of scale. I could not read the Barberis-Peruani full text to check whether their nematic regime looks like a network.
- **citations**: ["Barberis & Peruani, PRL 117, 248001 (2016)", "Negi, Winkler & Gompper, Soft Matter 18, 6167 (2022)", "Self-trapping of active particles with nonreciprocal interactions in disordered media, PRE 110, 064602 (2024)", "Jones, Characteristics of pattern formation and evolution in approximations of Physarum transport networks, Artificial Life 16(2) (2010)"]

### N_blocked_r3 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Filament/nematic-lane network state of the minimal cognitive flocking (Barberis-Peruani vision-cone, position-based attraction) model
- **behaviour_documented**: True
- **notes**: Interconnecting stream/nematic-band networks are documented for vision-cone position-attraction particles; honeycomb-like junction networks appear in iABP mixtures; similar networks appear in Physarum-type forward-sensing agents. The twist is quantitative only: D_f ~1.5, void area ~R^3, mostly 3-way junctions. Mostly 3-way junctions is generic for 2D networks. D_f comes from barely one decade of scale. The R^3 exponent comes from only a factor-of-2 range of R and is dimensionally suggestive of a second length (turning radius or persistence length); worth checking over a wider R range before treating it as meaningful. Barberis-Peruani full text not verified.
- **citations**: ["Barberis & Peruani, PRL 117, 248001 (2016)", "Collective behavior based on topological vision, Phys. Rev. Research 7, 023234 (2025)", "Negi, Winkler, Gompper et al., Binary mixtures of intelligent active Brownian particles with visual perception, New J. Phys. (2025), doi:10.1088/1367-2630/ae140c", "Negi, Winkler, Gompper, Soft Matter 18, 6167 (2022)", "Jones, Artificial Life 16, 127-153 (2010)"]
