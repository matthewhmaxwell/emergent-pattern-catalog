# planted:N05

## Target paper (ground truth; the gate never saw this)
- Title: Phase ordering kinetics in Light-Heavy-Vacancy model: unusual coarsening dynamics
- arXiv: 2609.14541 (v1 2026-09-13); authors: Chandradip Khamrai, Sakuntala Chatterjee
- Behaviour the authors report as new: In a 1D two-species driven lattice gas on a fluctuating landscape, early ordered domains disintegrate and new ones form, giving two coarsening exponents (0.5 then 0.34) and a non-monotonic landscape length scale.
- Authors' novelty statement: "Unlike conventional phase-ordering systems, where coarsening proceeds through the formation and merger of ordered domains, here we find that domains formed at early times become unsustainable at later times."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: a ring of 16,416 sites; between neighbouring sites is a bond tilted up or down (half of each, random), defining a height profile. Each site holds an A particle, a B particle or nothing (one third each, random). Rules, random sequential updates: (1) a particle swaps with a neighbouring site of different content; a move taking A downhill or B uphill succeeds with probability 0.9, the opposite move with 0.1; (2) a local peak (up-bond then down-bond) and a dip swap: a peak holding A becomes a dip with probability 0.5+b and back with 0.5-b, b=-0.4 (so A pushes ground up); a dip holding B becomes a peak with probability 0.5+b' and back with 0.5-b', b'=0.3; empty sites flip either way with 0.5. Track domain size over time.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Same rules with b'=0.18 (late growth exponent about 0.41).

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - early_growth_exponent = 0.5  (power-law exponent of domain size growth before t~1000)
  - late_growth_exponent = 0.34  (power-law exponent of domain size growth at late times)
  - crossover_time = 1000  (time at which early domains dissolve and growth slows)
  - late_macroscopic_hills = 1  (macroscopic hills in the height profile in the late regime)

This behaviour appears across 2 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Coupled particle-landscape (light-heavy) model of Chakraborty-Chatterjee-Barma: strong phase separation / compact-cluster ordered phases with algebraic coarsening (not FDPO, which only occurs on the separatrix)
- **behaviour_documented**: True
- **notes**: Agrees with first pass on the model family (KNOWN-VARIANT), but corrects the label: generic couplings (b=-0.4, b'=0.3) put the system in an ordered phase-separated phase, not FDPO, which lives only on the order-disorder separatrix. Unverified twist: the explicit early t^0.5 -> late t^~0.34 crossover near t~1000 with early-domain dissolution; full texts were not accessible. The late exponent shifts with b' (0.34 vs 0.41), suggesting effective/crossover exponents rather than a new class. Mapping of these parameters to the compact-cluster phase is inferred, not verified.
- **citations**: ["Chakraborty, Pal, Chatterjee, Barma, PRE 93, 050102(R) (2016)", "Chakraborty, Chatterjee, Barma, PRE 96, 022127 (2017)", "Chakraborty, Chatterjee, Barma, PRE 96, 022128 (2017)", "Mahapatra, Ramola, Barma, Phys. Rev. Research 2, 043279 (2020)"]

### N_blocked_r2 (first pass): verdict KNOWN
- **standard_name**: Coupled driven systems / light-heavy model of particles on a fluctuating surface (Chakraborty-Chatterjee-Barma); fluctuation-dominated phase ordering
- **model_family**: coupled driven lattice gas: two particle species + vacancies advected by a fluctuating (tilt-field) landscape that they modify (nonequilibrium statistical mechanics / phase separation)
- **notes**: Model and qualitative behaviour (phase separation, macroscopic landscape order, algebraic coarsening) clearly documented. Full texts were paywalled, so the specific two-stage exponents (0.5 -> ~0.34, crossover ~1000, 0.41 at b'=0.18) and the exact parameter values (b=-0.4, b'=0.3) were not checked against the papers. The late exponent may be a finite-time effective value. If the specific two-stage exponent crossover matters, check the 2016 PRE Rapid Comm. and 2017 Dynamic-properties paper full text.
- **citations**: ["Chakraborty, Pal, Chatterjee, Barma 2016, Phys. Rev. E 93, 050102(R) - macroscopic clusters, phase-separated species, ordered+disordered landscape, algebraic coarsening times", "Chakraborty, Chatterjee, Barma 2017, Phys. Rev. E 96, 022127 - static ordered phases: species phase separation, landscape with macroscopic regions of distinct slope", "Chakraborty, Chatterjee, Barma 2017, Phys. Rev. E 96, 022128 - dynamics of the ordered phases", "Phys. Rev. Research 2, 043279 (2020) - light-heavy model, coarsening length scale and early-time dynamics", "Das & Barma 2000, Phys. Rev. Lett. 85, 1602 - sliding particles on fluctuating surface, coarsening length ~ t^{1/z}", "Das, Barma, Majumdar 2001, Phys. Rev. E 64, 046126 - fluctuation-dominated phase ordering on stochastic surfaces"]

### N_blocked_r3 (first pass): verdict KNOWN
- **standard_name**: Coupled two-species sliding particles on a fluctuating landscape (light-heavy model): complete phase separation with ordered/disordered landscape coexistence and coarsening
- **model_family**: Nonequilibrium statistical mechanics: driven two-species lattice gas coupled to a fluctuating landscape (tilt fields); Lahiri-Ramaswamy / Barma-Chatterjee coupled-landscape (light-heavy) model
- **notes**: Model identity and the ordered-landscape phase are well supported. The specific exponents (0.5 early, 0.34 late, 0.41 for b'=0.18) and the crossover at t~1000 were not verified because the full texts were paywalled. Rate mapping to the published parameters was not checked. If the two-stage exponent values matter, read PRE 96 022128 and the PRE 93 050102 supplement.
- **citations**: ["Chakraborty, Pal, Chatterjee, Barma 2016, Phys. Rev. E 93, 050102(R) - macroscopic clustering, species phase-separated, landscape ordered regions coexisting with disordered region, algebraic coarsening timescales", "Chakraborty, Chatterjee, Barma 2017, Phys. Rev. E 96, 022127 - static properties of the ordered phases of the same model", "Chakraborty, Chatterjee, Barma 2017, Phys. Rev. E 96, 022128 - dynamic properties; clusters move on exponential timescale, landscape on power-law timescale", "Mahapatra, Ramola, Barma 2020, Phys. Rev. Research 2, 043279 - light-heavy model early-time and coarsening dynamics", "Das, Barma, Majumdar 2001, Phys. Rev. E 64, 046126 - single-species sliding particles on fluctuating surfaces, coarsening (distant relative)", "Ramaswamy, Barma, Das, Basu 2002, Phase Transitions 75, 363 - two-species lattice model with linear instability (Lahiri-Ramaswamy lineage)"]
