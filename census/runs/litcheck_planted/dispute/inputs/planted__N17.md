# planted:N17

## Target paper (ground truth; the gate never saw this)
- Title: Coexistence phase of the two-dimensional Vicsek model
- arXiv: 2610.26561 (v1 2026-10-19); authors: Yu Duan, Hugues Chaté, Benoît Mahault, Xia-qing Shi
- Behaviour the authors report as new: In the coexistence phase of the 2D Vicsek model, side-by-side lanes of band trains travelling in opposite directions (and fat/inverse bands) form at particle level; long lanes buckle and the ordered band states are metastable.
- Authors' novelty statement: "Here we report on the diversity of band patterns observed in its coexistence phase –some of them never observed before–, but we also argue that most disappear in the infinite-size limit, leaving nevertheless a different picture than anticipated previously."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: self-propelled point particles in a periodic 2D box (square of side 1200, or rectangular), density 1 per unit area. Rules, discrete time steps: (1) each particle sets its heading to the average heading of all particles within distance 1 (itself included), plus a uniformly random angle in [-eta*pi, eta*pi], with noise eta between 0.28 and 0.35; (2) it moves 0.5 along its new heading. Start from random positions and headings, or from alternating parallel strips of opposite headings. Run up to 10^8 steps and watch the dense travelling stripes that form.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Continuous-time version (sine alignment plus Gaussian angular noise): same phenomenology.
  - eta=0.29, density 1, long rectangular box: very long lanes buckle and merge.

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - parallel_lane_count = 7  (side-by-side lanes of stripe trains travelling in opposite directions)
  - random_start_success_share = 0.042  (share of random starts ending in opposite-travelling lanes (1 of 24))
  - noise_window_low = 0.28  (lower noise of window where opposite-travelling lanes form)
  - noise_window_high = 0.35  (upper noise of window where opposite-travelling lanes form)
  - lane_persistence_steps = 100000000.0  (time steps lanes persisted without disappearing)

This behaviour appears across 1 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Vicsek model band (liquid-gas coexistence / microphase-separated) regime; counter-propagating lane ('sandwich') state of metastable flocks
- **behaviour_documented**: True
- **notes**: Bands, cross-sea and the metastability of flocks to counter-propagating droplets are documented. Persistent opposite-polarity lanes are documented in discrete-symmetry (active Potts 'sandwich state') and two-species/bidirectional variants, but not found for the single-species continuous-angle standard VM. The twist is a multi-lane (7) antiparallel band-train state lasting 1e8 steps that is reached rarely (1 of 24) from random starts. Plausibly matters as a long-lived metastable state of the coexistence regime. Must rule out pinning by periodic box axes: test tilted or incommensurate boxes and system-size scaling. Full texts and supplementary videos of the large-system VM papers were not read.
- **citations**: ["Kürsten & Ihle, PRL 125, 188003 (2020) - cross sea phase in standard Vicsek model", "Chatterjee et al., PRE 112, 064115 (2025) - sandwich state: persistent opposite-polarity lane in active Potts flocks", "Karmakar, Chatterjee, Mangeat, Rieger, Noh, New J. Phys. 26, 043023 (2024) - discretized Vicsek model, counter-propagating bands, cross-sea to band transition near eta~0.35", "Chatterjee et al., PRE 107, 024607 (2023) - two-species Vicsek model antiparallel flocking (APF) bands", "Benvegnen et al., PRL 131, 218301 (2023) - flocks metastable to counter-propagating droplets", "Lane formation in a bidirectional Vicsek-type model (Zenodo 23263766, 2026) - lanes absent in standard VM"]

### N_blocked_r2 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Vicsek-model travelling bands (microphase separation) arranged as counter-propagating lanes; single-species analogue of the antiparallel-flocking (APF) band state
- **behaviour_documented**: True
- **notes**: Bands are textbook. Counter-propagating band lanes are documented only in models with two species, two preferred directions or anti-alignment. The active Potts and bidirectional-model papers say lanes are absent in the standard Vicsek model. The twist: lanes of transverse-band trains travelling in opposite directions, persisting for 1e8 steps in the single-species metric Vicsek model. Most likely a geometry-protected metastable state: 1 of 24 random starts, strip starts build it in, and periodic boundaries make lane reversal rare. Checks needed: box size and aspect ratio, starts not aligned with the box. Full texts of Chaté 2008, Solon 2015 and Caussin 2014 not read.
- **citations**: ["Kürsten & Ihle, PRL 125, 188003 (2020) - band and cross-sea phases of the standard Vicsek model", "Chatterjee, Mangeat et al., PRE 107, 024607 (2023) - two-species Vicsek model, APF state with counter-moving bands", "Mangeat, Chatterjee, Paul, Rieger, PRE 102, 042601 (2020) - band-to-lane transition, stated absent in Vicsek model", "PRE 89, 012918 (2014) - constant-speed noisy bands unstable in head-on collisions", "Discretized Vicsek model, New J. Phys. (2024) 10.1088/1367-2630/ad3ea7 - counter-propagating dense bands", "Bidirectional Vicsek-type model lane phase, Zenodo 23263766 (2026)", "arXiv:2412.10577 - lanes from anti-alignment"]

### N_blocked_r3 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Cluster-lane pattern: parallel lanes of polar band/cluster trains moving in opposite directions (overall nematic order) in the band-forming regime of polar-alignment (Vicsek-class) active matter
- **behaviour_documented**: True
- **notes**: Same verdict class as the first pass, but it corrects the first pass's claim that the phenomenology was not found: opposite-travelling cluster-lane trains are documented as a stable solution in the kinetic theory of Vicsek-like polar alignment (Thüroff et al. 2014). Twist: seen in particle simulations of the standard metric Vicsek model as a rare (1/24 random starts), long-lived (1e8 steps) state, mostly reached from strip starts; no particle-simulation report found. Possible box-size locking of the 7 lanes and the parameter regime should be checked. Full texts of Chaté 2008, Solon 2015 and Kürsten-Ihle not read.
- **citations**: ["Thüroff, Weber & Frey, Phys. Rev. X 4, 041030 (2014) - stable cluster-lane limit cycle in the BDG Boltzmann equation", "Kürsten & Ihle, Phys. Rev. Lett. 125, 188003 (2020) - cross-sea phase in the standard Vicsek model (related non-band pattern)", "Mangeat, Chatterjee, Paul & Rieger, Phys. Rev. E 102, 042601 (2020) - band-to-lane transition in the active Potts model, stated to be absent in Vicsek", "Stability of discrete-symmetry flocks: sandwich state, Phys. Rev. E (2025), doi 10.1103/1r19-ryx9"]
