# planted:N14

## Target paper (ground truth; the gate never saw this)
- Title: Networks of triblock ellipses: from kagome order to disordered porous networks
- arXiv: 2610.08522 (v1 2026-10-06); authors: Susanne Wagner, Gerhard Kahl, Carina Karner
- Behaviour the authors report as new: Hard ellipses with two polar sticky patches self-assemble into a kagome lattice that nucleates in one step directly from the gas (no chain intermediate), and into disordered porous networks at other conditions.
- Authors' novelty statement: "This is – to the best of our knowledge – the first observation of a kagome lattice assembled from patchy ellipses."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: a 2D surface exchanging particles with a reservoir (grand-canonical Monte Carlo; about 1,300 particles at close packing). Particles are hard elliptical disks with axis ratio 2 that cannot overlap; each carries two sticky patches, one at each tip of its long axis (patch half-angle 45 degrees, reach 0.17 particle diameters, bond energy -1), allowing up to two bonds per patch. Rules: Monte Carlo moves (single and cluster translations and rotations, insertions, deletions) at temperature T and chemical potential mu. Start from a dilute gas (density 0.004) quenched to T=0.145 with mu between -5 and 1 (full grid: T 0.055-0.16, mu -7 to 5). Track bonds and closed rings of bonded particles.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Lower temperatures 0.07-0.115 at mu=-5: the same open lattice nucleates directly from the gas, with more defects.

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - six_ring_order = 0.8  (global order of six-particle rings (1 = perfect lattice of triangles))
  - fraction_in_rings = 1.0  (fraction of particles belonging to closed rings)
  - triangle_bonds_after_assembly = 340  (three-patch bonding triangles present after equilibration)
  - pair_bond_peak_low_temperature = 90  (peak count of two-patch pair bonds during nucleation, T=0.07)
  - pair_triangle_onset_lag = 0  (delay between rises of pair and triangle bond counts)

This behaviour appears across 1 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Kagome lattice self-assembly of two-patch (triblock-Janus / Kern-Frenkel) particles in the two-bonds-per-patch limit
- **behaviour_documented**: True
- **notes**: Kagome lattice of ellipse centres (three tips meet at each triangle junction, six-particle rings) has the same connectivity as the triblock-Janus kagome; the elliptical core changes geometry but not bonding topology. Twist: at deep quench (T=0.07) pair and triangle bonds rise at the same time (zero lag) and the gas goes straight to the lattice, whereas the literature mostly reports two-step nucleation via a disordered network. This is probably expected kinetics at deep supercooling and/or limited by sampling resolution, so it likely does not matter. No paper with tip-patched hard ellipses of aspect ratio 2 was found. fraction_in_rings=1 and six_ring_order=0.8 indicate a genuine ordered open lattice, not a gel or artifact.
- **citations**: ["Chen, Bae & Granick, Nature 469, 381 (2011) - Directed self-assembly of a colloidal kagome lattice", "Romano & Sciortino, Soft Matter 7, 5799 (2011) - 2D assembly of triblock Janus particles in the two bond per patch limit", "Fejer et al., Soft Matter 11, 6663 (2015) - Design of a Kagome lattice from soft anisotropic particles", "Eslami et al., J. Chem. Theory Comput. 14 (2018) - Self-Assembly Mechanisms of Triblock Janus Particles (two-step kagome nucleation)", "Eslami, J. Chem. Theory Comput. 20, 6858 (2024) - GCMC two/three/four-patch 2D lattices", "Chapela et al., Soft Matter 10, 9167 (2014) - kagome from vibrating patchy particles"]

### N_blocked_r2 (strict pass): verdict KNOWN
- **standard_name**: Kagome lattice self-assembly of two-patch (triblock-Janus-type) particles in the two-bonds-per-patch limit
- **behaviour_documented**: True
- **notes**: Mild overturn: the first pass's 'variant' label rested on rule differences (ellipse shape, GCMC), not on behaviour. With ellipse tip patches, the particles sit on honeycomb edges and their centres form kagome, so the network is the same. The zero-lag simultaneous rise of pair and triangle bonds is not explicitly documented, but it is a minor kinetic detail that probably depends on sampling. The structure was inferred from ring and triangle counts, not seen; six-ring order 0.8 suggests defects or polycrystallinity. Check that the high-mu runs have not turned into the denser hexagonal phase that Romano and Sciortino report at high pressure.
- **citations**: ["Romano & Sciortino, Soft Matter 7, 5799 (2011), doi:10.1039/c0sm01494j", "Chen, Bae & Granick, Nature 469, 381 (2011), doi:10.1038/nature09713", "Design of a Kagome lattice from soft anisotropic particles, Soft Matter 11, 6663 (2015)", "Partially Bonded Crystals, ACS Nano 19, 5146 (2025)", "Elliptic patchy particles (co-vertex patches), Mol. Phys. (2019), doi:10.1080/00268976.2019.1592254"]

### N_blocked_r3 (strict pass): verdict KNOWN
- **standard_name**: Kagome lattice self-assembly of two-patch (triblock-Janus-like) particles in the two-bonds-per-patch limit, via two-step nucleation through a disordered network
- **behaviour_documented**: True
- **notes**: Mild overturn: the first pass said VARIANT; the behaviour itself is documented. Ellipses with three-tip junctions form a honeycomb network of particles whose centres sit on a kagome lattice. The elliptical hard core in place of a repulsive middle band is a rule difference, not new behaviour. Pairs and triangles appearing together at deep quench matches the Granick experiment (triangles and strings form first, then defective networks). Possible open point, not tested here: elongated cores may suppress the high-pressure hexagonal phase. Check snapshots for open hexagonal holes to rule out six-ring order coming from a dense phase at high chemical potential. Sources checked via abstracts and summaries, not full texts.
- **citations**: ["Romano & Sciortino, Soft Matter 7, 5799 (2011), doi:10.1039/C0SM01494J", "Chen, Bae & Granick, Nature 469, 381 (2011), Directed self-assembly of a colloidal Kagome lattice (OSTI 1875134)", "Eslami et al., J. Chem. Theory Comput. (doi:10.1021/acs.jctc.1c01116), Self-Assembly of Model Triblock Janus Colloidal Particles in Two Dimensions", "Design of a Kagome lattice from soft anisotropic particles, Soft Matter 11, 6663 (2015)", "Self-assembly and time-dependent control of active and passive triblock Janus colloids (PubMed 40709776)"]
