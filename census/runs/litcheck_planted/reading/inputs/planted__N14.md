# planted:N14

**Reported behavior (what the authors report as new):** Hard ellipses with two polar sticky patches self-assemble into a kagome lattice that nucleates in one step directly from the gas (no chain intermediate), and into disordered porous networks at other conditions.

**Rule on the card:** World: a 2D surface exchanging particles with a reservoir (grand-canonical Monte Carlo; about 1,300 particles at close packing). Particles are hard elliptical disks with axis ratio 2 that cannot overlap; each carries two sticky patches, one at each tip of its long axis (patch half-angle 45 degrees, reach 0.17 particle diameters, bond energy -1), allowing up to two bonds per patch. Rules: Monte Carlo moves (single and cluster translations and rotations, insertions, deletions) at temperature T and chemical potential mu. Start from a dilute gas (density 0.004) quenched to T=0.145 with mu between -5 and 1 (full grid: T 0.055-0.16, mu -7 to 5). Track bonds and closed rings of bonded particles.

**Observed features on the card:**
  - six_ring_order = 0.8 (global order of six-particle rings (1 = perfect lattice of triangles))
  - fraction_in_rings = 1.0 (fraction of particles belonging to closed rings)
  - triangle_bonds_after_assembly = 340 (three-patch bonding triangles present after equilibration)
  - pair_bond_peak_low_temperature = 90 (peak count of two-patch pair bonds during nucleation, T=0.07)
  - pair_triangle_onset_lag = 0 (delay between rises of pair and triangle bond counts)

## Reply R10

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Kagome lattice self-assembly of two-patch (triblock-Janus / Kern-Frenkel) particles in the two-bonds-per-patch limit", "notes": "Kagome lattice of ellipse centres (three tips meet at each triangle junction, six-particle rings) has the same connectivity as the triblock-Janus kagome; the elliptical core changes geometry but not bonding topology. Twist: at deep quench (T=0.07) pair and triangle bonds rise at the same time (zero lag) and the gas goes straight to the lattice, whereas the literature mostly reports two-step nucleation via a disordered network. This is probably expected kinetics at deep supercooling and/or limited by sampling resolution, so it likely does not matter. No paper with tip-patched hard ellipses of aspect ratio 2 was found. fraction_in_rings=1 and six_ring_order=0.8 indicate a genuine ordered open lattice, not a gel or artifact."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium confidence. This agrees with the first pass.** The main behaviour is already well documented: two-patch particles where each patch can bond to at most two others assemble into a kagome lattice. The elliptical shape and the zero lag between pair and triangle bonds are the only possibly new parts. I don't think either is a real discovery.

**Why the main behaviour is known**
- **Experiment.** Chen, Bae and Granick made spheres that are sticky at the two poles and repel around the middle ("triblock Janus" particles). They self-assemble into an intricate kagome-patterned open crystalline lattice (Nature 469, 381, 2011). The authors tracked the ordering using triangular bonding (blue triangles), the same triangle motif your system counts. They counted particles with four nearest neighboring bonds as having local kagome order.
- **Simulation with the same rule structure.** Romano and Sciortino used a two-patch model with a short-range attraction mimicking hydrophobic interactions and a patch width that allows for a maximum of two contacts per patch. They showed that it is able to reproduce the observed crystallisation pathway in the Kagome structure. They also found the kagome lattice is stable at low temperature and low pressure. That fits your low-T, grand-canonical setup.
- **Simple models with other shapes.** Kagome structures also appear for soft ellipsoid-based particles (Fejer et al. 2015), which can form energetically stabilized Kagome structures. They also appear for vibrating patchy particles (Chapela et al. 2014). A grand-canonical study of two-patch particles extended to 3 and 4 patches also exists. It reports kagome nucleation barriers of a similar size and a two-step mechanism (Eslami, JCTC 2024).
- **Ellipse-core patchy models.** These exist as well, e.g. Kahl's group's two dimensional model for elliptic (i.e. shape anisotropic) patchy colloids. In that work the patches sit on the short axis, not the tips. I found no paper with patches at the tips of an ellipse with aspect ratio 2.

**The twist I couldn't find described**
At T=0.07, pair bonds and triangle bonds start rising at the same time, and the system goes straight from gas to lattice. The literature mostly reports two steps instead: a disordered network forms first, then reorganises. The connection of colloidal particles through their patches initially generates a disordered liquid network. Subsequently, orientational rearrangements in the liquid precursors lead to the formation of ordered nuclei.

I doubt this difference matters, for three reasons:
1. At T=0.07, about half the target temperature, the system is deeply quenched. Bonds form almost as fast as particles meet, so pairs and triangles would naturally appear together.
2. A pair-bond count and a triangle count overlap, and a lag of 0 may simply be below the sampling resolution.
3. The two-step picture is about the order parameter (when the network becomes ordered), not about when bond types first appear. The observed metrics don't check it either way.

The elliptical shape doesn't change how the particles connect. Three tips meet at each triangle junction, each particle links two junctions, and the particle centres form a kagome (trihexagonal) lattice.

**Genuine pattern, not an artifact.** Every particle sits in a closed ring, six-ring order is 0.8, and about 340 triangles survive equilibration. That points to a real, ordered open network, not a frozen gel or random bonding. One caveat: the order value of 0.8 leaves room for defects or grain boundaries. I also couldn't confirm the zero-lag result against the time traces of Chen et al.

## Reply R31

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Kagome lattice self-assembly from two-patch (triblock Janus-like) particles in the two-bonds-per-patch limit", "notes": "Kagome identification is inferred from patch geometry (two patches, two bonds per patch, triangle bonds, six-rings, 45 deg half-angle), not stated in the rule description. No source found for hard ellipses (AR 2) with tip patches under GCMC; only abstracts/snippets were read, so the quantitative observables were not compared against any paper. Direct gas-to-kagome nucleation with defects at low T and the zero pair/triangle lag are not explicitly confirmed in the sources found."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The behaviour is a kagome lattice self-assembling from two-patch particles in the two-bonds-per-patch limit. This is well documented for patchy disks and spheres. I found no source that uses hard ellipses with axis ratio 2 and tip patches in a grand-canonical simulation.

**1. Model family.** Patchy-colloid self-assembly: a Kern–Frenkel-type two-patch model on a 2D surface, simulated by Monte Carlo. This is statistical physics and soft matter, not cellular automata or opinion dynamics.

**What the numbers say.** This is my own reading of the observables, not something a source states.
- Each patch allows two bonds, so three particles can meet in a bonded triangle. The reported "three-patch bonding triangles" are therefore the triangles of a kagome lattice.
- Kagome has hexagonal rings of six particles sitting on a triangular lattice, which fits the six_ring_order of 0.8 and the fraction_in_rings of 1.0.
- The 45° patch half-angle fits kagome geometry. A neighbouring bond lies about 30° off the long axis, which the 45° patch allows.
- It would not fit a honeycomb layout of ellipses as edges. That layout needs about 60° off-axis bonds. I therefore ruled it out.

**2. Earliest clear sources (checked against search results).**
- **Experiment:** Chen, Bae & Granick, "Directed self-assembly of a colloidal kagome lattice," Nature 469, 381–384 (2011).
  - Colloidal spheres with a simple pattern of hydrophobic domains self-assemble into a kagome lattice.
  - Particles first cluster into kinetically favoured triangles, strings or a mix of the two.
  - This is the same triangle-first pathway as the "triangle bonds rise with no lag after pair bonds" feature, though the paper reports it from experiment rather than simulation.
- **Closest simulation:** Romano & Sciortino, "Two dimensional assembly of triblock Janus particles into crystal phases in the two bond per patch limit," Soft Matter 7, 5799 (2011).
  - The model uses a short-range attraction and a patch width that allows at most two contacts per patch.
  - It reproduces the kagome crystallisation pathway, and kagome is stable at low temperature and low pressure.
  - Their own phrasing of "two bonds per patch" matches your rule exactly.
- **Review, entropy:** Mao, Chen & Granick, "Entropy favours open colloidal lattices," Nature Materials 12, 217–222 (2013).
  - Cates' News & Views accompanies it (Nature Mater. 12, 179–180, 2013).
  - Together they explain why the open kagome lattice is stable. Kagome and hexagonal lattices share the same rotational entropy, but close-packed lattices have lower vibrational entropy.
- **Review, general:** the ScienceDirect topic overview and the Comptes Rendus Chimie account, "Patchy colloidal particles for programmed self-assembly." I could not confirm the second one's author list from the search output.
- **Related:** Chapela et al., Soft Matter 10, 9167 (2014), on kagome from vibrating 2D patchy particles. Particles with extra repulsive sites form kagome but connect as entangled webs at low temperature. This resembles your "more defects at low T" variant.
- **Nucleation mechanism:** J. Chem. Theory Comput. 20, 6858 (2024), on 3- and 4-patch particles in 2D. It uses grand canonical simulations. It reports a two-step nucleation: a dense amorphous network first, then reorientation into crystalline nuclei. I did not find the authors in the results. Its stated nucleation barriers are close to the kagome values.

**3. What differs from the literature, and whether it matters.**
- The particle shape is an ellipse with axis ratio 2 and tip patches, not a sphere or disk with polar caps. Chen et al.'s particles had elongated caps, so patch shape and anisotropy were already part of the experiment. Shape elongation changes steric packing and may shift the stability window. It is unlikely to change the lattice identity, because the patch geometry sets that.
- The setting is grand-canonical Monte Carlo with a quench from a dilute gas. Romano & Sciortino used constant-pressure phase-diagram and free-energy calculations, and Chen et al. used experiments. The nucleation-from-gas protocol is a different route to the same structure.
- The low-T result (the same open lattice nucleating directly with more defects) matches the usual kinetic-trapping picture. I did not find a source that specifically describes direct gas-to-kagome nucleation with the pair-then-triangle timing, so I treat it as likely but unverified.

**

## Reply R52

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Kagome lattice self-assembly of triblock-Janus-like two-patch particles in the two-bond-per-patch limit (here with elliptical hard cores)", "notes": "Identifying the lattice as kagome is inferred from the rule structure (two polar patches, two bonds per patch, triangles plus six-rings), not stated in the prompt. No paper using ellipses with tip patches was found. Zero pair-to-triangle lag and direct gas-to-lattice nucleation at low T were not checked against the literature. The literature mainly reports two-step nucleation via an amorphous network, so this is a possible small difference. Several authors and years were not verifiable from search snippets."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** Patchy particles that form a kagome lattice are well documented. I found no source that uses hard ellipses with tip patches, so that exact shape is the variant.

**1. Model family.** This is a patchy-colloid self-assembly model in 2D, in the Kern–Frenkel style: hard-core particles with short-range square-well patches. The ensemble is grand-canonical Monte Carlo with cluster moves. The patch width allows two bonds per patch, which is the "two-bond-per-patch limit" of two-patch particles.

**2. What the observed behaviour is.** The reported features are a six-ring order of 0.8, every particle in a closed ring, and many bonded triangles. Together these describe a kagome lattice, which is a network of corner-sharing triangles and hexagons. A closed six-ring is the kagome hexagon. A triangle is three mutually bonded particles, which arises when each polar patch binds two neighbours. Kagome is the standard result for two-patch particles with this valence.

**3. Sources I checked.**
- **Experiment, earliest clear source.** Chen, Bae and Granick, "Directed self-assembly of a colloidal kagome lattice," Nature 469, 381 (2011). Spheres that attract at two polar regions and repel at the middle band assemble into a kagome lattice. Their design uses two patches at opposite poles, each wide enough in the plane to allow coordination with two nearest neighbours. The elongated caps in their experiment also allow two nearest neighbours only when the long patch axes of neighbouring particles are parallel. I read this in the Granick-group PDF and the OSTI copy.
- **Simulation of the equivalent model.** Romano and Sciortino, "Two dimensional assembly of triblock Janus particles into crystal phases in the two bond per patch limit," Soft Matter 7, 5799 (2011). The patch width allows a maximum of two contacts per patch. Free-energy calculations show kagome is stable at low temperature and low pressure. At higher pressure it transforms into a hexagonal lattice with alternating attractive and repulsive bands.
- **Reviews and perspectives.**
  - Romano and Sciortino, "Patchy from the bottom up," Nature Materials 10, 171 (2011).
  - Annual-review chapter "Toward Design Rules of Directional Janus Colloidal Assembly" (Northwestern, csml.northwestern.edu reprint). I did not verify its authors or year.
- **Nucleation pathway.** Eslami et al., Small 2024 (Wiley, "Self-Assembly Pathways of Triblock Janus Particles into 3D Open Lattices"), states that nucleation of a kagome lattice from solution of triblock Janus particles obeys a two-step mechanism. That means a dense amorphous network forms first, followed by crystalline nuclei. It also quotes a barrier of 5–10 k_BT.
- **Closest methodological match.** "Hierarchical self-assembly of patchy colloidal platelets," Soft Matter 16, 2774 (2020). It uses 2D rhombi with two patches, simulated with GCMC, single-particle moves and cluster moves, so the method matches. Its search snippet did not show authors, so I am not citing them. Related work is "Partially Bonded Crystals," ACS Nano (2025), which also uses rhombi and kagome-type lattices.
- **Lower-temperature, more-defective result.** The statement that lower T gives the same open lattice with more defects is the standard kinetic-trapping picture. A related result is the zigzag trapping reported for one-patch particles at strong attraction.

**4. What differs, and whether it matters.**
- The particle shape is an ellipse with axis ratio 2 and tip patches. I found no paper that uses exactly this, and one more search is the only check I could not fully make. Elongation matches the elongated caps in the Granick experiment, and rhombi are the closest 2D anisotropic-shape analogue. I don't expect the shape to change the outcome. The tip patches with a 45° half-angle and two bonds per patch fix the local geometry that gives kagome.
- The "pair_triangle_onset_lag = 0" and "direct nucleation from the gas" at low T may differ from the amorphous-network-first, two-step picture. That picture was reported near the kagome critical region. At low T and low mu, near-simultaneous growth of pairs and triangles is plausible, but I found no source that tests it. This is a small possible novelty, not a gap in the main finding.
- The quantitative values (0.8, 340, 90) are specific to this simulation and are not in the literature. I did not look for them.

**5. Confidence: medium.** The identification of kagome from two-bond-per-patch polar patches is solid, and the experiment, simulation and review sources all agree. Confidence is only medium because I found no paper that uses ellipses. I also did not directly check whether a patchy-ellipse paper reports the direct-nucleation timing.

**
