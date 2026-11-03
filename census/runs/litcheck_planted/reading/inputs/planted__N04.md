# planted:N04

**Reported behavior (what the authors report as new):** A 2D solid with transverse (odd) pair forces melts by defect-pair proliferation at any temperature above a threshold activity, and recrystallises at higher activity even at temperatures where the passive crystal is liquid.

**Rule on the card:** World: 10,000 slightly elongated particles (ellipse half-axes 0.55 and 0.65) in a periodic 2D box at density 0.87 (lattice spacing about 1.15). Rules, overdamped Brownian dynamics at temperature T: (1) each pair interacts through an orientation-dependent Lennard-Jones potential of strength e_p (repulsive core, weak attraction, depending on both particles' orientations, which can rotate); (2) each pair also feels a non-conservative force with Lennard-Jones magnitude and strength e_a, directed perpendicular to the line joining them, always in the same rotational sense; (3) thermal noise acts on positions and orientations, orientations relaxing faster. Sweep the ratio e_a/e_p from 0 to 33 and T from 0 to 6; measure six-fold bond order and defect counts.

**Observed features on the card:**
  - activity_disordered_at_zero_temperature = 2.0 (force ratio e_a/e_p where solid is disordered even at T=0)
  - activity_of_restored_order = 33.0 (force ratio e_a/e_p where crystalline order returns)
  - temperature_of_restored_order = 6.0 (temperature at which that high-activity solid stays ordered)
  - passive_partial_melting_temperature = 4.25 (temperature where the passive solid is only six-fold oriented)

## Reply R42

**JSON fields:** {"verdict": "KNOWN-VARIANT", "behaviour_documented": true, "undocumented_phenomenology": "", "standard_name": "Odd/transverse-force (chiral Lennard-Jones) melting of a 2D crystal into rotating hexatic patches, on top of KTHNY two-step passive melting", "notes": "Melting of the solid by transverse forces (disordered at e_a/e_p~2 even at T=0) and the passive hexatic at T~4.25 are documented. Twist: crystalline order returning at e_a/e_p~33 and T=6, above the passive hexatic temperature, was not found in the literature. It is suspicious, probably caused by huge per-step forces in overdamped integration (timestep), a rigidly rotating or jammed state, or an order-parameter artefact. Treat it as unverified until checked with a smaller timestep, different system sizes, mean-squared displacement and overlap checks. The particle anisotropy (ellipses with half-axes 0.55/0.65) seems to change details only, not the phenomenology."}

**Reply text:**

**Verdict: KNOWN-VARIANT.** The main behaviour is documented: a transverse (odd) pair force melts a 2D crystal. The one new-looking feature is that order comes back at e_a/e_p ≈ 33 and T = 6. I couldn't find that described, and I think it is more likely a numerical artefact than a real effect. I partly agree with the first pass, but I put less trust in the order-restoration result than it did.

**What is documented**
- **The closest match is the "chiral Lennard-Jones fluid".** This is Brownian particles with Lennard-Jones forces plus a transverse pair force (Caporusso, Gonnella & Levis, PRL 132, 168201, 2024). The abstract says "Chirality can also break the solid phase, giving rise to a dense fluid made of rotating hexatic patches." That matches the observation that the solid is disordered even at T=0 once e_a/e_p ≈ 2. I only read the abstract.
- **Transverse forces are known to destabilise crystals.** Transverse forces organize colloidal spinners into odd elastic crystals crisscrossed by motile dislocations, and both odd and Hall stresses can destabilize odd elastic crystals (Bililign et al., Nature Physics). Braverman et al. (PNAS 2025) use pairwise transverse forces with a shared conservative repulsion. They find odd crystals break up through interacting moving dislocations, and note that transverse interactions cannot compete with elastic forces to fission isolated topological defects. I did not check either paper's author list beyond what the search results showed.
- **Order that drops and then returns as chiral activity grows has been seen in other systems.** Granular spinner experiments show a chiral activity-mediated reentrant melting transition at a fixed areal spinner density. Simulations of chiral active Brownian particles find that the fluid-crystal phase diagram shows a re-entrant melting transition as a function of the radius of the circles.
- **The passive two-step melting is textbook.** The passive solid at T ≈ 4.25 that keeps only six-fold orientation looks like the KTHNY hexatic phase. In the 2D Lennard-Jones system, long-range "orientational" order persists up to a higher temperature than the loss of shear rigidity. The slightly elongated, orientation-dependent particles change details only, not the phenomenology.
- **The odd-force model is in use for Lennard-Jones systems.** It has been applied to a two-dimensional Kob-Andersen mixture, where Brownian particles interact via the Lennard-Jones potential and nonconservative odd forces.

**The twist and why I doubt it**
The undocumented part is crystalline order returning at a transverse/conservative force ratio of 33 and T = 6. That is well above the passive hexatic temperature of 4.25. I found no report of large odd forces re-crystallising a thermal 2D solid like this. The reentrance results above are about spinner fractions or orbit radii, not about pushing the transverse force very high.

It looks like an artefact for these reasons:
1. A transverse force 33 times the conservative one means very large forces per step in overdamped integration. If the timestep wasn't reduced, the dynamics are probably under-resolved. Large odd forces can also make crystals linearly unstable.
2. A solid that is ordered at T=6 when its passive version is partly melted at 4.25 is physically odd, since noise should disorder more, not less.
3. Plausible causes include a timestep or integrator problem, particles being pushed into a rigidly rotating or jammed configuration, or the six-fold measure picking up a frozen state.

To check it: shrink the timestep, test whether results depend on system size, look at the mean-squared displacement and particle overlaps, and check whether the "ordered" state is just a rigidly rotating lattice.

Overall, the melting and hexatic behaviour is a genuine, known emergent pattern. Only the order-restoration corner looks suspect.

**

## Reply R47

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "chiral Lennard-Jones solid: transverse-force (odd) melting of a 2D crystal, with chirality-induced reentrant ordering", "notes": "Only abstracts and snippets were read. The activity-disordered solid is covered by Caporusso et al. 2024 for disks. No source was found for the ellipse variant, or for a high-activity solid (e_a/e_p ~33) that stays ordered up to T=6, above the passive partial-melting temperature of 4.25. The re-entrant restoration is likely related to known chiral re-entrance but I could not confirm it for this exact model."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium-low confidence.** The closest published system is the chiral Lennard-Jones fluid of Brownian disks. Two things are not confirmed: the ellipse-shaped particles, and the return of order at very large e_a/e_p. I read abstracts and snippets only, not full texts.

**1. Model family.** This is chiral active matter with "odd" interactions, a branch of active matter and non-equilibrium statistical physics. Each pair feels a Lennard-Jones force plus a non-conservative force perpendicular to the line joining the pair, with a fixed sense of rotation. It models spinning colloids. The particles are slightly anisotropic, with an orientation-dependent Gay-Berne-like potential. Melting is measured by hexatic (six-fold) order ψ6 and defect counts, as in KTHNY-type 2D melting.

**2. Closest literature.**
- **Caporusso, Gonnella and Levis (2024).** *Phase Coexistence and Edge Currents in the Chiral Lennard-Jones Fluid*, Phys. Rev. Lett. 132, 168201. It uses Brownian disks with a Lennard-Jones potential plus a nonconservative transverse force. Chirality can also break the solid phase, giving rise to a dense fluid made of rotating hexatic patches. This is the earliest clear source I found for the effect your simulation shows at the start of the sweep: a crystal disordered by transverse forces. Its thesis summary says chirality modifies both liquid-gas coexistence and fluid-solid transitions.
- **Review.** *Modeling chiral active particles: from circular motion to odd interactions*, J. Stat. Mech. (IOP, 2026). It describes effective transverse forces, acting perpendicular to the line connecting the centers of two interacting particles, referred to as odd interactions.
- **Related: Melting of nonreciprocal solids** (PNAS, 2025; the DOI is 10.1073/pnas.2412993122). It treats odd/transverse forces and contrasts them with nonreciprocal melting. It says transverse interactions cannot compete with elastic forces to fission isolated topological defects, so odd crystals fragment by a different route.
- **Related: Reentrant melting of scarred odd crystals by self-shear** (Nature Communications, 2026; PubMed 41545359). It reports a chiral activity-mediated reentrant melting transition in spinner crystals. The tuning parameter is the spinner handedness ratio, and the mechanism is edge shear in confined systems, so it is not the same as yours.
- **Related: Crystallization of chiral active Brownian particles at low densities** (Jeong et al., Phys. Rev. Research, 2025). It finds a re-entrant melting transition as a function of the radius of the circles. These are circle-swimmers, not transverse pair forces.

**3. What is and isn't covered.**
- **Covered:** a passive solid that melts in two steps, a crystal broken by transverse forces at e_a/e_p ≈ 2 even at T = 0, and order returning at large chirality, which resembles known re-entrance in chiral systems.
- **Not found, so treat as unconfirmed:**
  - A source with ellipse-shaped or Gay-Berne particles combined with a transverse force.
  - A source showing the high-activity solid (e_a/e_p ≈ 33) staying ordered up to T = 6, which is above the passive partial-melting temperature of 4.25.
- **Likely mechanism:** very strong transverse forces probably make the lattice rotate or flow coherently. That would give re-entrant order, as in the reentrant-melting and crystallization papers above. I haven't verified this for pair-force systems.
- **Does the shape matter?** Elongation adds orientational degrees of freedom, so it probably shifts the thresholds. It may change the quantitative picture, though not the qualitative one.

**
