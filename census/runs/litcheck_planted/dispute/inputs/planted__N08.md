# planted:N08

## Target paper (ground truth; the gate never saw this)
- Title: Phase Separation in a system of Brownian inertial rotors
- arXiv: 2609.21556 (v1 2026-09-18); authors: Lucio Mauro Carenza, Pasquale Digregorio, Massimiliano Semeraro, Antonio Suma, Ignacio Pagonabarraga, Giuseppe Gonnella
- Behaviour the authors report as new: Purely repulsive inertial dumbbells driven by constant torques phase-separate into a dilute liquid and a dense hexatic phase (binodal, spinodal, t^1/3 coarsening), with rigidly rotating clusters; inertia is essential.
- Authors' novelty statement: "We uncover an emergent phase separation in a system of inertial dumbbells put into persistent rotation by a constant applied torque, breaking chiral symmetry, while interacting with each other through repulsive hard-core interactions."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: about 2.1 million dumbbells (two identical beads of diameter 1 joined by a stiff spring) in a square periodic 2D box at area fraction 0.65. Rules, Langevin dynamics with mass: (1) every pair of beads repels on contact (purely repulsive, no attraction); (2) each dumbbell feels a constant torque, equal and opposite forces perpendicular to its axis on its two beads, all spinning the same way, strength set by Peclet number 5 (swept); (3) beads have mass and weak friction (momentum persists; Prandtl number 0.005) plus thermal noise at temperature 0.01. Start homogeneous and watch the local density field.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Other torque strengths above Peclet number about 1, at area fractions inside the coexistence window.

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - onset_torque_peclet = 1.0  (torque strength above which the uniform state splits apart)
  - domain_growth_exponent = 0.333  (exponent of domain size growth with time after splitting)
  - coexisting_states = 2  (coexisting states: dilute liquid and dense six-fold ordered state)
  - example_area_fraction = 0.65  (global area fraction of the illustrated split state)

This behaviour appears across 1 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Rotation-induced phase separation (RIPS) in inertial chiral spinner fluids; with dumbbells, liquid-hexatic coexistence
- **behaviour_documented**: True
- **notes**: Torque-only inertial repulsive spinners phase-separating out of a uniform state is documented as RIPS (disks). Twist not found in the literature: the dense phase here is hexatic/six-fold ordered (RIPS reports a dense chiral liquid with packing near 1 and spin near 0), and coarsening keeps going with exponent 1/3 instead of the arrested phase separation RIPS reports in some regime. The hexatic order probably comes from the passive dumbbell liquid-hexatic transition, so it is plausible but likely minor. 1/3 is the generic conserved-coarsening value. The Pe~1 onset was not checked against a readable source. Artifact checks: low temperature 0.01 plus weak friction means strong heating from the spinning, so rule out a jammed dense state and fit the exponent over a long time range.
- **citations**: ["Digregorio, Pagonabarraga, Vega Reyes, Phys. Rev. Lett. 136, 218301 (2026), doi:10.1103/dbgp-pqsh", "Caporusso, Gonnella, Levis, Phys. Rev. Lett. 132, 168201 (2024)", "Transverse Self-Propulsion Enhances the Aggregation of Active Dumbbells, Entropy 27, 692 (2025), doi:10.3390/e27070692", "Caporusso et al., Soft Matter 20, 4208 (2024)"]

### N_blocked_r2 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Rotation-induced phase separation (RIPS) in an inertial chiral spinner fluid; dense phase hexatic as in dumbbell liquid-hexatic coexistence
- **behaviour_documented**: True
- **notes**: Same verdict as the first pass, with stronger evidence. van Zuiden et al. 2016 simulated essentially the same particles (inertial repulsive dimers with force-dipole torque) but reported crystal, liquid and jammed phases, not torque-driven dilute/dense coexistence with coarsening. Disk RIPS (2026) gives the mechanism and predicts it should occur in chiral fluids generally. Twists not found in the literature: (a) torque turns into translation through dumbbell shape, with no friction between particles; (b) the dense phase is hexatic, not a dense chiral liquid with spin near zero; (c) uninterrupted t^(1/3) coarsening, while disk RIPS shows arrest in some regimes. Things to check: whether 0.333 is assumed or fitted over a short window, whether the Pe=1 onset is just a sweep grid point, and whether torque-shifted effective temperature produces equilibrium-like liquid-hexatic coexistence instead of RIPS. Only abstracts and snippets were read.
- **citations**: ["Digregorio, Pagonabarraga, Vega Reyes, Phase Separation in a Chiral Active Fluid of Inertial Self-Spinning Disks, PRL 136, 218301 (2026), doi:10.1103/dbgp-pqsh", "van Zuiden, Paulose, Irvine, Bartolo, Vitelli, Spatiotemporal order and emergent edge currents in active spinner materials, PNAS 113, 12919 (2016)", "Transverse Self-Propulsion Enhances the Aggregation of Active Dumbbells, Entropy 27, 692 (2025), doi:10.3390/e27070692", "Gonnella, Marenduzzo, Suma, Tiribocchi, Motility-induced phase separation and coarsening in active matter, C. R. Physique (2015)"]

### N_blocked_r3 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Rotation-induced phase separation (RIPS) in inertial chiral spinner fluids, with liquid-hexatic coexistence as in active/passive dumbbell systems
- **behaviour_documented**: True
- **notes**: Dilute-liquid / dense-hexatic splitting driven only by spin, with inertia, is documented (RIPS for frictional spinner disks; hexatic dense phase for dumbbells). Twist: rotation is turned into motion by the dumbbell shape instead of contact friction, and coarsening continues as t^(1/3), while RIPS disks and colloidal spinners show arrested coarsening. The continued coarsening matters if it holds at late times. Onset Pe~1 and the 1/3 exponent were not verified against any source; 1/3 is the generic exponent for conserved-density coarsening. Checks to run: Pe=0 should be uniform at phi=0.65 (passive dumbbells have liquid-hexatic coexistence), and measure kinetic temperatures per phase because T=0.01 is small next to torque-driven heating.
- **citations**: ["Phase Separation in a Chiral Active Fluid of Inertial Self-Spinning Disks, PRL 136, 218301 (2026), doi:10.1103/dbgp-pqsh, arXiv:2504.08533", "Digregorio et al., Transverse Self-Propulsion Enhances the Aggregation of Active Dumbbells, Entropy 27, 692 (2025)", "Cugliandolo, Digregorio, Gonnella, Suma, Phase Coexistence in Two-Dimensional Passive and Active Dumbbell Systems, PRL 119, 268002 (2017)", "Caporusso et al., Phase separation kinetics and cluster dynamics in two-dimensional active dumbbell systems, Soft Matter 20, 4208 (2024)", "Massana-Cid et al., Arrested phase separation in chiral fluids of colloidal spinners, PRR 3, L042021 (2021)", "Translational and rotational temperature difference in coexisting phases of inertial active dumbbells, J. Chem. Phys. 164, 151101"]
