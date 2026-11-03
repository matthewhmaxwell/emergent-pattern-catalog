# planted:N16

## Target paper (ground truth; the gate never saw this)
- Title: Self-Avoidance Governs Re-entrant Coil-Globule Transition of Inertial Active Brownian Ring Polymers
- arXiv: 2610.16188 (v1 2026-10-12); authors: Rakesh Palariya, Arindam Panda, Roland G. Winkler, Sunil P Singh
- Behaviour the authors report as new: A self-avoiding ring of active Brownian beads with rotational inertia collapses into a dense globule (over tenfold smaller) at intermediate activity with outward-polarised propulsion, then re-swells at high activity; translational inertia opposes collapse.
- Authors' novelty statement: "Using Langevin dynamics simulations of an active Brownian ring polymer with translational and rotational inertia, we uncover an inertia-driven re-entrant coil–globule transition in the absence of any attractive interactions."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: one closed loop of 2048 beads (a chain joined end to end) in open 3D space. Rules, Langevin dynamics with inertia: (1) consecutive beads are linked by springs and all beads repel at contact (no attraction), so the chain cannot cross itself; (2) each bead pushes itself with constant force along its own heading (activity: Peclet number Pe, swept); (3) headings diffuse randomly but with rotational inertia (reduced moment of inertia J swept 0.001-100, e.g. J=5); (4) beads have translational inertia (reduced mass M=0.3), friction and thermal noise. Measure loop size (radius of gyration) versus Pe, J and bead number.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Same rules with J=30: very sharp switch back to an open coil within a narrow Pe window.
  - Same rules with J=10: compact loop size nearly independent of bead number.

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - size_reduction_factor = 10  (squared loop size shrinks by over this factor versus passive)
  - activity_of_max_compaction = 17  (Peclet number giving the most compact loop)
  - reswelling_activity = 40  (Peclet number above which the loop swells back open)
  - size_vs_length_exponent_j10 = 0  (exponent of loop size versus bead number at J=10 (near zero))

This behaviour appears across 3 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r1 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Motility-induced collapse of an underdamped (inertial) active Brownian particle polymer; non-monotonic Rg(Pe) of active rings
- **behaviour_documented**: True
- **notes**: Agrees with the first-pass KNOWN-VARIANT label but corrects its note: varying rotational inertia J in an inertial self-propelled-bead chain IS documented (Buglakov et al. 2025, linear chains; collapse at lower Pe for large J and m, MIPS-like core-plus-loops state). Remaining differences: (1) closed-loop topology (minor); (2) loop size independent of bead number (exponent ~0) at J=10, which is not documented and is physically implausible for a non-crossing chain at finite density (a globule gives ~1/3). This most likely comes from beads overlapping or the chain passing through itself at high inertia/Pe, unconverged runs, or a narrow fitting range. Verify with minimum non-neighbour bead distances, a check that the loop's knot type is unchanged, time-step sensitivity and convergence before treating it as real. The collapse and reswelling themselves look genuine, not trivial.
- **citations**: ["Buglakov, Lelecova, Chertovich, Phys. Chem. Chem. Phys. 27, 1691 (2025), doi:10.1039/d4cp03731f; ChemRxiv doi:10.26434/chemrxiv-2024-sv81t (underdamped ABP chain, m and J varied, collapse without attraction)", "Role of rotational inertia for collective phenomena in active matter (PubMed 36200385) (rotational inertia promotes MIPS)", "A computer simulation study of a chiral active ring polymer (PubMed 39513442) (non-monotonic Rg vs active force for ABP ring)", "Effects of Knotting on the Collapse of Active Ring Polymers, Macromolecules 2025, doi:10.1021/acs.macromol.5c02097", "Local Polar and Long-Range Isotropic Activity Assisted Swelling and Collapse Dynamics of an Active Ring Polymer, Macromolecules 2023, doi:10.1021/acs.macromol.3c00384"]

### N_blocked_r2 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Motility-induced collapse and re-swelling of an active Brownian particle polymer with translational and rotational inertia (here with ring topology)
- **behaviour_documented**: True
- **notes**: Collapse then re-swelling with activity in an inertial ABP chain is documented for linear chains with parameters close to these (Pe_min ~15 vs observed ~17). Twists: ring topology (likely minor); sharp J=30 re-opening (inertia shifting the collapse is documented qualitatively); size-vs-N exponent ~0 at J=10, which is physically implausible for a chain whose beads cannot overlap (a dense globule should give at least ~1/3). The exponent most likely reflects a narrow N range, runs not at steady state, or Rg dominated by loops; treat it as a possible artifact rather than a discovery. PCCP full text was not opened; Buglakov numbers come from preprint snippets.
- **citations**: ["Buglakov, Lelecova, Chertovich, Phys. Chem. Chem. Phys. 27, 1691 (2025), doi:10.1039/D4CP03731F (linear inertial ABP chain: collapse without attraction, Rg minimum ~Pe 15 at N=1024, re-swelling at high Pe, J and m dependence)", "ChemRxiv preprint 'Monomer activity solely induces collapse of polymer chain in a good solvent' (Buglakov et al. 2024)", "Locatelli et al., Phys. Rev. Lett. 126, 097801 (2021) (active ring collapse for long rings, tangential activity)", "Rep. Prog. Phys. review doi:10.1088/1361-6633/ad3e11 (ABP rings behave similarly to linear chains)"]

### N_blocked_r3 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Motility-induced collapse (and high-activity re-swelling) of an underdamped active Brownian polymer; non-monotonic Rg of active Brownian rings
- **behaviour_documented**: True
- **notes**: Collapse, re-swelling and inertia dependence are documented for the linear underdamped ABP chain; non-monotonic Rg is documented for overdamped ABP rings. The ring + inertia combination appears unstudied, but that is a rule-combination gap, not new phenomenology. Not found in the literature: a size-vs-N exponent of ~0 at J=10. It is physically suspicious (a dense globule would give ~1/3) and should be checked for incomplete equilibration, a narrow N range, or a measurement artifact before it is treated as real. The sharp J-dependent re-opening is a quantitative refinement of documented inertia effects. Full texts were not read.
- **citations**: ["Buglakov, Lelecova, Chertovich, Phys. Chem. Chem. Phys. 27, 1691 (2025), doi:10.1039/D4CP03731F (ChemRxiv 10.26434/chemrxiv-2024-sv81t) - underdamped ABP chain with m and J, collapse without attraction, re-swelling, inertia lowers collapse Pe", "Behavior of Active Polymer Knots, Macromolecules 58, 11229 (2025) - non-monotonic Rg vs Pe for unknotted ABP rings", "Locatelli, Bianco, Malgaretti, PRL 126, 097801 (2021) - activity-induced collapse of long active rings (tangential activity)", "A computer simulation study of a chiral active ring polymer (2024, PubMed 39513442) - non-monotonic Rg of 2D active Brownian ring"]
