# planted:N04

## Target paper (ground truth; the gate never saw this)
- Title: Topological instability and reentrant crystallization in active solids
- arXiv: 2609.09500 (v1 2026-09-08); authors: Zory Davoyan, Aditya Jha, Anton Souslov
- Behaviour the authors report as new: A 2D solid with transverse (odd) pair forces melts by defect-pair proliferation at any temperature above a threshold activity, and recrystallises at higher activity even at temperatures where the passive crystal is liquid.
- Authors' novelty statement: "However, at a threshold value of non-reciprocity, we discover a qualitatively distinct mechanism for 2D melting, which we call topological instability."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: 10,000 slightly elongated particles (ellipse half-axes 0.55 and 0.65) in a periodic 2D box at density 0.87 (lattice spacing about 1.15). Rules, overdamped Brownian dynamics at temperature T: (1) each pair interacts through an orientation-dependent Lennard-Jones potential of strength e_p (repulsive core, weak attraction, depending on both particles' orientations, which can rotate); (2) each pair also feels a non-conservative force with Lennard-Jones magnitude and strength e_a, directed perpendicular to the line joining them, always in the same rotational sense; (3) thermal noise acts on positions and orientations, orientations relaxing faster. Sweep the ratio e_a/e_p from 0 to 33 and T from 0 to 6; measure six-fold bond order and defect counts.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  (none distinct)

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - activity_disordered_at_zero_temperature = 2.0  (force ratio e_a/e_p where solid is disordered even at T=0)
  - activity_of_restored_order = 33.0  (force ratio e_a/e_p where crystalline order returns)
  - temperature_of_restored_order = 6.0  (temperature at which that high-activity solid stays ordered)
  - passive_partial_melting_temperature = 4.25  (temperature where the passive solid is only six-fold oriented)

This behaviour appears across 1 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r2 (strict pass): verdict KNOWN-VARIANT
- **standard_name**: Odd/transverse-force (chiral Lennard-Jones) melting of a 2D crystal into rotating hexatic patches, on top of KTHNY two-step passive melting
- **behaviour_documented**: True
- **notes**: Melting of the solid by transverse forces (disordered at e_a/e_p~2 even at T=0) and the passive hexatic at T~4.25 are documented. Twist: crystalline order returning at e_a/e_p~33 and T=6, above the passive hexatic temperature, was not found in the literature. It is suspicious, probably caused by huge per-step forces in overdamped integration (timestep), a rigidly rotating or jammed state, or an order-parameter artefact. Treat it as unverified until checked with a smaller timestep, different system sizes, mean-squared displacement and overlap checks. The particle anisotropy (ellipses with half-axes 0.55/0.65) seems to change details only, not the phenomenology.
- **citations**: ["Caporusso, Gonnella, Levis, Phys. Rev. Lett. 132, 168201 (2024) - Phase coexistence and edge currents in the chiral Lennard-Jones fluid", "Bililign et al., Nature Physics (2021) - Motile dislocations knead odd crystals into whorls", "Braverman et al., PNAS 122 (2025) doi:10.1073/pnas.2412993122 - Melting of nonreciprocal solids", "Reentrant melting of scarred odd crystals by self-shear, Nat. Commun. (2026) s41467-026-68510-4", "Crystallization of chiral active Brownian particles at low densities, Phys. Rev. Research (2025) doi:10.1103/yntx-ss9r", "Frenkel & McTague, PRL 42, 1632 (1979) - 2D LJ hexatic"]
