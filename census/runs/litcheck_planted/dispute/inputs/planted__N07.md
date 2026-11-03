# planted:N07

## Target paper (ground truth; the gate never saw this)
- Title: Odd diffusion and power-law correlations in chiral mass-transport processes
- arXiv: 2609.16821 (v1 2026-09-15); authors: Koushik Das, Animesh Hazra, Punyabrata Pradhan
- Behaviour the authors report as new: Isotropic mass-transfer rules with a fixed rotational sense (L-shaped chiral moves) produce r^-4 power-law density correlations in the steady state, with long-wavelength fluctuations minimised at an optimal chirality and a square-root cusp near maximal odd diffusion.
- Authors' novelty statement: "However, the broken symmetries give rise to fundamentally new fluctuation properties."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: an 800x800 square lattice with wrap-around edges; each site holds a continuous mass (mean 1 per site). Rules, repeated: (1) pick a random site and one of the four directions at random; (2) the site keeps 40% of its mass plus a random share of the rest and sends the remainder one step in that direction; (3) the receiving site pools it with its own mass, keeps 60% plus a random share of the rest, and passes the remainder one step to its left (a 90-degree counter-clockwise turn), so mass moves along L-shaped paths that always turn the same way. Rates are identical for all four directions and total mass is conserved. Measure steady-state density-density correlations and long-wavelength fluctuations.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Same rules in a closed box with reflecting walls: a steady mass current circulates along the edges.

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - correlation_decay_exponent = 4  (density correlations decay as distance to minus this power)
  - xy_correlation_ratio = 1.0  (ratio of correlations along x and along y (isotropic))
  - fluctuation_minima_vs_turning = 1  (optimal turning strengths minimising long-wavelength density fluctuations)
  - cusp_exponent = 0.5  (square-root cusp exponent of fluctuations near maximal turning strength)

This behaviour appears across 1 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r2 (first pass): verdict KNOWN
- **standard_name**: Generic power-law (~1/r^(d+2)) density correlations from multidirectional (chiral) hopping in conserved mass transport; related to center-of-mass-conserving hyperuniform dynamics
- **model_family**: statistical physics: nonequilibrium conserved-mass transport processes (multidirectional/chiral hopping), hyperuniformity and generic scale invariance
- **notes**: r^-4 decay and isotropy are matched by Hazra et al. 2025 (d+2 with d=2). The chiral section was seen only as an abstract. The optimal-turning minima, the square-root cusp exponent and the closed-box edge current were not verified in readable text. If they are not in that paper they could be a KNOWN-VARIANT detail, but they most likely follow from the mobility-tensor amplitude. Check the paper's chiral section directly.
- **citations**: ["Hazra et al. 2025, Phys Rev E 112(4):044130 (doi 10.1103/sn1k-k8ng), 'Generic power laws in higher-dimensional lattice models with multidirectional hopping' - density-density correlations decay as 1/r^(d+2) (r^-4 in 2D), amplitude set by the Onsager/mobility matrix, can be short-ranged in a restricted regime", "ICTS programme abstract (mftheory) - multidirectional hopping covers center-of-mass-conserving and chiral dynamics, both producing algebraic correlations (abstract only)", "J. Stat. Mech. 2025, 'Hyperuniformity in mass transport processes with center-of-mass conservation: some exact results' (doi 10.1088/1742-5468/ada88c) - related hyperuniform mass transport with S(q)~q^2; authors not confirmed", "Lyman thesis, Virginia Tech - anisotropic dynamics plus a conservation law gives power-law correlations (generic scale invariance background)"]
