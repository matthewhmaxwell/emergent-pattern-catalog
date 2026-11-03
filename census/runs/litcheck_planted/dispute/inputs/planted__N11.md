# planted:N11

## Target paper (ground truth; the gate never saw this)
- Title: Higher-Order Competition as a Minimal Mechanism for Spatial Pattern Diversity in Population Dynamics
- arXiv: 2609.29725 (v1 2026-09-24); authors: David Pinto-Ramos, Anudeep Surendran, Ricardo Martinez-Garcia
- Behaviour the authors report as new: Adding a three-body competition term to a purely competitive spatial birth-death model breaks the hexagonal-spots-only restriction, giving stripes and gaps and an unusual spots-stripes-spots sequence.
- Authors' novelty statement: "However, previous studies of the standard two-dimensional nonlocal FKPP equation reported only hexagonal arrays of population clusters."

## Card the gate saw
```
THE RULE SYSTEM (plain words):
World: individuals at continuous positions in a square periodic domain, starting with 4,000 placed at random. Rules (continuous-time birth and death): (1) each individual gives birth at rate 1920; the offspring lands at a Gaussian offset (standard deviation 0.05); (2) each individual dies at rate 192 + (density of others within radius 1.25) + g times (density within radius 0.0625, squared: a three-body crowding term). Densities are counts divided by disk area. Only competitive (harmful) interactions, no attraction or help. Sweep g from 3.3e-6 to 3.9e-5 and look at the steady spatial pattern.

OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):
  - Same rules with g = 2.35e-5 (spots stretching into stripes).
  - Death-to-birth ratio near 0.55: as g rises, spots, then stripes, then spots again.

WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):
  - anisotropy_at_stripe_onset = 1.0  (pattern anisotropy index jumps from 0 to this value)
  - g_for_hexagonal_spots = 3.27e-06  (three-body strength where clusters form a hexagonal spot array)
  - g_for_stripes = 3.14e-05  (three-body strength where clusters merge into parallel stripes)
  - g_for_gaps = 3.93e-05  (three-body strength where empty holes sit in a populated background)
  - reentry_death_birth_ratio = 0.55  (death/birth ratio where spots reappear beyond the stripes)

This behaviour appears across 4 distinct rule programs of this structure.
```

## Gate replies that called the behaviour documented (valid, blocked-search runs)

### N_blocked_r3 (strict pass): verdict KNOWN
- **standard_name**: Competition-induced (nonlocal logistic / Brownian-bug) Turing-like pattern formation: hexagonal spots -> stripes/labyrinths -> gaps, with spots returning at higher mortality
- **behaviour_documented**: True
- **notes**: Changes the label from variant to known: the only difference is the rule-level control parameter (three-body short-range crowding), not the behaviour. Spots probably exist at g=0 from the pairwise top-hat competition (the spot threshold sits at the bottom of the sweep), so g only drives spots->stripes->gaps by raising the filled area, as with a cubic saturation term. One competition-only variant with competition in the death rate (Martínez-García 2014) did not produce gaps; the squared-density term is a plausible reason gaps appear here, and that mechanism is known from Swift-Hohenberg-type models. The 0.55 death/birth ratio is unclear (realised steady-state ratio should be about 1) and the anisotropy index looks saturated; neither looks like an artifact that undermines the patterns.
- **citations**: ["Hernández-García & López, Phys. Rev. E 70, 016216 (2004)", "Ramos, López, Hernández-García & Muñoz, Phys. Rev. E 77, 011116 (2008)", "Heinsalu, Hernández-García & López, Phys. Rev. E 85, 041105 (2012)", "Martínez-García, Calabrese, Hernández-García & López, Geophys. Res. Lett. 40, 6143 (2013)", "Martínez-García et al., Phil. Trans. R. Soc. A (2014), Minimal mechanisms for vegetation patterns in semiarid regions (PMC4186217)", "Flow spatial structure determines pattern instabilities in nonlocal models of population dynamics, Communications Physics (2025)", "Vegetation pattern formation induced by local growth outpacing susceptibility to non-local competition, J. Theor. Biol. (2026)", "Localised labyrinthine patterns in ecosystems, Sci. Rep. (2021)"]
