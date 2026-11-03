# planted:N11

**Reported behavior (what the authors report as new):** Adding a three-body competition term to a purely competitive spatial birth-death model breaks the hexagonal-spots-only restriction, giving stripes and gaps and an unusual spots-stripes-spots sequence.

**Rule on the card:** World: individuals at continuous positions in a square periodic domain, starting with 4,000 placed at random. Rules (continuous-time birth and death): (1) each individual gives birth at rate 1920; the offspring lands at a Gaussian offset (standard deviation 0.05); (2) each individual dies at rate 192 + (density of others within radius 1.25) + g times (density within radius 0.0625, squared: a three-body crowding term). Densities are counts divided by disk area. Only competitive (harmful) interactions, no attraction or help. Sweep g from 3.3e-6 to 3.9e-5 and look at the steady spatial pattern.

**Observed features on the card:**
  - anisotropy_at_stripe_onset = 1.0 (pattern anisotropy index jumps from 0 to this value)
  - g_for_hexagonal_spots = 3.27e-06 (three-body strength where clusters form a hexagonal spot array)
  - g_for_stripes = 3.14e-05 (three-body strength where clusters merge into parallel stripes)
  - g_for_gaps = 3.93e-05 (three-body strength where empty holes sit in a populated background)
  - reentry_death_birth_ratio = 0.55 (death/birth ratio where spots reappear beyond the stripes)

## Reply R36

**JSON fields:** {"verdict": "KNOWN-VARIANT", "standard_name": "Competition-induced (nonlocal competition) Turing-like pattern formation: hexagonal spots -> stripes -> gaps in a nonlocal logistic / spatial birth-death model", "notes": "Pairwise nonlocal-competition patterns (spots/stripes/gaps) are well documented. No source found for a short-range three-body (density-squared) death term as the control parameter, nor for the specific re-entrant spots-stripes-spots sequence vs death/birth ratio; the latter is consistent with generic hexagon/stripe/inverted-hexagon theory but this was not found stated for this model. Did not read full texts of the Physica D or Communications Physics papers; numeric g values cannot be checked against the literature."}

**Reply text:**

**Verdict: KNOWN-VARIANT, medium confidence.** The model family and the main pattern sequence are well documented. I found no source for the specific three-body mechanism or the re-entrant spots.

**1. Model family**
- This is a spatial birth–death point process, an individual-based model in continuous space. It has Gaussian dispersal and a death rate that rises with the density of neighbours within a set radius. It is the spatial logistic model with nonlocal competition (Bolker–Pacala / Dieckmann–Law). A short-range term quadratic in density is added.
- Its continuum limit is a nonlocal Fisher–KPP-type equation. Pattern formation here is a finite-wavelength (Turing-like) instability driven only by competition, with no attraction or facilitation.

**2. What the literature describes**
- **Earliest clear source.** Hernández-García & López (2004) introduced a particle model with birth and death rates that depend on the number of neighbours. It shows a quasistationary pattern with many clusters of particles arranged periodically in space. They also derived the macroscopic density equation. A linear stability analysis shows a finite-wavelength instability that leads to pattern formation. Verified: Phys. Rev. E 70, 016216, published 30 July 2004.
- **Competition alone is enough.** Martínez-García, Calabrese, Hernández-García & López (2013) showed that patterns appear despite the absence of short-range facilitation, provided competition is strong enough. They also noted that a death term that increases with nonlocal density through a competition kernel also gives rise to pattern formation. Verified: Geophysical Research Letters, 2013.
- **Review.** Martínez-García et al. (2014), "Minimal mechanisms for vegetation patterns in semiarid regions," Phil. Trans. R. Soc. A 372:20140068. It presents a series of models that yield spatial patterns when finite-range competition is the only driving force. Its figures show both stripes and spots.
- **Pattern sequence.** Later work states that non-local vegetation models show spots, labyrinths, stripes and gaps in two dimensions. Spots with hexagonal symmetry are described as typical of competition-induced patterns. The spots → stripes → gaps order is also the classic Turing/Swift–Hohenberg sequence. The ordering of hexagons, stripes and inverted hexagons depends on the sign and size of the quadratic nonlinearity.
- **Individual-based versus continuum.** Surendran, Pinto-Ramos, Menezes & Martínez-García (2025, Physica D 477:134703) compare individual-based simulations of such models with spatial-moment and density-field approximations. I confirmed the title, authors, year and venue from the search results. I did not read the full text.

**3. What differs, and whether it matters**
- I found no source that uses a short-range three-body (density-squared) death term together with long-range pairwise competition. I also found none that sweeps its strength to get the spot → stripe → gap sequence.
  - Mechanistically this looks like a higher-order saturating nonlinearity added to the standard nonlocal-competition instability. That would likely change which planform wins (hexagons, stripes or inverted hexagons) but not the underlying instability. This is my inference, not something I found stated.
- The re-entrant sequence (spots → stripes → spots as the death/birth ratio varies near 0.55) is not documented for this model. It is qualitatively consistent with hexagon → stripe → inverted-hexagon transitions in amplitude-equation theory. I did not find that stated for this model, and I did not confirm whether the re-entrant spots are inverted spots (gaps).
- The specific numbers (g values, the 0.55 ratio, the anisotropy jump from 0 to 1) are properties of this parameterisation. The literature gives no numbers to compare them against.

**
