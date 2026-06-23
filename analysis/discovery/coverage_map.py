"""Ring 1 — dimensional-coverage map of the 32 catalog patterns.

Maps each catalog pattern onto the 11-dimension emergence ontology
(docs/ontology_v0_4.md) and computes the EMPTY cells — ontology values (and
selected value-combinations) that no catalogued pattern occupies. These empty
cells are the grounded, "periodic-table" discovery targets for Ring 1: regions
the ontology says *could* host an emergent pattern but the catalog doesn't cover.

The per-pattern classification below is a STRUCTURAL classification grounded in
the canonical models (Schelling, Vicsek, Kuramoto, sandpile, ...) and the model
metadata. It is cross-validated against the ontology's own stated gap analysis
(docs/ontology_v0_4.md §Dimensional Gap Analysis) and paper §5.8.
"""
from __future__ import annotations

from collections import Counter
from typing import Dict, List

# Ontology value sets (docs/ontology_v0_4.md). Abbreviated codes used below.
DIMS: Dict[str, List[str]] = {
    "scale":     ["local", "meso", "global"],
    "temporal":  ["steady", "transient", "oscillatory", "propagating"],
    "interact":  ["direct", "stigmergic", "field", "none_entropy"],
    "substrate": ["continuous", "lattice", "fixed_net", "evolving_net", "well_mixed"],
    "homog":     ["identical", "heterogeneous"],
    "goal":      ["goalless", "implicit", "explicit"],
    "feedback":  ["positive", "negative", "mixed", "nontransitive"],
    "memory":    ["memoryless", "local", "env_trace"],
    "conflict":  ["cooperative", "competitive", "mixed_motive", "nontransitive"],
    "driving":   ["autonomous", "ext_forced", "periodic", "field_coupled"],
    "update":    ["motion", "state", "reproduction", "resource_exchange"],
}
ORDER = list(DIMS.keys())

# Per-pattern classification: pX -> tuple in ORDER order.
P = {
 "P1_schelling":      ("meso","steady","direct","lattice","heterogeneous","implicit","positive","local","cooperative","autonomous","motion"),
 "P2_mips":           ("meso","steady","direct","continuous","identical","goalless","positive","memoryless","cooperative","autonomous","motion"),
 "P3_turing":         ("meso","steady","field","continuous","identical","goalless","mixed","memoryless","cooperative","autonomous","state"),
 "P4_territoriality": ("meso","steady","stigmergic","continuous","identical","implicit","negative","env_trace","competitive","autonomous","motion"),
 "P5_flocking":       ("global","steady","direct","continuous","identical","goalless","positive","memoryless","cooperative","autonomous","motion"),
 "P6_milling":        ("global","steady","direct","continuous","identical","goalless","mixed","memoryless","cooperative","autonomous","motion"),
 "P7_lanes":          ("meso","steady","direct","continuous","heterogeneous","implicit","mixed","memoryless","mixed_motive","autonomous","motion"),
 "P8_traffic":        ("meso","propagating","direct","lattice","identical","implicit","negative","memoryless","competitive","autonomous","motion"),
 "P9_sync":           ("global","oscillatory","direct","fixed_net","heterogeneous","goalless","positive","local","cooperative","autonomous","state"),
 "P10_chimera":       ("global","oscillatory","direct","fixed_net","identical","goalless","mixed","local","cooperative","autonomous","state"),
 "P11_lotka_volterra":("meso","oscillatory","direct","lattice","heterogeneous","goalless","mixed","memoryless","competitive","autonomous","reproduction"),
 "P12_rps":           ("meso","propagating","direct","lattice","heterogeneous","goalless","nontransitive","memoryless","nontransitive","autonomous","reproduction"),
 "P13_excitable":     ("meso","propagating","direct","lattice","identical","goalless","mixed","local","cooperative","autonomous","state"),
 "P14_soc_sandpile":  ("meso","transient","direct","lattice","identical","goalless","positive","local","cooperative","ext_forced","state"),
 "P15_gol":           ("meso","propagating","direct","lattice","identical","goalless","mixed","memoryless","cooperative","autonomous","state"),
 "P16_hopfield":      ("global","transient","direct","fixed_net","identical","explicit","mixed","local","cooperative","ext_forced","state"),
 "P17_gradient_sense":("meso","steady","field","continuous","identical","implicit","mixed","memoryless","cooperative","field_coupled","motion"),
 "P18_voter":         ("meso","transient","direct","lattice","identical","goalless","positive","memoryless","cooperative","autonomous","state"),
 "P19_leadership":    ("global","steady","direct","continuous","heterogeneous","implicit","positive","local","cooperative","autonomous","motion"),
 "P20_quorum":        ("global","transient","stigmergic","well_mixed","identical","implicit","positive","env_trace","cooperative","autonomous","state"),
 "P21_polarization":  ("meso","steady","direct","fixed_net","identical","goalless","positive","local","mixed_motive","autonomous","state"),
 "P22_cascade":       ("global","propagating","direct","fixed_net","identical","implicit","positive","local","cooperative","autonomous","state"),
 "P23_minority_game": ("global","steady","direct","well_mixed","heterogeneous","explicit","negative","local","competitive","autonomous","state"),
 "P24_homeostasis":   ("local","steady","direct","well_mixed","identical","explicit","negative","local","cooperative","ext_forced","state"),
 "P25_canalization":  ("global","transient","direct","fixed_net","identical","implicit","mixed","local","cooperative","autonomous","state"),
 "P26_stoch_resonance":("local","oscillatory","field","well_mixed","identical","goalless","positive","memoryless","cooperative","periodic","state"),
 "P27_spatial_recip": ("meso","steady","direct","lattice","identical","implicit","mixed","local","mixed_motive","autonomous","reproduction"),
 "P28_wealth_condens":("global","steady","direct","well_mixed","identical","goalless","positive","local","competitive","autonomous","resource_exchange"),
 "P29_trail_network": ("meso","steady","stigmergic","continuous","identical","implicit","positive","env_trace","cooperative","autonomous","motion"),
 "P30_autopoiesis":   ("meso","steady","direct","continuous","heterogeneous","implicit","mixed","local","cooperative","autonomous","state"),
 "P31_delayed_grat":  ("global","transient","direct","well_mixed","identical","explicit","mixed","local","mixed_motive","autonomous","state"),
 "P32_specialization":("global","steady","direct","well_mixed","heterogeneous","implicit","negative","local","cooperative","ext_forced","state"),
}

# Plausible 2-way combination cells to probe (physically realizable, untested).
COMBOS = [
    ("substrate", "evolving_net", None, None),          # whole substrate value empty
    ("interact", "none_entropy", None, None),           # entropy/depletion-driven ordering
    ("driving", "ext_forced", "feedback", "nontransitive"),  # driven cyclic dominance
    ("update", "resource_exchange", "homog", "heterogeneous"),  # heterogeneous-endowment exchange
    ("interact", "field", "conflict", "competitive"),   # field-mediated competition
    ("memory", "env_trace", "update", "state"),         # stigmergic state-only (no motion)
    ("substrate", "fixed_net", "update", "reproduction"),  # evolutionary games on networks
]


def main() -> None:
    n = len(P)
    print(f"# Ring-1 dimensional-coverage map  ({n} catalog patterns x 11 dimensions)\n")
    print("## First-order coverage: occupancy per ontology value (EMPTY = gap)\n")
    gaps_1st = []
    for d in ORDER:
        col = [P[p][ORDER.index(d)] for p in P]
        cnt = Counter(col)
        cells = []
        for v in DIMS[d]:
            c = cnt.get(v, 0)
            cells.append(f"{v}={c}" + ("  <-- EMPTY" if c == 0 else ""))
            if c == 0:
                gaps_1st.append((d, v))
        print(f"- **{d}**: " + " | ".join(cells))
    print("\n## EMPTY first-order cells (single-value gaps):")
    for d, v in gaps_1st:
        print(f"  - {d} = {v}")
    print("\n## Targeted combination cells (empty + plausible):")
    for d1, v1, d2, v2 in COMBOS:
        if d2 is None:
            hits = [p for p in P if P[p][ORDER.index(d1)] == v1]
            label = f"{d1}={v1}"
        else:
            hits = [p for p in P if P[p][ORDER.index(d1)] == v1 and P[p][ORDER.index(d2)] == v2]
            label = f"{d1}={v1} AND {d2}={v2}"
        status = "EMPTY <--" if not hits else f"occupied by {hits}"
        print(f"  - {label}: {status}")


if __name__ == "__main__":
    main()
