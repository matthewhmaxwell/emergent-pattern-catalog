"""Ring 0 discovery library — canonical emergent phenomena ABSENT from the 32-pattern
catalog.

Each entry is a well-established emergent phenomenon (cited) that is mechanistically
distinct from every catalog pattern. Run through the instrument by `sweep.py`; the
EXPECTED read is EMERGENT-UNCLASSIFIED (emergence detected, NO false MATCH) with a
named firing channel reported. Strong, clearly-distinct candidates graduate to a
full Phase-2a discriminating detector (P33+) — the same bar the catalog was held to.

This is the deliberately-curated ("rooted in the literature") arm of the discovery
sweep: the goal here is not novelty-to-science but novelty-to-catalog — it proves the
find -> verify -> catalog pipeline end-to-end and banks real catalog growth.

NULL controls are included so each sweep re-confirms specificity (NO-EMERGENCE) in
the same run, alongside the candidates.
"""
from __future__ import annotations

from typing import Any, Dict, List

# The 4 seed systems already validated as out-of-catalog in the T2c novelty arm.
from analysis.t2c_systems import (
    nov_dla, nov_keller_segel, nov_eden,
    null_spatial_noise, null_random_walk, null_uncoupled_phases, null_frozen_noise,
)
# Active nematic now uses the canonical faithful field model (the P33 model) rather
# than the original rough agent sketch — it is truly apolar with a sustained ±1/2
# defect gas, so it demonstrates the closed loop (EMERGENT-UNCLASSIFIED -> MATCH P33).
from epc.models.active_nematic import active_nematic_field

# New Ring-0 systems (this module).
from analysis.discovery.ring0_new import (
    cahn_hilliard_spinodal, gray_scott_selfrep, swarmalators,
    langton_ant, kauffman_rbn_critical,
)


# Each entry: name, family (descriptive emergence family), ref (literature anchor),
# build (seed -> (history, metadata)), stochastic (vary seed across replicates).
CANDIDATES: List[Dict[str, Any]] = [
    # --- seed set (T2c-validated out-of-catalog) ---
    {"name": "dla_fractal", "family": "fractal aggregation",
     "ref": "Witten & Sander 1981, PRL 47:1400",
     "build": nov_dla, "stochastic": True},
    {"name": "keller_segel_chemotaxis", "family": "chemotactic collapse",
     "ref": "Keller & Segel 1970, J Theor Biol 26:399",
     "build": nov_keller_segel, "stochastic": True},
    {"name": "active_nematic", "family": "nematic / topological defects",
     "ref": "Doostmohammadi et al. 2018, Nat Commun 9:3246",
     "build": active_nematic_field, "stochastic": True},
    {"name": "eden_kpz_interface", "family": "interface roughening (KPZ)",
     "ref": "Eden 1961; Kardar-Parisi-Zhang 1986, PRL 56:889",
     "build": nov_eden, "stochastic": True},
    # --- new Ring-0 systems ---
    {"name": "cahn_hilliard_spinodal", "family": "conserved phase separation",
     "ref": "Cahn & Hilliard 1958, J Chem Phys 28:258",
     "build": cahn_hilliard_spinodal, "stochastic": True},
    {"name": "gray_scott_selfrep", "family": "self-replicating spots",
     "ref": "Pearson 1993, Science 261:189",
     "build": gray_scott_selfrep, "stochastic": True},
    {"name": "swarmalators", "family": "space-phase coupling (sync+swarm)",
     "ref": "O'Keeffe, Hong & Strogatz 2017, Nat Commun 8:1504",
     "build": swarmalators, "stochastic": True},
    {"name": "langton_ant", "family": "emergent order from chaos (highway)",
     "ref": "Langton 1986, Physica D 22:120",
     "build": langton_ant, "stochastic": False},  # deterministic; seed ignored
    {"name": "kauffman_rbn_critical", "family": "edge-of-chaos Boolean dynamics",
     "ref": "Kauffman 1969, J Theor Biol 22:437",
     "build": kauffman_rbn_critical, "stochastic": True},
]

NULL_CONTROLS: List[Dict[str, Any]] = [
    {"name": "null_spatial_noise", "family": "null", "ref": None,
     "build": null_spatial_noise, "stochastic": True},
    {"name": "null_random_walk", "family": "null", "ref": None,
     "build": null_random_walk, "stochastic": True},
    {"name": "null_uncoupled_phases", "family": "null", "ref": None,
     "build": null_uncoupled_phases, "stochastic": True},
    {"name": "null_frozen_noise", "family": "null", "ref": None,
     "build": null_frozen_noise, "stochastic": True},
]

ALL: List[Dict[str, Any]] = CANDIDATES + NULL_CONTROLS
