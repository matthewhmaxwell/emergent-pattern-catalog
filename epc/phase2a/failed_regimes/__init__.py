"""Class C — pattern-specific failed-regime substrates for the Phase-2a panel.

Each config exports:

- ``CONFIG``: a dict with keys ``substrate_id``, ``format`` (``"grid"`` or
  ``"phases"``), ``description``, and ``regimes`` (list of 10 dicts each
  with ``label`` + ``params`` + ``seed``).
- ``build_substrate(regime: dict) -> list[dict]``: produces the substrate
  trajectory for a single regime entry.

The harness reads ``CONFIG["regimes"]`` and calls ``build_substrate`` for
each regime to feed the detector.
"""

from epc.phase2a.failed_regimes import p18_voter, p9_kuramoto, p12_rps, p13_gh, p11_lotka_volterra
from epc.phase2a.failed_regimes import p5_vicsek, p2_active_brownian, p6_dorsogna
from epc.phase2a.failed_regimes import p8_nagel_schreckenberg, p10_chimera

REGISTRY = {
    "P18": p18_voter,
    "P9": p9_kuramoto,
    "P12": p12_rps,
    "P13": p13_gh,
    "P11": p11_lotka_volterra,
    "P5": p5_vicsek,
    "P2": p2_active_brownian,
    "P6": p6_dorsogna,
    "P8": p8_nagel_schreckenberg,
    "P10": p10_chimera,
}

__all__ = [
    "REGISTRY",
    "p18_voter", "p9_kuramoto", "p12_rps", "p13_gh", "p11_lotka_volterra",
    "p5_vicsek", "p2_active_brownian", "p6_dorsogna",
    "p8_nagel_schreckenberg", "p10_chimera",
]
