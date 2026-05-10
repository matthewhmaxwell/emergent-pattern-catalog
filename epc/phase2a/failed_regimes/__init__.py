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

from epc.phase2a.failed_regimes import p18_voter, p9_kuramoto

REGISTRY = {
    "P18": p18_voter,
    "P9": p9_kuramoto,
}

__all__ = ["REGISTRY", "p18_voter", "p9_kuramoto"]
