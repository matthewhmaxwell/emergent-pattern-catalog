"""Ring 1 — gap-targeted systems (constructed from empty ontology cells).

The evolving-network gap (substrate=evolving_net) is now canonical in
`epc/models/adaptive_network.py` (P34). Re-exported here so the Ring-1 discovery
sweep + coverage docs keep a stable import. `adaptive_voter(phi=0.99)` fragments
(the P34 positive); lower phi (e.g. 0.7) develops degree heterogeneity instead.
"""
from __future__ import annotations

from epc.models.adaptive_network import adaptive_voter, null_random_rewire  # noqa: F401
