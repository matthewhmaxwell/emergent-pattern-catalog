"""P15 (GoL) — Class C declared N/A under Phase-2a spec v1.1.

Conway's Game of Life is a deterministic CA: the canonical positive
(R-pentomino at L=64) is a fixed initial condition with no parameter
that gates whether the pattern emerges. Per v1.1 spec §"Class C N/A
list", P15 has no parameter regime that suppresses the canonical
persistent-computation behavior, so Class C is N/A.

The panel runs Classes A and B only; the PASS criterion adjusts
proportionally per the v1.1 spec.
"""

from __future__ import annotations

from typing import Any, Dict, List


CONFIG: Dict[str, Any] = {
    "substrate_id": "P15_gol_class_c",
    "format": "grid",
    "status": "N/A",
    "n_a_reason": (
        "Game of Life is fully deterministic; the canonical positive "
        "(R-pentomino at L=64) is a fixed initial condition with no "
        "parameter that gates pattern emergence (v1.1 spec §Class C N/A list)."
    ),
    "regimes": [],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    raise RuntimeError(
        "P15 Class C is declared N/A under v1.1; build_substrate must not be called."
    )
