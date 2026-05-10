"""P18 (voter) — Class C declared N/A under Phase-2a spec v1.1.

The voter model on a finite lattice always reaches consensus eventually:
no parameter regime suppresses the canonical pattern from non-trivial
initial conditions. v1.0's Sprint 30 prototype (this file pre-revision)
used high-bias initial conditions (init_fraction ∈ [0.93, 0.999]) as a
proxy and found ~60% of those substrates still trip the detector — they
are true positives in disguise rather than failed-regime negatives.

Per v1.1 spec §"Class C N/A list", P18 joins P15 (GoL) and P31 (Zhang
sorting) as patterns whose Class C is N/A. The panel runs Classes A
and B (and B' supplements) only; the PASS criterion adjusts proportionally.

The v1.0 init-fraction config is preserved in
``docs/archive/phase2a_panel_spec_v1_0.md`` and
``analysis/outputs/archive/v1_0/p18_phase2a_panel.json`` for reference.
"""

from __future__ import annotations

from typing import Any, Dict, List


CONFIG: Dict[str, Any] = {
    "substrate_id": "P18_voter_class_c",
    "format": "grid",
    "status": "N/A",
    "n_a_reason": (
        "Voter model on a finite lattice always reaches consensus eventually; "
        "no parameter regime suppresses the canonical pattern from non-trivial "
        "initial conditions (v1.1 spec §Class C N/A list)."
    ),
    "regimes": [],  # empty by design — Class C N/A
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Should never be called when ``CONFIG['status'] == 'N/A'``.

    Kept for interface symmetry with other failed-regime modules. The panel
    runner checks ``CONFIG['status']`` and skips Class C when N/A.
    """
    raise RuntimeError(
        "P18 Class C is declared N/A under v1.1; build_substrate must not be called. "
        "If the panel runner reaches this code path, the Class C N/A handling has regressed."
    )
