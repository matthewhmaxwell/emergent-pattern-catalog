"""P3 (Gray-Scott) — Class C declared N/A under Phase-2a spec v1.1.

The Turing-pattern window in Gray-Scott (F, k) parameter space is a
continuous 2D region without a well-characterised binary sub-threshold
boundary analogous to Schelling's tolerance threshold or Kuramoto's K_c.

Parameter regimes outside the Turing window (e.g. uniform decay at
F=0.10, k=0.10; or uniform high at F=0.01, k=0.02) produce homogeneous
fields whose field_std ≈ 0 — these are rejected by P3's prerequisite
check (field_std >= 0.01), not by the Turing-wavelength detector itself.
Class C tests are designed to probe near-threshold detector behaviour,
not trivial prerequisite failures.

Per v1.1 spec §"Class C N/A list": patterns whose 'failed regime' is
defined as initial conditions (or parameter choices) where the pattern
simply does not form — and where the detector correctly rejects via a
substrate-level prerequisite rather than pattern-level discrimination —
qualify for Class C N/A.  P3 meets this criterion.  The panel runs
Classes A and B only; Class C is omitted.
"""

from __future__ import annotations

from typing import Any, Dict, List


CONFIG: Dict[str, Any] = {
    "substrate_id": "P3_gray_scott_non_turing",
    "format": "field",
    "status": "N/A",
    "n_a_reason": (
        "P3 Class C declared N/A: Gray-Scott non-Turing regimes (outside the "
        "F-k Turing window) produce homogeneous fields rejected by P3's "
        "field_std prerequisite, not by Turing-wavelength discrimination. "
        "No binary sub-threshold boundary comparable to Schelling/Kuramoto "
        "exists in Gray-Scott parameter space (v1.1 spec §Class C N/A list)."
    ),
    "regimes": [],  # empty by design — Class C N/A
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Should never be called when ``CONFIG['status'] == 'N/A'``.

    The panel runner checks ``CONFIG['status']`` and skips Class C when N/A.
    If called directly, raises RuntimeError to surface accidental usage.
    """
    raise RuntimeError(
        "P3 Class C is declared N/A under v1.1; build_substrate must not be "
        "called. Check panel runner logic."
    )
