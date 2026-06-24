"""Ring-2 model-free novelty tripwire — the Tier-2-special BRIDGE.

Combines the lens-agnostic complexity measures (MPR statistical complexity + Psi_CE
causal emergence) with the named-lens classification, and flags the one quadrant that
matters for the genuine-novelty search:

    COMPLEX (real structure beyond noise) AND UNCLASSIFIED (no named lens or
    catalogued detector recognises it).

That is the signal that the instrument is looking at structure it has no category for
-- a Ring-2 novelty lead, and the trigger to invent a Tier-3 lens (docs/validation_
rebuild/ring2_strategy.md). The model-free measures do NOT name the structure; they
only assert it is there. The named layer (generic_emergence + optional catalogued
detectors) says whether the instrument already has a word for it.

Thresholds are NULL-CALIBRATED: the model-free measures must exceed what the null
systems (noise / random walk / random graph) produce. Re-derive them with
analysis/ring2/stage0_bridge_validation.py if the null set changes.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from epc.phase2a.info_channels import micro_macro, mpr_complexity, psi_ce_best
from epc.phase2a.emergence import generic_emergence

# Null-calibrated thresholds (Stage-0 baseline, 2026-06-24; analysis/ring2/
# stage0_bridge_validation.py). On the null set (noise / random walk / random graph)
# the model-free measures topped out at C=0.127, psi=-2.07; floors are set a margin
# above. Validated QUIET on the known corpus: 0/3 nulls tripped, 17/17 known-emergent
# classified, 0 complex-but-unclassified residual blind spots.
C_THR = 0.16      # MPR statistical complexity floor (max-null C 0.127 + margin)
PSI_THR = 0.05    # Psi_CE causal-emergence floor (max-null psi << 0; floored at 0.05)
EM_THR = 0.50     # named-lens generic_emergence score that counts as "classified"


def model_free_complexity(history: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Lens-agnostic 'is there structure beyond noise here?' — returns C, psi, and a
    bool. Names nothing. Returns is_complex=False when the history carries no usable
    per-component series (micro_macro could not extract one)."""
    micro, cands = micro_macro(history)
    if micro is None or cands is None:
        return {"C": float("nan"), "psi": float("nan"), "is_complex": False, "macro_feat": None}
    mpr = mpr_complexity(cands.get("mean"))
    psi, feat = psi_ce_best(micro, cands)
    C = float(mpr["C"])
    psi = float(psi) if psi == psi else 0.0          # NaN -> 0
    return {"C": C, "psi": psi, "macro_feat": feat,
            "is_complex": bool(C > C_THR or psi > PSI_THR)}


def novelty_tripwire(history: List[Dict[str, Any]], metadata: Optional[Dict[str, Any]] = None,
                     detector_fns: Optional[Dict[str, Any]] = None,
                     battery: Optional[Any] = None, em_thr: float = EM_THR) -> Dict[str, Any]:
    """Triage one observation for the Ring-2 search.

    Returns a dict; the decisive field is `tripped` (complex AND unclassified). When
    detector_fns is supplied, a catalogued-detector MATCH also counts as classified
    (the stronger check used in the live search loop); without it, classification is
    by the named generic-emergence lenses alone (the fast baseline)."""
    mf = model_free_complexity(history)
    em = generic_emergence(history)
    em_score = float(em.get("score", 0.0)); em_kind = em.get("kind")

    named_classified = em_score >= em_thr
    matched = None
    if detector_fns is not None:
        from analysis.battery_profile import profile_observation
        prof = profile_observation(history, metadata, battery=battery,
                                   detector_fns=detector_fns, match_min_tier="confirmation")
        if prof["verdict"]["verdict"] == "MATCH":
            matched = prof["verdict"]["pattern_id"]
    classified = bool(named_classified or matched is not None)
    tripped = bool(mf["is_complex"] and not classified)

    return {"tripped": tripped, "is_complex": mf["is_complex"], "classified": classified,
            "matched": matched, "C": round(mf["C"], 4) if mf["C"] == mf["C"] else None,
            "psi": round(mf["psi"], 4), "em_score": round(em_score, 4), "em_kind": em_kind,
            "macro_feat": mf["macro_feat"],
            "reason": ("complex+unclassified (novelty lead)" if tripped
                       else f"classified ({matched or em_kind})" if classified
                       else "below complexity floor (trivial/null)")}
