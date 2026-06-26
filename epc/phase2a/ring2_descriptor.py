"""Ring-2 descriptor — the per-observation fingerprint for the genuine-novelty search.

Runs the named emergence indicator + the model-free bridge + the four ADMITTED Tier-2
lenses (structure_factor, persistent_homology, graph_structure, directed_info_flow) into
ONE labeled feature dict. Each lens self-guards on substrate, returning None when it does
not apply, so the fingerprint is substrate-aware and sparse: `lenses_fired` records the
coordinate subspace this observation actually occupies. This is the vector the Stage-2
search uses to (a) classify via the named lenses, (b) flag novelty via the tripwire, and
(c) cluster leads in lens space.

Deferred lenses (fractal_dimension, recurrence) are intentionally NOT wired in — only
lenses that earned their place contribute to the live descriptor. directed_info_flow is
fed the micro_macro component series and is interpreted per its documented scope
(coupling magnitude reliable; directionality = drive, not spatial propagation).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from epc.phase2a.emergence import generic_emergence
from epc.phase2a.novelty_tripwire import model_free_complexity, novelty_tripwire
from epc.phase2a.info_channels import micro_macro
from epc.metrics.structure_factor import structure_factor_peak
from epc.metrics.persistent_homology import persistent_homology
from epc.metrics.graph_structure import graph_structure
from epc.metrics.directed_info_flow import directed_transfer_entropy

# Stable key order so descriptors vectorise alignably across observations.
SCHEMA: List[str] = [
    "em_score", "em_kind", "mf_C", "mf_struct", "mf_psi", "mf_complex",
    "sk_peak", "h1_max", "ph_components",
    "degree_cv", "modularity", "clustering",
    "dte_mean_te", "dte_directionality",
    "tripped", "classified",
]


def _safe(fn, *a, **k):
    try:
        return fn(*a, **k)
    except Exception:
        return None


def _r(x):
    try:
        return round(float(x), 4) if x is not None and x == x else None
    except Exception:
        return None


def ring2_descriptor(history: List[Dict[str, Any]],
                     metadata: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    f: Dict[str, Any] = {k: None for k in SCHEMA}
    fired: List[str] = []

    em = _safe(generic_emergence, history) or {}
    f["em_score"] = _r(em.get("score", 0.0)); f["em_kind"] = em.get("kind")

    mf = _safe(model_free_complexity, history) or {}
    f["mf_C"] = _r(mf.get("C")); f["mf_struct"] = _r(mf.get("struct"))
    f["mf_psi"] = _r(mf.get("psi")); f["mf_complex"] = mf.get("is_complex")

    sf = _safe(structure_factor_peak, history)
    if sf:
        f["sk_peak"] = sf["sk_peak"]; fired.append("structure_factor")
    ph = _safe(persistent_homology, history)
    if ph:
        f["h1_max"] = ph["h1_max"]; f["ph_components"] = ph.get("n_components")
        fired.append("persistent_homology")
    gs = _safe(graph_structure, history)
    if gs:
        f["degree_cv"] = gs["degree_cv"]; f["modularity"] = gs["modularity"]
        f["clustering"] = gs["clustering"]; fired.append("graph_structure")

    mm = _safe(micro_macro, history)
    micro = mm[0] if isinstance(mm, tuple) else None
    if micro is not None:
        dte = _safe(directed_transfer_entropy, micro)
        if dte:
            f["dte_mean_te"] = dte["mean_te"]; f["dte_directionality"] = dte["directionality"]
            fired.append("directed_info_flow")

    tw = _safe(novelty_tripwire, history, metadata) or {}
    f["tripped"] = tw.get("tripped"); f["classified"] = tw.get("classified")

    f["lenses_fired"] = fired
    return f
