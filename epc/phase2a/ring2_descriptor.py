"""Ring-2 descriptor — the per-observation fingerprint for the genuine-novelty search.

Runs the named emergence indicator + the model-free bridge + the ADMITTED Tier-2 lens battery
into ONE labeled feature dict. Each lens self-guards on substrate (and on series length),
returning None when it does not apply, so the fingerprint is substrate-aware and sparse:
`lenses_fired` records the coordinate subspace this observation actually occupies. This is the
vector the Stage-2 search uses to (a) classify via the named lenses, (b) flag novelty via the
tripwire, and (c) cluster leads in lens space.

ADMITTED battery (21 lenses): the original 6 (structure_factor, persistent_homology[pos+field],
graph_structure, directed_info_flow, fractal_dimension[scoped], vorticity) + the comprehensive
build-out 15 (defect_census, pattern_symmetry, velocity_order, metastability, coarsening,
anomalous_transport, sync_texture, chaos_dimension, compressibility, hyperuniformity,
cross_substrate, multifractal, critical_slowing, nonequilibrium_current, extreme_value).

INTEGRITY: every battery lens is a COORDINATE (fingerprint axis), never a classifier. The
tripwire's `classified` consults ONLY generic_emergence + detector matches, NOT these lenses, so
adding them cannot move leads -- they only enrich fingerprints. Deferred/retired lenses
(recurrence -> subsumed by chaos_dimension; correlation-length -> FSS harness; DFA/ageing/PAC/
hierarchical -> documented LOW defers) are intentionally NOT wired.
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
from epc.metrics.fractal_dimension import fractal_dimension
from epc.metrics.vorticity import rotational_order
# --- comprehensive build-out (15 new coordinate lenses) ---
from epc.metrics.defect_census import defect_census
from epc.metrics.pattern_symmetry import pattern_symmetry
from epc.metrics.velocity_order import velocity_order
from epc.metrics.metastability import metastability
from epc.metrics.coarsening import coarsening
from epc.metrics.anomalous_transport import anomalous_transport
from epc.metrics.sync_texture import sync_texture
from epc.metrics.chaos_dimension import chaos_dimension
from epc.metrics.compressibility import compressibility
from epc.metrics.hyperuniformity import hyperuniformity
from epc.metrics.cross_substrate import cross_substrate
from epc.metrics.multifractal import multifractal
from epc.metrics.critical_slowing import critical_slowing
from epc.metrics.nonequilibrium_current import nonequilibrium_current
from epc.metrics.extreme_value import extreme_value

# Stable key order so descriptors vectorise alignably across observations.
SCHEMA: List[str] = [
    "em_score", "em_kind", "mf_C", "mf_struct", "mf_psi", "mf_complex",
    "sk_peak", "h1_max", "ph_components", "field_loop_area", "field_loops",
    "fractal_dim", "lacunarity", "ang_mom",
    "degree_cv", "modularity", "clustering",
    "dte_mean_te", "dte_directionality",
    # comprehensive build-out coordinates
    "defect_density", "defect_coherence",
    "lattice_symmetry", "azimuthal_concentration",
    "velcorr_length", "polarization", "mill_strength",
    "n_macrostates", "dwell_cv",
    "growth_exponent", "coarsen_fit",
    "msd_exponent", "ergodicity_break",
    "local_R_spread", "global_sync",
    "zero_one_K", "determinism",
    "lz_ratio",
    "number_variance_exponent",
    "cross_coherence",
    "mf_width",
    "ews_ar1_rise", "ews_var_rise",
    "circulation",
    "extremal_index",
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
        if ph.get("kind") == "field":
            f["field_loop_area"] = ph["field_loop_area"]; f["field_loops"] = ph["field_loops"]
        else:
            f["h1_max"] = ph["h1_max"]; f["ph_components"] = ph.get("n_components")
        fired.append("persistent_homology")
    fd = _safe(fractal_dimension, history)
    if fd:
        f["fractal_dim"] = fd["fractal_dim"]; f["lacunarity"] = fd.get("lacunarity")
        fired.append("fractal_dimension")
    vo = _safe(rotational_order, history)
    if vo:
        f["ang_mom"] = vo["ang_mom"]; fired.append("vorticity")
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

    # --- comprehensive build-out battery (each self-guards -> None off-substrate/too short) ---
    dc = _safe(defect_census, history)
    if dc:
        f["defect_density"] = dc["defect_density"]; f["defect_coherence"] = dc.get("coherence")
        fired.append("defect_census")
    ps = _safe(pattern_symmetry, history)
    if ps:
        f["lattice_symmetry"] = ps["lattice_symmetry"]; f["azimuthal_concentration"] = ps["azimuthal_concentration"]
        fired.append("pattern_symmetry")
    vlo = _safe(velocity_order, history)
    if vlo:
        f["velcorr_length"] = vlo["velcorr_length"]; f["polarization"] = vlo["polarization"]
        f["mill_strength"] = vlo["mill_strength"]; fired.append("velocity_order")
    msb = _safe(metastability, history)
    if msb:
        f["n_macrostates"] = msb["n_macrostates"]; f["dwell_cv"] = msb["dwell_cv"]
        fired.append("metastability")
    co = _safe(coarsening, history)
    if co:
        f["growth_exponent"] = co["growth_exponent"]; f["coarsen_fit"] = co["fit_quality"]
        fired.append("coarsening")
    at = _safe(anomalous_transport, history)
    if at:
        f["msd_exponent"] = at["msd_exponent"]; f["ergodicity_break"] = at["ergodicity_break"]
        fired.append("anomalous_transport")
    stx = _safe(sync_texture, history)
    if stx:
        f["local_R_spread"] = stx["local_R_spread"]; f["global_sync"] = stx["global_sync"]
        fired.append("sync_texture")
    cd = _safe(chaos_dimension, history)
    if cd:
        f["zero_one_K"] = cd["zero_one_K"]; f["determinism"] = cd["determinism"]
        fired.append("chaos_dimension")
    cz = _safe(compressibility, history)
    if cz:
        f["lz_ratio"] = cz["lz_ratio"]; fired.append("compressibility")
    hu = _safe(hyperuniformity, history)
    if hu:
        f["number_variance_exponent"] = hu["number_variance_exponent"]; fired.append("hyperuniformity")
    cs = _safe(cross_substrate, history)
    if cs:
        f["cross_coherence"] = cs["cross_coherence"]; fired.append("cross_substrate")
    mfr = _safe(multifractal, history)
    if mfr:
        f["mf_width"] = mfr["mf_width"]; fired.append("multifractal")
    csd = _safe(critical_slowing, history)
    if csd:
        f["ews_ar1_rise"] = csd["ews_ar1_rise"]; f["ews_var_rise"] = csd["ews_var_rise"]
        fired.append("critical_slowing")
    nq = _safe(nonequilibrium_current, history)
    if nq:
        f["circulation"] = nq["circulation"]; fired.append("nonequilibrium_current")
    ev = _safe(extreme_value, history)
    if ev:
        f["extremal_index"] = ev["extremal_index"]; fired.append("extreme_value")

    tw = _safe(novelty_tripwire, history, metadata) or {}
    f["tripped"] = tw.get("tripped"); f["classified"] = tw.get("classified")

    f["lenses_fired"] = fired
    return f
