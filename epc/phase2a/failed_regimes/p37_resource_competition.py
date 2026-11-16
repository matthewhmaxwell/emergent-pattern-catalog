"""P37 (field-mediated resource competition, Huisman-Weissing oscillation) -- Class C.

The target's own chemostat model (``epc.models.resource_competition``: Liebig
growth on shared depletable resources, D=0.25, r=1, m=0.25) with the competitive
structure that drives sustained oscillation removed, or with the canonical
Huisman-Weissing (2001) parameter set starved. Ground truth (long-run fixed point
or exclusion, from the model's own trajectories, not the detector):
analysis/validation_rebuild/p33_p37_class_c_groundtruth.py.

  competitive_exclusion   strict dominance hierarchy, 5 species / 3 resources ->
                          1 survivor at a fixed point (standalone, reused)
  stable_coexistence      5 decoupled diagonal niches -> all 5 coexist at a FIXED
                          POINT, late CV ~1e-16 (standalone, reused; the key control)
  single_resource         one limiting resource -> R* rule, 1 survivor (standalone)
  hw_uniform_consumption  added: the canonical HW K and supply, but every species
                          consumes every resource at the same (mean) rate, i.e. the
                          intransitive consumption loop is removed -> exclusion to
                          1 species at a fixed point
  hw_low_supply           added: canonical HW K and C with all supplies S_j = 1 ->
                          exclusion to 1 species at a fixed point

Rejected candidate (ground truth says the pattern still forms): HW with high
mortality m=0.6 keeps all 5 species oscillating (late CV ~1), so it is not a
non-forming regime and is not used.

Two seeds per regime (200-209), distinct from the positive seeds (0-4).
"""
from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

_LABELS = ("competitive_exclusion", "stable_coexistence", "single_resource",
           "hw_uniform_consumption", "hw_low_supply")

CONFIG: Dict[str, Any] = {
    "substrate_id": "P37_resource_competition_class_c",
    "format": "abundance_resource_timeseries",
    "status": "populated",
    "regimes": [
        {"label": lab, "params": {"t_max": 3000.0, "dt": 0.02, "n_frames": 80}, "seed": 200 + 2 * i + k}
        for i, lab in enumerate(_LABELS) for k in (0, 1)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    from epc.models import resource_competition as rc
    lab = regime["label"]
    seed = regime["seed"]
    if lab == "competitive_exclusion":
        return rc.neg_competitive_exclusion(seed)[0]
    if lab == "stable_coexistence":
        return rc.neg_stable_coexistence(seed)[0]
    if lab == "single_resource":
        return rc.neg_single_resource(seed)[0]
    p = regime["params"]
    rng = np.random.default_rng(seed * 7919 + 1)        # same as resource_competition()
    N0 = 0.1 + 0.02 * rng.random(5)
    K, C, S = rc._HW_K.copy(), rc._HW_C.copy(), rc._HW_S.copy()
    r, m = np.full(5, 1.0), np.full(5, 0.25)
    if lab == "hw_uniform_consumption":
        C = np.full_like(C, float(C.mean()))
    elif lab == "hw_low_supply":
        S = np.ones(5)
    else:
        raise KeyError(lab)
    return rc._integrate(N0, S.copy(), K, C, S, r, m, 0.25, p["t_max"], p["dt"], p["n_frames"])
