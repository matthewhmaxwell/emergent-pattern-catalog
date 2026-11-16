"""P34 (adaptive-network fission) -- Class C failed regimes.

The target's own model, ``epc.models.adaptive_network.adaptive_voter`` (N=200,
mean degree 4, 400k updates, the positive's protocol), at rewiring
probabilities where the network does NOT fission, plus the two mechanism-disabled
variants from the standalone validation. "Does not fission" is ground truth from
the final graph (giant component, components >= 10% of N, opinion purity), not
from the detector: analysis/validation_rebuild/p33_p37_class_c_groundtruth.py.
For phi <= 0.80 every run (2 seeds each, phi = 0..0.80) ends in opinion
consensus on ONE connected component (giant = 1.00); fission appears only at
phi >= 0.90 (giant 0.54-0.68, two pure-opinion components). This places the
fragmentation transition between 0.8 and 0.9 for this implementation.

  fixed_network_voter      no rewiring at all (coevolution disabled; standalone
                           neg_fixed_network_voter, reused unchanged)
  opinion_blind_rewiring   rewiring that ignores opinions (the topology moves but
                           cannot sort by opinion; standalone null_random_rewire)
  subcritical phi=0.50     coevolution present, below the transition -> consensus
  subcritical phi=0.80     coevolution present, just below the transition

Two seeds per regime (200-207), distinct from the positive seeds (0-4).
"""
from __future__ import annotations

from typing import Any, Dict, List

_AV = {"N": 200, "mean_deg": 4, "max_steps": 400000, "n_frames": 80}

CONFIG: Dict[str, Any] = {
    "substrate_id": "P34_adaptive_network_class_c",
    "format": "adjacency_timeseries",
    "status": "populated",
    "regimes": [
        {"label": "fixed_network_voter", "params": {"builder": "neg_fixed_network_voter"}, "seed": 200},
        {"label": "fixed_network_voter", "params": {"builder": "neg_fixed_network_voter"}, "seed": 201},
        {"label": "opinion_blind_rewiring", "params": {"builder": "null_random_rewire"}, "seed": 202},
        {"label": "opinion_blind_rewiring", "params": {"builder": "null_random_rewire"}, "seed": 203},
        {"label": "subcritical_phi=0.50", "params": dict(_AV, phi=0.50), "seed": 204},
        {"label": "subcritical_phi=0.50", "params": dict(_AV, phi=0.50), "seed": 205},
        {"label": "subcritical_phi=0.80", "params": dict(_AV, phi=0.80), "seed": 206},
        {"label": "subcritical_phi=0.80", "params": dict(_AV, phi=0.80), "seed": 207},
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    from epc.models import adaptive_network as an
    p = regime["params"]
    b = p.get("builder")
    if b == "neg_fixed_network_voter":
        return an.neg_fixed_network_voter(regime["seed"])[0]
    if b == "null_random_rewire":
        return an.null_random_rewire(regime["seed"])[0]
    return an.adaptive_voter(regime["seed"], N=p["N"], mean_deg=p["mean_deg"], phi=p["phi"],
                             max_steps=p["max_steps"], n_frames=p["n_frames"])[0]
