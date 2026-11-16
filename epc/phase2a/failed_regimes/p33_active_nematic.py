"""P33 (active nematic, +-1/2 defect gas) -- Class C failed regimes.

The target's own model, ``epc.models.active_nematic.active_nematic_field``, at
settings where the pattern (a locally ordered director field carrying a
SUSTAINED +-1/2 defect gas) does not form. Chosen from the model's own physical
observables (local nematic order S_loc, half-defect density over short and long
runs), not from the detector verdict:
analysis/validation_rebuild/p33_p37_class_c_groundtruth.py ->
analysis/outputs/p33_p37_class_c_groundtruth.json.

  uniform_nematic eta=0.05, aligned start   ordered nematic, ZERO defects at 400
                                            and 4000 steps (standalone negative
                                            neg_uniform_nematic, reused unchanged)
  ordered_nematic eta=0.10, aligned start   same at higher activity: S_loc 0.96,
                                            zero defects at 400 and 4000 steps
  isotropic eta=0.60 / eta=0.80, random     activity noise destroys local order:
                                            S_loc 0.19 / 0.16 vs isotropic
                                            baseline 0.15 (the disordered phase)

Excluded on ground truth (documented, not run): low activity from a RANDOM start
(eta 0-0.10) is NOT a non-forming regime in this model -- the +-1/2 defects do not
anneal (half-defect density 0.012-0.145 still at 4000 steps, S_loc 0.29-0.68), so
the observation contains the pattern's signature. eta 0.30-0.45 is a crossover
(S_loc 0.27-0.35), neither clearly ordered nor isotropic.

Two seeds per regime (200-207), distinct from the positive seeds (0-4).
"""
from __future__ import annotations

from typing import Any, Dict, List

CONFIG: Dict[str, Any] = {
    "substrate_id": "P33_active_nematic_class_c",
    "format": "orientation_field",
    "status": "populated",
    "regimes": [
        {"label": "uniform_nematic_eta=0.05_aligned",
         "params": {"builder": "neg_uniform_nematic"}, "seed": 200},
        {"label": "uniform_nematic_eta=0.05_aligned",
         "params": {"builder": "neg_uniform_nematic"}, "seed": 201},
        {"label": "ordered_nematic_eta=0.10_aligned",
         "params": {"eta": 0.10, "init_mode": "aligned", "G": 64, "n_steps": 400, "n_frames": 60},
         "seed": 202},
        {"label": "ordered_nematic_eta=0.10_aligned",
         "params": {"eta": 0.10, "init_mode": "aligned", "G": 64, "n_steps": 400, "n_frames": 60},
         "seed": 203},
        {"label": "isotropic_eta=0.60",
         "params": {"eta": 0.60, "init_mode": "random", "G": 64, "n_steps": 400, "n_frames": 60},
         "seed": 204},
        {"label": "isotropic_eta=0.60",
         "params": {"eta": 0.60, "init_mode": "random", "G": 64, "n_steps": 400, "n_frames": 60},
         "seed": 205},
        {"label": "isotropic_eta=0.80",
         "params": {"eta": 0.80, "init_mode": "random", "G": 64, "n_steps": 400, "n_frames": 60},
         "seed": 206},
        {"label": "isotropic_eta=0.80",
         "params": {"eta": 0.80, "init_mode": "random", "G": 64, "n_steps": 400, "n_frames": 60},
         "seed": 207},
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    from epc.models.active_nematic import active_nematic_field, neg_uniform_nematic
    p = regime["params"]
    if p.get("builder") == "neg_uniform_nematic":
        return neg_uniform_nematic(regime["seed"])[0]
    return active_nematic_field(regime["seed"], G=p["G"], eta=p["eta"], n_steps=p["n_steps"],
                                n_frames=p["n_frames"], init_mode=p["init_mode"])[0]
