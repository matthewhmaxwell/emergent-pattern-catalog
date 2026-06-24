"""Metrics for field-mediated resource competition (P37).

Operates on a trajectory of frames each carrying 'abundances' (species, shape (n,))
and 'resources' (shape (k,)). The signature of Huisman-Weissing competition is
SUSTAINED non-equilibrium dynamics with multi-species coexistence: abundances keep
oscillating (do not reach a fixed point) while several species persist together.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def coexistence_oscillation(history: List[Dict[str, Any]], persist_frac: float = 0.01,
                            late_frac: float = 2.0 / 3.0) -> Optional[Dict[str, float]]:
    """Sustained-oscillation + coexistence summary over the late window.

    Returns None if the frames do not carry both 'abundances' and 'resources'
    (substrate mismatch). Otherwise:
      osc_cv      : median per-species temporal coefficient of variation over the
                    late window (0 = settled fixed point; >0 = still oscillating)
      n_persist   : number of species above persist_frac * global-peak in the late mean
      n_resources : number of limiting resources (the competition 'field')
      n_species   : total species
    """
    ab = [np.asarray(f["abundances"], dtype=float) for f in history
          if isinstance(f, dict) and "abundances" in f]
    res = [np.asarray(f["resources"], dtype=float) for f in history
           if isinstance(f, dict) and "resources" in f]
    if len(ab) < 6 or not res:
        return None
    A = np.array(ab)                                   # (T, n)
    n_res = int(len(res[0]))
    late = A[int(len(A) * late_frac):]
    peak = float(A.max()) or 1.0
    persist_mask = late.mean(axis=0) > persist_frac * peak
    n_persist = int(persist_mask.sum())
    cvs = []
    for i in np.where(persist_mask)[0]:
        s = late[:, i]
        mu = float(s.mean())
        if mu > 1e-9:
            cvs.append(float(s.std() / mu))
    osc = float(np.median(cvs)) if cvs else 0.0
    return {"osc_cv": osc, "n_persist": float(n_persist),
            "n_resources": float(n_res), "n_species": float(A.shape[1])}
