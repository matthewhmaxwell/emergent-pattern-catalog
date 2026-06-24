"""Heterogeneous-agent kinetic wealth exchange — Pattern P36.

Ring-1 combination cell: update = resource_exchange AND homog = heterogeneous.

N agents hold a conserved scalar wealth and exchange it pairwise with a
Chakraborti-Chakrabarti (CC) saving propensity lambda_i per agent. With
HETEROGENEOUS lambda_i ~ U[0,1) the steady-state wealth distribution develops a
POWER-LAW (Pareto) tail, P(w) ~ w^-2 (Chatterjee, Chakrabarti & Manna 2004). The
Pareto tail is the SPECIFIC signature of agent heterogeneity, and it separates
this cell from every wealth look-alike already in the catalog / literature:

  - homogeneous (single) lambda  -> Gamma distribution, no power-law tail
  - lambda = 0 exchange          -> Boltzmann-Gibbs exponential (Dragulescu-Yakovenko)
  - pure Yard-Sale (P28)         -> condensation, Gini -> 1 (one agent holds ~all)
  - equal / lognormal            -> no power-law tail

References
----------
Chatterjee, A., Chakrabarti, B.K., Manna, S.S. (2004). Pareto law in a kinetic
    model of market with random saving propensity. Physica A 335, 155-163.
Chakraborti, A., Chakrabarti, B.K. (2000). Eur. Phys. J. B 17, 167.
Dragulescu, A., Yakovenko, V.M. (2000). Eur. Phys. J. B 17, 723.

Updates use the standard disjoint-pair "sweep" vectorization (each sweep is a
random perfect matching of the N agents; all pairs exchanged at once). This is
distribution-equivalent to sequential random pairing and converges to the same
Pareto law, but runs in O(n_sweeps) vectorized steps.
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

Built = Tuple[List[Dict[str, Any]], Dict[str, Any]]


def _cc_sweeps(seed: int, N: int, lam: np.ndarray, n_sweeps: int,
               n_frames: int, w0: float = 1.0) -> List[Dict[str, Any]]:
    """CC kinetic exchange with per-agent saving-propensity vector lam (N,)."""
    rng = np.random.default_rng(seed)
    w = np.full(N, float(w0))
    Nh = N // 2
    snap = max(1, n_sweeps // n_frames)
    frames: List[Dict[str, Any]] = [{"wealth": w.copy(), "step": 0}]
    for s in range(1, n_sweeps + 1):
        perm = rng.permutation(N)
        a = perm[:Nh]; b = perm[Nh:2 * Nh]
        eps = rng.random(Nh)
        unsaved = (1.0 - lam[a]) * w[a] + (1.0 - lam[b]) * w[b]
        w[a] = lam[a] * w[a] + eps * unsaved
        w[b] = lam[b] * w[b] + (1.0 - eps) * unsaved
        if s % snap == 0:
            frames.append({"wealth": w.copy(), "step": s})
    if frames[-1]["step"] != n_sweeps:
        frames.append({"wealth": w.copy(), "step": n_sweeps})
    return frames


def kinetic_exchange(seed: int = 0, N: int = 3000, hetero: bool = True,
                     lam_fixed: float = 0.5, n_sweeps: int = 5000,
                     n_frames: int = 40) -> Built:
    """Heterogeneous (positive) vs homogeneous (control) CC kinetic exchange.

    hetero=True : per-agent lambda_i ~ U[0, 1) -> Pareto power-law tail (P36).
    hetero=False: a single lambda = lam_fixed  -> Gamma (lam>0) / exponential (lam=0).
    """
    rng = np.random.default_rng(seed * 7919 + 1)
    if hetero:
        lam = np.clip(rng.uniform(0.0, 1.0, size=N), 0.0, 0.9999)
    else:
        lam = np.full(N, float(lam_fixed))
    frames = _cc_sweeps(seed, N, lam, n_sweeps, n_frames)
    return frames, {"model": "kinetic_exchange", "N": N, "hetero": hetero,
                    "lam_fixed": (None if hetero else lam_fixed),
                    "has_conserved_resource": True, "has_pairwise_exchange": True,
                    "heterogeneous_agents": hetero}


# ---- negatives: wealth look-alikes the Pareto-tail discriminator must reject ----

def neg_uniform_saving(seed: int = 0) -> Built:
    """Homogeneous saving propensity -> Gamma distribution (NO power-law tail).
    THE key discriminator: identical exchange rule, heterogeneity removed."""
    return kinetic_exchange(seed, hetero=False, lam_fixed=0.5)


def neg_boltzmann(seed: int = 0) -> Built:
    """lambda = 0 exchange -> Boltzmann-Gibbs exponential (Dragulescu-Yakovenko)."""
    return kinetic_exchange(seed, hetero=False, lam_fixed=0.0)


def neg_condensation(seed: int = 0, N: int = 3000, n_sweeps: int = 8000,
                     f: float = 0.1, n_frames: int = 40) -> Built:
    """Pure Yard-Sale multiplicative min-stake -> condensation (P28).
    High Gini, but a single agent holds ~all -> NOT a Pareto tail."""
    rng = np.random.default_rng(seed * 104729 + 3)
    w = np.full(N, 1.0); Nh = N // 2
    snap = max(1, n_sweeps // n_frames)
    frames: List[Dict[str, Any]] = [{"wealth": w.copy(), "step": 0}]
    for s in range(1, n_sweeps + 1):
        perm = rng.permutation(N)
        a = perm[:Nh]; b = perm[Nh:2 * Nh]
        stake = f * np.minimum(w[a], w[b])
        win_a = rng.random(Nh) < 0.5
        w[a] += np.where(win_a, stake, -stake)
        w[b] += np.where(win_a, -stake, stake)
        if s % snap == 0:
            frames.append({"wealth": w.copy(), "step": s})
    frames.append({"wealth": w.copy(), "step": n_sweeps})
    return frames, {"model": "yard_sale_condensation", "N": N}


def neg_equal(seed: int = 0, N: int = 3000) -> Built:
    """Near-equal wealth (tiny conserved jitter) -> no inequality, no tail."""
    rng = np.random.default_rng(seed * 13 + 5)
    w = np.full(N, 1.0)
    frames = [{"wealth": w.copy(), "step": 0}]
    for fr in range(1, 41):
        w = np.clip(w + rng.normal(0, 1e-4, N), 1e-9, None)
        w *= N / w.sum()
        frames.append({"wealth": w.copy(), "step": fr})
    return frames, {"model": "equal_wealth", "N": N}


def neg_lognormal(seed: int = 0, N: int = 3000) -> Built:
    """Static lognormal wealth: heavy-ish but exponential-class tail (NOT power law)."""
    rng = np.random.default_rng(seed * 31 + 7)
    frames = []
    for fr in range(20):
        w = rng.lognormal(mean=0.0, sigma=1.0, size=N)
        w *= N / w.sum()
        frames.append({"wealth": w.copy(), "step": fr})
    return frames, {"model": "lognormal_static", "N": N}
