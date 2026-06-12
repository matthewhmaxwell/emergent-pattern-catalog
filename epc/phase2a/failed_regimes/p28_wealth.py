"""P28 (wealth condensation) — Class C failed regimes + lookalikes.

Each regime is a wealth-distribution time series (frames carrying a 1-D
``wealth`` array) that should NOT reach DEFINITIVE P28 condensation:

  1. saving_propensity (Chakraborti-Chakrabarti 2000): conserved pairwise
     exchange WITH a saving propensity lambda>0 -> the distribution
     stabilizes at a Gamma with BOUNDED Gini (~0.3). Rejected at screening
     (Gini below the super-Boltzmann floor).
  2. redistribution (wealth tax): conserved exchange WITH a redistribution
     term chi>0 -> Gini plateaus at a bounded value (~0.5-0.6). Rejected at
     screening (Gini floor) / confirmation (no monotonic growth to extremes).
  3. bouchaud_mezard: multiplicative growth (geometric) with weak linear
     redistribution -> high Gini power-law inequality, but produced by a
     NON-CONSERVED multiplicative process (total wealth drifts). This is the
     key lookalike: it LOOKS like condensation (high Gini) but the mechanism
     is different. Rejected at screening by the CONSERVATION gate.
  4. boltzmann_exchange: additive symmetric conserved exchange -> the
     Boltzmann-Gibbs exponential equilibrium (Gini ~0.5). Rejected at
     screening (Gini below the super-Boltzmann floor); this is the null the
     yard-sale condensation must beat.

Reference: Dragulescu & Yakovenko 2000 (Boltzmann-Gibbs); Chakraborti &
Chakrabarti 2000 (saving propensity); Bouchaud & Mezard 2000 (multiplicative
redistribution / power law).
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


CONFIG: Dict[str, Any] = {
    "substrate_id": "P28_wealth_condensation_class_c",
    "format": "scalar_wealth",
    "status": "populated",
    "regimes": [
        {"label": "saving_propensity",
         "params": {"n_agents": 1000, "n_frames": 100, "n_tx_per_frame": 2000,
                    "lambda_save": 0.5, "chi": 0.0}, "seed": 7},
        {"label": "redistribution",
         "params": {"n_agents": 1000, "n_frames": 100, "n_tx_per_frame": 2000,
                    "lambda_save": 0.0, "chi": 0.05, "redistribute_every": 1000}, "seed": 7},
        {"label": "bouchaud_mezard",
         "params": {"n_agents": 1000, "n_frames": 100, "sigma": 0.12,
                    "redistribution_J": 0.02}, "seed": 7},
        {"label": "boltzmann_exchange",
         "params": {"n_agents": 1000, "n_frames": 60, "n_tx_per_frame": 4000},
         "seed": 7},
    ],
}


def _yardsale_history(p: Dict[str, Any], seed: int) -> List[Dict[str, Any]]:
    from epc.models.yard_sale import YardSale
    m = YardSale(
        n_agents=p["n_agents"], f=0.5,
        lambda_save=p.get("lambda_save", 0.0),
        chi=p.get("chi", 0.0),
        redistribute_every=p.get("redistribute_every", 1000),
        init_mode="equal", seed=seed,
    )
    s0 = m.setup()
    # Record the t=0 equal state so condensation EMERGENCE (Gini grows from ~0)
    # is visible to the detector.
    frames: List[Dict[str, Any]] = [
        {"wealth": np.asarray(s0["wealth"], dtype=np.float64).copy()}
    ]
    for _ in range(p["n_frames"]):
        s = m.step(n_transactions=p["n_tx_per_frame"])
        frames.append({"wealth": np.asarray(s["wealth"], dtype=np.float64).copy()})
    return frames


def _boltzmann_exchange_history(p: Dict[str, Any], seed: int) -> List[Dict[str, Any]]:
    """True ADDITIVE symmetric conserved exchange (Dragulescu-Yakovenko 2000).

    Pick two agents, pool their wealth, split it uniformly at random. Total is
    conserved; the equilibrium is the Boltzmann-Gibbs EXPONENTIAL distribution
    (Gini ~ 0.5) — emerges from equality but PLATEAUS well below the
    super-Boltzmann condensation floor. This is the null the multiplicative
    yard-sale condensation must beat, and a faithful conserved-but-not-
    condensing lookalike.
    """
    rng = np.random.default_rng(seed)
    N = p["n_agents"]
    w = np.full(N, 1.0, dtype=np.float64)
    frames: List[Dict[str, Any]] = [{"wealth": w.copy()}]
    n_tx = p["n_tx_per_frame"]
    for _ in range(p["n_frames"]):
        for _ in range(n_tx):
            i = int(rng.integers(0, N))
            j = int(rng.integers(0, N))
            if i == j:
                continue
            pool = w[i] + w[j]
            frac = rng.random()
            w[i] = frac * pool
            w[j] = (1.0 - frac) * pool
        frames.append({"wealth": w.copy()})
    return frames


def _bouchaud_mezard_history(p: Dict[str, Any], seed: int) -> List[Dict[str, Any]]:
    """Multiplicative growth with weak linear redistribution (NON-conserved).

    w_i <- w_i * exp(sigma * xi - sigma^2/2) + J*(<w> - w_i)

    The multiplicative kick is mean-preserving in expectation but variance
    accumulates -> heavy-tailed (power-law-ish) inequality. With only weak
    redistribution the realised total wealth DRIFTS (not exactly conserved),
    which is exactly what separates this mechanism from yard-sale condensation.
    """
    rng = np.random.default_rng(seed)
    N = p["n_agents"]
    sigma = p["sigma"]
    J = p["redistribution_J"]
    w = np.full(N, 1.0, dtype=np.float64)
    frames: List[Dict[str, Any]] = [{"wealth": w.copy()}]
    steps_per_frame = 50
    for _ in range(p["n_frames"] - 1):
        for _ in range(steps_per_frame):
            xi = rng.normal(0.0, 1.0, size=N)
            w = w * np.exp(sigma * xi - 0.5 * sigma * sigma)
            w = w + J * (w.mean() - w)
            w = np.clip(w, 1e-12, None)
        frames.append({"wealth": w.copy()})
    return frames


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    label = regime["label"]
    p = regime["params"]
    seed = regime["seed"]
    if label == "bouchaud_mezard":
        return _bouchaud_mezard_history(p, seed)
    if label == "boltzmann_exchange":
        return _boltzmann_exchange_history(p, seed)
    # saving / redistribution are conserved yard-sale (multiplicative) variants
    return _yardsale_history(p, seed)
