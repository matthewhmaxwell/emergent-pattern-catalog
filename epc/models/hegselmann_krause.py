"""Hegselmann-Krause bounded-confidence opinion dynamics — canonical model for P21.

Reference:
    Hegselmann, R. & Krause, U. (2002). "Opinion dynamics and bounded confidence:
    models, analysis, and simulation." Journal of Artificial Societies and Social
    Simulation, 5(3), 2.

Update rule (synchronous):
    x_i(t+1) = (1/|N_i(t)|) Σ_{j ∈ N_i(t)} x_j(t)

    where N_i(t) = {j : |x_i(t) - x_j(t)| ≤ ε} (including i itself).

    ε (epsilon): confidence bound. Agent i only considers opinions within
    distance ε of its own.

Key results (Hegselmann & Krause 2002):
    ε ≥ 0.5:  consensus (single cluster)
    ε ≈ 0.2:  2-3 clusters (polarization)
    ε ≈ 0.1:  many clusters (fragmentation)
    ε → 0:    frozen (no interaction)

    Critical threshold: ε_c ≈ 0.25-0.27 for consensus/polarization transition
    (exact value depends on N and initial distribution).

State history keys:
    opinions   : np.ndarray (N,) — agent opinions in [0, 1]
    step       : int — interaction round
    n_clusters : int — number of distinct opinion clusters
    variance   : float — opinion variance
    dip_stat   : float — Hartigan's dip statistic (multimodality)
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np


class HegselmannKrauseModel:
    """Hegselmann-Krause bounded-confidence opinion dynamics (2002).

    Parameters
    ----------
    n_agents : int
        Number of agents.
    epsilon : float
        Confidence bound. Agents only average with others within ε.
    init_mode : str
        'uniform' (U[0,1]), 'gaussian' (N(0.5, 0.2) clipped),
        'bimodal' (mixture of two Gaussians).
    convergence_tol : float
        Stop early if max opinion change < tol.
    seed : int
        Random seed.
    """

    def __init__(
        self,
        n_agents: int = 500,
        epsilon: float = 0.2,
        init_mode: str = "uniform",
        convergence_tol: float = 1e-8,
        seed: int = 42,
    ):
        self.n_agents = n_agents
        self.epsilon = epsilon
        self.init_mode = init_mode
        self.convergence_tol = convergence_tol
        self.seed = seed

        self.opinions: np.ndarray = np.empty(0)
        self.rng = np.random.default_rng(seed)
        self._step_count: int = 0
        self._converged: bool = False

    def setup(self) -> dict[str, Any]:
        """Initialize opinions."""
        self.rng = np.random.default_rng(self.seed)
        self._step_count = 0
        self._converged = False
        N = self.n_agents

        if self.init_mode == "uniform":
            self.opinions = self.rng.uniform(0, 1, size=N)
        elif self.init_mode == "gaussian":
            self.opinions = np.clip(self.rng.normal(0.5, 0.2, size=N), 0, 1)
        elif self.init_mode == "bimodal":
            half = N // 2
            left = self.rng.normal(0.3, 0.05, size=half)
            right = self.rng.normal(0.7, 0.05, size=N - half)
            self.opinions = np.clip(np.concatenate([left, right]), 0, 1)
            self.rng.shuffle(self.opinions)
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode}")

        return self._state_dict()

    def step(self) -> dict[str, Any]:
        """One synchronous update round."""
        if self._converged:
            self._step_count += 1
            return self._state_dict()

        N = self.n_agents
        eps = self.epsilon
        x = self.opinions

        # Compute pairwise distances and find neighbors within epsilon
        # For moderate N (≤5000), direct approach is fast enough
        new_x = np.empty(N)
        for i in range(N):
            diffs = np.abs(x - x[i])
            neighbors = diffs <= eps  # includes self
            new_x[i] = x[neighbors].mean()

        max_change = np.max(np.abs(new_x - self.opinions))
        self.opinions = new_x
        self._step_count += 1

        if max_change < self.convergence_tol:
            self._converged = True

        return self._state_dict()

    def run(self, n_steps: int = 500) -> list[dict[str, Any]]:
        """Run for n_steps or until convergence."""
        history = [self.setup()]
        for _ in range(n_steps):
            history.append(self.step())
            if self._converged:
                break
        return history

    def _state_dict(self) -> dict[str, Any]:
        """Build state dictionary."""
        x = self.opinions
        n_clusters = _count_clusters(x, gap=self.epsilon / 2)
        var = float(np.var(x))

        return {
            "opinions": x.copy(),
            "step": self._step_count,
            "n_clusters": n_clusters,
            "variance": var,
            "converged": self._converged,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata."""
        return {
            "model": "hegselmann_krause",
            "epsilon": self.epsilon,
            "n_agents": self.n_agents,
            "init_mode": self.init_mode,
            "substrate": "opinion_space",
        }


def _count_clusters(opinions: np.ndarray, gap: float = 0.05) -> int:
    """Count opinion clusters by gap detection.

    Sort opinions, then count groups separated by gaps > threshold.
    """
    if len(opinions) == 0:
        return 0
    sorted_x = np.sort(opinions)
    gaps = np.diff(sorted_x) > gap
    return int(np.sum(gaps)) + 1
