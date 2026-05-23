"""Class B' synthetic structured-but-non-pattern substrates (Phase-2a v1.1).

When a pattern's substrate type has fewer than 3 catalog mates the panel
supplements Class B with substrate-typed structured non-pattern substrates.
Each builder returns the substrate in the detector's *native* consumed
format (oscillator → phases dict; network → grid dict) so no cross-format
adapter is needed — the whole point of the v1.1 substrate-typed redesign.

Contracts:
- Each builder takes a positional ``seed: int`` and any number of keyword
  arguments with sensible defaults for fast in-test invocation.
- Each builder returns ``list[dict]`` matching the format used by Class A
  generators in :mod:`epc.phase2a.synthetic`.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List

import numpy as np


# --- Oscillator supplements -------------------------------------------------

def incoherent_phases(
    seed: int,
    *,
    n: int = 300,
    n_steps: int = 600,
    cadence: int = 10,
    **_: Any,
) -> List[Dict[str, Any]]:
    """Uniform random phases redrawn each step. Order parameter r ≈ 1/√n.

    No coupling, no temporal autocorrelation. A correctly-specific
    synchronization detector should not fire.
    """
    rng = np.random.default_rng(seed)
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        theta = rng.uniform(0.0, 2.0 * np.pi, n)
        r = float(np.abs(np.mean(np.exp(1j * theta))))
        out.append({"theta": theta, "r": r, "step": t * cadence})
    return out


def subcritical_kuramoto(
    seed: int,
    *,
    K: float = 0.1,
    n: int = 200,
    n_steps: int = 6000,
    cadence: int = 10,
    gamma: float = 0.5,
    dt: float = 0.05,
    **_: Any,
) -> List[Dict[str, Any]]:
    """Real Kuramoto integrator at K well below K_c (default K=0.1, K_c=1.0).

    Oscillators do not synchronize. Useful as a Class B' supplement for
    non-Kuramoto oscillator detectors (e.g. P10 chimera). Overlap with P9's
    Class C is intentional in v1.1 — the supplements are substrate-typed
    structured negatives, not detector-specific.
    """
    from epc.models.kuramoto import KuramotoModel, KuramotoParams

    model = KuramotoModel(KuramotoParams(N=n, K=K, gamma=gamma, dt=dt, seed=seed))
    return model.run(n_steps=n_steps, record_every=cadence)


# --- Network supplements ----------------------------------------------------

def random_graph_evolution(
    seed: int,
    *,
    rows: int = 64,
    cols: int = 64,
    n_steps: int = 200,
    **_: Any,
) -> List[Dict[str, Any]]:
    """Cells as i.i.d. Bernoulli(0.5) at every timestep — no opinion update.

    Same shape as a voter grid but no neighborhood interaction. P18 should not
    fire because there is no coarsening dynamic.
    """
    rng = np.random.default_rng(seed)
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        grid = rng.integers(0, 2, size=(rows, cols), dtype=np.int8)
        out.append({"grid": grid, "grid_dims": (rows, cols), "step": t})
    return out


def network_random_walks(
    seed: int,
    *,
    rows: int = 64,
    cols: int = 64,
    n_steps: int = 200,
    **_: Any,
) -> List[Dict[str, Any]]:
    """Each cell independently flips with p=0.5 each step. No graph coupling.

    Time-correlated noise at the cell level but no spatial structure.
    Distinguishable from random_graph_evolution by per-cell autocorrelation.
    """
    rng = np.random.default_rng(seed)
    state = rng.integers(0, 2, size=(rows, cols), dtype=np.int8)
    flips = rng.integers(0, 2, size=(n_steps, rows, cols), dtype=np.int8)
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        state = np.bitwise_xor(state, flips[t])
        out.append({"grid": state.copy(), "grid_dims": (rows, cols), "step": t})
    return out


# --- Continuous-2d supplements ----------------------------------------------

def uncorrelated_random_walks(
    seed: int,
    *,
    n: int = 200,
    box_size: float = 16.0,
    speed: float = 0.03,
    n_steps: int = 200,
    **_: Any,
) -> List[Dict[str, Any]]:
    """N particles with i.i.d. uniform headings redrawn each step — no coupling.

    Polarization r ≈ 1/√N. Same output shape as Vicsek/ABP history.
    A correctly-specific flocking or milling detector should not fire.
    """
    rng = np.random.default_rng(seed)
    positions = rng.uniform(0.0, box_size, size=(n, 2))
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        headings = rng.uniform(-np.pi, np.pi, size=n)
        velocities = speed * np.column_stack([np.cos(headings), np.sin(headings)])
        positions = (positions + velocities) % box_size
        out.append({
            "positions": positions.copy(),
            "headings": headings.copy(),
            "velocities": velocities.copy(),
            "speeds": np.full(n, speed, dtype=np.float64),
            "step": t,
        })
    return out


def independent_brownian_motion(
    seed: int,
    *,
    n: int = 200,
    box_size: float = 16.0,
    sigma: float = 0.1,
    n_steps: int = 200,
    **_: Any,
) -> List[Dict[str, Any]]:
    """N independent Brownian particles — no interactions.

    Each particle takes an i.i.d. Gaussian displacement per step.
    Tests against trivial diffusive motion that could spuriously trigger
    collective-motion detectors.
    """
    rng = np.random.default_rng(seed)
    positions = rng.uniform(0.0, box_size, size=(n, 2))
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        displacements = rng.normal(0.0, sigma, size=(n, 2))
        positions = (positions + displacements) % box_size
        speeds = np.linalg.norm(displacements, axis=1)
        headings = np.arctan2(displacements[:, 1], displacements[:, 0])
        out.append({
            "positions": positions.copy(),
            "headings": headings.copy(),
            "velocities": displacements.copy(),
            "speeds": speeds,
            "step": t,
        })
    return out


# --- Registry ---------------------------------------------------------------

SUPPLEMENTS_BY_SUBSTRATE_TYPE: Dict[str, List[str]] = {
    "oscillator": ["incoherent_phases", "subcritical_kuramoto"],
    "network": ["random_graph_evolution", "network_random_walks"],
    "continuous_2d": ["uncorrelated_random_walks", "independent_brownian_motion"],
}

SUPPLEMENT_BUILDERS: Dict[str, Callable[..., List[Dict[str, Any]]]] = {
    "incoherent_phases": incoherent_phases,
    "subcritical_kuramoto": subcritical_kuramoto,
    "random_graph_evolution": random_graph_evolution,
    "network_random_walks": network_random_walks,
    "uncorrelated_random_walks": uncorrelated_random_walks,
    "independent_brownian_motion": independent_brownian_motion,
}
