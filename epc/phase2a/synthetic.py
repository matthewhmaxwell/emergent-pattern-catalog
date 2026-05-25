"""Class A — synthetic null substrates for the Phase-2a panel.

Ten generators per the spec. Each returns a substrate in one of two
detector formats:

- ``"grid"``  → list[dict] each with ``{'grid': np.ndarray (H, W) int8,
  'grid_dims': (H, W)}``
- ``"phases"`` → list[dict] each with ``{'theta': np.ndarray (N,),
  'r': float, 'step': int}``

All generators are deterministic functions of (seed, format, **params)
so panel runs are reproducible. Two of the ten (``permutation_shuffled``
and ``time_shuffled``) require the detector's canonical positive as
input — these accept it via the ``positive`` kwarg.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List

import numpy as np


GRID_DEFAULT_SHAPE = (64, 64)
GRID_DEFAULT_STEPS = 200
PHASES_DEFAULT_N = 300
PHASES_DEFAULT_STEPS = 600
PHASES_DEFAULT_CADENCE = 10  # matches Kuramoto record_every=10 — keeps n_T_osc comparable
AVALANCHES_DEFAULT_N = 300   # default count of synthetic avalanche-size samples
AVALANCHES_DEFAULT_MAX_SIZE = 100


def _wrap_avalanches(sizes: np.ndarray) -> List[Dict[str, Any]]:
    """Wrap an avalanche-size array into a one-element history dict.

    P14's detector consumes a flat array of avalanche sizes rather than a
    state-history list. We package it as a single-step "history" so the
    panel runner's iteration model still works.
    """
    return [{"avalanche_sizes": sizes.astype(np.int64), "step": 0}]


def _grid_history_from_array(arr: np.ndarray, cadence: int = 1) -> List[Dict[str, Any]]:
    """Wrap a (T, H, W) array into the standard grid-history format."""
    T, H, W = arr.shape
    return [
        {"grid": arr[t].astype(np.int8), "grid_dims": (H, W), "step": t * cadence}
        for t in range(T)
    ]


def _phases_history_from_array(theta_t: np.ndarray, cadence: int = PHASES_DEFAULT_CADENCE) -> List[Dict[str, Any]]:
    """Wrap a (T, N) array of phases into the standard phases-history format.

    ``cadence`` sets the synthetic recording interval (``step *= cadence``)
    so that detectors which compute ``n_T_osc`` from consecutive ``step``
    deltas see a comparable observation density to the canonical positive.
    """
    T, _ = theta_t.shape
    out: List[Dict[str, Any]] = []
    for t in range(T):
        theta = theta_t[t]
        r = float(np.abs(np.mean(np.exp(1j * theta))))
        out.append({"theta": theta.copy(), "r": r, "step": t * cadence})
    return out


# --- Generators (10) ---------------------------------------------------------

def random_uniform_field(
    format: str,
    seed: int,
    *,
    shape: tuple = GRID_DEFAULT_SHAPE,
    n_steps: int = GRID_DEFAULT_STEPS,
    n: int = PHASES_DEFAULT_N,
    n_avalanches: int = AVALANCHES_DEFAULT_N,
    max_size: int = AVALANCHES_DEFAULT_MAX_SIZE,
    **_: Any,
) -> List[Dict[str, Any]]:
    """1. Random uniform field — i.i.d. uniform on the natural range."""
    rng = np.random.default_rng(seed)
    if format == "grid":
        arr = rng.integers(0, 2, size=(n_steps, *shape), dtype=np.int8)
        return _grid_history_from_array(arr)
    if format == "phases":
        theta_t = rng.uniform(0.0, 2.0 * np.pi, size=(n_steps, n))
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        return _wrap_avalanches(rng.integers(1, max_size + 1, size=n_avalanches))
    raise ValueError(f"unknown format: {format}")


def random_gaussian_field(
    format: str,
    seed: int,
    *,
    shape: tuple = GRID_DEFAULT_SHAPE,
    n_steps: int = GRID_DEFAULT_STEPS,
    n: int = PHASES_DEFAULT_N,
    n_avalanches: int = AVALANCHES_DEFAULT_N,
    **_: Any,
) -> List[Dict[str, Any]]:
    """2. Random Gaussian field — i.i.d. Gaussian thresholded to substrate range."""
    rng = np.random.default_rng(seed)
    if format == "grid":
        z = rng.standard_normal(size=(n_steps, *shape))
        arr = (z > 0.0).astype(np.int8)
        return _grid_history_from_array(arr)
    if format == "phases":
        z = rng.standard_normal(size=(n_steps, n))
        theta_t = np.mod(z * np.pi, 2.0 * np.pi)
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        sizes = np.clip(rng.normal(50, 15, size=n_avalanches), 1, None).astype(np.int64)
        return _wrap_avalanches(sizes)
    raise ValueError(f"unknown format: {format}")


def random_binary_field(
    format: str,
    seed: int,
    *,
    shape: tuple = GRID_DEFAULT_SHAPE,
    n_steps: int = GRID_DEFAULT_STEPS,
    n: int = PHASES_DEFAULT_N,
    n_avalanches: int = AVALANCHES_DEFAULT_N,
    **_: Any,
) -> List[Dict[str, Any]]:
    """3. Random binary field — i.i.d. Bernoulli(0.5)."""
    rng = np.random.default_rng(seed)
    if format == "grid":
        arr = rng.integers(0, 2, size=(n_steps, *shape), dtype=np.int8)
        return _grid_history_from_array(arr)
    if format == "phases":
        # Bernoulli(0.5) on {0, π} — anti-aligned coin flips.
        bits = rng.integers(0, 2, size=(n_steps, n))
        theta_t = bits.astype(float) * np.pi
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        # Avalanche sizes ∈ {1, 2}: trivially not power-law.
        return _wrap_avalanches(rng.integers(1, 3, size=n_avalanches))
    raise ValueError(f"unknown format: {format}")


def spatial_white_noise_series(
    format: str,
    seed: int,
    *,
    shape: tuple = GRID_DEFAULT_SHAPE,
    n_steps: int = GRID_DEFAULT_STEPS,
    n: int = PHASES_DEFAULT_N,
    n_avalanches: int = AVALANCHES_DEFAULT_N,
    **_: Any,
) -> List[Dict[str, Any]]:
    """4. Spatial white noise time series — independent random fields per step.

    Distinguished from random_uniform by intent: this stresses the
    no-temporal-autocorrelation case where every snapshot is a fresh draw.
    """
    rng = np.random.default_rng(seed)
    if format == "grid":
        arr = rng.integers(0, 2, size=(n_steps, *shape), dtype=np.int8)
        return _grid_history_from_array(arr)
    if format == "phases":
        theta_t = rng.uniform(0.0, 2.0 * np.pi, size=(n_steps, n))
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        # Exponential-distributed sizes — the canonical null for SOC tests.
        sizes = (rng.exponential(scale=20.0, size=n_avalanches) + 1).astype(np.int64)
        return _wrap_avalanches(sizes)
    raise ValueError(f"unknown format: {format}")


def temporal_white_noise_per_cell(
    format: str,
    seed: int,
    *,
    shape: tuple = GRID_DEFAULT_SHAPE,
    n_steps: int = GRID_DEFAULT_STEPS,
    n: int = PHASES_DEFAULT_N,
    n_avalanches: int = AVALANCHES_DEFAULT_N,
    **_: Any,
) -> List[Dict[str, Any]]:
    """5. Temporal white noise per cell — each cell evolves as independent walk.

    Every cell flips with i.i.d. probability 0.5 each step, no spatial coupling.
    """
    rng = np.random.default_rng(seed)
    if format == "grid":
        # Random walk in {0, 1}: at each step, each cell has p=0.5 to flip.
        H, W = shape
        flips = rng.integers(0, 2, size=(n_steps, H, W), dtype=np.int8)
        state = rng.integers(0, 2, size=shape, dtype=np.int8)
        out = np.empty((n_steps, H, W), dtype=np.int8)
        for t in range(n_steps):
            state = np.bitwise_xor(state, flips[t])
            out[t] = state
        return _grid_history_from_array(out)
    if format == "phases":
        # Each oscillator's phase performs a small random walk.
        steps = rng.normal(0.0, 0.3, size=(n_steps, n))
        theta = rng.uniform(0.0, 2.0 * np.pi, size=n)
        theta_t = np.empty((n_steps, n))
        for t in range(n_steps):
            theta = np.mod(theta + steps[t], 2.0 * np.pi)
            theta_t[t] = theta
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        # Poisson-distributed sizes (memoryless count process).
        sizes = (rng.poisson(lam=20.0, size=n_avalanches) + 1).astype(np.int64)
        return _wrap_avalanches(sizes)
    raise ValueError(f"unknown format: {format}")


def permutation_shuffled_positive(
    format: str,
    seed: int,
    *,
    positive: List[Dict[str, Any]] | None = None,
    **_: Any,
) -> List[Dict[str, Any]]:
    """6. Permutation-shuffled positive — final state with cells permuted.

    Destroys spatial structure but preserves the marginal value distribution.
    The shuffled state is repeated across the same number of timesteps as the
    original positive so the trajectory shape is preserved.
    """
    if positive is None:
        raise ValueError("permutation_shuffled_positive requires `positive` kwarg")
    rng = np.random.default_rng(seed)
    if format == "grid":
        snap = positive[-1]
        if "grid" in snap:
            last = np.asarray(snap["grid"])
        elif "field" in snap:
            # Continuous-field positive (e.g. Gray-Scott for P3): binarise at
            # median so the permutation generator can produce a valid grid
            # substrate. P3 will reject the resulting binary grid at its
            # substrate prerequisite (no 'field' key), correctly giving TN.
            f = np.asarray(snap["field"], dtype=float)
            last = (f > float(np.median(f))).astype(np.int8)
        else:
            raise KeyError(
                "permutation_shuffled: positive[-1] has neither 'grid' nor 'field'"
            )
        H, W = last.shape
        flat = last.flatten()
        rng.shuffle(flat)
        shuffled = flat.reshape(H, W).astype(np.int8)
        T = len(positive)
        arr = np.broadcast_to(shuffled, (T, H, W)).copy()
        return _grid_history_from_array(arr)
    if format == "phases":
        last_theta = positive[-1]["theta"].copy()
        rng.shuffle(last_theta)
        T = len(positive)
        theta_t = np.broadcast_to(last_theta, (T, last_theta.size)).copy()
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        # Permute the canonical positive's avalanche-size array.
        # NOTE: permutation preserves the marginal distribution → degenerate
        # case (same as oscillator C-class-a-degenerate from Sprint 32).
        sizes = positive[0]["avalanche_sizes"].copy()
        rng.shuffle(sizes)
        return _wrap_avalanches(sizes)
    raise ValueError(f"unknown format: {format}")


def time_shuffled_positive(
    format: str,
    seed: int,
    *,
    positive: List[Dict[str, Any]] | None = None,
    **_: Any,
) -> List[Dict[str, Any]]:
    """7. Time-shuffled positive — full trajectory with timesteps shuffled.

    Preserves spatial structure of each frame but destroys temporal ordering.
    """
    if positive is None:
        raise ValueError("time_shuffled_positive requires `positive` kwarg")
    rng = np.random.default_rng(seed)
    if format == "avalanches":
        # For avalanches the "positive" is a single-element history with a
        # 1-D sizes array; shuffling the array order is the natural analog.
        sizes = positive[0]["avalanche_sizes"].copy()
        rng.shuffle(sizes)
        return _wrap_avalanches(sizes)
    indices = np.arange(len(positive))
    rng.shuffle(indices)
    return [dict(positive[i], step=t) for t, i in enumerate(indices)]


def constant_field(
    format: str,
    seed: int,
    *,
    shape: tuple = GRID_DEFAULT_SHAPE,
    n_steps: int = GRID_DEFAULT_STEPS,
    n: int = PHASES_DEFAULT_N,
    n_avalanches: int = AVALANCHES_DEFAULT_N,
    value: int = 1,
    **_: Any,
) -> List[Dict[str, Any]]:
    """8. Constant field — all cells fixed to one value.

    Catches detectors that fire on triviality.
    """
    if format == "grid":
        arr = np.full((n_steps, *shape), value, dtype=np.int8)
        return _grid_history_from_array(arr)
    if format == "phases":
        theta_t = np.zeros((n_steps, n))
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        return _wrap_avalanches(np.full(n_avalanches, max(value, 1), dtype=np.int64))
    raise ValueError(f"unknown format: {format}")


def linear_gradient_field(
    format: str,
    seed: int,
    *,
    shape: tuple = GRID_DEFAULT_SHAPE,
    n_steps: int = GRID_DEFAULT_STEPS,
    n: int = PHASES_DEFAULT_N,
    n_avalanches: int = AVALANCHES_DEFAULT_N,
    **_: Any,
) -> List[Dict[str, Any]]:
    """9. Linear gradient — smooth monotonic spatial gradient (no emergence)."""
    if format == "grid":
        H, W = shape
        col = np.arange(W) / max(1, W - 1)
        gradient = np.broadcast_to(col, (H, W))
        arr_static = (gradient > 0.5).astype(np.int8)
        arr = np.broadcast_to(arr_static, (n_steps, H, W)).copy()
        return _grid_history_from_array(arr)
    if format == "phases":
        theta_static = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
        theta_t = np.broadcast_to(theta_static, (n_steps, n)).copy()
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        # Arithmetic progression — strictly monotonic, not power-law.
        return _wrap_avalanches(np.arange(1, n_avalanches + 1, dtype=np.int64))
    raise ValueError(f"unknown format: {format}")


def periodic_checkerboard(
    format: str,
    seed: int,
    *,
    shape: tuple = GRID_DEFAULT_SHAPE,
    n_steps: int = GRID_DEFAULT_STEPS,
    n: int = PHASES_DEFAULT_N,
    n_avalanches: int = AVALANCHES_DEFAULT_N,
    **_: Any,
) -> List[Dict[str, Any]]:
    """10. Periodic boundary checkerboard — alternating values, deterministic."""
    if format == "grid":
        H, W = shape
        rows = np.arange(H)[:, None]
        cols = np.arange(W)[None, :]
        cb = ((rows + cols) % 2).astype(np.int8)
        arr = np.broadcast_to(cb, (n_steps, H, W)).copy()
        return _grid_history_from_array(arr)
    if format == "phases":
        idx = np.arange(n)
        theta_static = (idx % 2).astype(float) * np.pi
        theta_t = np.broadcast_to(theta_static, (n_steps, n)).copy()
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        # Alternating two values — bimodal, definitely not power-law.
        sizes = np.tile([5, 10], n_avalanches // 2 + 1)[:n_avalanches].astype(np.int64)
        return _wrap_avalanches(sizes)
    raise ValueError(f"unknown format: {format}")


SYNTHETIC_GENERATORS: Dict[str, Callable[..., List[Dict[str, Any]]]] = {
    "random_uniform": random_uniform_field,
    "random_gaussian": random_gaussian_field,
    "random_binary": random_binary_field,
    "spatial_white_noise": spatial_white_noise_series,
    "temporal_white_noise": temporal_white_noise_per_cell,
    "permutation_shuffled": permutation_shuffled_positive,
    "time_shuffled": time_shuffled_positive,
    "constant": constant_field,
    "linear_gradient": linear_gradient_field,
    "checkerboard": periodic_checkerboard,
}
