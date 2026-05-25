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


# --- Lattice_1d supplements -------------------------------------------------

def independent_lane_traffic(
    seed: int,
    *,
    road_length: int = 100,
    n_cars: int = 30,
    n_steps: int = 200,
    **_: Any,
) -> List[Dict[str, Any]]:
    """N cars on a 1-D ring, each advancing 1 cell per step — no gap lookahead.

    Each car moves forward unconditionally (no slowing-down rule, no
    randomization). There is no jam formation: stopped_fraction = 0 always.
    A correctly-specific traffic-jam detector should not fire.

    Returns a list of step snapshots with ``array`` (binary occupancy of
    length ``road_length``) and ``step`` keys, matching the lattice_1d
    sequence shape contract.
    """
    rng = np.random.default_rng(seed)
    positions = rng.choice(road_length, size=n_cars, replace=False).astype(np.int32)
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        positions = (positions + 1) % road_length
        occupancy = np.zeros(road_length, dtype=np.int8)
        occupancy[positions] = 1
        out.append({"array": occupancy.copy(), "step": t})
    return out


def reverse_sorted_sequence(
    seed: int,
    *,
    n: int = 64,
    n_steps: int = 200,
    noise_scale: float = 0.5,
    **_: Any,
) -> List[Dict[str, Any]]:
    """Sorted-descending integer sequence with small per-step perturbations.

    Initial state: ``[n, n-1, ..., 1]`` (fully reverse-sorted). Each step
    adds i.i.d. integer noise in ``[-noise_scale, +noise_scale]`` (rounded)
    — no swap dynamics. Tests that detectors for emergent ordering don't fire
    on an already-sorted (imposed) sequence.

    Returns a list of step snapshots with ``array`` (int32 values of length
    ``n``) and ``step`` keys, matching the lattice_1d sequence shape contract.
    """
    rng = np.random.default_rng(seed)
    base = np.arange(n, 0, -1, dtype=np.float32)  # [n, n-1, ..., 1]
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        noise = rng.uniform(-noise_scale, noise_scale, size=n).astype(np.float32)
        arr = np.round(base + noise).astype(np.int32)
        out.append({"array": arr.copy(), "step": t})
    return out


# --- Lattice_2d_continuous supplements -------------------------------------
# Sprint 43: added to give P3 (Gray-Scott) a non-empty Class B pool.
# P3's substrate type is lattice_2d_continuous; Gray-Scott is the only
# registered catalog model of that type (0 Class B catalog mates).
# These supplements produce field_continuous histories — dicts with a
# 'field' key (2D float32 array) — matching the P3 detector's native
# consumption format so no cross-format adapter is needed.


def smooth_random_field(
    seed: int,
    *,
    rows: int = 64,
    cols: int = 64,
    n_steps: int = 100,
    **_: Any,
) -> List[Dict[str, Any]]:
    """Pure i.i.d. white Gaussian noise normalised to [0, 1].

    Each step is an independent draw of Gaussian white noise normalised to
    [0, 1].  The radial FFT spectrum is flat (no sharp peak): peak_to_mean
    ≈ 1.25 < 5.0 (P3 screening threshold).  P3 should reject at screening.

    NOTE: a spatially-correlated (Gaussian-filtered) variant was prototyped but
    produces a spectral peak at the filter's correlation length that falsely
    passes P3's screening threshold.  Pure white noise is the correct null.

    Produces kind='field_continuous' output matching ``_gen_p3_gray_scott``'s
    shape contract: each step is ``{"field": np.ndarray(rows, cols), "step": t}``.
    """
    rng = np.random.default_rng(seed)
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        noise = rng.standard_normal((rows, cols)).astype(np.float32)
        lo, hi = float(noise.min()), float(noise.max())
        denom = max(hi - lo, 1e-9)
        field = ((noise - lo) / denom).astype(np.float32)
        out.append({"field": field, "step": t})
    return out


def sinusoidal_traveling_wave(
    seed: int,
    *,
    rows: int = 64,
    cols: int = 64,
    n_steps: int = 100,
    period: float = 50.0,
    **_: Any,
) -> List[Dict[str, Any]]:
    """Spatially-uniform temporal oscillation — every cell has the same value.

    ``field[t, i, j] = 0.5 + 0.5 * sin(2π * t / period)``

    The field is spatially homogeneous at every snapshot: field_std = 0.
    P3's prerequisite check (field_std >= 0.01) rejects the substrate before
    the Turing-wavelength pipeline is reached.

    NOTE: a spatially-varying traveling wave (with explicit wavelength) was
    prototyped first but produced a sharp FFT peak that falsely passed P3's
    screening threshold.  Spatially-uniform oscillation is the correct null.

    The ``seed`` parameter is accepted but unused (deterministic by construction).
    """
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        value = float(0.5 + 0.5 * np.sin(2.0 * np.pi * t / period))
        field = np.full((rows, cols), value, dtype=np.float32)
        out.append({"field": field, "step": t})
    return out


# --- Registry ---------------------------------------------------------------

SUPPLEMENTS_BY_SUBSTRATE_TYPE: Dict[str, List[str]] = {
    "oscillator": ["incoherent_phases", "subcritical_kuramoto"],
    "network": ["random_graph_evolution", "network_random_walks"],
    "continuous_2d": ["uncorrelated_random_walks", "independent_brownian_motion"],
    "lattice_1d": ["independent_lane_traffic", "reverse_sorted_sequence"],
    "lattice_2d_continuous": ["smooth_random_field", "sinusoidal_traveling_wave"],
}

SUPPLEMENT_BUILDERS: Dict[str, Callable[..., List[Dict[str, Any]]]] = {
    "incoherent_phases": incoherent_phases,
    "subcritical_kuramoto": subcritical_kuramoto,
    "random_graph_evolution": random_graph_evolution,
    "network_random_walks": network_random_walks,
    "uncorrelated_random_walks": uncorrelated_random_walks,
    "independent_brownian_motion": independent_brownian_motion,
    "independent_lane_traffic": independent_lane_traffic,
    "reverse_sorted_sequence": reverse_sorted_sequence,
    "smooth_random_field": smooth_random_field,
    "sinusoidal_traveling_wave": sinusoidal_traveling_wave,
}
