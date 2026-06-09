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
PARTICLES_DEFAULT_N = 200
PARTICLES_DEFAULT_BOX = 16.0
PARTICLES_DEFAULT_SPEED = 0.03


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


SEQUENCE_DEFAULT_N = 200  # default lattice_1d road/array length
OPINIONS_DEFAULT_N = 100  # default number of agents for opinions format

# Scalar timeseries defaults (P24 homeostasis panel)
SCALAR_TS_DT = 0.1                   # integration timestep
SCALAR_TS_SETPOINT = 10.0            # regulated-variable setpoint
SCALAR_TS_PERT_AMPLITUDE = 5.0       # perturbation amplitude
SCALAR_TS_PERT_ONSET_FRAC = 0.25     # perturbation onset at 25% of trajectory

# Choice timeseries defaults (P23 anti-coordination panel)
CHOICE_TS_N_AGENTS = 101             # default agent count (odd for MG symmetry)
CHOICE_TS_DEFAULT_STEPS = 500        # default number of rounds

# Noise-sweep timeseries defaults (P26 stochastic resonance panel)
NOISE_SWEEP_N_LEVELS = 10            # number of noise levels in sweep
NOISE_SWEEP_STEPS_PER_LEVEL = 2000   # timesteps per noise level
NOISE_SWEEP_SIGNAL_FREQ = 0.005      # signal frequency (Hz)
NOISE_SWEEP_SIGNAL_AMP = 1.0         # signal amplitude
NOISE_SWEEP_DT = 0.01                # integration timestep

# State-vector (attractor-network) defaults (P16 associative memory panel)
STATE_VECTOR_DEFAULT_N = 100         # default number of neurons
STATE_VECTOR_DEFAULT_P = 5           # default number of stored patterns
STATE_VECTOR_DEFAULT_TRIALS = 5      # default number of retrieval trials
STATE_VECTOR_DEFAULT_STEPS_PER_TRIAL = 10  # default steps per trial

# Canalization bundle defaults (P25 equifinality panel)
CANAL_DEFAULT_N_ICS = 20             # default number of initial conditions
CANAL_DEFAULT_N_DIMS = 10            # default state-space dimensionality
CANAL_DEFAULT_N_STEPS = 200          # default steps per trajectory


def _noise_sweep_history(
    x_per_level: List[np.ndarray],
    noise_levels: np.ndarray,
    dt: float = NOISE_SWEEP_DT,
    signal_freq: float = NOISE_SWEEP_SIGNAL_FREQ,
    signal_amp: float = NOISE_SWEEP_SIGNAL_AMP,
) -> List[Dict[str, Any]]:
    """Wrap per-level x arrays into noise-sweep-timeseries history for P26.

    Parameters
    ----------
    x_per_level : list of 1-D arrays, one per noise level
    noise_levels : 1-D array of noise level values
    """
    out: List[Dict[str, Any]] = []
    for nl_idx, (x_arr, D) in enumerate(zip(x_per_level, noise_levels)):
        for step_i, x_val in enumerate(x_arr):
            t = step_i * dt
            signal = signal_amp * np.sin(2.0 * np.pi * signal_freq * t)
            out.append({
                'time': t,
                'x': float(x_val),
                'signal': float(signal),
                'noise_level': float(D),
                'noise_level_idx': nl_idx,
                'step': step_i,
            })
    return out


def _default_null_noise_levels(n_levels: int = NOISE_SWEEP_N_LEVELS) -> np.ndarray:
    """Default noise levels for null noise-sweep substrates."""
    return np.linspace(0.0, 10.0, n_levels)


def _noise_sweep_pure_noise(
    rng: np.random.Generator,
    n_levels: int = NOISE_SWEEP_N_LEVELS,
    n_steps: int = NOISE_SWEEP_STEPS_PER_LEVEL,
    scale_func: str = "linear",
) -> List[Dict[str, Any]]:
    """Generate noise-sweep substrate where x = noise only (no signal correlation).

    At all noise levels, x is i.i.d. noise → coherent response ≈ 0 everywhere
    → no inverted-U → P26 rejects at screening.
    """
    noise_levels = _default_null_noise_levels(n_levels)
    x_per_level = []
    for D in noise_levels:
        if scale_func == "linear":
            x = rng.normal(0.0, max(D, 0.01), size=n_steps)
        elif scale_func == "uniform":
            x = rng.uniform(-max(D, 0.01), max(D, 0.01), size=n_steps)
        else:
            x = rng.normal(0.0, 1.0, size=n_steps)
        x_per_level.append(x)
    return _noise_sweep_history(x_per_level, noise_levels)


def _noise_sweep_monotone(
    rng: np.random.Generator,
    n_levels: int = NOISE_SWEEP_N_LEVELS,
    n_steps: int = NOISE_SWEEP_STEPS_PER_LEVEL,
    direction: str = "increasing",
) -> List[Dict[str, Any]]:
    """Generate noise-sweep substrate with monotone performance vs noise.

    x = amplitude * signal + noise, where amplitude increases (or decreases)
    monotonically with noise level → no interior peak → P26 rejects.
    """
    noise_levels = _default_null_noise_levels(n_levels)
    x_per_level = []
    dt = NOISE_SWEEP_DT
    for i, D in enumerate(noise_levels):
        if direction == "increasing":
            amp = (i + 1) / n_levels * 5.0
        else:
            amp = (n_levels - i) / n_levels * 5.0
        x = np.zeros(n_steps)
        for step_i in range(n_steps):
            t = step_i * dt
            signal = NOISE_SWEEP_SIGNAL_AMP * np.sin(
                2.0 * np.pi * NOISE_SWEEP_SIGNAL_FREQ * t
            )
            x[step_i] = amp * signal + rng.normal(0.0, max(D, 0.01))
        x_per_level.append(x)
    return _noise_sweep_history(x_per_level, noise_levels)


def _choice_ts_history(
    attendance: np.ndarray,
    n_agents: int = CHOICE_TS_N_AGENTS,
) -> List[Dict[str, Any]]:
    """Wrap an attendance array into choice_timeseries history format."""
    return [
        {'round': int(t), 'attendance': int(a), 'n_agents': n_agents}
        for t, a in enumerate(attendance)
    ]


def _state_vector_null_history(
    rng: np.random.Generator,
    N: int = STATE_VECTOR_DEFAULT_N,
    P: int = STATE_VECTOR_DEFAULT_P,
    n_trials: int = STATE_VECTOR_DEFAULT_TRIALS,
    n_steps_per_trial: int = STATE_VECTOR_DEFAULT_STEPS_PER_TRIAL,
) -> List[Dict[str, Any]]:
    """Null state-vector history: random ±1 states with no associative memory.

    Generates a format-compatible attractor-network history where states are
    i.i.d. random (no convergence to stored patterns). The P16 detector should
    reject: completion accuracy ≈ chance overlap √(1/N).
    """
    stored_patterns = rng.choice([-1, 1], size=(P, N)).astype(np.int8)
    history: List[Dict[str, Any]] = []
    for trial in range(n_trials):
        cue_idx = trial % P
        for step in range(n_steps_per_trial):
            state = rng.choice([-1, 1], size=N).astype(np.int8)
            overlap = float(
                np.dot(stored_patterns[cue_idx].astype(np.int32),
                       state.astype(np.int32)) / N
            )
            history.append({
                'state': state,
                'step': step,
                'trial': trial,
                'cue_pattern_idx': cue_idx,
                'overlap': overlap,
                'stored_patterns': stored_patterns.copy(),
                'converged': step == n_steps_per_trial - 1,
            })
    return history


def _canalization_bundle_null(
    rng: np.random.Generator,
    n_ics: int = CANAL_DEFAULT_N_ICS,
    n_dims: int = CANAL_DEFAULT_N_DIMS,
    n_steps: int = CANAL_DEFAULT_N_STEPS,
    dynamics: str = "random_walk",
) -> List[Dict[str, Any]]:
    """Null canalization bundle: non-converging multi-IC trajectories.

    Generates observation-bundle-format histories where ICs do NOT converge
    to a common target. P25 should reject at screening (convergence variance
    ratio ≥ 0.1) or at observation-bundle extraction (missing keys → no detect).

    Parameters
    ----------
    dynamics : str
        "random_walk" — pure diffusion (default).
        "constant" — ICs stay in place (no dynamics, spread preserved).
        "divergent" — ICs actively repel from target.
    """
    target = rng.standard_normal(n_dims)
    history: List[Dict[str, Any]] = []
    for trial in range(n_ics):
        x = rng.standard_normal(n_dims) * 5.0
        ic = x.copy()
        for step in range(n_steps):
            dist = float(np.linalg.norm(x - target))
            history.append({
                'state': x.copy(),
                'step': step,
                'trial': trial,
                'ic': ic.copy(),
                'target': target.copy(),
                'distance_to_target': dist,
                'converged': dist < 0.1,
            })
            if dynamics == "random_walk":
                x = x + rng.standard_normal(n_dims) * 0.1
            elif dynamics == "divergent":
                delta = x - target
                norm = np.linalg.norm(delta) + 1e-12
                x = x + 0.05 * delta / norm + rng.standard_normal(n_dims) * 0.05
            # "constant" → x unchanged
    return history


def _opinions_history_from_array(ops: np.ndarray) -> List[Dict[str, Any]]:
    """Wrap a (T, N) float array into the standard opinions-history format.

    Each step dict has an ``opinions`` key (1-D float array ∈ [0, 1] of
    length N).  These are the null substrates for P21's dip-test-based
    polarization detector — they should be unimodal so the detector rejects.
    """
    T, N = ops.shape
    return [{"opinions": ops[t].astype(np.float64), "step": t} for t in range(T)]


def _scalar_ts_history(
    x_arr: np.ndarray,
    n_steps: int,
    dt: float = SCALAR_TS_DT,
    setpoint: float = SCALAR_TS_SETPOINT,
    pert_amplitude: float = SCALAR_TS_PERT_AMPLITUDE,
    onset_frac: float = SCALAR_TS_PERT_ONSET_FRAC,
) -> List[Dict[str, Any]]:
    """Wrap a 1-D x array into scalar-timeseries history format for P24.

    Perturbation is a sustained step function starting at ``onset_frac`` of
    the trajectory, matching the ProportionalHomeostat convention.
    """
    onset_step = int(n_steps * onset_frac)
    out: List[Dict[str, Any]] = []
    for t in range(min(n_steps, len(x_arr))):
        pert = pert_amplitude if t >= onset_step else 0.0
        out.append({
            "time": t * dt,
            "x": float(x_arr[t]),
            "setpoint": setpoint,
            "perturbation": pert,
            "deviation": float(x_arr[t] - setpoint),
            "step": t,
        })
    return out


def _scalar_ts_uncontrolled(
    rng: np.random.Generator,
    n_steps: int,
    noise_std: float = 0.5,
    dt: float = SCALAR_TS_DT,
    setpoint: float = SCALAR_TS_SETPOINT,
    pert_amplitude: float = SCALAR_TS_PERT_AMPLITUDE,
    onset_frac: float = SCALAR_TS_PERT_ONSET_FRAC,
) -> np.ndarray:
    """Generate an uncontrolled scalar trajectory: x drifts under perturbation.

    dx = perturbation(t) * dt + noise — no feedback correction.
    Deviation from setpoint grows linearly post-onset → growth_ratio >> 2.0.
    """
    onset_step = int(n_steps * onset_frac)
    x = np.full(n_steps, setpoint)
    for t in range(1, n_steps):
        pert = pert_amplitude if t >= onset_step else 0.0
        noise = rng.normal(0, noise_std) * np.sqrt(dt)
        x[t] = x[t - 1] + pert * dt + noise
    return x


def _sequence_history_from_array(arrays: np.ndarray) -> List[Dict[str, Any]]:
    """Wrap a (T, N) integer array into the standard sequence-history format.

    Each step dict has an ``array`` key (1-D integer array of length N).
    These histories lack a ``velocities`` key so the P8 traffic-jam detector
    rejects them at prerequisites — the intended behaviour for null substrates.
    """
    T, N = arrays.shape
    return [{"array": arrays[t].copy(), "step": t} for t in range(T)]


def _particles_history_from_random(
    rng: np.random.Generator,
    n: int,
    n_steps: int,
    box_size: float = PARTICLES_DEFAULT_BOX,
    speed: float = PARTICLES_DEFAULT_SPEED,
) -> List[Dict[str, Any]]:
    """Null particle substrate: random positions + random headings, no interactions."""
    positions = rng.uniform(0.0, box_size, (n, 2))
    out: List[Dict[str, Any]] = []
    for t in range(n_steps):
        headings = rng.uniform(-np.pi, np.pi, n)
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
    if format == "particles":
        return _particles_history_from_random(rng, n, n_steps)
    if format == "sequence":
        arr = rng.integers(0, 2, size=(n_steps, n), dtype=np.int8)
        return _sequence_history_from_array(arr)
    if format == "opinions":
        # Uniform opinions ∈ [0,1] — unimodal, P21 dip-test should reject.
        ops = rng.uniform(0.0, 1.0, size=(n_steps, OPINIONS_DEFAULT_N))
        return _opinions_history_from_array(ops)
    if format == "scalar_timeseries":
        # Uncontrolled drift under perturbation + uniform noise.
        x = _scalar_ts_uncontrolled(rng, n_steps, noise_std=0.5)
        return _scalar_ts_history(x, n_steps)
    if format == "noise_sweep":
        return _noise_sweep_pure_noise(rng, scale_func="uniform")
    if format == "choice_timeseries":
        # i.i.d. Binomial(N, 0.5) attendance — random-choice null.
        att = rng.binomial(CHOICE_TS_N_AGENTS, 0.5, size=n_steps)
        return _choice_ts_history(att)
    if format == "state_vector":
        return _state_vector_null_history(rng)
    if format == "canalization_bundle":
        return _canalization_bundle_null(rng, dynamics="random_walk")
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
    if format == "particles":
        return _particles_history_from_random(rng, n, n_steps)
    if format == "sequence":
        z = rng.standard_normal(size=(n_steps, n))
        arr = (z > 0.0).astype(np.int8)
        return _sequence_history_from_array(arr)
    if format == "opinions":
        # Gaussian opinions clipped to [0,1] — unimodal centred at 0.5.
        z = rng.standard_normal(size=(n_steps, OPINIONS_DEFAULT_N))
        ops = np.clip(z * 0.2 + 0.5, 0.0, 1.0)
        return _opinions_history_from_array(ops)
    if format == "scalar_timeseries":
        # Uncontrolled drift under perturbation + Gaussian noise.
        x = _scalar_ts_uncontrolled(rng, n_steps, noise_std=1.0)
        return _scalar_ts_history(x, n_steps)
    if format == "noise_sweep":
        return _noise_sweep_pure_noise(rng, scale_func="linear")
    if format == "choice_timeseries":
        # Gaussian-rounded Binomial attendance — random-choice null with heavier tails.
        att = np.clip(rng.normal(CHOICE_TS_N_AGENTS / 2, np.sqrt(CHOICE_TS_N_AGENTS) / 2,
                                  size=n_steps).round().astype(int), 0, CHOICE_TS_N_AGENTS)
        return _choice_ts_history(att)
    if format == "state_vector":
        return _state_vector_null_history(rng)
    if format == "canalization_bundle":
        return _canalization_bundle_null(rng, dynamics="random_walk")
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
    if format == "particles":
        return _particles_history_from_random(rng, n, n_steps)
    if format == "sequence":
        arr = rng.integers(0, 2, size=(n_steps, n), dtype=np.int8)
        return _sequence_history_from_array(arr)
    if format == "opinions":
        # Continuous uniform opinions — NOT binary {0,1} which would be bimodal.
        # "Binary" here means the grid representation; opinions are continuous [0,1].
        ops = rng.uniform(0.0, 1.0, size=(n_steps, OPINIONS_DEFAULT_N))
        return _opinions_history_from_array(ops)
    if format == "scalar_timeseries":
        # Uncontrolled drift with small noise — binary-like steps.
        x = _scalar_ts_uncontrolled(rng, n_steps, noise_std=0.3)
        return _scalar_ts_history(x, n_steps)
    if format == "noise_sweep":
        # Binary: x ∈ {-1, +1} i.i.d. at all noise levels → no signal correlation.
        noise_levels = _default_null_noise_levels()
        x_per_level = [rng.choice([-1.0, 1.0], size=NOISE_SWEEP_STEPS_PER_LEVEL)
                        for _ in noise_levels]
        return _noise_sweep_history(x_per_level, noise_levels)
    if format == "choice_timeseries":
        # Bernoulli(0.5) binary attendance — random-choice null.
        att = rng.binomial(CHOICE_TS_N_AGENTS, 0.5, size=n_steps)
        return _choice_ts_history(att)
    if format == "state_vector":
        return _state_vector_null_history(rng)
    if format == "canalization_bundle":
        return _canalization_bundle_null(rng, dynamics="random_walk")
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
    if format == "particles":
        return _particles_history_from_random(rng, n, n_steps)
    if format == "sequence":
        arr = rng.integers(0, 2, size=(n_steps, n), dtype=np.int8)
        return _sequence_history_from_array(arr)
    if format == "opinions":
        # Fresh i.i.d. uniform opinions each step — no temporal autocorrelation.
        ops = rng.uniform(0.0, 1.0, size=(n_steps, OPINIONS_DEFAULT_N))
        return _opinions_history_from_array(ops)
    if format == "scalar_timeseries":
        # i.i.d. uniform x each step (no temporal autocorrelation) + perturbation.
        x = _scalar_ts_uncontrolled(rng, n_steps, noise_std=2.0)
        return _scalar_ts_history(x, n_steps)
    if format == "noise_sweep":
        # Fresh i.i.d. Gaussian x each step, same variance at all noise levels.
        noise_levels = _default_null_noise_levels()
        x_per_level = [rng.normal(0.0, 1.0, size=NOISE_SWEEP_STEPS_PER_LEVEL)
                        for _ in noise_levels]
        return _noise_sweep_history(x_per_level, noise_levels)
    if format == "choice_timeseries":
        # i.i.d. Binomial attendance redrawn each step — no temporal autocorrelation.
        att = rng.binomial(CHOICE_TS_N_AGENTS, 0.5, size=n_steps)
        return _choice_ts_history(att)
    if format == "state_vector":
        return _state_vector_null_history(rng)
    if format == "canalization_bundle":
        return _canalization_bundle_null(rng, dynamics="random_walk")
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
    if format == "particles":
        return _particles_history_from_random(rng, n, n_steps)
    if format == "sequence":
        flips = rng.integers(0, 2, size=(n_steps, n), dtype=np.int8)
        state = rng.integers(0, 2, size=n, dtype=np.int8)
        out = np.empty((n_steps, n), dtype=np.int8)
        for t in range(n_steps):
            state = np.bitwise_xor(state, flips[t])
            out[t] = state
        return _sequence_history_from_array(out)
    if format == "opinions":
        # Each agent's opinion performs an independent random walk clamped to [0,1].
        ops = np.empty((n_steps, OPINIONS_DEFAULT_N))
        state_op = rng.uniform(0.0, 1.0, size=OPINIONS_DEFAULT_N)
        for t in range(n_steps):
            state_op = np.clip(state_op + rng.normal(0.0, 0.05, OPINIONS_DEFAULT_N), 0.0, 1.0)
            ops[t] = state_op
        return _opinions_history_from_array(ops)
    if format == "scalar_timeseries":
        # Random walk under perturbation — noise-dominated drift.
        x = _scalar_ts_uncontrolled(rng, n_steps, noise_std=0.8)
        return _scalar_ts_history(x, n_steps)
    if format == "noise_sweep":
        # Random walk x at each noise level — noise-dominated, no signal correlation.
        noise_levels = _default_null_noise_levels()
        x_per_level = []
        for D in noise_levels:
            x = np.cumsum(rng.normal(0.0, max(D, 0.01) * 0.01,
                                      size=NOISE_SWEEP_STEPS_PER_LEVEL))
            x_per_level.append(x)
        return _noise_sweep_history(x_per_level, noise_levels)
    if format == "choice_timeseries":
        # Random-walk attendance with independent per-step noise.
        att = np.clip(
            np.cumsum(rng.choice([-1, 0, 1], size=n_steps)) + CHOICE_TS_N_AGENTS // 2,
            0, CHOICE_TS_N_AGENTS,
        ).astype(int)
        return _choice_ts_history(att)
    if format == "state_vector":
        return _state_vector_null_history(rng)
    if format == "canalization_bundle":
        return _canalization_bundle_null(rng, dynamics="divergent")
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
    if format == "particles":
        # Particle permutation: shuffle headings in the final frame, keep positions.
        # Destroys velocity alignment while preserving the spatial configuration.
        # NOTE: for density-based detectors (P2 MIPS), this is degenerate by
        # construction — local density is unchanged — carry-forward C-p2-perm-shuffled-fp.
        snap = positive[-1]
        positions = np.asarray(snap["positions"], dtype=np.float64).copy()
        speeds_arr = np.asarray(snap.get("speeds", np.full(len(positions), PARTICLES_DEFAULT_SPEED)))
        headings_last = np.asarray(snap["headings"], dtype=np.float64).copy()
        rng.shuffle(headings_last)
        velocities = np.column_stack([
            speeds_arr * np.cos(headings_last),
            speeds_arr * np.sin(headings_last),
        ])
        T = len(positive)
        return [
            {
                "positions": positions.copy(),
                "headings": headings_last.copy(),
                "velocities": velocities.copy(),
                "speeds": speeds_arr.copy(),
                "step": t,
            }
            for t in range(T)
        ]
    if format == "sequence":
        # Sequence permutation: shuffle velocities (or array) from the final frame.
        # For NS-format positives (velocities key present), this preserves
        # stopped_fraction (permutation-invariant) so P8 may still detect at
        # SCREENING tier — expected FP per brief (permutation_invariant flag
        # is missing from detector_invariance.py; carry-forward C-p8-perm-shuffled-fp).
        snap = positive[-1]
        if "velocities" in snap:
            vels = np.asarray(snap["velocities"]).copy()
            rng.shuffle(vels)
            T = len(positive)
            return [dict(positive[t], velocities=vels.copy(), step=t) for t in range(T)]
        # Fallback: shuffle the 'array' key if present, else produce zeros.
        arr_last = np.asarray(snap.get("array", np.zeros(SEQUENCE_DEFAULT_N, dtype=np.int8))).copy()
        rng.shuffle(arr_last)
        T = len(positive)
        return [{"array": arr_last.copy(), "step": t} for t in range(T)]
    if format == "opinions":
        # Permute the final state's opinion values (preserves distribution).
        # For HK bimodal positives this IS degenerate: shuffled opinions keep
        # the same bimodal distribution → dip test fires → expected FP.
        # Carry-forward: C-p21-perm-shuffled-fp.
        snap = positive[-1]
        last_ops = np.asarray(snap.get("opinions", np.zeros(OPINIONS_DEFAULT_N))).copy()
        rng.shuffle(last_ops)
        T = len(positive)
        return [{"opinions": last_ops.copy(), "step": t} for t in range(T)]
    if format == "scalar_timeseries":
        # Scalar timeseries has a single variable — permutation is a no-op.
        # Return the positive unchanged (degenerate-by-construction, will be
        # skipped by invariance flag permutation_invariant=True).
        return list(positive)
    if format == "noise_sweep":
        # Noise-sweep has a single scalar x — permutation is a no-op.
        # Degenerate-by-construction (will be skipped by permutation_invariant=True).
        return list(positive)
    if format == "choice_timeseries":
        # Permute agent indices — for attendance-only time series this is a no-op
        # (σ²/N is permutation-invariant). Degenerate-by-construction; will be
        # skipped by permutation_invariant=True.
        return list(positive)
    if format == "state_vector":
        # Consistent permutation of neuron indices across state AND stored_patterns.
        # Overlap = (1/N)|Σ ξ_{π(i)} s_{π(i)}| = (1/N)|Σ ξ_i s_i| — preserved.
        # Degenerate-by-construction; will be skipped by permutation_invariant=True.
        perm = rng.permutation(len(positive[0]['state']))
        out: List[Dict[str, Any]] = []
        for h in positive:
            new_h = dict(h)
            new_h['state'] = np.asarray(h['state'])[perm].copy()
            new_h['stored_patterns'] = np.asarray(h['stored_patterns'])[:, perm].copy()
            out.append(new_h)
        return out
    if format == "canalization_bundle":
        # Permute trial indices in the positive — convergence variance ratio
        # is invariant under IC relabelling → degenerate-by-construction.
        # Will be skipped by permutation_invariant=True.
        return list(positive)
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
    if format == "scalar_timeseries":
        # Shuffle timesteps of the scalar timeseries positive. Preserves
        # deviation distribution but destroys temporal recovery structure.
        indices = np.arange(len(positive))
        rng.shuffle(indices)
        return [dict(positive[i], step=t, time=t * SCALAR_TS_DT)
                for t, i in enumerate(indices)]
    if format == "noise_sweep":
        # Shuffle timesteps across ALL noise levels. Each dict retains its
        # original noise_level/noise_level_idx, so the detector still groups
        # correctly — but within each level, x and signal are decorrelated
        # because they come from random timesteps → coherent response ≈ 0.
        indices = np.arange(len(positive))
        rng.shuffle(indices)
        return [dict(positive[i], step=t, time=t * NOISE_SWEEP_DT)
                for t, i in enumerate(indices)]
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
    rng = np.random.default_rng(seed)
    if format == "grid":
        arr = np.full((n_steps, *shape), value, dtype=np.int8)
        return _grid_history_from_array(arr)
    if format == "phases":
        theta_t = np.zeros((n_steps, n))
        return _phases_history_from_array(theta_t)
    if format == "avalanches":
        return _wrap_avalanches(np.full(n_avalanches, max(value, 1), dtype=np.int64))
    if format == "particles":
        # Constant field: all particles move in the same direction (θ=0).
        # Produces φ=1.0 — a known degenerate substrate for flocking detectors
        # (analogous to the P9 constant_field trivial-sync carry-forward).
        positions = rng.uniform(0.0, PARTICLES_DEFAULT_BOX, (n, 2))
        headings = np.zeros(n, dtype=np.float64)
        velocities = PARTICLES_DEFAULT_SPEED * np.column_stack([
            np.cos(headings), np.sin(headings)
        ])
        out = []
        for t in range(n_steps):
            positions = (positions + velocities) % PARTICLES_DEFAULT_BOX
            out.append({
                "positions": positions.copy(),
                "headings": headings.copy(),
                "velocities": velocities.copy(),
                "speeds": np.full(n, PARTICLES_DEFAULT_SPEED, dtype=np.float64),
                "step": t,
            })
        return out
    if format == "sequence":
        arr = np.full((n_steps, n), value, dtype=np.int8)
        return _sequence_history_from_array(arr)
    if format == "opinions":
        # All agents share a single opinion at 0.5 — trivially unimodal.
        ops = np.full((n_steps, OPINIONS_DEFAULT_N), 0.5, dtype=np.float64)
        return _opinions_history_from_array(ops)
    if format == "scalar_timeseries":
        # x = setpoint for all time, perturbation = 0 → no perturbation onset.
        # Detector rejects at "No perturbation detected" prerequisite.
        x = np.full(n_steps, SCALAR_TS_SETPOINT)
        return _scalar_ts_history(x, n_steps, pert_amplitude=0.0)
    if format == "noise_sweep":
        # Constant x = 0 at all noise levels → coherent response = 0 → screening fails.
        noise_levels = _default_null_noise_levels()
        x_per_level = [np.zeros(NOISE_SWEEP_STEPS_PER_LEVEL) for _ in noise_levels]
        return _noise_sweep_history(x_per_level, noise_levels)
    if format == "choice_timeseries":
        # All agents choose the same side → attendance = N/2 every round.
        # Variance = 0, mean_ratio = 1.0 — screening passes but confirmation
        # fails (variance at zero is below baseline, but trivially so).
        att = np.full(n_steps, CHOICE_TS_N_AGENTS // 2, dtype=int)
        return _choice_ts_history(att)
    if format == "state_vector":
        # All neurons fixed to +1 — trivially not content-addressable.
        N = STATE_VECTOR_DEFAULT_N
        P = STATE_VECTOR_DEFAULT_P
        stored_patterns = rng.choice([-1, 1], size=(P, N)).astype(np.int8)
        state = np.ones(N, dtype=np.int8)
        history: List[Dict[str, Any]] = []
        for trial in range(STATE_VECTOR_DEFAULT_TRIALS):
            for step in range(STATE_VECTOR_DEFAULT_STEPS_PER_TRIAL):
                overlap = float(np.dot(stored_patterns[trial % P].astype(np.int32),
                                       state.astype(np.int32)) / N)
                history.append({
                    'state': state.copy(), 'step': step, 'trial': trial,
                    'cue_pattern_idx': trial % P, 'overlap': overlap,
                    'stored_patterns': stored_patterns.copy(), 'converged': True,
                })
        return history
    if format == "canalization_bundle":
        return _canalization_bundle_null(rng, dynamics="constant")
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
    rng = np.random.default_rng(seed)
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
    if format == "particles":
        # Linear gradient: headings spread uniformly over [0, 2π) → φ≈0.
        positions = rng.uniform(0.0, PARTICLES_DEFAULT_BOX, (n, 2))
        headings_static = np.linspace(-np.pi, np.pi, n, endpoint=False)
        velocities = PARTICLES_DEFAULT_SPEED * np.column_stack([
            np.cos(headings_static), np.sin(headings_static)
        ])
        out = []
        for t in range(n_steps):
            positions = (positions + velocities) % PARTICLES_DEFAULT_BOX
            out.append({
                "positions": positions.copy(),
                "headings": headings_static.copy(),
                "velocities": velocities.copy(),
                "speeds": np.full(n, PARTICLES_DEFAULT_SPEED, dtype=np.float64),
                "step": t,
            })
        return out
    if format == "sequence":
        arr_static = np.arange(n, dtype=np.int32) % 2
        arr = np.broadcast_to(arr_static.astype(np.int8), (n_steps, n)).copy()
        return _sequence_history_from_array(arr)
    if format == "opinions":
        # Opinions linearly spaced from 0 to 1 — uniform distribution, unimodal.
        ops_static = np.linspace(0.0, 1.0, OPINIONS_DEFAULT_N, endpoint=True)
        ops = np.broadcast_to(ops_static, (n_steps, OPINIONS_DEFAULT_N)).copy()
        return _opinions_history_from_array(ops)
    if format == "scalar_timeseries":
        # Uncontrolled monotonic drift under perturbation — x drifts at
        # perturbation rate post-onset. growth_ratio >> 2.0 → screening fails.
        x = _scalar_ts_uncontrolled(rng, n_steps, noise_std=0.1)
        return _scalar_ts_history(x, n_steps)
    if format == "noise_sweep":
        # Monotonically increasing performance: x = k * signal at each level,
        # where k increases with noise → no interior peak → screening fails.
        return _noise_sweep_monotone(rng, direction="increasing")
    if format == "choice_timeseries":
        # Linear ramp attendance from 0 to N — monotonic drift, not anti-coordination.
        att = np.linspace(0, CHOICE_TS_N_AGENTS, n_steps).round().astype(int)
        return _choice_ts_history(att)
    if format == "state_vector":
        return _state_vector_null_history(rng)
    if format == "canalization_bundle":
        return _canalization_bundle_null(rng, dynamics="divergent")
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
    rng = np.random.default_rng(seed)
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
    if format == "particles":
        # Checkerboard: alternating headings 0 and π → φ≈0 (antiparallel).
        positions = rng.uniform(0.0, PARTICLES_DEFAULT_BOX, (n, 2))
        headings_static = np.where(np.arange(n) % 2 == 0, 0.0, np.pi)
        velocities = PARTICLES_DEFAULT_SPEED * np.column_stack([
            np.cos(headings_static), np.sin(headings_static)
        ])
        out = []
        for t in range(n_steps):
            positions = (positions + velocities) % PARTICLES_DEFAULT_BOX
            out.append({
                "positions": positions.copy(),
                "headings": headings_static.copy(),
                "velocities": velocities.copy(),
                "speeds": np.full(n, PARTICLES_DEFAULT_SPEED, dtype=np.float64),
                "step": t,
            })
        return out
    if format == "sequence":
        arr_static = (np.arange(n) % 2).astype(np.int8)
        arr = np.broadcast_to(arr_static, (n_steps, n)).copy()
        return _sequence_history_from_array(arr)
    if format == "opinions":
        # Opinions alternating between two monotone halves (0→0.5 then 0.5→1).
        # Produces a flat uniform distribution, not bimodal — P21 should reject.
        half = OPINIONS_DEFAULT_N // 2
        ops_static = np.concatenate([
            np.linspace(0.0, 0.5, half, endpoint=False),
            np.linspace(0.5, 1.0, OPINIONS_DEFAULT_N - half, endpoint=True),
        ])
        ops = np.broadcast_to(ops_static, (n_steps, OPINIONS_DEFAULT_N)).copy()
        return _opinions_history_from_array(ops)
    if format == "scalar_timeseries":
        # Uncontrolled drift with oscillatory noise — x drifts post-onset.
        # growth_ratio >> 2.0 → screening fails.
        x = _scalar_ts_uncontrolled(rng, n_steps, noise_std=0.2)
        return _scalar_ts_history(x, n_steps)
    if format == "noise_sweep":
        # Monotonically decreasing performance: x = k * signal at each level,
        # where k decreases with noise → no interior peak → screening fails.
        return _noise_sweep_monotone(rng, direction="decreasing")
    if format == "choice_timeseries":
        # Alternating high/low attendance — deterministic oscillation, variance > baseline.
        half = CHOICE_TS_N_AGENTS
        att = np.tile([0, half], n_steps // 2 + 1)[:n_steps].astype(int)
        return _choice_ts_history(att)
    if format == "state_vector":
        # Alternating ±1 blocks — structured but no memory.
        N = STATE_VECTOR_DEFAULT_N
        P = STATE_VECTOR_DEFAULT_P
        stored_patterns = rng.choice([-1, 1], size=(P, N)).astype(np.int8)
        state = np.array([1 if i % 2 == 0 else -1 for i in range(N)], dtype=np.int8)
        history: List[Dict[str, Any]] = []
        for trial in range(STATE_VECTOR_DEFAULT_TRIALS):
            for step in range(STATE_VECTOR_DEFAULT_STEPS_PER_TRIAL):
                overlap = float(np.dot(stored_patterns[trial % P].astype(np.int32),
                                       state.astype(np.int32)) / N)
                history.append({
                    'state': state.copy(), 'step': step, 'trial': trial,
                    'cue_pattern_idx': trial % P, 'overlap': overlap,
                    'stored_patterns': stored_patterns.copy(), 'converged': True,
                })
        return history
    if format == "canalization_bundle":
        return _canalization_bundle_null(rng, dynamics="constant")
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
