"""
Stochastic resonance models.

Two independent SR implementations for T1b cross-model generalization:

1. BistableDoubleWell — Overdamped particle in a symmetric double-well
   potential V(x) = -a*x²/2 + b*x⁴/4, driven by a weak subthreshold
   periodic signal A*sin(2πft) and Gaussian white noise of intensity D.
   At intermediate noise, the particle hops between wells in synchrony
   with the driving signal, producing a peak in output SNR.

2. ThresholdUnit — A simple threshold detector: output = 1 if
   (signal + noise > threshold) else 0. At intermediate noise,
   detection probability is maximized (noise occasionally pushes
   subthreshold signal above threshold without overwhelming it).

Both models sweep noise levels and return the output time series at
each level plus the driving signal, enabling the P26 detector to
identify the inverted-U SNR-vs-noise curve.

Reference: Collins, J. J., Chow, C. C., & Imhoff, T. T. (1995).
Stochastic resonance without tuning. Nature, 376(6537), 236–238.
Gammaitoni, L. et al. (1998). Stochastic resonance. Rev. Mod. Phys.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class DoubleWellParams:
    """Parameters for the bistable double-well SR model."""
    a: float = 4.0                # Quadratic coefficient (well depth)
    b: float = 1.0                # Quartic coefficient (confinement)
    signal_amplitude: float = 1.0  # Subthreshold periodic signal amplitude
    signal_frequency: float = 0.005  # Signal frequency (slow, in Hz)
    dt: float = 0.01              # Integration timestep
    n_steps: int = 40000          # Steps per noise level (2 signal periods)
    n_trials: int = 20            # Independent trials per noise level
    noise_levels: Optional[List[float]] = None  # Noise intensities to sweep
    seed: int = 42


@dataclass
class ThresholdUnitParams:
    """Parameters for the threshold-unit SR model."""
    threshold: float = 1.0        # Detection threshold
    signal_amplitude: float = 0.7  # Subthreshold (< threshold)
    signal_frequency: float = 0.05 # Signal frequency (Hz)
    dt: float = 0.01              # Timestep
    n_steps: int = 20000          # Steps per noise level
    noise_levels: Optional[List[float]] = None
    seed: int = 42


def _default_noise_levels() -> List[float]:
    """Default noise sweep: 15 levels from 0 to high noise.

    For the double-well with a=4, b=1 (barrier=4.0), the Kramers
    escape rate matches f_signal=0.005 at D ≈ 1.0–1.5. For the
    threshold unit with threshold=1.0 and amplitude=0.7, optimal
    noise is D ≈ 0.5–1.0. The sweep covers both regimes.
    """
    return [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5,
            2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 20.0]


class BistableDoubleWell:
    """Overdamped bistable double-well with periodic driving + noise.

    V(x) = -a*x²/2 + b*x⁴/4
    dx/dt = a*x - b*x³ + A*sin(2πft) + √(2D)*ξ(t)

    where D is noise intensity and ξ is unit white noise.

    The barrier height is a²/(4b). The signal must be subthreshold:
    A < a*√(a/(3b)) for the double-well to have two minima at all
    times (weak signal that cannot alone push x over the barrier).

    Parameters
    ----------
    params : DoubleWellParams
        Model parameters.
    """

    def __init__(self, params: Optional[DoubleWellParams] = None) -> None:
        if params is None:
            params = DoubleWellParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run noise-sweep simulation with multiple trials per level.

        For each noise level, runs n_trials independent simulations
        (each with a fresh random seed derived from the master seed)
        and concatenates their outputs. The detector averages across
        trials to get a clean performance estimate.

        Returns
        -------
        list[dict]
            One dict per timestep per trial per noise level, with keys:
            'time', 'x', 'signal', 'noise_level', 'step',
            'noise_level_idx'.
        """
        p = self.params
        noise_levels = p.noise_levels if p.noise_levels is not None else _default_noise_levels()

        history: List[Dict[str, Any]] = []
        well_pos = np.sqrt(p.a / p.b)

        for nl_idx, D in enumerate(noise_levels):
            for trial in range(p.n_trials):
                # Alternate starting well across trials
                x = -well_pos if trial % 2 == 0 else well_pos

                for step_i in range(p.n_steps):
                    t = step_i * p.dt
                    signal = p.signal_amplitude * np.sin(
                        2.0 * np.pi * p.signal_frequency * t
                    )

                    history.append({
                        'time': t,
                        'x': float(x),
                        'signal': float(signal),
                        'noise_level': float(D),
                        'noise_level_idx': nl_idx,
                        'step': step_i,
                    })

                    # Euler-Maruyama integration
                    drift = p.a * x - p.b * x**3 + signal
                    noise = 0.0
                    if D > 0:
                        noise = (np.sqrt(2.0 * D)
                                 * self.rng.normal()
                                 * np.sqrt(p.dt))
                    x = x + drift * p.dt + noise

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        p = self.params
        barrier = p.a**2 / (4.0 * p.b)
        return {
            'model_name': 'bistable_double_well',
            'model_class': 'stochastic_resonance',
            'a': p.a,
            'b': p.b,
            'signal_amplitude': p.signal_amplitude,
            'signal_frequency': p.signal_frequency,
            'barrier_height': barrier,
            'dt': p.dt,
            'n_steps': p.n_steps,
            'seed': p.seed,
            'substrate': 'noise_sweep_timeseries',
            'has_subthreshold_signal': True,
            'has_noise_sweep': True,
            'reference': 'Gammaitoni et al. (1998)',
        }


class ThresholdUnit:
    """Simple threshold detector for stochastic resonance.

    output = 1 if (signal(t) + noise(t)) > threshold else 0.

    The signal is subthreshold: amplitude < threshold. At intermediate
    noise intensity, the signal+noise sum crosses the threshold in
    synchrony with the signal, producing peak detection accuracy.

    This provides an independent SR implementation for T1b cross-model
    generalization testing.

    Parameters
    ----------
    params : ThresholdUnitParams
        Model parameters.
    """

    def __init__(self, params: Optional[ThresholdUnitParams] = None) -> None:
        if params is None:
            params = ThresholdUnitParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run noise-sweep simulation.

        Returns
        -------
        list[dict]
            One dict per timestep per noise level, with keys:
            'time', 'x' (binary output), 'signal', 'noise_level',
            'step', 'noise_level_idx'.
        """
        p = self.params
        noise_levels = p.noise_levels if p.noise_levels is not None else _default_noise_levels()

        history: List[Dict[str, Any]] = []

        for nl_idx, D in enumerate(noise_levels):
            for step_i in range(p.n_steps):
                t = step_i * p.dt
                signal = p.signal_amplitude * np.sin(2.0 * np.pi * p.signal_frequency * t)

                noise_val = 0.0
                if D > 0:
                    noise_val = D * self.rng.normal()

                # Threshold detection
                x = 1.0 if (signal + noise_val) > p.threshold else 0.0

                history.append({
                    'time': t,
                    'x': float(x),
                    'signal': float(signal),
                    'noise_level': float(D),
                    'noise_level_idx': nl_idx,
                    'step': step_i,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        p = self.params
        return {
            'model_name': 'threshold_unit',
            'model_class': 'stochastic_resonance',
            'threshold': p.threshold,
            'signal_amplitude': p.signal_amplitude,
            'signal_frequency': p.signal_frequency,
            'dt': p.dt,
            'n_steps': p.n_steps,
            'seed': p.seed,
            'substrate': 'noise_sweep_timeseries',
            'has_subthreshold_signal': p.signal_amplitude < p.threshold,
            'has_noise_sweep': True,
            'reference': 'Collins et al. (1995)',
        }


# ---------------------------------------------------------------------------
# Vectorised driver (repro v2, P26). Code addition only: the update rule is
# the one in BistableDoubleWell.simulate(), applied to all trials at once.
# Instead of one dict per step per trial it returns Fourier projections of
# x(t) over a measurement window, so a 3.8e5-step x 400-trial run fits in
# memory. simulate() and the catalog defaults are unchanged.
# ---------------------------------------------------------------------------

def double_well_fourier_driver(
    params: DoubleWellParams,
    D: float,
    n_trials: int,
    n_burn: int,
    n_meas: int,
    omegas: Optional[List[float]] = None,
    rng: Optional[np.random.Generator] = None,
    normals: Optional[np.ndarray] = None,
    chunk: int = 4096,
    return_trajectory: bool = False,
) -> Dict[str, Any]:
    """Integrate n_trials independent copies of the double well together.

    Update rule (identical to ``BistableDoubleWell.simulate``, Euler-Maruyama):
        t      = step_i * dt
        signal = A * sin(2 pi f t)
        drift  = a*x - b*x**3 + signal
        noise  = sqrt(2 D) * z * sqrt(dt),   z ~ N(0, 1)   (0 if D == 0)
        x      = x + drift*dt + noise
    with x recorded *before* the update at each step (as simulate() records
    it) and trial k started in the well -sqrt(a/b) (k even) or +sqrt(a/b)
    (k odd), as in simulate().

    Steps 0 .. n_burn-1 are burn-in; steps n_burn .. n_burn+n_meas-1 form the
    measurement window of length T = n_meas*dt. For every angular frequency
    w in ``omegas`` the driver returns, per trial,
        X_k(w) = sum_{window} x(t) exp(-i w t) dt
    (complex, shape (len(omegas), n_trials)). The caller forms
    c_k = (2/T) X_k(Omega).

    Normals come from ``normals`` (shape (n_burn+n_meas, n_trials), row =
    step) when given, else from ``rng.standard_normal`` drawn in chunks of
    ``chunk`` steps (row-major: step, then trial).
    """
    p = params
    n_total = int(n_burn) + int(n_meas)
    if omegas is None:
        omegas = [2.0 * np.pi * p.signal_frequency]
    omegas_arr = np.asarray(omegas, dtype=np.float64)
    if normals is not None:
        normals = np.asarray(normals, dtype=np.float64)
        if normals.shape != (n_total, n_trials):
            raise ValueError(
                f"normals shape {normals.shape} != {(n_total, n_trials)}"
            )
    elif rng is None:
        raise ValueError("need rng or normals")

    well_pos = np.sqrt(p.a / p.b)
    x = np.where(np.arange(n_trials) % 2 == 0, -well_pos, well_pos).astype(np.float64)

    s2D = np.sqrt(2.0 * D)
    sdt = np.sqrt(p.dt)
    proj_re = np.zeros((len(omegas_arr), n_trials), dtype=np.float64)
    proj_im = np.zeros((len(omegas_arr), n_trials), dtype=np.float64)
    traj = np.empty((n_total, n_trials), dtype=np.float64) if return_trajectory else None
    buf = np.empty((chunk, n_trials), dtype=np.float64)

    start = 0
    while start < n_total:
        m = min(chunk, n_total - start)
        steps = np.arange(start, start + m)
        t = steps * p.dt
        sig = p.signal_amplitude * np.sin(2.0 * np.pi * p.signal_frequency * t)
        if D > 0:
            z = normals[start:start + m] if normals is not None else rng.standard_normal((m, n_trials))
            noise = s2D * z * sdt
        else:
            noise = np.zeros((m, n_trials))
        for k in range(m):
            buf[k] = x
            drift = p.a * x - p.b * x**3 + sig[k]
            x = x + drift * p.dt + noise[k]
        if traj is not None:
            traj[start:start + m] = buf[:m]
        # Fourier projection over the part of this chunk inside the window
        lo = max(start, n_burn)
        hi = min(start + m, n_total)
        if hi > lo:
            tw = np.arange(lo, hi) * p.dt
            ph = np.outer(omegas_arr, tw)            # (n_omega, len)
            seg = buf[lo - start:hi - start]          # (len, n_trials)
            proj_re += (np.cos(ph) @ seg) * p.dt
            proj_im -= (np.sin(ph) @ seg) * p.dt
        start += m

    out: Dict[str, Any] = {
        "X": proj_re + 1j * proj_im,
        "omegas": omegas_arr,
        "T": n_meas * p.dt,
        "x_final": x,
    }
    if traj is not None:
        out["trajectory"] = traj
    return out
