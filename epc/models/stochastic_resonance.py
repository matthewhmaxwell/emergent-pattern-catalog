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
