"""
Minimal homeostatic regulation model.

A scalar regulated variable x is driven by:
  (i)  external perturbation schedule (step, sustained, pulse)
  (ii) active negative-feedback corrective response proportional to deviation
       from a set-point: dx/dt = -gain * (x - setpoint) + perturbation(t)

With gain=0 the system is uncontrolled (variable drifts under perturbation).
With gain>0 the system exhibits homeostatic regulation: bounded deviation,
finite recovery time, and deviation integral ≪ uncontrolled baseline.

Two controller variants are provided:
  - ProportionalHomeostat: dx/dt = -gain * (x - setpoint) + perturbation(t)
  - IntegralHomeostat: uses integral control (accumulates error over time)
The second variant supports T1b cross-model generalization testing.

Reference: Ashby, W. R. (1956). An Introduction to Cybernetics. Chapman & Hall.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class HomeostatParams:
    """Parameters for the proportional homeostat model."""
    setpoint: float = 10.0          # Target regulated variable value
    gain: float = 0.5               # Feedback gain (0 = no regulation)
    dt: float = 0.1                 # Integration timestep
    noise_std: float = 0.0          # Process noise (Gaussian)
    seed: int = 42                  # RNG seed
    x0: Optional[float] = None      # Initial value (None → setpoint)


@dataclass
class PerturbationSchedule:
    """External perturbation schedule."""
    onset: float = 50.0             # Time at which perturbation begins
    offset: Optional[float] = None  # Time at which perturbation ends (None → sustained)
    amplitude: float = 5.0          # Perturbation magnitude


def _perturbation_value(t: float, schedule: PerturbationSchedule) -> float:
    """Return perturbation value at time t."""
    if t < schedule.onset:
        return 0.0
    if schedule.offset is not None and t > schedule.offset:
        return 0.0
    return schedule.amplitude


class ProportionalHomeostat:
    """Proportional negative-feedback homeostat.

    dx/dt = -gain * (x - setpoint) + perturbation(t) + noise

    Parameters
    ----------
    params : HomeostatParams
        Model parameters.
    """

    def __init__(self, params: Optional[HomeostatParams] = None) -> None:
        if params is None:
            params = HomeostatParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)
        self.x = params.x0 if params.x0 is not None else params.setpoint

    def simulate(
        self,
        n_steps: int,
        schedule: Optional[PerturbationSchedule] = None,
    ) -> List[Dict[str, Any]]:
        """Run simulation and return state history.

        Parameters
        ----------
        n_steps : int
            Number of integration steps.
        schedule : PerturbationSchedule, optional
            Perturbation schedule. If None, no perturbation.

        Returns
        -------
        list[dict]
            State history with keys: 'time', 'x', 'setpoint',
            'perturbation', 'deviation', 'step'.
        """
        if schedule is None:
            schedule = PerturbationSchedule(onset=float('inf'))

        p = self.params
        x = self.x
        history: List[Dict[str, Any]] = []

        for step_i in range(n_steps):
            t = step_i * p.dt
            pert = _perturbation_value(t, schedule)
            deviation = x - p.setpoint

            history.append({
                'time': t,
                'x': float(x),
                'setpoint': p.setpoint,
                'perturbation': float(pert),
                'deviation': float(deviation),
                'step': step_i,
            })

            # Euler integration
            dxdt = -p.gain * deviation + pert
            noise = 0.0
            if p.noise_std > 0:
                noise = self.rng.normal(0, p.noise_std)
            x = x + dxdt * p.dt + noise * np.sqrt(p.dt)

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'proportional_homeostat',
            'model_class': 'homeostatic_regulation',
            'controller_type': 'proportional',
            'setpoint': self.params.setpoint,
            'gain': self.params.gain,
            'dt': self.params.dt,
            'noise_std': self.params.noise_std,
            'seed': self.params.seed,
            'substrate': 'scalar_timeseries',
            'has_active_feedback': self.params.gain > 0,
            'has_perturbation': True,
            'reference': 'Ashby (1956)',
        }


class IntegralHomeostat:
    """Integral-controller homeostat (PID with P=0, I=gain, D=0).

    dx/dt = perturbation(t) + control(t)
    control(t) = -gain * ∫(x - setpoint) dt

    The integral controller accumulates error over time and applies
    corrective action proportional to the accumulated error. This
    provides a second, independent homeostat implementation for T1b
    cross-model generalization testing.

    Parameters
    ----------
    params : HomeostatParams
        Model parameters. gain is the integral gain.
    """

    def __init__(self, params: Optional[HomeostatParams] = None) -> None:
        if params is None:
            params = HomeostatParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)
        self.x = params.x0 if params.x0 is not None else params.setpoint

    def simulate(
        self,
        n_steps: int,
        schedule: Optional[PerturbationSchedule] = None,
    ) -> List[Dict[str, Any]]:
        """Run simulation and return state history.

        Returns same observation bundle as ProportionalHomeostat.
        """
        if schedule is None:
            schedule = PerturbationSchedule(onset=float('inf'))

        p = self.params
        x = self.x
        error_integral = 0.0
        history: List[Dict[str, Any]] = []

        for step_i in range(n_steps):
            t = step_i * p.dt
            pert = _perturbation_value(t, schedule)
            deviation = x - p.setpoint

            history.append({
                'time': t,
                'x': float(x),
                'setpoint': p.setpoint,
                'perturbation': float(pert),
                'deviation': float(deviation),
                'step': step_i,
            })

            # Accumulate error integral
            error_integral += deviation * p.dt

            # Euler integration: dx/dt = perturbation + control
            control = -p.gain * error_integral
            dxdt = pert + control
            noise = 0.0
            if p.noise_std > 0:
                noise = self.rng.normal(0, p.noise_std)
            x = x + dxdt * p.dt + noise * np.sqrt(p.dt)

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'integral_homeostat',
            'model_class': 'homeostatic_regulation',
            'controller_type': 'integral',
            'setpoint': self.params.setpoint,
            'gain': self.params.gain,
            'dt': self.params.dt,
            'noise_std': self.params.noise_std,
            'seed': self.params.seed,
            'substrate': 'scalar_timeseries',
            'has_active_feedback': self.params.gain > 0,
            'has_perturbation': True,
            'reference': 'Ashby (1956)',
        }
