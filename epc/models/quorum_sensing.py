"""
Quorum sensing / threshold-activated collective response models.

Two independent implementations for T1b cross-model generalization:

1. AutoinducerQuorum — Mean-field bacterial quorum-sensing model: agents
   produce autoinducer at a constant rate; local concentration scales with
   agent density. When concentration crosses activation_threshold, agents
   switch ON and produce at an enhanced rate (positive feedback). The
   deactivation_threshold is LOWER than activation, creating hysteresis:
   the system turns on at high density but stays on at lower densities
   (because enhanced production keeps concentration above deactivation).
   The model sweeps density up then down to expose hysteresis.

2. FractionThresholdModel — Simple fraction-threshold model: agents
   switch ON when a density-amplified effective fraction exceeds a
   threshold. The threshold for switching OFF is lower (explicit
   hysteresis). Provides an independent implementation for T1b.

Both models return (density, concentration, collective_state, fraction_on)
time series for up-sweep and down-sweep of density.

Reference: Waters, C. M., & Bassler, B. L. (2005). Quorum sensing:
cell-to-cell communication in bacteria. Annual Review of Cell and
Developmental Biology, 21, 319–346.
Nealson, K. H., Platt, T., & Hastings, J. W. (1970). Cellular control
of the synthesis and activity of the bacterial luminescent system.
Journal of Bacteriology, 104(1), 313–322.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class AutoinducerParams:
    """Parameters for the autoinducer quorum-sensing model.

    At equilibrium with all OFF: C_eq = density * base_production / degradation.
    Switch ON when C_eq > activation_threshold →
      density_ON = activation_threshold * degradation / base_production.
    At equilibrium with all ON: C_eq = density * total_production / degradation.
      where total_production = base_production + enhanced_production.
    Switch OFF when C_eq < deactivation_threshold →
      density_OFF = deactivation_threshold * degradation / total_production.
    Hysteresis: density_OFF < density_ON (requires activation > deactivation
    and enhanced_production > 0).
    """
    base_production: float = 1.0     # Autoinducer production rate (OFF state)
    enhanced_production: float = 4.0  # ADDITIONAL production when ON (positive feedback)
    degradation_rate: float = 1.0    # Autoinducer decay rate
    activation_threshold: float = 1.5  # Concentration to switch ON
    deactivation_threshold: float = 1.0  # Concentration to switch OFF (< activation)
    noise_std: float = 0.02         # Agent-level switching noise
    dt: float = 0.05                # Integration timestep
    n_steps_per_density: int = 300  # Steps to equilibrate at each density
    n_density_levels: int = 40      # Number of density levels per sweep direction
    density_min: float = 0.1        # Min density
    density_max: float = 3.0        # Max density
    seed: int = 42


@dataclass
class FractionThresholdParams:
    """Parameters for the fraction-threshold model."""
    n_agents: int = 200              # Number of agents
    activation_fraction: float = 0.5  # Effective fraction needed to switch ON
    deactivation_fraction: float = 0.2  # Effective fraction to stay ON
    noise_std: float = 0.03          # Per-agent decision noise
    n_steps_per_density: int = 200   # Steps per density level
    n_density_levels: int = 40
    density_min: float = 0.1
    density_max: float = 3.0
    density_scaling: float = 1.0     # Density at which coupling = 1.0
    dt: float = 0.1
    seed: int = 42


def _density_sweep(
    n_levels: int,
    d_min: float,
    d_max: float,
) -> List[float]:
    """Generate up-then-down density sweep."""
    up = np.linspace(d_min, d_max, n_levels).tolist()
    down = np.linspace(d_max, d_min, n_levels).tolist()
    return up + down


class AutoinducerQuorum:
    """Mean-field autoinducer quorum-sensing model.

    Population-level ODE for autoinducer concentration C:
      dC/dt = density * [f_off * r_base + f_on * (r_base + r_enh)] - γ * C

    where f_on is the fraction of ON agents, r_base and r_enh are
    production rates, and γ is the degradation rate.

    At each density level the system equilibrates for n_steps_per_density
    timesteps. The density is swept up then down.

    Parameters
    ----------
    params : AutoinducerParams
        Model parameters.
    """

    def __init__(self, params: Optional[AutoinducerParams] = None) -> None:
        if params is None:
            params = AutoinducerParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run density-sweep simulation (up then down).

        Returns
        -------
        list[dict]
            One dict per timestep per density level, with keys:
            'density', 'concentration', 'collective_state' (0 or 1),
            'fraction_on', 'step', 'density_idx', 'sweep_direction'.
        """
        p = self.params
        densities = _density_sweep(p.n_density_levels, p.density_min, p.density_max)
        n_up = p.n_density_levels  # first half is up-sweep

        history: List[Dict[str, Any]] = []

        # Population state
        n_agents = 200  # fixed population, density controls production
        agent_on = np.zeros(n_agents, dtype=bool)
        concentration = 0.0

        for d_idx, density in enumerate(densities):
            direction = 'up' if d_idx < n_up else 'down'

            for step_i in range(p.n_steps_per_density):
                f_on = float(np.mean(agent_on))
                f_off = 1.0 - f_on

                # Production rate: density modulates total autoinducer output
                total_prod = density * (
                    f_off * p.base_production
                    + f_on * (p.base_production + p.enhanced_production)
                )
                # ODE: dC/dt = production - degradation
                d_conc = total_prod - p.degradation_rate * concentration
                noise = p.noise_std * self.rng.normal() * np.sqrt(p.dt)
                concentration = max(0.0, concentration + d_conc * p.dt + noise)

                # Agent switching (individual with noise)
                for i in range(n_agents):
                    agent_noise = p.noise_std * self.rng.normal()
                    eff_conc = concentration + agent_noise
                    if not agent_on[i]:
                        if eff_conc > p.activation_threshold:
                            agent_on[i] = True
                    else:
                        if eff_conc < p.deactivation_threshold:
                            agent_on[i] = False

                f_on = float(np.mean(agent_on))
                collective = 1 if f_on > 0.5 else 0

                history.append({
                    'density': float(density),
                    'concentration': float(concentration),
                    'collective_state': collective,
                    'fraction_on': f_on,
                    'step': step_i,
                    'density_idx': d_idx,
                    'sweep_direction': direction,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        p = self.params
        return {
            'model_name': 'autoinducer_quorum',
            'model_class': 'quorum_sensing',
            'activation_threshold': p.activation_threshold,
            'deactivation_threshold': p.deactivation_threshold,
            'base_production': p.base_production,
            'enhanced_production': p.enhanced_production,
            'degradation_rate': p.degradation_rate,
            'dt': p.dt,
            'seed': p.seed,
            'substrate': 'density_sweep_timeseries',
            'has_threshold_activation': True,
            'has_hysteresis': p.activation_threshold > p.deactivation_threshold,
            'has_positive_feedback': p.enhanced_production > 0,
            'reference': 'Waters & Bassler (2005)',
        }


class FractionThresholdModel:
    """Simple fraction-threshold quorum model.

    Agents switch ON when the density-amplified effective ON fraction
    exceeds activation_fraction, and switch OFF when it drops below
    deactivation_fraction. The separation between these thresholds
    creates hysteresis.

    Parameters
    ----------
    params : FractionThresholdParams
        Model parameters.
    """

    def __init__(self, params: Optional[FractionThresholdParams] = None) -> None:
        if params is None:
            params = FractionThresholdParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run density-sweep simulation.

        Returns
        -------
        list[dict]
            Same schema as AutoinducerQuorum.
        """
        p = self.params
        densities = _density_sweep(p.n_density_levels, p.density_min, p.density_max)
        n_up = p.n_density_levels

        history: List[Dict[str, Any]] = []
        agent_on = np.zeros(p.n_agents, dtype=bool)

        for d_idx, density in enumerate(densities):
            direction = 'up' if d_idx < n_up else 'down'
            # Density coupling: higher density amplifies social pressure
            coupling = min(density / p.density_scaling, 5.0)

            for step_i in range(p.n_steps_per_density):
                f_on = float(np.mean(agent_on))
                # Effective signal = social fraction * coupling + density seed
                # The density seed represents background signal that grows with
                # population density, bootstrapping the collective switch
                density_seed = 0.35 * density / p.density_scaling
                effective = f_on * coupling + density_seed

                # Agent switching with noise
                for i in range(p.n_agents):
                    noise = p.noise_std * self.rng.normal()
                    eff = effective + noise
                    if not agent_on[i]:
                        if eff > p.activation_fraction:
                            agent_on[i] = True
                    else:
                        if eff < p.deactivation_fraction:
                            agent_on[i] = False

                f_on = float(np.mean(agent_on))
                concentration = f_on * coupling  # proxy signal
                collective = 1 if f_on > 0.5 else 0

                history.append({
                    'density': float(density),
                    'concentration': float(concentration),
                    'collective_state': collective,
                    'fraction_on': f_on,
                    'step': step_i,
                    'density_idx': d_idx,
                    'sweep_direction': direction,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        p = self.params
        return {
            'model_name': 'fraction_threshold',
            'model_class': 'quorum_sensing',
            'n_agents': p.n_agents,
            'activation_fraction': p.activation_fraction,
            'deactivation_fraction': p.deactivation_fraction,
            'noise_std': p.noise_std,
            'density_scaling': p.density_scaling,
            'dt': p.dt,
            'seed': p.seed,
            'substrate': 'density_sweep_timeseries',
            'has_threshold_activation': True,
            'has_hysteresis': p.activation_fraction > p.deactivation_fraction,
            'has_positive_feedback': False,
            'reference': 'threshold-activation model',
        }


class GradedResponseModel:
    """Linear/graded response model (negative control).

    Collective response increases linearly with density — no sharp
    threshold, no hysteresis. Used as a negative control.

    Parameters
    ----------
    slope : float
        Response per unit density.
    noise_std : float
    n_steps_per_density : int
    n_density_levels : int
    density_min, density_max : float
    seed : int
    """

    def __init__(
        self,
        slope: float = 0.25,
        noise_std: float = 0.05,
        n_steps_per_density: int = 100,
        n_density_levels: int = 40,
        density_min: float = 0.1,
        density_max: float = 3.0,
        seed: int = 42,
    ) -> None:
        self.slope = slope
        self.noise_std = noise_std
        self.n_steps_per_density = n_steps_per_density
        self.n_density_levels = n_density_levels
        self.density_min = density_min
        self.density_max = density_max
        self.seed = seed
        self.rng = np.random.default_rng(seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run density-sweep with graded (linear) response.

        Returns
        -------
        list[dict]
            Same schema as AutoinducerQuorum.
        """
        densities = _density_sweep(
            self.n_density_levels, self.density_min, self.density_max,
        )
        n_up = self.n_density_levels

        history: List[Dict[str, Any]] = []

        for d_idx, density in enumerate(densities):
            direction = 'up' if d_idx < n_up else 'down'
            for step_i in range(self.n_steps_per_density):
                noise = self.noise_std * self.rng.normal()
                fraction_on = float(np.clip(self.slope * density + noise, 0.0, 1.0))
                concentration = density * fraction_on
                collective = 1 if fraction_on > 0.5 else 0

                history.append({
                    'density': float(density),
                    'concentration': float(concentration),
                    'collective_state': collective,
                    'fraction_on': fraction_on,
                    'step': step_i,
                    'density_idx': d_idx,
                    'sweep_direction': direction,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'graded_response',
            'model_class': 'graded_response',
            'slope': self.slope,
            'noise_std': self.noise_std,
            'seed': self.seed,
            'substrate': 'density_sweep_timeseries',
            'has_threshold_activation': False,
            'has_hysteresis': False,
            'has_positive_feedback': False,
            'reference': 'linear response (negative control)',
        }
