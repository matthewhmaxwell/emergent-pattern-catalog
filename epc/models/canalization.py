"""
Canalized restoration / equifinality model.

A dynamical-landscape model with a canalized target basin: many random
initial conditions flow (via gradient + nonlinear dynamics) to the same
target macrostate. Optionally, a perturbation can displace the state
mid-trajectory and it restores.

Two model variants:
  - CanalizedLandscape: multi-dimensional gradient flow on a double-well-
    like potential with a dominant canalized basin. The target attractor is
    a fixed point; diverse ICs converge to it via deterministic dynamics.
  - MultiBasinGRN: a Boolean-inspired gene regulatory network with a
    dominant basin of attraction. Provides T1b cross-model validation.

The model is NOT a trivial constant map (all dynamics set to zero / map
everything to a point regardless of dynamics). Convergence is driven by
genuine gradient flow with a finite relaxation time that depends on IC
distance from the target.

Reference: Waddington, C. H. (1957). The Strategy of the Genes. Allen & Unwin.
           Huang, S. et al. (2005). Cell Fates as High-Dimensional Attractor
           States of a Complex Gene Regulatory Network. PRL 94, 128701.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class CanalizedLandscapeParams:
    """Parameters for the canalized landscape model."""
    n_dims: int = 10               # State-space dimensionality
    target: Optional[np.ndarray] = None  # Target attractor (None → random)
    basin_strength: float = 2.0    # Gradient steepness toward target
    noise_std: float = 0.0         # Process noise
    dt: float = 0.05               # Integration timestep
    n_ics: int = 20                # Number of initial conditions to simulate
    ic_spread: float = 5.0         # IC sampling radius (uniform ball)
    perturbation_time: Optional[int] = None   # Step at which to perturb (None → no perturbation)
    perturbation_magnitude: float = 3.0       # Perturbation displacement magnitude
    n_steps: int = 200             # Steps per trajectory
    seed: int = 42


class CanalizedLandscape:
    """Multi-dimensional gradient-flow landscape with a canalized target basin.

    dx/dt = -basin_strength * (x - target) - nonlinear_restoring(x) + noise

    The nonlinear restoring term adds a quartic potential component that
    creates a genuine basin structure (not just a trivial linear sink).
    The target is a stable fixed point of the full dynamics.

    Parameters
    ----------
    params : CanalizedLandscapeParams
        Model parameters.
    """

    def __init__(self, params: Optional[CanalizedLandscapeParams] = None) -> None:
        if params is None:
            params = CanalizedLandscapeParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

        # Set target attractor
        if params.target is not None:
            self.target = np.asarray(params.target, dtype=float)
        else:
            self.target = self.rng.standard_normal(params.n_dims)

    def _nonlinear_force(self, x: np.ndarray) -> np.ndarray:
        """Nonlinear restoring force (quartic potential gradient).

        V_nl(x) = 0.1 * |x - target|^4 / n_dims^2
        F_nl = -dV_nl/dx = -0.4 * |x - target|^2 * (x - target) / n_dims^2

        This ensures the basin is genuinely nonlinear: trajectories from
        far-away ICs experience stronger restoring than near ones, and
        the convergence is NOT just exponential decay.
        """
        delta = x - self.target
        r2 = np.sum(delta ** 2) / (self.params.n_dims ** 2)
        return -0.4 * r2 * delta

    def simulate(self) -> List[Dict[str, Any]]:
        """Simulate trajectories from diverse ICs.

        Returns
        -------
        list[dict]
            State history with keys: 'state', 'step', 'trial', 'ic',
            'target', 'distance_to_target', 'converged'.
        """
        p = self.params
        rng = np.random.default_rng(p.seed)

        history: List[Dict[str, Any]] = []

        for trial in range(p.n_ics):
            # Sample IC uniformly in a ball of radius ic_spread
            direction = rng.standard_normal(p.n_dims)
            direction /= np.linalg.norm(direction) + 1e-12
            radius = p.ic_spread * rng.random() ** (1.0 / p.n_dims)
            x = self.target + direction * radius

            ic = x.copy()

            for step in range(p.n_steps):
                dist = float(np.linalg.norm(x - self.target))
                converged = dist < 0.1

                history.append({
                    'state': x.copy(),
                    'step': step,
                    'trial': trial,
                    'ic': ic.copy(),
                    'target': self.target.copy(),
                    'distance_to_target': dist,
                    'converged': converged,
                })

                # Apply perturbation mid-trajectory if requested
                if p.perturbation_time is not None and step == p.perturbation_time:
                    pert_dir = rng.standard_normal(p.n_dims)
                    pert_dir /= np.linalg.norm(pert_dir) + 1e-12
                    x = x + pert_dir * p.perturbation_magnitude

                # Dynamics: gradient + nonlinear restoring
                force_linear = -p.basin_strength * (x - self.target)
                force_nonlinear = self._nonlinear_force(x)
                dxdt = force_linear + force_nonlinear

                noise = np.zeros(p.n_dims)
                if p.noise_std > 0:
                    noise = rng.normal(0, p.noise_std, p.n_dims)

                x = x + dxdt * p.dt + noise * np.sqrt(p.dt)

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'canalized_landscape',
            'model_class': 'equifinality',
            'n_dims': self.params.n_dims,
            'basin_strength': self.params.basin_strength,
            'ic_spread': self.params.ic_spread,
            'n_ics': self.params.n_ics,
            'n_steps': self.params.n_steps,
            'dt': self.params.dt,
            'noise_std': self.params.noise_std,
            'seed': self.params.seed,
            'substrate': 'canalization_landscape',
            'has_canalized_basin': True,
            'has_perturbation_recovery': self.params.perturbation_time is not None,
            'perturbation_time': self.params.perturbation_time,
            'perturbation_magnitude': self.params.perturbation_magnitude,
            'reference': 'Waddington (1957)',
        }


class MultiBasinGRN:
    """Boolean-inspired GRN with a dominant canalized attractor.

    A continuous-valued GRN where gene states x_i ∈ [0, 1] evolve under
    sigmoidal regulatory interactions. The network is constructed so that
    one attractor basin dominates (captures most of IC space).

    dx_i/dt = -decay * x_i + sigmoid(W @ x + bias_i)

    where W is a weight matrix designed to have a dominant fixed point.

    This provides a biologically-motivated T1b cross-model for P25.
    """

    def __init__(
        self,
        n_genes: int = 10,
        n_ics: int = 20,
        n_steps: int = 200,
        dt: float = 0.05,
        decay: float = 1.0,
        seed: int = 42,
    ) -> None:
        self.n_genes = n_genes
        self.n_ics = n_ics
        self.n_steps = n_steps
        self.dt = dt
        self.decay = decay
        self.seed = seed

        rng = np.random.default_rng(seed)

        # Construct weight matrix with a dominant attractor
        # Strategy: create a target pattern and build Hebbian-like weights
        self.target_pattern = rng.choice([0.0, 1.0], size=n_genes)
        # Hebbian weights: W_ij = (2*xi_i - 1)(2*xi_j - 1) / n_genes
        signed = 2.0 * self.target_pattern - 1.0
        self.W = 2.0 * np.outer(signed, signed) / n_genes
        np.fill_diagonal(self.W, 0.0)  # No self-regulation

        # Bias toward the target pattern (strong bias ensures deep basin)
        self.bias = 3.0 * signed

    @staticmethod
    def _sigmoid(x: np.ndarray) -> np.ndarray:
        """Logistic sigmoid."""
        return 1.0 / (1.0 + np.exp(-np.clip(x, -20, 20)))

    def simulate(self) -> List[Dict[str, Any]]:
        """Simulate GRN from diverse random ICs.

        Returns
        -------
        list[dict]
            State history matching canalized landscape observation contract.
        """
        rng = np.random.default_rng(self.seed)
        history: List[Dict[str, Any]] = []

        for trial in range(self.n_ics):
            # Random IC: uniform in [0, 1]^n_genes
            x = rng.random(self.n_genes)
            ic = x.copy()

            for step in range(self.n_steps):
                dist = float(np.linalg.norm(x - self.target_pattern))
                converged = dist < 0.1

                history.append({
                    'state': x.copy(),
                    'step': step,
                    'trial': trial,
                    'ic': ic.copy(),
                    'target': self.target_pattern.copy(),
                    'distance_to_target': dist,
                    'converged': converged,
                })

                # GRN dynamics: dx/dt = -decay * x + sigmoid(Wx + bias)
                activation = self.W @ x + self.bias
                dxdt = -self.decay * x + self._sigmoid(activation)
                x = np.clip(x + dxdt * self.dt, 0.0, 1.0)

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'multi_basin_grn',
            'model_class': 'equifinality',
            'n_genes': self.n_genes,
            'n_ics': self.n_ics,
            'n_steps': self.n_steps,
            'dt': self.dt,
            'decay': self.decay,
            'seed': self.seed,
            'substrate': 'canalization_landscape',
            'has_canalized_basin': True,
            'has_perturbation_recovery': False,
            'reference': 'Huang et al. (2005)',
        }


class DiffusiveDynamics:
    """Negative control: diffusive dynamics where ICs do NOT converge.

    dx/dt = noise (pure random walk / diffusion)

    Each IC drifts independently; final states have variance comparable
    to or larger than initial-state variance. This is the primary
    negative control for P25.
    """

    def __init__(
        self,
        n_dims: int = 10,
        n_ics: int = 20,
        n_steps: int = 200,
        dt: float = 0.05,
        diffusion_std: float = 1.0,
        seed: int = 42,
    ) -> None:
        self.n_dims = n_dims
        self.n_ics = n_ics
        self.n_steps = n_steps
        self.dt = dt
        self.diffusion_std = diffusion_std
        self.seed = seed

    def simulate(self) -> List[Dict[str, Any]]:
        """Simulate diffusive trajectories."""
        rng = np.random.default_rng(self.seed)
        target = np.zeros(self.n_dims)  # Arbitrary reference point
        history: List[Dict[str, Any]] = []

        for trial in range(self.n_ics):
            x = rng.standard_normal(self.n_dims) * 2.0
            ic = x.copy()

            for step in range(self.n_steps):
                dist = float(np.linalg.norm(x - target))
                history.append({
                    'state': x.copy(),
                    'step': step,
                    'trial': trial,
                    'ic': ic.copy(),
                    'target': target.copy(),
                    'distance_to_target': dist,
                    'converged': False,
                })
                # Pure diffusion
                x = x + rng.normal(0, self.diffusion_std, self.n_dims) * np.sqrt(self.dt)

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'diffusive_dynamics',
            'model_class': 'diffusion',
            'n_dims': self.n_dims,
            'n_ics': self.n_ics,
            'n_steps': self.n_steps,
            'dt': self.dt,
            'seed': self.seed,
            'substrate': 'canalization_landscape',
            'has_canalized_basin': False,
            'has_perturbation_recovery': False,
        }


class TrivialCollapse:
    """Negative control: trivial constant map (everything → same point).

    The 'dynamics' immediately set x = constant regardless of state.
    This is NOT canalization — there's no dynamics-driven convergence,
    just a trivial collapse that the detector must reject.

    Detection gate: the detector checks that convergence takes finite
    time (relaxation time > 0), not instant collapse.
    """

    def __init__(
        self,
        n_dims: int = 10,
        n_ics: int = 20,
        n_steps: int = 200,
        collapse_point: Optional[np.ndarray] = None,
        seed: int = 42,
    ) -> None:
        self.n_dims = n_dims
        self.n_ics = n_ics
        self.n_steps = n_steps
        self.seed = seed

        if collapse_point is not None:
            self.collapse_point = np.asarray(collapse_point, dtype=float)
        else:
            self.collapse_point = np.zeros(n_dims)

    def simulate(self) -> List[Dict[str, Any]]:
        """Simulate trivial collapse trajectories."""
        rng = np.random.default_rng(self.seed)
        history: List[Dict[str, Any]] = []

        for trial in range(self.n_ics):
            ic = rng.standard_normal(self.n_dims) * 5.0

            for step in range(self.n_steps):
                # Instant collapse: state is always the collapse point
                x = self.collapse_point.copy()
                dist = 0.0 if step > 0 else float(np.linalg.norm(ic - self.collapse_point))

                history.append({
                    'state': x.copy(),
                    'step': step,
                    'trial': trial,
                    'ic': ic.copy(),
                    'target': self.collapse_point.copy(),
                    'distance_to_target': dist,
                    'converged': True,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'trivial_collapse',
            'model_class': 'trivial',
            'n_dims': self.n_dims,
            'n_ics': self.n_ics,
            'n_steps': self.n_steps,
            'seed': self.seed,
            'substrate': 'canalization_landscape',
            'has_canalized_basin': False,
            'has_perturbation_recovery': False,
        }
