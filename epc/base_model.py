"""Base model interface.

Every model subclasses BaseModel and implements:
- setup(): initialize state from parameters
- step(): advance one timestep, return state snapshot
- run(): execute full simulation, return state history
- get_metadata(): return model rules and parameters for metadata-required detectors

State history is a list of dicts with standardized keys:
    positions, velocities, types, grid, wealth, opinion, strategy, trail,
    scent, active, theta (phases), etc.

Not every model uses every key — each model documents which keys it provides.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

import numpy as np


class BaseModel(ABC):
    """Abstract base for all agent-based models.

    Parameters
    ----------
    seed : int
        Random seed for reproducibility.
    params : dict
        Model-specific parameters.
    """

    def __init__(self, seed: int = 42, **params: Any) -> None:
        self.seed = seed
        self.params = params
        self.rng = np.random.default_rng(seed)
        self._state_history: list[dict[str, Any]] = []
        self._step_count: int = 0
        self._is_setup: bool = False

    @abstractmethod
    def setup(self) -> dict[str, Any]:
        """Initialize model state.

        Returns the initial state snapshot (same format as step()).
        Must set self._is_setup = True.
        """
        ...

    @abstractmethod
    def step(self) -> dict[str, Any]:
        """Advance one timestep.

        Returns a state snapshot dict with standardized keys.
        Implementations should increment self._step_count.
        """
        ...

    @abstractmethod
    def get_metadata(self) -> dict[str, Any]:
        """Return model rules, parameters, and interaction structure.

        Used by metadata-required and metadata-assisted detectors.
        Should include at minimum: model_name, model_class, params,
        interaction_type, update_mode.
        """
        ...

    def run(self, n_steps: int, record_every: int = 1) -> list[dict[str, Any]]:
        """Execute full simulation.

        Parameters
        ----------
        n_steps : int
            Number of timesteps to run.
        record_every : int
            Record state every N steps (1 = every step).

        Returns
        -------
        list[dict]
            State history — list of state snapshots.
        """
        if not self._is_setup:
            initial_state = self.setup()
            self._state_history = [initial_state]

        for i in range(n_steps):
            state = self.step()
            if (i + 1) % record_every == 0:
                self._state_history.append(state)

        return self._state_history

    @property
    def state_history(self) -> list[dict[str, Any]]:
        """Access recorded state history."""
        return self._state_history

    @property
    def step_count(self) -> int:
        return self._step_count

    def is_converged(self) -> bool:
        """Check if model has reached a steady state.

        Default: False. Override for models with natural convergence criteria.
        """
        return False

    def get_timescale(self) -> float:
        """Return the system-intrinsic timescale τ in timesteps.

        Must be overridden by each model with the appropriate timescale
        from the scale-normalization convention:
        T_osc, T_prop, T_cross, T_traverse, T_sort, T_gen, T_relax.
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement get_timescale()"
        )
