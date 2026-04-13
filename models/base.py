"""Abstract base class for all minimal agent-based models."""

from abc import ABC, abstractmethod
from typing import Any, Callable

import numpy as np


class BaseModel(ABC):
    """Abstract base class for all minimal agent-based models.

    All models maintain a seeded RNG for reproducibility and record
    their full state history during execution.
    """

    def __init__(self, seed: int = 42):
        self.rng = np.random.default_rng(seed)
        self.history: list[dict[str, Any]] = []
        self.step_count: int = 0

    @abstractmethod
    def setup(self, **kwargs) -> None:
        """Initialize the model state."""
        ...

    @abstractmethod
    def step(self) -> None:
        """Execute one time step."""
        ...

    @abstractmethod
    def get_state(self) -> dict[str, Any]:
        """Return current state as a dictionary."""
        ...

    def run(
        self,
        max_steps: int = 10000,
        stop_condition: Callable[[dict[str, Any]], bool] | None = None,
    ) -> list[dict[str, Any]]:
        """Run the model, recording state history.

        Args:
            max_steps: Maximum number of steps to execute.
            stop_condition: Optional callable that receives state dict
                and returns True to stop early.

        Returns:
            List of state dictionaries (one per step, including initial state).
        """
        self.history = [self.get_state()]
        for _ in range(max_steps):
            self.step()
            self.step_count += 1
            state = self.get_state()
            self.history.append(state)
            if stop_condition and stop_condition(state):
                break
        return self.history
