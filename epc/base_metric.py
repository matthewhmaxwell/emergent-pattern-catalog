"""Base metric interface.

Metrics compute quantitative measurements from state histories.
Detectors consume metric outputs to make detection decisions.

The distinction: a metric computes Moran's I. A detector uses Moran's I
plus null models, secondaries, and exclusions to decide if P1 is present.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

import numpy as np


class BaseMetric(ABC):
    """Abstract base for all metrics.

    Parameters
    ----------
    name : str
        Human-readable metric name, e.g. 'morans_i'.
    """

    def __init__(self, name: str) -> None:
        self.name = name

    @abstractmethod
    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Compute metric from state history.

        Parameters
        ----------
        state_history : list[dict]
            Sequence of state snapshots.
        **kwargs
            Metric-specific parameters (e.g., window sizes, thresholds).

        Returns
        -------
        dict
            Named measurements. At minimum includes the metric's primary
            value under the key self.name.
        """
        ...

    @abstractmethod
    def required_keys(self) -> list[str]:
        """State history keys this metric needs, e.g. ['positions', 'types']."""
        ...

    def validate_history(self, state_history: list[dict[str, Any]]) -> list[str]:
        """Check that state history contains required keys.

        Returns list of missing keys (empty if all present).
        """
        if not state_history:
            return ["empty state history"]
        first = state_history[0]
        return [k for k in self.required_keys() if k not in first]

    def compute_timeseries(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> np.ndarray:
        """Compute metric value at each timestep.

        Default implementation calls compute() on single-frame histories.
        Override for efficiency when the metric supports incremental computation.

        Returns
        -------
        np.ndarray
            1D array of metric values, one per timestep.
        """
        values = []
        for i, state in enumerate(state_history):
            result = self.compute([state], **kwargs)
            values.append(result.get(self.name, np.nan))
        return np.array(values)
