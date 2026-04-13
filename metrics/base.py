"""Abstract base class for all detection metrics."""

from abc import ABC, abstractmethod
from typing import Any


class BaseMetric(ABC):
    """Abstract base class for all detection metrics.

    Each metric operates on a state history produced by a model run
    and returns a dictionary of named measurements.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable metric name."""
        ...

    @property
    @abstractmethod
    def pattern_id(self) -> str:
        """Pattern catalog ID this metric detects (e.g., 'P1', 'P2')."""
        ...

    @abstractmethod
    def compute(self, history: list[dict[str, Any]], **kwargs) -> dict[str, Any]:
        """Compute the metric from a state history.

        Args:
            history: List of state dictionaries from a model run.
            **kwargs: Additional metric-specific parameters.

        Returns:
            Dict with named measurements, e.g.:
            {'peak_aggregation': 0.72, 'time_to_peak': 42,
             'significant': True, 'p_value': 0.001}
        """
        ...

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(pattern={self.pattern_id})"
