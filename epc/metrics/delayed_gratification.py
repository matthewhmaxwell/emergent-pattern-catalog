"""P31 — Delayed gratification metrics (monotonicity-based).

Matches Zhang et al.'s definition: DG is computed from the monotonicity error
trajectory, measuring episodes where global sortedness temporarily worsens
then recovers. The DG value is the accumulated sum of (net_gain/backtrack_depth)
ratios across backtracking episodes.

For state histories from the sorting model, this uses per-swap monotonicity
error traces. For other systems, the "error" can be any objective function
where lower = better (distance from target, disorder metric, etc.).

Required state history keys (sorting model):
    monotonicity_error : int — per-state monotonicity error count
    OR: use compute_from_trace() with a pre-computed error trajectory
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_metric import BaseMetric
from epc.models.cell_view_sorting import compute_delayed_gratification, get_monotonicity_error


class DelayedGratification(BaseMetric):
    """Zhang's delayed gratification metric.

    Operates on a trajectory of error/objective values. Detects episodes
    where the error temporarily increases (backtracking) then decreases
    more than the increase (net gain from the detour).
    """

    def __init__(self) -> None:
        super().__init__(name="dg_index")

    def required_keys(self) -> list[str]:
        return ["monotonicity_error"]

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Compute DG from state history with monotonicity_error per state.

        For per-swap resolution, pass the full state history from
        model.run_to_completion(). For coarser resolution, pass sampled states.
        """
        errors = [s["monotonicity_error"] for s in state_history]
        return self.compute_from_trace(errors)

    def compute_from_trace(self, error_trace: list[int]) -> dict[str, Any]:
        """Compute DG from a pre-computed error trajectory.

        This is the primary entry point for sorting models that track
        monotonicity error per swap (via model.mono_error_trace).
        """
        dg = compute_delayed_gratification(error_trace)

        # Also compute trajectory statistics
        if len(error_trace) >= 2:
            initial_error = error_trace[0]
            final_error = error_trace[-1]
            max_error = max(error_trace)
            min_error = min(error_trace)
        else:
            initial_error = final_error = max_error = min_error = 0

        # Count backtracking episodes
        deduped = [error_trace[0]]
        for v in error_trace[1:]:
            if v != deduped[-1]:
                deduped.append(v)
        n_increases = sum(1 for i in range(1, len(deduped)) if deduped[i] > deduped[i - 1])

        return {
            "dg_index": dg,
            "initial_error": initial_error,
            "final_error": final_error,
            "max_error": max_error,
            "min_error": min_error,
            "trace_length": len(error_trace),
            "n_backtrack_steps": n_increases,
        }


class DGConditionComparison(BaseMetric):
    """Compare DG across conditions (e.g., different frozen cell counts).

    This tests Zhang's key finding: DG increases with more obstacles,
    indicating context-sensitive problem-solving rather than random noise.
    """

    def __init__(self) -> None:
        super().__init__(name="dg_condition_comparison")

    def required_keys(self) -> list[str]:
        return []

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Not meaningful for single run — use compare_conditions() instead."""
        return {}

    def compare_conditions(
        self,
        condition_traces: dict[str, list[list[int]]],
    ) -> dict[str, Any]:
        """Compare DG across named conditions.

        Parameters
        ----------
        condition_traces : dict[str, list[list[int]]]
            Maps condition name (e.g., '0_frozen', '1_frozen') to a list
            of error traces (one per trial).

        Returns
        -------
        dict with per-condition DG statistics and trend analysis.
        """
        results = {}
        means = []
        condition_names = []

        for name, traces in condition_traces.items():
            dg_values = [compute_delayed_gratification(t) for t in traces]
            arr = np.array(dg_values)
            results[name] = {
                "dg_mean": float(arr.mean()),
                "dg_std": float(arr.std()),
                "dg_median": float(np.median(arr)),
                "n_trials": len(traces),
            }
            means.append(float(np.median(arr)))
            condition_names.append(name)

        # Check for monotonic increase
        is_increasing = all(means[i] <= means[i + 1] for i in range(len(means) - 1))
        total_increase = means[-1] - means[0] if len(means) >= 2 else 0.0

        results["trend"] = {
            "monotonically_increasing": is_increasing,
            "total_increase": total_increase,
            "condition_order": condition_names,
            "median_values": means,
        }

        return results
