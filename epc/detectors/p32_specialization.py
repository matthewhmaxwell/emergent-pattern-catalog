"""P32 Emergent Specialization / Division of Labor detector.

Detects emergent role differentiation where initially identical agents
spontaneously specialize into distinct functional roles that increase
collective efficiency, without pre-assigned task allocation.

Primary metrics:
  - Per-agent behavioral entropy decline over time
  - Efficiency gain vs non-specialized baseline
  - Task-switching frequency decline

Three-tier detection:
  - Screening: mean per-agent task entropy decreases over time
  - Confirmation: agents settle into DISTINCT roles (low per-agent entropy
    AND coverage of all tasks across the population)
  - Definitive: collective efficiency exceeds non-specialized baseline

Distinctness from P23 (anti-coordination): P32 is role-FIXED specialization;
P23 is continuous re-balancing. From P19 (leadership): P32 is a stable task
identity; P19 leadership is transient.

Reference:
  Bonabeau, E., Theraulaz, G. & Deneubourg, J.-L. (1996).
  Quantitative study of the fixed threshold model for the regulation
  of division of labour in insect societies.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType


def extract_observation_bundle(history: list[dict[str, Any]]) -> dict[str, Any]:
    """T1a adapter: extract per-agent task-assignment time series.

    Parameters
    ----------
    history : list[dict]
        State history. Each dict must have 'task_assignments' (ndarray of
        int, -1 for idle), 'n_agents' (int), and 'n_tasks' (int).

    Returns
    -------
    dict with keys:
        task_assignments : ndarray(T, N) — task index per agent per step
        n_agents : int
        n_tasks : int
        steps : ndarray(T,)
    """
    task_assignments = np.array([s["task_assignments"] for s in history])
    steps = np.array([s.get("step", i) for i, s in enumerate(history)])
    n_agents = history[0]["n_agents"]
    n_tasks = history[0]["n_tasks"]
    return {
        "task_assignments": task_assignments,
        "n_agents": n_agents,
        "n_tasks": n_tasks,
        "steps": steps,
    }


def _per_agent_entropy_window(
    assignments: np.ndarray,
    n_tasks: int,
    start: int,
    end: int,
) -> np.ndarray:
    """Compute per-agent Shannon entropy over a window of timesteps.

    Parameters
    ----------
    assignments : ndarray(T, N)
    n_tasks : int
    start, end : int
        Window slice [start:end].

    Returns
    -------
    ndarray(N,) — entropy per agent in bits
    """
    window = assignments[start:end]
    n_agents = window.shape[1]
    counts = np.zeros((n_agents, n_tasks), dtype=np.float64)

    for t in range(window.shape[0]):
        for i in range(n_agents):
            task = window[t, i]
            if task >= 0:
                counts[i, task] += 1.0

    totals = counts.sum(axis=1, keepdims=True)
    totals = np.where(totals > 0, totals, 1.0)
    probs = counts / totals

    log_probs = np.where(probs > 0, np.log2(probs), 0.0)
    entropy = -np.sum(probs * log_probs, axis=1)
    return entropy


def _task_coverage(assignments: np.ndarray, n_tasks: int, start: int, end: int) -> float:
    """Fraction of tasks that have at least one agent assigned in the window."""
    window = assignments[start:end]
    tasks_seen = set()
    for t in range(window.shape[0]):
        for task in window[t]:
            if task >= 0:
                tasks_seen.add(int(task))
    return len(tasks_seen) / n_tasks


def _switching_frequency(assignments: np.ndarray, start: int, end: int) -> float:
    """Mean per-agent task-switching frequency in the window."""
    window = assignments[start:end]
    if window.shape[0] < 2:
        return 0.0
    n_agents = window.shape[1]
    switches = 0
    transitions = 0
    for i in range(n_agents):
        for t in range(1, window.shape[0]):
            if window[t, i] >= 0 and window[t - 1, i] >= 0:
                transitions += 1
                if window[t, i] != window[t - 1, i]:
                    switches += 1
    if transitions == 0:
        return 0.0
    return switches / transitions


class P32SpecializationDetector(BaseDetector):
    """Detector for P32: Emergent specialization / division of labor.

    Parameters
    ----------
    n_permutations : int
        Number of permutations for the null model.
    seed : int
        Random seed.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        seed: int = 42,
    ) -> None:
        super().__init__(
            pattern_id="P32",
            excluded_patterns=["P23"],  # anti-coordination: continuous re-balancing
            allowed_co_occurrences=[],
            observable_scope="state_history_only",
        )
        self.n_permutations = n_permutations
        self._seed = seed

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """Timescale: number of agents (one reinforcement cycle)."""
        if state_history and "n_agents" in state_history[0]:
            return float(state_history[0]["n_agents"])
        return float(len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        """Check minimum run length and required keys."""
        warnings = super()._validate_prerequisites(state_history, timescale)
        if not state_history:
            warnings.append("empty state history")
            return warnings
        if "task_assignments" not in state_history[0]:
            warnings.append("missing 'task_assignments' key in state history")
        if "n_tasks" not in state_history[0]:
            warnings.append("missing 'n_tasks' key in state history")
        return warnings

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Compute entropy decline: early vs late per-agent task entropy."""
        bundle = extract_observation_bundle(state_history)
        assignments = bundle["task_assignments"]
        n_tasks = bundle["n_tasks"]
        T = assignments.shape[0]

        window = max(T // 4, 10)
        early_entropy = _per_agent_entropy_window(assignments, n_tasks, 0, window)
        late_entropy = _per_agent_entropy_window(assignments, n_tasks, T - window, T)

        mean_early = float(np.mean(early_entropy))
        mean_late = float(np.mean(late_entropy))
        entropy_decline = mean_early - mean_late

        # Task coverage in late window
        late_coverage = _task_coverage(assignments, n_tasks, T - window, T)

        # Switching frequency decline
        early_switch = _switching_frequency(assignments, 0, window)
        late_switch = _switching_frequency(assignments, T - window, T)
        switch_decline = early_switch - late_switch

        return {
            "entropy_decline": entropy_decline,
            "mean_early_entropy": mean_early,
            "mean_late_entropy": mean_late,
            "late_coverage": late_coverage,
            "early_switch_freq": early_switch,
            "late_switch_freq": late_switch,
            "switch_decline": switch_decline,
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening: per-agent entropy decreases over time WITH maintained coverage.

        Content prerequisite (Sprint 91): P32 requires entropy decline AND
        maintained population-level task coverage. Single-task collapse
        (all agents do one task) shows entropy decline but is NOT division
        of labor — it lacks distinct functional roles. Requiring
        late_coverage ≥ 0.5 gates out collapse-to-one-task scenarios.
        Literature basis: Bonabeau et al. 1996 define division of labor as
        allocation of DIFFERENT tasks across individuals.
        """
        entropy_decline = primary_result.get("entropy_decline", 0.0)
        late_coverage = primary_result.get("late_coverage", 0.0)
        return entropy_decline > 0.1 and late_coverage >= 0.5

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Compute secondary metrics: role distinctness + efficiency."""
        bundle = extract_observation_bundle(state_history)
        assignments = bundle["task_assignments"]
        n_tasks = bundle["n_tasks"]
        n_agents = bundle["n_agents"]
        T = assignments.shape[0]

        # Per-agent dominant task in late window
        window = max(T // 4, 10)
        late_assignments = assignments[T - window:T]
        dominant_tasks = np.zeros(n_agents, dtype=np.int64)
        for i in range(n_agents):
            agent_tasks = late_assignments[:, i]
            active = agent_tasks[agent_tasks >= 0]
            if len(active) > 0:
                dominant_tasks[i] = int(np.bincount(active, minlength=n_tasks).argmax())
            else:
                dominant_tasks[i] = -1

        # Role diversity: how many distinct dominant tasks exist?
        unique_roles = set(dominant_tasks[dominant_tasks >= 0])
        role_diversity = len(unique_roles) / n_tasks if n_tasks > 0 else 0.0

        # Check for single-task collapse: all agents have same dominant task
        single_task_collapse = len(unique_roles) <= 1

        # Efficiency: fraction of late timesteps where all tasks are covered
        late_efficiency = 0
        for t in range(T - window, T):
            tasks_covered = set(assignments[t][assignments[t] >= 0])
            if len(tasks_covered) >= n_tasks:
                late_efficiency += 1
        late_efficiency = late_efficiency / window if window > 0 else 0.0

        # Baseline efficiency (early window)
        early_efficiency = 0
        for t in range(min(window, T)):
            tasks_covered = set(assignments[t][assignments[t] >= 0])
            if len(tasks_covered) >= n_tasks:
                early_efficiency += 1
        early_efficiency = early_efficiency / window if window > 0 else 0.0

        efficiency_gain = late_efficiency - early_efficiency

        return {
            "role_diversity": role_diversity,
            "single_task_collapse": single_task_collapse,
            "late_efficiency": late_efficiency,
            "early_efficiency": early_efficiency,
            "efficiency_gain": efficiency_gain,
            "n_unique_roles": len(unique_roles),
        }

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Null model: time-shuffle per-agent assignments.

        Shuffling timestamps within each agent's history destroys any
        temporal trend in specialization while preserving marginal
        task-frequency distributions.
        """
        bundle = extract_observation_bundle(state_history)
        assignments = bundle["task_assignments"]
        n_tasks = bundle["n_tasks"]
        T = assignments.shape[0]
        n_agents = assignments.shape[1]

        rng = np.random.default_rng(self._seed)
        window = max(T // 4, 10)
        observed_decline = primary_result.get("entropy_decline", 0.0)

        null_declines = []
        for _ in range(self.n_permutations):
            shuffled = assignments.copy()
            for i in range(n_agents):
                rng.shuffle(shuffled[:, i])

            early_ent = _per_agent_entropy_window(shuffled, n_tasks, 0, window)
            late_ent = _per_agent_entropy_window(shuffled, n_tasks, T - window, T)
            null_decline = float(np.mean(early_ent) - np.mean(late_ent))
            null_declines.append(null_decline)

        null_declines = np.array(null_declines)
        p_value = float(np.mean(null_declines >= observed_decline))
        if p_value == 0:
            p_value = 1.0 / (self.n_permutations + 1)

        return (
            p_value,
            NullType.SHUFFLE,
            {
                "mean": float(null_declines.mean()),
                "std": float(null_declines.std()),
            },
        )

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Confirmation: distinct roles + task coverage + null rejection.

        Requires:
        - Null p < 0.01
        - Role diversity: agents cover multiple tasks (not single-task collapse)
        - Late coverage of all tasks across the population
        - Low late per-agent entropy
        """
        if null_p >= 0.01:
            return False
        if secondary_result.get("single_task_collapse", True):
            return False
        if secondary_result.get("role_diversity", 0.0) < 0.5:
            return False
        if primary_result.get("late_coverage", 0.0) < 0.8:
            return False
        if primary_result.get("mean_late_entropy", 1.0) > 1.0:
            return False
        return True

    def _check_definitive(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        null_type: NullType,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> bool:
        """Definitive: efficiency gain + strong entropy decline + strong null.

        Requires:
        - Null p < 0.005
        - Efficiency gain > 0 (late > early or late efficiency > 0.5)
        - Entropy decline > 0.3 bits
        - Role diversity >= 0.67 (at least 2/3 of tasks covered by distinct agents)
        """
        if null_p >= 0.005:
            return False
        if primary_result.get("entropy_decline", 0.0) < 0.3:
            return False
        late_eff = secondary_result.get("late_efficiency", 0.0)
        eff_gain = secondary_result.get("efficiency_gain", 0.0)
        if late_eff < 0.3 and eff_gain <= 0.0:
            return False
        if secondary_result.get("role_diversity", 0.0) < 0.67:
            return False
        return True

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Exclude P23 (anti-coordination / continuous re-balancing).

        P32 specialization is role-FIXED: late switching frequency should be
        low. P23 anti-coordination is continuous re-balancing: agents keep
        changing choices. If late switching frequency is high (> 0.5),
        P23 is not excluded.
        """
        bundle = extract_observation_bundle(state_history)
        assignments = bundle["task_assignments"]
        T = assignments.shape[0]
        window = max(T // 4, 10)

        late_switch = _switching_frequency(assignments, T - window, T)

        checked = ["P23"]
        if late_switch < 0.3:
            # Low switching → agents are specialized, not re-balancing → exclude P23
            results = {"P23": "excluded"}
        else:
            results = {"P23": "inconclusive"}

        return checked, results

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        """All secondary metrics must pass."""
        return (
            not secondary_result.get("single_task_collapse", True)
            and secondary_result.get("role_diversity", 0.0) >= 0.5
            and secondary_result.get("late_efficiency", 0.0) > 0.0
        )

    def _check_finite_size_robustness(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> bool:
        """Finite-size check: run length >= 10τ."""
        return len(state_history) >= 10 * timescale
