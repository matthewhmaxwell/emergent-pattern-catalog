"""Response-threshold model for emergent division of labor.

Implements the Bonabeau, Theraulaz & Deneubourg (1996) response-threshold
model: identical agents each have per-task response thresholds. Task demand
generates stimulus; an agent performs a task when its stimulus exceeds its
threshold. Thresholds reinforce with performance (performing a task lowers
the agent's threshold for that task). Over time, agents differentiate into
consistent functional roles.

Reference:
  Bonabeau, E., Theraulaz, G. & Deneubourg, J.-L. (1996).
  Quantitative study of the fixed threshold model for the regulation
  of division of labour in insect societies.
  Proceedings of the Royal Society B, 263(1376), 1565–1569.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_model import BaseModel


class ResponseThresholdModel(BaseModel):
    """Response-threshold model for emergent specialization.

    Parameters
    ----------
    n_agents : int
        Number of identical agents.
    n_tasks : int
        Number of task types.
    reinforcement_rate : float
        Rate at which thresholds decrease for performed tasks (ξ).
        Higher → faster specialization.
    forgetting_rate : float
        Rate at which thresholds increase for unperformed tasks (φ).
    stimulus_rate : float
        Base rate at which task stimulus accumulates per timestep.
    initial_threshold : float
        Starting threshold for all agents on all tasks.
    seed : int
        Random seed.
    """

    def __init__(
        self,
        n_agents: int = 20,
        n_tasks: int = 3,
        reinforcement_rate: float = 0.05,
        forgetting_rate: float = 0.01,
        stimulus_rate: float = 0.1,
        initial_threshold: float = 0.5,
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)
        self.n_agents = n_agents
        self.n_tasks = n_tasks
        self.reinforcement_rate = reinforcement_rate
        self.forgetting_rate = forgetting_rate
        self.stimulus_rate = stimulus_rate
        self.initial_threshold = initial_threshold

        # State arrays (initialized in setup)
        self._thresholds: np.ndarray | None = None
        self._stimulus: np.ndarray | None = None
        self._task_assignments: np.ndarray | None = None  # last assignment per agent

    def setup(self) -> dict[str, Any]:
        """Initialize all agents with identical thresholds."""
        self._thresholds = np.full(
            (self.n_agents, self.n_tasks),
            self.initial_threshold,
            dtype=np.float64,
        )
        self._stimulus = np.full(self.n_tasks, 0.5, dtype=np.float64)
        self._task_assignments = np.full(self.n_agents, -1, dtype=np.int64)
        self._step_count = 0
        self._is_setup = True
        return self._snapshot()

    def step(self) -> dict[str, Any]:
        """Advance one timestep.

        Each agent independently evaluates all tasks. The probability of
        responding to task j is  s_j^2 / (s_j^2 + θ_{i,j}^2)  where s_j
        is the stimulus for task j and θ_{i,j} is agent i's threshold for
        task j (Bonabeau et al. 1996, Eq. 1).

        If an agent responds to a task:
          - Its threshold for that task DECREASES by reinforcement_rate
        If an agent does NOT respond to a task:
          - Its threshold for that task INCREASES by forgetting_rate

        Stimulus accumulates at stimulus_rate per timestep for each task,
        reduced by the number of agents performing it.
        """
        assert self._thresholds is not None

        # Compute response probabilities: s^2 / (s^2 + θ^2)
        s2 = self._stimulus ** 2  # (n_tasks,)
        th2 = self._thresholds ** 2  # (n_agents, n_tasks)
        denom = s2[np.newaxis, :] + th2
        denom = np.where(denom > 0, denom, 1e-12)
        prob = s2[np.newaxis, :] / denom  # (n_agents, n_tasks)

        # Each agent decides independently for each task
        rand = self.rng.random((self.n_agents, self.n_tasks))
        responds = rand < prob  # boolean (n_agents, n_tasks)

        # Assign each agent to at most one task (highest-probability task
        # among those it responds to; if none, idle)
        assignments = np.full(self.n_agents, -1, dtype=np.int64)
        for i in range(self.n_agents):
            active_tasks = np.where(responds[i])[0]
            if len(active_tasks) > 0:
                # Pick the task with highest response probability
                best = active_tasks[np.argmax(prob[i, active_tasks])]
                assignments[i] = best

        self._task_assignments = assignments

        # Threshold reinforcement / forgetting
        for i in range(self.n_agents):
            if assignments[i] >= 0:
                task = assignments[i]
                # Decrease threshold for performed task
                self._thresholds[i, task] = max(
                    0.01,
                    self._thresholds[i, task] - self.reinforcement_rate,
                )
                # Increase threshold for all other tasks
                for j in range(self.n_tasks):
                    if j != task:
                        self._thresholds[i, j] = min(
                            1.0,
                            self._thresholds[i, j] + self.forgetting_rate,
                        )
            else:
                # Idle agent: all thresholds drift up
                self._thresholds[i, :] = np.minimum(
                    1.0,
                    self._thresholds[i, :] + self.forgetting_rate,
                )

        # Stimulus dynamics: accumulate, reduce by workers
        workers_per_task = np.zeros(self.n_tasks, dtype=np.float64)
        for i in range(self.n_agents):
            if assignments[i] >= 0:
                workers_per_task[assignments[i]] += 1.0

        self._stimulus = np.clip(
            self._stimulus + self.stimulus_rate - 0.1 * workers_per_task,
            0.0,
            1.0,
        )

        self._step_count += 1
        return self._snapshot()

    def _snapshot(self) -> dict[str, Any]:
        """Return current state as a dict."""
        assert self._thresholds is not None
        assert self._task_assignments is not None

        # Per-agent task entropy (over assignment history is computed externally;
        # here we return instantaneous data)
        return {
            "step": self._step_count,
            "task_assignments": self._task_assignments.copy(),
            "thresholds": self._thresholds.copy(),
            "stimulus": self._stimulus.copy(),
            "n_agents": self.n_agents,
            "n_tasks": self.n_tasks,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata."""
        return {
            "model_name": "division_of_labor",
            "model_class": "task_allocation",
            "n_agents": self.n_agents,
            "n_tasks": self.n_tasks,
            "reinforcement_rate": self.reinforcement_rate,
            "forgetting_rate": self.forgetting_rate,
            "stimulus_rate": self.stimulus_rate,
            "service_rate": 0.1,   # per-worker demand reduction (stimulus dynamics)
            "stimulus_max": 1.0,   # per-task stimulus clip ceiling
            "initial_threshold": self.initial_threshold,
            "interaction_type": "threshold_response",
            "update_mode": "synchronous",
            "has_threshold_reinforcement": True,
            "has_forgetting": True,
            "seed": self.seed,
            "reference": "Bonabeau, Theraulaz & Deneubourg 1996",
        }

    def get_timescale(self) -> float:
        """Characteristic timescale: 1 / reinforcement_rate."""
        return 1.0 / max(self.reinforcement_rate, 1e-6)


class NoReinforcementModel(BaseModel):
    """Control model: identical response-threshold model without reinforcement.

    Agents respond to stimuli but thresholds never change — no specialization
    emerges. Used as a negative control.
    """

    def __init__(
        self,
        n_agents: int = 20,
        n_tasks: int = 3,
        stimulus_rate: float = 0.1,
        initial_threshold: float = 0.5,
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)
        self.n_agents = n_agents
        self.n_tasks = n_tasks
        self.stimulus_rate = stimulus_rate
        self.initial_threshold = initial_threshold

        self._thresholds: np.ndarray | None = None
        self._stimulus: np.ndarray | None = None
        self._task_assignments: np.ndarray | None = None

    def setup(self) -> dict[str, Any]:
        """Initialize identical agents."""
        self._thresholds = np.full(
            (self.n_agents, self.n_tasks),
            self.initial_threshold,
            dtype=np.float64,
        )
        self._stimulus = np.full(self.n_tasks, 0.5, dtype=np.float64)
        self._task_assignments = np.full(self.n_agents, -1, dtype=np.int64)
        self._step_count = 0
        self._is_setup = True
        return self._snapshot()

    def step(self) -> dict[str, Any]:
        """Step without any threshold reinforcement."""
        assert self._thresholds is not None

        s2 = self._stimulus ** 2
        th2 = self._thresholds ** 2
        denom = s2[np.newaxis, :] + th2
        denom = np.where(denom > 0, denom, 1e-12)
        prob = s2[np.newaxis, :] / denom

        rand = self.rng.random((self.n_agents, self.n_tasks))
        responds = rand < prob

        assignments = np.full(self.n_agents, -1, dtype=np.int64)
        for i in range(self.n_agents):
            active_tasks = np.where(responds[i])[0]
            if len(active_tasks) > 0:
                best = active_tasks[np.argmax(prob[i, active_tasks])]
                assignments[i] = best

        self._task_assignments = assignments

        # No threshold update — this is the control

        # Stimulus dynamics
        workers_per_task = np.zeros(self.n_tasks, dtype=np.float64)
        for i in range(self.n_agents):
            if assignments[i] >= 0:
                workers_per_task[assignments[i]] += 1.0

        self._stimulus = np.clip(
            self._stimulus + self.stimulus_rate - 0.1 * workers_per_task,
            0.0,
            1.0,
        )

        self._step_count += 1
        return self._snapshot()

    def _snapshot(self) -> dict[str, Any]:
        """Return current state."""
        assert self._thresholds is not None
        assert self._task_assignments is not None
        return {
            "step": self._step_count,
            "task_assignments": self._task_assignments.copy(),
            "thresholds": self._thresholds.copy(),
            "stimulus": self._stimulus.copy(),
            "n_agents": self.n_agents,
            "n_tasks": self.n_tasks,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata."""
        return {
            "model_name": "division_of_labor_no_reinforcement",
            "model_class": "task_allocation",
            "n_agents": self.n_agents,
            "n_tasks": self.n_tasks,
            "reinforcement_rate": 0.0,
            "forgetting_rate": 0.0,
            "stimulus_rate": self.stimulus_rate,
            "initial_threshold": self.initial_threshold,
            "interaction_type": "threshold_response",
            "update_mode": "synchronous",
            "has_threshold_reinforcement": False,
            "has_forgetting": False,
            "seed": self.seed,
            "reference": "Bonabeau, Theraulaz & Deneubourg 1996 (control)",
        }

    def get_timescale(self) -> float:
        """No characteristic timescale for the control."""
        return float(self.n_agents)


def compute_per_agent_entropy(history: list[dict[str, Any]], n_tasks: int) -> np.ndarray:
    """Compute per-agent task entropy from a state history.

    Parameters
    ----------
    history : list[dict]
        State history from ResponseThresholdModel or NoReinforcementModel.
    n_tasks : int
        Number of task types.

    Returns
    -------
    ndarray(n_agents,)
        Shannon entropy of each agent's task distribution over the run.
    """
    if not history:
        return np.array([])

    n_agents = history[0]["n_agents"]
    counts = np.zeros((n_agents, n_tasks), dtype=np.float64)

    for state in history:
        assignments = state["task_assignments"]
        for i in range(n_agents):
            if assignments[i] >= 0:
                counts[i, assignments[i]] += 1.0

    # Normalize to probabilities
    totals = counts.sum(axis=1, keepdims=True)
    totals = np.where(totals > 0, totals, 1.0)
    probs = counts / totals

    # Shannon entropy: -sum(p * log(p))
    log_probs = np.where(probs > 0, np.log2(probs), 0.0)
    entropy = -np.sum(probs * log_probs, axis=1)

    return entropy


def compute_windowed_entropy(
    history: list[dict[str, Any]],
    n_tasks: int,
    window_size: int = 50,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-agent entropy in sliding windows.

    Returns
    -------
    (early_entropy, late_entropy) : tuple of ndarray(n_agents,)
        Mean per-agent entropy in the first and last windows.
    """
    if len(history) < 2 * window_size:
        # Fall back to halves
        mid = len(history) // 2
        early = compute_per_agent_entropy(history[:mid], n_tasks)
        late = compute_per_agent_entropy(history[mid:], n_tasks)
        return early, late

    early = compute_per_agent_entropy(history[:window_size], n_tasks)
    late = compute_per_agent_entropy(history[-window_size:], n_tasks)
    return early, late


def compute_efficiency(
    history: list[dict[str, Any]],
    n_tasks: int,
) -> float:
    """Compute collective efficiency: fraction of timesteps where all tasks
    have at least one worker assigned.

    Parameters
    ----------
    history : list[dict]
        State history.
    n_tasks : int
        Number of task types.

    Returns
    -------
    float
        Fraction of timesteps with full task coverage.
    """
    if not history:
        return 0.0

    covered_steps = 0
    for state in history:
        assignments = state["task_assignments"]
        tasks_covered = set(assignments[assignments >= 0])
        if len(tasks_covered) >= n_tasks:
            covered_steps += 1

    return covered_steps / len(history)
