"""Cell-view sorting model — faithful reimplementation of Zhang et al. (2024).

Reference: https://github.com/Zhangtaining/cell_research
Paper: "Classical sorting algorithms as a model of morphogenesis"
       Zhang, Goldstein & Levin, Adaptive Behavior 33, 25-54.

Key architectural choices matching the reference:
- 1D array (default 100 elements), not 2D grid
- Cell-view: each cell makes autonomous local decisions
- Concurrent execution simulated via randomized activation order
- Frozen cells: damaged elements that refuse to swap
- Chimeric arrays: cells with different algotypes mixed together
- DG metric: monotonicity-based (global sortedness backtracking), not distance-based

State history keys produced:
    array       : list[int] — current array state (values at positions)
    cell_types  : list[str] — algotype per position ('Bubble'/'Insertion'/'Selection')
    frozen_mask : list[bool] — which positions contain frozen cells
    step        : int — round counter
    swap_count  : int — total swaps so far
    monotonicity_error : int — count of adjacent inversions
    sortedness  : float — fraction of elements in correct sorted position
    n           : int — array size
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_model import BaseModel


class CellViewSorting(BaseModel):
    """Faithful reimplementation of Zhang et al. cell-view sorting.

    Parameters
    ----------
    n : int
        Array size (default 100, matching paper).
    algorithm : str
        'bubble', 'insertion', 'selection', or 'chimera'.
    chimera_mix : dict, optional
        For chimeric arrays: {'Bubble': 0.5, 'Selection': 0.5} etc.
        Fractions of each algotype. Only used when algorithm='chimera'.
    n_frozen : int
        Number of frozen (damaged) cells.
    frozen_type : str
        'movable' (others can push it) or 'immovable'.
    seed : int
        Random seed.
    """

    def __init__(
        self,
        n: int = 100,
        algorithm: str = "bubble",
        chimera_mix: dict[str, float] | None = None,
        n_frozen: int = 0,
        frozen_type: str = "movable",
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)
        self.n = n
        self.algorithm = algorithm
        self.chimera_mix = chimera_mix
        self.n_frozen = n_frozen
        self.frozen_type = frozen_type

        # Cell state
        self._values: np.ndarray | None = None
        self._cell_types: list[str] | None = None
        self._frozen: np.ndarray | None = None
        self._ideal_pos: np.ndarray | None = None

        self._converged = False
        self._swap_count = 0

    def setup(self) -> dict[str, Any]:
        values = list(range(self.n))
        self.rng.shuffle(values)
        self._values = np.array(values, dtype=int)

        if self.algorithm == "chimera" and self.chimera_mix:
            self._cell_types = self._assign_chimera_types()
        else:
            algo_name = {"bubble": "Bubble", "insertion": "Insertion",
                         "selection": "Selection"}[self.algorithm]
            self._cell_types = [algo_name] * self.n

        self._frozen = np.zeros(self.n, dtype=bool)
        if self.n_frozen > 0:
            frozen_indices = self.rng.choice(self.n, size=self.n_frozen, replace=False)
            self._frozen[frozen_indices] = True

        self._ideal_pos = np.zeros(self.n, dtype=int)

        self._step_count = 0
        self._swap_count = 0
        self._converged = False
        self._is_setup = True

        # Per-swap traces for DG and aggregation (matching Zhang's per-swap recording)
        self._mono_error_trace: list[int] = [get_monotonicity_error(self._values.tolist())]
        self._aggregation_trace: list[float] = [get_aggregation_value(self._cell_types)]

        return self._snapshot()

    def step(self) -> dict[str, Any]:
        """One round: activate all non-frozen cells in random order."""
        if self._converged:
            return self._snapshot()

        active_indices = np.where(~self._frozen)[0]
        self.rng.shuffle(active_indices)

        any_swap = False
        for i in active_indices:
            swapped = self._cell_action(int(i))
            if swapped:
                any_swap = True
                self._swap_count += 1

        self._step_count += 1

        if not any_swap and self._no_cell_should_move():
            self._converged = True

        return self._snapshot()

    def get_metadata(self) -> dict[str, Any]:
        return {
            "model_name": "zhang_cell_view_sorting",
            "model_class": "sorting",
            "algorithm": self.algorithm,
            "n": self.n,
            "n_frozen": self.n_frozen,
            "frozen_type": self.frozen_type,
            "chimera_mix": self.chimera_mix,
            "interaction_type": "direct",
            "update_mode": "motion_transport",
            "seed": self.seed,
            "reference": "Zhang, Goldstein & Levin (2024)",
        }

    def is_converged(self) -> bool:
        return self._converged

    def get_timescale(self) -> float:
        if self._converged:
            return float(self._step_count)
        return float(self.n * self.n)

    def run_to_completion(self, max_rounds: int | None = None) -> list[dict[str, Any]]:
        if not self._is_setup:
            initial = self.setup()
            self._state_history = [initial]

        if max_rounds is None:
            max_rounds = self.n * self.n * 2

        for _ in range(max_rounds):
            state = self.step()
            self._state_history.append(state)
            if self._converged:
                break

        return self._state_history

    # --- Cell actions (matching Zhang's cell implementations) ---

    def _cell_action(self, i: int) -> bool:
        cell_type = self._cell_types[i]
        if cell_type == "Bubble":
            return self._bubble_action(i)
        elif cell_type == "Insertion":
            return self._insertion_action(i)
        elif cell_type == "Selection":
            return self._selection_action(i)
        return False

    def _bubble_action(self, i: int) -> bool:
        """Bubble cell: randomly check left or right, swap if out of order."""
        check_right = self.rng.random() < 0.5
        if check_right:
            target = i + 1
            if target >= self.n:
                return False
            if self._frozen[target] and self.frozen_type == "immovable":
                return False
            if self._values[i] > self._values[target]:
                return self._do_swap(i, target)
        else:
            target = i - 1
            if target < 0:
                return False
            if self._frozen[target] and self.frozen_type == "immovable":
                return False
            if self._values[i] < self._values[target]:
                return self._do_swap(i, target)
        return False

    def _insertion_action(self, i: int) -> bool:
        """Insertion cell: swap left if left portion is sorted and left neighbor is larger."""
        if i == 0:
            return False
        if not self._is_left_sorted(i):
            return False
        target = i - 1
        if self._frozen[target] and self.frozen_type == "immovable":
            return False
        if self._values[i] < self._values[target]:
            return self._do_swap(i, target)
        return False

    def _selection_action(self, i: int) -> bool:
        """Selection cell: try to swap with ideal position (leftmost unsorted)."""
        ideal = self._ideal_pos[i]
        if ideal >= self.n or ideal == i:
            return False
        if self._frozen[ideal]:
            self._ideal_pos[i] = min(ideal + 1, self.n - 1)
            return False
        if self._values[i] < self._values[ideal]:
            return self._do_swap(i, ideal)
        else:
            self._ideal_pos[i] = min(ideal + 1, self.n - 1)
            return False

    def _do_swap(self, a: int, b: int) -> bool:
        if self._frozen[a]:
            return False
        self._values[a], self._values[b] = self._values[b], self._values[a]
        self._cell_types[a], self._cell_types[b] = self._cell_types[b], self._cell_types[a]
        self._frozen[a], self._frozen[b] = self._frozen[b], self._frozen[a]
        self._ideal_pos[a], self._ideal_pos[b] = self._ideal_pos[b], self._ideal_pos[a]
        # Record per-swap traces (matching Zhang's StatusProbe.record_sorting_step)
        self._mono_error_trace.append(get_monotonicity_error(self._values.tolist()))
        self._aggregation_trace.append(get_aggregation_value(self._cell_types))
        return True

    def _is_left_sorted(self, pos: int) -> bool:
        prev = -1
        for j in range(0, pos):
            if self._frozen[j]:
                prev = -1
                continue
            if self._values[j] < prev:
                return False
            prev = self._values[j]
        return True

    def _no_cell_should_move(self) -> bool:
        for i in range(self.n):
            if self._frozen[i]:
                continue
            ct = self._cell_types[i]
            if ct == "Bubble":
                if i > 0 and self._values[i] < self._values[i - 1]:
                    return False
                if i < self.n - 1 and self._values[i] > self._values[i + 1]:
                    return False
            elif ct == "Insertion":
                if i > 0 and self._is_left_sorted(i) and self._values[i] < self._values[i - 1]:
                    return False
            elif ct == "Selection":
                ideal = self._ideal_pos[i]
                if ideal < self.n and ideal != i:
                    if not self._frozen[ideal] and self._values[i] < self._values[ideal]:
                        return False
        return True

    def _assign_chimera_types(self) -> list[str]:
        types = []
        type_names = list(self.chimera_mix.keys())
        fractions = list(self.chimera_mix.values())
        cum_fracs = np.cumsum(fractions)
        for _ in range(self.n):
            r = self.rng.random()
            assigned = type_names[-1]
            for j, cf in enumerate(cum_fracs):
                if r < cf:
                    assigned = type_names[j]
                    break
            types.append(assigned)
        return types

    def _snapshot(self) -> dict[str, Any]:
        arr = self._values.tolist()
        return {
            "array": arr,
            "cell_types": list(self._cell_types),
            "frozen_mask": self._frozen.tolist(),
            "step": self._step_count,
            "swap_count": self._swap_count,
            "monotonicity_error": get_monotonicity_error(arr),
            "sortedness": get_sortedness(arr),
            "n": self.n,
        }

    @property
    def mono_error_trace(self) -> list[int]:
        """Per-swap monotonicity error trace for DG computation."""
        return self._mono_error_trace

    @property
    def aggregation_trace(self) -> list[float]:
        """Per-swap aggregation value trace for chimeric analysis."""
        return self._aggregation_trace

    def compute_dg(self) -> float:
        """Compute delayed gratification from per-swap monotonicity trace."""
        return compute_delayed_gratification(self._mono_error_trace)


# --- Zhang's exact metrics ---

def get_monotonicity_error(arr: list[int]) -> int:
    """Count adjacent inversions. Matches Zhang's get_monotonicity()."""
    errors = 0
    for i in range(1, len(arr)):
        if arr[i] < arr[i - 1]:
            errors += 1
    return errors


def get_sortedness(arr: list[int]) -> float:
    """Fraction of elements in their correct sorted position."""
    expected = sorted(arr)
    correct = sum(1 for a, b in zip(arr, expected) if a == b)
    return correct / len(arr) if arr else 0.0


def get_aggregation_value(cell_types: list[str]) -> float:
    """Fraction of neighbor pairs sharing the same algotype.

    Matches Zhang's aggregation metric for chimeric arrays.
    """
    if len(cell_types) <= 1:
        return 1.0
    same = sum(1 for i in range(len(cell_types) - 1) if cell_types[i] == cell_types[i + 1])
    return same / (len(cell_types) - 1)


def compute_delayed_gratification(monotonicity_errors: list[int]) -> float:
    """Compute DG from monotonicity error trajectory.

    Direct port of Zhang's avg_wandering_range() + get_discrepency_arr().
    Reference: cell_research/analysis/delay_gratification_analysis.py

    IMPORTANT: This function has an early return that returns the accumulated
    sum (not average) when it reaches the end of the runs or hits two
    consecutive positive runs. This matches Zhang's code exactly.
    """
    if len(monotonicity_errors) < 2:
        return 0.0

    # --- dedup: remove consecutive duplicates ---
    deduped = [monotonicity_errors[0]]
    for v in monotonicity_errors[1:]:
        if v != deduped[-1]:
            deduped.append(v)

    if len(deduped) < 2:
        return 0.0

    # --- get_discrepency_arr: diffs then group same-sign runs ---
    temp = [deduped[i] - deduped[i - 1] for i in range(1, len(deduped))]

    if not temp:
        return 0.0

    runs = [temp[0]]
    for i in range(1, len(temp)):
        if temp[i - 1] * temp[i] > 0:
            runs[-1] += temp[i]
        else:
            runs.append(temp[i])

    # --- get_first_pos: first positive run index ---
    first_pos = len(runs)
    for idx in range(len(runs)):
        if runs[idx] > 0:
            first_pos = idx
            break

    if first_pos >= len(runs) - 1:
        return 0.0

    # --- avg_wandering_range: walk paired runs ---
    res = 0.0
    n = 0
    i = first_pos

    while i < len(runs):
        j = i + 1
        # Zhang's early return: if j out of bounds or next run is positive
        if not (j < len(runs)) or runs[j] > 0:
            return res  # accumulated SUM, not average — matches Zhang exactly

        change = (-(runs[j] + runs[i])) / runs[i]
        res += change
        n += 1
        i += 2

    if n == 0:
        return 0.0
    return res / n  # average — only reached if all runs are perfectly paired
