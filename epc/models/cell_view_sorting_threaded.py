"""Cell-view sorting — threaded implementation matching Zhang et al. exactly.

This implements Zhang's actual architecture: each cell is a Python thread
running a tight while-loop, acquiring a shared lock to check neighbors and
swap. The OS thread scheduler determines which cell acts next, creating
genuine contention that produces the context-sensitive DG behavior around
frozen cells.

Reference: https://github.com/Zhangtaining/cell_research
Paper: Zhang, Goldstein & Levin (2024), Adaptive Behavior 33, 25-54.

NOTE: This is nondeterministic by design (thread scheduling). Each run
of the same initial condition may produce different results. Run multiple
trials and report statistics, matching Zhang's methodology (100 trials).
"""

from __future__ import annotations

import random
import threading
import time
from enum import Enum
from typing import Any

from epc.base_model import BaseModel


class CellStatus(Enum):
    ACTIVE = 1
    SLEEP = 2
    MOVING = 4
    INACTIVE = 5
    FREEZE = 7


class StatusProbe:
    """Records per-swap state snapshots, matching Zhang's StatusProbe."""

    def __init__(self) -> None:
        self.sorting_steps: list[list[int]] = []
        self.cell_type_steps: list[list[str]] = []
        self.swap_count: int = 0
        self.compare_and_swap_count: int = 0
        self.frozen_swap_attempts: int = 0

    def record_swap(self, cells: list[ThreadedCell]) -> None:
        self.swap_count += 1
        self.sorting_steps.append([c.value for c in cells])
        self.cell_type_steps.append([c.cell_type for c in cells])

    def record_compare(self) -> None:
        self.compare_and_swap_count += 1

    def count_frozen_attempt(self) -> None:
        self.frozen_swap_attempts += 1


class ThreadedCell(threading.Thread):
    """Base threaded cell. Each cell runs as its own thread."""

    def __init__(
        self,
        value: int,
        position: int,
        cells: list[ThreadedCell],
        lock: threading.Lock,
        probe: StatusProbe,
        cell_type: str = "Bubble",
    ) -> None:
        threading.Thread.__init__(self, daemon=True)
        self.value = value
        self.current_position = position
        self.cells = cells
        self.lock = lock
        self.probe = probe
        self.cell_type = cell_type
        self.status = CellStatus.ACTIVE
        self.previous_status = CellStatus.ACTIVE
        self.tried_to_swap_with_frozen = False

    def set_freeze(self) -> None:
        self.status = CellStatus.FREEZE
        self.previous_status = CellStatus.FREEZE

    def swap(self, target_pos: int) -> None:
        """Swap this cell with the cell at target_pos."""
        if self.status == CellStatus.FREEZE:
            if not self.tried_to_swap_with_frozen:
                self.probe.count_frozen_attempt()
                self.tried_to_swap_with_frozen = True
            return

        target_cell = self.cells[target_pos]
        self.tried_to_swap_with_frozen = False
        target_cell.tried_to_swap_with_frozen = False

        self.status = CellStatus.MOVING
        target_cell.status = CellStatus.MOVING

        # Swap in the shared cells list
        self.cells[self.current_position] = target_cell
        self.cells[target_pos] = self

        # Update positions
        target_cell.current_position = self.current_position
        self.current_position = target_pos

        # Immediate completion (disable_visualization=True in Zhang)
        self.status = self.previous_status
        target_cell.status = target_cell.previous_status

        self.probe.record_swap(self.cells)

    def should_move(self) -> bool:
        raise NotImplementedError

    def move(self) -> None:
        raise NotImplementedError

    def run(self) -> None:
        while self.status != CellStatus.INACTIVE:
            self.move()


class BubbleSortCell(ThreadedCell):
    """Bubble sort cell — randomly checks left or right neighbor."""

    def __init__(self, value, position, cells, lock, probe):
        super().__init__(value, position, cells, lock, probe, "Bubble")

    def should_move(self) -> bool:
        pos = self.current_position
        n = len(self.cells)
        smaller_left = False
        if pos > 0:
            left = self.cells[pos - 1]
            smaller_left = (self.value < left.value and left.status == CellStatus.ACTIVE)
        bigger_right = False
        if pos < n - 1:
            right = self.cells[pos + 1]
            bigger_right = (self.value > right.value and right.status == CellStatus.ACTIVE)
        return smaller_left or bigger_right

    def move(self) -> None:
        self.lock.acquire()
        try:
            if self.should_move():
                self.probe.record_compare()

            pos = self.current_position
            n = len(self.cells)
            check_right = random.random() < 0.5

            if check_right:
                target = pos + 1
                if target < n and self.status == CellStatus.ACTIVE:
                    target_cell = self.cells[target]
                    if (target_cell.status == CellStatus.ACTIVE or
                            target_cell.status == CellStatus.FREEZE):
                        if self.value > target_cell.value:
                            self.swap(target)
            else:
                target = pos - 1
                if target >= 0 and self.status == CellStatus.ACTIVE:
                    target_cell = self.cells[target]
                    if (target_cell.status == CellStatus.ACTIVE or
                            target_cell.status == CellStatus.FREEZE):
                        if self.value < target_cell.value:
                            self.swap(target)
        finally:
            self.lock.release()


class InsertionSortCell(ThreadedCell):
    """Insertion sort cell — moves left if left portion is sorted."""

    def __init__(self, value, position, cells, lock, probe):
        super().__init__(value, position, cells, lock, probe, "Insertion")

    def _is_left_sorted(self) -> bool:
        prev = -1
        for j in range(0, self.current_position):
            c = self.cells[j]
            if c.status == CellStatus.FREEZE:
                prev = -1
                continue
            if c.value < prev:
                return False
            prev = c.value
        return True

    def should_move(self) -> bool:
        pos = self.current_position
        if pos == 0:
            return False
        if not self._is_left_sorted():
            return False
        left = self.cells[pos - 1]
        return (self.value < left.value and left.status == CellStatus.ACTIVE)

    def move(self) -> None:
        self.lock.acquire()
        try:
            if not self._is_left_sorted():
                return

            if self.should_move():
                self.probe.record_compare()

            pos = self.current_position
            target = pos - 1
            if target >= 0 and self.status == CellStatus.ACTIVE:
                target_cell = self.cells[target]
                if (target_cell.status == CellStatus.ACTIVE or
                        target_cell.status == CellStatus.FREEZE):
                    if self.value < target_cell.value:
                        self.swap(target)
        finally:
            self.lock.release()


class SelectionSortCell(ThreadedCell):
    """Selection sort cell — swaps with ideal position (leftmost unsorted)."""

    def __init__(self, value, position, cells, lock, probe):
        super().__init__(value, position, cells, lock, probe, "Selection")
        self.ideal_position = 0

    def should_move(self) -> bool:
        return (self.current_position != self.ideal_position and
                0 <= self.ideal_position < len(self.cells))

    def move(self) -> None:
        self.lock.acquire()
        try:
            if self.should_move():
                self.probe.record_compare()

            ideal = self.ideal_position
            n = len(self.cells)
            if ideal < 0 or ideal >= n or ideal == self.current_position:
                return
            if self.status != CellStatus.ACTIVE:
                return

            target_cell = self.cells[ideal]

            if target_cell.status == CellStatus.FREEZE:
                self.ideal_position = min(ideal + 1, n - 1)
                if self.value < target_cell.value:
                    self.swap(ideal)  # counts as frozen attempt
                return

            if target_cell.status == CellStatus.ACTIVE:
                if self.value < target_cell.value:
                    self.swap(ideal)
                else:
                    self.ideal_position = min(ideal + 1, n - 1)
        finally:
            self.lock.release()


# --- Experiment runner ---

_CELL_CLASSES = {
    "bubble": BubbleSortCell,
    "insertion": InsertionSortCell,
    "selection": SelectionSortCell,
}


class GroupMonitor(threading.Thread):
    """Matches Zhang's CellGroup thread behavior.

    Periodically acquires the shared lock to check if the array is sorted.
    This creates additional lock contention that affects cell timing patterns.
    Zhang's CellGroup does this every ~50ms with time.sleep(0.05).
    """

    def __init__(self, cells: list[ThreadedCell], lock: threading.Lock) -> None:
        threading.Thread.__init__(self, daemon=True)
        self.cells = cells
        self.lock = lock
        self.active = True

    def run(self) -> None:
        while self.active:
            self.lock.acquire()
            try:
                # Check if sorted (matching CellGroup.is_group_sorted)
                sorted_check = True
                prev = self.cells[0]
                for c in self.cells[1:]:
                    if c.status == CellStatus.MOVING or c.value < prev.value:
                        sorted_check = False
                        break
                    prev = c
            finally:
                self.lock.release()
            time.sleep(0.05)  # Matching Zhang's CellGroup sleep


def _no_cells_should_move(cells: list[ThreadedCell]) -> bool:
    for c in cells:
        if c.status == CellStatus.SLEEP:
            return False
        if c.status == CellStatus.ACTIVE and c.should_move():
            return False
    return True


def run_threaded_experiment(
    n: int = 100,
    algorithm: str = "bubble",
    n_frozen: int = 0,
    poll_interval: float = 0.5,
    timeout: float = 60.0,
) -> dict[str, Any]:
    """Run one threaded sorting experiment matching Zhang's methodology.

    Parameters
    ----------
    n : int
        Array size.
    algorithm : str
        'bubble', 'insertion', or 'selection'.
    n_frozen : int
        Number of frozen cells.
    poll_interval : float
        Seconds between convergence checks.
    timeout : float
        Maximum wall-clock seconds before aborting.

    Returns
    -------
    dict with:
        sorting_steps : list[list[int]] — per-swap array snapshots
        cell_type_steps : list[list[str]] — per-swap algotype snapshots
        swap_count : int
        compare_count : int
        frozen_attempts : int
        converged : bool
        final_array : list[int]
    """
    values = list(range(n))
    random.shuffle(values)

    lock = threading.Lock()
    probe = StatusProbe()

    CellClass = _CELL_CLASSES[algorithm]
    cells: list[ThreadedCell] = []
    for i in range(n):
        cell = CellClass(values[i], i, cells, lock, probe)
        cells.append(cell)

    # Freeze random cells
    if n_frozen > 0:
        frozen_indices = random.sample(range(n), n_frozen)
        for idx in frozen_indices:
            cells[idx].set_freeze()

    # Start all non-frozen cells
    lock.acquire()
    for cell in cells:
        if cell.status != CellStatus.FREEZE:
            cell.start()
    # Start group monitor (matches Zhang's CellGroup thread)
    monitor = GroupMonitor(cells, lock)
    monitor.start()
    lock.release()

    # Poll for convergence
    start_time = time.time()
    converged = False
    while time.time() - start_time < timeout:
        time.sleep(poll_interval)
        if _no_cells_should_move(cells):
            converged = True
            break

    # Kill all threads
    lock.acquire()
    for c in cells:
        c.status = CellStatus.INACTIVE
    monitor.active = False
    lock.release()

    # Wait for threads to finish
    for c in cells:
        if c.is_alive():
            c.join(timeout=2.0)

    return {
        "sorting_steps": probe.sorting_steps,
        "cell_type_steps": probe.cell_type_steps,
        "swap_count": probe.swap_count,
        "compare_count": probe.compare_and_swap_count,
        "frozen_attempts": probe.frozen_swap_attempts,
        "converged": converged,
        "final_array": [c.value for c in cells],
    }


def run_chimeric_experiment(
    n: int = 100,
    mix: dict[str, float] | None = None,
    n_frozen: int = 0,
    poll_interval: float = 0.5,
    timeout: float = 60.0,
) -> dict[str, Any]:
    """Run one chimeric (mixed algotype) threaded experiment.

    Parameters
    ----------
    mix : dict
        e.g. {'bubble': 0.5, 'selection': 0.5}
    """
    if mix is None:
        mix = {"bubble": 0.5, "selection": 0.5}

    values = list(range(n))
    random.shuffle(values)

    lock = threading.Lock()
    probe = StatusProbe()

    # Assign algotypes
    type_names = list(mix.keys())
    fractions = list(mix.values())
    cum_fracs = []
    s = 0
    for f in fractions:
        s += f
        cum_fracs.append(s)

    cells: list[ThreadedCell] = []
    for i in range(n):
        r = random.random()
        algo = type_names[-1]
        for j, cf in enumerate(cum_fracs):
            if r < cf:
                algo = type_names[j]
                break
        CellClass = _CELL_CLASSES[algo]
        cell = CellClass(values[i], i, cells, lock, probe)
        cells.append(cell)

    if n_frozen > 0:
        frozen_indices = random.sample(range(n), n_frozen)
        for idx in frozen_indices:
            cells[idx].set_freeze()

    lock.acquire()
    for cell in cells:
        if cell.status != CellStatus.FREEZE:
            cell.start()
    lock.release()

    start_time = time.time()
    converged = False
    while time.time() - start_time < timeout:
        time.sleep(poll_interval)
        if _no_cells_should_move(cells):
            converged = True
            break

    lock.acquire()
    for c in cells:
        c.status = CellStatus.INACTIVE
    lock.release()

    for c in cells:
        if c.is_alive():
            c.join(timeout=2.0)

    return {
        "sorting_steps": probe.sorting_steps,
        "cell_type_steps": probe.cell_type_steps,
        "swap_count": probe.swap_count,
        "compare_count": probe.compare_and_swap_count,
        "frozen_attempts": probe.frozen_swap_attempts,
        "converged": converged,
        "final_array": [c.value for c in cells],
        "final_cell_types": [c.cell_type for c in cells],
    }
