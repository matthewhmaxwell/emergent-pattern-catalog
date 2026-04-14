"""P13 — Excitable wave metrics.

Computes wavefront detection, propagation speed, spiral tip tracking,
and wave source persistence for excitable media.

Designed for Greenberg-Hastings and similar excitable CAs with discrete
states: resting (0), excited (1), refractory (2+).

Supported state history format:
    grid      : np.ndarray (rows, cols) — cell states
    grid_dims : (rows, cols)
    n_states  : int — κ (number of states)
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_metric import BaseMetric


class WavefrontSpeed(BaseMetric):
    """Measure wavefront propagation speed from excitable media state history.

    Wavefronts are identified as resting→excited transitions (state 0 → 1).
    Speed is estimated from the spatial displacement of wavefront centroids
    between consecutive timesteps.

    For GH automata with θ=1 and VN neighborhood, theoretical speed is
    exactly 1 cell/step. Speed CV < 0.15 indicates stable propagation.
    """

    def __init__(self) -> None:
        super().__init__(name="wavefront_speed")

    def required_keys(self) -> list[str]:
        return ["grid", "grid_dims"]

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Compute wavefront speed statistics over the trajectory.

        Parameters
        ----------
        state_history : list[dict]
            Full state trajectory.

        Returns
        -------
        dict with:
            wavefront_speed_mean : float — mean speed (cells/step)
            wavefront_speed_std  : float — std of speed
            wavefront_speed_cv   : float — coefficient of variation
            wavefront_count_mean : float — mean number of wavefront cells
            speeds               : np.ndarray — per-step speed estimates
        """
        if len(state_history) < 3:
            return {
                "wavefront_speed_mean": 0.0,
                "wavefront_speed_std": 0.0,
                "wavefront_speed_cv": float("inf"),
                "wavefront_count_mean": 0.0,
                "speeds": np.array([]),
            }

        speeds = []
        wavefront_counts = []
        prev_wf_centroid = None

        for t in range(1, len(state_history)):
            grid_curr = state_history[t]["grid"]
            grid_prev = state_history[t - 1]["grid"]

            # Wavefront: cells that are excited now and were resting before
            wavefront = (grid_curr == 1) & (grid_prev == 0)
            wf_count = int(wavefront.sum())
            wavefront_counts.append(wf_count)

            if wf_count == 0:
                prev_wf_centroid = None
                continue

            # Centroid of wavefront
            wf_positions = np.argwhere(wavefront)
            centroid = wf_positions.mean(axis=0)

            if prev_wf_centroid is not None:
                # Speed = displacement of centroid
                # Handle periodic boundary wrapping
                rows, cols = state_history[t]["grid_dims"]
                dr = centroid[0] - prev_wf_centroid[0]
                dc = centroid[1] - prev_wf_centroid[1]
                # Minimum image convention for periodic BC
                if abs(dr) > rows / 2:
                    dr = dr - np.sign(dr) * rows
                if abs(dc) > cols / 2:
                    dc = dc - np.sign(dc) * cols
                speed = np.sqrt(dr**2 + dc**2)
                speeds.append(speed)

            prev_wf_centroid = centroid

        speeds = np.array(speeds)
        wf_counts = np.array(wavefront_counts)

        if len(speeds) == 0:
            return {
                "wavefront_speed_mean": 0.0,
                "wavefront_speed_std": 0.0,
                "wavefront_speed_cv": float("inf"),
                "wavefront_count_mean": float(wf_counts.mean()) if len(wf_counts) > 0 else 0.0,
                "speeds": speeds,
            }

        mean_speed = float(speeds.mean())
        std_speed = float(speeds.std())
        cv_speed = std_speed / mean_speed if mean_speed > 0 else float("inf")

        return {
            "wavefront_speed_mean": mean_speed,
            "wavefront_speed_std": std_speed,
            "wavefront_speed_cv": cv_speed,
            "wavefront_count_mean": float(wf_counts.mean()),
            "speeds": speeds,
        }


class WavefrontSpeedLocal(BaseMetric):
    """Per-cell wavefront speed via inter-excitation interval.

    For each cell, measure the time between successive excitations
    (inter-excitation interval = temporal period T). The spatial period
    (wavelength λ) is κ — the number of states in the color wheel,
    since the repeating spatial unit is one full cycle
    (..., 0, 1, 2, ..., κ-1, 0, 1, 2, ...). Phase speed = λ / T.

    For a GH automaton with Moore neighborhood and θ=1, T = κ in
    spiral waves, giving speed = κ/κ = 1.0 cells/step — matching
    the true wavefront propagation speed.

    Note: for Von Neumann spirals, T may be slightly larger than κ
    due to the spiral geometry (wavefront zigzags), giving apparent
    speed < 1.0. This is a real geometric effect, not a bug.

    More robust than centroid-based speed for multi-spiral systems
    where centroids move erratically.
    """

    def __init__(self) -> None:
        super().__init__(name="wavefront_speed_local")

    def required_keys(self) -> list[str]:
        return ["grid", "grid_dims"]

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Compute local wavefront speed from inter-excitation intervals."""
        if len(state_history) < 10:
            return {
                "wavefront_speed_mean": 0.0,
                "wavefront_speed_std": 0.0,
                "wavefront_speed_cv": float("inf"),
            }

        rows, cols = state_history[0]["grid_dims"]
        n_states = state_history[0].get("n_states", 3)
        wavelength = n_states  # spatial period: one full cycle (0, 1, ..., κ-1)

        # Track excitation times for a sample of cells
        n_sample = min(500, rows * cols)
        rng = np.random.default_rng(0)
        sample_indices = rng.choice(rows * cols, size=n_sample, replace=False)
        sample_rc = [(idx // cols, idx % cols) for idx in sample_indices]

        # Collect excitation times
        excitation_times: dict[int, list[int]] = {i: [] for i in range(n_sample)}

        for t, state in enumerate(state_history):
            grid = state["grid"]
            for i, (r, c) in enumerate(sample_rc):
                if grid[r, c] == 1:
                    excitation_times[i].append(t)

        # Compute inter-excitation intervals → local periods
        all_speeds = []
        for i in range(n_sample):
            times = excitation_times[i]
            if len(times) < 2:
                continue
            intervals = np.diff(times)
            # Filter out consecutive excitations (artifact of recording)
            intervals = intervals[intervals > 1]
            if len(intervals) == 0:
                continue
            # speed = wavelength / period
            local_speeds = wavelength / intervals
            all_speeds.extend(local_speeds.tolist())

        speeds = np.array(all_speeds)
        if len(speeds) == 0:
            return {
                "wavefront_speed_mean": 0.0,
                "wavefront_speed_std": 0.0,
                "wavefront_speed_cv": float("inf"),
            }

        mean_s = float(speeds.mean())
        std_s = float(speeds.std())
        cv_s = std_s / mean_s if mean_s > 0 else float("inf")

        return {
            "wavefront_speed_mean": mean_s,
            "wavefront_speed_std": std_s,
            "wavefront_speed_cv": cv_s,
        }


class SpiralTipDetector(BaseMetric):
    """Detect spiral wave tips via topological charge.

    In an excitable medium with κ states arranged on a color wheel,
    a spiral tip is a point where all κ states meet. The topological
    charge Q is computed by summing phase differences around each 2×2
    plaquette. Q = ±1 indicates a spiral tip; Q = 0 indicates no tip.

    For periodic media, the total topological charge sums to 0.
    """

    def __init__(self) -> None:
        super().__init__(name="spiral_tips")

    def required_keys(self) -> list[str]:
        return ["grid", "grid_dims", "n_states"]

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Count spiral tips at each timestep.

        Returns
        -------
        dict with:
            tip_count_mean : float — mean number of tips over trajectory
            tip_count_std  : float
            tip_positions  : list[list[tuple]] — tip (r,c) at each recorded step
            net_charge     : list[int] — net topological charge at each step
            max_tip_persistence : int — longest streak a tip region is active
        """
        t_start = kwargs.get("t_start", 0)
        t_end = kwargs.get("t_end", len(state_history))

        tip_counts = []
        net_charges = []
        all_tip_positions = []

        for t in range(t_start, t_end):
            state = state_history[t]
            grid = state["grid"]
            n_states = state["n_states"]
            rows, cols = state["grid_dims"]

            tips, charges = self._find_tips(grid, rows, cols, n_states)
            tip_counts.append(len(tips))
            net_charges.append(sum(charges))
            all_tip_positions.append(tips)

        tc = np.array(tip_counts)

        return {
            "tip_count_mean": float(tc.mean()) if len(tc) > 0 else 0.0,
            "tip_count_std": float(tc.std()) if len(tc) > 0 else 0.0,
            "tip_positions": all_tip_positions,
            "net_charge": net_charges,
            "tip_count_timeseries": tip_counts,
        }

    @staticmethod
    def _find_tips(
        grid: np.ndarray,
        rows: int,
        cols: int,
        n_states: int,
    ) -> tuple[list[tuple[int, int]], list[int]]:
        """Find spiral tips via topological charge on 2×2 plaquettes.

        The phase at each cell is φ = 2π × state / κ. The topological
        charge of a plaquette is Q = (1/2π) Σ Δφ around the plaquette,
        where Δφ is the phase difference wrapped to [-π, π].

        Returns (tip_positions, tip_charges).
        """
        tips = []
        charges = []

        for r in range(rows):
            for c in range(cols):
                # 2×2 plaquette: (r,c), (r,c+1), (r+1,c+1), (r+1,c)
                # Periodic boundary
                r1 = (r + 1) % rows
                c1 = (c + 1) % cols

                states = [
                    grid[r, c],
                    grid[r, c1],
                    grid[r1, c1],
                    grid[r1, c],
                ]

                # Phase = 2π × state / κ
                phases = [2 * np.pi * s / n_states for s in states]

                # Sum phase differences around plaquette
                total_diff = 0.0
                for i in range(4):
                    diff = phases[(i + 1) % 4] - phases[i]
                    # Wrap to [-π, π]
                    diff = (diff + np.pi) % (2 * np.pi) - np.pi
                    total_diff += diff

                # Topological charge
                charge = round(total_diff / (2 * np.pi))
                if charge != 0:
                    tips.append((r, c))
                    charges.append(charge)

        return tips, charges


class WavePersistence(BaseMetric):
    """Measure persistence of wave activity over time.

    Checks whether excitable wave dynamics are sustained (not transient).
    For P13 screening, we need persistent wavefronts for ≥ 5 × T_prop.
    """

    def __init__(self) -> None:
        super().__init__(name="wave_persistence")

    def required_keys(self) -> list[str]:
        return ["grid", "grid_dims"]

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Compute wave persistence metrics.

        Returns
        -------
        dict with:
            active_fraction      : float — fraction of timesteps with activity
            activity_density_mean: float — mean excited fraction
            activity_density_std : float
            longest_active_streak: int — longest run of consecutive active steps
            died_out             : bool — whether activity went to zero permanently
        """
        activities = []
        for state in state_history:
            grid = state["grid"]
            excited = (grid == 1).sum()
            total = grid.size
            activities.append(excited / total)

        activities = np.array(activities)
        active = activities > 0

        # Longest active streak
        streaks = []
        current_streak = 0
        for a in active:
            if a:
                current_streak += 1
            else:
                if current_streak > 0:
                    streaks.append(current_streak)
                current_streak = 0
        if current_streak > 0:
            streaks.append(current_streak)

        longest = max(streaks) if streaks else 0

        # Check if died out (last 10% all dead)
        late_start = max(0, int(len(activities) * 0.9))
        died_out = all(activities[late_start:] == 0) if late_start < len(activities) else True

        return {
            "active_fraction": float(active.mean()),
            "activity_density_mean": float(activities.mean()),
            "activity_density_std": float(activities.std()),
            "longest_active_streak": longest,
            "died_out": died_out,
        }
