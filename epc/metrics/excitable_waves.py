"""Excitable wave metrics for P13 detection.

Metrics for detecting and characterizing excitable spiral and target waves
in cellular automata with rest/excited/refractory states.

- WavefrontSpeed: propagation speed and coefficient of variation
- SpiralTipDetector: topological charge method for spiral tip counting
- WaveSourceCount: persistent excitation sources
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_metric import BaseMetric


class WavefrontSpeed(BaseMetric):
    """Measure wavefront propagation speed via center-of-mass tracking.

    Wavefronts are cells that transition rest(0) → excited(1) between
    consecutive frames. Speed is the frame-to-frame displacement of the
    wavefront center of mass.
    """

    def __init__(self) -> None:
        super().__init__(name="wavefront_speed")

    def required_keys(self) -> list[str]:
        return ["grid", "grid_dims"]

    def compute(
        self, state_history: list[dict[str, Any]], **kwargs: Any
    ) -> dict[str, Any]:
        min_wavefront_cells: int = kwargs.get("min_wavefront_cells", 5)
        speeds: list[float] = []
        prev_com: np.ndarray | None = None
        persistence_frames = 0
        total_wavefront_frames = 0

        for t in range(1, len(state_history)):
            prev_grid = state_history[t - 1]["grid"]
            curr_grid = state_history[t]["grid"]

            # Wavefront: cells that are excited now and were resting before
            wavefront = (curr_grid == 1) & (prev_grid == 0)
            n_front = int(wavefront.sum())

            if n_front < min_wavefront_cells:
                prev_com = None
                continue

            total_wavefront_frames += 1
            # Center of mass of wavefront cells
            rows, cols = np.where(wavefront)
            com = np.array([rows.mean(), cols.mean()])

            if prev_com is not None:
                # Periodic distance (use grid_dims)
                dims = np.array(state_history[t]["grid_dims"], dtype=float)
                delta = com - prev_com
                # Minimum image convention
                delta = delta - dims * np.round(delta / dims)
                speed = float(np.linalg.norm(delta))
                speeds.append(speed)

            prev_com = com

        if not speeds:
            return {
                "wavefront_speed": self.name,
                "wavefront_speed_mean": 0.0,
                "wavefront_speed_std": 0.0,
                "wavefront_speed_cv": float("inf"),
                "n_wavefront_frames": total_wavefront_frames,
                "n_speed_measurements": 0,
            }

        arr = np.array(speeds)
        mean_speed = float(arr.mean())
        std_speed = float(arr.std())
        cv = std_speed / mean_speed if mean_speed > 0 else float("inf")

        return {
            self.name: mean_speed,
            "wavefront_speed_mean": mean_speed,
            "wavefront_speed_std": std_speed,
            "wavefront_speed_cv": cv,
            "n_wavefront_frames": total_wavefront_frames,
            "n_speed_measurements": len(speeds),
        }


class SpiralTipDetector(BaseMetric):
    """Detect spiral wave tips via topological charge on 2x2 plaquettes.

    Maps discrete cell states to phases: phase = 2π·state / n_states.
    A spiral tip exists where the winding number around a 2x2 plaquette
    is ±1 (topological charge ±1).

    Reference: Barkley, D. (1991). A model for fast computer simulation
    of waves in excitable media. Physica D.
    """

    def __init__(self) -> None:
        super().__init__(name="spiral_tip")

    def required_keys(self) -> list[str]:
        return ["grid", "grid_dims"]

    def compute(
        self, state_history: list[dict[str, Any]], **kwargs: Any
    ) -> dict[str, Any]:
        n_states: int = kwargs.get("n_states", None)
        max_track_distance: float = kwargs.get("max_track_distance", 5.0)

        tip_counts: list[int] = []
        tip_positions_per_frame: list[list[tuple[int, int]]] = []

        for state in state_history:
            grid = state["grid"]
            if n_states is None:
                n_states = int(grid.max()) + 1
            tips = self._find_tips(grid, n_states)
            tip_counts.append(len(tips))
            tip_positions_per_frame.append(tips)

        # Track tips across frames and count rotations
        max_rotations = 0.0
        tip_persistence = 0
        if any(tip_counts):
            max_rotations, tip_persistence = self._track_tips(
                tip_positions_per_frame,
                state_history[0]["grid_dims"],
                max_track_distance,
            )

        net_charge = 0
        if tip_positions_per_frame and tip_positions_per_frame[-1]:
            grid = state_history[-1]["grid"]
            if n_states is None:
                n_states = int(grid.max()) + 1
            net_charge = self._net_charge(grid, n_states)

        return {
            self.name: float(max(tip_counts)) if tip_counts else 0.0,
            "n_spiral_tips_max": int(max(tip_counts)) if tip_counts else 0,
            "n_spiral_tips_mean": float(np.mean(tip_counts)) if tip_counts else 0.0,
            "max_rotations": max_rotations,
            "tip_persistence_frames": tip_persistence,
            "topological_charge_net": net_charge,
        }

    @staticmethod
    def _find_tips(
        grid: np.ndarray, n_states: int
    ) -> list[tuple[int, int]]:
        """Find spiral tips via topological charge on 2x2 plaquettes."""
        if n_states < 3:
            return []

        rows, cols = grid.shape
        phase = 2.0 * np.pi * grid.astype(float) / n_states
        tips = []

        for r in range(rows):
            for c in range(cols):
                # 2x2 plaquette with periodic wrap
                r1 = (r + 1) % rows
                c1 = (c + 1) % cols
                corners = [
                    phase[r, c],
                    phase[r, c1],
                    phase[r1, c1],
                    phase[r1, c],
                ]
                # Compute winding number
                winding = 0.0
                for i in range(4):
                    diff = corners[(i + 1) % 4] - corners[i]
                    # Wrap to (-π, π]
                    diff = (diff + np.pi) % (2 * np.pi) - np.pi
                    winding += diff

                charge = winding / (2 * np.pi)
                if abs(charge) > 0.4:  # Should be ±1 for a tip
                    tips.append((r, c))

        return tips

    @staticmethod
    def _net_charge(grid: np.ndarray, n_states: int) -> int:
        """Compute net topological charge (should be 0 for periodic BC)."""
        tips = SpiralTipDetector._find_tips(grid, n_states)
        if not tips:
            return 0
        rows, cols = grid.shape
        phase = 2.0 * np.pi * grid.astype(float) / n_states
        total = 0.0
        for r, c in tips:
            r1 = (r + 1) % rows
            c1 = (c + 1) % cols
            corners = [phase[r, c], phase[r, c1], phase[r1, c1], phase[r1, c]]
            winding = 0.0
            for i in range(4):
                diff = corners[(i + 1) % 4] - corners[i]
                diff = (diff + np.pi) % (2 * np.pi) - np.pi
                winding += diff
            total += winding / (2 * np.pi)
        return int(round(total))

    @staticmethod
    def _track_tips(
        tips_per_frame: list[list[tuple[int, int]]],
        grid_dims: tuple[int, int],
        max_distance: float,
    ) -> tuple[float, int]:
        """Track tips across frames via nearest-neighbor linking.

        Returns (max_rotations, max_persistence_frames).
        Rotation estimated as cumulative angular displacement / 2π.
        """
        if not any(tips_per_frame):
            return 0.0, 0

        # Simple nearest-neighbor tracker
        tracks: list[list[tuple[int, int]]] = []  # each track = list of positions
        active_track_indices: list[int] = []

        rows, cols = grid_dims

        for frame_tips in tips_per_frame:
            used = set()
            new_active: list[int] = []

            for ti in active_track_indices:
                if not frame_tips:
                    break
                last_pos = tracks[ti][-1]
                # Find closest tip
                best_dist = float("inf")
                best_j = -1
                for j, tp in enumerate(frame_tips):
                    if j in used:
                        continue
                    dr = abs(tp[0] - last_pos[0])
                    dc = abs(tp[1] - last_pos[1])
                    dr = min(dr, rows - dr)
                    dc = min(dc, cols - dc)
                    d = (dr**2 + dc**2) ** 0.5
                    if d < best_dist:
                        best_dist = d
                        best_j = j
                if best_j >= 0 and best_dist <= max_distance:
                    tracks[ti].append(frame_tips[best_j])
                    used.add(best_j)
                    new_active.append(ti)

            # Start new tracks for unmatched tips
            for j, tp in enumerate(frame_tips):
                if j not in used:
                    tracks.append([tp])
                    new_active.append(len(tracks) - 1)

            active_track_indices = new_active

        if not tracks:
            return 0.0, 0

        # Estimate rotations per track: cumulative angular displacement / 2π
        max_rot = 0.0
        max_persist = 0
        for track in tracks:
            if len(track) < 3:
                continue
            max_persist = max(max_persist, len(track))
            # Angular displacement from center of track
            pts = np.array(track, dtype=float)
            center = pts.mean(axis=0)
            angles = np.arctan2(pts[:, 0] - center[0], pts[:, 1] - center[1])
            # Cumulative angular displacement
            diffs = np.diff(angles)
            diffs = (diffs + np.pi) % (2 * np.pi) - np.pi
            total_angle = float(np.abs(diffs).sum())
            rotations = total_angle / (2 * np.pi)
            max_rot = max(max_rot, rotations)

        return max_rot, max_persist


class WaveSourceCount(BaseMetric):
    """Count persistent wave emission sources.

    A wave source is a site that transitions rest(0) → excited(1)
    without any excited neighbors in the previous frame (spontaneous
    excitation or spiral tip emission).
    """

    def __init__(self) -> None:
        super().__init__(name="wave_source_count")

    def required_keys(self) -> list[str]:
        return ["grid", "grid_dims"]

    def compute(
        self, state_history: list[dict[str, Any]], **kwargs: Any
    ) -> dict[str, Any]:
        min_emissions: int = kwargs.get("min_emissions", 3)
        from scipy.signal import convolve2d

        moore = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], dtype=np.int32)
        source_counts: dict[tuple[int, int], int] = {}

        for t in range(1, len(state_history)):
            prev_grid = state_history[t - 1]["grid"]
            curr_grid = state_history[t]["grid"]

            # Newly excited cells
            new_excited = (curr_grid == 1) & (prev_grid == 0)
            if not new_excited.any():
                continue

            # Count excited neighbors in previous frame
            prev_excited = (prev_grid == 1).astype(np.int32)
            n_excited_neighbors = convolve2d(
                prev_excited, moore, mode="same", boundary="wrap"
            )

            # Spontaneous: new_excited AND no excited neighbors previously
            spontaneous = new_excited & (n_excited_neighbors == 0)
            rows, cols = np.where(spontaneous)
            for r, c in zip(rows, cols):
                key = (int(r), int(c))
                source_counts[key] = source_counts.get(key, 0) + 1

        # Persistent sources: emitted >= min_emissions times
        persistent = {k: v for k, v in source_counts.items() if v >= min_emissions}

        return {
            self.name: len(source_counts),
            "wave_source_count": len(source_counts),
            "persistent_sources": len(persistent),
            "source_positions": list(persistent.keys()),
            "max_emissions": max(source_counts.values()) if source_counts else 0,
        }
