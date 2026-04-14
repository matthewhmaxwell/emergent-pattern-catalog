"""Vicsek flocking model.

N self-propelled particles in a 2D periodic box. Each agent moves at
constant speed and aligns its heading with neighbors within an interaction
radius, plus noise. Exhibits an order-disorder phase transition.

Optional attraction_strength parameter biases heading toward the group
center of mass, enabling milling / vortex formation (P6).

Canonical model for P5 (translational alignment / flocking) and,
with attraction, P6 (milling / vortex formation).

Reference: Vicsek, T. et al. (1995). Novel type of phase transition in
a system of self-driven particles. Physical Review Letters, 75(6), 1226.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from epc.base_model import BaseModel


def _circular_mean(angles: np.ndarray) -> float:
    """Compute circular mean of angles in radians."""
    return float(np.arctan2(np.mean(np.sin(angles)), np.mean(np.cos(angles))))


def _periodic_com(positions: np.ndarray, box_size: float) -> np.ndarray:
    """Compute center of mass in periodic domain using circular mean trick."""
    com = np.empty(2)
    for d in range(2):
        theta = 2.0 * np.pi * positions[:, d] / box_size
        com[d] = box_size * np.arctan2(
            np.mean(np.sin(theta)), np.mean(np.cos(theta))
        ) / (2.0 * np.pi)
        # Wrap to [0, box_size)
        com[d] = com[d] % box_size
    return com


class VicsekModel(BaseModel):
    """Vicsek self-propelled particle model with optional attraction.

    Parameters
    ----------
    n_agents : int
        Number of particles.
    box_size : float
        Side length of the periodic square domain.
    speed : float
        Constant speed v₀ for all agents.
    noise : float
        Noise amplitude η. Heading perturbation drawn from
        Uniform(-η·π, η·π). Set η=0 for deterministic alignment.
    interaction_radius : float
        Agents within this distance interact (metric interaction).
    attraction_strength : float
        Bias toward group center of mass (0 = standard Vicsek).
        Heading update becomes weighted average of alignment heading
        and direction-to-COM. Enables milling at moderate values (~0.3).
    """

    def __init__(
        self,
        n_agents: int = 300,
        box_size: float = 25.0,
        speed: float = 0.03,
        noise: float = 0.1,
        interaction_radius: float = 1.0,
        attraction_strength: float = 0.0,
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)
        self.n_agents = n_agents
        self.box_size = box_size
        self.speed = speed
        self.noise = noise
        self.interaction_radius = interaction_radius
        self.attraction_strength = attraction_strength
        self._positions: np.ndarray = np.empty((n_agents, 2))
        self._headings: np.ndarray = np.empty(n_agents)

    def setup(self) -> dict[str, Any]:
        self._positions = self.rng.uniform(
            0, self.box_size, size=(self.n_agents, 2)
        )
        self._headings = self.rng.uniform(0, 2 * np.pi, size=self.n_agents)
        self._step_count = 0
        self._is_setup = True
        return self._snapshot()

    def step(self) -> dict[str, Any]:
        L = self.box_size
        R = self.interaction_radius

        # Build periodic KD-tree
        tree = cKDTree(self._positions, boxsize=L)
        neighbor_lists = tree.query_ball_tree(tree, r=R)

        new_headings = np.empty(self.n_agents)

        for i, neighbors in enumerate(neighbor_lists):
            neighbor_headings = self._headings[neighbors]
            align_heading = np.arctan2(
                np.mean(np.sin(neighbor_headings)),
                np.mean(np.cos(neighbor_headings)),
            )

            if self.attraction_strength > 0:
                # Bias toward COM
                com = _periodic_com(self._positions, L)
                # Periodic displacement from agent to COM
                dx = com[0] - self._positions[i, 0]
                dy = com[1] - self._positions[i, 1]
                # Minimum-image convention
                dx = dx - L * np.round(dx / L)
                dy = dy - L * np.round(dy / L)
                attract_heading = np.arctan2(dy, dx)

                # Weighted circular average
                w_align = 1.0 - self.attraction_strength
                w_attract = self.attraction_strength
                sx = w_align * np.sin(align_heading) + w_attract * np.sin(attract_heading)
                cx = w_align * np.cos(align_heading) + w_attract * np.cos(attract_heading)
                target_heading = np.arctan2(sx, cx)
            else:
                target_heading = align_heading

            # Add noise
            noise_term = self.rng.uniform(-self.noise * np.pi, self.noise * np.pi)
            new_headings[i] = target_heading + noise_term

        self._headings = new_headings

        # Update positions
        self._positions[:, 0] += self.speed * np.cos(self._headings)
        self._positions[:, 1] += self.speed * np.sin(self._headings)
        # Periodic wrap
        self._positions %= L

        self._step_count += 1
        return self._snapshot()

    def _snapshot(self) -> dict[str, Any]:
        velocities = np.column_stack([
            self.speed * np.cos(self._headings),
            self.speed * np.sin(self._headings),
        ])
        return {
            "positions": self._positions.copy(),
            "velocities": velocities,
            "headings": self._headings.copy(),
            "step": self._step_count,
            "n_agents": self.n_agents,
            "box_size": self.box_size,
        }

    def get_metadata(self) -> dict[str, Any]:
        return {
            "model_name": "vicsek_flocking",
            "model_class": "active_matter",
            "n_agents": self.n_agents,
            "box_size": self.box_size,
            "speed": self.speed,
            "noise": self.noise,
            "interaction_radius": self.interaction_radius,
            "attraction_strength": self.attraction_strength,
            "interaction_type": "metric_radius",
            "update_mode": "synchronous",
            "params": {
                "n_agents": self.n_agents,
                "box_size": self.box_size,
                "speed": self.speed,
                "noise": self.noise,
                "interaction_radius": self.interaction_radius,
                "attraction_strength": self.attraction_strength,
            },
            "seed": self.seed,
            "reference": "Vicsek et al. (1995)",
        }

    def get_timescale(self) -> float:
        """T_cross = box_size / speed (domain crossing time)."""
        if self.speed <= 0:
            return float(self.box_size)
        return self.box_size / self.speed
