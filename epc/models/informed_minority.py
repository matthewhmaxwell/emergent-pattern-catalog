"""
Informed Minority Flocking Model — Couzin et al. (2005).

Reference:
    Couzin, I. D., Krause, J., Franks, N. R., & Levin, S. A. (2005).
    Effective leadership and decision-making in animal groups on the move.
    Nature, 433(7025), 513–516. DOI: 10.1038/nature03236

A Vicsek-style flock where an informed fraction ρ has an additional
weak bias toward a fixed preferred direction. Even a small informed
minority (ρ ~ 0.05–0.10) can steer the entire group. The key signature
is **influence asymmetry**: transfer entropy from informed → naive
exceeds the reverse direction, and group heading aligns with the
preferred direction disproportionately to ρ.

Update rule (synchronous):
    For naive agents:
        θ_i(t+1) = arctan2(⟨sin θ⟩_neighbors, ⟨cos θ⟩_neighbors) + η·ξ_i
    For informed agents:
        θ_i(t+1) = (1-ω) · ⟨θ⟩_neighbors + ω · θ_preferred + η·ξ_i
        (implemented via weighted circular mean)

    x_i(t+1) = x_i(t) + v₀·cos(θ_i(t+1))·Δt
    y_i(t+1) = y_i(t) + v₀·sin(θ_i(t+1))·Δt

Periodic boundary conditions on [0, L) × [0, L).

Parameters from Couzin (2005):
    N = 200, ρ ∈ [0, 0.5], ω = 0.3 (bias weight),
    preferred direction = 0 (East)
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Any, Optional
from scipy.spatial import cKDTree


class InformedMinorityModel:
    """Informed-minority flocking model (Couzin et al. 2005).

    Vicsek-style alignment + an informed fraction with directional bias.
    The informed agents have a weak preference for a fixed direction,
    while naive agents align purely with neighbors.

    State dict keys:
        positions: (N, 2) float64 — particle positions in [0, L)²
        velocities: (N, 2) float64 — particle velocities
        headings: (N,) float64 — heading angles θ in [−π, π)
        informed_mask: (N,) bool — True for informed agents
        step: int — current timestep
    """

    def __init__(
        self,
        n_particles: int = 200,
        box_size: float = 10.0,
        speed: float = 0.03,
        noise: float = 0.1,
        interaction_radius: float = 1.0,
        dt: float = 1.0,
        informed_fraction: float = 0.1,
        bias_weight: float = 0.3,
        preferred_direction: float = 0.0,
        init_mode: str = "random",
        seed: Optional[int] = None,
    ):
        """Initialize informed-minority model.

        Args:
            n_particles: Number of self-propelled particles (N).
            box_size: Side length of periodic square domain (L).
            speed: Constant particle speed (v₀).
            noise: Noise strength η. Angular perturbation Uniform[-η/2, η/2].
            interaction_radius: Metric interaction radius (r).
            dt: Timestep size (Δt).
            informed_fraction: Fraction ρ of agents that are informed (0 to 1).
            bias_weight: Weight ω for preferred direction in informed agents.
            preferred_direction: Preferred heading angle (radians).
            init_mode: One of 'random', 'aligned'.
            seed: Random seed for reproducibility.
        """
        self.n_particles = n_particles
        self.box_size = box_size
        self.speed = speed
        self.noise = noise
        self.interaction_radius = interaction_radius
        self.dt = dt
        self.informed_fraction = informed_fraction
        self.bias_weight = bias_weight
        self.preferred_direction = preferred_direction
        self.init_mode = init_mode
        self.seed = seed

        # Derived
        self.density = n_particles / box_size**2
        self.n_informed = max(0, int(round(n_particles * informed_fraction)))

        # State (set in setup())
        self.positions: NDArray[np.float64] = np.empty(0)
        self.headings: NDArray[np.float64] = np.empty(0)
        self.informed_mask: NDArray[np.bool_] = np.empty(0, dtype=bool)
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self._step_count: int = 0

    def setup(self) -> dict[str, Any]:
        """Initialize particle positions, headings, and informed labels.

        Returns:
            Initial state dict.
        """
        self.rng = np.random.default_rng(self.seed)
        self._step_count = 0
        N = self.n_particles
        L = self.box_size

        if self.init_mode == "random":
            self.positions = self.rng.uniform(0, L, size=(N, 2))
            self.headings = self.rng.uniform(-np.pi, np.pi, size=N)
        elif self.init_mode == "aligned":
            self.positions = self.rng.uniform(0, L, size=(N, 2))
            self.headings = np.zeros(N)
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode!r}")

        # Assign informed agents: first n_informed agents are informed
        self.informed_mask = np.zeros(N, dtype=bool)
        self.informed_mask[:self.n_informed] = True

        return self.get_state()

    def step(self) -> dict[str, Any]:
        """Execute one synchronous update step.

        Naive agents: pure Vicsek alignment with neighbors.
        Informed agents: weighted circular mean of neighbor alignment
        and preferred direction.

        Returns:
            State dict after this step.
        """
        N = self.n_particles
        L = self.box_size
        r = self.interaction_radius

        # --- Neighbor finding with periodic boundary conditions ---
        tree = cKDTree(self.positions, boxsize=L)
        W = tree.sparse_distance_matrix(tree, r, output_type="coo_matrix")

        from scipy.sparse import csr_matrix
        adjacency = csr_matrix(
            (np.ones(W.nnz, dtype=np.float64), (W.row, W.col)),
            shape=(N, N),
        )

        # --- Compute mean heading of neighbors (circular mean) ---
        sin_h = np.sin(self.headings)
        cos_h = np.cos(self.headings)

        sum_sin = adjacency @ sin_h
        sum_cos = adjacency @ cos_h

        n_neighbors = np.array(adjacency.sum(axis=1)).ravel()
        n_neighbors = np.maximum(n_neighbors, 1)
        mean_sin = sum_sin / n_neighbors
        mean_cos = sum_cos / n_neighbors

        # Social heading from neighbor alignment
        social_heading = np.arctan2(mean_sin, mean_cos)

        # --- Combine social + preferred direction for informed agents ---
        new_headings = social_heading.copy()

        if self.n_informed > 0:
            omega = self.bias_weight
            informed = self.informed_mask

            # Weighted circular mean: (1-ω) * social + ω * preferred
            social_sin = np.sin(social_heading[informed])
            social_cos = np.cos(social_heading[informed])
            pref_sin = np.sin(self.preferred_direction)
            pref_cos = np.cos(self.preferred_direction)

            combined_sin = (1 - omega) * social_sin + omega * pref_sin
            combined_cos = (1 - omega) * social_cos + omega * pref_cos

            new_headings[informed] = np.arctan2(combined_sin, combined_cos)

        # --- Add angular noise ---
        noise_vals = self.rng.uniform(-self.noise / 2, self.noise / 2, size=N)
        new_headings += noise_vals

        # --- Wrap headings to [-π, π) ---
        self.headings = (new_headings + np.pi) % (2 * np.pi) - np.pi

        # --- Update positions ---
        vx = self.speed * np.cos(self.headings) * self.dt
        vy = self.speed * np.sin(self.headings) * self.dt
        self.positions[:, 0] += vx
        self.positions[:, 1] += vy

        # --- Periodic boundary conditions ---
        self.positions %= L

        self._step_count += 1
        return self.get_state()

    def run(self, n_steps: int, record_interval: int = 1) -> list[dict[str, Any]]:
        """Run simulation for n_steps, recording state at intervals.

        Args:
            n_steps: Number of timesteps to run.
            record_interval: Record state every this many steps.

        Returns:
            List of state dicts (including initial state at index 0).
        """
        history = [self.setup()]
        for t in range(1, n_steps + 1):
            self.step()
            if t % record_interval == 0:
                history.append(self.get_state())
        return history

    def get_state(self) -> dict[str, Any]:
        """Return current state as a dict with numpy arrays.

        Returns:
            Dict with keys: positions, velocities, headings,
            informed_mask, step.
        """
        velocities = np.column_stack([
            self.speed * np.cos(self.headings),
            self.speed * np.sin(self.headings),
        ])
        return {
            "positions": self.positions.copy(),
            "velocities": velocities,
            "headings": self.headings.copy(),
            "informed_mask": self.informed_mask.copy(),
            "step": self._step_count,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata for detector use."""
        return {
            "model_family": "informed_minority",
            "model_name": "Informed Minority Flocking (Couzin 2005)",
            "model_class": "informed_minority",
            "n_particles": self.n_particles,
            "box_size": self.box_size,
            "speed": self.speed,
            "noise": self.noise,
            "interaction_radius": self.interaction_radius,
            "dt": self.dt,
            "density": self.density,
            "informed_fraction": self.informed_fraction,
            "n_informed": self.n_informed,
            "bias_weight": self.bias_weight,
            "preferred_direction": self.preferred_direction,
            "init_mode": self.init_mode,
            "seed": self.seed,
            "space_type": "continuous_2d_periodic",
            "dynamics_type": "first_order_heading",
            "has_alignment_rule": True,
            "has_attraction_rule": False,
            "has_informed_minority": True,
        }

    def get_timescale(self) -> dict[str, float]:
        """Return characteristic timescales."""
        return {
            "T_cross": self.box_size / self.speed,
            "T_align": self.interaction_radius / self.speed,
        }


def group_accuracy(state: dict[str, Any], preferred_direction: float) -> float:
    """Compute group directional accuracy relative to preferred direction.

    accuracy = cos(θ_group - θ_preferred) where θ_group is the mean heading.
    Returns value in [-1, 1], where 1 = perfect alignment with preferred.

    Args:
        state: State dict with 'headings' key.
        preferred_direction: The target direction (radians).

    Returns:
        Accuracy ∈ [-1, 1].
    """
    headings = state["headings"]
    mean_vx = np.mean(np.cos(headings))
    mean_vy = np.mean(np.sin(headings))
    group_heading = np.arctan2(mean_vy, mean_vx)
    return float(np.cos(group_heading - preferred_direction))


def polarization(state: dict[str, Any]) -> float:
    """Compute polarization order parameter φ.

    φ = |(1/N) Σ (cos θ_i, sin θ_i)|. φ = 1 perfect alignment, φ ≈ 1/√N disordered.
    """
    headings = state["headings"]
    mean_vx = np.mean(np.cos(headings))
    mean_vy = np.mean(np.sin(headings))
    return float(np.sqrt(mean_vx**2 + mean_vy**2))
