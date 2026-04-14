"""
Standard Vicsek Model (SVM) — Vicsek et al. 1995.

Reference:
    Vicsek, T., Czirók, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995).
    Novel type of phase transition in a system of self-driven particles.
    Physical Review Letters, 75(6), 1226–1229.
    DOI: 10.1103/PhysRevLett.75.1226

Update rule (synchronous, angular noise variant):
    θ_i(t+1) = arctan2(⟨sin θ⟩_neighbors, ⟨cos θ⟩_neighbors) + η·ξ_i
    x_i(t+1) = x_i(t) + v₀·cos(θ_i(t+1))·Δt
    y_i(t+1) = y_i(t) + v₀·sin(θ_i(t+1))·Δt

where ξ_i ~ Uniform[-1/2, 1/2] and neighbors are all particles j with
|r_i - r_j| ≤ r (metric interaction, including self).

Periodic boundary conditions on [0, L) × [0, L).

Key parameters from original paper (Figure 1):
    N = 300, L = 7 (ρ = N/L² ≈ 6.12), v₀ = 0.03, r = 1, Δt = 1

Validation targets:
    1. Phase transition: φ → 1 at low η, φ → 0 at high η
    2. Disordered baseline: φ ≈ 1/√N
    3. Critical scaling: |v_a| ~ (η_c - η)^β, β ≈ 0.45
    4. Qualitative snapshots matching paper Figure 1(a-d)
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Any, Optional
from scipy.spatial import cKDTree


class VicsekModel:
    """Standard Vicsek Model (1995) for self-propelled particles.

    Implements the exact angular-noise variant from the original paper.
    All particles move at constant speed v0, align with local neighbors
    within metric radius r, subject to angular noise of strength eta.

    State dict keys:
        positions: (N, 2) float64 — particle positions in [0, L)²
        velocities: (N, 2) float64 — particle velocities (v0 * [cos θ, sin θ])
        headings: (N,) float64 — heading angles θ in [−π, π)
        step: int — current timestep

    Init modes:
        random: uniform random positions and headings
        aligned: uniform positions, all headings = 0 (test perfect order)
        half_aligned: half θ=0, half θ=π (test symmetry breaking)
        single_cluster: particles in central L/4 box, random headings
    """

    def __init__(
        self,
        n_particles: int = 300,
        box_size: float = 7.0,
        speed: float = 0.03,
        noise: float = 0.1,
        interaction_radius: float = 1.0,
        dt: float = 1.0,
        init_mode: str = "random",
        seed: Optional[int] = None,
    ):
        """Initialize Vicsek model.

        Args:
            n_particles: Number of self-propelled particles (N).
            box_size: Side length of periodic square domain (L).
            speed: Constant particle speed (v₀).
            noise: Noise strength η ∈ [0, 2π]. Angular perturbation
                   is uniform in [-η/2, η/2].
            interaction_radius: Metric interaction radius (r).
            dt: Timestep size (Δt).
            init_mode: One of 'random', 'aligned', 'half_aligned',
                       'single_cluster'.
            seed: Random seed for reproducibility.
        """
        self.n_particles = n_particles
        self.box_size = box_size
        self.speed = speed
        self.noise = noise
        self.interaction_radius = interaction_radius
        self.dt = dt
        self.init_mode = init_mode
        self.seed = seed

        # Derived
        self.density = n_particles / box_size**2

        # State (set in setup())
        self.positions: NDArray[np.float64] = np.empty(0)
        self.headings: NDArray[np.float64] = np.empty(0)
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self._step_count: int = 0

    def setup(self) -> dict[str, Any]:
        """Initialize particle positions and headings.

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

        elif self.init_mode == "half_aligned":
            self.positions = self.rng.uniform(0, L, size=(N, 2))
            self.headings = np.where(
                np.arange(N) < N // 2, 0.0, np.pi
            )

        elif self.init_mode == "single_cluster":
            center = L / 2
            spread = L / 4
            self.positions = self.rng.uniform(
                center - spread / 2, center + spread / 2, size=(N, 2)
            )
            self.positions %= L  # wrap
            self.headings = self.rng.uniform(-np.pi, np.pi, size=N)

        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode!r}")

        return self.get_state()

    def step(self) -> dict[str, Any]:
        """Execute one synchronous Vicsek update step.

        The update is simultaneous: all particles compute their new heading
        from the CURRENT configuration, then all positions are updated.

        Uses vectorized sparse-matrix approach for neighbor averaging.

        Returns:
            State dict after this step.
        """
        N = self.n_particles
        L = self.box_size
        r = self.interaction_radius

        # --- Neighbor finding with periodic boundary conditions ---
        # Build sparse neighbor matrix: W[i,j] = 1 if |r_i - r_j| ≤ r
        tree = cKDTree(self.positions, boxsize=L)
        W = tree.sparse_distance_matrix(tree, r, output_type="coo_matrix")

        # Convert to CSR for efficient row operations
        # W.data are distances; we just need the adjacency structure
        from scipy.sparse import csr_matrix
        adjacency = csr_matrix(
            (np.ones(W.nnz, dtype=np.float64), (W.row, W.col)),
            shape=(N, N),
        )

        # --- Compute mean heading of neighbors (circular mean) ---
        # Sum sin(θ) and cos(θ) over neighbors using sparse matrix-vector product
        sin_h = np.sin(self.headings)
        cos_h = np.cos(self.headings)

        sum_sin = adjacency @ sin_h   # (N,) sum of sin(θ_j) for neighbors of i
        sum_cos = adjacency @ cos_h   # (N,) sum of cos(θ_j) for neighbors of i

        # Divide by neighbor count to get mean
        n_neighbors = np.array(adjacency.sum(axis=1)).ravel()
        n_neighbors = np.maximum(n_neighbors, 1)  # safety (shouldn't be 0)
        mean_sin = sum_sin / n_neighbors
        mean_cos = sum_cos / n_neighbors

        new_headings = np.arctan2(mean_sin, mean_cos)

        # --- Add angular noise: Δθ ~ Uniform[-η/2, η/2] ---
        noise = self.rng.uniform(
            -self.noise / 2, self.noise / 2, size=N
        )
        new_headings += noise

        # --- Update headings (wrap to [-π, π) for cleanliness) ---
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
            Dict with keys: positions, velocities, headings, step,
            plus metadata keys (n_particles, box_size, etc.).
        """
        velocities = np.column_stack([
            self.speed * np.cos(self.headings),
            self.speed * np.sin(self.headings),
        ])
        return {
            "positions": self.positions.copy(),
            "velocities": velocities,
            "headings": self.headings.copy(),
            "step": self._step_count,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata for detector use."""
        return {
            "model_family": "vicsek",
            "model_name": "Standard Vicsek Model (1995)",
            "n_particles": self.n_particles,
            "box_size": self.box_size,
            "speed": self.speed,
            "noise": self.noise,
            "interaction_radius": self.interaction_radius,
            "dt": self.dt,
            "density": self.density,
            "init_mode": self.init_mode,
            "seed": self.seed,
            "space_type": "continuous_2d_periodic",
            "dynamics_type": "first_order_heading",
        }

    def get_timescale(self) -> dict[str, float]:
        """Return characteristic timescales for this model.

        T_cross = L / v₀ — time for a particle to cross the domain.
        T_align = r / v₀ — time for a particle to traverse one
                  interaction radius (local mixing time).
        """
        return {
            "T_cross": self.box_size / self.speed,
            "T_align": self.interaction_radius / self.speed,
        }


def polarization(state: dict[str, Any]) -> float:
    """Compute instantaneous polarization order parameter φ.

    φ(t) = |(1/N) Σ v̂_i(t)| = |(1/N) Σ (cos θ_i, sin θ_i)|

    φ = 1 means perfect alignment (all same direction).
    φ ≈ 1/√N means disordered (random headings).

    Args:
        state: State dict with 'headings' key.

    Returns:
        Polarization φ ∈ [0, 1].
    """
    headings = state["headings"]
    N = len(headings)
    mean_vx = np.mean(np.cos(headings))
    mean_vy = np.mean(np.sin(headings))
    return np.sqrt(mean_vx**2 + mean_vy**2)


def angular_momentum(state: dict[str, Any], box_size: float) -> float:
    """Compute angular momentum order parameter L.

    L(t) = (1/N) Σ (r̂_i × v̂_i)

    where r̂_i is the unit vector from center of mass to particle i,
    and v̂_i is the unit heading vector.

    |L| near 1 means milling (rotation), |L| near 0 means no rotation.

    Args:
        state: State dict with 'positions' and 'headings' keys.
        box_size: Domain size for periodic unwrapping.

    Returns:
        Angular momentum L ∈ [-1, 1].
    """
    positions = state["positions"]
    headings = state["headings"]
    N = len(headings)

    # Center of mass (periodic-aware)
    # Use circular mean of positions for periodic BC
    theta_x = 2 * np.pi * positions[:, 0] / box_size
    theta_y = 2 * np.pi * positions[:, 1] / box_size
    com_x = box_size / (2 * np.pi) * np.arctan2(
        np.mean(np.sin(theta_x)), np.mean(np.cos(theta_x))
    ) % box_size
    com_y = box_size / (2 * np.pi) * np.arctan2(
        np.mean(np.sin(theta_y)), np.mean(np.cos(theta_y))
    ) % box_size

    # Displacement from COM (minimum image)
    dx = positions[:, 0] - com_x
    dy = positions[:, 1] - com_y
    dx -= box_size * np.round(dx / box_size)
    dy -= box_size * np.round(dy / box_size)

    # Unit radial vectors
    dist = np.sqrt(dx**2 + dy**2)
    dist = np.maximum(dist, 1e-12)  # avoid division by zero
    rx_hat = dx / dist
    ry_hat = dy / dist

    # Unit heading vectors
    vx_hat = np.cos(headings)
    vy_hat = np.sin(headings)

    # Cross product (2D): r̂ × v̂ = rx*vy - ry*vx
    cross = rx_hat * vy_hat - ry_hat * vx_hat

    return np.mean(cross)


def group_speed_ratio(state: dict[str, Any]) -> float:
    """Compute group speed ratio R = |V_cm| / ⟨|v_i|⟩.

    For constant-speed Vicsek, ⟨|v_i|⟩ = v₀, so R = |V_cm| / v₀.
    R ≈ φ for constant-speed models. R > 0.5 indicates net transport.

    Args:
        state: State dict with 'velocities' key.

    Returns:
        Group speed ratio R ∈ [0, 1].
    """
    velocities = state["velocities"]
    v_cm = np.mean(velocities, axis=0)
    v_cm_mag = np.linalg.norm(v_cm)
    mean_speed = np.mean(np.linalg.norm(velocities, axis=1))
    if mean_speed < 1e-12:
        return 0.0
    return v_cm_mag / mean_speed
