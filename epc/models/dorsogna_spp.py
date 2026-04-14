"""
D'Orsogna Self-Propelled Particle Model with Morse Interaction.

Reference:
    D'Orsogna, M.R., Chuang, Y.-L., Bertozzi, A.L., & Chayes, L.S. (2006).
    Self-propelled particles with soft-core interactions: patterns, stability,
    and collapse. Physical Review Letters, 96(10), 104302.
    DOI: 10.1103/PhysRevLett.96.104302

    Carrillo, J.A., D'Orsogna, M.R., & Panferov, V. (2009).
    Double milling in self-propelled swarms from kinetic theory.
    Kinetic and Related Models, 2(2), 363–378.

Equations of motion (2D):
    dx_i/dt = v_i
    dv_i/dt = (α - β|v_i|²) v_i - Σ_{j≠i} ∇U(|x_i - x_j|)

    U(r) = -C_a exp(-r/l_a) + C_r exp(-r/l_r)   (Morse potential)
    ∇U(r) = U'(r) r̂ = [(C_a/l_a)exp(-r/l_a) - (C_r/l_r)exp(-r/l_r)] r̂

Self-propulsion: Rayleigh friction with equilibrium speed v₀ = √(α/β).
    - α drives acceleration, β|v|² provides speed-dependent damping.

Published milling parameters (Carrillo et al. 2009, Fig 3.1):
    N=100, C_a=0.5, C_r=1.0, l_a=3.0, l_r=0.5, α=1.0, β=0.5
    Random IC → single rotating mill.

Validation targets:
    1. Single mill: |L| → 1, φ → 0, R → 0
    2. Equilibrium speed: |v_i| → √(α/β) = √2 ≈ 1.414
    3. Ring-like radial density profile
    4. Stable rotation for extended time
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Any, Optional


class DOrsognaSPPModel:
    """D'Orsogna Self-Propelled Particle Model with Morse interaction.

    Second-order dynamics: positions and velocities are both dynamic.
    Particles self-propel via Rayleigh friction and interact via
    pairwise Morse potential (short-range repulsion, long-range attraction).

    Produces milling (P6), flocking, or collapse depending on parameters.

    State dict keys:
        positions: (N, 2) float64 — particle positions (open space)
        velocities: (N, 2) float64 — particle velocities
        headings: (N,) float64 — heading angles derived from velocities
        step: int — current timestep

    Init modes:
        random: uniform positions in [-R, R]², random heading, speed v₀
        ring: particles on a circle, tangential velocities (starts milling)
        flock: clustered positions, aligned velocities
    """

    def __init__(
        self,
        n_particles: int = 100,
        C_a: float = 0.5,
        C_r: float = 1.0,
        l_a: float = 3.0,
        l_r: float = 0.5,
        alpha: float = 1.0,
        beta: float = 0.5,
        dt: float = 0.01,
        init_mode: str = "random",
        init_radius: float = 5.0,
        seed: Optional[int] = None,
    ):
        """Initialize D'Orsogna SPP model.

        Args:
            n_particles: Number of self-propelled particles (N).
            C_a: Attractive potential strength.
            C_r: Repulsive potential strength.
            l_a: Attractive potential range.
            l_r: Repulsive potential range.
            alpha: Self-propulsion coefficient.
            beta: Speed damping coefficient. Equilibrium speed = √(α/β).
            dt: Integration timestep (RK4).
            init_mode: One of 'random', 'ring', 'flock'.
            init_radius: Characteristic initial spread radius.
            seed: Random seed for reproducibility.
        """
        self.n_particles = n_particles
        self.C_a = C_a
        self.C_r = C_r
        self.l_a = l_a
        self.l_r = l_r
        self.alpha = alpha
        self.beta = beta
        self.dt = dt
        self.init_mode = init_mode
        self.init_radius = init_radius
        self.seed = seed

        # Derived
        self.v_eq = np.sqrt(alpha / beta)  # equilibrium speed

        # State (set in setup())
        self.positions: NDArray[np.float64] = np.empty(0)
        self.velocities: NDArray[np.float64] = np.empty(0)
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self._step_count: int = 0

    def setup(self) -> dict[str, Any]:
        """Initialize particle positions and velocities."""
        self.rng = np.random.default_rng(self.seed)
        self._step_count = 0
        N = self.n_particles
        R = self.init_radius

        if self.init_mode == "random":
            # Random positions in square, random directions at equilibrium speed
            self.positions = self.rng.uniform(-R, R, size=(N, 2))
            angles = self.rng.uniform(-np.pi, np.pi, size=N)
            self.velocities = self.v_eq * np.column_stack([
                np.cos(angles), np.sin(angles)
            ])

        elif self.init_mode == "ring":
            # Particles on a circle with tangential velocities
            angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
            self.positions = R * np.column_stack([
                np.cos(angles), np.sin(angles)
            ])
            # Tangential velocities (counterclockwise)
            self.velocities = self.v_eq * np.column_stack([
                -np.sin(angles), np.cos(angles)
            ])

        elif self.init_mode == "flock":
            # Clustered positions, aligned velocities
            self.positions = self.rng.normal(0, R / 3, size=(N, 2))
            self.velocities = self.v_eq * np.column_stack([
                np.ones(N), np.zeros(N)
            ])

        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode!r}")

        return self.get_state()

    def _compute_forces(
        self,
        positions: NDArray[np.float64],
        velocities: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Compute total acceleration for all particles.

        Returns:
            (N, 2) acceleration array.
        """
        N = len(positions)

        # --- Self-propulsion (Rayleigh friction) ---
        speeds_sq = np.sum(velocities**2, axis=1, keepdims=True)  # (N, 1)
        a_self = (self.alpha - self.beta * speeds_sq) * velocities  # (N, 2)

        # --- Morse interaction forces ---
        # Pairwise displacement vectors: dx[i,j] = x_i - x_j
        dx = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]  # (N, N, 2)

        # Pairwise distances
        dist = np.sqrt(np.sum(dx**2, axis=2))  # (N, N)
        dist_safe = np.maximum(dist, 1e-10)  # avoid division by zero

        # Morse force magnitude: U'(r) = (C_a/l_a)exp(-r/l_a) - (C_r/l_r)exp(-r/l_r)
        # Force on i from j: F_ij = -U'(r_ij) * r̂_ij = -U'(r_ij) * (x_i-x_j)/r_ij
        U_prime = (
            (self.C_a / self.l_a) * np.exp(-dist_safe / self.l_a)
            - (self.C_r / self.l_r) * np.exp(-dist_safe / self.l_r)
        )  # (N, N)

        # Zero self-interaction
        np.fill_diagonal(U_prime, 0.0)

        # Unit displacement vectors
        r_hat = dx / dist_safe[:, :, np.newaxis]  # (N, N, 2)

        # Force on particle i: -Σ_j U'(r_ij) r̂_ij
        a_interaction = -np.sum(
            U_prime[:, :, np.newaxis] * r_hat,
            axis=1,
        )  # (N, 2)

        return a_self + a_interaction

    def step(self) -> dict[str, Any]:
        """Execute one RK4 integration step.

        Returns:
            State dict after this step.
        """
        dt = self.dt
        x = self.positions
        v = self.velocities

        # RK4 for dx/dt = v, dv/dt = f(x, v)
        k1_x = v
        k1_v = self._compute_forces(x, v)

        k2_x = v + 0.5 * dt * k1_v
        k2_v = self._compute_forces(x + 0.5 * dt * k1_x, v + 0.5 * dt * k1_v)

        k3_x = v + 0.5 * dt * k2_v
        k3_v = self._compute_forces(x + 0.5 * dt * k2_x, v + 0.5 * dt * k2_v)

        k4_x = v + dt * k3_v
        k4_v = self._compute_forces(x + dt * k3_x, v + dt * k3_v)

        self.positions = x + (dt / 6) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x)
        self.velocities = v + (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v)

        self._step_count += 1
        return self.get_state()

    def run(
        self,
        n_steps: int,
        record_interval: int = 1,
    ) -> list[dict[str, Any]]:
        """Run simulation for n_steps, recording at intervals.

        Args:
            n_steps: Number of integration steps.
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
        """Return current state dict."""
        speeds = np.linalg.norm(self.velocities, axis=1)
        headings = np.arctan2(self.velocities[:, 1], self.velocities[:, 0])
        return {
            "positions": self.positions.copy(),
            "velocities": self.velocities.copy(),
            "headings": headings,
            "speeds": speeds,
            "step": self._step_count,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata."""
        return {
            "model_family": "dorsogna_spp",
            "model_name": "D'Orsogna SPP with Morse potential (2006)",
            "n_particles": self.n_particles,
            "C_a": self.C_a,
            "C_r": self.C_r,
            "l_a": self.l_a,
            "l_r": self.l_r,
            "alpha": self.alpha,
            "beta": self.beta,
            "v_eq": self.v_eq,
            "dt": self.dt,
            "init_mode": self.init_mode,
            "init_radius": self.init_radius,
            "seed": self.seed,
            "space_type": "continuous_2d_open",
            "dynamics_type": "second_order_newtonian",
        }

    def get_timescale(self) -> dict[str, float]:
        """Return characteristic timescales.

        T_cross: estimated time to traverse the group diameter.
        T_orbit: estimated orbital period for milling (if applicable).
        """
        # Rough estimate: group diameter ~ 2*init_radius, speed ~ v_eq
        T_cross = 2 * self.init_radius / self.v_eq if self.v_eq > 0 else 100.0
        return {
            "T_cross": T_cross,
            "dt": self.dt,
            "steps_per_T_cross": T_cross / self.dt,
        }
