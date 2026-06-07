"""
Counterflow Social-Force Model — Lane Formation (P7).

Reference:
    Helbing, D. & Molnár, P. (1995). Social force model for pedestrian
    dynamics. Physical Review E, 51(5), 4282–4286.
    DOI: 10.1103/PhysRevE.51.4282

Additional reference for lane order parameter:
    Nowak, S. & Schadschneider, A. (2012). Quantitative analysis of
    pedestrian counterflow in a cellular automaton model. Physical Review
    E, 85(6), 066128.

Model description:
    Two populations of agents (N/2 each) move in opposite directions (+x
    and -x) through a 2D corridor with periodic boundaries along the flow
    axis (x) and reflecting walls on the lateral axis (y). Agents
    experience a driving force toward their desired velocity and a
    repulsive social force from nearby agents.

Update rule (first-order overdamped, Euler integration):
    v_i(t+dt) = v_desired_i + F_social_i / tau
    x_i(t+dt) = x_i(t) + v_i(t+dt) * dt

    F_social_i = sum_{j != i} A * exp(-d_ij / B) * n_ij

where d_ij is the distance between agents i and j, n_ij is the unit
vector pointing from j to i, A is the repulsion amplitude, B is the
repulsion range.

Key parameters (Helbing-Molnár-style):
    N = 200, corridor W x H = 20.0 x 4.0, v0 = 1.0, A = 5.0, B = 0.3,
    dt = 0.05, tau = 0.5 (relaxation time)

Validation targets:
    1. Lane order parameter (Nowak-Schadschneider) rises from ~0 to ~0.8+
    2. Head-on encounter rate decreases as lanes form
    3. Throughput (mean speed along desired direction) increases vs mixed
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Any, Optional


class LaneFormationModel:
    """Counterflow social-force model for lane formation (P7).

    Two populations with opposing desired directions (+x, -x) in a
    periodic corridor. Social force repulsion drives spontaneous lateral
    segregation into same-direction lanes.

    State dict keys:
        positions: (N, 2) float64 — agent positions
        velocities: (N, 2) float64 — agent velocities
        labels: (N,) int — population label (0 = +x, 1 = -x)
        step: int — current timestep

    Init modes:
        random: uniform random positions, velocities initialized to desired
        mixed: alternating population positions (maximizes initial mixing)
    """

    def __init__(
        self,
        n_agents: int = 200,
        corridor_width: float = 20.0,
        corridor_height: float = 4.0,
        desired_speed: float = 1.0,
        repulsion_amplitude: float = 5.0,
        repulsion_range: float = 0.3,
        tau: float = 0.5,
        dt: float = 0.05,
        wall_repulsion_amplitude: float = 10.0,
        wall_repulsion_range: float = 0.2,
        interaction_radius: float = 2.0,
        init_mode: str = "random",
        seed: Optional[int] = None,
    ):
        """Initialize lane formation model.

        Args:
            n_agents: Total number of agents (split evenly into two populations).
            corridor_width: Length along flow axis (periodic).
            corridor_height: Width across flow axis (reflecting walls at y=0, y=H).
            desired_speed: Magnitude of desired velocity |v0|.
            repulsion_amplitude: Social force strength A.
            repulsion_range: Social force decay length B.
            tau: Relaxation time for velocity adaptation.
            dt: Integration timestep.
            wall_repulsion_amplitude: Wall repulsion strength.
            wall_repulsion_range: Wall repulsion decay length.
            interaction_radius: Cutoff for social force computation.
            init_mode: One of 'random', 'mixed'.
            seed: Random seed for reproducibility.
        """
        self.n_agents = n_agents
        self.corridor_width = corridor_width
        self.corridor_height = corridor_height
        self.desired_speed = desired_speed
        self.repulsion_amplitude = repulsion_amplitude
        self.repulsion_range = repulsion_range
        self.tau = tau
        self.dt = dt
        self.wall_repulsion_amplitude = wall_repulsion_amplitude
        self.wall_repulsion_range = wall_repulsion_range
        self.interaction_radius = interaction_radius
        self.init_mode = init_mode
        self.seed = seed

        # State
        self.positions: NDArray[np.float64] = np.empty(0)
        self.velocities: NDArray[np.float64] = np.empty(0)
        self.labels: NDArray[np.int64] = np.empty(0, dtype=np.int64)
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self._step_count: int = 0

    def setup(self) -> dict[str, Any]:
        """Initialize agent positions, velocities, and population labels.

        Returns:
            Initial state dict.
        """
        self.rng = np.random.default_rng(self.seed)
        self._step_count = 0
        N = self.n_agents
        W = self.corridor_width
        H = self.corridor_height

        # Population labels: first half = +x (label 0), second half = -x (label 1)
        self.labels = np.zeros(N, dtype=np.int64)
        self.labels[N // 2:] = 1

        if self.init_mode == "random":
            self.positions = np.column_stack([
                self.rng.uniform(0, W, size=N),
                self.rng.uniform(0, H, size=N),
            ])
        elif self.init_mode == "mixed":
            # Place agents in a grid-like pattern for maximum initial mixing
            self.positions = np.column_stack([
                self.rng.uniform(0, W, size=N),
                self.rng.uniform(0, H, size=N),
            ])
            # Shuffle labels to ensure mixing (labels already assigned by index)
            self.rng.shuffle(self.labels)
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode!r}")

        # Initialize velocities to desired direction
        self.velocities = np.zeros((N, 2), dtype=np.float64)
        self.velocities[:, 0] = np.where(
            self.labels == 0, self.desired_speed, -self.desired_speed
        )

        return self.get_state()

    def _desired_velocity(self) -> NDArray[np.float64]:
        """Compute desired velocity for each agent based on population label."""
        v_desired = np.zeros((self.n_agents, 2), dtype=np.float64)
        v_desired[:, 0] = np.where(
            self.labels == 0, self.desired_speed, -self.desired_speed
        )
        return v_desired

    def _compute_social_force(self) -> NDArray[np.float64]:
        """Compute pairwise social repulsion forces.

        Uses periodic boundary in x (flow axis) and direct distance in y.
        Cutoff at interaction_radius for efficiency.
        """
        N = self.n_agents
        W = self.corridor_width
        r_cut = self.interaction_radius
        A = self.repulsion_amplitude
        B = self.repulsion_range

        forces = np.zeros((N, 2), dtype=np.float64)

        # Compute pairwise displacement with periodic x
        # dx_ij = x_i - x_j (with periodic wrapping in x)
        pos = self.positions
        for i in range(N):
            dx = pos[i, 0] - pos[:, 0]
            dy = pos[i, 1] - pos[:, 1]

            # Periodic boundary in x
            dx = dx - W * np.round(dx / W)

            dist = np.sqrt(dx**2 + dy**2)
            dist[i] = np.inf  # exclude self

            # Apply cutoff
            mask = dist < r_cut
            if not np.any(mask):
                continue

            d = dist[mask]
            # Unit vector from j to i
            nx = dx[mask] / d
            ny = dy[mask] / d

            # Social force magnitude: A * exp(-d / B)
            f_mag = A * np.exp(-d / B)

            forces[i, 0] += np.sum(f_mag * nx)
            forces[i, 1] += np.sum(f_mag * ny)

        return forces

    def _compute_social_force_vectorized(self) -> NDArray[np.float64]:
        """Vectorized social force computation for better performance."""
        N = self.n_agents
        W = self.corridor_width
        r_cut = self.interaction_radius
        A = self.repulsion_amplitude
        B = self.repulsion_range

        pos = self.positions

        # Pairwise displacements (i - j)
        dx = pos[:, 0:1] - pos[:, 0]  # (N, N)
        dy = pos[:, 1:2] - pos[:, 1]  # (N, N)

        # Periodic boundary in x
        dx = dx - W * np.round(dx / W)

        dist = np.sqrt(dx**2 + dy**2)
        np.fill_diagonal(dist, np.inf)

        # Cutoff mask
        mask = dist < r_cut

        # Unit vectors (safe division)
        safe_dist = np.where(mask, dist, 1.0)
        nx = np.where(mask, dx / safe_dist, 0.0)
        ny = np.where(mask, dy / safe_dist, 0.0)

        # Force magnitudes
        f_mag = np.where(mask, A * np.exp(-dist / B), 0.0)

        # Sum forces on each agent
        forces = np.column_stack([
            np.sum(f_mag * nx, axis=1),
            np.sum(f_mag * ny, axis=1),
        ])

        return forces

    def _compute_wall_force(self) -> NDArray[np.float64]:
        """Compute repulsive force from corridor walls (y=0 and y=H)."""
        N = self.n_agents
        H = self.corridor_height
        A_w = self.wall_repulsion_amplitude
        B_w = self.wall_repulsion_range

        forces = np.zeros((N, 2), dtype=np.float64)

        # Distance from bottom wall (y=0)
        d_bottom = self.positions[:, 1]
        f_bottom = A_w * np.exp(-d_bottom / B_w)
        forces[:, 1] += f_bottom

        # Distance from top wall (y=H)
        d_top = H - self.positions[:, 1]
        f_top = A_w * np.exp(-d_top / B_w)
        forces[:, 1] -= f_top

        return forces

    def step(self) -> dict[str, Any]:
        """Execute one integration step.

        Updates velocities via driving + social + wall forces, then updates
        positions with periodic x and clamped y.

        Returns:
            State dict after this step.
        """
        W = self.corridor_width
        H = self.corridor_height

        # Compute forces
        v_desired = self._desired_velocity()
        f_social = self._compute_social_force_vectorized()
        f_wall = self._compute_wall_force()

        # Velocity update: v = v_desired + (f_social + f_wall) / tau
        # (overdamped first-order: velocity relaxes toward desired + force/tau)
        self.velocities = v_desired + (f_social + f_wall) * (self.dt / self.tau)

        # Position update
        self.positions += self.velocities * self.dt

        # Periodic boundary in x
        self.positions[:, 0] %= W

        # Reflecting walls in y: clamp to [0, H]
        self.positions[:, 1] = np.clip(self.positions[:, 1], 0.0, H)

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
            Dict with keys: positions, velocities, labels, step.
        """
        return {
            "positions": self.positions.copy(),
            "velocities": self.velocities.copy(),
            "labels": self.labels.copy(),
            "step": self._step_count,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata for detector use."""
        return {
            "model_family": "lane_formation",
            "model_name": "Counterflow Social-Force (Helbing-Molnár 1995)",
            "model_class": "counterflow",
            "n_agents": self.n_agents,
            "corridor_width": self.corridor_width,
            "corridor_height": self.corridor_height,
            "desired_speed": self.desired_speed,
            "repulsion_amplitude": self.repulsion_amplitude,
            "repulsion_range": self.repulsion_range,
            "tau": self.tau,
            "dt": self.dt,
            "interaction_radius": self.interaction_radius,
            "density": self.n_agents / (self.corridor_width * self.corridor_height),
            "substrate_type": "continuous_2d",
            "has_counterflow": True,
            "has_alignment_rule": False,
            "has_attraction_rule": False,
            "n_populations": 2,
        }
