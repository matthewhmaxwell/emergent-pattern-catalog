"""
Collective Sensing Model — Distributed Gradient Detection (P17).

Reference:
    Berdahl, A., Torney, C. J., Ioannou, C. C., Faria, J. J., & Couzin, I. D.
    (2013). Emergent sensing of complex environments by mobile animal groups.
    Science, 339(6119), 574–576. DOI: 10.1126/science.1225883

    Camley, B. A., Zimmermann, J., Levine, H., & Rappel, W.-J. (2016).
    Collective gradient sensing and chemotaxis: modeling and recent developments.
    J. Phys.: Condens. Matter, 28, 253001.

Model description:
    Agents perform correlated random walks in a 2D periodic domain containing
    a scalar environmental field (Gaussian resource peak). The collective
    sensing mechanism from Berdahl (2013):

    1. Each agent senses the local field with additive Gaussian noise.
    2. Agents that sense higher quality SLOW DOWN (speed modulation).
    3. Social coupling: agents are attracted toward the group centroid.

    The speed-modulation mechanism creates a bias: agents near the peak spend
    more time there (slower), so the group centroid drifts toward the peak.
    For a single agent, the noisy sensing means it cannot reliably detect
    the gradient (signal is buried in noise). For a group, independent noise
    samples are effectively averaged through the centroid attraction, giving
    the group much better gradient-following ability.

    The defining P17 signature: group chemotactic index (CI) RISES with group
    size N, while isolated individuals (N=1) remain near chance.

Update rule:
    1. Sense: q_i = field(x_i) + N(0, σ_sense)
    2. Speed: v_i = v_max * (1 - α * clamp(q_i / A, 0, 1))
    3. Heading: θ_i += social_torque_i + N(0, σ_turn)
    4. Position: x_i += v_i * [cos(θ_i), sin(θ_i)] * dt (periodic BC)

Key insight for N-scaling:
    The social attraction toward centroid means each agent's EFFECTIVE field
    estimate is a weighted average of all agents' estimates. For N agents with
    iid noise of variance σ², the effective sensing noise scales as σ/√N.
    This is why CI rises with N — the signal-to-noise ratio improves with
    group size.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Any, Optional


class CollectiveSensingModel:
    """Berdahl-style collective gradient sensing model (P17).

    Agents in a periodic 2D domain sense a Gaussian resource field with noise,
    slow down in high-field regions, and are attracted toward the group centroid.
    The collective speed-modulation mechanism produces emergent gradient climbing
    that improves with group size.

    State dict keys:
        positions: (N, 2) float64 — agent positions
        headings: (N,) float64 — heading angles in radians
        speeds: (N,) float64 — current speed of each agent
        field_samples: (N,) float64 — noisy field value sensed by each agent
        step: int — current timestep
    """

    def __init__(
        self,
        n_agents: int = 50,
        box_size: float = 20.0,
        v_max: float = 0.3,
        turn_noise: float = 0.3,
        sensing_noise: float = 1.5,
        alpha: float = 0.9,
        social_strength: float = 0.15,
        field_center: Optional[tuple[float, float]] = None,
        field_sigma: float = 4.0,
        field_amplitude: float = 1.0,
        dt: float = 1.0,
        init_mode: str = "offset",
        seed: Optional[int] = None,
    ):
        """Initialize collective sensing model.

        Args:
            n_agents: Number of agents in the group.
            box_size: Side length of square periodic domain.
            v_max: Maximum speed (at zero field).
            turn_noise: Std of heading noise per step (radians).
            sensing_noise: Std of additive Gaussian noise on field sensing.
            alpha: Speed modulation depth (0=no modulation, 1=stop at peak).
            social_strength: Strength of attraction toward group centroid (0-1).
            field_center: (x, y) center of the Gaussian field.
            field_sigma: Gaussian field width.
            field_amplitude: Peak value of the field.
            dt: Integration timestep.
            init_mode: 'random', 'offset' (clustered away from peak).
            seed: Random seed.
        """
        self.n_agents = n_agents
        self.box_size = box_size
        self.v_max = v_max
        self.turn_noise = turn_noise
        self.sensing_noise = sensing_noise
        self.alpha = alpha
        self.social_strength = social_strength
        self.field_center = field_center if field_center is not None else (box_size / 2, box_size / 2)
        self.field_sigma = field_sigma
        self.field_amplitude = field_amplitude
        self.dt = dt
        self.init_mode = init_mode
        self.seed = seed

        # State
        self.positions: NDArray[np.float64] = np.empty(0)
        self.headings: NDArray[np.float64] = np.empty(0)
        self.speeds: NDArray[np.float64] = np.empty(0)
        self.field_samples: NDArray[np.float64] = np.empty(0)
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self._step_count: int = 0

    def _evaluate_field(self, positions: NDArray[np.float64]) -> NDArray[np.float64]:
        """Evaluate the scalar field at positions (periodic distance to center)."""
        cx, cy = self.field_center
        L = self.box_size
        dx = positions[:, 0] - cx
        dy = positions[:, 1] - cy
        dx = dx - L * np.round(dx / L)
        dy = dy - L * np.round(dy / L)
        r_sq = dx**2 + dy**2
        return self.field_amplitude * np.exp(-r_sq / (2 * self.field_sigma**2))

    def _periodic_com(self) -> NDArray[np.float64]:
        """Compute center of mass using circular mean (handles periodicity)."""
        L = self.box_size
        theta_x = 2 * np.pi * self.positions[:, 0] / L
        theta_y = 2 * np.pi * self.positions[:, 1] / L
        com_x = L * np.arctan2(np.mean(np.sin(theta_x)), np.mean(np.cos(theta_x))) / (2 * np.pi)
        com_y = L * np.arctan2(np.mean(np.sin(theta_y)), np.mean(np.cos(theta_y))) / (2 * np.pi)
        return np.array([(com_x % L), (com_y % L)])

    def setup(self) -> dict[str, Any]:
        """Initialize agent positions and headings."""
        self.rng = np.random.default_rng(self.seed)
        self._step_count = 0
        N = self.n_agents
        L = self.box_size

        if self.init_mode == "random":
            self.positions = np.column_stack([
                self.rng.uniform(0, L, size=N),
                self.rng.uniform(0, L, size=N),
            ])
        elif self.init_mode == "offset":
            # Start agents clustered ~L/3 away from peak
            cx, cy = self.field_center
            start_x = (cx + L * 0.33) % L
            start_y = cy  # same y as peak to make gradient clear
            spread = min(2.0, L * 0.1)
            self.positions = np.column_stack([
                (self.rng.normal(start_x, spread, size=N)) % L,
                (self.rng.normal(start_y, spread, size=N)) % L,
            ])
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode!r}")

        self.headings = self.rng.uniform(0, 2 * np.pi, size=N)
        self.speeds = np.full(N, self.v_max)
        self.field_samples = self._evaluate_field(self.positions)

        return self.get_state()

    def step(self) -> dict[str, Any]:
        """Advance simulation by one timestep."""
        N = self.n_agents
        L = self.box_size

        # 1. Sense local field + noise
        true_field = self._evaluate_field(self.positions)
        noise = self.rng.normal(0, self.sensing_noise, size=N)
        sensed = true_field + noise
        self.field_samples = sensed

        # 2. Speed modulation: slow down in high-quality regions
        # Clamp sensed to [0, amplitude] before normalizing
        q_norm = np.clip(sensed / self.field_amplitude, 0.0, 1.0)
        self.speeds = self.v_max * (1.0 - self.alpha * q_norm)

        # 3. Social attraction toward group centroid
        social_torque = np.zeros(N)
        if N > 1:
            com = self._periodic_com()
            # Periodic vector from each agent to centroid
            to_com_x = com[0] - self.positions[:, 0]
            to_com_y = com[1] - self.positions[:, 1]
            to_com_x = to_com_x - L * np.round(to_com_x / L)
            to_com_y = to_com_y - L * np.round(to_com_y / L)

            # Angle toward centroid
            target_heading = np.arctan2(to_com_y, to_com_x)
            angle_diff = target_heading - self.headings
            angle_diff = (angle_diff + np.pi) % (2 * np.pi) - np.pi
            social_torque = self.social_strength * angle_diff

        # 4. Update headings with noise + social torque
        heading_noise = self.rng.normal(0, self.turn_noise, size=N)
        self.headings = (self.headings + heading_noise + social_torque) % (2 * np.pi)

        # 5. Move
        vx = self.speeds * np.cos(self.headings)
        vy = self.speeds * np.sin(self.headings)
        self.positions[:, 0] = (self.positions[:, 0] + vx * self.dt) % L
        self.positions[:, 1] = (self.positions[:, 1] + vy * self.dt) % L

        self._step_count += 1
        return self.get_state()

    def get_state(self) -> dict[str, Any]:
        """Return current state snapshot."""
        return {
            "positions": self.positions.copy(),
            "headings": self.headings.copy(),
            "speeds": self.speeds.copy(),
            "field_samples": self.field_samples.copy(),
            "step": self._step_count,
        }

    def run(self, n_steps: int = 500, record_interval: int = 1) -> list[dict[str, Any]]:
        """Run simulation for n_steps, recording every record_interval steps.

        Args:
            n_steps: Total number of simulation steps.
            record_interval: Record state every this many steps.

        Returns:
            List of state dicts.
        """
        history = []
        state = self.setup()
        history.append(state)

        for t in range(1, n_steps + 1):
            state = self.step()
            if t % record_interval == 0:
                history.append(state)

        return history

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata for detector dispatch."""
        return {
            "model_name": "collective_sensing",
            "model_class": "collective_sensing",
            "n_agents": self.n_agents,
            "box_size": self.box_size,
            "v_max": self.v_max,
            "turn_noise": self.turn_noise,
            "sensing_noise": self.sensing_noise,
            "alpha": self.alpha,
            "social_strength": self.social_strength,
            "field_center": self.field_center,
            "field_sigma": self.field_sigma,
            "field_amplitude": self.field_amplitude,
            "dt": self.dt,
            "substrate_type": "continuous_2d",
            "has_alignment_rule": False,
            "has_attraction_rule": True,
            "has_field_sensing": True,
            "has_speed_modulation": True,
            "has_turning_modulation": False,
        }


def compute_chemotactic_index(
    history: list[dict[str, Any]],
    field_center: tuple[float, float],
    box_size: float,
) -> float:
    """Compute chemotactic index: net approach to field peak.

    CI = (d_initial - d_final) / d_initial where d is the periodic distance
    from the group center of mass to the field peak.

    Args:
        history: List of state dicts with 'positions' key.
        field_center: (x, y) location of the field peak.
        box_size: Domain size for periodic distance computation.

    Returns:
        Chemotactic index in [-1, 1]. Positive = approaching peak.
    """
    if len(history) < 2:
        return 0.0

    cx, cy = field_center
    L = box_size

    def periodic_com(positions: NDArray) -> tuple[float, float]:
        """Center of mass via circular mean."""
        theta_x = 2 * np.pi * positions[:, 0] / L
        theta_y = 2 * np.pi * positions[:, 1] / L
        com_x = L * np.arctan2(np.mean(np.sin(theta_x)), np.mean(np.cos(theta_x))) / (2 * np.pi)
        com_y = L * np.arctan2(np.mean(np.sin(theta_y)), np.mean(np.cos(theta_y))) / (2 * np.pi)
        return (com_x % L, com_y % L)

    def periodic_dist(x1: float, y1: float, x2: float, y2: float) -> float:
        dx = x1 - x2
        dy = y1 - y2
        dx = dx - L * round(dx / L)
        dy = dy - L * round(dy / L)
        return float(np.sqrt(dx**2 + dy**2))

    # Initial distance (first 10% of history)
    early_end = max(2, len(history) // 10)
    early_dists = []
    for s in history[:early_end]:
        com = periodic_com(s["positions"])
        early_dists.append(periodic_dist(com[0], com[1], cx, cy))
    d_init = float(np.mean(early_dists))

    # Final distance (last 25% of history)
    late_start = len(history) * 3 // 4
    late_dists = []
    for s in history[late_start:]:
        com = periodic_com(s["positions"])
        late_dists.append(periodic_dist(com[0], com[1], cx, cy))
    d_final = float(np.mean(late_dists))

    if d_init < 1e-8:
        return 1.0 if d_final < 1e-8 else 0.0

    ci = (d_init - d_final) / d_init
    return float(np.clip(ci, -1.0, 1.0))
