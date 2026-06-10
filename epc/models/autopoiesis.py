"""
Minimal SCL-style Autopoiesis Model.

Reference:
    Varela, F.G., Maturana, H.R., & Uribe, R. (1974).
    Autopoiesis: the organization of living systems, its characterization
    and a model. Biosystems, 5(4), 187–196.

    McMullin, B. & Varela, F.J. (1997). Rediscovering computational
    autopoiesis. In P. Husbands & I. Harvey (Eds.), Proc. ECAL 1997,
    pp. 38–47.

Model description:
    A 2D particle system with three types: Substrate (S), Catalyst (C),
    and Link (L). In a continuous periodic domain:

    - Substrate particles diffuse freely (Brownian).
    - Catalyst particles are stationary/slow-diffusing; they convert
      nearby substrate into link particles.
    - Link particles are attracted toward catalysts at a preferred radius
      (membrane_equilibrium_radius), forming a ring/shell.
    - Links also attract each other (chain formation).
    - Links decay stochastically back to substrate.
    - The membrane (ring of links around catalysts) is maintained by
      ongoing production from incoming substrate.

State dict keys:
    positions: (N, 2) float64 — particle positions
    types: (N,) int — 0=substrate, 1=catalyst, 2=link
    bonds: (N_links, 2) int — pairs of bonded link-particle indices
    step: int — current timestep
    n_particles: int — total particle count
    box_size: float — domain size
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Any, Optional


class AutopoiesisModel:
    """Minimal SCL autopoiesis model in 2D periodic domain.

    Catalyst particles convert nearby substrate into link particles.
    Links are attracted to catalysts at a preferred radius (forming a ring),
    and attract each other (chain formation). Links decay stochastically,
    requiring ongoing production to maintain the boundary (self-repair).

    Parameters
    ----------
    n_substrate : int
        Initial number of substrate particles.
    n_catalyst : int
        Number of catalyst particles (fixed).
    box_size : float
        Side length of periodic square domain.
    production_rate : float
        Probability per step that a substrate within production_radius
        of a catalyst converts to link.
    decay_rate : float
        Probability per step that a link decays back to substrate.
    production_radius : float
        Distance within which catalyst converts substrate to link.
    membrane_equilibrium_radius : float
        Preferred distance of links from nearest catalyst (ring radius).
    catalyst_link_attraction : float
        Strength of radial spring toward membrane_equilibrium_radius.
    link_attraction : float
        Strength of attraction between nearby links (chain cohesion).
    substrate_diffusion : float
        Diffusion coefficient for substrate particles.
    link_diffusion : float
        Diffusion coefficient for link particles (much lower).
    catalyst_diffusion : float
        Diffusion coefficient for catalyst particles (very low).
    repulsion_radius : float
        Short-range repulsion distance (steric).
    repulsion_strength : float
        Strength of short-range repulsion.
    dt : float
        Integration timestep.
    seed : int
        Random seed.
    """

    TYPE_SUBSTRATE = 0
    TYPE_CATALYST = 1
    TYPE_LINK = 2

    def __init__(
        self,
        n_substrate: int = 100,
        n_catalyst: int = 3,
        box_size: float = 20.0,
        production_rate: float = 0.15,
        decay_rate: float = 0.01,
        production_radius: float = 3.0,
        membrane_equilibrium_radius: float = 3.0,
        catalyst_link_attraction: float = 0.5,
        link_attraction: float = 0.3,
        substrate_diffusion: float = 0.3,
        link_diffusion: float = 0.01,
        catalyst_diffusion: float = 0.002,
        repulsion_radius: float = 0.4,
        repulsion_strength: float = 1.5,
        dt: float = 1.0,
        seed: int = 42,
    ) -> None:
        self.n_substrate = n_substrate
        self.n_catalyst = n_catalyst
        self.n_total = n_substrate + n_catalyst
        self.box_size = box_size
        self.production_rate = production_rate
        self.decay_rate = decay_rate
        self.production_radius = production_radius
        self.membrane_equilibrium_radius = membrane_equilibrium_radius
        self.catalyst_link_attraction = catalyst_link_attraction
        self.link_attraction = link_attraction
        self.substrate_diffusion = substrate_diffusion
        self.link_diffusion = link_diffusion
        self.catalyst_diffusion = catalyst_diffusion
        self.repulsion_radius = repulsion_radius
        self.repulsion_strength = repulsion_strength
        self.dt = dt
        self.seed = seed
        self.rng = np.random.default_rng(seed)

        self.positions: NDArray[np.float64] = np.zeros((self.n_total, 2))
        self.types: NDArray[np.int32] = np.zeros(self.n_total, dtype=np.int32)
        self._step_count = 0
        self._is_setup = False

    def setup(self) -> dict[str, Any]:
        """Initialize particle positions and types."""
        center = self.box_size / 2.0
        # Place catalyst particles near center
        self.positions[:self.n_catalyst] = center + self.rng.uniform(
            -1.0, 1.0, size=(self.n_catalyst, 2)
        )
        self.types[:self.n_catalyst] = self.TYPE_CATALYST

        # Place substrate particles randomly
        self.positions[self.n_catalyst:] = self.rng.uniform(
            0, self.box_size, size=(self.n_substrate, 2)
        )
        self.types[self.n_catalyst:] = self.TYPE_SUBSTRATE

        self._is_setup = True
        self._step_count = 0
        return self._get_state()

    def step(self) -> dict[str, Any]:
        """Advance one timestep."""
        self._step_count += 1
        self._do_production()
        self._do_decay()
        self._move_particles()
        return self._get_state()

    def _periodic_delta(self, p1: NDArray, p2: NDArray) -> NDArray:
        """Periodic displacement from p1 to p2."""
        delta = p2 - p1
        delta -= self.box_size * np.round(delta / self.box_size)
        return delta

    def _do_production(self) -> None:
        """Convert substrate particles near catalysts to links."""
        cat_idx = np.where(self.types == self.TYPE_CATALYST)[0]
        sub_idx = np.where(self.types == self.TYPE_SUBSTRATE)[0]
        if len(cat_idx) == 0 or len(sub_idx) == 0:
            return

        cat_pos = self.positions[cat_idx]
        for si in sub_idx:
            delta = self.positions[si] - cat_pos
            delta -= self.box_size * np.round(delta / self.box_size)
            dists = np.sqrt(np.sum(delta**2, axis=1))
            if np.min(dists) < self.production_radius:
                if self.rng.random() < self.production_rate:
                    self.types[si] = self.TYPE_LINK

    def _do_decay(self) -> None:
        """Links stochastically decay back to substrate."""
        link_idx = np.where(self.types == self.TYPE_LINK)[0]
        for li in link_idx:
            if self.rng.random() < self.decay_rate:
                self.types[li] = self.TYPE_SUBSTRATE

    def _move_particles(self) -> None:
        """Compute forces and update positions."""
        N = self.n_total
        forces = np.zeros((N, 2), dtype=np.float64)
        cat_idx = np.where(self.types == self.TYPE_CATALYST)[0]
        link_idx = np.where(self.types == self.TYPE_LINK)[0]

        # Catalyst-link radial spring: attract links to membrane_equilibrium_radius
        if len(cat_idx) > 0 and len(link_idx) > 0:
            cat_pos = self.positions[cat_idx]
            for li in link_idx:
                # Find nearest catalyst
                delta = self._periodic_delta(cat_pos, self.positions[li:li+1].repeat(len(cat_idx), axis=0))
                dists = np.sqrt(np.sum(delta**2, axis=1))
                nearest = np.argmin(dists)
                d = dists[nearest]
                if d > 0.01:
                    # Radial spring toward equilibrium radius
                    direction = delta[nearest] / d  # unit vector from catalyst to link
                    displacement = d - self.membrane_equilibrium_radius
                    f = -self.catalyst_link_attraction * displacement * direction
                    forces[li] += f

        # Link-link attraction (nearby links attract for chain cohesion)
        if len(link_idx) > 1:
            link_pos = self.positions[link_idx]
            for i in range(len(link_idx)):
                for j in range(i + 1, len(link_idx)):
                    delta = self._periodic_delta(link_pos[i], link_pos[j])
                    d = np.sqrt(np.sum(delta**2))
                    if 0.01 < d < self.membrane_equilibrium_radius:
                        f = self.link_attraction * delta / d
                        forces[link_idx[i]] += f
                        forces[link_idx[j]] -= f

        # Short-range repulsion (all pairs within repulsion_radius)
        # Only compute for nearby particles to keep it manageable
        for i in range(N):
            for j in range(i + 1, N):
                delta = self._periodic_delta(self.positions[i], self.positions[j])
                d = np.sqrt(np.sum(delta**2))
                if 0.01 < d < self.repulsion_radius:
                    f_mag = self.repulsion_strength * (self.repulsion_radius - d) / d
                    forces[i] -= f_mag * delta / d
                    forces[j] += f_mag * delta / d

        # Apply forces + diffusion
        for i in range(N):
            t = self.types[i]
            if t == self.TYPE_SUBSTRATE:
                D = self.substrate_diffusion
            elif t == self.TYPE_CATALYST:
                D = self.catalyst_diffusion
            else:
                D = self.link_diffusion

            noise = self.rng.normal(0, np.sqrt(2 * D * self.dt), size=2)
            self.positions[i] += forces[i] * self.dt + noise

        self.positions %= self.box_size

    def _get_state(self) -> dict[str, Any]:
        """Return current state snapshot."""
        # Compute bonds on the fly for state output
        link_idx = np.where(self.types == self.TYPE_LINK)[0]
        bonds = []
        if len(link_idx) > 1:
            link_pos = self.positions[link_idx]
            bond_r = self.membrane_equilibrium_radius * 0.6
            for i in range(len(link_idx)):
                for j in range(i + 1, len(link_idx)):
                    delta = self._periodic_delta(link_pos[i], link_pos[j])
                    if np.sqrt(np.sum(delta**2)) < bond_r:
                        bonds.append((int(link_idx[i]), int(link_idx[j])))

        return {
            'positions': self.positions.copy(),
            'types': self.types.copy(),
            'bonds': bonds,
            'step': self._step_count,
            'n_particles': self.n_total,
            'box_size': self.box_size,
        }

    def run(self, n_steps: int, record_every: int = 1) -> list[dict[str, Any]]:
        """Execute full simulation."""
        if not self._is_setup:
            initial = self.setup()
            history = [initial]
        else:
            history = [self._get_state()]

        for i in range(n_steps):
            state = self.step()
            if (i + 1) % record_every == 0:
                history.append(state)

        return history

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'autopoiesis_scl',
            'model_class': 'autopoiesis',
            'model_family': 'continuous_2d',
            'n_particles': self.n_total,
            'n_substrate_initial': self.n_substrate,
            'n_catalyst': self.n_catalyst,
            'box_size': self.box_size,
            'production_rate': self.production_rate,
            'decay_rate': self.decay_rate,
            'production_radius': self.production_radius,
            'membrane_equilibrium_radius': self.membrane_equilibrium_radius,
            'catalyst_link_attraction': self.catalyst_link_attraction,
            'link_attraction': self.link_attraction,
            'substrate_diffusion': self.substrate_diffusion,
            'link_diffusion': self.link_diffusion,
            'catalyst_diffusion': self.catalyst_diffusion,
            'repulsion_radius': self.repulsion_radius,
            'repulsion_strength': self.repulsion_strength,
            'dt': self.dt,
            'seed': self.seed,
            'space_type': 'continuous_2d_periodic',
            'dynamics_type': 'first_order_overdamped',
            'has_production': True,
            'has_decay': True,
            'has_bonding': True,
        }

    def get_timescale(self) -> float:
        """Return characteristic timescale (production time)."""
        return 1.0 / self.production_rate if self.production_rate > 0 else 100.0

    def breach_membrane(self, fraction: float = 0.1) -> None:
        """Remove a fraction of link particles (convert to substrate) for self-repair test."""
        link_idx = np.where(self.types == self.TYPE_LINK)[0]
        n_remove = max(1, int(len(link_idx) * fraction))
        if len(link_idx) == 0:
            return
        remove_idx = self.rng.choice(link_idx, size=min(n_remove, len(link_idx)),
                                     replace=False)
        for ri in remove_idx:
            self.types[ri] = self.TYPE_SUBSTRATE


class NonBondingParticleModel:
    """Control model: particles that do NOT form bonds or membranes.

    All particles diffuse freely with no production, decay, or bonding.
    Used as negative control for P30.
    """

    def __init__(
        self,
        n_particles: int = 100,
        box_size: float = 20.0,
        diffusion: float = 0.3,
        dt: float = 1.0,
        seed: int = 42,
    ) -> None:
        self.n_particles = n_particles
        self.box_size = box_size
        self.diffusion = diffusion
        self.dt = dt
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.positions = np.zeros((n_particles, 2))
        self.types = np.zeros(n_particles, dtype=np.int32)
        self._step_count = 0
        self._is_setup = False

    def setup(self) -> dict[str, Any]:
        self.positions = self.rng.uniform(0, self.box_size,
                                         size=(self.n_particles, 2))
        self.types[:] = 0
        self._is_setup = True
        self._step_count = 0
        return self._get_state()

    def step(self) -> dict[str, Any]:
        self._step_count += 1
        noise = self.rng.normal(0, np.sqrt(2 * self.diffusion * self.dt),
                                size=(self.n_particles, 2))
        self.positions += noise
        self.positions %= self.box_size
        return self._get_state()

    def _get_state(self) -> dict[str, Any]:
        return {
            'positions': self.positions.copy(),
            'types': self.types.copy(),
            'bonds': [],
            'step': self._step_count,
            'n_particles': self.n_particles,
            'box_size': self.box_size,
        }

    def run(self, n_steps: int, record_every: int = 1) -> list[dict[str, Any]]:
        if not self._is_setup:
            history = [self.setup()]
        else:
            history = [self._get_state()]
        for i in range(n_steps):
            state = self.step()
            if (i + 1) % record_every == 0:
                history.append(state)
        return history

    def get_metadata(self) -> dict[str, Any]:
        return {
            'model_name': 'non_bonding_particles',
            'model_class': 'brownian',
            'model_family': 'continuous_2d',
            'n_particles': self.n_particles,
            'box_size': self.box_size,
            'diffusion': self.diffusion,
            'dt': self.dt,
            'seed': self.seed,
            'space_type': 'continuous_2d_periodic',
            'has_production': False,
            'has_decay': False,
            'has_bonding': False,
        }


class DenseClusterModel:
    """Control model: particles that form a dense cluster but no closed boundary.

    Particles attract each other (like P1 aggregation) but do NOT form
    a topologically closed membrane with inside/outside distinction.
    Used as negative control for P30 vs P1.
    """

    def __init__(
        self,
        n_particles: int = 100,
        box_size: float = 20.0,
        attraction_strength: float = 0.3,
        attraction_radius: float = 3.0,
        diffusion: float = 0.1,
        repulsion_radius: float = 0.3,
        repulsion_strength: float = 1.0,
        dt: float = 1.0,
        seed: int = 42,
    ) -> None:
        self.n_particles = n_particles
        self.box_size = box_size
        self.attraction_strength = attraction_strength
        self.attraction_radius = attraction_radius
        self.diffusion = diffusion
        self.repulsion_radius = repulsion_radius
        self.repulsion_strength = repulsion_strength
        self.dt = dt
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.positions = np.zeros((n_particles, 2))
        self.types = np.zeros(n_particles, dtype=np.int32)
        self._step_count = 0
        self._is_setup = False

    def setup(self) -> dict[str, Any]:
        self.positions = self.rng.uniform(0, self.box_size,
                                         size=(self.n_particles, 2))
        self.types[:] = 0
        self._is_setup = True
        self._step_count = 0
        return self._get_state()

    def step(self) -> dict[str, Any]:
        self._step_count += 1
        N = self.n_particles
        forces = np.zeros((N, 2))

        for i in range(N):
            for j in range(i + 1, N):
                delta = self.positions[j] - self.positions[i]
                delta -= self.box_size * np.round(delta / self.box_size)
                dist = np.sqrt(np.sum(delta**2))
                if dist < 0.01:
                    continue
                if dist < self.repulsion_radius:
                    f_mag = self.repulsion_strength * (self.repulsion_radius - dist)
                    forces[i] -= f_mag * delta / dist
                    forces[j] += f_mag * delta / dist
                elif dist < self.attraction_radius:
                    f_mag = self.attraction_strength
                    forces[i] += f_mag * delta / dist
                    forces[j] -= f_mag * delta / dist

        noise = self.rng.normal(0, np.sqrt(2 * self.diffusion * self.dt),
                                size=(N, 2))
        self.positions += forces * self.dt + noise
        self.positions %= self.box_size
        return self._get_state()

    def _get_state(self) -> dict[str, Any]:
        return {
            'positions': self.positions.copy(),
            'types': self.types.copy(),
            'bonds': [],
            'step': self._step_count,
            'n_particles': self.n_particles,
            'box_size': self.box_size,
        }

    def run(self, n_steps: int, record_every: int = 1) -> list[dict[str, Any]]:
        if not self._is_setup:
            history = [self.setup()]
        else:
            history = [self._get_state()]
        for i in range(n_steps):
            state = self.step()
            if (i + 1) % record_every == 0:
                history.append(state)
        return history

    def get_metadata(self) -> dict[str, Any]:
        return {
            'model_name': 'dense_cluster',
            'model_class': 'aggregation',
            'model_family': 'continuous_2d',
            'n_particles': self.n_particles,
            'box_size': self.box_size,
            'attraction_strength': self.attraction_strength,
            'attraction_radius': self.attraction_radius,
            'diffusion': self.diffusion,
            'dt': self.dt,
            'seed': self.seed,
            'space_type': 'continuous_2d_periodic',
            'has_production': False,
            'has_decay': False,
            'has_bonding': False,
        }
