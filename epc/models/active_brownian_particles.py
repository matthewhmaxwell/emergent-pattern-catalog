"""
Active Brownian Particles (ABP) with density-dependent self-propulsion.

Reference:
    Fily, Y. & Marchetti, M. C. (2012).
    Athermal Phase Separation of Self-Propelled Particles with No Alignment.
    Physical Review Letters, 108(23), 235702.
    DOI: 10.1103/PhysRevLett.108.235702

Equations of motion (overdamped Langevin, no translational noise):

    dr_i/dt = v(rho_i) * e_i(theta_i)
    dtheta_i/dt = sqrt(2 D_r) * xi_i(t)

where e_i = (cos theta_i, sin theta_i), xi_i is unit Gaussian white noise,
and the density-dependent speed is the Fily-Marchetti linear law

    v(rho) = v_0 * max(0, 1 - rho / rho_star).

The local density rho_i at particle i is estimated by counting neighbors
inside a coarse-graining disk of radius r_cg (default sigma = 1.0, the
particle diameter) and dividing by the disk area. This matches Fily &
Marchetti's measurement-scale-agnostic mean-field ansatz; the detector
computes rho_i independently on its own r_cg.

Periodic boundary conditions on [0, L) x [0, L).

KEY CONTRAST WITH VICSEK AND D'ORSOGNA (continuous_2d neighbors):

  Vicsek (1995): ALIGNMENT rule. Particles average neighbor headings.
    Constant speed v_0. No density-dependent speed. Ordered state shows
    bulk polarization, NOT density bimodality.

  D'Orsogna (2006): ATTRACTION + REPULSION rule (Morse potential). Speed
    is NOT constant (second-order Newtonian), but clustering here is
    driven by the explicit -grad U_att pairwise force, not by kinetic
    slowdown.

  ABP (Fily-Marchetti 2012): NEITHER alignment NOR attraction. The only
    inter-particle effect is the density-dependent speed reduction
    v(rho). Clustering emerges purely from the kinetic mechanism:
    particles slow where density is already high, which feeds back to
    raise density further (positive-feedback amplification of density
    fluctuations). This is the defining character of pattern P2 (MIPS).

The get_metadata() output includes explicit boolean flags
``has_alignment_rule = False`` and ``has_attraction_rule = False``
together with ``has_density_dependent_speed = True`` so that detectors
can assert the mechanistic absence-of-attraction required for a DEFINITIVE
P2 call. See Architecture Decision 43.

Observables exposed in state snapshots:
  - positions: (N, 2) float64
  - velocities: (N, 2) float64  (= speeds * [cos theta, sin theta])
  - headings: (N,) float64 in [-pi, pi)
  - speeds: (N,) float64  (magnitude; VARIABLE, unlike Vicsek)
  - step: int

Canonical replication parameters (Fily-Marchetti 2012 Fig 1 regime):
  N = 400, box_size L such that phi = N * pi * sigma^2 / 4 / L^2 = 0.6
  (packing fraction), D_r = 1e-3, v_0 = 2.0e-1 to 5.0e-1 (Pe = v_0 /
  (D_r * sigma); our native unit sigma = 1 makes Pe = v_0 / D_r).

  Phase-separated regime: Pe >= 50 and packing fraction phi in [0.3, 0.7].

API surface (same shape as Vicsek, NOT BaseModel subclass, matching
Vicsek's standalone pattern).
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import cKDTree


class ActiveBrownianParticles:
    """Fily-Marchetti (2012) active Brownian particles with v(rho) law.

    Parameters
    ----------
    n_particles : int
        Number of self-propelled particles (N).
    box_size : float
        Side length of periodic square domain (L). Particle diameter
        sigma is fixed at 1.0; the control parameter is the packing
        fraction phi = N * pi * sigma^2 / 4 / L^2.
    v0 : float
        Reference self-propulsion speed (v_0).
    rho_star : float
        Critical density at which v -> 0 in the Fily-Marchetti linear
        slowdown law. Default 10.0 per unit area. Matches the mean-field
        ansatz v(rho) = v_0 * (1 - rho / rho_star).
    D_r : float
        Rotational diffusion constant. Persistence length L_p = v_0 / D_r.
        Peclet number Pe = v_0 / (D_r * sigma) = v_0 / D_r at sigma = 1.
    r_cg : float
        Coarse-graining radius for local density estimation. Default 1.0
        (particle diameter). Density at particle i is estimated as
        (number of neighbors within r_cg, INCLUDING self) /
        (pi * r_cg^2).
    dt : float
        Integration timestep.
    init_mode : str
        One of 'uniform', 'lattice', 'clump'. 'uniform' = random uniform
        positions, random headings; 'lattice' = triangular lattice with
        random headings (clean Vicsek-style IC); 'clump' = all particles
        in central L/4 box, random headings (stress test for mixing).
    seed : int
        Random seed for reproducibility.
    """

    def __init__(
        self,
        n_particles: int = 400,
        box_size: float = 20.0,
        v0: float = 0.3,
        rho_star: float = 10.0,
        D_r: float = 3e-3,
        r_cg: float = 1.0,
        dt: float = 0.1,
        init_mode: str = "uniform",
        seed: Optional[int] = None,
    ) -> None:
        if n_particles < 10:
            raise ValueError("n_particles must be >= 10 for meaningful density statistics")
        if box_size <= 0:
            raise ValueError("box_size must be positive")
        if v0 < 0:
            raise ValueError("v0 must be non-negative")
        if rho_star <= 0:
            raise ValueError("rho_star must be positive")
        if D_r < 0:
            raise ValueError("D_r must be non-negative")
        if r_cg <= 0:
            raise ValueError("r_cg must be positive")
        if dt <= 0:
            raise ValueError("dt must be positive")
        if init_mode not in ("uniform", "lattice", "clump"):
            raise ValueError(f"unknown init_mode: {init_mode!r}")

        self.n_particles = n_particles
        self.box_size = box_size
        self.v0 = v0
        self.rho_star = rho_star
        self.D_r = D_r
        self.r_cg = r_cg
        self.dt = dt
        self.init_mode = init_mode
        self.seed = seed

        # Derived
        # "sigma" = particle diameter is fixed at 1.0; packing fraction
        # phi = N * pi / 4 / L^2
        self.sigma = 1.0
        self.packing_fraction = (
            n_particles * np.pi * self.sigma**2 / 4.0 / box_size**2
        )
        # Peclet number per Fily-Marchetti convention
        if D_r > 0:
            self.peclet = v0 / (D_r * self.sigma)
        else:
            self.peclet = float("inf")

        # State
        self.positions: NDArray[np.float64] = np.empty((0, 2))
        self.headings: NDArray[np.float64] = np.empty(0)
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self._step_count: int = 0
        self._is_setup: bool = False

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def setup(self) -> dict[str, Any]:
        """Initialize positions and headings; return initial state."""
        self.rng = np.random.default_rng(self.seed)
        self._step_count = 0
        N = self.n_particles
        L = self.box_size

        if self.init_mode == "uniform":
            self.positions = self.rng.uniform(0, L, size=(N, 2))
        elif self.init_mode == "lattice":
            # Triangular lattice; pick nearest grid that fits N particles
            nx = int(np.ceil(np.sqrt(N)))
            ny = int(np.ceil(N / nx))
            ax = L / nx
            ay = L / ny
            xs = np.arange(nx) * ax + ax / 2
            ys = np.arange(ny) * ay + ay / 2
            XX, YY = np.meshgrid(xs, ys, indexing="ij")
            # Stagger every other row for triangular lattice
            XX[:, 1::2] += ax / 2
            XX %= L
            pts = np.column_stack([XX.ravel(), YY.ravel()])
            self.positions = pts[:N]
        elif self.init_mode == "clump":
            center = L / 2
            spread = L / 4
            self.positions = self.rng.uniform(
                center - spread / 2, center + spread / 2, size=(N, 2)
            )
            self.positions %= L
        else:  # pragma: no cover - guarded in __init__
            raise ValueError(f"unknown init_mode: {self.init_mode!r}")

        self.headings = self.rng.uniform(-np.pi, np.pi, size=N)
        self._is_setup = True
        return self.get_state()

    def step(self) -> dict[str, Any]:
        """Advance one dt under Fily-Marchetti ABP dynamics.

        Order of operations (Euler-Maruyama):
          1. Compute local densities rho_i at current positions.
          2. Compute speeds v_i = v0 * max(0, 1 - rho_i / rho_star).
          3. Update positions: r_i <- r_i + v_i * e(theta_i) * dt.
          4. Update headings: theta_i <- theta_i + sqrt(2 D_r dt) * eta,
             eta ~ N(0, 1).
          5. Wrap positions into [0, L) (periodic).
        """
        N = self.n_particles
        L = self.box_size

        # --- Local density via cKDTree with periodic BC ---
        tree = cKDTree(self.positions, boxsize=L)
        # count_neighbors with eps=r_cg INCLUDES self (distance 0)
        counts = tree.query_ball_point(self.positions, r=self.r_cg, return_length=True)
        counts = np.asarray(counts, dtype=np.float64)
        area = np.pi * self.r_cg**2
        rho_local = counts / area  # (N,)

        # --- Density-dependent speed (Fily-Marchetti linear law) ---
        v_local = self.v0 * np.maximum(0.0, 1.0 - rho_local / self.rho_star)

        # --- Propulsion step ---
        cos_t = np.cos(self.headings)
        sin_t = np.sin(self.headings)
        self.positions[:, 0] += v_local * cos_t * self.dt
        self.positions[:, 1] += v_local * sin_t * self.dt

        # --- Rotational diffusion step ---
        eta = self.rng.standard_normal(N)
        self.headings += np.sqrt(2.0 * self.D_r * self.dt) * eta
        # Wrap to [-pi, pi)
        self.headings = (self.headings + np.pi) % (2 * np.pi) - np.pi

        # --- Periodic BC on positions ---
        self.positions %= L

        self._step_count += 1
        return self.get_state()

    def run(self, n_steps: int, record_interval: int = 1) -> list[dict[str, Any]]:
        """Run simulation; record state every `record_interval` steps.

        Returns a list of state dicts, index 0 = initial state (after
        setup), subsequent entries at each record interval.
        """
        if not self._is_setup:
            history = [self.setup()]
        else:
            history = [self.get_state()]
        for t in range(1, n_steps + 1):
            self.step()
            if t % record_interval == 0:
                history.append(self.get_state())
        return history

    # ------------------------------------------------------------------
    # Observables / metadata
    # ------------------------------------------------------------------

    def get_state(self) -> dict[str, Any]:
        """Snapshot of current configuration.

        Contains both geometric observables (positions, headings) and
        derived instantaneous kinematic observables (velocities, speeds).
        The ``speeds`` field is NON-CONSTANT across particles, which is
        the empirical signature of density-dependent propulsion; this
        differs from Vicsek, where speeds == v0 for all particles at
        all times.
        """
        # Recompute instantaneous densities and speeds for the snapshot.
        # (step() already computed these but did not cache them; doing
        # it again here keeps get_state() idempotent and side-effect-free.)
        L = self.box_size
        tree = cKDTree(self.positions, boxsize=L)
        counts = tree.query_ball_point(
            self.positions, r=self.r_cg, return_length=True
        )
        counts = np.asarray(counts, dtype=np.float64)
        area = np.pi * self.r_cg**2
        rho_local = counts / area
        v_local = self.v0 * np.maximum(0.0, 1.0 - rho_local / self.rho_star)

        cos_t = np.cos(self.headings)
        sin_t = np.sin(self.headings)
        velocities = np.column_stack([v_local * cos_t, v_local * sin_t])
        return {
            "positions": self.positions.copy(),
            "velocities": velocities,
            "headings": self.headings.copy(),
            "speeds": v_local.copy(),
            "local_density": rho_local.copy(),
            "step": self._step_count,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Full metadata including mechanistic rule flags.

        The flags ``has_alignment_rule``, ``has_attraction_rule``, and
        ``has_density_dependent_speed`` encode the P2 mechanistic
        signature (Decision 43). Detectors can promote an empirical P2
        signal to DEFINITIVE only when attraction and alignment are
        both declared absent AND density-dependent speed is declared
        present.
        """
        return {
            "model_family": "active_brownian_particles",
            "model_name": "Fily-Marchetti Active Brownian Particles (2012)",
            "model_class": "self_propelled_particles",
            "n_particles": self.n_particles,
            "box_size": self.box_size,
            "v0": self.v0,
            "rho_star": self.rho_star,
            "D_r": self.D_r,
            "r_cg": self.r_cg,
            "dt": self.dt,
            "init_mode": self.init_mode,
            "seed": self.seed,
            "sigma": self.sigma,
            "packing_fraction": self.packing_fraction,
            "peclet": self.peclet,
            "space_type": "continuous_2d_periodic",
            "dynamics_type": "overdamped_langevin",
            # Mechanistic flags for P2 detector (Decision 43)
            "has_alignment_rule": False,
            "has_attraction_rule": False,
            "has_density_dependent_speed": True,
            "interaction_type": "kinetic_density_feedback",
        }

    def get_timescale(self) -> dict[str, float]:
        """Characteristic timescales for the ABP model.

        T_rot = 1 / D_r : rotational correlation time (persistence time).
        T_cross = L / v0 : ballistic crossing time at full speed.
        L_p = v0 / D_r   : persistence length (in units of sigma).
        """
        T_rot = 1.0 / self.D_r if self.D_r > 0 else float("inf")
        T_cross = self.box_size / self.v0 if self.v0 > 0 else float("inf")
        L_p = self.v0 / self.D_r if self.D_r > 0 else float("inf")
        return {
            "T_rot": T_rot,
            "T_cross": T_cross,
            "L_p": L_p,
            "dt": self.dt,
            "steps_per_T_rot": T_rot / self.dt,
        }
