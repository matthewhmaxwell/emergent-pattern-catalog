"""Non-local Kuramoto ring — chimera states on an oscillator ring.

References
----------
Abrams, D. M. & Strogatz, S. H. (2004). Chimera states for coupled
    oscillators. Physical Review Letters 93, 174102.
    doi:10.1103/PhysRevLett.93.174102

Kuramoto, Y. & Battogtokh, D. (2002). Coexistence of coherence and
    incoherence in non-locally coupled phase oscillators. Nonlinear
    Phenomena in Complex Systems 5, 380–385. arXiv:cond-mat/0210694

Dynamics
--------

N identical phase oscillators (ω_i ≡ 0) sit at positions
x_i = 2π·i/N on a ring (periodic 1D domain). They evolve as

    dθ_i/dt = -(1/N) Σ_j G(x_i - x_j) sin(θ_i - θ_j + α)

with coupling kernel (Abrams-Strogatz, "cosine")

    G(x) = (1/(2π)) (1 + A cos(x)),   A ∈ [0, 1]

or (Kuramoto-Battogtokh, "step")

    G(x) = (1/(4π·r)) · 1[|x| < 2π·r],    r ∈ (0, 1/2]

and phase lag α = π/2 - β with β small and positive.

The reversible limit β → 0 has a manifold of fixed points; chimera states
emerge as a stable bistable branch for (A, β) in a ribbon near the
reversible manifold. For A = 0.995 and β in roughly [0.02, 0.20], the
chimera basin coexists with the full-sync basin — which initial condition
you pick determines which attractor the system reaches.

CANONICAL CHIMERA POSITIVE (Sprint 18)

  A = 0.995, β = 0.05, N = 128, T ≥ 50, dt = 0.025, record_dt = 1.0,
  init_mode = "asymmetric_gaussian" (amp = 2.0, sigma = 0.5),
  seed = 0 or 42.

With this IC the chimera basin is reliably reached. The paper's β = 0.18
is narrower — only roughly 3/6 random seeds land in the chimera basin
there. β = 0.05 is the empirically robust choice for testing.

See ADR 51 in REPLICATION_NOTES.md.

IMPLEMENTATION

Unlike the existing ``KuramotoModel``, we CANNOT use the mean-field O(N)
reformulation. Non-local coupling requires the explicit O(N²) kernel
sum. This model runs in O(N² · n_steps) — about 20 s for N = 128,
T = 50 on CPU. Larger N benefits from vectorised broadcast.

The existing all-to-all ``KuramotoModel`` is unchanged by this file; the
only cross-file change is adding ``has_nonlocal_coupling = False`` to
``KuramotoModel.get_metadata()`` for symmetry with this file's
``has_nonlocal_coupling = True``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
from numpy.typing import NDArray


@dataclass
class KuramotoNonlocalParams:
    """Parameters for the non-local Kuramoto ring."""

    N: int = 128                       # Number of oscillators on the ring
    A: float = 0.995                   # Cosine-kernel amplitude (in [0, 1])
    beta: float = 0.05                 # Phase lag: α = π/2 - β
    dt: float = 0.025                  # Integration timestep
    record_dt: float = 1.0             # Wall-time between recorded frames
    kernel: str = "cosine"             # "cosine" (Abrams-Strogatz) or "step"
    step_kernel_radius: float = 0.2    # For step kernel: r ∈ (0, 0.5]
    seed: int = 42                     # RNG seed
    init_mode: str = "asymmetric_gaussian"  # see _init_theta
    asym_amp: float = 2.0              # amplitude for asymmetric_gaussian IC
    asym_sigma: float = 0.5            # width for asymmetric_gaussian IC


class KuramotoNonlocal:
    """Non-local Kuramoto ring supporting chimera states.

    Follows the standalone-model convention used by Vicsek, ABP, Yard-Sale:
    not a BaseModel subclass; provides ``run()``, ``get_metadata()``, and
    ``get_timescale()`` directly.

    Parameters
    ----------
    params : KuramotoNonlocalParams, optional
        If omitted, defaults to the canonical chimera positive
        (N=128, A=0.995, β=0.05, asymmetric_gaussian IC, seed=42).

    Alternatively, pass individual kwargs; they override the defaults.

    Attributes
    ----------
    theta : (N,) current phases, mod 2π.
    x : (N,) ring positions in [0, 2π).
    G : (N, N) precomputed coupling kernel.
    alpha : float, phase lag = π/2 - β.

    Examples
    --------
    >>> m = KuramotoNonlocal()                   # canonical chimera positive
    >>> history = m.run(n_frames=50)             # 50 recorded frames
    >>> len(history), history[0]["theta"].shape
    (50, (128,))
    """

    def __init__(
        self,
        params: Optional[KuramotoNonlocalParams] = None,
        *,
        N: Optional[int] = None,
        A: Optional[float] = None,
        beta: Optional[float] = None,
        dt: Optional[float] = None,
        record_dt: Optional[float] = None,
        kernel: Optional[str] = None,
        step_kernel_radius: Optional[float] = None,
        seed: Optional[int] = None,
        init_mode: Optional[str] = None,
        asym_amp: Optional[float] = None,
        asym_sigma: Optional[float] = None,
    ) -> None:
        if params is None:
            params = KuramotoNonlocalParams()
        # kwarg overrides
        if N is not None: params.N = int(N)
        if A is not None: params.A = float(A)
        if beta is not None: params.beta = float(beta)
        if dt is not None: params.dt = float(dt)
        if record_dt is not None: params.record_dt = float(record_dt)
        if kernel is not None: params.kernel = str(kernel)
        if step_kernel_radius is not None:
            params.step_kernel_radius = float(step_kernel_radius)
        if seed is not None: params.seed = int(seed)
        if init_mode is not None: params.init_mode = str(init_mode)
        if asym_amp is not None: params.asym_amp = float(asym_amp)
        if asym_sigma is not None: params.asym_sigma = float(asym_sigma)

        if params.N < 16:
            raise ValueError(f"N must be >= 16 for windowing, got {params.N}")
        if not (0.0 <= params.A <= 1.0):
            raise ValueError(f"A must be in [0, 1], got {params.A}")
        if params.dt <= 0:
            raise ValueError(f"dt must be > 0, got {params.dt}")
        if params.record_dt <= 0:
            raise ValueError(f"record_dt must be > 0, got {params.record_dt}")
        if params.kernel not in {"cosine", "step"}:
            raise ValueError(
                f"kernel must be 'cosine' or 'step', got {params.kernel}"
            )

        self.params = params
        self.rng = np.random.default_rng(params.seed)
        self.alpha: float = float(np.pi / 2.0 - params.beta)

        # Ring positions
        self.x: NDArray[np.float64] = (
            2 * np.pi * np.arange(params.N) / params.N
        )

        # Signed pairwise ring separation, wrapped to [-π, π]
        dx = self.x[:, None] - self.x[None, :]
        dx = np.mod(dx + np.pi, 2 * np.pi) - np.pi

        # Coupling kernel G[i, j]
        if params.kernel == "cosine":
            self.G: NDArray[np.float64] = (
                (1.0 / (2 * np.pi)) * (1.0 + params.A * np.cos(dx))
            )
        else:  # "step"
            r = float(params.step_kernel_radius)
            if not (0.0 < r <= 0.5):
                raise ValueError(
                    f"step_kernel_radius must be in (0, 0.5], got {r}"
                )
            # Uniform inside |dx| < 2π r, zero outside, normalized so row
            # sum ~ 1 (matches the cosine kernel's integrated weight).
            G = np.where(np.abs(dx) < 2 * np.pi * r, 1.0, 0.0)
            row_sum = G.sum(axis=1, keepdims=True)
            row_sum = np.where(row_sum > 0, row_sum, 1.0)
            self.G = G / row_sum * (params.N / (2 * np.pi))
            # This normalization keeps the deriv scale consistent with the
            # cosine kernel case (both evaluate to (1/N) · ~1 · sin(...)).

        # Initial phases
        self.theta: NDArray[np.float64] = self._init_theta()

    # ------------------------------------------------------------------
    # Initial conditions
    # ------------------------------------------------------------------

    def _init_theta(self) -> NDArray[np.float64]:
        p = self.params
        if p.init_mode == "asymmetric_gaussian":
            # Centered at x = π (opposite side of ring from x = 0).
            envelope = p.asym_amp * np.exp(
                -((self.x / np.pi - 1.0) ** 2) / (2 * p.asym_sigma ** 2)
            )
            return np.mod(envelope * self.rng.standard_normal(p.N), 2 * np.pi)
        if p.init_mode == "uniform":
            return self.rng.uniform(0, 2 * np.pi, p.N)
        if p.init_mode == "coherent":
            # Small noise around zero — falls into full-sync basin typically
            return np.mod(0.01 * self.rng.standard_normal(p.N), 2 * np.pi)
        if p.init_mode == "paper":
            # Abrams-Strogatz 2004 Fig. 2 IC: localized perturbation with
            # much narrower envelope. Less robust to chimera than the
            # asymmetric_gaussian default; retained for replication.
            envelope = 6.0 * np.exp(-30 * ((self.x / np.pi - 1.0) ** 2))
            return np.mod(envelope * self.rng.standard_normal(p.N) / 1.5,
                          2 * np.pi)
        raise ValueError(f"Unknown init_mode: {p.init_mode}")

    # ------------------------------------------------------------------
    # RK4 integration
    # ------------------------------------------------------------------

    def _derivatives(self, theta: NDArray[np.float64]) -> NDArray[np.float64]:
        """Time derivative of phases.

        Implements  dθ_i/dt = -Σ_j G(x_i - x_j) sin(θ_i - θ_j + α)
        where G already carries the ``1/(2π)`` normalization (for the
        cosine kernel) or the row-normalized weight (for the step kernel).
        This matches the Abrams-Strogatz (2004) discretization convention
        in which the Riemann measure ``dy = 2π/N`` is absorbed into an
        O(1) effective timescale — i.e., integration times in this module
        are not expressed in the PDE's natural time units. The canonical
        chimera at A=0.995, β=0.05, N=128 emerges at T ≈ 50 in these
        units (see Phase 1 calibration in REPLICATION_NOTES Sprint 18).
        """
        diff = theta[:, None] - theta[None, :] + self.alpha
        return -np.sum(self.G * np.sin(diff), axis=1)

    def _rk4_step(
        self, theta: NDArray[np.float64], dt: float
    ) -> NDArray[np.float64]:
        k1 = self._derivatives(theta)
        k2 = self._derivatives(theta + 0.5 * dt * k1)
        k3 = self._derivatives(theta + 0.5 * dt * k2)
        k4 = self._derivatives(theta + dt * k3)
        return theta + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

    # ------------------------------------------------------------------
    # Public run / observables / metadata
    # ------------------------------------------------------------------

    def run(
        self,
        n_frames: int = 50,
        equilibration: float = 0.0,
    ) -> list[dict[str, Any]]:
        """Integrate the ring and record `n_frames` snapshots.

        Parameters
        ----------
        n_frames : int
            Number of snapshots to record. Total integrated time is
            ``n_frames * record_dt + equilibration``.
        equilibration : float
            Wall-time to integrate before recording starts. Default 0.

        Returns
        -------
        list[dict] of state snapshots with keys
            'theta'      : (N,) phases mod 2π
            'r'          : global order parameter magnitude
            'psi'        : global mean phase
            'positions'  : (N,) ring positions in [0, 2π) (constant)
            'step'       : int, cumulative integration step count
            't'          : float, integrated time (post-equilibration)
        """
        p = self.params
        dt = p.dt
        theta = self.theta.copy()
        steps_per_frame = max(1, int(round(p.record_dt / dt)))

        # Equilibration
        n_eq = max(0, int(round(equilibration / dt)))
        for _ in range(n_eq):
            theta = self._rk4_step(theta, dt)

        history: list[dict[str, Any]] = []
        step_count = n_eq
        for frame in range(n_frames):
            for _ in range(steps_per_frame):
                theta = self._rk4_step(theta, dt)
                step_count += 1
            theta_mod = np.mod(theta, 2 * np.pi)
            z = np.mean(np.exp(1j * theta_mod))
            history.append({
                "theta": theta_mod.copy(),
                "r": float(np.abs(z)),
                "psi": float(np.angle(z)),
                "positions": self.x.copy(),
                "step": int(step_count),
                "t": float(step_count * dt),
            })

        self.theta = theta
        return history

    def get_metadata(self) -> dict[str, Any]:
        """Metadata for detector dispatch and mechanistic gating."""
        p = self.params
        return {
            "model_family": "kuramoto_nonlocal",
            "model_name": "Abrams-Strogatz non-local Kuramoto ring (2004)",
            "model_class": "oscillator_ring",
            "N": p.N,
            "A": p.A,
            "beta": p.beta,
            "alpha": self.alpha,
            "dt": p.dt,
            "record_dt": p.record_dt,
            "kernel": p.kernel,
            "step_kernel_radius": p.step_kernel_radius,
            "seed": p.seed,
            "init_mode": p.init_mode,
            # Substrate / dynamics classification
            "substrate_type": "oscillator",
            "space_type": "ring_1d",
            "dynamics_type": "deterministic_phase_ode",
            # Mechanistic flags for P10 detector
            "has_nonlocal_coupling": True,
            "has_frequency_heterogeneity": False,  # identical ω = 0
            "coupling_kernel": p.kernel,
            "freq_dist": "identical",
        }

    def get_timescale(self) -> dict[str, float]:
        """Characteristic timescales of the ring.

        For identical oscillators (ω = 0), there is no intrinsic
        oscillation period. The characteristic timescale is set by the
        inverse of the largest eigenvalue of the linearized system around
        the synchronous manifold, which is O(1) in our units. We use
        ``T_relax = 1 / max(A, β)`` as a heuristic.
        """
        p = self.params
        T_relax = 1.0 / max(p.A, p.beta, 0.05)
        return {
            "T_relax": float(T_relax),
            "T_record": float(p.record_dt),
        }
