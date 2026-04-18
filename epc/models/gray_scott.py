"""Gray-Scott reaction-diffusion on a 2D periodic lattice.

Primary reference:
  Pearson, J. E. (1993). "Complex Patterns in a Simple System."
  Science 261, 189-192.
  https://doi.org/10.1126/science.261.5118.189

Related references:
  Turing, A. M. (1952). "The Chemical Basis of Morphogenesis."
  Phil. Trans. R. Soc. Lond. B 237, 37-72.
  — Original paper on reaction-diffusion instabilities as a mechanism
  for biological pattern formation.

  Gray, P. & Scott, S. K. (1983). "Autocatalytic reactions in the
  isothermal, continuous stirred tank reactor: isolas and other forms
  of multistability." Chem. Eng. Sci. 38, 29-43.
  — Original introduction of the well-mixed form of the model.

The Gray-Scott reaction-diffusion system:

    ∂u/∂t = D_u ∇²u − u v² + F (1 − u)
    ∂v/∂t = D_v ∇²v + u v² − (F + k) v

where u, v are continuous concentration fields on [0, 1], D_u and D_v
are diffusion coefficients, F is the feed rate, and k is the kill rate.

Discretization: 5-point Laplacian on an integer grid with periodic
boundaries, forward-Euler integration with step dt = 1.0. Standard
grid-scale coefficients: D_u = 0.16, D_v = 0.08. With these values the
Turing window spans roughly F ∈ [0.010, 0.080], k ∈ [0.040, 0.070].

Initial condition ("Pearson seed"): u = 1 everywhere, v = 0 everywhere,
except a square patch (size ~ N/5) in the center where u = 0.50,
v = 0.25. Small iid Gaussian noise (std 0.02) breaks symmetry.

State history keys per snapshot:
    field          : np.ndarray (rows, cols), dtype float — the v-field
    field_u        : np.ndarray (rows, cols), dtype float — the u-field
    grid_dims      : tuple (rows, cols)
    step           : int — timestep
    v_mean         : float — spatial mean of v
    v_std          : float — spatial std of v
    u_mean         : float — spatial mean of u
    u_std          : float — spatial std of u
    activity       : float — v_std (alias for detector code)

Canonical regimes (empirically characterized, Sprint 13, N=128, T≥4000):
  Labyrinth        (F=0.037, k=0.060): λ ≈ 12.8 px, peak/mean ≈ 23
  Spots            (F=0.030, k=0.062): λ ≈ 11.6 px, peak/mean ≈ 24
  Long-wavelength  (F=0.062, k=0.0609): needs N≥256 to resolve
                    ("spots" in Pearson's original nomenclature; at N=128
                    gives wavelength ≈ N/2 as a domain-size artifact).
  Uniform decay    (F=0.100, k=0.100): v → 0, no pattern.

WAVELENGTH STABILIZATION:
After the Pearson IC, the system undergoes an initial transient where
the central perturbation invades the surrounding u=1 medium. The Turing
wavelength is selected in ≈ 2000–4000 timesteps at N=128 (characterization
table, Sprint 13). Detectors should either:
  (a) use the full final state after ≥ 4000 steps, OR
  (b) use time-averaged snapshots from the last quarter of the run.

P3 boundary:
Gray-Scott produces a stationary periodic spatial pattern with a
characteristic wavelength: this is the canonical positive for P3
(spatial pattern selection by reaction-diffusion instability).

Discrimination from P1 (aggregation): P1 measures spatial autocorrelation
without a wavelength signature. A Schelling-segregation field shows high
Moran's I but no periodic FFT peak; Gray-Scott shows a sharp FFT peak at
the Turing wavenumber. P1's exclusion of P3 is already coded (FFT peak
check on the final grid). Gray-Scott provides the canonical positive for
P3 and thereby activates the P1×P3 discrimination test.

Discrimination from RPS (cyclic dominance): RPS at low mobility produces
rotating spirals whose raw integer grid shows a large FFT peak too.
Discrimination is by STATE TYPE: Gray-Scott's field is a continuous
(float) variable with many thousands of distinct values per snapshot;
RPS has ≤ 4 discrete states. P3's prerequisite `n_unique_values >= 50`
cleanly rejects all discrete-state models.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_model import BaseModel


class GrayScott(BaseModel):
    """Gray-Scott reaction-diffusion on a 2D periodic lattice.

    Parameters
    ----------
    rows : int
        Grid height (default 128).
    cols : int
        Grid width (default 128).
    feed_rate : float
        F, the feed rate for u. Default 0.037 (labyrinth regime).
    kill_rate : float
        k, the kill rate for v. Default 0.060 (labyrinth regime).
    Du : float
        Diffusion coefficient for u. Default 0.16 (grid-scale).
    Dv : float
        Diffusion coefficient for v. Default 0.08 (grid-scale).
    dt : float
        Forward-Euler step size. Default 1.0 with grid-scale D's.
    init_mode : str
        'pearson_seed' (default): small central patch of u=0.5, v=0.25,
        surrounded by u=1, v=0. Tiny iid noise added.
        'custom': user supplies (u0, v0) via init_fields.
    init_patch_size : int
        For 'pearson_seed': half-width of the seed patch as a fraction of
        grid side length. Default: 10 → patch occupies rows/10 × cols/10.
    init_noise : float
        Std-dev of Gaussian noise added to the IC. Default 0.02.
    init_fields : tuple[np.ndarray, np.ndarray] or None
        For 'custom' mode: (u0, v0) user-supplied initial fields.
    seed : int
        RNG seed.
    """

    def __init__(
        self,
        rows: int = 128,
        cols: int = 128,
        feed_rate: float = 0.037,
        kill_rate: float = 0.060,
        Du: float = 0.16,
        Dv: float = 0.08,
        dt: float = 1.0,
        init_mode: str = "pearson_seed",
        init_patch_size: int = 10,
        init_noise: float = 0.02,
        init_fields: tuple[np.ndarray, np.ndarray] | None = None,
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)
        if rows < 16 or cols < 16:
            raise ValueError("rows and cols must be >= 16")
        if not 0.0 < feed_rate < 1.0:
            raise ValueError(f"feed_rate must be in (0, 1), got {feed_rate}")
        if not 0.0 < kill_rate < 1.0:
            raise ValueError(f"kill_rate must be in (0, 1), got {kill_rate}")
        if Du <= 0.0 or Dv <= 0.0:
            raise ValueError("diffusion coefficients must be positive")
        if dt <= 0.0:
            raise ValueError("dt must be positive")
        # Simple CFL-like guard for 5-point Laplacian: dt * max_D <= 0.25
        if dt * max(Du, Dv) > 0.25:
            raise ValueError(
                f"dt * max(Du, Dv) = {dt * max(Du, Dv):.3f} exceeds stability "
                f"bound 0.25 for 5-point Laplacian. Reduce dt or diffusion."
            )
        if init_mode not in ("pearson_seed", "custom"):
            raise ValueError(f"init_mode must be 'pearson_seed' or 'custom'")
        if init_patch_size < 2:
            raise ValueError("init_patch_size must be >= 2")
        if init_noise < 0:
            raise ValueError("init_noise must be >= 0")

        self.rows = rows
        self.cols = cols
        self.feed_rate = feed_rate
        self.kill_rate = kill_rate
        self.Du = Du
        self.Dv = Dv
        self.dt = dt
        self.init_mode = init_mode
        self.init_patch_size = init_patch_size
        self.init_noise = init_noise
        self.init_fields = init_fields

        self._u: np.ndarray | None = None
        self._v: np.ndarray | None = None

    def setup(self) -> dict[str, Any]:
        if self.init_mode == "pearson_seed":
            self._u, self._v = self._init_pearson_seed()
        elif self.init_mode == "custom":
            if self.init_fields is None:
                raise ValueError("init_fields required for 'custom' mode")
            u0, v0 = self.init_fields
            if u0.shape != (self.rows, self.cols) or v0.shape != (self.rows, self.cols):
                raise ValueError(
                    f"init_fields shapes must be ({self.rows}, {self.cols})"
                )
            self._u = u0.astype(float).copy()
            self._v = v0.astype(float).copy()
        else:
            raise ValueError(f"unknown init_mode: {self.init_mode}")

        self._is_setup = True
        self._step_count = 0
        snap = self._snapshot()
        self._state_history = [snap]
        return snap

    def step(self) -> dict[str, Any]:
        if not self._is_setup:
            self.setup()

        u, v = self._u, self._v
        lu = self._laplacian5(u)
        lv = self._laplacian5(v)
        uvv = u * v * v

        u_new = u + self.dt * (self.Du * lu - uvv + self.feed_rate * (1.0 - u))
        v_new = v + self.dt * (self.Dv * lv + uvv - (self.feed_rate + self.kill_rate) * v)

        # Stability clamp. Concentrations are not supposed to leave [0, 1] but
        # finite-precision drift is occasionally possible; clamp avoids NaN
        # runaway if user supplies a very unstable (F, k) pair.
        np.clip(u_new, 0.0, 1.0, out=u_new)
        np.clip(v_new, 0.0, 1.0, out=v_new)

        self._u = u_new
        self._v = v_new
        self._step_count += 1
        return self._snapshot()

    def get_metadata(self) -> dict[str, Any]:
        return {
            "model_name": "gray_scott",
            "model_class": "reaction_diffusion",  # deliberately not "ca" / "excitable"
            "rows": self.rows,
            "cols": self.cols,
            "n_states": None,  # continuous field — n_states not meaningful
            "state_dtype": "float",
            "feed_rate": self.feed_rate,
            "kill_rate": self.kill_rate,
            "Du": self.Du,
            "Dv": self.Dv,
            "dt": self.dt,
            "init_mode": self.init_mode,
            "init_patch_size": self.init_patch_size,
            "init_noise": self.init_noise,
            "boundary": "periodic",
            "interaction_type": "local_diffusion_reaction",
            "update_mode": "synchronous",
            "seed": self.seed,
            "reference": (
                "Pearson (1993), Science 261, 189-192; Gray & Scott (1983); "
                "Turing (1952)."
            ),
        }

    def get_timescale(self) -> float:
        """System-intrinsic timescale.

        For Turing-wavelength selection, the relevant timescale is the
        diffusive crossing time across one wavelength:
          T_wave ≈ λ² / (4 D_v)
        Empirical Sprint 13 characterization: λ ≈ 12 px at (F=0.037,
        k=0.060), so T_wave ≈ 144 / 0.32 ≈ 450 steps. The full pattern
        takes ~ 10 × T_wave to stabilize from the Pearson seed.
        """
        # Conservative estimate: characteristic wavelength ~ sqrt(Dv / (F+k))
        # times a prefactor of order 10. For grid-scale (Dv=0.08, F+k ~ 0.1)
        # this gives ~ 10 * sqrt(0.08/0.1) ~ 9 pixels, crossing time
        # ~ 9^2 / (4*0.08) ~ 250 steps.
        if self._v is None:
            # Estimate from params
            wavelength_est = 10.0 * np.sqrt(
                self.Dv / (self.feed_rate + self.kill_rate + 1e-12)
            )
        else:
            wavelength_est = 12.0  # empirical default for labyrinth regime
        t_wave = (wavelength_est ** 2) / (4.0 * self.Dv)
        return float(t_wave)

    def is_converged(self) -> bool:
        """Converged when field change per step is negligible.

        Used only for early-termination guidance; not called by run()
        unless explicitly checked by the caller.
        """
        if self._v is None or len(self._state_history) < 2:
            return False
        prev_v = self._state_history[-2].get("field")
        if prev_v is None:
            return False
        # L2 difference per cell
        diff = float(np.mean((self._v - prev_v) ** 2))
        return diff < 1e-10

    # --- Initialization ---

    def _init_pearson_seed(self) -> tuple[np.ndarray, np.ndarray]:
        u = np.ones((self.rows, self.cols), dtype=float)
        v = np.zeros((self.rows, self.cols), dtype=float)
        r = max(2, min(self.rows, self.cols) // self.init_patch_size)
        cr, cc = self.rows // 2, self.cols // 2
        r0, r1 = cr - r, cr + r
        c0, c1 = cc - r, cc + r
        u[r0:r1, c0:c1] = 0.50
        v[r0:r1, c0:c1] = 0.25
        if self.init_noise > 0:
            u += self.init_noise * self.rng.standard_normal(u.shape)
            v += self.init_noise * self.rng.standard_normal(v.shape)
            np.clip(u, 0.0, 1.0, out=u)
            np.clip(v, 0.0, 1.0, out=v)
        return u, v

    # --- Laplacian (5-point stencil with periodic BC) ---

    @staticmethod
    def _laplacian5(arr: np.ndarray) -> np.ndarray:
        return (
            np.roll(arr, 1, axis=0) + np.roll(arr, -1, axis=0)
            + np.roll(arr, 1, axis=1) + np.roll(arr, -1, axis=1)
            - 4.0 * arr
        )

    # --- Snapshot ---

    def _snapshot(self) -> dict[str, Any]:
        v = self._v.copy()
        u = self._u.copy()
        return {
            "field": v,                     # canonical continuous-field key
            "field_u": u,                   # optional second field
            "grid_dims": (self.rows, self.cols),
            "step": self._step_count,
            "v_mean": float(v.mean()),
            "v_std": float(v.std()),
            "u_mean": float(u.mean()),
            "u_std": float(u.std()),
            "activity": float(v.std()),     # alias convenience
        }
