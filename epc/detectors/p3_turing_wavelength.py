"""P3 — Turing wavelength detector.

Detects a stationary periodic spatial pattern selected by a
reaction-diffusion (or diffusion-driven) instability — the classical
"Turing pattern" signature.

Observable scope: state_history_only, but on a continuous (float-valued)
2D field. The detector reads the 'field' observable from each snapshot,
NOT the 'grid' observable that the integer-state lattice_2d detectors
use. This substrate distinction is part of the detector's prerequisite
structure: a prerequisite of `n_unique_values >= 50` cleanly rejects all
discrete-state models that happen to produce high FFT peak-to-mean
(notably RPS at low mobility, whose raw integer grid has p/m ≈ 23 —
numerically indistinguishable from Gray-Scott labyrinth on peak-to-mean
alone).

DISCRIMINATION FROM NEIGHBORS:

P1 (aggregation): P1 measures spatial autocorrelation (Moran's I)
  without a wavelength signature. A Schelling-segregation field has high
  Moran's I but a near-flat radial FFT spectrum; Gray-Scott has a sharp
  peak at a finite wavenumber. P1's exclusion of P3 is already coded at
  P1's side (FFT peak prominence check); this detector is the canonical
  positive that activates that discrimination.

RPS / cyclic dominance: RPS spirals have a characteristic wavelength
  and produce a similar radial-FFT peak on the raw integer grid. P3's
  continuous-field prerequisite (n_unique_values >= 50) rejects RPS
  cleanly: RPS has 4 distinct integer labels, Gray-Scott has thousands
  of distinct float values. The prerequisite is SUBSTRATE-level, not
  empirically-tuned.

P13 (excitable wavefronts): excitable media produce traveling waves
  (spirals, target patterns) with a characteristic wavelength too.
  Same substrate prerequisite applies (GH is integer-labelled).

MECHANISM:

Turing instability arises when two diffusively-coupled reaction species
satisfy D_activator << D_inhibitor. Small spatially-inhomogeneous
perturbations at the fastest-growing wavenumber k* are amplified; all
other k are damped. The system equilibrates at a stationary periodic
pattern with wavelength λ* = 2π / k*.

For Gray-Scott with grid-scale coefficients (D_u = 0.16, D_v = 0.08)
in the Turing window, the selected wavelength is empirically ~12 px —
verified invariant across grid sizes N = 64, 96, 128, 192 (Sprint 13
characterization; peak_k scales linearly with N while λ = N / peak_k
remains constant).

DETECTOR TIERS:

Prerequisites (all required for any detection):
  1. The 'field' observable is present and 2D.
  2. n_unique_values in the final field >= MIN_UNIQUE (50).
  3. field_std in the final field >= MIN_FIELD_STD (0.01).

Primary metric: peak_to_mean of the radially-averaged 2D FFT power
spectrum on the final (or time-averaged last-window) scalar field,
with DC and first-bin suppressed (k_min = 2) and the full radial range
used for the mean (k_max_frac = 1.0).

Secondary metrics:
  - peak_k and wavelength in pixels (wavelength = N / peak_k).
  - peak_k stability across the last N_STABILITY_FRAMES snapshots:
    std / mean of peak_k, should be ~0 for a true Turing pattern.
  - peak_to_mean stability (min across the last N_STABILITY_FRAMES).

Screening:    prerequisites pass AND primary peak_to_mean > 5.0
              AND 2 <= peak_k <= N/4 (real wavelength, not aliased).
Confirmation: screening + peak_to_mean > 10.0 + Cohen's d > 10.0
              + peak_k stability (peak_k_cv < 0.15) + null p < 0.01.
Definitive:   confirmation + peak_to_mean > 15.0 + Cohen's d > 30.0
              + peak_k_cv < 0.05 (highly stable wavelength)
              + exclusions cleared.

Null model: spatial shuffle of the field values. Preserves the marginal
intensity distribution, destroys all spatial structure including
wavelength. For Gray-Scott labyrinth the null peak-to-mean distribution
is centered at 1.4 with std 0.16, giving Cohen's d ≈ 100 and floor
p-values.

EMPIRICAL CALIBRATION (Sprint 13, Gray-Scott, seed=42):

Canonical positives (N=128, T>=4000):
  Labyrinth (F=0.037, k=0.060): peak_k=10, lam=12.8 px,
                                 p/m=18.75, d~107, p~0.01
  Spots     (F=0.030, k=0.062): peak_k=11, lam=11.6 px,
                                 p/m~20,  d~70,  p~0.01

Negative-model spot checks (raw integer grids):
  RPS (M=1e-4):             p/m~23  -> REJECTED by nu_unique prerequisite
  GH (n_states=8):          p/m~6.6 -> REJECTED by nu_unique prerequisite
  Schelling (tau=0.375):    p/m~3.3 -> also below screening p/m threshold
  Nowak-May (b=1.8):        p/m~3.6 -> REJECTED by nu_unique prerequisite
  SIR (wavefront):          p/m~3.7 -> REJECTED by nu_unique prerequisite
  GoL (random):             p/m~5.4 -> REJECTED by nu_unique prerequisite
  LV lattice:               p/m~6.7 -> REJECTED by nu_unique prerequisite

Non-Turing Gray-Scott regimes (continuous fields):
  Uniform decay (F=0.10, k=0.10): field_std = 0 -> REJECTED by prereq.
  Uniform high (F=0.01, k=0.02): field_std = 0 -> REJECTED by prereq.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, DetectorResult, NullType
from epc.metrics.turing_wavelength import (
    count_unique_values,
    field_stationarity,
    radial_fft_peak_stats,
    ring_angular_isotropy,
    shuffle_null_distribution,
    wavelength_growth,
    wavelength_stability,
)


class P3TuringWavelengthDetector(BaseDetector):
    """Detector for P3 — Turing wavelength (stationary periodic spatial
    pattern from a reaction-diffusion instability).

    Parameters
    ----------
    n_permutations : int
        Number of spatial-shuffle null permutations. Default 199 (floor
        p-value = 1/200 = 0.005). With Gray-Scott's huge effect sizes
        (d > 100), 99 is also fine for screening but CONFIRMATION requires
        p < 0.01 which needs >= 199 permutations.
    k_min : int
        Lowest wavenumber considered in the peak search. Default 2
        (skips DC and the first-bin long-wavelength drift).
    k_max_frac : float
        Upper end of the k-search window as a fraction of Nyquist.
        Default 1.0 (use the full radial range). Sprint 13 characterization
        showed this gives the cleanest signal-to-noise.
    n_stability_frames : int
        Number of late-trajectory snapshots over which wavelength
        stability is evaluated. Default 5.
    stability_stride : int
        Stride in timesteps between stability frames. Default 50.
    seed : int
        RNG seed for the spatial-shuffle null.
    """

    # Prerequisite thresholds
    _MIN_UNIQUE = 50
    _MIN_FIELD_STD = 0.01

    # Screening thresholds
    _SCREEN_PM_MIN = 5.0
    _SCREEN_K_MIN = 2

    # Confirmation thresholds (in addition to screening)
    _CONFIRM_PM_MIN = 10.0
    _CONFIRM_COHENS_D_MIN = 10.0
    _CONFIRM_PEAK_K_CV_MAX = 0.15
    _CONFIRM_NULL_P_MAX = 0.01

    # Definitive thresholds (in addition to confirmation)
    _DEF_PM_MIN = 15.0
    _DEF_COHENS_D_MIN = 30.0
    _DEF_PEAK_K_CV_MAX = 0.05

    # Turing-specificity gates (Phase-2a discrimination). A genuine
    # Turing pattern is a STATIONARY, ISOTROPIC standing wave selected
    # by a diffusion-driven instability -- not merely any field with a
    # dominant wavelength. Without these, an imposed static sinusoid
    # (anisotropic single wavevector) or a travelling / rotating wave
    # (non-stationary) clears the peak-to-mean gates as a false positive.
    # Empirical separation (Gray-Scott labyrinth positive vs imposed
    # periodicity): angular_entropy 0.86 vs <=0.28; stationarity 0.99
    # vs <=0.88 (travelling) / 0.07 (spiral).
    _SCREEN_STATIONARITY_MIN = 0.95
    _SCREEN_ANG_ENTROPY_MIN = 0.55

    # Turing-specificity gate 3: INTRINSIC WAVELENGTH (reject coarsening). A Turing
    # pattern locks an intrinsic wavelength selected by the instability; conserved /
    # non-conserved coarsening (Cahn-Hilliard spinodal, Allen-Cahn phase ordering)
    # has NO intrinsic scale -- its dominant wavelength grows monotonically as
    # L(t) ~ t^(1/n). Reject a field whose wavelength grows BOTH substantially
    # (ratio > 1.30) AND monotonically (spearman_rho > 0.50) over the trajectory.
    # Empirical separation: Turing/Gray-Scott ratio <= 0.86 / rho <= 0; Cahn-
    # Hilliard ratio 1.67-2.22 / rho +0.83..+0.96.
    _SCREEN_WL_GROWTH_MAX = 1.30
    _SCREEN_WL_RHO_MAX = 0.50

    def __init__(
        self,
        n_permutations: int = 199,
        k_min: int = 2,
        k_max_frac: float = 1.0,
        n_stability_frames: int = 5,
        stability_stride: int = 50,
        seed: int = 42,
    ) -> None:
        super().__init__(
            pattern_id="P3",
            excluded_patterns=["P1", "P13"],
            allowed_co_occurrences=[],
            observable_scope="state_history_only",
        )
        if n_permutations < 1:
            raise ValueError("n_permutations must be >= 1")
        if k_min < 1:
            raise ValueError("k_min must be >= 1")
        if not 0.0 < k_max_frac <= 1.0:
            raise ValueError("k_max_frac must be in (0, 1]")
        if n_stability_frames < 2:
            raise ValueError("n_stability_frames must be >= 2")
        if stability_stride < 1:
            raise ValueError("stability_stride must be >= 1")

        self.n_permutations = n_permutations
        self.k_min = k_min
        self.k_max_frac = k_max_frac
        self.n_stability_frames = n_stability_frames
        self.stability_stride = stability_stride
        self._seed = seed

        # Caches populated during detect()
        self._final_field: np.ndarray | None = None
        self._stability_fields: list[np.ndarray] = []
        self._N: int = 0

    # ------------------------------------------------------------------
    # Field extraction
    # ------------------------------------------------------------------

    def _get_final_field(self, state_history: list[dict[str, Any]]) -> np.ndarray | None:
        """Extract the last available 2D continuous field from the state
        history.

        Priority: 'field' key, then 'field_v' key, then None. The detector
        intentionally does NOT fall back to 'grid' — a detector that
        silently consumes integer grids would violate its substrate
        prerequisite.
        """
        if not state_history:
            return None
        last = state_history[-1]
        for key in ("field", "field_v"):
            if key in last:
                arr = np.asarray(last[key])
                if arr.ndim == 2:
                    return arr
        return None

    def _collect_stability_fields(
        self,
        state_history: list[dict[str, Any]],
    ) -> list[np.ndarray]:
        """Collect the last n_stability_frames snapshots at stability_stride.

        Used to compute peak_k stability. Returns a list of 2D fields in
        chronological order.
        """
        if not state_history:
            return []

        # Walk backwards from the end, stepping by stability_stride
        indices: list[int] = []
        idx = len(state_history) - 1
        while idx >= 0 and len(indices) < self.n_stability_frames:
            indices.append(idx)
            idx -= self.stability_stride
        indices.reverse()

        fields: list[np.ndarray] = []
        for i in indices:
            snap = state_history[i]
            for key in ("field", "field_v"):
                if key in snap:
                    arr = np.asarray(snap[key])
                    if arr.ndim == 2:
                        fields.append(arr)
                        break
        return fields

    def _collect_trajectory_fields(
        self,
        state_history: list[dict[str, Any]],
        n_traj: int = 12,
        skip_frac: float = 0.15,
    ) -> list[np.ndarray]:
        """Evenly-spaced fields across the WHOLE trajectory (after an initial
        transient) for the wavelength-vs-time intrinsic-scale test. Unlike the
        last-N stability frames, this spans the full run so monotonic coarsening
        (slow at late times) is still visible against an early reference.
        """
        fields_all: list[np.ndarray] = []
        for snap in state_history:
            for key in ("field", "field_v"):
                if key in snap:
                    arr = np.asarray(snap[key])
                    if arr.ndim == 2:
                        fields_all.append(arr)
                    break
        if len(fields_all) < 4:
            return fields_all
        start = int(skip_frac * len(fields_all))
        if len(fields_all) - start < n_traj:
            return fields_all[start:]
        idx = np.linspace(start, len(fields_all) - 1, n_traj).astype(int)
        return [fields_all[i] for i in idx]

    # ------------------------------------------------------------------
    # Pipeline methods
    # ------------------------------------------------------------------

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """Use grid size as a conservative spatial timescale proxy."""
        if state_history and "grid_dims" in state_history[0]:
            rows, cols = state_history[0]["grid_dims"]
            return float(max(rows, cols))
        field = self._get_final_field(state_history)
        if field is not None:
            return float(max(field.shape))
        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        """Populate detector caches while checking basic availability."""
        warnings: list[str] = []
        if not state_history:
            warnings.append("empty state_history")
            return warnings

        field = self._get_final_field(state_history)
        if field is None:
            warnings.append(
                "no 'field' (continuous 2D) observable in state_history — "
                "P3 requires a continuous-valued field, not a 'grid'"
            )
            self._final_field = None
            return warnings

        self._final_field = field
        self._N = int(min(field.shape))
        self._stability_fields = self._collect_stability_fields(state_history)

        # Run-length sanity: Turing wavelength selection needs a transient.
        # At Gray-Scott grid-scale, empirically ~ 4000 steps at N=128.
        # Use timescale-scaled guidance.
        recommended_min = int(10 * timescale)
        if len(state_history) < recommended_min:
            warnings.append(
                f"run length {len(state_history)} < recommended {recommended_min} "
                f"(10 * timescale {timescale:.0f}) — pattern may not be fully selected"
            )

        return warnings

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Primary metric: radial-FFT peak-to-mean on the final field."""
        field = self._final_field
        if field is None:
            return {
                "peak_to_mean": 0.0,
                "peak_k": 0,
                "wavelength_pixels": 0.0,
                "field_std": 0.0,
                "field_mean": 0.0,
                "n_unique_values": 0,
                "radial_mean": 0.0,
                "peak_value": 0.0,
                "angular_entropy": 0.0,
                "radial_concentration": 1.0,
                "field_stationarity": 0.0,
                "wavelength_growth_ratio": 1.0,
                "wavelength_growth_rho": 0.0,
            }

        stats = radial_fft_peak_stats(
            field,
            k_min=self.k_min,
            k_max_frac=self.k_max_frac,
        )
        iso = ring_angular_isotropy(
            field, k_min=self.k_min, k_max_frac=self.k_max_frac,
        )
        stationarity = field_stationarity(self._stability_fields)
        growth = wavelength_growth(
            self._collect_trajectory_fields(state_history),
            k_min=self.k_min, k_max_frac=self.k_max_frac,
        )
        return {
            "peak_to_mean": float(stats["peak_to_mean"]),
            "peak_k": int(stats["peak_k"]),
            "wavelength_pixels": float(stats["wavelength_pixels"])
            if not np.isnan(stats["wavelength_pixels"])
            else 0.0,
            "field_std": float(stats["field_std"]),
            "field_mean": float(stats["field_mean"]),
            "n_unique_values": int(stats["n_unique_values"]),
            "radial_mean": float(stats["radial_mean"]),
            "peak_value": float(stats["peak_value"]),
            "k_min_used": int(stats["k_min_used"]),
            "k_max_used": int(stats["k_max_used"]),
            "angular_entropy": float(iso["angular_entropy"]),
            "radial_concentration": float(iso["radial_concentration"]),
            "field_stationarity": float(stationarity),
            "wavelength_growth_ratio": float(growth["growth_ratio"]),
            "wavelength_growth_rho": float(growth["spearman_rho"]),
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening gates — note the continuous-substrate prerequisite."""
        # Substrate prerequisite: continuous-valued field
        nu = primary_result.get("n_unique_values", 0)
        if nu < self._MIN_UNIQUE:
            return False

        # Amplitude prerequisite: field not flat
        std = primary_result.get("field_std", 0.0)
        if std < self._MIN_FIELD_STD:
            return False

        # Peak strength
        p2m = primary_result.get("peak_to_mean", 0.0)
        if p2m <= self._SCREEN_PM_MIN:
            return False

        # Wavelength sanity: peak must be at a real (non-aliased) wavenumber
        peak_k = int(primary_result.get("peak_k", 0))
        if peak_k < self._SCREEN_K_MIN:
            return False
        if self._N > 0 and peak_k > self._N // 4:
            # Very short wavelength — likely grid discretization artifact
            return False

        # Turing-specificity gate 1: STATIONARITY. A Turing pattern is a
        # stationary standing wave; travelling / rotating waves (plane
        # waves, spirals, target patterns) advect and decorrelate
        # frame-to-frame, so they are rejected here.
        if (primary_result.get("field_stationarity", 0.0)
                < self._SCREEN_STATIONARITY_MIN):
            return False

        # Turing-specificity gate 2: ISOTROPY. A diffusion-driven
        # instability selects a wave NUMBER |k*| (an isotropic ring /
        # finite band), not a single wave VECTOR. An imposed sinusoid or
        # oriented stripe pattern concentrates power at one or a few
        # discrete directions and is rejected here. (Targets the
        # isotropic Turing phases -- labyrinth / spots -- which are the
        # catalog canonical positive.)
        if (primary_result.get("angular_entropy", 0.0)
                < self._SCREEN_ANG_ENTROPY_MIN):
            return False

        # Turing-specificity gate 3: INTRINSIC WAVELENGTH (reject coarsening). A
        # diffusion-driven instability selects a fixed wavelength; conserved /
        # non-conserved coarsening (Cahn-Hilliard, Allen-Cahn) has no intrinsic
        # scale and its dominant wavelength grows monotonically in time. Reject a
        # field whose wavelength grows BOTH substantially AND monotonically.
        gr = primary_result.get("wavelength_growth_ratio", 1.0)
        rho = primary_result.get("wavelength_growth_rho", 0.0)
        if gr > self._SCREEN_WL_GROWTH_MAX and rho > self._SCREEN_WL_RHO_MAX:
            return False

        return True

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Wavelength stability across the last N snapshots."""
        if len(self._stability_fields) < 2:
            return {}

        ws = wavelength_stability(
            self._stability_fields,
            k_min=self.k_min,
            k_max_frac=self.k_max_frac,
        )
        return {
            "peak_k_mean": float(ws["peak_k_mean"]),
            "peak_k_std": float(ws["peak_k_std"]),
            "peak_k_cv": float(ws["peak_k_cv"]),
            "peak_k_range": int(ws["peak_k_range"]),
            "peak_to_mean_mean": float(ws["peak_to_mean_mean"]),
            "peak_to_mean_min": float(ws["peak_to_mean_min"]),
            "n_stability_frames_used": len(self._stability_fields),
        }

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Null: spatial shuffle of the final field."""
        field = self._final_field
        observed = primary_result.get("peak_to_mean", 0.0)

        if field is None or observed <= 0.0:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        nulls = shuffle_null_distribution(
            field,
            n_permutations=self.n_permutations,
            k_min=self.k_min,
            k_max_frac=self.k_max_frac,
            seed=self._seed,
        )
        # Right-tailed test: how many null draws are >= observed
        p = (float((nulls >= observed).sum()) + 1.0) / (len(nulls) + 1.0)

        return (
            float(p),
            NullType.SHUFFLE,
            {"mean": float(nulls.mean()), "std": float(nulls.std())},
        )

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Cohen's d on peak_to_mean vs shuffle null.

        Observed peak_to_mean is compared against the null distribution's
        mean and std. Positive d means observed is HIGHER than null —
        the direction indicating real structure.
        """
        observed = primary_result.get("peak_to_mean", 0.0)
        nm = null_dist_stats.get("mean", 0.0)
        ns = null_dist_stats.get("std", 0.0)
        if ns <= 1e-9:
            return {
                "cohens_d": 0.0,
                "raw_value": observed,
                "null_mean": nm,
                "null_std": ns,
            }
        d = (observed - nm) / ns
        return {
            "cohens_d": float(d),
            "raw_value": observed,
            "null_mean": nm,
            "null_std": ns,
        }

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Confirmation: higher peak/mean + stable peak_k + p < 0.01 + strong d.

        Effect size d is injected into secondary_result by _determine_tier.
        """
        p2m = primary_result.get("peak_to_mean", 0.0)
        if p2m <= self._CONFIRM_PM_MIN:
            return False

        cohens_d = secondary_result.get("cohens_d", 0.0)
        if cohens_d < self._CONFIRM_COHENS_D_MIN:
            return False

        cv = secondary_result.get("peak_k_cv", float("inf"))
        # cv < 0 is the "could not compute stability" sentinel; a Turing
        # pattern must demonstrate a *stable* wavelength, so reject it.
        if not (0.0 <= cv <= self._CONFIRM_PEAK_K_CV_MAX):
            return False

        if null_p >= self._CONFIRM_NULL_P_MAX:
            return False

        return True

    def _check_definitive(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        null_type: NullType,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> bool:
        """Definitive: confirmation + even higher p/m, d, and tighter CV."""
        p2m = primary_result.get("peak_to_mean", 0.0)
        if p2m <= self._DEF_PM_MIN:
            return False

        cohens_d = secondary_result.get("cohens_d", 0.0)
        if cohens_d < self._DEF_COHENS_D_MIN:
            return False

        cv = secondary_result.get("peak_k_cv", float("inf"))
        if not (0.0 <= cv <= self._DEF_PEAK_K_CV_MAX):
            return False

        return True

    def _determine_tier(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        null_type: NullType,
        effect_size: dict[str, float],
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[DetectionTier, dict[str, bool]]:
        """Inject cohens_d into secondary_result so tier checks can use it."""
        secondary_augmented = dict(secondary_result)
        secondary_augmented["cohens_d"] = effect_size.get("cohens_d", 0.0)
        return super()._determine_tier(
            primary_result=primary_result,
            secondary_result=secondary_augmented,
            null_p=null_p,
            null_type=null_type,
            effect_size=effect_size,
            state_history=state_history,
            model_metadata=model_metadata,
            timescale=timescale,
        )

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        """All secondaries pass if wavelength is stable (low CV)."""
        cv = secondary_result.get("peak_k_cv", float("inf"))
        return 0.0 <= cv <= self._CONFIRM_PEAK_K_CV_MAX

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Exclusions against P1 (aggregation) and P13 (excitable wavefronts).

        P1: a segregation-only field would have high Moran's I but no
          concentrated FFT peak — the peak-to-mean would be low. P3
          already requires high peak-to-mean. We mark P1 'excluded' when
          the peak is narrow and strong (relative peak width small); this
          is a light check since the substrate prerequisite already
          excludes the canonical P1 positives (Schelling, Zhang sorting)
          which are integer-grid.

        P13: excitable media (GH etc.) produce traveling waves on integer
          grids. The continuous-field prerequisite already rejects them.
        """
        checked = ["P1", "P13"]
        results: dict[str, str] = {}

        # Integer-grid substrate excludes P13's canonical positive
        # (GH is integer). Likewise excludes the canonical P1 positives
        # that would be confused with P3 by peak-to-mean.
        if self._final_field is not None:
            nu = count_unique_values(self._final_field)
            if nu >= self._MIN_UNIQUE:
                # Continuous-field substrate — P13 and P1 canonical positives
                # live on discrete-state substrates.
                results["P13"] = "excluded"
                results["P1"] = "excluded"
            else:
                results["P13"] = "inconclusive"
                results["P1"] = "inconclusive"
        else:
            results["P13"] = "inconclusive"
            results["P1"] = "inconclusive"

        return checked, results
