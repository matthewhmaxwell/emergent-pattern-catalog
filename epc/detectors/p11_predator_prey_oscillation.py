"""P11 — Predator-prey (bilateral) oscillation detector.

Detects the signature of a bilateral two-species predator-prey
relationship: strong anti-phase oscillation between two species'
densities, mediated by a shared empty reservoir (so their sum is
not trivially conserved). Canonical positive model: Lotka-Volterra
lattice (Mobilia-Georgiev-Täuber 2007).

Observable scope: state_history_only (species labels on a grid).

DISCRIMINATION FROM NEIGHBORS:

P12 (cyclic dominance): needs >= 3 competing species in a cyclic
  (intransitive) structure. P11 requires exactly 2 species — any
  system with 3 or more candidate species fails P11's prerequisite.
  The two detectors are mutually exclusive by n_species.

P9 (synchronization): intra-species coherence (e.g., Kuramoto phase
  alignment). P11 is an INTER-species anti-phase signal. P11's primary
  metric (rho_anti < 0) is irrelevant to P9 (which measures order
  parameters > 0), and vice versa.

P1 (aggregation): can co-occur with P11, since the spatial clusters
  that predators and prey form during pursuit-evasion waves produce
  nonzero Moran's I. The Sprint 11 characterization measured LV density
  Moran's I ≈ 0.28 (predator) / 0.60 (prey). LV × P1 is expected at
  screening level.

MECHANISM:

The lattice LV model (A = predator, B = prey, E = empty):
  A → E     at rate μ        (predator death)
  B + E → B + B  at rate σ   (prey birth into empty)
  A + B → A + A  at rate λ   (predation with offspring replacing prey)

Resonant amplification of demographic noise in finite-size lattices
produces long-lived erratic oscillations with period ~50-400 generations
(depending on L and rates). Global species density amplitudes scale
as O(1/√L²) but the power-spectrum peak-to-mean ratio remains large
(typically > 15) throughout the coexistence phase.

Predator density LAGS prey density, in the sense that the cross-
correlation Pearson(prey(t), predator(t+tau)) has its MINIMUM value
(most anti-correlated) at a nonzero lag |tau| ~ 10-20 generations
(roughly half an oscillation period shifted by quarter-period phase
lag). At lag 0 the correlation is moderately negative; the minimum
deepens at small positive or negative lag, reflecting the phase-shifted
anti-phase relationship.

DETECTOR TIERS:

All tiers require prerequisites:
  - Exactly 2 non-empty species on the grid (from state_history)
  - Both species have nontrivial time variance (std > 0.005)
  - Species sum has nontrivial variance (std > 0.005) — this catches
    strictly-conserved 2-species systems (Nowak-May: coop + defect = 1)

Primary metric: rho_anti = min_{|tau| >= 5} Pearson(A(t), B(t+tau))
Secondary: FFT peak-to-mean of species B density; |tau_anti|; rho at lag 0.

Screening:    prerequisites pass AND rho_anti < -0.3 AND |tau_anti| >= 5
              AND fft_peak_to_mean > 6.
Confirmation: screening AND rho_anti < -0.5 AND null p < 0.05.
Definitive:   confirmation AND rho_anti < -0.7 AND fft_peak_to_mean > 12
              AND null effect size |cohens_d| > 1.0
              AND P12 exclusion cleared (n_species != 3 → auto-pass).

Null model: circular-shift of one species' time series. Preserves each
series' marginal, autocorrelation and power spectrum; destroys the
phase relationship between the two series. This is a SHUFFLE-type null.

The shuffle null is INTENTIONALLY STRONG for oscillating systems: it
preserves the slow-mode autocorrelation that makes LV and RPS look
similar in cross-correlation, so the p-value alone does not separate
them. Separation is achieved by the rho_anti MAGNITUDE (LV |rho| ~ 0.7-0.9)
and the n_species prerequisite (RPS has 3+ species and fails).

EMPIRICAL CALIBRATION (Sprint 11):

Canonical positive (LV L=100 lambda=4 sigma=mu=1, seed=42, 1500 steps,
burn-in 100): rho_anti ~= -0.86, fft_peak_to_mean ~= 25, null p ~= 0.07,
cohens_d ~= -2.0.

Seed stability (seeds 42, 7, 123): rho_anti ∈ [-0.86, -0.72],
fft_peak_to_mean ∈ [13.6, 21.6].

Negative-model spot checks:
  Schelling (tau=0.375):         prereq fail (species_std = 0.000).
  Nowak-May (b=1.8):             prereq fail (total_std = 0.000, strict
                                 conservation).
  SIR (inf=0.3, rec=0.1):        rho_anti ~= -0.35 (screening fail).
  RPS species pairs (M=1e-4):    n_species prereq fails (3 species).
  White noise:                   rho_anti ~= -0.08 (screening fail).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, DetectorResult, NullType, compute_confidence
from epc.metrics.predator_prey_crosscorr import (
    circular_shift_null,
    extract_species_fractions,
    fft_peak_to_mean,
    predator_prey_rho_anti,
    species_time_series_variance_check,
)


class P11PredatorPreyDetector(BaseDetector):
    """Detector for P11 — bilateral predator-prey oscillation.

    Parameters
    ----------
    n_permutations : int
        Number of circular-shift null permutations (default 199, giving
        p-resolution of ~0.005).
    max_lag : int
        Maximum absolute lag to search for rho_anti. Default 150 — enough
        to capture LV oscillation periods seen in the characterization.
    min_abs_lag : int
        Minimum |lag| considered when searching for rho_anti. Excludes
        near-instantaneous anti-correlation (which would be a conservation
        artifact). Default 5, based on LV characterization (|tau_anti|
        observed in [11, 18]).
    burn_in : int
        Number of initial timesteps to discard before computing metrics.
        Default 100 — covers the initial deterministic-like transient of
        LV without wasting much of shorter trajectories.
    species_state_hint : tuple[int, int] or None
        Optional (prey_state, predator_state) integer labels if known
        in advance (e.g., from model_metadata). If None, the detector
        auto-detects the two most-common non-zero states observed in
        the trajectory.
    seed : int
        RNG seed for the null model.
    """

    # Screening thresholds
    _SCREEN_RHO_MAX = -0.3        # rho_anti must be < this (more negative)
    _SCREEN_ABS_TAU_MIN = 5
    _SCREEN_FFT_MIN = 6.0

    # Confirmation thresholds (in addition to screening)
    _CONFIRM_RHO_MAX = -0.5
    _CONFIRM_COHENS_D_MAX = -1.5  # observed rho more negative than null mean by >=1.5 SD

    # Definitive thresholds (in addition to confirmation)
    _DEF_RHO_MAX = -0.7
    _DEF_FFT_MIN = 12.0
    _DEF_COHENS_D_MAX = -1.5  # same as confirmation; deepens separation via rho+fft
    # Note: null_p is reported as a diagnostic but NOT used as a tier gate.
    # Rationale (Sprint 11 empirics): the circular-shift null preserves each
    # series' autocorrelation, so on LV the null frequently produces its own
    # extreme anti-correlations (p ~ 0.05-0.15 even when observed d ~ -2).
    # Separation LV-vs-noise is achieved by rho_anti magnitude (|LV|~0.8,
    # |noise|~0.1), and separation LV-vs-RPS by the n_species prerequisite.

    # Prerequisite variance floors
    _MIN_SPECIES_STD = 0.005
    _MIN_TOTAL_STD = 0.005

    def __init__(
        self,
        n_permutations: int = 199,
        max_lag: int = 150,
        min_abs_lag: int = 5,
        burn_in: int = 100,
        species_state_hint: tuple[int, int] | None = None,
        seed: int = 42,
    ) -> None:
        super().__init__(
            pattern_id="P11",
            excluded_patterns=["P9", "P12"],
            allowed_co_occurrences=["P1"],
            observable_scope="state_history_only",
        )
        if n_permutations < 1:
            raise ValueError("n_permutations must be >= 1")
        if max_lag < 1:
            raise ValueError("max_lag must be >= 1")
        if min_abs_lag < 0:
            raise ValueError("min_abs_lag must be >= 0")
        if burn_in < 0:
            raise ValueError("burn_in must be >= 0")

        self.n_permutations = n_permutations
        self.max_lag = max_lag
        self.min_abs_lag = min_abs_lag
        self.burn_in = burn_in
        self.species_state_hint = species_state_hint
        self._seed = seed

        # Cache of inferred species states (populated during detect())
        self._resolved_species: tuple[int, int] | None = None
        self._n_unique_species_observed: int = 0

    # ------------------------------------------------------------------
    # Species identification
    # ------------------------------------------------------------------

    def _resolve_species_states(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> tuple[int, int] | None:
        """Determine which grid integer labels correspond to the two species.

        Priority:
          1. self.species_state_hint if provided.
          2. model_metadata hints: if metadata has explicit keys
             'prey_state' / 'predator_state' (or similar), use them.
          3. Auto-detect: take the two most-populous nonzero integer
             labels observed across the trajectory. Record the total
             count of distinct nonzero labels as a prerequisite stat.

        Returns (prey_state, predator_state), or None if cannot resolve
        (e.g., fewer than 2 nonzero species observed).

        The AUTO-DETECTION does not depend on which label is actually
        prey vs predator: rho_anti is symmetric in its two arguments
        (x <-> y at lag tau equals y <-> x at lag -tau). The labels
        are purely identifiers.
        """
        if self.species_state_hint is not None:
            return self.species_state_hint

        # Inspect metadata for explicit hints
        if model_metadata:
            # Common metadata keys
            if (
                isinstance(model_metadata.get("prey_state"), int)
                and isinstance(model_metadata.get("predator_state"), int)
            ):
                return (
                    int(model_metadata["prey_state"]),
                    int(model_metadata["predator_state"]),
                )

        # Auto-detect: count occurrences of each nonzero label
        counts: dict[int, int] = {}
        for s in state_history:
            grid = s.get("grid")
            if grid is None:
                continue
            vals, cnts = np.unique(grid, return_counts=True)
            for v, c in zip(vals, cnts):
                iv = int(v)
                if iv == 0:
                    continue
                counts[iv] = counts.get(iv, 0) + int(c)

        self._n_unique_species_observed = len(counts)

        if len(counts) < 2:
            return None
        if len(counts) > 2:
            # P11 requires exactly 2 species; >= 3 means the prerequisite
            # will fail. But still return the top 2 for diagnostics.
            sorted_keys = sorted(counts.keys(), key=lambda k: -counts[k])
            return (sorted_keys[0], sorted_keys[1])

        # Exactly 2 species
        keys = sorted(counts.keys())
        return (keys[0], keys[1])

    # ------------------------------------------------------------------
    # BaseDetector overrides
    # ------------------------------------------------------------------

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """Timescale: LV pursuit-evasion period ~ max(rows, cols).

        Consistent with other lattice_2d detectors; a finer estimate
        would come from the FFT dominant period, but we don't have
        that before _compute_primary runs.
        """
        if state_history and "grid_dims" in state_history[0]:
            rows, cols = state_history[0]["grid_dims"]
            return float(max(rows, cols))
        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        warnings = super()._validate_prerequisites(state_history, timescale)
        if not state_history:
            warnings.append("empty state_history")
            return warnings

        s0 = state_history[0]
        if "grid" not in s0:
            warnings.append("no 'grid' key in state_history — P11 needs grids")
            return warnings

        # Run length sanity: need enough samples post-burn-in for a stable
        # rho_anti (rule of thumb: at least 3x max_lag + burn_in)
        needed = self.burn_in + 3 * self.max_lag
        if len(state_history) < needed:
            warnings.append(
                f"run length {len(state_history)} < recommended {needed} "
                f"(burn_in {self.burn_in} + 3*max_lag {3*self.max_lag})"
            )

        return warnings

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Primary metric: rho_anti + FFT peak-to-mean + variance checks."""
        # Resolve species labels
        species = self._resolve_species_states(state_history, None)
        self._resolved_species = species

        if species is None:
            return {
                "rho_anti": 0.0,
                "tau_anti": 0,
                "rho_at_zero_lag": 0.0,
                "fft_peak_to_mean": 0.0,
                "dominant_period": -1.0,
                "prey_std": 0.0,
                "predator_std": 0.0,
                "total_std": 0.0,
                "n_unique_species_observed": self._n_unique_species_observed,
                "n_samples": 0,
            }

        prey_state, pred_state = species

        # Extract time series
        prey_full, pred_full = extract_species_fractions(
            state_history,
            prey_state=prey_state,
            predator_state=pred_state,
        )

        # Apply burn-in
        bi = min(self.burn_in, max(0, len(prey_full) - 10))
        prey_ss = prey_full[bi:]
        pred_ss = pred_full[bi:]

        # Variance prerequisites (these are part of primary metric so that
        # the detector can report them even if screening fails)
        var_stats = species_time_series_variance_check(prey_ss, pred_ss)

        # rho_anti: anti-correlation at nonzero lag
        rho_stats = predator_prey_rho_anti(
            prey_ss, pred_ss,
            max_lag=self.max_lag,
            min_abs_lag=self.min_abs_lag,
        )

        # FFT peak-to-mean on predator density (the more informative of
        # the two — predator densities are typically lower-amplitude and
        # therefore have cleaner relative spectra)
        fft_stats = fft_peak_to_mean(pred_ss)

        return {
            **rho_stats,
            **fft_stats,
            **var_stats,
            "n_unique_species_observed": self._n_unique_species_observed,
            "prey_state": int(prey_state),
            "predator_state": int(pred_state),
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening gates:

        1. n_unique_species_observed == 2 (exactly 2 non-zero species).
        2. prey_std and predator_std both >= _MIN_SPECIES_STD.
        3. total_std (std of prey+predator) >= _MIN_TOTAL_STD.
        4. rho_anti < _SCREEN_RHO_MAX (= -0.3).
        5. |tau_anti| >= _SCREEN_ABS_TAU_MIN.
        6. fft_peak_to_mean > _SCREEN_FFT_MIN.
        """
        n_sp = primary_result.get("n_unique_species_observed", 0)
        if n_sp != 2:
            return False

        if primary_result.get("prey_std", 0.0) < self._MIN_SPECIES_STD:
            return False
        if primary_result.get("predator_std", 0.0) < self._MIN_SPECIES_STD:
            return False
        if primary_result.get("total_std", 0.0) < self._MIN_TOTAL_STD:
            return False

        rho_anti = primary_result.get("rho_anti", 0.0)
        if rho_anti >= self._SCREEN_RHO_MAX:
            return False

        tau_anti = primary_result.get("tau_anti", 0)
        if abs(int(tau_anti)) < self._SCREEN_ABS_TAU_MIN:
            return False

        p2m = primary_result.get("fft_peak_to_mean", 0.0)
        if p2m <= self._SCREEN_FFT_MIN:
            return False

        return True

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Secondary metrics:
          - prey FFT peak-to-mean (separate from the primary predator p2m)
          - period agreement (prey dominant period vs predator)
          - stability: half-and-half robustness of rho_anti sign
        """
        if self._resolved_species is None:
            return {}

        prey_state, pred_state = self._resolved_species
        prey_full, pred_full = extract_species_fractions(
            state_history,
            prey_state=prey_state,
            predator_state=pred_state,
        )

        bi = min(self.burn_in, max(0, len(prey_full) - 10))
        prey_ss = prey_full[bi:]
        pred_ss = pred_full[bi:]

        if len(prey_ss) < 30:
            return {}

        # Prey FFT
        prey_fft = fft_peak_to_mean(prey_ss)

        # Half-trajectory robustness: compute rho_anti on first and second
        # halves and check both are negative.
        h = len(prey_ss) // 2
        r1 = predator_prey_rho_anti(
            prey_ss[:h], pred_ss[:h],
            max_lag=min(self.max_lag, max(10, h // 3)),
            min_abs_lag=self.min_abs_lag,
        )
        r2 = predator_prey_rho_anti(
            prey_ss[h:], pred_ss[h:],
            max_lag=min(self.max_lag, max(10, (len(prey_ss) - h) // 3)),
            min_abs_lag=self.min_abs_lag,
        )
        halves_agree = (
            r1.get("rho_anti", 0.0) < self._SCREEN_RHO_MAX
            and r2.get("rho_anti", 0.0) < self._SCREEN_RHO_MAX
        )

        return {
            "prey_fft_peak_to_mean": prey_fft["fft_peak_to_mean"],
            "prey_dominant_period": prey_fft["dominant_period"],
            "rho_anti_first_half": r1["rho_anti"],
            "rho_anti_second_half": r2["rho_anti"],
            "halves_agree": bool(halves_agree),
        }

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Null model: circular shift of predator time series.

        Preserves each series' autocorrelation and power spectrum; destroys
        cross-correlation. See module docstring for design rationale.
        """
        if self._resolved_species is None:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        prey_state, pred_state = self._resolved_species
        prey_full, pred_full = extract_species_fractions(
            state_history,
            prey_state=prey_state,
            predator_state=pred_state,
        )

        bi = min(self.burn_in, max(0, len(prey_full) - 10))
        prey_ss = prey_full[bi:]
        pred_ss = pred_full[bi:]

        if len(prey_ss) < 3 * self.min_abs_lag + 3:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        null = circular_shift_null(
            prey_ss, pred_ss,
            n_permutations=self.n_permutations,
            max_lag=self.max_lag,
            min_abs_lag=self.min_abs_lag,
            seed=self._seed,
        )

        return (
            float(null["p_value"]),
            NullType.SHUFFLE,
            {"mean": null["null_mean"], "std": null["null_std"]},
        )

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Cohen's d on rho_anti vs null distribution.

        Observed more NEGATIVE than null → d negative. We take the absolute
        value for tier checks.
        """
        observed = primary_result.get("rho_anti", 0.0)
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
        """Confirmation: rho_anti <= -0.5 AND cohens_d <= -1.5 AND halves agree.

        Note: null_p is intentionally NOT a gate. See class docstring on
        _DEF_COHENS_D_MAX — the circular-shift null preserves autocorrelation
        so p-values are not a clean signal-vs-noise separator. The real
        separator is rho_anti magnitude (|LV|~0.8 vs |noise|~0.1).
        """
        rho = primary_result.get("rho_anti", 0.0)
        if rho > self._CONFIRM_RHO_MAX:
            return False
        # Effect size is injected into secondary_result by _determine_tier
        cohens_d = secondary_result.get("cohens_d", 0.0)
        if cohens_d > self._CONFIRM_COHENS_D_MAX:
            return False
        # Halves-agree is a robustness secondary; fine to require.
        if not secondary_result.get("halves_agree", False):
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
        """Definitive: confirmation + rho_anti <= -0.7 + fft_p2m > 12
        + cohens_d <= -2.0 + P12/P9 exclusions cleared (auto-pass via
        n_species == 2 prerequisite).
        """
        rho = primary_result.get("rho_anti", 0.0)
        if rho > self._DEF_RHO_MAX:
            return False
        p2m = primary_result.get("fft_peak_to_mean", 0.0)
        if p2m <= self._DEF_FFT_MIN:
            return False
        cohens_d = secondary_result.get("cohens_d", 0.0)
        if cohens_d > self._DEF_COHENS_D_MAX:
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
        """Inject effect_size cohens_d into secondary_result so tier checks
        have access, then delegate to the base implementation."""
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
        return bool(secondary_result.get("halves_agree", False))

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Check P9 and P12 exclusions.

        P12 (cyclic dominance): requires >= 3 species. P11's prerequisite
          n_unique_species_observed == 2 auto-clears this.

        P9 (synchronization): intra-species coherence — not meaningfully
          applicable to our two-species density time series. We mark
          it 'excluded' if primary n_species == 2 (there's no phase
          variable to synchronize), 'inconclusive' otherwise.
        """
        checked = ["P9", "P12"]
        results: dict[str, str] = {}

        n_sp = self._n_unique_species_observed

        if n_sp == 2:
            results["P12"] = "excluded"
            results["P9"] = "excluded"
        elif n_sp >= 3:
            # Shouldn't reach here if screening passed, but defensively:
            results["P12"] = "not_excluded"
            results["P9"] = "inconclusive"
        else:
            results["P12"] = "inconclusive"
            results["P9"] = "inconclusive"

        return checked, results
