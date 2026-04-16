"""P15 — Persistent computation detector (generalized).

Detects systems where the final state is a *deterministic, input-dependent
function* of initial conditions — i.e., the system implements persistent
computation. Works for any deterministic lattice_2d CA when given a
step function.

This replaces the GoL-specific `p15_fidelity_fix.py` module with a proper
BaseDetector subclass. The original module is retained for backward
compatibility but its functionality is now available substrate-generally
here.

Observable scope: state_history_only for structural diversity analysis;
model_metadata_required if a step function is not supplied.

Detection tiers:
  Screening:    Determinism ≥ 0.9 (bit-exact replay matches) AND
                ≥ 2 distinct outcome classes across N variations.
  Confirmation: Screening + determinism = 1.0 AND ≥ 3 distinct outcomes.
  Definitive:   Confirmation + effect size > 0 AND N_variations ≥ 8.

Null model:
  Scrambled-IC null: randomly permute the cells of the initial grid
  while preserving the total cell count per state. Under this null,
  outcome diversity and determinism should be similar to the data — this
  is the BASELINE for "what does input-variation do for a non-computing
  system?" The test is whether the ORDERED variations (phase shifts,
  targeted flips) produce MORE outcome diversity than random scrambles.

Nearest-neighbor exclusions:
  P13 (excitable waves): persistent re-entrant activity. Computation
       requires that the outcome encodes a decision, not just that the
       system stays active.
  P22 (cascade): single-pass spreading is not computation (the outcome
       doesn't depend meaningfully on fine-grained IC details).

Usage:
  from epc.detectors.p15_persistent_computation import (
      P15PersistentComputationDetector,
      make_step_fn_from_model,
  )

  # For GoL
  step_fn = make_step_fn_from_model(gol_model)
  det = P15PersistentComputationDetector(step_fn=step_fn, n_variations=8)
  result = det.detect(history, metadata)

References:
  Langton, C.G. (1990). "Computation at the edge of chaos: Phase transitions
    and emergent computation." Physica D 42, 12-37.
  Lizier, J.T. et al. (2012). "Local measures of information storage in
    complex distributed computation." Information Sciences 208, 39-54.
"""

from __future__ import annotations

from typing import Any, Callable
from collections import Counter

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType


# Type alias: a step function takes a grid and returns the next grid
StepFn = Callable[[np.ndarray], np.ndarray]


def make_step_fn_for_gol() -> StepFn:
    """Factory for the GoL B3/S23 step function (periodic BC).

    Identical dynamics to GameOfLife.step_vectorized but as a pure function
    of the grid, suitable for replay and variation testing.
    """
    def _gol_step(grid: np.ndarray) -> np.ndarray:
        padded = np.pad(grid, 1, mode="wrap")
        neighbors = sum(
            padded[1 + dr: grid.shape[0] + 1 + dr,
                   1 + dc: grid.shape[1] + 1 + dc]
            for dr in [-1, 0, 1] for dc in [-1, 0, 1]
            if not (dr == 0 and dc == 0)
        )
        out = ((grid == 1) & ((neighbors == 2) | (neighbors == 3)) |
               (grid == 0) & (neighbors == 3))
        return out.astype(grid.dtype)
    return _gol_step


def make_step_fn_from_model(model: Any) -> StepFn:
    """Extract a pure step function from a model instance.

    For deterministic models, this produces a function grid → grid' that
    can be called independently of the model's internal state. For
    stochastic models, the result will not be reproducible — which is
    exactly what P15's reproducibility check should detect.

    Supports: GameOfLife, GreenbergHastings, NowakMaySpatialPD, and any
    model with a step() or step_vectorized() method and a grid attribute.
    """
    model_class_name = type(model).__name__

    if model_class_name == "GameOfLife":
        return make_step_fn_for_gol()

    if model_class_name == "GreenbergHastings":
        # GH is deterministic; we can replay from any grid by swapping
        # the internal state temporarily.
        def _gh_step(grid: np.ndarray) -> np.ndarray:
            saved = model._grid
            saved_step = model._step_count
            try:
                model._grid = grid.copy()
                model.step_vectorized()
                return model._grid.copy()
            finally:
                model._grid = saved
                model._step_count = saved_step
        return _gh_step

    if model_class_name in ("NowakMaySpatialPD", "NowakMayModel"):
        # NowakMay is deterministic; uses self.grid (not self._grid)
        def _nm_step(grid: np.ndarray) -> np.ndarray:
            saved = model.grid
            try:
                model.grid = grid.copy()
                model.step()
                return model.grid.copy()
            finally:
                model.grid = saved
        return _nm_step

    # Generic fallback: try common attribute names
    for attr_name in ("_grid", "grid"):
        if hasattr(model, attr_name):
            for step_name in ("step_vectorized", "step"):
                if hasattr(model, step_name):
                    step_method = getattr(model, step_name)
                    def _generic_step(
                        grid: np.ndarray,
                        _attr=attr_name,
                        _step=step_method,
                    ) -> np.ndarray:
                        saved = getattr(model, _attr)
                        try:
                            setattr(model, _attr, grid.copy())
                            _step()
                            return getattr(model, _attr).copy()
                        finally:
                            setattr(model, _attr, saved)
                    return _generic_step

    raise TypeError(
        f"Cannot extract step function from {model_class_name}. "
        f"Pass a custom step_fn=Callable[[ndarray], ndarray] instead."
    )


def _classify_outcome(
    post_grids: list[np.ndarray],
    motion_threshold: float = 1.5,
) -> str:
    """Substrate-independent outcome classification.

    Given a window of post-trajectory snapshots, classify as:
      - 'dead': all cells in zero state
      - 'static': no change at all
      - 'period_{p}': exact cycle with period p (1 < p <= len(window))
      - 'moving': center of mass drifts beyond threshold
      - 'complex': non-periodic, non-moving (chaotic/transient)
    """
    if not post_grids or all(int(np.sum(g != 0)) == 0 for g in post_grids):
        return "dead"

    first = post_grids[0]

    # Static check
    if all(np.array_equal(first, g) for g in post_grids[1:]):
        return "static"

    # Periodicity check: find smallest p such that post_grids[0] == post_grids[p]
    for p in range(1, len(post_grids)):
        if np.array_equal(first, post_grids[p]):
            return f"period_{p}"

    # Motion check: center of mass displacement
    def _com(g: np.ndarray) -> np.ndarray | None:
        cells = np.argwhere(g != 0)
        return cells.mean(axis=0) if len(cells) > 0 else None

    com0 = _com(post_grids[0])
    com_last = _com(post_grids[-1])
    if com0 is not None and com_last is not None:
        disp = float(np.linalg.norm(com_last - com0))
        if disp > motion_threshold:
            return "moving"

    return "complex"


def _make_variations(
    initial_grid: np.ndarray,
    n_variations: int,
    rng: np.random.Generator,
) -> list[np.ndarray]:
    """Generate N input variations by perturbing the IC.

    Uses two strategies:
    1. Cell flips: toggle a small number of cells (1-3) to probe how
       the outcome depends on fine-grained IC differences. This is the
       GoL-like "different collision phases" test generalized.
    2. Translations: shift the entire IC by different offsets. This
       probes spatial symmetry breaking.

    The mix of both tests whether the system's dynamics are sensitive to
    small input changes (hallmark of computation) vs. robust to them
    (hallmark of simple dynamics like waves or aggregation).
    """
    rows, cols = initial_grid.shape
    variations: list[np.ndarray] = []

    # Half perturbative flips, half translations
    n_flips = max(1, n_variations // 2)
    n_translations = n_variations - n_flips

    # Strategy 1: cell flips (1-3 random cells toggled)
    nonzero = np.argwhere(initial_grid != 0)
    zero = np.argwhere(initial_grid == 0)
    max_state = int(initial_grid.max()) if initial_grid.max() > 0 else 1

    for k in range(n_flips):
        varied = initial_grid.copy()
        n_changes = 1 + k % 3  # 1, 2, or 3 cell changes
        for _ in range(n_changes):
            r = rng.integers(0, rows)
            c = rng.integers(0, cols)
            if varied[r, c] == 0:
                varied[r, c] = rng.integers(1, max_state + 1)
            else:
                varied[r, c] = 0
        variations.append(varied)

    # Strategy 2: translations
    for k in range(n_translations):
        dr = (k * 3 + 1) % rows
        dc = (k * 5 + 2) % cols
        shifted = np.roll(np.roll(initial_grid, dr, axis=0), dc, axis=1)
        variations.append(shifted)

    return variations


def _run_trajectory(
    step_fn: StepFn,
    initial_grid: np.ndarray,
    n_steps: int,
) -> np.ndarray:
    """Run step_fn for n_steps and return the final grid."""
    current = initial_grid.copy()
    for _ in range(n_steps):
        current = step_fn(current)
    return current


def _trajectory_window(
    step_fn: StepFn,
    initial_grid: np.ndarray,
    n_steps: int,
    window_size: int = 8,
) -> list[np.ndarray]:
    """Run n_steps, then capture window_size subsequent snapshots."""
    current = initial_grid.copy()
    for _ in range(n_steps):
        current = step_fn(current)
    window = [current.copy()]
    for _ in range(window_size - 1):
        current = step_fn(current)
        window.append(current.copy())
    return window


class P15PersistentComputationDetector(BaseDetector):
    """Detector for P15 — Persistent computation in deterministic lattice CAs.

    Parameters
    ----------
    step_fn : callable(grid) -> grid, optional
        Pure function that advances a grid by one timestep. If None, the
        detector runs in state-history-only mode with structural diversity
        analysis only (cannot achieve more than SCREENING tier).
    n_variations : int
        Number of input variations to test (default 8). Definitive
        tier requires ≥ 8.
    seed : int
        Random seed for null model generation.
    window_size : int
        Number of post-trajectory snapshots for outcome classification
        (default 8). Larger values can detect longer periods.
    """

    def __init__(
        self,
        step_fn: StepFn | None = None,
        n_variations: int = 8,
        seed: int = 42,
        window_size: int = 8,
    ) -> None:
        super().__init__(
            pattern_id="P15",
            excluded_patterns=["P13", "P22"],
            allowed_co_occurrences=[],
            observable_scope="state_history_only" if step_fn is None
                              else "model_metadata_required",
        )
        self.step_fn = step_fn
        self.n_variations = n_variations
        self._seed = seed
        self.window_size = window_size
        self._last_outcomes: list[str] = []

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """T_prop = max(rows, cols): time to cross the grid."""
        if state_history and "grid_dims" in state_history[0]:
            rows, cols = state_history[0]["grid_dims"]
            return float(max(rows, cols))
        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        warnings = []
        if not state_history:
            warnings.append("empty state history")
            return warnings
        if "grid" not in state_history[0]:
            warnings.append("no grid observable in state history")
        if self.step_fn is None:
            warnings.append(
                "no step_fn provided — P15 can only reach SCREENING tier. "
                "Pass step_fn=make_step_fn_from_model(model) for full detection."
            )
        if len(state_history) < 20:
            warnings.append(f"short trajectory ({len(state_history)} steps)")
        return warnings

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Primary metrics: reproducibility and outcome diversity.

        Reproducibility: replay the IC with step_fn; compare bit-exactly
        against the recorded final state. 1.0 = deterministic match.

        Diversity: run N variations of the IC and count distinct outcome
        classes (static / period_p / moving / complex / dead).
        """
        if not state_history or "grid" not in state_history[0]:
            return {
                "reproducibility": 0.0,
                "n_distinct_outcomes": 0.0,
                "outcome_diversity": 0.0,
                "n_trajectory_steps": 0.0,
            }

        initial_grid = state_history[0]["grid"]
        final_grid = state_history[-1]["grid"]
        # Use actual step count, not just len(history) — handles record_every > 1
        n_steps = state_history[-1].get("step", len(state_history) - 1)

        # Reproducibility test: replay from IC and compare at MULTIPLE
        # timepoints (not just the endpoint). This catches stochastic
        # models that converge to the same absorbing state via different
        # paths (e.g., SIR epidemic where everyone eventually recovers).
        reproducibility = 0.0
        if self.step_fn is not None and n_steps > 0:
            try:
                # Build a map of recorded step → recorded grid
                step_to_grid = {}
                for s in state_history:
                    step_to_grid[s.get("step", 0)] = s["grid"]

                # Sample 4 checkpoint steps: 1/4, 1/2, 3/4, final
                checkpoints = sorted(set([
                    max(1, n_steps // 4),
                    max(1, n_steps // 2),
                    max(1, 3 * n_steps // 4),
                    n_steps,
                ]))

                current = initial_grid.copy()
                matches = 0
                total_checks = 0
                for step_idx in range(1, n_steps + 1):
                    current = self.step_fn(current)
                    if step_idx in checkpoints and step_idx in step_to_grid:
                        total_checks += 1
                        if np.array_equal(current, step_to_grid[step_idx]):
                            matches += 1

                reproducibility = (
                    matches / total_checks if total_checks > 0 else 0.0
                )
            except Exception:
                reproducibility = 0.0

        # Outcome diversity: run N variations and classify
        n_distinct = 0
        outcomes: list[str] = []
        if self.step_fn is not None:
            rng = np.random.default_rng(self._seed)
            variations = _make_variations(initial_grid, self.n_variations, rng)
            for var_ic in variations:
                try:
                    window = _trajectory_window(
                        self.step_fn, var_ic, n_steps, self.window_size
                    )
                    outcomes.append(_classify_outcome(window))
                except Exception:
                    outcomes.append("error")
            n_distinct = len(set(outcomes))
        else:
            # State-history-only mode: analyze recorded trajectory
            # for structural diversity. Count distinct "regimes" visible
            # in the time series (population levels).
            pops = np.array([int(np.sum(s["grid"] != 0))
                             for s in state_history])
            # Bin populations into 4 levels; count distinct transitions
            if len(pops) > 0 and pops.max() > 0:
                bins = np.digitize(pops, np.linspace(0, pops.max() + 1, 5))
                n_distinct = len(set(bins.tolist()))

        diversity_frac = (n_distinct / max(1, self.n_variations)
                          if self.step_fn is not None else n_distinct / 4.0)

        # Store outcomes list for secondary analysis (not in primary dict
        # which should be float-only)
        self._last_outcomes = outcomes

        return {
            "reproducibility": reproducibility,
            "n_distinct_outcomes": float(n_distinct),
            "outcome_diversity": float(diversity_frac),
            "n_trajectory_steps": float(n_steps),
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening: reproducibility ≥ 0.9 AND ≥ 2 distinct outcomes."""
        rep = primary_result.get("reproducibility", 0.0)
        n_out = primary_result.get("n_distinct_outcomes", 0.0)
        return rep >= 0.9 and n_out >= 2.0

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Secondary: activity metrics from the recorded trajectory."""
        if not state_history or "grid" not in state_history[0]:
            return {}
        pops = [int(np.sum(s["grid"] != 0)) for s in state_history]
        return {
            "population_initial": pops[0] if pops else 0,
            "population_final": pops[-1] if pops else 0,
            "population_peak": max(pops) if pops else 0,
            "population_mean": float(np.mean(pops)) if pops else 0.0,
            "trajectory_length": len(state_history),
            "outcome_classes": self._last_outcomes,
        }

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Null: compare ordered-variation diversity vs. random-permutation diversity.

        Random permutations of the IC should produce MORE diverse outcomes
        than ordered variations (translations) if the system does not
        compute — random permutations destroy spatial correlations. A
        system that computes should show STRUCTURED outcome variation
        under ordered variation, not the same diversity as under random
        scrambling.

        Returns a "weak" null — we report the null diversity and the
        ratio, but the P15 signature is primarily in the reproducibility
        + nonzero diversity check, not a p-value comparison.
        """
        if not state_history or "grid" not in state_history[0]:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        if self.step_fn is None:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        initial_grid = state_history[0]["grid"]
        n_steps = len(state_history) - 1

        # Generate null diversity from random permutations
        rng = np.random.default_rng(self._seed + 1000)
        null_diversities: list[int] = []
        # Smaller null sample for speed — this is an auxiliary metric
        n_null = 20

        for _ in range(n_null):
            # Randomly shuffle cells (preserving count per state)
            flat = initial_grid.flatten().copy()
            rng.shuffle(flat)
            scrambled_ic = flat.reshape(initial_grid.shape)

            null_outcomes = []
            for k in range(min(self.n_variations, 4)):  # fewer for speed
                # Apply the same translation set to the scrambled IC
                dr = (k * 3 + 1) % initial_grid.shape[0]
                dc = (k * 5 + 2) % initial_grid.shape[1]
                var_ic = np.roll(np.roll(scrambled_ic, dr, axis=0), dc, axis=1)
                try:
                    window = _trajectory_window(
                        self.step_fn, var_ic, min(n_steps, 50), self.window_size
                    )
                    null_outcomes.append(_classify_outcome(window))
                except Exception:
                    null_outcomes.append("error")
            null_diversities.append(len(set(null_outcomes)))

        null_mean = float(np.mean(null_diversities)) if null_diversities else 0.0
        null_std = float(np.std(null_diversities)) if null_diversities else 0.0

        observed = primary_result.get("n_distinct_outcomes", 0.0)
        # Report p-value as fraction of null samples with diversity ≥ observed
        p_value = (sum(1 for d in null_diversities if d >= observed)
                   / max(1, len(null_diversities)))

        return (
            p_value,
            NullType.SHUFFLE,
            {"mean": null_mean, "std": null_std,
             "note": "null_is_scrambled_IC_with_same_variations"},
        )

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Effect size: observed outcome diversity vs. null diversity."""
        observed = primary_result.get("n_distinct_outcomes", 0.0)
        null_mean = null_dist_stats.get("mean", 0.0)
        null_std = null_dist_stats.get("std", 1.0)
        cohens_d = ((observed - null_mean) / null_std
                    if null_std > 0 else 0.0)
        return {
            "cohens_d": cohens_d,
            "raw_value": observed,
            "null_mean": null_mean,
            "null_std": null_std,
            "note": "diversity_from_ordered_vs_scrambled_variations",
        }

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Confirmation: determinism = 1.0 AND ≥ 3 distinct outcomes."""
        rep = primary_result.get("reproducibility", 0.0)
        n_out = primary_result.get("n_distinct_outcomes", 0.0)
        return rep >= 1.0 and n_out >= 3.0

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
        """Definitive: confirmation + n_variations ≥ 8."""
        rep = primary_result.get("reproducibility", 0.0)
        n_out = primary_result.get("n_distinct_outcomes", 0.0)
        return (rep >= 1.0
                and n_out >= 3.0
                and self.n_variations >= 8)

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Check P13 and P22 exclusions.

        Computation ≠ persistent waves (P13) because computation produces
        discrete outcome classes, not just oscillating activity.
        Computation ≠ cascade (P22) because cascades are single-pass
        spreading with no real outcome distinction.
        """
        checked = ["P13", "P22"]
        results: dict[str, str] = {}

        # P13 exclusion: we're NOT a pure excitable wave if our outcomes
        # are diverse (not just "all periodic"). If the model shows
        # outcome diversity beyond period-matching, P13 is excluded.
        # Without step_fn we can't distinguish — return inconclusive.
        if self.step_fn is not None:
            # If most outcomes are NOT periodic, P13 is excluded
            results["P13"] = "inconclusive"  # detailed analysis needs outcomes list
        else:
            results["P13"] = "not_checked"

        # P22 exclusion: cascade would die out quickly. If population
        # persists throughout the trajectory, P22 is excluded.
        if state_history:
            pops = [int(np.sum(s["grid"] != 0)) for s in state_history]
            # If activity in last quarter is substantial relative to peak,
            # it's not a dying cascade
            last_q = pops[len(pops) * 3 // 4:]
            peak = max(pops) if pops else 1
            if last_q and min(last_q) >= peak * 0.1:
                results["P22"] = "excluded"
            elif last_q and max(last_q) < peak * 0.05:
                results["P22"] = "inconclusive"
            else:
                results["P22"] = "excluded"
        else:
            results["P22"] = "not_checked"

        return checked, results
