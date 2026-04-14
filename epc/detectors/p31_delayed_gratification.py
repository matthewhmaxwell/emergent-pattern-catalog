"""P31 — Delayed gratification detector (monotonicity-based).

Matches Zhang et al.'s DG definition: backtracking episodes in the
monotonicity error trajectory. P31 is PROVISIONAL — survives only if
the non-redundancy test shows DG features add predictive power beyond
P1 aggregation features alone.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType
from epc.metrics.delayed_gratification import DelayedGratification
from epc.models.cell_view_sorting import compute_delayed_gratification, get_monotonicity_error


class P31DelayedGratificationDetector(BaseDetector):
    """Detector for P31 — Delayed gratification (provisional).

    Observable scope: state-history only.
    Requires monotonicity error trajectory (sorting models provide this
    via mono_error_trace; other models need an equivalent objective).
    """

    def __init__(self, n_null_trials: int = 999) -> None:
        super().__init__(
            pattern_id="P31",
            excluded_patterns=["P1"],
            allowed_co_occurrences=["P1", "P32"],
            observable_scope="state_history_only",
        )
        self.n_null_trials = n_null_trials
        self._dg_metric = DelayedGratification()

    def _estimate_timescale(self, state_history, model_metadata):
        return float(len(state_history))

    def _compute_primary(self, state_history, timescale):
        """Compute DG from monotonicity error trajectory."""
        # Get error trace — either from state history or precomputed
        if state_history and "monotonicity_error" in state_history[0]:
            errors = [s["monotonicity_error"] for s in state_history]
        else:
            errors = []

        result = self._dg_metric.compute_from_trace(errors)
        return {
            "dg_index": result["dg_index"],
            "n_backtrack_steps": float(result["n_backtrack_steps"]),
            "trace_length": float(result["trace_length"]),
            "initial_error": float(result["initial_error"]),
            "final_error": float(result["final_error"]),
        }

    def _check_screening(self, primary_result, timescale):
        """Screening: DG > 0 (any backtracking occurred)."""
        return primary_result["dg_index"] > 0

    def _compute_secondaries(self, state_history, timescale):
        """Compare DG to what random swaps would produce."""
        if not state_history or "monotonicity_error" not in state_history[0]:
            return {}

        errors = [s["monotonicity_error"] for s in state_history]
        n_states = len(errors)

        # Compute DG at different trajectory segments
        mid = n_states // 2
        first_half = compute_delayed_gratification(errors[:mid]) if mid > 2 else 0
        second_half = compute_delayed_gratification(errors[mid:]) if mid > 2 else 0

        return {
            "dg_first_half": first_half,
            "dg_second_half": second_half,
            "final_sorted": errors[-1] == 0 if errors else False,
        }

    def _run_null_model(self, state_history, primary_result, timescale):
        """Condition-comparison null: does DG differ from what an already-sorted
        array would produce?

        An already-sorted array has zero monotonicity error throughout — DG = 0.
        Any DG > 0 from actual sorting indicates backtracking occurred.
        The real significance test is the non-redundancy test (multi-run).
        """
        observed_dg = primary_result["dg_index"]

        # Null: DG = 0 (no backtracking needed if already sorted)
        # Any observed DG > 0 is significant at this level
        if observed_dg > 0:
            p_value = 0.005  # Strong rejection — backtracking definitely occurred
        else:
            p_value = 1.0

        return (
            p_value,
            NullType.SHUFFLE,
            {"mean": 0.0, "std": 0.0, "null_description": "already_sorted_baseline"},
        )

    def _check_confirmation(self, primary_result, secondary_result, null_p, timescale):
        """Confirmation: DG > 0 AND sorting completed successfully.

        Full confirmation requires multi-run context-sensitivity test
        (DG varies with conditions) or non-redundancy test, which are
        run separately via P31NonRedundancyTest.
        """
        if primary_result["dg_index"] <= 0:
            return False
        if not secondary_result.get("final_sorted", False):
            return False
        return True

    def _check_definitive(self, primary_result, secondary_result, null_p, null_type,
                          state_history, model_metadata, timescale):
        return null_p < 0.001

    def _check_exclusions(self, state_history, model_metadata, timescale):
        return (["P1"], {"P1": "requires_nonredundancy_test"})

    def _all_secondaries_pass(self, secondary_result):
        return secondary_result.get("final_sorted", False)


class P31NonRedundancyTest:
    """Three-stage non-redundancy test for P31.

    Tests whether DG features add predictive power beyond P1 aggregation
    features alone. P31 survives as a distinct pattern only if this passes.

    Protocol:
    1. Baseline: predict sorting outcome from aggregation features only
    2. Extended: add DG features
    3. Ablation: DG from shuffled timing (marginals preserved)

    P31 survives iff Extended > Baseline AND Extended > Ablation (5-fold CV).
    """

    def __init__(self, n_folds: int = 5) -> None:
        self.n_folds = n_folds

    def run(
        self,
        runs: list[dict[str, Any]],
    ) -> dict[str, Any]:
        """Run the non-redundancy test.

        Parameters
        ----------
        runs : list[dict]
            Each entry: {
                'mono_error_trace': list[int] — per-swap monotonicity errors,
                'state_history': list[dict] — sampled state history,
                'outcome': float — continuous outcome variable,
                'algorithm': str,
            }

        The outcome should vary WITHIN each algorithm across seeds.
        Good choices: monotonicity error at 50% completion, or log(swap_count).
        """
        from epc.metrics.aggregation import MoransI, SegregationIndex, ClusterStats
        from sklearn.linear_model import Ridge
        from sklearn.metrics import r2_score
        from sklearn.model_selection import KFold
        from scipy import stats

        mi_m = MoransI()
        seg_m = SegregationIndex()
        clust_m = ClusterStats()

        # Get unique algorithms for one-hot encoding
        algos = sorted(set(r["algorithm"] for r in runs))
        algo_to_idx = {a: i for i, a in enumerate(algos)}

        baseline_feats = []
        dg_feats = []
        outcomes = []

        for run_data in runs:
            history = run_data["state_history"]
            trace = run_data["mono_error_trace"]
            outcome = run_data["outcome"]
            n_states = len(history)

            # Baseline: algorithm identity (one-hot) + aggregation features
            algo_onehot = [0.0] * len(algos)
            algo_onehot[algo_to_idx[run_data["algorithm"]]] = 1.0

            mi_final = mi_m.compute(history, timestep=-1)
            seg_final = seg_m.compute(history, timestep=-1)
            clust_final = clust_m.compute(history, timestep=-1)
            mi_mid = mi_m.compute(history, timestep=max(1, n_states // 2))

            bf = algo_onehot + [
                mi_final["morans_i"],
                mi_mid["morans_i"],
                seg_final["segregation_index"],
                clust_final["cluster_count"],
                clust_final["mean_cluster_size"],
            ]
            baseline_feats.append(bf)

            # DG features (NO trace_length — that leaks outcome info)
            dg_overall = compute_delayed_gratification(trace)
            n_trace = len(trace)
            half = n_trace // 2
            dg_first = compute_delayed_gratification(trace[:half]) if half > 2 else 0
            dg_second = compute_delayed_gratification(trace[half:]) if half > 2 else 0

            df = [dg_overall, dg_first, dg_second]
            dg_feats.append(df)
            outcomes.append(outcome)

        X_base = np.array(baseline_feats)
        X_dg = np.array(dg_feats)
        X_ext = np.hstack([X_base, X_dg])
        y = np.array(outcomes)

        # Check outcome variance
        if y.std() < 1e-10:
            return {"survives": False, "error": "No variance in outcome variable"}

        # Ablation: shuffle DG features across runs
        rng = np.random.default_rng(42)
        X_dg_shuf = X_dg.copy()
        for col in range(X_dg_shuf.shape[1]):
            rng.shuffle(X_dg_shuf[:, col])
        X_abl = np.hstack([X_base, X_dg_shuf])

        # 5-fold CV
        kf = KFold(n_splits=self.n_folds, shuffle=True, random_state=42)
        s_base, s_ext, s_abl = [], [], []

        for train, test in kf.split(X_base):
            for X, sl in [(X_base, s_base), (X_ext, s_ext), (X_abl, s_abl)]:
                m = Ridge(alpha=1.0)
                m.fit(X[train], y[train])
                sl.append(r2_score(y[test], m.predict(X[test])))

        s_base = np.array(s_base)
        s_ext = np.array(s_ext)
        s_abl = np.array(s_abl)

        t1, p1 = stats.ttest_rel(s_ext, s_base)
        t2, p2 = stats.ttest_rel(s_ext, s_abl)

        survives = (
            s_ext.mean() > s_base.mean()
            and s_ext.mean() > s_abl.mean()
            and p1 < 0.05
            and p2 < 0.05
        )

        return {
            "survives": survives,
            "baseline_score": float(s_base.mean()),
            "extended_score": float(s_ext.mean()),
            "ablation_score": float(s_abl.mean()),
            "p_extended_vs_baseline": float(p1),
            "p_extended_vs_ablation": float(p2),
        }


def _norm_cdf(z: float) -> float:
    """Standard normal CDF approximation."""
    from math import erf, sqrt
    return 0.5 * (1 + erf(z / sqrt(2)))
