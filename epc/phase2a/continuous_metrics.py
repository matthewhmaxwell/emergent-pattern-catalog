"""Per-pattern continuous discriminating scalar for the faithful Phase-2a effect size.

The panel's primary Cohen's d (``cohens_d_positive_vs_panel``) is computed over the
DISCRETE tier-confidence score, which is undefined (NaN, verdict
``TNR-PASS-EFFECT-UNDEFINED``) at perfect separation — a correctly-discriminating
detector gives positives a constant high confidence and rejected negatives a
constant 0, so there is no within-group spread to standardise by.

This registry names, per pattern, the CONTINUOUS canonical discriminating
statistic the detector actually computes (the quantity its screening gate
thresholds and/or feeds to its effect-size computation — the literature's metric
for the phenomenon). The panel then also reports a faithful continuous effect
size = Cohen's d of that scalar, positive seeds vs the pooled negative panel,
ORIENTED so a positive d means the positive separates correctly.

Identified by per-pattern source analysis of each detector's _compute_primary /
screening gate / _compute_effect_size (validation-rebuild, 2026-06-14), with
P8/P14/P19/P21 corrected after verifying the scalar is the one the gate actually
keys on (the automated pass had picked a secondary or the discrete score).

Entry fields:
  kind      "dict"   -> result.primary_metric[key]
            "field"  -> getattr(result, key)
            "nested" -> dotted getattr path (e.g. "fit.lr_vs_exponential"), None-guarded
  direction +1 if a HIGHER value means the pattern is present, -1 if LOWER
  null      value to substitute when the scalar is absent (detector rejected the
            negative before computing it); None -> exclude that item from the d
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


CONTINUOUS_METRIC: Dict[str, Dict[str, Any]] = {
    "P1":  {"kind": "dict",   "key": "morans_i",                   "direction": +1, "null": 0.0},
    "P2":  {"kind": "dict",   "key": "two_phase_score",            "direction": +1, "null": 0.0},
    "P3":  {"kind": "derived","key": "turing_discriminant",        "direction": +1, "null": 0.0},
    "P4":  {"kind": "dict",   "key": "avoidance_ratio",            "direction": -1, "null": 1.0},
    "P5":  {"kind": "dict",   "key": "polarization_mean",          "direction": +1, "null": 0.0},
    "P6":  {"kind": "dict",   "key": "angular_momentum_abs_mean",  "direction": +1, "null": 0.0},
    "P7":  {"kind": "dict",   "key": "lane_order_parameter_mean",  "direction": +1, "null": 0.0},
    "P8":  {"kind": "dict",   "key": "stopped_fraction",           "direction": +1, "null": 0.0},
    "P9":  {"kind": "field",  "key": "r_mean",                     "direction": +1, "null": 0.055},
    "P10": {"kind": "dict",   "key": "pos_vel_ac",                 "direction": +1, "null": 0.0},
    "P11": {"kind": "dict",   "key": "rho_anti",                   "direction": -1, "null": 0.0},
    "P12": {"kind": "dict",   "key": "intransitivity_score",       "direction": +1, "null": 0.0},
    "P13": {"kind": "dict",   "key": "wavefront_speed_cv",         "direction": -1, "null": None},
    "P14": {"kind": "derived","key": "slope_consistency",          "direction": -1, "null": None},
    "P15": {"kind": "dict",   "key": "n_distinct_outcomes",        "direction": +1, "null": 0.0},
    "P16": {"kind": "dict",   "key": "mean_completion_accuracy",   "direction": +1, "null": 0.0},
    "P17": {"kind": "dict",   "key": "chemotactic_index",          "direction": +1, "null": 0.0},
    "P18": {"kind": "dict",   "key": "moran_spearman_early",       "direction": +1, "null": 0.0},
    "P19": {"kind": "dict",   "key": "group_accuracy_mean",        "direction": +1, "null": 0.0},
    "P20": {"kind": "dict",   "key": "step_r2",                    "direction": +1, "null": 0.0},
    "P21": {"kind": "field",  "key": "between_cluster_distance",   "direction": +1, "null": 0.0},
    "P22": {"kind": "dict",   "key": "moran_i_time",               "direction": +1, "null": 0.0},
    "P23": {"kind": "dict",   "key": "scaled_variance",            "direction": -1, "null": 0.25},
    "P24": {"kind": "dict",   "key": "deviation_integral",         "direction": -1, "null": None},
    "P25": {"kind": "dict",   "key": "convergence_variance_ratio", "direction": -1, "null": 1.0},
    "P26": {"kind": "dict",   "key": "peak_performance",           "direction": +1, "null": 0.0},
    "P27": {"kind": "field",  "key": "coop_fraction",              "direction": +1, "null": 0.0},
    "P28": {"kind": "dict",   "key": "gini_final",                 "direction": +1, "null": 0.0},
    "P29": {"kind": "dict",   "key": "weight_dist_corr",           "direction": +1, "null": 0.0},
    "P30": {"kind": "dict",   "key": "radial_cv",                  "direction": -1, "null": 1.0},
    "P32": {"kind": "dict",   "key": "entropy_decline",            "direction": +1, "null": 0.0},
    "P33": {"kind": "dict",   "key": "coherent_half_defect_density","direction": +1, "null": 0.0},
    "P34": {"kind": "dict",   "key": "modularity_gain",            "direction": +1, "null": 0.0},
    "P35": {"kind": "dict",   "key": "psi6_gain",                  "direction": +1, "null": 0.0},
    "P36": {"kind": "dict",   "key": "pareto_advantage",           "direction": +1, "null": 0.0},
    "P37": {"kind": "dict",   "key": "coexistence_oscillation",    "direction": +1, "null": 0.0},
}


# Derived metrics: a few patterns discriminate on a quantity the detector computes
# but does not store as a single primary_metric key. These reproduce the EXACT
# scalar the screening gate keys on.
def _p3_turing_discriminant(result):
    """Binding-gate margin for the 2-D Turing criterion (stationarity x isotropy).

    P3's discrimination is genuinely two-dimensional: a static sinusoid passes
    stationarity but fails isotropy; a spiral passes isotropy but fails
    stationarity. The single faithful scalar is the minimum of the two normalised
    gate margins -- high only when BOTH gates are cleared (a true Turing pattern).
    """
    pm = getattr(result, "primary_metric", None)
    if not isinstance(pm, dict):
        return None
    ae = pm.get("angular_entropy")
    fs = pm.get("field_stationarity")
    if ae is None or fs is None:
        return None
    return min(float(ae) / 0.55, float(fs) / 0.95)


def _p14_slope_consistency(result):
    """|tau_MLE - tau_logbin| -- the SOC slope-consistency the screening gate uses.

    A genuine scale-free law agrees across the two slope estimators (positives
    |delta| <= 0.14); a dissipative-cutoff distribution diverges (>= 0.28).
    """
    fit = getattr(result, "fit", None)
    if fit is None:
        return None
    return abs(float(fit.tau) - float(fit.tau_logbin))


DERIVED_FN = {
    "P3": _p3_turing_discriminant,
    "P14": _p14_slope_consistency,
}
_DERIVED_BY_KEY = {
    "turing_discriminant": _p3_turing_discriminant,
    "slope_consistency": _p14_slope_consistency,
}


def _raw(result: Any, entry: Dict[str, Any]) -> Optional[float]:
    kind = entry["kind"]
    key = entry["key"]
    try:
        if kind == "derived":
            fn = _DERIVED_BY_KEY.get(key)
            return fn(result) if fn is not None else None
        if kind == "dict":
            pm = getattr(result, "primary_metric", None)
            return pm.get(key) if isinstance(pm, dict) else None
        if kind == "field":
            return getattr(result, key, None)
        if kind == "nested":
            cur = result
            for part in key.split("."):
                if cur is None:
                    return None
                cur = getattr(cur, part, None)
            return cur
    except Exception:
        return None
    return None


def extract_continuous(result: Any, entry: Optional[Dict[str, Any]]) -> Optional[float]:
    """Continuous discriminating scalar for one detector result, or None if absent."""
    if entry is None:
        return None
    v = _raw(result, entry)
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if np.isfinite(f) else None


def continuous_cohens_d(
    pos_vals: List[Optional[float]],
    neg_vals: List[Optional[float]],
    entry: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    """Oriented Cohen's d of the continuous canonical metric, positives vs negatives.

    Absent values are filled with the entry's null baseline when one is defined,
    else excluded. d is oriented by ``direction`` so a positive value means the
    canonical positive separates correctly. Returns d=None when fewer than two
    usable values remain in a group, or when both groups are constant (the
    continuous metric is itself saturated — perfect-but-degenerate separation).
    """
    if entry is None:
        return {"d": None, "metric": None, "direction": None, "n_pos": 0, "n_neg": 0, "saturated": False}
    null = entry.get("null")
    direction = entry["direction"]

    def fill(vals: List[Optional[float]]) -> List[float]:
        out: List[float] = []
        for v in vals:
            if v is None:
                if null is not None:
                    out.append(float(null))
            else:
                out.append(float(v))
        return out

    pos = fill(pos_vals)
    neg = fill(neg_vals)
    metric = entry["key"]
    dir_str = "higher" if direction > 0 else "lower"
    base = {"metric": metric, "direction": dir_str, "n_pos": len(pos), "n_neg": len(neg),
            "pos_values": pos, "neg_values": neg}
    if len(pos) < 2 or len(neg) < 2:
        return {**base, "d": None, "saturated": False}
    pa = np.asarray(pos, dtype=float)
    na = np.asarray(neg, dtype=float)
    var_pos = float(np.var(pa, ddof=1))
    var_neg = float(np.var(na, ddof=1))
    pooled = ((pa.size - 1) * var_pos + (na.size - 1) * var_neg) / (pa.size + na.size - 2)
    mean_diff = float(pa.mean()) - float(na.mean())
    if pooled <= 1e-12:
        # The continuous metric is itself saturated (positives constant, negatives
        # constant): separation is perfect but the standardized d is undefined,
        # exactly as for the discrete score. Report None + saturated=True +
        # whether the means are separated in the correct direction.
        return {**base, "d": None, "saturated": True,
                "separated": bool(direction * mean_diff > 0 and abs(mean_diff) > 1e-9)}
    d = direction * mean_diff / float(np.sqrt(pooled))
    return {**base, "d": float(d), "saturated": False}
