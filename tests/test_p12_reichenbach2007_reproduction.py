"""Tests pinning the Sprint 59 near-M_c dense sweep reproduction of Reichenbach 2007 Fig. 2c.

These tests read the pre-computed Sprint 59 output.
Run `analysis/reproductions/p12_reichenbach2007.py` first if the JSON is absent.

Sprint 54 (narrow sweep, slope=0.366) and Sprint 58 (wide sweep, slope=0.107) are
superseded by Sprint 59 (near-M_c dense sweep, 7 points in [2e-4, 5e-4], 30 seeds).
The output JSON at `analysis/outputs/p12_reichenbach2007_reproduction.json` must
carry sprint=59.
"""
import json
import pathlib

import pytest

OUTPUT_PATH = pathlib.Path("analysis/outputs/p12_reichenbach2007_reproduction.json")


def _load() -> dict:
    """Load the reproduction JSON; skip if absent."""
    if not OUTPUT_PATH.exists():
        pytest.skip(f"Output not found: {OUTPUT_PATH}. Run the reproduction script first.")
    with open(OUTPUT_PATH) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# B.1 — test_sprint59_loglog_slope_in_band (slow)
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_sprint59_loglog_slope_in_band() -> None:
    """Sprint 59 log-log slope must be in [0.40, 0.60] with R² ≥ 0.90."""
    d = _load()
    assert d["sprint"] == 59, f"Expected sprint=59, got {d['sprint']}"
    assert d["n_fit_points"] >= 4, f"n_fit_points={d['n_fit_points']} < 4"
    assert d["log_log_slope"] is not None, "log_log_slope is None"
    assert 0.40 <= d["log_log_slope"] <= 0.60, (
        f"log_log_slope={d['log_log_slope']:.4f} outside [0.40, 0.60]"
    )
    assert d["r_squared"] >= 0.90, (
        f"r_squared={d['r_squared']:.4f} below target 0.90"
    )
    assert d["overall_pass"] is True, "overall_pass is not True"


# ---------------------------------------------------------------------------
# B.2 — test_sprint59_near_mc_relative_errors (slow)
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_sprint59_near_mc_relative_errors() -> None:
    """Per-point relative errors < 25% for non-extinction points with n_valid ≥ 15."""
    d = _load()
    assert d["sprint"] == 59, f"Expected sprint=59, got {d['sprint']}"
    qualifying = []
    for p in d["per_point"]:
        if p["near_extinction"] is False and p["n_valid"] >= 15:
            qualifying.append(p)
            assert p["relative_error"] < 0.25, (
                f"M={p['M']:.2e}: relative_error={p['relative_error']:.4f} ≥ 0.25"
            )
    assert len(qualifying) >= 4, (
        f"Only {len(qualifying)} points with near_extinction=False and n_valid≥15; need ≥4"
    )


# ---------------------------------------------------------------------------
# B.3 — test_sprint59_json_schema (NOT slow)
# ---------------------------------------------------------------------------
def test_sprint59_json_schema() -> None:
    """Output JSON must have Sprint 59 schema with all required keys."""
    d = _load()
    assert d["sprint"] == 59, f"Expected sprint=59, got {d['sprint']}"

    # Top-level required keys
    required_top = [
        "sprint", "m_c", "M_values", "fit_M_values", "n_fit_points",
        "per_point", "log_log_slope", "r_squared", "tolerance_pass",
        "r_squared_pass", "n_fit_points_pass", "overall_pass",
    ]
    for key in required_top:
        assert key in d, f"Missing top-level key: {key}"

    # 7 per_point entries
    assert len(d["per_point"]) == 7, (
        f"Expected 7 per_point entries, got {len(d['per_point'])}"
    )

    # Per-point required keys
    required_pp = [
        "M", "lambda_mean", "lambda_sem", "lambda_formula",
        "relative_error", "n_valid", "n_seeds", "near_extinction",
    ]
    for i, p in enumerate(d["per_point"]):
        for key in required_pp:
            assert key in p, f"per_point[{i}] missing key: {key}"

    # n_seeds check
    assert d["per_point"][0]["n_seeds"] == 30, (
        f"Expected n_seeds=30, got {d['per_point'][0]['n_seeds']}"
    )
