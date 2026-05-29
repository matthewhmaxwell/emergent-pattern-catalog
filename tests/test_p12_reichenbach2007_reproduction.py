"""Tests pinning the Sprint 58 wide-sweep reproduction of Reichenbach 2007 Fig. 2c.

These tests read the pre-computed Sprint 58 output.
Run `analysis/reproductions/p12_reichenbach2007.py` first if the JSON is absent.

Sprint 54 (narrow sweep, slope=0.366) is superseded by Sprint 58 (wide sweep,
7 points, 15 seeds each). The output JSON at
`analysis/outputs/p12_reichenbach2007_reproduction.json` must carry sprint=58.
"""
import json
import os
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
# B.1 — Sprint 58 log-log slope and R² in acceptance band
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_sprint58_loglog_slope_in_band() -> None:
    """Slope must lie in [0.40, 0.60] and R² ≥ 0.90 for dim1 closure."""
    result = _load()
    assert result["sprint"] == 58, f"Expected sprint=58, got {result['sprint']}"
    slope = result["log_log_slope"]
    r2 = result["r_squared"]
    assert 0.40 <= slope <= 0.60, (
        f"log_log_slope={slope:.4f} outside acceptance band [0.40, 0.60]. "
        f"Sprint 58 wide-sweep ({result.get('note', '')})."
    )
    assert r2 >= 0.90, f"R²={r2:.4f} below target 0.90."
    assert result["overall_pass"] is True, (
        f"overall_pass=False (slope={slope:.4f}, R²={r2:.4f})."
    )


# ---------------------------------------------------------------------------
# B.2 — Per-point relative errors and valid-seed counts
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_sprint58_per_point_relative_error() -> None:
    """Each of the 7 M points must have rel_error < 0.25 and n_valid >= 10."""
    result = _load()
    assert result["sprint"] == 58
    per_point = result["per_point"]
    assert len(per_point) == 7, f"Expected 7 per-point entries, got {len(per_point)}"
    failures = []
    for p in per_point:
        M = p["M"]
        rel_err = p.get("relative_error")
        n_valid = p["n_valid"]
        if rel_err is None or rel_err >= 0.25:
            failures.append(
                f"M={M:.2e}: relative_error={rel_err} (must be < 0.25)"
            )
        if n_valid < 10:
            failures.append(
                f"M={M:.2e}: n_valid={n_valid} < 10 (minimum valid-seed requirement)"
            )
    assert not failures, "Per-point failures:\n" + "\n".join(failures)


# ---------------------------------------------------------------------------
# B.3 — Sprint 54 output superseded (non-slow: only reads JSON)
# ---------------------------------------------------------------------------
def test_sprint54_output_superseded() -> None:
    """Output JSON must carry sprint=58, confirming Sprint 54 result is superseded."""
    result = _load()
    assert result["sprint"] == 58, (
        f"Output JSON has sprint={result['sprint']}, expected 58. "
        "Sprint 54 narrow-sweep result has not been overwritten."
    )
    assert "per_point" in result, "Missing 'per_point' key in output JSON."
    assert len(result["per_point"]) == 7, (
        f"Expected 7 per-point entries (Sprint 58 wide sweep), got {len(result['per_point'])}."
    )
    assert "log_log_slope" in result, "Missing 'log_log_slope' key (Sprint 54 used 'fit_slope')."
    assert "r_squared" in result, "Missing 'r_squared' key."
    assert "overall_pass" in result, "Missing 'overall_pass' key."
