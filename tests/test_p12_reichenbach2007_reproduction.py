"""Tests pinning the Sprint 63 L=200 FFT-ring reproduction of Reichenbach 2007 Fig. 2c.

These tests read the pre-computed Sprint 63 output.
Run `analysis/reproductions/p12_reichenbach2007.py` first if the JSON is absent.

Sprint history: Sprint 54 (slope=0.366), Sprint 58 (slope=0.107), Sprint 59
(slope=0.244), Sprint 63 (slope=0.161, L=200 FFT-ring, accepted limitation).
The output JSON at `analysis/outputs/p12_reichenbach2007_reproduction.json` must
carry sprint=63.
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
# B.1 — test_sprint63_json_schema (NOT slow)
# ---------------------------------------------------------------------------
def test_sprint63_json_schema() -> None:
    """Output JSON must have Sprint 63 schema with all required keys."""
    d = _load()
    assert d["sprint"] == 63, f"Expected sprint=63, got {d['sprint']}"

    # Top-level required keys
    required_top = [
        "sprint", "L", "estimator", "m_c", "M_values", "n_fit_points",
        "per_point", "log_log_slope", "r_squared", "tolerance_pass",
        "r_squared_pass", "overall_pass", "n_seeds", "T_eq",
        "sprint_history",
    ]
    for key in required_top:
        assert key in d, f"Missing top-level key: {key}"

    # L=200 and FFT estimator
    assert d["L"] == 200, f"Expected L=200, got {d['L']}"
    assert d["estimator"] == "fft_ring_peak", f"Expected fft_ring_peak, got {d['estimator']}"

    # 5 per_point entries
    assert len(d["per_point"]) == 5, (
        f"Expected 5 per_point entries, got {len(d['per_point'])}"
    )

    # Per-point required keys
    required_pp = ["M", "lambda_mean", "lambda_sem", "n_valid", "n_seeds"]
    for i, p in enumerate(d["per_point"]):
        for key in required_pp:
            assert key in p, f"per_point[{i}] missing key: {key}"

    # n_seeds check
    assert d["per_point"][0]["n_seeds"] == 15, (
        f"Expected n_seeds=15, got {d['per_point'][0]['n_seeds']}"
    )


# ---------------------------------------------------------------------------
# B.2 — test_sprint63_all_seeds_valid (NOT slow)
# ---------------------------------------------------------------------------
def test_sprint63_all_seeds_valid() -> None:
    """All seeds should produce valid wavelength measurements at L=200."""
    d = _load()
    assert d["sprint"] == 63
    for p in d["per_point"]:
        assert p["n_valid"] == d["n_seeds"], (
            f"M={p['M']:.2e}: n_valid={p['n_valid']} < n_seeds={d['n_seeds']}"
        )


# ---------------------------------------------------------------------------
# B.3 — test_sprint63_accepted_limitation (NOT slow)
# ---------------------------------------------------------------------------
def test_sprint63_accepted_limitation() -> None:
    """Sprint 63 is an accepted limitation — slope outside band is expected."""
    d = _load()
    assert d["sprint"] == 63
    # The slope is outside the tolerance band — this is the documented outcome
    assert d["log_log_slope"] is not None, "log_log_slope should not be None"
    assert d["overall_pass"] is False, (
        "Sprint 63 is an accepted limitation; overall_pass should be False"
    )
    # Sprint history should record all four attempts
    assert len(d["sprint_history"]) >= 3, (
        f"Expected ≥3 sprint_history entries, got {len(d['sprint_history'])}"
    )
