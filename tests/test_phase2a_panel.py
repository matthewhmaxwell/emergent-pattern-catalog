"""Tests for the Phase-2a panel harness.

Fast suite: deterministic synthetic generators, hand-checked TNR + Cohen's d
math, catalog loaders smoke (small sizes), JSON output schema.

The actual prototype panel runs (P18, P9) live in
``analysis/run_phase2a_panel.py`` and are out of the fast test scope —
the test here covers the harness, not the prototype results.
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List

import numpy as np
import pytest

from epc.phase2a import PANEL_VERSION
from epc.phase2a import catalog as catalog_mod
from epc.phase2a import synthetic as synth_mod
from epc.phase2a.failed_regimes import p18_voter as p18_failed
from epc.phase2a.failed_regimes import p9_kuramoto as p9_failed
from epc.phase2a.panel import (
    cohens_d,
    compute_tnr,
    overall_verdict,
    run_panel,
)


# --- Spec-compliance metadata ----------------------------------------------

def test_panel_version_constant_is_one_zero():
    assert PANEL_VERSION == "1.0"


def test_synthetic_generators_count_is_ten():
    """Spec Class A: exactly ten synthetic substrates."""
    assert len(synth_mod.SYNTHETIC_GENERATORS) == 10


def test_catalog_ids_count_is_ten():
    """Spec Class B: exactly ten fixed catalog substrates."""
    assert len(catalog_mod.CATALOG_IDS_FIXED) == 10


def test_catalog_self_replacement_uses_fallback():
    """Self-replacement: P9_kuramoto replaced by P12 fallback when D_i tests P9."""
    ids = catalog_mod.catalog_ids_for_pattern("P9")
    assert len(ids) == 10
    assert "P9_kuramoto" not in ids
    assert catalog_mod.FALLBACK_ID in ids


def test_catalog_self_replacement_p18():
    """P18 panel does not include P18_voter; fallback substituted."""
    ids = catalog_mod.catalog_ids_for_pattern("P18")
    assert len(ids) == 10
    assert "P18_voter" not in ids
    assert catalog_mod.FALLBACK_ID in ids


# --- Synthetic generator determinism ---------------------------------------

@pytest.mark.parametrize("name", list(synth_mod.SYNTHETIC_GENERATORS.keys()))
def test_synthetic_generators_deterministic_grid(name):
    """Same seed → byte-identical grid output for every Class A generator."""
    gen = synth_mod.SYNTHETIC_GENERATORS[name]
    # Two generators need a positive trajectory injected; build a tiny one.
    pos = synth_mod.random_binary_field("grid", 0, n_steps=8, shape=(4, 4))
    kwargs: Dict[str, Any] = {"n_steps": 8, "shape": (4, 4)}
    if name in ("permutation_shuffled", "time_shuffled"):
        kwargs = {"positive": pos}
    a = gen("grid", 7, **kwargs)
    b = gen("grid", 7, **kwargs)
    assert len(a) == len(b)
    for sa, sb in zip(a, b):
        assert np.array_equal(sa["grid"], sb["grid"])


def test_synthetic_generator_seeds_differ():
    """Different seed → different output (sanity)."""
    g = synth_mod.random_uniform_field
    a = g("grid", 1, n_steps=4, shape=(8, 8))
    b = g("grid", 2, n_steps=4, shape=(8, 8))
    assert not np.array_equal(a[0]["grid"], b[0]["grid"])


def test_synthetic_phases_format_returns_phases_dict():
    """Phases-format generators emit dicts with `theta` and `r`."""
    h = synth_mod.random_uniform_field("phases", 0, n_steps=5, n=20)
    assert all("theta" in s and "r" in s and "step" in s for s in h)
    assert h[0]["theta"].shape == (20,)


def test_constant_field_grid_is_constant():
    h = synth_mod.constant_field("grid", 0, n_steps=5, shape=(4, 4), value=1)
    assert all(np.all(s["grid"] == 1) for s in h)


def test_periodic_checkerboard_alternates():
    h = synth_mod.periodic_checkerboard("grid", 0, n_steps=2, shape=(4, 4))
    g = h[0]["grid"]
    # Adjacent cells differ along both axes.
    assert (g[0, 0] != g[0, 1]) and (g[0, 0] != g[1, 0])


def test_phases_synthetic_uses_cadence_step():
    """Synthetic phases output has step *= PHASES_DEFAULT_CADENCE so detectors
    that compute n_T_osc from consecutive `step` deltas see the right cadence."""
    h = synth_mod.random_uniform_field("phases", 0, n_steps=5, n=10)
    record_interval = h[1]["step"] - h[0]["step"]
    assert record_interval == synth_mod.PHASES_DEFAULT_CADENCE


# --- TNR + Cohen's d math --------------------------------------------------

def test_tnr_hand_constructed_28_of_30():
    """28 negatives correctly rejected (detected=False), 2 false positives → 28/30 = 0.9333."""
    verdicts = [False] * 28 + [True] * 2
    assert compute_tnr(verdicts) == pytest.approx(28.0 / 30.0, rel=1e-9)


def test_tnr_all_correct():
    assert compute_tnr([False] * 10) == 1.0


def test_tnr_all_wrong():
    assert compute_tnr([True] * 10) == 0.0


def test_cohens_d_hand_computed():
    """Cohen's d for two distributions with known means and equal variances.

    Pos: mean=2.0, var=1.0 (sample); Neg: mean=0.0, var=1.0 (sample).
    Pooled var = 1.0, d = (2 − 0) / 1 = 2.0.
    """
    pos = [1.0, 2.0, 3.0]   # mean=2, sample var = 1.0
    neg = [-1.0, 0.0, 1.0]  # mean=0, sample var = 1.0
    d = cohens_d(pos, neg)
    assert d == pytest.approx(2.0, rel=1e-9)


def test_overall_verdict_pass():
    """≥0.95 overall AND every class ≥0.90 → PASS."""
    assert overall_verdict(0.97, {"a": 0.95, "b": 1.0, "c": 0.95}) == "PASS"


def test_overall_verdict_pass_with_weakness():
    """≥0.95 overall but some class <0.90 → PASS-with-weakness."""
    assert overall_verdict(0.97, {"a": 0.95, "b": 0.80, "c": 1.0}) == "PASS-with-weakness"


def test_overall_verdict_partial():
    """<0.95 overall → PARTIAL."""
    assert overall_verdict(0.93, {"a": 1.0, "b": 1.0, "c": 1.0}) == "PARTIAL"


# --- Catalog loaders -------------------------------------------------------

def test_catalog_loader_p18_voter_grid_format():
    """Catalog loader returns grid-format trajectory for P18 voter substrate."""
    h = catalog_mod.load_catalog_substrate_for_format(
        "P18_voter", "grid", target_steps=20, target_shape=(8, 8),
    )
    assert len(h) == 20
    assert "grid" in h[0]


def test_catalog_loader_p9_kuramoto_phases_format():
    """Catalog loader returns phases-format trajectory for P9 kuramoto substrate."""
    h = catalog_mod.load_catalog_substrate_for_format(
        "P9_kuramoto", "phases", target_steps=20, target_n=50,
    )
    assert len(h) == 20
    assert "theta" in h[0]
    assert h[0]["theta"].shape == (50,)


def test_catalog_loader_unknown_substrate_raises():
    with pytest.raises(KeyError):
        catalog_mod.load_native_substrate("P99_does_not_exist")


# --- Failed-regime registry ------------------------------------------------

def test_p18_failed_regime_config_has_ten_regimes():
    assert len(p18_failed.CONFIG["regimes"]) == 10
    assert p18_failed.CONFIG["format"] == "grid"


def test_p9_failed_regime_config_has_ten_regimes():
    assert len(p9_failed.CONFIG["regimes"]) == 10
    assert p9_failed.CONFIG["format"] == "phases"


def test_p9_failed_regime_k_values_in_subcritical_range():
    """All 10 regimes use K ≤ 0.5 (sub-critical for K_c = 1.0 at γ=0.5)."""
    K_c = 1.0
    K_values = [r["params"]["K"] for r in p9_failed.CONFIG["regimes"]]
    assert all(k <= 0.5 * K_c + 1e-9 for k in K_values)
    assert all(k >= 0.05 * K_c - 1e-9 for k in K_values)


# --- End-to-end harness with a stub detector -------------------------------

class _StubResult:
    def __init__(self, detected: bool, confidence: float):
        self.detected = detected
        self.confidence = confidence
        self.tier = type("Tier", (), {"name": "SCREENING"})()


def _stub_detector_always_reject(history, metadata=None):
    return _StubResult(detected=False, confidence=0.0)


def _stub_detector_always_fire(history, metadata=None):
    return _StubResult(detected=True, confidence=0.9)


class _MinimalFailedRegime:
    """Build a single trivial substrate so the panel ends quickly in tests."""
    CONFIG = {
        "substrate_id": "stub",
        "format": "grid",
        "description": "test stub",
        "regimes": [
            {"label": f"stub_{i}", "params": {}, "seed": i}
            for i in range(10)
        ],
    }

    @staticmethod
    def build_substrate(regime):
        return synth_mod.random_binary_field("grid", regime["seed"], n_steps=10, shape=(4, 4))


def test_panel_writes_json_at_expected_path(tmp_path):
    """Harness writes a JSON file at the requested path with the spec schema."""
    out = str(tmp_path / "stub_panel.json")
    pos = [synth_mod.random_binary_field("grid", s, n_steps=10, shape=(4, 4)) for s in (0, 1)]
    summary = run_panel(
        _stub_detector_always_reject,
        pattern_id="STUB",
        detector_format="grid",
        canonical_positive_runs=pos,
        canonical_metadata={},
        failed_regime_module=_MinimalFailedRegime,
        output_path=out,
        target_steps=10,
        target_shape=(4, 4),
        verbose=False,
    )
    assert os.path.isfile(out)
    with open(out) as f:
        on_disk = json.load(f)
    # Schema check: top-level keys.
    for k in ("pattern_id", "head_commit", "panel_version", "canonical_positive",
              "synthetic", "catalog", "failed_regime", "summary"):
        assert k in on_disk, f"missing key: {k}"
    # Summary keys.
    for k in ("overall_tnr", "synthetic_tnr", "catalog_tnr", "failed_regime_tnr",
              "cohens_d_positive_vs_panel", "verdict"):
        assert k in on_disk["summary"]
    # Always-reject stub: TNR = 1.0 across all classes.
    assert on_disk["summary"]["overall_tnr"] == 1.0
    assert on_disk["summary"]["verdict"] == "PASS"
    # And in-memory mirror matches.
    assert summary["summary"]["overall_tnr"] == 1.0


def test_panel_summary_with_always_fire_stub(tmp_path):
    """Always-fire stub → TNR 0.0, verdict PARTIAL."""
    out = str(tmp_path / "stub_fire_panel.json")
    pos = [synth_mod.random_binary_field("grid", s, n_steps=10, shape=(4, 4)) for s in (0, 1)]
    summary = run_panel(
        _stub_detector_always_fire,
        pattern_id="STUB",
        detector_format="grid",
        canonical_positive_runs=pos,
        canonical_metadata={},
        failed_regime_module=_MinimalFailedRegime,
        output_path=out,
        target_steps=10,
        target_shape=(4, 4),
        verbose=False,
    )
    assert summary["summary"]["overall_tnr"] == 0.0
    assert summary["summary"]["verdict"] == "PARTIAL"
    assert summary["summary"]["n_negatives"] == 30
