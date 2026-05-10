"""Tests for the Phase-2a panel harness (v1.1).

Fast suite: deterministic synthetic generators, hand-checked TNR + Cohen's d
math, v1.1 verdict labels, substrate-typed Class B selection
(``class_b_for_pattern``), Class C N/A handling, catalog loaders smoke,
JSON output schema.

The actual prototype panel runs (P18, P9) live in
``analysis/run_phase2a_panel.py`` and are out of the fast test scope.
"""

from __future__ import annotations

import json
import os
from typing import Any, Dict, List

import numpy as np
import pytest

from epc.phase2a import PANEL_VERSION
from epc.phase2a import catalog as catalog_mod
from epc.phase2a import structured as structured_mod
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

def test_panel_version_constant_is_one_one():
    assert PANEL_VERSION == "1.1"


def test_synthetic_generators_count_is_ten():
    """Spec Class A: exactly ten synthetic substrates."""
    assert len(synth_mod.SYNTHETIC_GENERATORS) == 10


def test_catalog_ids_count_is_ten():
    """Legacy v1.0 fixed Class B list still has 10 entries (kept for reference)."""
    assert len(catalog_mod.CATALOG_IDS_FIXED) == 10


def test_legacy_catalog_self_replacement_p9():
    """v1.0 helper ``catalog_ids_for_pattern`` still works for backward compatibility."""
    ids = catalog_mod.catalog_ids_for_pattern("P9")
    assert len(ids) == 10
    assert "P9_kuramoto" not in ids
    assert catalog_mod.FALLBACK_ID in ids


# --- Synthetic generator determinism ---------------------------------------

@pytest.mark.parametrize("name", list(synth_mod.SYNTHETIC_GENERATORS.keys()))
def test_synthetic_generators_deterministic_grid(name):
    """Same seed → byte-identical grid output for every Class A generator."""
    gen = synth_mod.SYNTHETIC_GENERATORS[name]
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
    assert (g[0, 0] != g[0, 1]) and (g[0, 0] != g[1, 0])


def test_phases_synthetic_uses_cadence_step():
    """Synthetic phases output has step *= PHASES_DEFAULT_CADENCE."""
    h = synth_mod.random_uniform_field("phases", 0, n_steps=5, n=10)
    record_interval = h[1]["step"] - h[0]["step"]
    assert record_interval == synth_mod.PHASES_DEFAULT_CADENCE


# --- TNR + Cohen's d math --------------------------------------------------

def test_tnr_hand_constructed_28_of_30():
    """28 negatives correctly rejected, 2 false positives → 28/30 = 0.9333."""
    verdicts = [False] * 28 + [True] * 2
    assert compute_tnr(verdicts) == pytest.approx(28.0 / 30.0, rel=1e-9)


def test_tnr_all_correct():
    assert compute_tnr([False] * 10) == 1.0


def test_tnr_all_wrong():
    assert compute_tnr([True] * 10) == 0.0


def test_cohens_d_hand_computed():
    """Pos: mean=2.0, var=1.0; Neg: mean=0.0, var=1.0. Pooled var=1.0, d=2.0."""
    pos = [1.0, 2.0, 3.0]
    neg = [-1.0, 0.0, 1.0]
    assert cohens_d(pos, neg) == pytest.approx(2.0, rel=1e-9)


# --- v1.1 verdict labels ---------------------------------------------------

def test_v1_1_verdict_pass():
    """≥0.95 overall AND every gating class ≥0.90 AND d ≥1.0 → PASS."""
    assert overall_verdict(0.97, {"a": 0.95, "b": 1.0}, cohens_d_value=2.0) == "PASS"


def test_v1_1_verdict_pass_with_weakness():
    """≥0.95 overall AND d ≥1.0 but a gating class <0.90 → PASS-with-weakness."""
    assert overall_verdict(0.97, {"a": 0.95, "b": 0.80}, cohens_d_value=2.0) == "PASS-with-weakness"


def test_v1_1_verdict_partial():
    """<0.95 overall AND d ≥0.5 → PARTIAL."""
    assert overall_verdict(0.80, {"a": 0.85, "b": 0.85}, cohens_d_value=1.5) == "PARTIAL"


def test_v1_1_verdict_fail():
    """<0.95 overall AND d <0.5 → FAIL."""
    assert overall_verdict(0.50, {"a": 0.40, "b": 0.40}, cohens_d_value=0.2) == "FAIL"


def test_v1_1_verdict_partial_when_tnr_passes_but_d_fails_pass_threshold():
    """TNR clears 0.95 but d <1.0 (still ≥0.5) → PARTIAL (criterion 4 not met)."""
    assert overall_verdict(0.97, {"a": 0.95, "b": 1.0}, cohens_d_value=0.7) == "PARTIAL"


# --- Substrate-typed Class B (v1.1) ----------------------------------------

def test_class_b_p9_oscillator_one_mate_two_supplements():
    """P9 (oscillator): one catalog-mate (P10_chimera) + two B' supplements."""
    r = catalog_mod.class_b_for_pattern("P9")
    assert r["substrate_type"] == "oscillator"
    assert r["catalog_mates"] == ["P10_chimera"]
    assert len(r["synthetic_supplements"]) >= 2
    assert "incoherent_phases" in r["synthetic_supplements"]
    assert "subcritical_kuramoto" in r["synthetic_supplements"]


def test_class_b_p18_network_override():
    """P18 reclassified to network per v1.1 spec; mate is P21_hegselmann_krause + 2 supps."""
    r = catalog_mod.class_b_for_pattern("P18")
    assert r["substrate_type"] == "network"
    assert r["catalog_mates"] == ["P21_hegselmann_krause"]
    assert len(r["synthetic_supplements"]) >= 2


def test_class_b_p21_network_override():
    """P21 reclassified to network per v1.1 spec; mate is P18_voter + 2 supps."""
    r = catalog_mod.class_b_for_pattern("P21")
    assert r["substrate_type"] == "network"
    assert r["catalog_mates"] == ["P18_voter"]


def test_class_b_p1_lattice_2d_seven_or_more_mates_no_supplements():
    """P1 (lattice_2d): ≥7 catalog-mates from same substrate type, no supplements."""
    r = catalog_mod.class_b_for_pattern("P1")
    assert r["substrate_type"] == "lattice_2d"
    assert len(r["catalog_mates"]) >= 7
    assert "P1_schelling" not in r["catalog_mates"]
    assert r["synthetic_supplements"] == []


def test_class_b_unknown_pattern_returns_empty():
    """Unknown pattern_id → empty composition, no error."""
    r = catalog_mod.class_b_for_pattern("PX_unknown")
    assert r["substrate_type"] is None
    assert r["catalog_mates"] == []
    assert r["synthetic_supplements"] == []


def test_substrate_type_overrides_p18_and_p21():
    """Sanity: the v1.1 reclassification for voter and HK is in effect."""
    assert catalog_mod.SUBSTRATE_TYPE_BY_PATTERN["P18"] == "network"
    assert catalog_mod.SUBSTRATE_TYPE_BY_PATTERN["P21"] == "network"


def test_supplement_builders_deterministic():
    """Class B' supplements are deterministic functions of seed."""
    a = structured_mod.incoherent_phases(0, n=10, n_steps=4, cadence=1)
    b = structured_mod.incoherent_phases(0, n=10, n_steps=4, cadence=1)
    assert all(np.array_equal(sa["theta"], sb["theta"]) for sa, sb in zip(a, b))


def test_random_graph_evolution_is_grid_format():
    """Network supplement returns grid-format substrate ready for grid detectors."""
    h = structured_mod.random_graph_evolution(0, rows=4, cols=4, n_steps=3)
    assert all("grid" in s for s in h)
    assert h[0]["grid"].shape == (4, 4)


# --- Catalog loaders -------------------------------------------------------

def test_catalog_loader_p18_voter_grid_format():
    h = catalog_mod.load_catalog_substrate_for_format(
        "P18_voter", "grid", target_steps=20, target_shape=(8, 8),
    )
    assert len(h) == 20
    assert "grid" in h[0]


def test_catalog_loader_p9_kuramoto_phases_format():
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

def test_p18_failed_regime_is_class_c_n_a():
    """v1.1: P18 declares Class C N/A — voter has no parameter that suppresses consensus."""
    assert p18_failed.CONFIG.get("status") == "N/A"
    assert "n_a_reason" in p18_failed.CONFIG
    assert p18_failed.CONFIG["regimes"] == []


def test_p18_build_substrate_raises_when_called():
    """Defensive: build_substrate must not be reachable when Class C is N/A."""
    with pytest.raises(RuntimeError):
        p18_failed.build_substrate({"label": "x", "params": {}, "seed": 0})


def test_p9_failed_regime_config_has_ten_regimes():
    """v1.1 Class C populated for P9 (sub-K_c regimes)."""
    assert p9_failed.CONFIG.get("status") != "N/A"
    assert len(p9_failed.CONFIG["regimes"]) == 10
    assert p9_failed.CONFIG["format"] == "phases"


def test_p9_failed_regime_k_values_in_subcritical_range():
    K_c = 1.0
    K_values = [r["params"]["K"] for r in p9_failed.CONFIG["regimes"]]
    assert all(k <= 0.5 * K_c + 1e-9 for k in K_values)
    assert all(k >= 0.05 * K_c - 1e-9 for k in K_values)


# --- End-to-end harness with stub detectors --------------------------------

class _StubResult:
    def __init__(self, detected: bool, confidence: float):
        self.detected = detected
        self.confidence = confidence
        self.tier = type("Tier", (), {"name": "SCREENING"})()


class _DiscriminatingStub:
    """Stub that returns positive verdict for the first ``n_positive`` calls.

    The panel runner calls the detector first against the canonical positive
    seeds, then against the negative substrates. With ``n_positive`` matching
    the positive-seed count we get a clean separation: positives always fire,
    negatives never fire — Cohen's d → +∞ and verdict → PASS.
    """

    def __init__(self, n_positive: int):
        self.n_positive = n_positive
        self.calls = 0

    def __call__(self, history, metadata=None):
        idx = self.calls
        self.calls += 1
        if idx < self.n_positive:
            return _StubResult(detected=True, confidence=0.9)
        return _StubResult(detected=False, confidence=0.0)


def _stub_detector_always_fire(history, metadata=None):
    return _StubResult(detected=True, confidence=0.9)


class _MinimalFailedRegime:
    """Minimal Class C config with 10 trivial regimes for end-to-end smoke."""
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


class _NaFailedRegime:
    """Minimal Class C config with N/A status."""
    CONFIG = {
        "substrate_id": "stub_na",
        "format": "grid",
        "status": "N/A",
        "n_a_reason": "test reason",
        "regimes": [],
    }

    @staticmethod
    def build_substrate(regime):
        raise RuntimeError("must not be called")


def test_panel_writes_json_with_v1_1_schema(tmp_path):
    """Harness writes JSON containing the v1.1 schema fields."""
    out = str(tmp_path / "stub_panel.json")
    pos = [synth_mod.random_binary_field("grid", s, n_steps=10, shape=(4, 4)) for s in (0, 1)]
    summary = run_panel(
        _DiscriminatingStub(n_positive=2),
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
    # v1.1 schema additions.
    assert on_disk["panel_version"] == "1.1"
    assert "class_b_composition" in on_disk
    assert set(on_disk["class_b_composition"].keys()) >= {"catalog_mates", "synthetic_supplements", "substrate_type"}
    assert "class_c_status" in on_disk
    assert on_disk["class_c_status"] == "populated"
    # Discriminating stub: TNR=1.0 (negs all rejected), positives all fire → d very large → PASS.
    assert on_disk["summary"]["overall_tnr"] == 1.0
    assert on_disk["summary"]["verdict"] == "PASS"
    # Class B is empty for unknown STUB pattern (substrate_type is None) → n_negatives = 10 (Class A) + 0 (Class B) + 10 (Class C) = 20.
    assert on_disk["summary"]["n_negatives"] == 20


def test_panel_class_c_n_a_skipped_in_summary(tmp_path):
    """When failed_regime_module declares status='N/A', Class C is skipped."""
    out = str(tmp_path / "stub_na_panel.json")
    pos = [synth_mod.random_binary_field("grid", s, n_steps=10, shape=(4, 4)) for s in (0, 1)]
    summary = run_panel(
        _DiscriminatingStub(n_positive=2),
        pattern_id="STUB",
        detector_format="grid",
        canonical_positive_runs=pos,
        canonical_metadata={},
        failed_regime_module=_NaFailedRegime,
        output_path=out,
        target_steps=10,
        target_shape=(4, 4),
        verbose=False,
    )
    assert summary["class_c_status"] == "N/A"
    assert summary["class_c_n_a_reason"] == "test reason"
    assert summary["failed_regime"] == []
    # n_negatives = Class A (10) + Class B (0 for STUB) + Class C N/A (0) = 10.
    assert summary["summary"]["n_negatives"] == 10
    assert summary["summary"]["failed_regime_tnr"] is None


def test_panel_always_fire_stub_yields_fail_verdict(tmp_path):
    """Always-fire stub: TNR=0, no signal → FAIL under v1.1 verdict labels."""
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
    # All inputs (positive + negative) score the same constant → no discriminating signal → FAIL.
    assert summary["summary"]["verdict"] == "FAIL"
