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

def test_panel_version_constant_is_current():
    """Sprint 34: PANEL_VERSION is the active spec version."""
    # Updated each spec-revision sprint. test_panel_version_constant_is_one_two
    # added in Sprint 34 has the explicit "1.2" check.
    assert PANEL_VERSION in ("1.2",)


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


def test_p21_generator_deterministic():
    """Same seed → byte-identical opinion trajectory (closes C-p21-generator)."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS["P21_hegselmann_krause"])
    a = catalog_mod._gen_p21_hegselmann_krause(p)
    b = catalog_mod._gen_p21_hegselmann_krause(p)
    assert a["kind"] == "opinions"
    assert np.array_equal(a["opinions"], b["opinions"])


def test_p21_canonical_positive_is_fragmented():
    """ε=0.2 should produce a multimodal final opinion distribution (not consensus)."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS["P21_hegselmann_krause"])
    native = catalog_mod._gen_p21_hegselmann_krause(p)
    final = native["opinions"][-1]
    sorted_op = np.sort(final)
    n_clusters = int(np.sum(np.diff(sorted_op) > 0.1)) + 1
    assert n_clusters >= 2, f"expected ≥2 clusters at ε=0.2, got {n_clusters}"
    # And not trivial single-spike consensus.
    assert (final.max() - final.min()) > 0.2, f"opinion spread too small: {final.max() - final.min()}"


def test_p21_native_kind_is_opinions():
    """P21 generator's output 'kind' field marks it as opinion-vector substrate."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS["P21_hegselmann_krause"])
    native = catalog_mod._gen_p21_hegselmann_krause(p)
    assert native["kind"] == "opinions"
    assert "n_agents" in native and native["n_agents"] == p["n_agents"]
    assert "epsilon" in native and native["epsilon"] == p["epsilon"]


def test_p18_class_b_now_includes_p21_generator():
    """After Sprint 32 P21 generator landed, P18's network catalog-mate is loadable."""
    r = catalog_mod.class_b_for_pattern("P18")
    assert "P21_hegselmann_krause" in r["catalog_mates"]
    # And the substrate id is now backed by a real generator.
    assert "P21_hegselmann_krause" in catalog_mod._GENERATORS


# Sprint 33 — three new lattice_2d generators (P11, P13, P22) needed for P14/P15 Class B.

@pytest.mark.parametrize("sid", ["P11_lotka_volterra", "P13_greenberg_hastings", "P22_sir_epidemic"])
def test_lattice_2d_generator_registered(sid):
    """Sprint 33: each new lattice_2d catalog-mate has a generator in _GENERATORS."""
    assert sid in catalog_mod._GENERATORS
    assert sid in catalog_mod.SUBSTRATE_PARAMS


@pytest.mark.parametrize("sid", ["P11_lotka_volterra", "P13_greenberg_hastings", "P22_sir_epidemic"])
def test_lattice_2d_generator_returns_grid_categorical(sid):
    """Each new generator returns kind=grid_categorical with stacked int8 grids."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS[sid])
    native = catalog_mod._GENERATORS[sid](p)
    assert native["kind"] == "grid_categorical"
    assert native["grids"].dtype == np.int8
    assert native["grids"].ndim == 3  # (T, rows, cols)


def test_p14_class_b_fully_loadable_after_sprint_33():
    """All P14 catalog-mates now have generators (no missing)."""
    r = catalog_mod.class_b_for_pattern("P14")
    missing = [m for m in r["catalog_mates"] if m not in catalog_mod._GENERATORS]
    assert missing == [], f"P14 catalog-mates without generators: {missing}"


def test_p15_class_b_fully_loadable_after_sprint_33():
    """All P15 catalog-mates now have generators (no missing)."""
    r = catalog_mod.class_b_for_pattern("P15")
    missing = [m for m in r["catalog_mates"] if m not in catalog_mod._GENERATORS]
    assert missing == [], f"P15 catalog-mates without generators: {missing}"


# Sprint 34 — primary-metric invariance flags (Phase-2a v1.2 §Change 1+2)

from epc.phase2a.detector_invariance import (
    DETECTOR_INVARIANCE_FLAGS,
    get_flags,
    is_permutation_invariant,
    is_time_shuffle_invariant,
)
from epc.phase2a.panel import (
    SKIP_VERDICT,
    PERM_SHUFFLED_SUBSTRATE,
    TIME_SHUFFLED_SUBSTRATE,
)


def test_panel_version_constant_is_one_two():
    """Sprint 34: PANEL_VERSION bumped to 1.2."""
    from epc.phase2a import PANEL_VERSION
    assert PANEL_VERSION == "1.2"


@pytest.mark.parametrize("pid, perm, time", [
    ("P9",  True,  True),
    ("P1",  False, False),
    ("P14", True,  True),
    ("P15", False, False),
    ("P5",  True,  False),
    ("P11", False, True),
    ("P18", True,  True),
    ("P21", True,  True),
    ("P28", True,  True),
    ("P27", False, True),   # provisional per spec Change 2
])
def test_invariance_flags_match_v1_2_spec_table(pid, perm, time):
    """Spec v1.2 §Change 2 authoritative flag values per pattern."""
    f = get_flags(pid)
    assert f.permutation_invariant is perm, f"{pid} perm_inv: spec={perm} got={f.permutation_invariant}"
    assert f.time_shuffle_invariant is time, f"{pid} time_inv: spec={time} got={f.time_shuffle_invariant}"


def test_invariance_flags_unknown_pattern_returns_safe_default():
    """Unknown pattern → both False (conservative; runs every substrate)."""
    f = get_flags("PX_unknown")
    assert f.permutation_invariant is False
    assert f.time_shuffle_invariant is False


def test_invariance_lookup_helpers():
    """Module exports get_flags + helper booleans."""
    assert is_permutation_invariant("P9") is True
    assert is_time_shuffle_invariant("P9") is True
    assert is_permutation_invariant("P1") is False


# --- Harness skip-logic + JSON schema (v1.2) -------------------------------


class _CountingStub:
    """Discriminating stub that ALSO counts how many positive vs negative calls it gets."""
    def __init__(self, n_positive: int):
        self.n_positive = n_positive
        self.calls = 0

    def __call__(self, history, metadata=None):
        idx = self.calls
        self.calls += 1
        if idx < self.n_positive:
            return _StubResult(detected=True, confidence=0.9)
        return _StubResult(detected=False, confidence=0.0)


def _stub_failed_module_minimal():
    return _MinimalFailedRegime


def test_panel_skips_permutation_shuffled_for_perm_invariant_detector(tmp_path):
    """v1.2 §Change 1: detector with perm_inv=True → permutation_shuffled is skipped."""
    out = str(tmp_path / "perm_inv_panel.json")
    pos = [synth_mod.random_binary_field("grid", s, n_steps=10, shape=(4, 4)) for s in (0, 1)]
    summary = run_panel(
        _CountingStub(n_positive=2),
        pattern_id="P9",   # perm_inv=True, time_inv=True
        detector_format="grid",
        canonical_positive_runs=pos,
        canonical_metadata={},
        failed_regime_module=_MinimalFailedRegime,
        output_path=out,
        target_steps=10,
        target_shape=(4, 4),
        verbose=False,
    )
    syn = summary["synthetic"]
    # Both permutation_shuffled and time_shuffled should be skipped for P9.
    skipped = [s for s in syn if s["verdict"] == SKIP_VERDICT]
    assert len(skipped) == 2, f"expected 2 skips for P9, got {len(skipped)}"
    skipped_ids = {s["substrate"] for s in skipped}
    assert skipped_ids == {PERM_SHUFFLED_SUBSTRATE, TIME_SHUFFLED_SUBSTRATE}
    # Skip reasons must reference the correct flag.
    for s in skipped:
        if s["substrate"] == PERM_SHUFFLED_SUBSTRATE:
            assert s["skip_reason"] == "primary_metric_permutation_invariant"
        else:
            assert s["skip_reason"] == "primary_metric_time_shuffle_invariant"


def test_panel_runs_all_class_a_for_non_invariant_detector(tmp_path):
    """v1.2 §Change 1: detector with both flags False → all 10 Class A substrates run."""
    out = str(tmp_path / "non_inv_panel.json")
    pos = [synth_mod.random_binary_field("grid", s, n_steps=10, shape=(4, 4)) for s in (0, 1)]
    summary = run_panel(
        _CountingStub(n_positive=2),
        pattern_id="P1",   # perm_inv=False, time_inv=False
        detector_format="grid",
        canonical_positive_runs=pos,
        canonical_metadata={},
        failed_regime_module=_MinimalFailedRegime,
        output_path=out,
        target_steps=10,
        target_shape=(4, 4),
        verbose=False,
    )
    skipped = [s for s in summary["synthetic"] if s["verdict"] == SKIP_VERDICT]
    assert skipped == []
    assert summary["summary"]["class_a_size_total"] == 10
    assert summary["summary"]["class_a_size_evaluated"] == 10


def test_panel_partial_skip_for_perm_only_invariant(monkeypatch, tmp_path):
    """perm_inv=True + time_inv=False → only permutation_shuffled is skipped.

    Uses a temporary stub pattern (no implemented catalog mates) so the panel
    runner doesn't try to load unimplemented continuous_2d generators when
    exercising P5's class B.
    """
    from epc.phase2a.detector_invariance import (
        DETECTOR_INVARIANCE_FLAGS, InvarianceFlags,
    )
    monkeypatch.setitem(
        DETECTOR_INVARIANCE_FLAGS,
        "STUB_PERM_ONLY",
        InvarianceFlags(True, False, "stub", rationale="test only"),
    )

    out = str(tmp_path / "perm_only_panel.json")
    pos = [synth_mod.random_binary_field("grid", s, n_steps=10, shape=(4, 4)) for s in (0, 1)]
    summary = run_panel(
        _CountingStub(n_positive=2),
        pattern_id="STUB_PERM_ONLY",
        detector_format="grid",
        canonical_positive_runs=pos,
        canonical_metadata={},
        failed_regime_module=_MinimalFailedRegime,
        output_path=out,
        target_steps=10,
        target_shape=(4, 4),
        verbose=False,
    )
    skipped = [s for s in summary["synthetic"] if s["verdict"] == SKIP_VERDICT]
    assert len(skipped) == 1
    assert skipped[0]["substrate"] == PERM_SHUFFLED_SUBSTRATE
    assert summary["summary"]["class_a_size_evaluated"] == 9


def test_panel_tnr_uses_evaluated_denominator_not_total():
    """v1.2 §Harness output: TNR is computed over evaluated Class A only.

    Hand-constructed: 8 evaluated (2 skipped), 7 correctly rejected, 1 false
    positive → TNR_syn = 7/8 = 0.875, NOT 7/10 = 0.7.
    """
    # Directly compute compute_tnr on the 8-element evaluated list.
    evaluated = [False] * 7 + [True] * 1
    tnr = compute_tnr(evaluated)
    assert tnr == pytest.approx(7.0 / 8.0, rel=1e-9)


def test_panel_json_includes_class_a_size_fields(tmp_path):
    """Spec v1.2 §"Harness output": summary block includes class_a_size_{total,evaluated}."""
    out = str(tmp_path / "schema_panel.json")
    pos = [synth_mod.random_binary_field("grid", s, n_steps=10, shape=(4, 4)) for s in (0, 1)]
    run_panel(
        _CountingStub(n_positive=2),
        pattern_id="P14",   # both flags True → 2 skipped
        detector_format="grid",
        canonical_positive_runs=pos,
        canonical_metadata={},
        failed_regime_module=_MinimalFailedRegime,
        output_path=out,
        target_steps=10,
        target_shape=(4, 4),
        verbose=False,
    )
    with open(out) as f:
        on_disk = json.load(f)
    assert on_disk["panel_version"] == "1.2"
    assert on_disk["summary"]["class_a_size_total"] == 10
    assert on_disk["summary"]["class_a_size_evaluated"] == 8
    # Schema: detector_invariance metadata at top level.
    assert "detector_invariance" in on_disk
    assert on_disk["detector_invariance"]["permutation_invariant"] is True
    assert on_disk["detector_invariance"]["time_shuffle_invariant"] is True
    # Skipped substrates have the canonical skip_reason format.
    skipped = [s for s in on_disk["synthetic"] if s["verdict"] == SKIP_VERDICT]
    assert all("skip_reason" in s for s in skipped)


def test_p27_time_shuffle_flag_is_provisional_per_spec():
    """P27's time_shuffle_invariant=True is provisional per v1.2 §Change 2 + carry-forward
    C-p27-time-shuffle-invariance. The flag value is True per the spec table; the
    *test* exists to make the provisional status explicit in code.
    """
    f = get_flags("P27")
    assert f.permutation_invariant is False
    assert f.time_shuffle_invariant is True
    assert "Provisional" in f.rationale or "C-p27-time-shuffle-invariance" in f.rationale


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
    # v1.1+ schema fields (still present in v1.2).
    assert on_disk["panel_version"] == "1.2"
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


# --- Sprint 37: continuous_2d catalog generators (P2 ABP + P6 D'Orsogna) ---

def test_gen_p2_active_brownian_deterministic():
    """Same seed → byte-identical ABP trajectory."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS["P2_abp"])
    a = catalog_mod._gen_p2_active_brownian(p)
    b = catalog_mod._gen_p2_active_brownian(p)
    assert a["kind"] == "particles"
    assert np.array_equal(a["headings"], b["headings"])
    assert np.array_equal(a["positions"], b["positions"])


def test_gen_p6_dorsogna_deterministic():
    """Same seed → byte-identical D'Orsogna trajectory."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS["P6_dorsogna"])
    a = catalog_mod._gen_p6_dorsogna(p)
    b = catalog_mod._gen_p6_dorsogna(p)
    assert a["kind"] == "particles"
    assert np.array_equal(a["headings"], b["headings"])
    assert np.array_equal(a["positions"], b["positions"])


def test_gen_p2_active_brownian_output_format():
    """P2 generator returns kind='particles' schema with correct shapes."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS["P2_abp"])
    native = catalog_mod._gen_p2_active_brownian(p)
    assert native["kind"] == "particles"
    assert native["headings"].ndim == 2           # (T, N)
    assert native["positions"].ndim == 3          # (T, N, 2)
    assert native["headings"].shape[0] == native["positions"].shape[0]
    assert native["positions"].shape[2] == 2
    assert "box_size" in native and native["box_size"] > 0


def test_gen_p6_dorsogna_output_format():
    """P6 generator returns kind='particles' schema with correct shapes."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS["P6_dorsogna"])
    native = catalog_mod._gen_p6_dorsogna(p)
    assert native["kind"] == "particles"
    assert native["headings"].ndim == 2           # (T, N)
    assert native["positions"].ndim == 3          # (T, N, 2)
    assert native["headings"].shape[0] == native["positions"].shape[0]
    assert native["positions"].shape[2] == 2
    assert "box_size" in native and native["box_size"] > 0


def test_class_b_p5_contains_p2_abp_and_p6_dorsogna():
    """Sprint 37: class_b_for_pattern('P5') includes both new continuous_2d mates."""
    r = catalog_mod.class_b_for_pattern("P5")
    assert r["substrate_type"] == "continuous_2d"
    assert "P2_abp" in r["catalog_mates"]
    assert "P6_dorsogna" in r["catalog_mates"]


def test_continuous_2d_supplements_registered():
    """Sprint 37: SUPPLEMENTS_BY_SUBSTRATE_TYPE has both new continuous_2d builders."""
    sup = structured_mod.SUPPLEMENTS_BY_SUBSTRATE_TYPE.get("continuous_2d", [])
    assert "uncorrelated_random_walks" in sup
    assert "independent_brownian_motion" in sup
    # Builders callable and in SUPPLEMENT_BUILDERS
    assert "uncorrelated_random_walks" in structured_mod.SUPPLEMENT_BUILDERS
    assert "independent_brownian_motion" in structured_mod.SUPPLEMENT_BUILDERS


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


# --- Sprint 38: P8 NS generator + lattice_1d supplements -------------------

def test_gen_p8_nagel_schreckenberg_deterministic():
    """Sprint 38: same seed → byte-identical NS occupancy trajectory."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS["P8_nagel_schreckenberg"])
    a = catalog_mod._gen_p8_nagel_schreckenberg(p)
    b = catalog_mod._gen_p8_nagel_schreckenberg(p)
    assert a["kind"] == "sequence"
    assert np.array_equal(a["arrays"], b["arrays"])


def test_gen_p8_nagel_schreckenberg_output_shape():
    """Sprint 38: P8 generator returns kind='sequence' with shape (T, road_length)."""
    p = dict(catalog_mod.SUBSTRATE_PARAMS["P8_nagel_schreckenberg"])
    native = catalog_mod._gen_p8_nagel_schreckenberg(p)
    assert native["kind"] == "sequence"
    arrays = native["arrays"]
    assert arrays.ndim == 2
    assert arrays.shape[1] == p["road_length"]
    assert arrays.shape[0] > 0


def test_class_b_p31_includes_p8_nagel_schreckenberg():
    """Sprint 38: class_b_for_pattern('P31') catalog_mates includes P8_nagel_schreckenberg."""
    r = catalog_mod.class_b_for_pattern("P31")
    assert r["substrate_type"] == "lattice_1d"
    assert "P8_nagel_schreckenberg" in r["catalog_mates"]


def test_lattice_1d_supplements_registered():
    """Sprint 38: SUPPLEMENTS_BY_SUBSTRATE_TYPE['lattice_1d'] has both new builders."""
    sup = structured_mod.SUPPLEMENTS_BY_SUBSTRATE_TYPE.get("lattice_1d", [])
    assert "independent_lane_traffic" in sup
    assert "reverse_sorted_sequence" in sup
    assert "independent_lane_traffic" in structured_mod.SUPPLEMENT_BUILDERS
    assert "reverse_sorted_sequence" in structured_mod.SUPPLEMENT_BUILDERS
