"""Orchestration layer tests — substrate-aware detector dispatch.

Tests the compatibility matrix: 5 substrate types, 11 models × 10 detectors,
24 compatible pairs. Verifies substrate filtering, observable guards, canonical
positive/negative assignments, and that Sprint 5 models (Nowak-May, HK) are
correctly registered.

Architecture decision #25 (updated Sprint 5 verification).
"""

import pytest


# ---------------------------------------------------------------------------
# Import orchestration (adjust if module path differs)
# ---------------------------------------------------------------------------
from epc.orchestration import (
    build_default_registry,
    Substrate,
)


@pytest.fixture(scope="module")
def registry():
    """Build the default registry once for all tests."""
    return build_default_registry()


# ---------------------------------------------------------------------------
# Substrate type counts
# ---------------------------------------------------------------------------

class TestSubstrateCounts:
    """Verify the 5 substrate types are represented."""

    def test_lattice_1d_models_exist(self, registry):
        models = [m for m in registry.models.values()
                  if m.substrate == Substrate.LATTICE_1D]
        assert len(models) >= 1, "Need at least 1 lattice_1d model (sorting)"

    def test_lattice_2d_models_exist(self, registry):
        models = [m for m in registry.models.values()
                  if m.substrate == Substrate.LATTICE_2D]
        # GH, GoL, Schelling, Nowak-May
        assert len(models) >= 4, f"Expected ≥4 lattice_2d models, got {len(models)}"

    def test_continuous_2d_models_exist(self, registry):
        models = [m for m in registry.models.values()
                  if m.substrate == Substrate.CONTINUOUS_2D]
        # Vicsek, D'Orsogna
        assert len(models) >= 2, f"Expected ≥2 continuous_2d models, got {len(models)}"

    def test_oscillator_models_exist(self, registry):
        models = [m for m in registry.models.values()
                  if m.substrate == Substrate.OSCILLATOR]
        assert len(models) >= 1, "Need at least 1 oscillator model (Kuramoto)"

    def test_opinion_space_models_exist(self, registry):
        models = [m for m in registry.models.values()
                  if m.substrate == Substrate.OPINION_SPACE]
        assert len(models) >= 1, "Need at least 1 opinion_space model (HK)"

    def test_five_substrate_types_total(self, registry):
        substrates = {m.substrate for m in registry.models.values()}
        assert len(substrates) == 5, f"Expected 5 substrate types, got {substrates}"


# ---------------------------------------------------------------------------
# Model and detector counts
# ---------------------------------------------------------------------------

class TestRegistryCounts:
    """Verify correct number of models and detectors registered."""

    def test_model_count(self, registry):
        # 11 model configurations (some models appear in multiple configs)
        assert len(registry.models) >= 11, \
            f"Expected ≥11 models, got {len(registry.models)}"

    def test_detector_count(self, registry):
        # 10 detectors (P1, P5, P6, P9, P13, P14, P15, P21, P27, P31)
        assert len(registry.detectors) >= 10, \
            f"Expected ≥10 detectors, got {len(registry.detectors)}"


# ---------------------------------------------------------------------------
# Compatibility matrix
# ---------------------------------------------------------------------------

class TestCompatibility:
    """Verify the 24 compatible pairs and substrate filtering."""

    def test_total_compatible_pairs(self, registry):
        compatible = []
        for m_name, m_cfg in registry.models.items():
            for d_name, d_cfg in registry.detectors.items():
                if d_cfg.required_observables.issubset(m_cfg.provides):
                    compatible.append((m_name, d_name))
        assert len(compatible) == 24, \
            f"Expected 24 compatible pairs, got {len(compatible)}: {compatible}"

    def test_substrate_mismatch_count(self, registry):
        total = len(registry.models) * len(registry.detectors)
        compatible = sum(
            1 for m in registry.models.values()
            for d in registry.detectors.values()
            if d.required_observables.issubset(m.provides)
        )
        mismatches = total - compatible
        assert mismatches > 0, "Should have substrate mismatches"


# ---------------------------------------------------------------------------
# Canonical positive pairs (each detector detects on its canonical model)
# ---------------------------------------------------------------------------

class TestCanonicalPairs:
    """Verify that canonical model-detector pairs are compatible."""

    @pytest.mark.parametrize("model_key,detector_key", [
        # Each detector's canonical positive model
        ("schelling", "P1"),
        ("vicsek_ordered", "P5"),
        ("dorsogna_mill", "P6"),
        ("kuramoto_sync", "P9"),
        ("gh_spiral", "P13"),
        ("btw_sandpile", "P14"),
        ("gol_rpent", "P15"),
        ("hk_polar", "P21"),
        ("nowak_may", "P27"),
    ])
    def test_canonical_pair_compatible(self, registry, model_key, detector_key):
        """Each canonical model must be substrate-compatible with its detector."""
        # Find model (name matching may be approximate)
        model = None
        for name, cfg in registry.models.items():
            if model_key in name.lower().replace("-", "_").replace(" ", "_"):
                model = cfg
                break
        assert model is not None, f"Model '{model_key}' not found in registry"

        detector = None
        for name, cfg in registry.detectors.items():
            if detector_key.lower() in name.lower():
                detector = cfg
                break
        assert detector is not None, f"Detector '{detector_key}' not found in registry"

        assert detector.required_observables.issubset(model.provides), \
            f"{model_key} does not provide {detector.required_observables - model.provides} for {detector_key}"


# ---------------------------------------------------------------------------
# Observable filtering — P27 uses coop_fraction, not grid
# ---------------------------------------------------------------------------

class TestObservableGuards:
    """Verify fine-grained observable filtering beyond substrate type."""

    def test_p27_requires_coop_fraction(self, registry):
        """P27 must require coop_fraction to avoid matching other lattice_2d models."""
        p27 = None
        for name, cfg in registry.detectors.items():
            if "p27" in name.lower() or "27" in name:
                p27 = cfg
                break
        assert p27 is not None, "P27 not found in registry"
        assert "coop_fraction" in p27.required_observables, \
            f"P27 should require coop_fraction, has {p27.required_observables}"

    def test_p27_not_compatible_with_gh(self, registry):
        """GH (lattice_2d) should NOT match P27 despite same substrate."""
        gh = None
        p27 = None
        for name, cfg in registry.models.items():
            if "gh" in name.lower() or "greenberg" in name.lower():
                gh = cfg
                break
        for name, cfg in registry.detectors.items():
            if "p27" in name.lower() or "27" in name:
                p27 = cfg
                break
        assert gh is not None and p27 is not None
        assert not p27.required_observables.issubset(gh.provides), \
            "P27 should NOT be compatible with GH"

    def test_p27_compatible_with_nowak_may(self, registry):
        """Nowak-May (lattice_2d + coop_fraction) SHOULD match P27."""
        nm = None
        p27 = None
        for name, cfg in registry.models.items():
            if "nowak" in name.lower():
                nm = cfg
                break
        for name, cfg in registry.detectors.items():
            if "p27" in name.lower() or "27" in name:
                p27 = cfg
                break
        assert nm is not None and p27 is not None
        assert p27.required_observables.issubset(nm.provides), \
            f"Nowak-May should provide {p27.required_observables}"

    def test_p21_compatible_with_hk(self, registry):
        """HK (opinion_space) should match P21."""
        hk = None
        p21 = None
        for name, cfg in registry.models.items():
            if "hk" in name.lower() or "hegselmann" in name.lower():
                hk = cfg
                break
        for name, cfg in registry.detectors.items():
            if "p21" in name.lower() or "21" in name:
                p21 = cfg
                break
        assert hk is not None and p21 is not None
        assert p21.required_observables.issubset(hk.provides), \
            f"HK should provide {p21.required_observables}"
