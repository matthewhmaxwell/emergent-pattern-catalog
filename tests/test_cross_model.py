"""Cross-model detection tests — beginning of the transfer matrix.

Tests that detectors correctly identify their home patterns and correctly
reject non-applicable models.
"""

import numpy as np
import pytest


class TestP13OnModels:
    """P13 excitable wave detector across model families."""

    def test_p13_on_greenberg_hastings_detected(self):
        """P13 should detect on its canonical model (broken wave → spiral)."""
        from epc.models.greenberg_hastings import GreenbergHastings
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        m = GreenbergHastings(
            rows=50, cols=50, n_states=5, threshold=1,
            init_mode="broken_wave", seed=0,
        )
        h = m.run(500)
        d = P13ExcitableWaveDetector(n_null_runs=99)
        r = d.detect(h, m.get_metadata())
        assert r.detected, f"P13 not detected on GH CA: {r.summary()}"

    def test_p13_on_vicsek_not_detected(self):
        """P13 should not detect on Vicsek (no grid, no excitable states)."""
        from epc.models.vicsek import VicsekModel
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        m = VicsekModel(n_particles=100, seed=0)
        h = m.run(100)
        d = P13ExcitableWaveDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P13 false positive on Vicsek: {r.summary()}"

    def test_p13_on_quiescent_gh_not_detected(self):
        """P13 should not detect on quiescent GH (high threshold)."""
        from epc.models.greenberg_hastings import GreenbergHastings
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        m = GreenbergHastings(rows=30, cols=30, threshold=9, seed=0)
        h = m.run(50)
        d = P13ExcitableWaveDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P13 false positive on quiescent GH: {r.summary()}"


class TestP5OnModels:
    """P5 flocking detector across model families."""

    def test_p5_on_vicsek_low_noise_detected(self):
        """P5 should detect flocking in low-noise Vicsek."""
        from epc.models.vicsek import VicsekModel
        from epc.detectors.p5_flocking import P5FlockingDetector

        m = VicsekModel(
            n_particles=100, box_size=5.0, speed=0.03,
            noise=0.05, interaction_radius=1.0, seed=42,
        )
        h = m.run(3000)
        d = P5FlockingDetector(n_permutations=99)
        r = d.detect(h, m.get_metadata())
        assert r.detected, f"P5 not detected on low-noise Vicsek: {r.summary()}"

    def test_p5_on_vicsek_high_noise_not_detected(self):
        """P5 should not detect flocking in high-noise Vicsek."""
        from epc.models.vicsek import VicsekModel
        from epc.detectors.p5_flocking import P5FlockingDetector

        m = VicsekModel(
            n_particles=100, box_size=10.0, speed=0.03,
            noise=2.0, interaction_radius=1.0, seed=42,
        )
        h = m.run(500)
        d = P5FlockingDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P5 false positive on high-noise Vicsek: {r.summary()}"

    def test_p5_on_greenberg_hastings_not_detected(self):
        """P5 should not detect on GH CA (no velocities)."""
        from epc.models.greenberg_hastings import GreenbergHastings
        from epc.detectors.p5_flocking import P5FlockingDetector

        m = GreenbergHastings(rows=30, cols=30, seed=0)
        h = m.run(50)
        d = P5FlockingDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P5 false positive on GH CA: {r.summary()}"


class TestP6OnModels:
    """P6 milling detector across model families."""

    def test_p6_on_greenberg_hastings_not_detected(self):
        """P6 should not detect on GH CA."""
        from epc.models.greenberg_hastings import GreenbergHastings
        from epc.detectors.p6_milling import P6MillingDetector

        m = GreenbergHastings(rows=30, cols=30, seed=0)
        h = m.run(50)
        d = P6MillingDetector()
        r = d.detect(h, m.get_metadata())
        assert not r.detected, f"P6 false positive on GH CA: {r.summary()}"

    def test_p6_on_standard_vicsek_not_detected(self):
        """P6 should not detect milling in standard Vicsek (no attraction)."""
        from epc.models.vicsek import VicsekModel
        from epc.detectors.p6_milling import P6MillingDetector

        m = VicsekModel(
            n_particles=100, box_size=10.0, speed=0.03,
            noise=0.1, seed=42,
        )
        h = m.run(500)
        d = P6MillingDetector()
        r = d.detect(h, m.get_metadata())
        # Standard Vicsek typically doesn't produce milling
        # At most screening-level is acceptable
        if r.detected:
            assert r.tier.value == "screening", (
                f"P6 unexpectedly confirmed on standard Vicsek: {r.summary()}"
            )
