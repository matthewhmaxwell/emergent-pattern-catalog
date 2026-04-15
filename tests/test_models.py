"""Unit tests for model implementations."""

import numpy as np
import pytest


# --- Greenberg-Hastings CA ---

class TestGreenbergHastingsCA:
    def test_setup_and_state_keys(self):
        from epc.models.greenberg_hastings import GreenbergHastings
        m = GreenbergHastings(rows=20, cols=20, seed=0)
        state = m.setup()
        assert "grid" in state
        assert "step" in state
        assert "excited_count" in state
        assert "refractory_count" in state
        assert "resting_count" in state
        assert "grid_dims" in state
        assert state["grid_dims"] == (20, 20)

    def test_state_counts_sum_to_grid_size(self):
        from epc.models.greenberg_hastings import GreenbergHastings
        m = GreenbergHastings(rows=30, cols=30, seed=42)
        h = m.run(50)
        for state in h:
            total = state["excited_count"] + state["refractory_count"] + state["resting_count"]
            assert total == 30 * 30, f"State counts don't sum to grid size: {total}"

    def test_high_threshold_quiescence(self):
        """With threshold > 8 (max Moore neighbors), no rest cell can excite."""
        from epc.models.greenberg_hastings import GreenbergHastings
        m = GreenbergHastings(rows=20, cols=20, threshold=9, init_density=0.05, seed=0)
        h = m.run(30)
        final = h[-1]
        assert final["excited_count"] == 0

    def test_broken_wave_sustains_activity(self):
        """Broken wave init should produce sustained spiral activity."""
        from epc.models.greenberg_hastings import GreenbergHastings
        m = GreenbergHastings(
            rows=50, cols=50, n_states=5, threshold=1,
            init_mode="broken_wave", seed=0,
        )
        h = m.run(200)
        late_excited = [s["excited_count"] for s in h[-50:]]
        assert max(late_excited) > 0, "Spiral activity died"

    def test_convergence_detection(self):
        from epc.models.greenberg_hastings import GreenbergHastings
        # threshold=9 means no propagation; need enough steps for all
        # refractory cells to decay
        m = GreenbergHastings(rows=10, cols=10, threshold=9, init_density=0.05, n_states=5, seed=0)
        m.run(50)
        assert m.is_converged()

    def test_timescale(self):
        from epc.models.greenberg_hastings import GreenbergHastings
        m = GreenbergHastings(rows=50, cols=50)
        ts = m.get_timescale()
        assert ts == 50.0

    def test_metadata(self):
        from epc.models.greenberg_hastings import GreenbergHastings
        m = GreenbergHastings(rows=40, cols=40, n_states=7, threshold=2)
        md = m.get_metadata()
        assert md["model_name"] == "greenberg_hastings"
        assert md["model_class"] == "excitable_ca"
        assert md["rows"] == 40
        assert md["cols"] == 40
        assert md["n_states"] == 7
        assert md["threshold"] == 2


# --- Vicsek Model ---

class TestVicsekModel:
    def test_setup_and_state_keys(self):
        from epc.models.vicsek import VicsekModel
        m = VicsekModel(n_particles=50, seed=0)
        state = m.setup()
        assert "positions" in state
        assert "velocities" in state
        assert "headings" in state
        assert "step" in state
        assert state["positions"].shape == (50, 2)
        assert state["velocities"].shape == (50, 2)

    def test_positions_stay_in_box(self):
        from epc.models.vicsek import VicsekModel
        m = VicsekModel(n_particles=100, box_size=10.0, speed=0.5, seed=42)
        h = m.run(200)
        for state in h:
            pos = state["positions"]
            assert np.all(pos >= 0), "Position below 0"
            assert np.all(pos < 10.0), "Position above box_size"

    def test_low_noise_produces_flocking(self):
        """Low noise should produce high polarization."""
        from epc.models.vicsek import VicsekModel
        m = VicsekModel(
            n_particles=100, box_size=5.0, speed=0.03,
            noise=0.05, interaction_radius=1.0, seed=42,
        )
        h = m.run(2000)
        late_states = h[-len(h) // 5:]
        phis = []
        for s in late_states:
            v = s["velocities"]
            speeds = np.linalg.norm(v, axis=1, keepdims=True)
            speeds = np.where(speeds == 0, 1, speeds)
            v_hat = v / speeds
            phi = float(np.linalg.norm(v_hat.mean(axis=0)))
            phis.append(phi)
        mean_phi = np.mean(phis)
        assert mean_phi > 0.6, f"Expected flocking (phi > 0.6), got {mean_phi:.3f}"

    def test_high_noise_produces_disorder(self):
        """High noise should produce low polarization."""
        from epc.models.vicsek import VicsekModel
        # Use low density (N=100, L=25 → ρ=0.16) and high noise for
        # definitive disorder. Validated test uses η=2π, N=300.
        m = VicsekModel(
            n_particles=100, box_size=25.0, speed=0.03,
            noise=5.0, interaction_radius=1.0, seed=42,
        )
        h = m.run(500)
        late_states = h[-100:]
        phis = []
        for s in late_states:
            v = s["velocities"]
            speeds = np.linalg.norm(v, axis=1, keepdims=True)
            speeds = np.where(speeds == 0, 1, speeds)
            v_hat = v / speeds
            phi = float(np.linalg.norm(v_hat.mean(axis=0)))
            phis.append(phi)
        mean_phi = np.mean(phis)
        assert mean_phi < 0.4, f"Expected disorder (phi < 0.4), got {mean_phi:.3f}"

    def test_constant_speed(self):
        """All agents should maintain constant speed."""
        from epc.models.vicsek import VicsekModel
        m = VicsekModel(n_particles=50, speed=0.5, seed=0)
        h = m.run(100)
        for state in h:
            speeds = np.linalg.norm(state["velocities"], axis=1)
            np.testing.assert_allclose(speeds, 0.5, atol=1e-10)

    def test_timescale(self):
        from epc.models.vicsek import VicsekModel
        m = VicsekModel(box_size=25.0, speed=0.5)
        ts = m.get_timescale()
        assert isinstance(ts, dict)
        assert ts["T_cross"] == 50.0

    def test_metadata(self):
        from epc.models.vicsek import VicsekModel
        m = VicsekModel(n_particles=200, noise=0.3)
        md = m.get_metadata()
        assert md["model_family"] == "vicsek"
        assert md["n_particles"] == 200
        assert md["noise"] == 0.3
