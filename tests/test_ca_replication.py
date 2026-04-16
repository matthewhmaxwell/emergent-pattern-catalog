"""Quantitative replication tests for excitable-media and computation CAs.

Tests verify that our Greenberg-Hastings and Game of Life implementations
reproduce PUBLISHED QUANTITATIVE RESULTS from the canonical literature.
Each test targets a specific numerical claim from the source paper or
LifeWiki reference, not just qualitative behavior.

References:
  Greenberg, J.M. & Hastings, S.P. (1978). "Spatial patterns for discrete
  models of diffusion in excitable media." SIAM J. Applied Math 34, 515-523.

  Fisch, R., Gravner, J. & Griffeath, D. (1991). "Threshold-range scaling
  of excitable cellular automata." Statistics and Computing 1, 23-39.

  Gardner, M. (1970). "The fantastic combinations of John Conway's new
  solitaire game 'life'." Scientific American 223(4), 120-123.

  LifeWiki (conwaylife.com): canonical patterns and trajectories.
"""

import numpy as np

from epc.models.greenberg_hastings import GreenbergHastings
from epc.models.game_of_life import GameOfLife


# ============================================================
# Greenberg-Hastings Quantitative Replication
# ============================================================


class TestGHWaveSpeed:
    """Wave propagation speed in GH must equal 1.0 cells/step (L1) in VN
    neighborhood and sqrt(2) cells/step (L2) in Moore neighborhood.

    This is a direct consequence of the synchronous update rule: a single
    wavefront cell excites all neighbors at the next timestep, so the
    front advances by exactly 1 step along the shortest neighbor
    connection.
    """

    def _measure_speed(self, model, n_steps, fit_start=2, fit_end=25):
        """Measure wavefront expansion speed by linear fit on max radius."""
        history = model.run(n_steps, record_every=1)
        center = np.array([model.rows // 2, model.cols // 2], dtype=float)
        radii = []
        for t, state in enumerate(history):
            touched = np.argwhere(state["grid"] != 0)
            if len(touched) > 0:
                dists = np.sqrt(np.sum((touched - center) ** 2, axis=1))
                radii.append((t, float(dists.max())))

        if len(radii) < fit_start + 5:
            return None, None
        data = np.array([r for r in radii if fit_start <= r[0] <= fit_end])
        if len(data) < 5:
            return None, None
        ts, rs = data[:, 0], data[:, 1]
        slope, intercept = np.polyfit(ts, rs, 1)
        predicted = slope * ts + intercept
        ss_res = np.sum((rs - predicted) ** 2)
        ss_tot = np.sum((rs - np.mean(rs)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return slope, r2

    def test_von_neumann_wave_speed(self):
        """VN neighborhood: wave expands at exactly 1.0 cell/step."""
        m = GreenbergHastings(
            rows=101, cols=101, n_states=3, threshold=1,
            neighborhood="von_neumann", init_mode="single_seed",
            boundary="fixed", seed=42,
        )
        speed, r2 = self._measure_speed(m, 30)
        assert speed is not None, "Could not measure speed"
        assert abs(speed - 1.0) < 0.01, f"VN wave speed = {speed:.4f}, expected 1.0"
        assert r2 > 0.999, f"VN speed linearity R²={r2:.4f}, expected > 0.999"

    def test_moore_wave_speed(self):
        """Moore neighborhood: wave expands at sqrt(2) ≈ 1.414 cells/step
        along diagonals (L2 max distance)."""
        m = GreenbergHastings(
            rows=101, cols=101, n_states=3, threshold=1,
            neighborhood="moore", init_mode="single_seed",
            boundary="fixed", seed=42,
        )
        speed, r2 = self._measure_speed(m, 30)
        assert speed is not None
        expected = np.sqrt(2)
        assert abs(speed - expected) < 0.02, \
            f"Moore wave speed = {speed:.4f}, expected {expected:.4f}"
        assert r2 > 0.999, f"Moore speed linearity R²={r2:.4f}"

    def test_wave_speed_independent_of_kappa(self):
        """Wave speed does not depend on κ (number of states).

        The refractory length only affects how long cells wait to be
        re-excitable — it doesn't affect the leading edge speed.
        """
        speeds = []
        for kappa in [3, 5, 10]:
            m = GreenbergHastings(
                rows=101, cols=101, n_states=kappa, threshold=1,
                neighborhood="moore", init_mode="single_seed",
                boundary="fixed", seed=42,
            )
            speed, r2 = self._measure_speed(m, 30)
            speeds.append((kappa, speed, r2))
            assert r2 > 0.999

        # All three speeds should be essentially equal
        spread = max(s[1] for s in speeds) - min(s[1] for s in speeds)
        assert spread < 0.01, \
            f"Wave speeds vary with κ: {speeds} (spread={spread:.4f})"


class TestGHThresholdRule:
    """Threshold controls whether an excited cell can trigger neighbors.

    A single excited cell has only 1 excited neighbor to anyone. With
    θ ≥ 2, a lone spark cannot propagate.
    """

    def test_single_seed_propagates_at_theta_1(self):
        """θ=1: single seed propagates (has 1 excited neighbor ≥ 1)."""
        m = GreenbergHastings(
            rows=51, cols=51, n_states=3, threshold=1,
            neighborhood="von_neumann", init_mode="single_seed",
            boundary="fixed", seed=42,
        )
        h = m.run(20)
        # After 20 steps, the wave should have spread significantly
        activity = sum(s["excited_count"] for s in h)
        assert activity > 50, "Single seed should propagate at θ=1"

    def test_single_seed_fails_at_high_threshold(self):
        """θ ≥ 2: single seed cannot propagate (only 1 excited neighbor)."""
        for nbhd in ["von_neumann", "moore"]:
            for theta in [2, 3, 4]:
                m = GreenbergHastings(
                    rows=51, cols=51, n_states=3, threshold=theta,
                    neighborhood=nbhd, init_mode="single_seed",
                    boundary="fixed", seed=42,
                )
                h = m.run(20)
                # The initial cell cycles through excited→refractory→rest
                # but no new cells should ever become excited
                total_activations = sum(s["wavefront_count"] for s in h)
                assert total_activations == 0, \
                    f"{nbhd} θ={theta}: had {total_activations} activations, " \
                    f"should be 0 (lone seed has only 1 neighbor)"

    def test_all_resting_stays_dead(self):
        """Empty grid with all cells resting must remain empty forever."""
        m = GreenbergHastings(
            rows=20, cols=20, n_states=5, threshold=1,
            neighborhood="moore", init_mode="custom",
            init_grid=np.zeros((20, 20), dtype=int), seed=42,
        )
        h = m.run(50)
        for s in h:
            assert s["excited_count"] == 0
            assert s["refractory_count"] == 0


class TestGHSpiralPersistence:
    """GH spirals from broken-wave IC must persist indefinitely (periodic
    bdry) and show a stable fundamental period.

    The spiral period at a fixed point far from the core equals κ+1 for
    the minimal spiral (κ=3): excitation cycles through all κ states,
    but must also wait one step for a neighbor to re-excite it.
    """

    def test_broken_wave_activity_persists(self):
        """Activity must remain > 0 for at least 3×T_prop steps."""
        m = GreenbergHastings(
            rows=60, cols=60, n_states=3, threshold=1,
            neighborhood="von_neumann", init_mode="broken_wave",
            boundary="periodic", seed=42,
        )
        h = m.run(200, record_every=1)
        # After any transient, activity in last 50 steps should be > 0
        late_activity = [s["excited_count"] for s in h[-50:]]
        assert min(late_activity) > 0, \
            f"Spiral activity died: min(excited in t=150..200) = {min(late_activity)}"

    def test_spiral_period_at_kappa_3(self):
        """κ=3, θ=1, VN: period at fixed point is exactly 4.

        This is the classical minimal GH spiral. The +1 over κ comes from
        the requirement that a neighbor be in the excited state one step
        before the cell itself can be re-excited — adding one "waiting"
        step to the bare state cycle.
        """
        m = GreenbergHastings(
            rows=60, cols=60, n_states=3, threshold=1,
            neighborhood="von_neumann", init_mode="broken_wave",
            boundary="periodic", seed=42,
        )
        h = m.run(100, record_every=1)

        # Sample several interior points far from spiral core
        periods_observed = []
        for point in [(40, 40), (20, 20), (50, 10), (10, 50)]:
            series = np.array([s["grid"][point] for s in h])
            excited_events = np.where(series == 1)[0]
            if len(excited_events) >= 4:
                # Use steady-state periods (skip first 2 for transient)
                ss = np.diff(excited_events[-5:])
                periods_observed.extend(ss.tolist())

        assert len(periods_observed) > 0, "No spiral re-excitations found"
        # All observed steady-state periods should be 4
        unique_periods = set(periods_observed)
        assert unique_periods == {4}, \
            f"Expected all periods = 4 for κ=3 GH spiral, got {unique_periods}"


class TestGHSelfOrganization:
    """Random IC → ordered spiral state (Greenberg-Hastings self-organization).

    Starting from a random configuration, the GH CA self-organizes into
    a collection of rotating spirals whose activity persists indefinitely.
    This is one of the key published claims.
    """

    def test_random_ic_produces_persistent_activity(self):
        """Across 10 seeds, random IC must produce persistent activity."""
        persistent = 0
        for seed in range(10):
            m = GreenbergHastings(
                rows=60, cols=60, n_states=6, threshold=1,
                neighborhood="moore", init_mode="random",
                init_density=0.5, boundary="periodic", seed=seed,
            )
            h = m.run(150, record_every=1)
            # Activity in both early and late windows
            early = np.mean([s["excited_count"] for s in h[:50]])
            late = np.mean([s["excited_count"] for s in h[100:]])
            if early > 0 and late > 0:
                persistent += 1
        assert persistent >= 9, \
            f"GH should self-organize in most runs; persistent in {persistent}/10"


# ============================================================
# Game of Life Quantitative Replication
# ============================================================


class TestGoLStillLifes:
    """Canonical still lifes must remain unchanged forever."""

    def _place_pattern(self, pattern_cells, grid_size=10, offset=(4, 4)):
        """Place a pattern on an empty grid at the given offset."""
        grid = np.zeros((grid_size, grid_size), dtype=int)
        r0, c0 = offset
        for dr, dc in pattern_cells:
            grid[r0 + dr, c0 + dc] = 1
        return grid

    def test_block_stable(self):
        """Block (2×2): stable for at least 20 steps."""
        block = [(0, 0), (0, 1), (1, 0), (1, 1)]
        grid = self._place_pattern(block)
        m = GameOfLife(rows=10, cols=10, init_mode="custom", init_grid=grid, seed=0)
        h = m.run(20)
        for i, state in enumerate(h):
            assert np.array_equal(state["grid"], h[0]["grid"]), \
                f"Block not stable at step {i}"

    def test_beehive_stable(self):
        """Beehive (6 cells): stable for at least 20 steps."""
        beehive = [(0, 1), (0, 2), (1, 0), (1, 3), (2, 1), (2, 2)]
        grid = self._place_pattern(beehive)
        m = GameOfLife(rows=10, cols=10, init_mode="custom", init_grid=grid, seed=0)
        h = m.run(20)
        for state in h:
            assert np.array_equal(state["grid"], h[0]["grid"])

    def test_loaf_stable(self):
        """Loaf (7 cells): stable for at least 20 steps."""
        loaf = [(0, 1), (0, 2), (1, 0), (1, 3), (2, 1), (2, 3), (3, 2)]
        grid = self._place_pattern(loaf)
        m = GameOfLife(rows=10, cols=10, init_mode="custom", init_grid=grid, seed=0)
        h = m.run(20)
        for state in h:
            assert np.array_equal(state["grid"], h[0]["grid"])


class TestGoLOscillators:
    """Canonical oscillators must have exactly the published periods."""

    def _place(self, cells, grid_size=20, offset=(8, 8)):
        grid = np.zeros((grid_size, grid_size), dtype=int)
        r0, c0 = offset
        for dr, dc in cells:
            grid[r0 + dr, c0 + dc] = 1
        return grid

    def test_blinker_period_2(self):
        """Blinker (period 2): horizontal 3-cell row ↔ vertical 3-cell column."""
        grid = self._place([(0, 0), (0, 1), (0, 2)])
        m = GameOfLife(rows=20, cols=20, init_mode="custom", init_grid=grid, seed=0)
        h = m.run(10)
        # t and t+2 equal, t and t+1 differ
        for t in range(8):
            assert np.array_equal(h[t]["grid"], h[t + 2]["grid"]), \
                f"Blinker not period-2 at step {t}"
            assert not np.array_equal(h[t]["grid"], h[t + 1]["grid"]), \
                f"Blinker has shorter period at step {t}"

    def test_toad_period_2(self):
        """Toad (period 2): 6-cell oscillator."""
        toad = [(0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2)]
        grid = self._place(toad)
        m = GameOfLife(rows=20, cols=20, init_mode="custom", init_grid=grid, seed=0)
        h = m.run(10)
        for t in range(8):
            assert np.array_equal(h[t]["grid"], h[t + 2]["grid"])

    def test_beacon_period_2(self):
        """Beacon (period 2): two blocks that flash."""
        beacon = [(0, 0), (0, 1), (1, 0), (1, 1),
                  (2, 2), (2, 3), (3, 2), (3, 3)]
        grid = self._place(beacon)
        m = GameOfLife(rows=20, cols=20, init_mode="custom", init_grid=grid, seed=0)
        h = m.run(10)
        for t in range(8):
            assert np.array_equal(h[t]["grid"], h[t + 2]["grid"])

    def test_pulsar_period_3(self):
        """Pulsar (period 3): 48-cell oscillator, one of the largest common ones."""
        def pulsar_cells():
            cells = []
            # Horizontal 3-cell bars
            for r in [0, 5, 7, 12]:
                for c in [2, 3, 4, 8, 9, 10]:
                    cells.append((r, c))
            # Vertical 3-cell bars
            for c in [0, 5, 7, 12]:
                for r in [2, 3, 4, 8, 9, 10]:
                    cells.append((r, c))
            return cells

        grid = self._place(pulsar_cells(), grid_size=30, offset=(8, 8))
        m = GameOfLife(rows=30, cols=30, init_mode="custom", init_grid=grid, seed=0)
        h = m.run(20)
        # t and t+3 equal
        for t in range(15):
            assert np.array_equal(h[t]["grid"], h[t + 3]["grid"]), \
                f"Pulsar not period-3 at step {t}"


class TestGoLSpaceships:
    """Spaceships must travel at published velocities."""

    def test_glider_speed_c_over_4_diagonal(self):
        """Glider: period 4, displaces 1 cell diagonally per period.

        Speed = c/4 where c = max speed 1 cell/step. So in 40 steps the
        center of mass should move exactly (10, 10) cells.
        """
        # Standard glider (bottom-right moving)
        glider_cells = [(0, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
        grid = np.zeros((50, 50), dtype=int)
        for dr, dc in glider_cells:
            grid[5 + dr, 5 + dc] = 1
        m = GameOfLife(rows=50, cols=50, init_mode="custom",
                       init_grid=grid, boundary="fixed", seed=0)
        h = m.run(40)

        # Center of mass at each recorded step
        def com(grid):
            cells = np.argwhere(grid == 1)
            return cells.mean(axis=0)

        com_0 = com(h[0]["grid"])
        com_40 = com(h[40]["grid"])
        dr = com_40[0] - com_0[0]
        dc = com_40[1] - com_0[1]

        # Expected: diagonal move of exactly 10 cells in each axis
        assert abs(dr - 10.0) < 0.01, f"Glider vertical displacement = {dr}, expected 10"
        assert abs(dc - 10.0) < 0.01, f"Glider horizontal displacement = {dc}, expected 10"

    def test_lwss_speed_c_over_2_orthogonal(self):
        """LWSS: period 4, displaces 2 cells per period in one axis.

        Speed = c/2. Over 20 steps = 5 periods, should move 10 cells.
        """
        m = GameOfLife(rows=50, cols=50, init_mode="lwss",
                       boundary="fixed", seed=0)
        h = m.run(20)

        def com(grid):
            cells = np.argwhere(grid == 1)
            return cells.mean(axis=0) if len(cells) > 0 else None

        com_0 = com(h[0]["grid"])
        com_20 = com(h[20]["grid"])
        assert com_0 is not None and com_20 is not None, "LWSS disappeared"

        displacement = np.linalg.norm(com_20 - com_0)
        # Should be approximately 10 cells (20 steps / 2 steps per cell)
        assert abs(displacement - 10.0) < 0.5, \
            f"LWSS displacement over 20 steps = {displacement:.2f}, expected 10"


class TestGoLRPentomino:
    """R-pentomino is Conway's canonical methuselah.

    Published LifeWiki results:
      - Peak population: 319 cells at generation 821
      - Stabilizes at generation 1103
      - Final stable population: 116 cells (including 6 escaping gliders)
    """

    def test_r_pentomino_peak_population(self):
        """Peak population should be 319 cells, occurring around step 821."""
        m = GameOfLife(rows=300, cols=300, init_mode="r_pentomino",
                       boundary="fixed", seed=0)
        h = m.run(900)
        pop = [s["alive_count"] for s in h]
        peak = max(pop)
        peak_step = int(np.argmax(pop))

        # Peak should be ~319 cells (small variation possible due to grid size)
        assert abs(peak - 319) <= 10, f"R-pentomino peak population = {peak}, expected ~319"
        # Peak should occur between step 800 and 850
        assert 800 <= peak_step <= 850, \
            f"R-pentomino peak at step {peak_step}, expected ~821"

    def test_r_pentomino_stabilizes_near_1103(self):
        """Population should stabilize around step 1103 with ~116 cells.

        On a fixed-boundary grid of finite size, some escaping gliders may
        annihilate at the boundary, so final count can be 110-116.
        """
        m = GameOfLife(rows=500, cols=500, init_mode="r_pentomino",
                       boundary="fixed", seed=0)
        h = m.run(1200)
        pop = [s["alive_count"] for s in h]

        # Check population stability between step 1103 and 1200
        late_pops = pop[1103:1200]
        spread = max(late_pops) - min(late_pops)
        # Should be "stable" modulo oscillator beats
        # (blinkers create period-2 oscillation in total count)
        assert spread <= 10, \
            f"Population not stable after step 1103: spread={spread}"

        # Published final count: 116. Allow 110-125 on finite grid.
        final = pop[-1]
        assert 105 <= final <= 125, \
            f"R-pentomino final population = {final}, expected ~116"
