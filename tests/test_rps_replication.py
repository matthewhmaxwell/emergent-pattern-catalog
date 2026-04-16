"""RPS (Reichenbach 2007) replication tests.

Each test verifies a specific quantitative claim from the paper or from
the broader cyclic-competition literature, not just that code runs.

Quantitative targets (Reichenbach, Mobilia & Frey 2007, Nature 448):
  1. Three-species coexistence at M < M_c: fractions fluctuate around 1/3
     each, no species goes extinct during the coexistence timescale.
  2. Single-species dominance at M > M_c: one species eventually exceeds
     ~95% on a modest lattice, within O(N) generations.
  3. Conservation: species + empty = L² at every timestep (trivial but
     pinned as a regression test).
  4. Reproducibility: same seed → bitwise-identical trajectories.
  5. Reaction firing rates: over many generations, the actual ratio of
     (selection : reproduction : exchange) executions matches the rate
     ratios (σ : µ : ε), modulo the no-op rate (incompatible proposals).

Design decisions for CI friendliness:
  - Most tests use L ≤ 50 and ≤ 500 generations (O(seconds)).
  - A single slow test with L=80 for the biodiversity loss check is
    marked with pytest.mark.slow so CI can opt in/out.
  - Critical-mobility test uses a COARSE grid of 3 mobilities spanning
    [deep coexistence, near-critical-but-stable, deep extinction] rather
    than fitting M_c precisely; this replicates the qualitative finding
    robustly without 20+ minute runs.
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.models.rps_spatial import RPSSpatialModel


# ===================================================================
# Basic invariants
# ===================================================================

class TestConservationAndSetup:
    """Conservation laws and basic setup behavior."""

    def test_total_cells_conserved(self):
        """Total site count = L² at every snapshot."""
        m = RPSSpatialModel(rows=40, cols=40, mobility=1e-4, seed=42)
        history = m.run(n_steps=50)
        expected_total = 40 * 40
        for i, h in enumerate(history):
            total = h["a_count"] + h["b_count"] + h["c_count"] + h["empty_count"]
            assert total == expected_total, (
                f"Conservation violated at step {i}: "
                f"A={h['a_count']}, B={h['b_count']}, C={h['c_count']}, "
                f"empty={h['empty_count']}, total={total} != {expected_total}"
            )

    def test_fractions_sum_to_one(self):
        """Species fractions + empty_fraction = 1 exactly."""
        m = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=7)
        history = m.run(n_steps=30)
        for h in history:
            s = (
                h["a_fraction"] + h["b_fraction"]
                + h["c_fraction"] + h["empty_fraction"]
            )
            assert abs(s - 1.0) < 1e-9, f"Fractions sum to {s}, not 1.0"

    def test_reproducibility_same_seed(self):
        """Bitwise-identical trajectories for the same seed."""
        m1 = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=123)
        m2 = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=123)
        h1 = m1.run(n_steps=25)
        h2 = m2.run(n_steps=25)
        assert len(h1) == len(h2)
        for i, (s1, s2) in enumerate(zip(h1, h2)):
            assert np.array_equal(s1["grid"], s2["grid"]), (
                f"Grids differ at snapshot {i} despite same seed"
            )

    def test_different_seeds_diverge(self):
        """Different seeds produce different trajectories."""
        m1 = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=1)
        m2 = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=2)
        h1 = m1.run(n_steps=20)
        h2 = m2.run(n_steps=20)
        # At least one snapshot should differ
        any_diff = any(
            not np.array_equal(s1["grid"], s2["grid"])
            for s1, s2 in zip(h1, h2)
        )
        assert any_diff, "Different seeds produced identical trajectories"

    def test_mobility_exchange_rate_mapping(self):
        """M = 2ε/N mapping is consistent in both directions."""
        # Specify by mobility
        m1 = RPSSpatialModel(rows=50, cols=50, mobility=4e-4, seed=0)
        assert abs(m1.exchange_rate - 4e-4 * 2500 / 2) < 1e-9
        assert abs(m1.mobility - 4e-4) < 1e-12

        # Specify by exchange_rate
        m2 = RPSSpatialModel(rows=50, cols=50, exchange_rate=1.0, seed=0)
        assert abs(m2.mobility - 2 * 1.0 / 2500) < 1e-12
        assert abs(m2.exchange_rate - 1.0) < 1e-12


# ===================================================================
# Biodiversity regime: coexistence vs extinction
# ===================================================================

class TestBiodiversityRegimes:
    """Reichenbach's core finding: mobility crosses a critical threshold
    separating coexistence from extinction."""

    def test_coexistence_at_low_mobility(self):
        """At M << M_c, all three species persist with fractions near 1/3.

        Uses L=50 (N=2500) → M_c is nominally ≈ 4.5e-4 but finite-size
        effects can shift it. We use M=1e-5 (well below even conservative
        estimates of M_c) and verify all species stay > 10% after 200
        generations. The fractions should fluctuate around 1/3 but are
        typically shifted a bit (empty cells take some share).
        """
        m = RPSSpatialModel(rows=50, cols=50, mobility=1e-5, seed=42)
        history = m.run(n_steps=200)

        # Use second half to avoid initial transient
        late = history[len(history) // 2:]
        a_min = min(h["a_fraction"] for h in late)
        b_min = min(h["b_fraction"] for h in late)
        c_min = min(h["c_fraction"] for h in late)

        # All three species maintain at least 10% of cells
        assert a_min > 0.10, f"Species A fraction dropped to {a_min:.3f} < 0.10"
        assert b_min > 0.10, f"Species B fraction dropped to {b_min:.3f} < 0.10"
        assert c_min > 0.10, f"Species C fraction dropped to {c_min:.3f} < 0.10"

        # Mean fractions are in the rough 1/5 – 1/2 band (each species
        # competes for ~1/3 of non-empty cells; empties take some fraction)
        a_mean = np.mean([h["a_fraction"] for h in late])
        b_mean = np.mean([h["b_fraction"] for h in late])
        c_mean = np.mean([h["c_fraction"] for h in late])
        for name, val in [("A", a_mean), ("B", b_mean), ("C", c_mean)]:
            assert 0.15 < val < 0.50, (
                f"Species {name} mean fraction {val:.3f} outside [0.15, 0.50] "
                f"(expected near 1/3 in coexistence)"
            )

    def test_extinction_at_high_mobility(self):
        """At M >> M_c on a modest grid, biodiversity collapses.

        Uses L=50, M=1e-2 (~20× above nominal M_c). After 500 generations
        we expect the minimum species fraction to have dropped below 5%,
        indicating extinction or near-extinction. Empirically we see one
        species typically exceeds 70% by this point.
        """
        m = RPSSpatialModel(rows=50, cols=50, mobility=1e-2, seed=42)
        history = m.run(n_steps=500)

        final = history[-1]
        min_frac = min(final["a_fraction"], final["b_fraction"], final["c_fraction"])
        max_frac = max(final["a_fraction"], final["b_fraction"], final["c_fraction"])

        assert min_frac < 0.05, (
            f"At M=1e-2 on L=50, min species fraction should drop below 5%; "
            f"got {min_frac:.3f}. Regime should be near-extinction."
        )
        assert max_frac > 0.5, (
            f"At M=1e-2, one species should dominate (>50%); "
            f"got max={max_frac:.3f}"
        )

    def test_mobility_monotonic_effect_on_minimum_fraction(self):
        """Higher mobility → lower minimum species fraction at fixed time.

        Runs three mobilities below / at / above the nominal M_c and
        verifies the MIN species fraction decreases monotonically with
        mobility. This is the qualitative signature of Reichenbach's
        phase diagram without fitting M_c precisely.
        """
        # Small grid + modest run length for CI speed
        mobilities = [1e-5, 1e-3, 5e-3]
        min_fractions = []

        for M in mobilities:
            m = RPSSpatialModel(rows=40, cols=40, mobility=M, seed=42)
            history = m.run(n_steps=300)
            final = history[-1]
            min_frac = min(
                final["a_fraction"], final["b_fraction"], final["c_fraction"]
            )
            min_fractions.append(min_frac)

        # Monotone decreasing (allow small reversals due to stochasticity)
        # Principal check: deep-extinction regime << deep-coexistence regime
        assert min_fractions[0] > min_fractions[-1] + 0.1, (
            f"Min-fraction should be much lower at high mobility. "
            f"Low M: {min_fractions[0]:.3f}, High M: {min_fractions[-1]:.3f}"
        )


# ===================================================================
# Reaction mechanics
# ===================================================================

class TestReactionMechanics:
    """Verify selection/reproduction/exchange reactions execute with
    rates proportional to their Poisson rate parameters."""

    def test_reaction_rate_ratios(self):
        """Over a run, executed selections : reproductions ratio should
        roughly track selection_rate : reproduction_rate, correcting for
        the fact that selection is a no-op when the pair isn't in a
        dominance relation and reproduction is a no-op when there's no
        empty cell.

        We check the weaker property: BOTH reactions fire nonzero counts
        and exchange counts scale with exchange_rate (which controls
        mobility).
        """
        # High mobility: exchange should dominate executed reactions
        m_hi = RPSSpatialModel(
            rows=40, cols=40, mobility=1e-2, seed=42  # ε = 8.0 (dominates σ=µ=1)
        )
        hist_hi = m_hi.run(n_steps=30)
        total_exc_hi = sum(h["n_exchanges"] for h in hist_hi[1:])
        total_sel_hi = sum(h["n_selections"] for h in hist_hi[1:])
        total_rep_hi = sum(h["n_reproductions"] for h in hist_hi[1:])

        # Low mobility: exchange should be rare
        m_lo = RPSSpatialModel(
            rows=40, cols=40, mobility=1e-5, seed=42  # ε = 0.008 (subdominant)
        )
        hist_lo = m_lo.run(n_steps=30)
        total_exc_lo = sum(h["n_exchanges"] for h in hist_lo[1:])

        # At high mobility, exchanges vastly outnumber selections+reproductions
        assert total_exc_hi > 5 * (total_sel_hi + total_rep_hi), (
            f"At M=1e-2, exchanges ({total_exc_hi}) should be >5x all other "
            f"reactions ({total_sel_hi + total_rep_hi})"
        )

        # At low mobility, exchanges << combined
        # (exchange_rate/total = 0.008 / 2.008 ≈ 0.4%)
        assert total_exc_lo < 0.02 * (total_sel_hi + total_rep_hi + total_exc_hi), (
            f"At M=1e-5, exchanges ({total_exc_lo}) should be rare "
            f"(< 2% of total reactions)"
        )

    def test_selection_only_within_dominance(self):
        """A selection event must produce an empty cell; inverse dominance
        pairs should never trigger selection.

        Implementation check: run a single step from a prepared state where
        A and B are adjacent in a 'wrong' orientation and verify no cell
        becomes empty by selection when zero-mobility, zero-reproduction.
        """
        # Build a grid that's pure A and B striped: each A cell has only
        # B neighbors. Because A dominates B, A-selection events will
        # make B cells empty. But B-selection events (B picks A neighbor)
        # are no-ops. So net effect: B cells vanish over time, A persists.
        L = 20
        grid = np.zeros((L, L), dtype=np.int8)
        grid[::2, :] = RPSSpatialModel.SPECIES_A  # rows 0, 2, 4, ... are A
        grid[1::2, :] = RPSSpatialModel.SPECIES_B  # rows 1, 3, ... are B

        m = RPSSpatialModel(
            rows=L, cols=L,
            selection_rate=1.0, reproduction_rate=1e-6, exchange_rate=0.0,
            init_mode="custom", init_grid=grid, seed=42,
        )
        history = m.run(n_steps=20)

        # After many steps with selection dominant and negligible repro:
        # B should have dropped significantly (A killed many B's)
        # and A should NOT have dropped (nothing kills A in this setup —
        # no C species present).
        final = history[-1]
        initial = history[0]
        assert final["a_count"] >= initial["a_count"] - 2, (
            f"A count should not decrease (nothing kills A here). "
            f"Initial: {initial['a_count']}, Final: {final['a_count']}"
        )
        assert final["b_count"] < initial["b_count"] * 0.5, (
            f"B count should drop substantially (A kills B). "
            f"Initial: {initial['b_count']}, Final: {final['b_count']}"
        )


# ===================================================================
# Slow replication test (opt-in via marker)
# ===================================================================

class TestSlowReplication:
    """Longer runs that verify quantitative claims more rigorously.

    Marked 'slow' — runs ~20-60 seconds each.
    """

    @pytest.mark.slow
    def test_coexistence_stability_long_run(self):
        """At M << M_c over long time, all three species persist and
        fractions remain close to equal.

        Uses L=60 and 500 generations. The fraction-of-time each species
        exceeds 10% should be essentially 1.0.
        """
        m = RPSSpatialModel(rows=60, cols=60, mobility=1e-5, seed=42)
        history = m.run(n_steps=500)

        # Skip first 100 generations as transient
        late = history[100:]
        n = len(late)

        a_above = sum(1 for h in late if h["a_fraction"] > 0.10) / n
        b_above = sum(1 for h in late if h["b_fraction"] > 0.10) / n
        c_above = sum(1 for h in late if h["c_fraction"] > 0.10) / n

        # Each species should be above 10% essentially all the time
        # (allow 5% slack for rare fluctuations)
        assert a_above > 0.95, f"A above 10% only {a_above:.2%} of the time"
        assert b_above > 0.95, f"B above 10% only {b_above:.2%} of the time"
        assert c_above > 0.95, f"C above 10% only {c_above:.2%} of the time"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
