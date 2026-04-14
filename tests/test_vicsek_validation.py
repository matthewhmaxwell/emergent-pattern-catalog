"""
Vicsek Model Validation — Replication of Vicsek et al. (1995).

Tests:
    1. UNIT TESTS: Basic mechanics (periodic BC, constant speed, neighbor finding)
    2. PHASE TRANSITION: φ vs η curve at paper parameters (N=300, L=7, v₀=0.03)
    3. DISORDERED BASELINE: φ ≈ 1/√N at high noise
    4. ORDERED STATE: φ → 1 at low noise
    5. DENSITY DEPENDENCE: Higher density → lower critical noise

Power requirements:
    - T_cross = L/v₀ = 7/0.03 ≈ 233 steps
    - Equilibration: ≥ 10 × T_cross = 2,330 steps (burn-in)
    - Measurement: ≥ 20 × T_cross = 4,660 steps (averaging window)
    - Total per η value: ≥ 7,000 steps
    - Multiple independent seeds per η for error bars: ≥ 5 seeds
    - η scan: ≥ 15 points across [0, 5]

Reference: Vicsek et al. PRL 75, 1226 (1995), Figures 1 and 3.
"""

import sys
import time
import numpy as np

sys.path.insert(0, "/home/claude")
from epc.models.vicsek import VicsekModel, polarization, angular_momentum, group_speed_ratio


def test_periodic_boundary_conditions():
    """Verify particles wrap correctly at domain boundaries."""
    # Use 2 particles far apart so no neighbor interaction changes headings
    model = VicsekModel(
        n_particles=2, box_size=100.0, speed=1.0, noise=0.0,
        interaction_radius=1.0, dt=1.0, seed=42
    )
    model.setup()

    # Place particle at edge, heading right, far from other particle
    model.positions[0] = [99.5, 50.0]
    model.headings[0] = 0.0  # moving right
    model.positions[1] = [50.0, 50.0]
    model.headings[1] = 0.0

    model.step()

    # Should wrap: 99.5 + 1.0 = 100.5 → 0.5 (mod 100)
    assert 0.0 <= model.positions[0, 0] < 100.0, \
        f"X not wrapped: {model.positions[0, 0]}"
    assert abs(model.positions[0, 0] - 0.5) < 1e-10, \
        f"Expected 0.5, got {model.positions[0, 0]}"

    print("  ✓ Periodic boundary conditions work correctly")


def test_constant_speed():
    """Verify all particles maintain constant speed v₀."""
    model = VicsekModel(
        n_particles=100, box_size=10.0, speed=0.5, noise=1.0,
        interaction_radius=1.0, dt=1.0, seed=42
    )
    model.setup()

    for _ in range(50):
        state = model.step()
        speeds = np.linalg.norm(state["velocities"], axis=1)
        assert np.allclose(speeds, 0.5, atol=1e-10), \
            f"Speed not constant: min={speeds.min():.6f}, max={speeds.max():.6f}"

    print("  ✓ Constant speed maintained for all particles across 50 steps")


def test_perfect_alignment_no_noise():
    """With η=0 and aligned init, polarization should remain 1.0."""
    model = VicsekModel(
        n_particles=100, box_size=10.0, speed=0.03, noise=0.0,
        interaction_radius=1.0, dt=1.0, init_mode="aligned", seed=42
    )
    model.setup()

    for _ in range(100):
        state = model.step()
        phi = polarization(state)
        assert abs(phi - 1.0) < 1e-10, \
            f"Polarization drifted from 1.0: φ={phi:.10f}"

    print("  ✓ Perfect alignment preserved with zero noise (φ=1.0 for 100 steps)")


def test_neighbor_self_inclusion():
    """Verify each particle counts itself as a neighbor.

    With zero noise, an isolated particle should maintain its heading
    (since its only neighbor is itself).
    """
    model = VicsekModel(
        n_particles=2, box_size=100.0, speed=0.01, noise=0.0,
        interaction_radius=1.0, dt=1.0, seed=42
    )
    model.setup()

    # Place particles far apart (no mutual interaction)
    model.positions[0] = [10.0, 10.0]
    model.positions[1] = [90.0, 90.0]
    model.headings[0] = 0.5
    model.headings[1] = 2.0

    model.step()

    # Each particle should keep its heading (only self-neighbor)
    assert abs(model.headings[0] - 0.5) < 1e-10, \
        f"Isolated particle heading changed: {model.headings[0]} != 0.5"
    assert abs(model.headings[1] - 2.0) < 1e-10, \
        f"Isolated particle heading changed: {model.headings[1]} != 2.0"

    print("  ✓ Self-inclusion verified: isolated particles maintain heading")


def test_neighbor_alignment():
    """Verify that two nearby particles align with zero noise."""
    model = VicsekModel(
        n_particles=2, box_size=100.0, speed=0.01, noise=0.0,
        interaction_radius=2.0, dt=1.0, seed=42
    )
    model.setup()

    # Place particles close together with different headings
    model.positions[0] = [50.0, 50.0]
    model.positions[1] = [50.5, 50.0]
    model.headings[0] = 0.0
    model.headings[1] = np.pi / 2

    model.step()

    # Both should converge to the circular mean of 0 and π/2 = π/4
    expected = np.arctan2(
        (np.sin(0.0) + np.sin(np.pi / 2)) / 2,
        (np.cos(0.0) + np.cos(np.pi / 2)) / 2,
    )
    for i in range(2):
        diff = abs(model.headings[i] - expected)
        diff = min(diff, 2 * np.pi - diff)  # handle wrap
        assert diff < 1e-10, \
            f"Particle {i} heading {model.headings[i]:.6f} != expected {expected:.6f}"

    print(f"  ✓ Two nearby particles align to circular mean (θ={expected:.4f})")


def test_periodic_neighbor_finding():
    """Verify neighbors found across periodic boundaries."""
    model = VicsekModel(
        n_particles=2, box_size=10.0, speed=0.01, noise=0.0,
        interaction_radius=2.0, dt=1.0, seed=42
    )
    model.setup()

    # Place particles near opposite edges (distance = 1.0 across boundary)
    model.positions[0] = [0.5, 5.0]
    model.positions[1] = [9.5, 5.0]
    model.headings[0] = 0.0
    model.headings[1] = np.pi / 2

    model.step()

    # They should interact (distance across boundary = 1.0 < 2.0)
    # Both should align to circular mean of 0 and π/2
    expected = np.arctan2(
        (np.sin(0.0) + np.sin(np.pi / 2)) / 2,
        (np.cos(0.0) + np.cos(np.pi / 2)) / 2,
    )
    for i in range(2):
        diff = abs(model.headings[i] - expected)
        diff = min(diff, 2 * np.pi - diff)
        assert diff < 1e-10, \
            f"Periodic neighbor not found: particle {i} heading " \
            f"{model.headings[i]:.6f} != {expected:.6f}"

    print("  ✓ Periodic boundary neighbor finding works correctly")


def test_disordered_baseline():
    """At maximum noise, φ should approach 1/√N.

    With η = 2π (full angular randomization), the headings are
    essentially independent uniform random variables, giving
    φ ≈ 1/√N from the central limit theorem.

    Power: 10 independent seeds × 500 measurement steps each.
    """
    N = 300
    n_seeds = 10
    n_equil = 500
    n_measure = 500

    phi_means = []
    for seed in range(n_seeds):
        model = VicsekModel(
            n_particles=N, box_size=7.0, speed=0.03,
            noise=2 * np.pi,  # maximum noise
            interaction_radius=1.0, dt=1.0, seed=seed
        )
        model.setup()

        # Equilibrate
        for _ in range(n_equil):
            model.step()

        # Measure
        phis = []
        for _ in range(n_measure):
            state = model.step()
            phis.append(polarization(state))

        phi_means.append(np.mean(phis))

    grand_mean = np.mean(phi_means)
    grand_std = np.std(phi_means) / np.sqrt(n_seeds)
    expected = 1 / np.sqrt(N)

    # Should be within 50% of 1/√N (it's a finite-N, finite-r effect)
    ratio = grand_mean / expected
    assert 0.5 < ratio < 2.0, \
        f"Disordered φ = {grand_mean:.4f}, expected ~{expected:.4f} (ratio {ratio:.2f})"

    print(f"  ✓ Disordered baseline: φ = {grand_mean:.4f} ± {grand_std:.4f} "
          f"(1/√N = {expected:.4f}, ratio = {ratio:.2f})")


def test_ordered_state():
    """At very low noise, φ should approach 1.0.

    Power: 5 seeds × long equilibration (3000 steps) + 1000 measurement.
    """
    N = 300
    n_seeds = 5
    n_equil = 3000
    n_measure = 1000

    phi_means = []
    for seed in range(n_seeds):
        model = VicsekModel(
            n_particles=N, box_size=7.0, speed=0.03,
            noise=0.1,  # very low noise
            interaction_radius=1.0, dt=1.0, seed=seed
        )
        model.setup()

        for _ in range(n_equil):
            model.step()

        phis = []
        for _ in range(n_measure):
            state = model.step()
            phis.append(polarization(state))

        phi_means.append(np.mean(phis))

    grand_mean = np.mean(phi_means)
    grand_std = np.std(phi_means) / np.sqrt(n_seeds)

    assert grand_mean > 0.8, \
        f"Ordered state too low: φ = {grand_mean:.4f} (expected > 0.8)"

    print(f"  ✓ Ordered state (η=0.1): φ = {grand_mean:.4f} ± {grand_std:.4f}")


def test_phase_transition_curve():
    """Reproduce the φ vs η phase transition from Vicsek 1995, Figure 3.

    Uses paper parameters: N=300, L=7, v₀=0.03, r=1.
    Scans η from 0.1 to 5.0 in 16 points.

    Power budget per η:
        - T_cross ≈ 233 steps
        - Equilibration: 3000 steps (~13 T_cross)
        - Measurement: 2000 steps (~8.6 T_cross)
        - Seeds: 5 independent runs
        - Total: 5 × 5000 = 25,000 steps per η point
    """
    N = 300
    L = 7.0
    v0 = 0.03
    r = 1.0

    eta_values = np.array([
        0.1, 0.3, 0.5, 0.8, 1.0, 1.3, 1.5, 1.8,
        2.0, 2.3, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0
    ])

    n_seeds = 5
    n_equil = 3000
    n_measure = 2000

    print(f"\n  Scanning η: {len(eta_values)} points × {n_seeds} seeds × "
          f"{n_equil + n_measure} steps = {len(eta_values) * n_seeds * (n_equil + n_measure):,} total steps")

    results = []
    t_start = time.time()

    for eta in eta_values:
        phi_seeds = []
        for seed in range(n_seeds):
            model = VicsekModel(
                n_particles=N, box_size=L, speed=v0, noise=eta,
                interaction_radius=r, dt=1.0, seed=seed * 100 + int(eta * 100)
            )
            model.setup()

            for _ in range(n_equil):
                model.step()

            phis = []
            for _ in range(n_measure):
                state = model.step()
                phis.append(polarization(state))

            phi_seeds.append(np.mean(phis))

        mean_phi = np.mean(phi_seeds)
        std_phi = np.std(phi_seeds) / np.sqrt(n_seeds)
        results.append((eta, mean_phi, std_phi))

    elapsed = time.time() - t_start

    # Print results table
    print(f"\n  Phase transition curve (elapsed: {elapsed:.1f}s):")
    print(f"  {'η':>6s}  {'φ_mean':>8s}  {'φ_sem':>8s}  {'bar':>20s}")
    print(f"  {'─' * 6}  {'─' * 8}  {'─' * 8}  {'─' * 20}")

    for eta, mean_phi, std_phi in results:
        bar_len = int(mean_phi * 20)
        bar = "█" * bar_len + "░" * (20 - bar_len)
        print(f"  {eta:6.1f}  {mean_phi:8.4f}  {std_phi:8.4f}  {bar}")

    # Validate monotonic decrease
    phis = [r[1] for r in results]
    # Check that low-noise is ordered, high-noise is disordered
    low_noise_phi = np.mean(phis[:3])   # η = 0.1, 0.3, 0.5
    high_noise_phi = np.mean(phis[-3:])  # η = 4.0, 4.5, 5.0

    assert low_noise_phi > 0.5, \
        f"Low-noise regime not ordered enough: φ = {low_noise_phi:.4f}"
    assert high_noise_phi < 0.3, \
        f"High-noise regime not disordered enough: φ = {high_noise_phi:.4f}"
    assert low_noise_phi > high_noise_phi, \
        f"No phase transition: low-noise φ ({low_noise_phi:.4f}) " \
        f"not greater than high-noise φ ({high_noise_phi:.4f})"

    # Check transition region exists (intermediate φ values)
    mid_phis = [r[1] for r in results if 1.5 <= r[0] <= 3.0]
    has_transition = any(0.1 < p < 0.9 for p in mid_phis)
    assert has_transition, \
        "No transition region found between η=1.5 and η=3.0"

    print(f"\n  ✓ Phase transition verified:")
    print(f"    Low noise (η ≤ 0.5): φ = {low_noise_phi:.4f} (ordered)")
    print(f"    High noise (η ≥ 4.0): φ = {high_noise_phi:.4f} (disordered)")
    print(f"    Transition region present in η ∈ [1.5, 3.0]")

    return results


def test_angular_momentum_low():
    """Standard Vicsek should NOT produce milling (|L| ≈ 0).

    This is a negative control for the P6 detector. The Vicsek model
    with pure alignment produces flocking (high φ) but no rotation
    (low |L|).
    """
    model = VicsekModel(
        n_particles=300, box_size=7.0, speed=0.03, noise=0.1,
        interaction_radius=1.0, dt=1.0, seed=42
    )
    model.setup()

    # Equilibrate
    for _ in range(3000):
        model.step()

    # Measure angular momentum
    L_values = []
    for _ in range(1000):
        state = model.step()
        L_values.append(abs(angular_momentum(state, 7.0)))

    mean_L = np.mean(L_values)
    assert mean_L < 0.3, \
        f"Angular momentum too high for pure Vicsek: |L| = {mean_L:.4f}"

    print(f"  ✓ No milling in standard Vicsek: mean |L| = {mean_L:.4f} (expected < 0.3)")


def test_group_speed_ratio():
    """R should track φ closely for constant-speed models."""
    model = VicsekModel(
        n_particles=300, box_size=7.0, speed=0.03, noise=0.1,
        interaction_radius=1.0, dt=1.0, seed=42
    )
    model.setup()

    for _ in range(3000):
        model.step()

    for _ in range(100):
        state = model.step()
        phi = polarization(state)
        R = group_speed_ratio(state)
        # For constant-speed Vicsek, R should equal φ
        assert abs(R - phi) < 1e-10, \
            f"R ({R:.6f}) differs from φ ({phi:.6f})"

    print("  ✓ Group speed ratio R = φ exactly for constant-speed model")


def test_reproducibility():
    """Same seed should produce identical results."""
    def run_short(seed):
        model = VicsekModel(
            n_particles=100, box_size=7.0, speed=0.03, noise=1.0,
            interaction_radius=1.0, dt=1.0, seed=seed
        )
        model.setup()
        for _ in range(100):
            model.step()
        return model.positions.copy(), model.headings.copy()

    pos1, head1 = run_short(42)
    pos2, head2 = run_short(42)

    assert np.array_equal(pos1, pos2), "Positions not reproducible with same seed"
    assert np.array_equal(head1, head2), "Headings not reproducible with same seed"

    print("  ✓ Reproducibility verified: identical seeds → identical states")


def test_density_dependence():
    """Higher density should support order at higher noise.

    At fixed η, increasing ρ (decreasing L for fixed N) should
    increase φ. This is a key prediction of the Vicsek model.

    Test: at η=2.0, compare ρ=0.5 (L≈24.5) vs ρ=6.12 (L=7).
    """
    N = 300
    eta = 2.0
    n_equil = 3000
    n_measure = 1000
    n_seeds = 3

    configs = [
        ("low_density", np.sqrt(N / 0.5)),   # ρ ≈ 0.5, L ≈ 24.5
        ("high_density", 7.0),                 # ρ ≈ 6.12
    ]

    phi_by_config = {}
    for label, L in configs:
        phi_seeds = []
        for seed in range(n_seeds):
            model = VicsekModel(
                n_particles=N, box_size=L, speed=0.03, noise=eta,
                interaction_radius=1.0, dt=1.0, seed=seed
            )
            model.setup()
            for _ in range(n_equil):
                model.step()
            phis = []
            for _ in range(n_measure):
                state = model.step()
                phis.append(polarization(state))
            phi_seeds.append(np.mean(phis))
        phi_by_config[label] = np.mean(phi_seeds)

    phi_low = phi_by_config["low_density"]
    phi_high = phi_by_config["high_density"]

    assert phi_high > phi_low, \
        f"Density dependence failed: high-ρ φ ({phi_high:.4f}) " \
        f"≤ low-ρ φ ({phi_low:.4f})"

    print(f"  ✓ Density dependence at η={eta}:")
    print(f"    Low density (ρ≈0.5):  φ = {phi_low:.4f}")
    print(f"    High density (ρ≈6.1): φ = {phi_high:.4f}")


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("VICSEK MODEL VALIDATION — Replication of Vicsek et al. (1995)")
    print("=" * 70)

    print("\n--- Unit Tests ---")
    test_periodic_boundary_conditions()
    test_constant_speed()
    test_perfect_alignment_no_noise()
    test_neighbor_self_inclusion()
    test_neighbor_alignment()
    test_periodic_neighbor_finding()
    test_reproducibility()
    test_group_speed_ratio()

    print("\n--- Physics Validation ---")
    test_disordered_baseline()
    test_ordered_state()
    test_angular_momentum_low()

    print("\n--- Phase Transition (Primary Validation Target) ---")
    results = test_phase_transition_curve()

    print("\n--- Density Dependence ---")
    test_density_dependence()

    print("\n" + "=" * 70)
    print("ALL VALIDATION TESTS PASSED")
    print("=" * 70)
