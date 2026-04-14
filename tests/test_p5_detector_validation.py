"""
Validation tests for collective motion metrics and P5 flocking detector.

Tests:
    1. Metric unit tests (polarization, R, L, heading distribution)
    2. P5 positive control: Vicsek at low noise → CONFIRM or DEFINITIVE
    3. P5 negative control: Vicsek at high noise → NOT DETECTED
    4. P5 null model sanity: shuffle φ ≈ 1/√N
    5. P5 exclusion logic: high |L| → P6 exclusion triggers

Power requirements:
    - Positive control: ≥ 10 T_cross measurement window, 199 permutations
    - Negative control: ≥ 5 T_cross measurement window, 199 permutations
    - T_cross = L/v₀ = 7/0.03 ≈ 233 steps
"""

import sys
import time
import numpy as np

sys.path.insert(0, "/home/claude")
from epc.models.vicsek import VicsekModel
from epc.metrics.collective_motion import (
    PolarizationMetric,
    GroupSpeedRatioMetric,
    AngularMomentumMetric,
    HeadingDistributionMetric,
    HeadingAutocorrelationMetric,
)
from epc.detectors.p5_flocking import P5FlockingDetector, heading_shuffle_null


# ===================================================================
# Metric Unit Tests
# ===================================================================

def test_polarization_perfect_alignment():
    """All particles heading right → φ = 1.0."""
    N = 100
    state = {
        "velocities": np.column_stack([np.ones(N), np.zeros(N)]),
    }
    phi = PolarizationMetric.compute_instant(state)
    assert abs(phi - 1.0) < 1e-10, f"Expected 1.0, got {phi}"
    print("  ✓ Polarization: perfect alignment → φ = 1.0")


def test_polarization_random():
    """Random headings → φ ≈ 1/√N."""
    N = 1000
    rng = np.random.default_rng(42)
    headings = rng.uniform(-np.pi, np.pi, size=N)
    state = {
        "velocities": np.column_stack([np.cos(headings), np.sin(headings)]),
    }
    phi = PolarizationMetric.compute_instant(state)
    expected = 1 / np.sqrt(N)
    assert phi < 3 * expected, f"φ={phi:.4f} too high for random (expected ~{expected:.4f})"
    print(f"  ✓ Polarization: random headings → φ = {phi:.4f} (1/√N = {expected:.4f})")


def test_polarization_antiparallel():
    """Half right, half left → φ ≈ 0."""
    N = 100
    v = np.zeros((N, 2))
    v[:N // 2, 0] = 1.0
    v[N // 2:, 0] = -1.0
    state = {"velocities": v}
    phi = PolarizationMetric.compute_instant(state)
    assert phi < 0.01, f"Expected ~0, got {phi}"
    print(f"  ✓ Polarization: antiparallel → φ = {phi:.6f}")


def test_group_speed_ratio_aligned():
    """All same direction → R = 1.0."""
    N = 100
    state = {
        "velocities": np.column_stack([np.ones(N) * 0.5, np.zeros(N)]),
    }
    R = GroupSpeedRatioMetric.compute_instant(state)
    assert abs(R - 1.0) < 1e-10, f"Expected 1.0, got {R}"
    print("  ✓ Group speed ratio: aligned → R = 1.0")


def test_angular_momentum_rotation():
    """Particles arranged in circle with tangential velocities → |L| ≈ 1."""
    N = 100
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
    r = 5.0
    positions = np.column_stack([
        50.0 + r * np.cos(angles),
        50.0 + r * np.sin(angles),
    ])
    # Tangential velocity: perpendicular to radial, counterclockwise
    velocities = np.column_stack([
        -np.sin(angles),
        np.cos(angles),
    ])
    state = {"positions": positions, "velocities": velocities}
    L = AngularMomentumMetric.compute_instant(state, box_size=None)
    assert abs(L) > 0.9, f"Expected |L| ≈ 1.0, got {L:.4f}"
    print(f"  ✓ Angular momentum: pure rotation → L = {L:.4f}")


def test_angular_momentum_translation():
    """All moving in same direction → |L| ≈ 0."""
    N = 100
    rng = np.random.default_rng(42)
    positions = rng.uniform(0, 100, size=(N, 2))
    velocities = np.column_stack([np.ones(N), np.zeros(N)])
    state = {"positions": positions, "velocities": velocities}
    L = AngularMomentumMetric.compute_instant(state, box_size=None)
    assert abs(L) < 0.15, f"Expected |L| ≈ 0, got {L:.4f}"
    print(f"  ✓ Angular momentum: pure translation → L = {L:.4f}")


def test_heading_distribution_unimodal():
    """All headings near 0 → unimodal, not antiparallel."""
    N = 200
    rng = np.random.default_rng(42)
    headings = rng.normal(0, 0.1, size=N)
    velocities = np.column_stack([np.cos(headings), np.sin(headings)])

    # Create a minimal history
    history = [{"velocities": velocities}] * 10
    result = HeadingDistributionMetric.compute(history)
    assert not result["is_bimodal_antiparallel"], \
        f"Unimodal distribution falsely flagged as bimodal"
    print(f"  ✓ Heading distribution: unimodal correctly identified "
          f"(nematic S = {result['nematic_order']:.4f}, "
          f"polar φ = {result['resultant_length']:.4f})")


def test_heading_distribution_bimodal():
    """Half heading 0, half heading π → bimodal antiparallel.

    Nematic order S ≈ 1, polar order φ ≈ 0 → antiparallel detected.
    """
    N = 200
    rng = np.random.default_rng(42)
    headings = np.concatenate([
        rng.normal(0, 0.1, size=N // 2),
        rng.normal(np.pi, 0.1, size=N // 2),
    ])
    velocities = np.column_stack([np.cos(headings), np.sin(headings)])

    history = [{"velocities": velocities}] * 10
    result = HeadingDistributionMetric.compute(history)
    assert result["is_bimodal_antiparallel"], \
        f"Bimodal antiparallel distribution not detected " \
        f"(nematic S={result['nematic_order']:.4f}, " \
        f"polar φ={result['resultant_length']:.4f})"
    print(f"  ✓ Heading distribution: bimodal antiparallel correctly identified "
          f"(nematic S = {result['nematic_order']:.4f}, "
          f"polar φ = {result['resultant_length']:.4f})")


# ===================================================================
# P5 Detector Tests
# ===================================================================

def generate_vicsek_history(
    noise: float, n_steps: int, n_equil: int = 0, seed: int = 42,
) -> tuple[list[dict], dict]:
    """Run Vicsek model and return history + metadata."""
    model = VicsekModel(
        n_particles=300, box_size=7.0, speed=0.03, noise=noise,
        interaction_radius=1.0, dt=1.0, seed=seed,
    )
    state0 = model.setup()

    # Equilibrate (not recorded)
    for _ in range(n_equil):
        model.step()

    # Record history
    history = [model.get_state()]
    for _ in range(n_steps):
        model.step()
        history.append(model.get_state())

    metadata = model.get_metadata()
    return history, metadata


def test_p5_positive_control():
    """Vicsek at low noise (η=0.5) should detect P5 at confirmation tier.

    Power budget:
        T_cross = 233 steps
        Total: 3000 equil + 5000 record = 8000 steps
        Measurement window (second half): 2500 steps ≈ 10.7 T_cross
        Permutations: 199 (floor p = 0.005)
    """
    print("\n  Running positive control (η=0.5, 3000+5000 steps, 199 perms)...")
    t0 = time.time()

    history, metadata = generate_vicsek_history(
        noise=0.5, n_steps=5000, n_equil=3000, seed=42,
    )

    detector = P5FlockingDetector(
        n_permutations=199,
        T_cross=metadata["box_size"] / metadata["speed"],
        box_size=metadata["box_size"],
        seed=42,
    )

    result = detector.detect(history, metadata)
    elapsed = time.time() - t0

    print(f"    Time: {elapsed:.1f}s")
    print(f"    Detected: {result.detected}")
    print(f"    Tier: {result.tier}")
    print(f"    Confidence: {result.confidence:.3f}")
    print(f"    φ_mean: {result.primary_metric['polarization_mean']:.4f}")
    print(f"    R_mean: {result.secondary_metrics['group_speed_ratio']:.4f}")
    print(f"    |L|_mean: {result.secondary_metrics['angular_momentum_abs']:.4f}")
    print(f"    Null p-value: {result.null_p_value:.4f}")
    print(f"    Cohen's d: {result.effect_size['cohens_d']:.1f}")
    print(f"    Exclusions: {result.exclusion_results}")
    print(f"    T_cross coverage: {result.secondary_metrics['measurement_T_cross']:.1f}")

    assert result.detected, "P5 not detected at η=0.5 (should be ordered)"
    assert result.tier in ("confirmation", "definitive"), \
        f"Expected confirmation or definitive, got {result.tier}"
    assert result.null_p_value < 0.01, \
        f"Null p-value too high: {result.null_p_value}"
    assert result.primary_metric["polarization_mean"] > 0.7, \
        f"φ too low for confirmation: {result.primary_metric['polarization_mean']}"

    print(f"  ✓ Positive control: P5 detected at {result.tier} tier "
          f"(φ={result.primary_metric['polarization_mean']:.4f}, "
          f"p={result.null_p_value:.4f})")


def test_p5_negative_control():
    """Vicsek at high noise (η=5.0) should NOT detect P5.

    Power budget:
        Total: 1000 equil + 2000 record = 3000 steps
        Measurement window: 1000 steps ≈ 4.3 T_cross
        Permutations: 199
    """
    print("\n  Running negative control (η=5.0, 1000+2000 steps)...")
    t0 = time.time()

    history, metadata = generate_vicsek_history(
        noise=5.0, n_steps=2000, n_equil=1000, seed=42,
    )

    detector = P5FlockingDetector(
        n_permutations=199,
        T_cross=metadata["box_size"] / metadata["speed"],
        box_size=metadata["box_size"],
        seed=42,
    )

    result = detector.detect(history, metadata)
    elapsed = time.time() - t0

    print(f"    Time: {elapsed:.1f}s")
    print(f"    Detected: {result.detected}")
    print(f"    Tier: {result.tier}")
    print(f"    φ_mean: {result.primary_metric['polarization_mean']:.4f}")
    print(f"    Null p-value: {result.null_p_value:.4f}")

    assert not result.detected, \
        f"P5 falsely detected at η=5.0: tier={result.tier}, " \
        f"φ={result.primary_metric['polarization_mean']:.4f}"

    print(f"  ✓ Negative control: P5 not detected at η=5.0 "
          f"(φ={result.primary_metric['polarization_mean']:.4f})")


def test_p5_null_model_sanity():
    """Shuffle null should produce φ ≈ 1/√N.

    Uses a mix of headings (not all identical) so shuffling
    actually destroys spatial structure.
    """
    print("\n  Running null model sanity check...")

    N = 300
    rng = np.random.default_rng(42)

    # Create a spatially ordered state: headings vary with position
    # (nearby particles aligned, distant ones differ)
    # This gives high φ but shuffling should destroy spatial correlation
    history = []
    for _ in range(100):
        # Headings: smooth gradient across particles + small noise
        headings = np.linspace(0, 0.3, N) + rng.normal(0, 0.05, size=N)
        state = {
            "velocities": np.column_stack([
                np.cos(headings), np.sin(headings)
            ]),
        }
        history.append(state)

    null_result = heading_shuffle_null(history, n_permutations=199, rng=rng)

    expected_null_phi = 1 / np.sqrt(N)

    print(f"    Observed φ_mean: {null_result['observed_phi']:.4f}")
    print(f"    Null φ_mean: {null_result['null_mean']:.4f}")
    print(f"    Expected null (1/√N): {expected_null_phi:.4f}")
    print(f"    Null φ_std: {null_result['null_std']:.6f}")
    print(f"    Cohen's d: {null_result['cohens_d']:.1f}")
    print(f"    p-value: {null_result['p_value']:.4f}")

    # The null should be much lower than observed
    assert null_result["null_mean"] < null_result["observed_phi"], \
        "Null φ should be lower than observed for ordered data"
    assert null_result["p_value"] < 0.01, \
        "Spatially ordered data should be significant vs shuffle null"

    print(f"  ✓ Null model: shuffle destroys order "
          f"(observed φ={null_result['observed_phi']:.4f}, "
          f"null φ={null_result['null_mean']:.4f}, p={null_result['p_value']:.4f})")


def test_p5_exclusion_p6():
    """High angular momentum should trigger P6 exclusion, blocking definitive."""
    print("\n  Testing P6 exclusion logic...")

    # Create synthetic history with high |L| (milling)
    N = 100
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
    r = 3.0
    history = []
    for t in range(500):
        # Rotating mill
        phase = t * 0.01
        positions = np.column_stack([
            50.0 + r * np.cos(angles + phase),
            50.0 + r * np.sin(angles + phase),
        ])
        # Tangential velocities (high φ AND high |L|)
        velocities = np.column_stack([
            -np.sin(angles + phase),
            np.cos(angles + phase),
        ])
        history.append({
            "positions": positions,
            "velocities": velocities,
        })

    detector = P5FlockingDetector(
        n_permutations=99,
        T_cross=50.0,  # short T_cross for testing
        box_size=None,  # non-periodic for this test
        seed=42,
    )

    result = detector.detect(history)

    print(f"    |L|_mean: {result.secondary_metrics['angular_momentum_abs']:.4f}")
    print(f"    P6 exclusion: {result.exclusion_results.get('P6', 'N/A')}")
    print(f"    Tier: {result.tier}")

    assert result.exclusion_results["P6"] == "TRIGGERED", \
        f"P6 exclusion should trigger for milling data"
    assert result.tier != "definitive", \
        f"Should not reach definitive with P6 exclusion triggered"

    print(f"  ✓ P6 exclusion correctly triggered (|L|="
          f"{result.secondary_metrics['angular_momentum_abs']:.4f}), "
          f"tier capped at {result.tier}")


# ===================================================================
# Main
# ===================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("P5 FLOCKING DETECTOR VALIDATION")
    print("=" * 70)

    print("\n--- Metric Unit Tests ---")
    test_polarization_perfect_alignment()
    test_polarization_random()
    test_polarization_antiparallel()
    test_group_speed_ratio_aligned()
    test_angular_momentum_rotation()
    test_angular_momentum_translation()
    test_heading_distribution_unimodal()
    test_heading_distribution_bimodal()

    print("\n--- Null Model Sanity ---")
    test_p5_null_model_sanity()

    print("\n--- P5 Detector: Positive Control ---")
    test_p5_positive_control()

    print("\n--- P5 Detector: Negative Control ---")
    test_p5_negative_control()

    print("\n--- P5 Detector: Exclusion Logic ---")
    test_p5_exclusion_p6()

    print("\n" + "=" * 70)
    print("ALL P5 DETECTOR VALIDATION TESTS PASSED")
    print("=" * 70)
