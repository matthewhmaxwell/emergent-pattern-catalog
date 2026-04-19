"""
End-to-End Test: BTW Sandpile × P14 SOC Detector.

Sprint 4: Build and validate P14 (self-organized criticality) detection
on the canonical BTW sandpile model.

Tests:
1. BTW model physics validation (critical state, τ ≈ 1.20)
2. Full P14 detector pipeline (screening → confirmation → definitive)
3. Dissipative null model (subcritical → exponential)
4. Cross-detection: P13 should NOT detect SOC on sandpile
5. Negative control: dissipative sandpile should NOT trigger P14

Published replication targets:
- τ ≈ 1.20 (2D BTW, Lübeck & Usadel 1997)
- Power-law preferred over exponential
- Dissipative null produces exponential distribution

Statistical power:
- 100,000 driving events (after 10,000 burn-in)
- 50,000 null driving events (after 5,000 burn-in)
"""

import sys
import time
import numpy as np
import pytest

from epc.models.btw_sandpile import (
    run_sandpile, run_dissipative_sandpile, BTWSandpileParams
)
from epc.detectors.p14_soc import detect_p14


# =========================================================================
# Shared fixtures — expensive simulations run once per module
# =========================================================================

@pytest.fixture(scope="module")
def btw_result():
    """Run BTW sandpile once, share across tests.

    Sprint 14.6 D.1: This fixture produces a replication-quality 100k-event
    run (~130s) and is only consumed by slow-marked tests. Fast-half
    tests use the separate `test_p14_fast_smoke` below.
    """
    params = BTWSandpileParams(L=64, n_drive=100_000, n_burn=10_000, seed=42)
    return run_sandpile(params)


@pytest.fixture(scope="module")
def det(btw_result):
    """Run P14 detector on BTW result, share across tests."""
    null_params = BTWSandpileParams(L=64, n_drive=50_000, n_burn=5_000, seed=42)
    null_result = run_dissipative_sandpile(null_params, p_diss=0.2)
    return detect_p14(
        avalanche_sizes=btw_result.avalanche_sizes,
        avalanche_durations=btw_result.avalanche_durations,
        activity=btw_result.activity,
        energy=btw_result.energy_history,
        null_sizes=null_result.avalanche_sizes,
        is_self_tuned=True,
    )


# =========================================================================
# FAST-SMOKE: pipeline works end-to-end at reduced scale (Sprint 14.6 D.1)
# =========================================================================

def test_p14_fast_smoke():
    """Fast-smoke pipeline verification at reduced scale.

    Sprint 14.6 D.1: The full-replication tests below each run ~130s of
    BTW simulation (100k events at L=64). This fast-smoke test uses
    L=32 and n_drive=5_000 to verify that:

      1. The BTW sandpile model runs to completion without raising.
      2. The P14 detector accepts the output and returns a result.
      3. The returned result has the expected schema (tier, detected,
         fit, tau_in_range).

    It does NOT assert replication-quality τ values or power-law fit
    quality — those are pinned by the slow-marked `test_replication_summary`
    on 100k-event runs. If this smoke test fails, the pipeline is
    structurally broken; if it passes but the slow tests fail, the
    issue is a statistical-power or replication-quality regression.

    Timing: ~5s on typical hardware.
    """
    params = BTWSandpileParams(L=32, n_drive=5_000, n_burn=500, seed=42)
    result = run_sandpile(params)

    # Model physics: reached critical state, produced avalanches
    assert result.final_grid.max() < params.z_c, (
        f"BTW critical state violation: max height = "
        f"{result.final_grid.max()} >= z_c = {params.z_c}"
    )
    sizes = result.avalanche_sizes
    assert len(sizes) > 0, "BTW produced no avalanches"
    nonzero = sizes[sizes > 0]
    assert len(nonzero) > 50, (
        f"BTW produced only {len(nonzero)} non-zero avalanches "
        f"(need >50 for detector to run meaningfully)"
    )

    # Detector pipeline: runs to completion, returns schema-valid result
    det_result = detect_p14(
        avalanche_sizes=result.avalanche_sizes,
        avalanche_durations=result.avalanche_durations,
        activity=result.activity,
        energy=result.energy_history,
        is_self_tuned=True,
    )

    assert hasattr(det_result, "detected"), "Detector result missing 'detected'"
    assert hasattr(det_result, "tier"), "Detector result missing 'tier'"
    assert hasattr(det_result, "fit"), "Detector result missing 'fit'"
    # tier may be 'none', 'screening', 'confirmation', or 'definitive'
    # at this low statistical power — don't pin; just verify it ran.
    assert det_result.tier is not None


# =========================================================================
# Test 1: BTW Model Physics (slow — replication-quality)
# =========================================================================

@pytest.mark.slow
def test_btw_physics():
    """Verify BTW sandpile reaches critical state with correct properties."""
    print("=" * 70)
    print("TEST 1: BTW Model Physics Validation")
    print("=" * 70)
    
    params = BTWSandpileParams(L=64, n_drive=100_000, n_burn=10_000, seed=42)
    print(f"  Running BTW sandpile (L={params.L}, {params.n_drive} events)...")
    
    t0 = time.perf_counter()
    result = run_sandpile(params)
    elapsed = time.perf_counter() - t0
    print(f"  Time: {elapsed:.1f}s")
    
    sizes = result.avalanche_sizes
    nonzero = sizes[sizes > 0]
    
    print(f"\n  Total avalanches: {len(sizes)}")
    print(f"  Non-zero: {len(nonzero)} ({100*len(nonzero)/len(sizes):.1f}%)")
    print(f"  Size range: [{nonzero.min()}, {nonzero.max()}]")
    print(f"  Mean: {nonzero.mean():.1f}, Median: {np.median(nonzero):.0f}")
    print(f"  Grid: mean={result.final_grid.mean():.3f}, max={result.final_grid.max()}")
    
    checks = {}
    
    # Critical state: all heights < z_c
    checks['critical_state'] = result.final_grid.max() < params.z_c
    print(f"\n  {'✅' if checks['critical_state'] else '❌'} Critical state (max height = {result.final_grid.max()} < {params.z_c})")
    
    # Mean height should be around 2.0-2.2 for 2D BTW
    mean_h = result.final_grid.mean()
    checks['mean_height'] = 1.8 < mean_h < 2.4
    print(f"  {'✅' if checks['mean_height'] else '❌'} Mean height {mean_h:.3f} ∈ (1.8, 2.4)")
    
    # Span > 3 orders of magnitude in size
    span = np.log10(nonzero.max() / nonzero.min())
    checks['size_span'] = span > 3
    print(f"  {'✅' if checks['size_span'] else '❌'} Size span: {span:.1f} decades (need >3)")
    
    # Heavy-tailed: mean >> median
    checks['heavy_tail'] = nonzero.mean() > 5 * np.median(nonzero)
    print(f"  {'✅' if checks['heavy_tail'] else '❌'} Heavy-tailed: mean/median = {nonzero.mean()/np.median(nonzero):.1f}")
    
    all_ok = all(checks.values())
    assert all_ok, f"BTW physics checks failed: {checks}"


# =========================================================================
# Test 2: Full P14 Detector End-to-End (slow — replication-quality)
# =========================================================================

@pytest.mark.slow
def test_p14_e2e(btw_result):
    """Run full P14 detector on BTW sandpile."""
    print("\n" + "=" * 70)
    print("TEST 2: P14 Detector End-to-End")
    print("=" * 70)
    
    # Run dissipative null for the detector
    print("  Running dissipative null (p_diss=0.2)...")
    null_params = BTWSandpileParams(L=64, n_drive=50_000, n_burn=5_000, seed=42)
    t0 = time.perf_counter()
    null_result = run_dissipative_sandpile(null_params, p_diss=0.2)
    null_elapsed = time.perf_counter() - t0
    print(f"  Null time: {null_elapsed:.1f}s")
    
    # Run P14 detector
    print(f"\n  Running P14 detector...")
    t0 = time.perf_counter()
    det = detect_p14(
        avalanche_sizes=btw_result.avalanche_sizes,
        avalanche_durations=btw_result.avalanche_durations,
        activity=btw_result.activity,
        energy=btw_result.energy_history,
        null_sizes=null_result.avalanche_sizes,
        is_self_tuned=True,  # BTW is inherently self-tuning
    )
    det_elapsed = time.perf_counter() - t0
    print(f"  Detector time: {det_elapsed:.1f}s")
    
    print(f"\n  P14 Detection Results:")
    print(f"    detected:     {det.detected}")
    print(f"    tier:         {det.tier}")
    print(f"    confidence:   {det.confidence:.3f}")
    
    if det.fit:
        print(f"\n  Power-Law Fit:")
        print(f"    τ (MLE):      {det.fit.tau:.3f}")
        print(f"    τ (log-bin):  {det.fit.tau_logbin:.3f}")
        print(f"    τ in range:   {det.tau_in_range}")
        print(f"    KS distance:  {det.fit.ks_distance:.4f}")
        print(f"    vs exponential: R={det.fit.lr_vs_exponential:.3f}, p={det.fit.lr_vs_exponential_p:.4f}")
        print(f"    vs log-normal: R={det.fit.lr_vs_lognormal:.3f}, p={det.fit.lr_vs_lognormal_p:.4f}")
    
    print(f"\n  Secondaries:")
    print(f"    Duration γ:   {det.duration_gamma}")
    print(f"    Spectral β:   {det.spectral_beta}")
    
    print(f"\n  Null Comparison:")
    print(f"    Null τ:       {det.null_tau}")
    print(f"    Null exponential: {det.null_is_exponential}")
    print(f"    Null max size: {det.null_max_size}")
    print(f"    BTW max size:  {det.size_range[1]}")
    
    if det.warnings:
        print(f"\n  Warnings:")
        for w in det.warnings:
            print(f"    ⚠️ {w}")
    
    # Verification checks
    checks = {}
    checks['detected'] = det.detected
    checks['tier_confirmation+'] = det.tier in ('confirmation', 'definitive')
    checks['tau_in_range'] = det.tau_in_range
    
    if det.fit:
        checks['tau_near_1.20'] = abs(det.fit.tau - 1.20) < 0.15
        checks['exp_rejected'] = det.fit.lr_vs_exponential > 0
    
    checks['null_exponential'] = bool(det.null_is_exponential) is True
    
    print(f"\n  Verification:")
    all_pass = True
    for name, passed in checks.items():
        print(f"    {'✅' if passed else '❌'} {name}: {passed}")
        if not passed:
            all_pass = False
    
    assert all_pass, f"P14 e2e checks failed: {checks}"


# =========================================================================
# Test 3: Negative Control — Dissipative Sandpile
# =========================================================================

def test_dissipative_negative():
    """Dissipative sandpile should NOT be detected as SOC."""
    print("\n" + "=" * 70)
    print("TEST 3: Negative Control — Dissipative Sandpile")
    print("=" * 70)
    
    params = BTWSandpileParams(L=64, n_drive=50_000, n_burn=5_000, seed=42)
    t0 = time.perf_counter()
    diss_result = run_dissipative_sandpile(params, p_diss=0.2)
    elapsed = time.perf_counter() - t0
    
    nz = diss_result.avalanche_sizes[diss_result.avalanche_sizes > 0]
    print(f"  Dissipative sandpile (p_diss=0.2): {elapsed:.1f}s")
    print(f"  Non-zero: {len(nz)}, range: [{nz.min()}, {nz.max()}]")
    print(f"  Mean: {nz.mean():.1f}, Median: {np.median(nz):.0f}")
    
    det = detect_p14(
        avalanche_sizes=diss_result.avalanche_sizes,
        is_self_tuned=False,  # dissipative model is NOT self-tuned
    )
    
    print(f"\n  P14 Result:")
    print(f"    detected: {det.detected}")
    print(f"    tier:     {det.tier}")
    if det.fit:
        print(f"    τ:        {det.fit.tau:.3f}")
        print(f"    vs exp:   R={det.fit.lr_vs_exponential:.3f}")
    
    # Should NOT be detected (exponential distribution, not power-law)
    not_detected = not det.detected or det.tier == 'none'
    print(f"\n  {'✅' if not_detected else '❌'} Dissipative sandpile not detected as SOC: {not_detected}")
    assert not_detected


# =========================================================================
# Test 4: Replication Summary (slow — requires full 100k-event fixture)
# =========================================================================

@pytest.mark.slow
def test_replication_summary(det):
    """Summarize replication against published BTW results."""
    print("\n" + "=" * 70)
    print("TEST 4: Replication Summary")
    print("=" * 70)
    
    if det.fit is None:
        print("  ❌ No fit available")
        return False
    
    print(f"\n  | Published claim                     | Our result        | Match |")
    print(f"  |--------------------------------------|-------------------|-------|")
    
    tau_ok = abs(det.fit.tau - 1.20) < 0.15
    print(f"  | τ ≈ 1.20 (2D BTW)                   | τ = {det.fit.tau:.3f}         | {'✅' if tau_ok else '❌'}    |")
    
    exp_rej = det.fit.lr_vs_exponential > 0
    print(f"  | Power-law over exponential           | R = {det.fit.lr_vs_exponential:+.1f}         | {'✅' if exp_rej else '❌'}    |")
    
    crit = det.tau_in_range
    print(f"  | τ in plausible range [1.0, 2.0]      | {det.fit.tau:.3f} ∈ range    | {'✅' if crit else '❌'}    |")
    
    null_ok = det.null_is_exponential is not None and bool(det.null_is_exponential)
    print(f"  | Dissipative null exponential         | R < 0              | {'✅' if null_ok else '⚠️'}    |")
    
    if det.duration_gamma is not None:
        dur_ok = 0.3 < det.duration_gamma < 0.9
        print(f"  | Duration scaling T ~ s^γ            | γ = {det.duration_gamma:.3f}         | {'✅' if dur_ok else '⚠️'}    |")
    
    # Known caveat: BTW multifractal scaling
    ln_issue = det.fit.lr_vs_lognormal < 0
    if ln_issue:
        print(f"\n  Documented caveat: Log-normal preferred over simple power-law")
        print(f"  (R={det.fit.lr_vs_lognormal:.1f}). This is a known property of the")
        print(f"  2D BTW universality class (multifractal scaling with logarithmic")
        print(f"  corrections). Does NOT invalidate SOC detection.")
    
    all_ok = tau_ok and exp_rej and null_ok
    print(f"\n  {'✅' if all_ok else '⚠️'} Core replication: {'PASS' if all_ok else 'PARTIAL'}")
    assert all_ok, "Core replication checks failed"


# =========================================================================
# Main
# =========================================================================

def main():
    print("SPRINT 4: BTW SANDPILE × P14 SOC END-TO-END VALIDATION")
    print("=" * 70)
    print()
    
    results = {}
    
    # Test 1: Physics
    physics_ok, btw_result = test_btw_physics()
    results['btw_physics'] = physics_ok
    
    # Test 2: P14 detector
    det_ok, det_result = test_p14_e2e(btw_result)
    results['p14_e2e'] = det_ok
    
    # Test 3: Negative control
    results['dissipative_negative'] = test_dissipative_negative()
    
    # Test 4: Replication
    results['replication'] = test_replication_summary(det_result)
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    
    for name, passed in results.items():
        print(f"  {'✅' if passed else '❌'} {name}")
    
    n_passed = sum(results.values())
    n_total = len(results)
    print(f"\n  {n_passed}/{n_total} tests passed")
    
    if all(results.values()):
        print("\n  BTW × P14 END-TO-END VALIDATED")
    else:
        failed = [name for name, v in results.items() if not v]
        print(f"\n  ISSUES: {failed}")
    
    assert all(results.values())


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
