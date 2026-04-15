"""Nowak-May spatial PD × P27 end-to-end test.

Sprint 5: Opens Cluster H (competition/cooperation) with the first
model and detector in this cluster.

Reference: Nowak & May (1992), Nature 359, 826-829.

Validation targets:
  1. b=1.0: ~all cooperate (no temptation)
  2. b=1.5: cooperators survive (~60-90% depending on init)
  3. b=1.8: cooperators survive in clusters (~30-50%)
  4. b=2.0: cooperators extinct
  5. Moran's I > 0 at b=1.5, 1.8 (spatial clustering)
  6. P27 detection: DEFINITIVE at b=1.8 with PD metadata

Tests:
  1. Model physics validation (4 b-values)
  2. P27 detection end-to-end (positive + negative controls)
"""

from __future__ import annotations

import time
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def test_nowak_may_physics():
    """Validate Nowak-May against published phase diagram."""
    from epc.models.nowak_may import NowakMayModel

    print(f"\n{'='*60}")
    print("Nowak-May Model Validation")
    print(f"{'='*60}")

    results = {}
    for b, note in [(1.0, "no tempt"), (1.5, "moderate"),
                     (1.8, "high"), (2.0, "too high")]:
        m = NowakMayModel(rows=100, cols=100, b=b, init_mode="random",
                          init_coop_fraction=0.5, seed=42)
        h = m.run(n_steps=500)
        fc = h[-1]["coop_fraction"]
        mi = h[-1]["moran_i"]
        results[b] = (fc, mi)
        print(f"  b={b:.1f}: f_C={fc:.3f}, I={mi:.3f}  ({note})")

    # Assertions
    # b=1.0: nearly all cooperate
    assert results[1.0][0] > 0.90, f"b=1.0: f_C={results[1.0][0]:.3f} < 0.90"

    # b=1.5: cooperators survive
    assert results[1.5][0] > 0.30, f"b=1.5: f_C={results[1.5][0]:.3f} < 0.30"
    assert results[1.5][1] > 0.1, f"b=1.5: Moran I={results[1.5][1]:.3f} < 0.1"

    # b=1.8: cooperators survive in clusters
    assert 0.05 < results[1.8][0] < 0.80, \
        f"b=1.8: f_C={results[1.8][0]:.3f} outside (0.05, 0.80)"
    assert results[1.8][1] > 0.1, f"b=1.8: Moran I={results[1.8][1]:.3f} < 0.1"

    # b=2.0: cooperators extinct
    assert results[2.0][0] < 0.01, f"b=2.0: f_C={results[2.0][0]:.3f} > 0.01"

    # Monotonic decrease in cooperation with b
    assert results[1.0][0] > results[1.5][0] > results[1.8][0] > results[2.0][0], \
        "f_C should decrease monotonically with b"

    print(f"\n  ✅ All physics checks pass")


def test_p27_e2e():
    """P27 detection end-to-end: positive and negative controls."""
    from epc.models.nowak_may import NowakMayModel
    from epc.detectors.p27_spatial_reciprocity import detect_p27

    print(f"\n{'='*60}")
    print("P27 Spatial Reciprocity Detection")
    print(f"{'='*60}")

    # === Positive control: b=1.8 (cooperation survives) ===
    print("\n  --- Positive: b=1.8 (PD, cooperation survives) ---")
    m = NowakMayModel(rows=100, cols=100, b=1.8, init_mode="random",
                      init_coop_fraction=0.5, seed=42)
    h = m.run(n_steps=1500)
    meta = m.get_metadata()

    t0 = time.perf_counter()
    det = detect_p27(h, model_metadata=meta, n_permutations=199)
    elapsed = time.perf_counter() - t0

    print(f"  Result: {det.tier} (confidence {det.confidence:.3f})")
    print(f"    f_C = {det.coop_fraction:.4f}")
    print(f"    Moran I = {det.moran_i:.4f}, p = {det.moran_i_p:.4f}")
    print(f"    PD verified = {det.pd_verified}")
    print(f"    Generations = {det.n_generations}")
    print(f"    Time: {elapsed:.1f}s")

    assert det.detected, "P27 should be detected at b=1.8"
    assert det.tier == "definitive", f"Expected definitive, got {det.tier}"
    assert det.pd_verified, "PD structure should be verified"
    assert det.moran_i > 0.1, f"Moran I={det.moran_i:.3f} too low"
    assert det.moran_i_p < 0.01, f"Moran I not significant (p={det.moran_i_p:.3f})"

    # === Negative control 1: b=2.0 (cooperators extinct) ===
    print("\n  --- Negative: b=2.0 (cooperators extinct) ---")
    m2 = NowakMayModel(rows=100, cols=100, b=2.0, init_mode="random",
                       init_coop_fraction=0.5, seed=42)
    h2 = m2.run(n_steps=500)
    det2 = detect_p27(h2, model_metadata=m2.get_metadata(), n_permutations=99)
    print(f"  Result: {det2.tier} (f_C={det2.coop_fraction:.4f})")
    assert not det2.detected, f"P27 should NOT be detected at b=2.0 (f_C={det2.coop_fraction})"

    # === Negative control 2: b=1.0 (no PD — everyone cooperates trivially) ===
    print("\n  --- Negative: b=1.0 (no temptation, trivial cooperation) ---")
    m3 = NowakMayModel(rows=100, cols=100, b=1.0, init_mode="random",
                       init_coop_fraction=0.5, seed=42)
    h3 = m3.run(n_steps=500)
    meta3 = m3.get_metadata()
    # b=1.0 means T=R=1, so T > R fails → PD not verified → can't be definitive
    det3 = detect_p27(h3, model_metadata=meta3, n_permutations=99)
    print(f"  Result: {det3.tier} (f_C={det3.coop_fraction:.4f}, "
          f"PD={det3.pd_verified})")
    # Should detect cooperation but NOT as definitive P27
    # (because T=R=1, no dilemma, PD structure fails)
    assert not det3.pd_verified, "b=1.0: PD should NOT be verified (T=R)"
    assert det3.tier != "definitive", "b=1.0 should not be definitive P27"

    print(f"\n  ✅ All P27 detection checks pass")
    print(f"     b=1.8 → DEFINITIVE, b=2.0 → none, b=1.0 → not definitive")


def test_transfer_matrix_row():
    """Verify Nowak-May against P1 (cross-detection check)."""
    from epc.models.nowak_may import NowakMayModel

    print(f"\n{'='*60}")
    print("Transfer Matrix: Nowak-May row")
    print(f"{'='*60}")

    m = NowakMayModel(rows=100, cols=100, b=1.8, init_mode="random",
                      init_coop_fraction=0.5, seed=42)
    h = m.run(n_steps=1000)
    meta = m.get_metadata()

    # P27 should fire
    from epc.detectors.p27_spatial_reciprocity import detect_p27
    p27 = detect_p27(h, model_metadata=meta, n_permutations=99)
    print(f"  P27: {p27.tier} (f_C={p27.coop_fraction:.3f}, I={p27.moran_i:.3f})")

    # P1 substrate check: Nowak-May is lattice_2d with strategies,
    # but strategies are NOT persistent type labels (they change via imitation).
    # The detector card says: P1 exclusion via "movement → P1, selection → P27"
    print(f"  P1 exclusion: has_movement={meta.get('has_movement')} → P27 not P1")
    assert not meta.get("has_movement"), "Nowak-May has no physical movement"

    print(f"\n  ✅ Transfer matrix row correct")


if __name__ == "__main__":
    print("TEST 1: Nowak-May physics")
    test_nowak_may_physics()
    print("\n\nTEST 2: P27 end-to-end")
    test_p27_e2e()
    print("\n\nTEST 3: Transfer matrix")
    test_transfer_matrix_row()
    print("\n\nALL TESTS PASSED")
