"""Hegselmann-Krause × P21 end-to-end test.

Sprint 5: Opens Cluster F (Decision-Making/Social Dynamics). 14th model, 11th detector.

Reference: Hegselmann & Krause (2002), JASSS 5(3).

Tests:
  1. Model physics: consensus/polarization/fragmentation phase diagram
  2. P21 detection: positive (ε=0.2) + negative controls (ε=0.5, ε→0)
"""

from __future__ import annotations

import time
import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def test_hk_physics():
    """Validate HK against published phase diagram."""
    from epc.models.hegselmann_krause import HegselmannKrauseModel

    print(f"\n{'='*60}")
    print("Hegselmann-Krause Model Validation")
    print(f"{'='*60}")

    configs = [
        (0.50, "consensus", 1, 1),
        (0.30, "few clusters", 1, 3),
        (0.20, "polarization", 2, 5),
        (0.10, "fragmentation", 3, 8),
        (0.05, "many clusters", 5, 15),
    ]

    results = {}
    for eps, note, min_cl, max_cl in configs:
        m = HegselmannKrauseModel(n_agents=500, epsilon=eps, init_mode="uniform", seed=42)
        h = m.run(n_steps=500)
        f = h[-1]
        results[eps] = f
        print(f"  ε={eps:.2f}: {f['n_clusters']:2d} clusters, "
              f"var={f['variance']:.4f}, steps={f['step']}, conv={f['converged']} ({note})")
        assert min_cl <= f["n_clusters"] <= max_cl, \
            f"ε={eps}: {f['n_clusters']} clusters outside [{min_cl}, {max_cl}]"

    # Cluster count should increase as epsilon decreases
    epsilons = [0.50, 0.30, 0.20, 0.10, 0.05]
    clusters = [results[e]["n_clusters"] for e in epsilons]
    assert all(clusters[i] <= clusters[i + 1] for i in range(len(clusters) - 1)), \
        f"Cluster count should be non-decreasing: {clusters}"

    # All should converge
    for eps in epsilons:
        assert results[eps]["converged"], f"ε={eps} did not converge"

    print(f"\n  ✅ All physics checks pass")


def test_p21_e2e():
    """P21 detection: positive and negative controls."""
    from epc.models.hegselmann_krause import HegselmannKrauseModel
    from epc.detectors.p21_polarization import detect_p21

    print(f"\n{'='*60}")
    print("P21 Polarization Detection")
    print(f"{'='*60}")

    # === Positive: ε=0.2 (polarization, 2-3 clusters from uniform IC) ===
    print("\n  --- Positive: ε=0.2 (polarization) ---")
    m = HegselmannKrauseModel(n_agents=500, epsilon=0.2, init_mode="uniform", seed=42)
    h = m.run(n_steps=500)
    meta = m.get_metadata()

    # Pad history to meet persistence requirement (model converges fast)
    # After convergence, state is frozen, so repeating final state is valid
    while len(h) < 200:
        h.append(h[-1].copy())

    t0 = time.perf_counter()
    det = detect_p21(h, model_metadata=meta, n_boot=1000)
    elapsed = time.perf_counter() - t0

    print(f"  Result: {det.tier} (confidence {det.confidence:.3f})")
    print(f"    Clusters: {det.n_clusters}")
    print(f"    Dip: stat={det.dip_stat:.4f}, p={det.dip_p:.4f}")
    print(f"    Variance: {det.variance:.4f}")
    print(f"    Persistence: {det.persistence_steps} steps")
    print(f"    From unimodal: {det.from_unimodal}")
    print(f"    Time: {elapsed:.1f}s")
    if det.warnings:
        for w in det.warnings:
            print(f"    ⚠️ {w}")

    assert det.detected, "P21 should be detected at ε=0.2"
    assert det.n_clusters >= 2, f"Expected ≥2 clusters, got {det.n_clusters}"
    assert det.dip_p < 0.01, f"Dip test not significant: p={det.dip_p}"
    assert det.from_unimodal, "Should recognize uniform IC as unimodal"

    # === Negative: ε=0.5 (consensus, 1 cluster) ===
    print("\n  --- Negative: ε=0.5 (consensus) ---")
    m2 = HegselmannKrauseModel(n_agents=500, epsilon=0.5, init_mode="uniform", seed=42)
    h2 = m2.run(n_steps=300)
    det2 = detect_p21(h2, model_metadata=m2.get_metadata(), n_boot=500)
    print(f"  Result: {det2.tier} (clusters={det2.n_clusters}, dip_p={det2.dip_p:.3f})")
    assert not det2.detected or det2.tier == "none", \
        f"P21 should NOT be detected at ε=0.5 (consensus)"

    # === Positive: ε=0.1 (fragmentation, many clusters) ===
    print("\n  --- Positive: ε=0.1 (fragmentation) ---")
    m3 = HegselmannKrauseModel(n_agents=500, epsilon=0.1, init_mode="uniform", seed=42)
    h3 = m3.run(n_steps=500)
    while len(h3) < 200:
        h3.append(h3[-1].copy())
    det3 = detect_p21(h3, model_metadata=m3.get_metadata(), n_boot=1000)
    print(f"  Result: {det3.tier} (clusters={det3.n_clusters}, dip_p={det3.dip_p:.4f})")
    assert det3.detected, "P21 should be detected at ε=0.1"
    assert det3.n_clusters >= 3, f"Expected ≥3 clusters at ε=0.1, got {det3.n_clusters}"

    print(f"\n  ✅ All P21 detection checks pass")
    print(f"     ε=0.2 → {det.tier}, ε=0.5 → {det2.tier}, ε=0.1 → {det3.tier}")


if __name__ == "__main__":
    print("TEST 1: HK physics")
    test_hk_physics()
    print("\n\nTEST 2: P21 end-to-end")
    test_p21_e2e()
    print("\n\nALL TESTS PASSED")
