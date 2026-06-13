"""P4 territoriality — movement-causality discrimination test.

Reproducible validation that the P4 detector distinguishes scent-mediated
territoriality from a SCENT-BLIND RANDOM WALK using a per-step
movement-causality metric, resolving the deferral finding that the legacy
cumulative-occupancy metrics (exclusivity / overlap / occupancy-scent
correlation) cannot tell them apart (a random walk scores equally or more
exclusive).

Metric: avoidance_ratio = foreign-scent at the CHOSEN cell / mean foreign-scent
over available neighbours, summed over moves. Territorial movement avoids
foreign scent (ratio << 1); a scent-blind walk is indifferent (~1). The two
runs share identical scent dynamics and deposition — only the MOVEMENT RULE
differs — so the metric isolates the causal dependence of movement on foreign
scent.

Run:  PYTHONPATH=. python3.12 analysis/validation_rebuild/p4_movement_causality_test.py
"""

from __future__ import annotations

from epc.models.territoriality import ScentMarkingModel, ScentMarkingParams
from epc.detectors.p4_territoriality import P4TerritorialityDetector


def bundle(seed: int, scent_blind: bool):
    p = ScentMarkingParams(
        n_agents=9, grid_size=48, deposition_rate=0.1, decay_rate=0.03,
        repulsion_strength=2.0, home_attraction=2.0, temperature=0.5,
        seed=42 + seed,
    )
    return ScentMarkingModel(p).simulate_movement_bundle(
        burnin=3000, window=1500, scent_blind=scent_blind)


def main() -> int:
    det = P4TerritorialityDetector(n_permutations=199, seed=42)
    failures = []

    print("=== TERRITORIAL positives (detected; DEFINITIVE on strong avoidance, CONFIRMATION on moderate) ===")
    for s in range(5):
        r = det.detect(bundle(s, scent_blind=False))
        ar = r.primary_metric.get("avoidance_ratio")
        print(f"  seed{s}: detected={r.detected!s:5s} tier={r.tier.value:11s} avoidance_ratio={ar:.3f} null_p={r.null_p_value:.4f}")
        if not (r.detected and r.tier.value in ("confirmation", "definitive")):
            failures.append(f"territorial seed{s} not detected as territorial (tier={r.tier.value})")

    print("=== SCENT-BLIND RANDOM WALK (must be rejected — same scent, no avoidance) ===")
    for s in range(5):
        r = det.detect(bundle(s, scent_blind=True))
        ar = r.primary_metric.get("avoidance_ratio")
        reason = r.warnings[0][:55] if r.warnings else ""
        print(f"  seed{s}: detected={r.detected!s:5s} tier={r.tier.value:11s} avoidance_ratio={ar:.3f}  {reason}")
        if r.detected:
            failures.append(f"random-walk seed{s} was DETECTED (false positive)")

    print()
    if failures:
        print("FAIL:")
        for f in failures:
            print("  -", f)
        return 1
    print("PASS: territorial movement avoids foreign scent (DEFINITIVE); a "
          "scent-blind random walk on identical scent dynamics is rejected.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
