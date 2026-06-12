"""P28 wealth-condensation — mechanism-from-trajectory discrimination test.

Reproducible validation that the P28 detector's tier is determined by the
TRAJECTORY (conserved resource + condensation emerging from near-equality +
monotonic super-Boltzmann growth), NOT by self-reported metadata booleans.

This is the honest panel-equivalent for P28 pending the full Phase-2a harness
wiring (which needs a ``scalar_wealth`` Class A synthetic format across the ten
generators — tracked separately). It directly exercises the audit's three
complaints:

  1. "mechanism gate reads 4 self-reported metadata booleans; flipping one flips
     CONFIRMATION->DEFINITIVE"  -> metadata-flip invariance check below.
  2. "add a Bouchaud-Mezard lookalike"  -> non-conserved multiplicative-growth
     lookalike, rejected by the conservation gate.
  3. "build a real panel"  -> positive vs 4 conserved/non-conserved lookalikes.

Run:  PYTHONPATH=. python3.12 analysis/validation_rebuild/p28_mechanism_discrimination.py
"""

from __future__ import annotations

import numpy as np

from epc.models.yard_sale import YardSale
from epc.detectors.p28_wealth_condensation import P28WealthCondensationDetector
from epc.phase2a.failed_regimes import p28_wealth as fr


def positive_history(seed: int, n_frames: int = 100, n_tx: int = 2000):
    """Pure yard-sale condensation (conserved, multiplicative), t=0 recorded."""
    m = YardSale(n_agents=1000, f=0.5, lambda_save=0.0, chi=0.0,
                 init_mode="equal", seed=seed)
    s0 = m.setup()
    frames = [{"wealth": np.asarray(s0["wealth"], dtype=np.float64).copy()}]
    for _ in range(n_frames):
        s = m.step(n_transactions=n_tx)
        frames.append({"wealth": np.asarray(s["wealth"], dtype=np.float64).copy()})
    return frames, m.get_metadata()


def main() -> int:
    det = P28WealthCondensationDetector()
    failures = []

    # 1) Positives reach DEFINITIVE
    for sd in range(1, 6):
        h, meta = positive_history(sd)
        r = det.detect(h, model_metadata=meta)
        ok = r.tier.value == "definitive"
        print(f"POSITIVE seed{sd}: tier={r.tier.value} gini={r.primary_metric['gini_final']:.3f} "
              f"g0={r.primary_metric['gini_initial']:.3f} drift={r.primary_metric['conservation_drift']:.4f}")
        if not ok:
            failures.append(f"positive seed{sd} tier={r.tier.value} != definitive")

    # 2) Lookalikes are rejected (detected=False) by a MEANINGFUL trajectory metric
    expected_reason = {
        "saving_propensity": "gini_below_screening_floor",
        "redistribution": "gini_below_screening_floor",
        "bouchaud_mezard": "resource_not_conserved",
        "boltzmann_exchange": "gini_below_screening_floor",
    }
    for reg in fr.CONFIG["regimes"]:
        h = fr.build_substrate(reg)
        r = det.detect(h, model_metadata=None)
        reason = r.primary_metric.get("screening_rejection_reason", "-")
        print(f"LOOKALIKE {reg['label']:18s}: detected={r.detected} reason={reason} "
              f"gini={r.primary_metric.get('gini_final', 0):.3f} drift={r.primary_metric.get('conservation_drift', -1):.4f}")
        if r.detected:
            failures.append(f"lookalike {reg['label']} was DETECTED (false positive)")
        if reason != expected_reason[reg["label"]]:
            failures.append(f"lookalike {reg['label']} reason={reason} != {expected_reason[reg['label']]}")

    # 3) Metadata-flip invariance: flipping the 4 booleans must NOT change the tier
    h, meta = positive_history(1)
    base = det.detect(h, model_metadata=meta).tier.value
    flipped = dict(meta)
    flipped.update(has_conserved_resource=False, has_multiplicative_stake=False,
                   has_saving_propensity=True, has_redistribution=True)
    flip = det.detect(h, model_metadata=flipped).tier.value
    nometa = det.detect(h, model_metadata=None).tier.value
    print(f"METADATA-FLIP: true={base} flipped={flip} nometa={nometa}")
    if not (base == flip == nometa == "definitive"):
        failures.append(f"metadata-flip not invariant: {base}/{flip}/{nometa}")

    # 4) Time-shuffle defeats it: shuffling frame order destroys emergence/monotonicity
    rng = np.random.default_rng(0)
    h, meta = positive_history(1)
    shuffled = [h[i] for i in rng.permutation(len(h))]
    r = det.detect(shuffled, model_metadata=meta)
    print(f"TIME-SHUFFLED positive: detected={r.detected} tier={r.tier.value} "
          f"reason={r.primary_metric.get('screening_rejection_reason', '-')}")
    if r.detected:
        failures.append("time-shuffled positive still detected (not temporally sensitive)")

    print()
    if failures:
        print("FAIL:")
        for f in failures:
            print("  -", f)
        return 1
    print("PASS: tier is trajectory-determined; lookalikes rejected; metadata-flip invariant; time-sensitive.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
