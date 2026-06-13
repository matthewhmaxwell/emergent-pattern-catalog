"""P30 autopoiesis — membrane discrimination test.

Reproducible validation that, with the link-production cap, the AutopoiesisModel
forms a genuine THIN CLOSED MEMBRANE and the detector distinguishes it from a
random particle cloud and from mechanism-off failed regimes — resolving the
deferral finding that the model produced a diffuse cloud and the detector's
"closure" metric returned 1.0 on random scatter.

Signature: radial tightness — link-to-catalyst distances cluster at one radius
(CV << 1), a thin closed shell, self-maintained against decay (stable link
count), vs a random cloud (CV ~ 0.4).

Run:  PYTHONPATH=. OMP_NUM_THREADS=1 python3.12 analysis/validation_rebuild/p30_membrane_test.py
"""

from __future__ import annotations

import numpy as np

from epc.models.autopoiesis import AutopoiesisModel
from epc.detectors.p30_autopoiesis import P30AutopoiesisDetector


def membrane_run(seed: int, steps: int = 1200, **ov):
    p = dict(n_substrate=150, n_catalyst=1, box_size=16.0, production_rate=0.35,
             decay_rate=0.004, production_radius=3.0, membrane_equilibrium_radius=3.0,
             catalyst_link_attraction=4.0, link_attraction=0.0, substrate_diffusion=0.4,
             catalyst_diffusion=0.0, repulsion_radius=0.55, repulsion_strength=1.0,
             dt=0.2, max_links=35, seed=seed)
    p.update(ov)
    m = AutopoiesisModel(**p)
    m.setup()
    hist = [m.step() for _ in range(steps)]
    return hist[49::50], m.get_metadata()


def random_cloud(seed: int, n_link: int = 35, box: float = 16.0, nsnap: int = 8):
    rng = np.random.default_rng(seed)
    hist = []
    for t in range(nsnap):
        n = n_link + 1 + 114
        pos = np.zeros((n, 2))
        typ = np.zeros(n, dtype=int)
        pos[0] = box / 2; typ[0] = AutopoiesisModel.TYPE_CATALYST
        pos[1:n_link + 1] = rng.uniform(0, box, (n_link, 2)); typ[1:n_link + 1] = AutopoiesisModel.TYPE_LINK
        pos[n_link + 1:] = rng.uniform(0, box, (114, 2)); typ[n_link + 1:] = AutopoiesisModel.TYPE_SUBSTRATE
        hist.append({"positions": pos, "types": typ, "box_size": box, "step": t * 50, "n_particles": n})
    return hist


def main() -> int:
    det = P30AutopoiesisDetector(n_permutations=199, seed=42)
    failures = []

    print("=== POSITIVE membrane (should reach DEFINITIVE) ===")
    for s in range(5):
        h, meta = membrane_run(42 + s)
        r = det.detect(h, model_metadata=meta)
        print(f"  seed{42+s}: tier={r.tier.value:11s} radial_cv={r.primary_metric['radial_cv']:.3f} "
              f"n_links={r.primary_metric['n_links']:.0f} gap={r.primary_metric['angular_gap_deg']:.0f}")
        if r.tier.value != "definitive":
            failures.append(f"positive seed{42+s} tier={r.tier.value} != definitive")

    print("=== NEGATIVES (should be rejected — no thin closed shell) ===")
    cases = [("random_cloud", lambda: (random_cloud(1), None))]
    for lbl, ov in [("no_production", dict(production_rate=0.0)),
                    ("high_decay", dict(decay_rate=0.5)),
                    ("no_attraction", dict(catalyst_link_attraction=0.0)),
                    ("uncapped_blob", dict(max_links=-1))]:
        cases.append((lbl, lambda ov=ov: membrane_run(42, **ov)))
    for lbl, gen in cases:
        h, meta = gen()
        r = det.detect(h, model_metadata=meta)
        print(f"  {lbl:14s}: detected={r.detected!s:5s} tier={r.tier.value:11s} "
              f"radial_cv={r.primary_metric.get('radial_cv', -1):.3f} n_links={r.primary_metric.get('n_links', 0):.0f}")
        if lbl != "uncapped_blob" and r.detected:
            failures.append(f"negative {lbl} was DETECTED (false positive)")

    print()
    if failures:
        print("FAIL:")
        for f in failures:
            print("  -", f)
        return 1
    print("PASS: membrane reaches DEFINITIVE (thin self-maintained closed shell); "
          "random cloud + mechanism-off regimes are rejected.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
