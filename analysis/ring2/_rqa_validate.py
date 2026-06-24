"""Does the recurrence-determinism lens earn its place? Run it over the probe corpus
(known-emergent + nulls), the 4 mundane leads, and a classified particle-life. A lens
earns admission if it SEPARATES emergent (high DET) from null (low DET) and adds a
discriminating axis (e.g., ranks the mundane churn leads below genuine structure)."""
import numpy as np
from analysis.blind_spot_probes import PROBES
from epc.metrics.recurrence import recurrence_determinism
from epc.models.particle_life import particle_life

rows = []
for probe in PROBES:
    try:
        p = probe(0)
        rows.append((probe.__name__, p.get("truth"), recurrence_determinism(p["history"])))
    except Exception as e:
        rows.append((probe.__name__, "err", None))
for i in [180, 270, 288, 298, 205]:
    rng = np.random.default_rng(i); K = int(rng.integers(3, 7)); F = rng.uniform(-1.0, 1.0, (K, K))
    h, _ = particle_life(i, K=K, F=F, N=200, steps=250)
    rows.append((f"pl_{i}", ("lead" if i != 205 else "pl-class"), recurrence_determinism(h)))

print(f"{'system':<24}{'truth':<12}{'DET':>7}{'LAM':>7}{'RR':>7}")
for nm, truth, r in sorted(rows, key=lambda x: -((x[2] or {}).get('DET', -1))):
    if r:
        print(f"{nm:<24}{str(truth):<12}{r['DET']:>7.3f}{r['LAM']:>7.3f}{r['RR']:>7.3f}")
    else:
        print(f"{nm:<24}{str(truth):<12}   (no trajectory)")

em = [r['DET'] for nm, t, r in rows if t == 'emergent' and r]
nl = [r['DET'] for nm, t, r in rows if t == 'null' and r]
ld = [r['DET'] for nm, t, r in rows if t == 'lead' and r]
print(f"\nemergent DET: mean {np.mean(em):.3f} range {min(em):.3f}-{max(em):.3f}  (n={len(em)})")
print(f"null     DET: mean {np.mean(nl):.3f} range {min(nl):.3f}-{max(nl):.3f}  (n={len(nl)})")
print(f"leads    DET: {[round(x,3) for x in ld]}")
