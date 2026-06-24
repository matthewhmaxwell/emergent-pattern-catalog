"""Does the structure-factor lens earn its place? Want: patterned/characteristic-scale
systems (stripes, lattice, clusters) HIGH sk_peak; gas / noise nulls LOW; mundane
churn leads LOW."""
import numpy as np
from analysis.blind_spot_probes import PROBES
from epc.metrics.structure_factor import structure_factor_peak
from epc.models.particle_life import particle_life

rows = []
for probe in PROBES:
    try:
        p = probe(0)
        rows.append((probe.__name__, p.get("truth"), structure_factor_peak(p["history"])))
    except Exception:
        rows.append((probe.__name__, "err", None))
for i in [180, 270, 288, 298, 205]:
    rng = np.random.default_rng(i); K = int(rng.integers(3, 7)); F = rng.uniform(-1.0, 1.0, (K, K))
    h, _ = particle_life(i, K=K, F=F, N=200, steps=250)
    rows.append((f"pl_{i}", ("lead" if i != 205 else "pl-class"), structure_factor_peak(h)))

print(f"{'system':<24}{'truth':<12}{'sk_peak':>9}")
for nm, truth, r in sorted(rows, key=lambda x: -((x[2] or {}).get('sk_peak', -1))):
    print(f"{nm:<24}{str(truth):<12}{r['sk_peak']:>9.2f}" if r else f"{nm:<24}{str(truth):<12}   (n/a)")

em = [r['sk_peak'] for nm, t, r in rows if t == 'emergent' and r]
nl = [r['sk_peak'] for nm, t, r in rows if t == 'null' and r]
ld = [r['sk_peak'] for nm, t, r in rows if t == 'lead' and r]
print(f"\nemergent sk_peak: mean {np.mean(em):.1f} range {min(em):.1f}-{max(em):.1f} (n={len(em)})")
print(f"null     sk_peak: {[round(x,1) for x in nl]}")
print(f"leads    sk_peak: {[round(x,1) for x in ld]}")
