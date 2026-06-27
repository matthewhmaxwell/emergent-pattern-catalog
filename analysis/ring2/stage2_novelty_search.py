"""Stage-2 novelty search: run the FULL ring2_descriptor over a generative-family sweep,
cluster observations by lens fingerprint, and surface BOTH (a) the tripwire's complex+unclassified
leads and (b) the FULL per-config fingerprint (all 21 admitted lenses) for downstream lens-space
outlier analysis -- the discovery mechanism the comprehensive battery enables.

A 'tripwire lead' is structure the em/MPR-C tripwire flags as complex+unclassified. The richer
battery does not change WHICH configs trip (the tripwire never consults the Ring-2 lenses) -- it
enriches each fingerprint so leads can be vetted across 21 axes AND so lens-space outliers (configs
structurally unusual across the battery, which em/C alone may miss) can be found separately.
Honest prior: genuine novelty is rare; the deliverable is an AUDITED FRONTIER.

Usage: stage2_novelty_search.py <N> <family> [seed_base]
  family: particle_life | lenia | reaction_diffusion | kuramoto
"""
import json
import sys
import signal
from collections import Counter

import numpy as np

from epc.phase2a.ring2_descriptor import ring2_descriptor

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
FAMILY = sys.argv[2] if len(sys.argv) > 2 else "particle_life"
SEED_BASE = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
PER_TIMEOUT = 180


class _TO(Exception):
    pass


signal.signal(signal.SIGALRM, lambda s, f: (_ for _ in ()).throw(_TO()))


def gen(i):
    rng = np.random.default_rng(SEED_BASE + i)
    if FAMILY == "lenia":
        from epc.models.lenia import lenia
        soup = bool(rng.random() < 0.2)
        nc = int(rng.integers(1, 7))
        mu = float(rng.uniform(0.10, 0.22)); sig = float(rng.uniform(0.010, 0.022))
        h = lenia(i, N=64, mu=mu, sigma=sig, steps=400, record=60, n_creatures=nc, soup=soup)
        return h, {"family": "lenia", "mu": mu, "sigma": sig, "nc": nc, "soup": soup}
    if FAMILY == "reaction_diffusion":
        from epc.models.reaction_diffusion import reaction_diffusion
        F = float(rng.uniform(0.010, 0.070)); k = float(rng.uniform(0.045, 0.070))
        h = reaction_diffusion(i, F=F, k=k, N=96, steps=8000, record=24)
        return h, {"family": "reaction_diffusion", "F": F, "k": k}
    if FAMILY == "kuramoto":
        from epc.models.kuramoto_lattice import kuramoto_lattice
        K = float(rng.uniform(0.1, 3.0)); osp = float(rng.uniform(0.0, 2.0))
        init = ["random", "spiral", "plane"][int(rng.integers(0, 3))]
        h = kuramoto_lattice(i, N=40, K=K, omega_spread=osp, init=init, steps=3000, record=120)
        return h, {"family": "kuramoto", "K": K, "omega_spread": osp, "init": init}
    from epc.models.particle_life import particle_life
    K = int(rng.integers(3, 9)); F = rng.uniform(-1.5, 1.5, (K, K))
    h, meta = particle_life(i, K=K, F=F, N=200, steps=250)
    return h, (meta or {})


rows = []
for i in range(N):
    try:
        signal.alarm(PER_TIMEOUT)
        hist, meta = gen(i)
        d = ring2_descriptor(hist, meta)
        signal.alarm(0)
        row = {"i": i, "fired": d["lenses_fired"]}
        row.update({k: v for k, v in d.items() if k != "lenses_fired"})  # FULL 21-lens fingerprint
        rows.append(row)
    except _TO:
        signal.alarm(0); rows.append({"i": i, "error": "timeout"})
    except Exception as e:
        signal.alarm(0); rows.append({"i": i, "error": repr(e)[:90]})
    if (i + 1) % 20 == 0:
        print(f"...{i + 1}/{N}", flush=True)

ok = [r for r in rows if "error" not in r]
errs = [r for r in rows if "error" in r]
clusters = Counter((tuple(sorted(r["fired"])), r["em_kind"]) for r in ok)
leads = [r for r in ok if r.get("tripped")]

print(f"\n=== Stage-2 {FAMILY} sweep: {len(ok)} ok / {len(errs)} errors (N={N}, seed_base={SEED_BASE}) ===")
print("fingerprint clusters  (#lenses_fired | em_kind) : count")
for (fired, kind), c in clusters.most_common():
    print(f"  {len(fired):>2} lenses  {str(kind):<18} {c}")
if errs:
    print(f"errors: {Counter(r['error'] for r in errs)}")
print(f"\nTRIPWIRE LEADS (complex + unclassified): {len(leads)} / {len(ok)}")
for r in leads:
    print(f"  #{r['i']} em={r.get('em_score')} C={r.get('mf_C')} struct={r.get('mf_struct')} "
          f"psi={r.get('mf_psi')} fired={len(r['fired'])} lenses")

out = f"analysis/ring2/stage2_results_{FAMILY}.json"
with open(out, "w") as fh:
    json.dump(rows, fh, indent=1, default=str)
print(f"\nsaved {out}  (full 21-lens fingerprints for lens-space outlier analysis)")
