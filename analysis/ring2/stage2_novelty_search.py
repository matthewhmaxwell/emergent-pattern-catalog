"""Stage-2 novelty search: run the FULL ring2_descriptor over a generative-family sweep,
cluster observations by lens fingerprint, and surface the tripwire's complex+unclassified
leads for vetting.

Upgrade over Stage-1 (which used the tripwire alone): the broadened battery means a 'lead'
is structure that NO admitted lens can name (checked against 5 lenses + the named emergence
channels), so lead PRECISION is higher. Honest prior (Stage-1): genuine novelty is rare and
leads are usually mundane churn; the deliverable is an AUDITED FRONTIER, not a manufactured
claim. Classification here uses the fast named-lens baseline; survivors get the full
catalogued-detector check in vetting.

Usage: stage2_novelty_search.py <N> <family>   (family: particle_life | lenia)
"""
import json
import sys
import signal
from collections import Counter

import numpy as np

from epc.phase2a.ring2_descriptor import ring2_descriptor

N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
FAMILY = sys.argv[2] if len(sys.argv) > 2 else "particle_life"
PER_TIMEOUT = 90


class _TO(Exception):
    pass


signal.signal(signal.SIGALRM, lambda s, f: (_ for _ in ()).throw(_TO()))


def gen(i):
    rng = np.random.default_rng(1000 + i)
    if FAMILY == "lenia":
        from epc.models.lenia import lenia
        nc = int(rng.integers(1, 6))
        mu = float(rng.uniform(0.10, 0.20)); sig = float(rng.uniform(0.012, 0.020))
        h = lenia(i, N=64, mu=mu, sigma=sig, steps=400, record=60, n_creatures=nc)
        return h, {"family": "lenia", "mu": mu, "sigma": sig, "nc": nc}
    from epc.models.particle_life import particle_life
    K = int(rng.integers(3, 7)); F = rng.uniform(-1, 1, (K, K))
    h, meta = particle_life(i, K=K, F=F, N=200, steps=250)
    return h, (meta or {})


rows = []
for i in range(N):
    try:
        signal.alarm(PER_TIMEOUT)
        hist, meta = gen(i)
        d = ring2_descriptor(hist, meta)
        signal.alarm(0)
        rows.append({"i": i, "fired": d["lenses_fired"], "em_kind": d["em_kind"],
                     "em": d["em_score"], "tripped": d["tripped"], "cx": d["mf_complex"],
                     "C": d.get("mf_C"), "struct": d.get("mf_struct"), "psi": d.get("mf_psi"),
                     "sk": d.get("sk_peak"), "h1": d.get("h1_max"),
                     "field_loops": d.get("field_loops"), "lacun": d.get("lacunarity"),
                     "dte_mte": d.get("dte_mean_te"), "dte_dir": d.get("dte_directionality")})
    except _TO:
        signal.alarm(0); rows.append({"i": i, "error": "timeout"})
    except Exception as e:
        signal.alarm(0); rows.append({"i": i, "error": repr(e)[:90]})
    if (i + 1) % 20 == 0:
        print(f"...{i + 1}/{N}", flush=True)

ok = [r for r in rows if "error" not in r]
errs = [r for r in rows if "error" in r]
clusters = Counter((tuple(sorted(r["fired"])), r["em_kind"]) for r in ok)
leads = [r for r in ok if r["tripped"]]

print(f"\n=== Stage-2 {FAMILY} sweep: {len(ok)} ok / {len(errs)} errors (N={N}) ===")
print("fingerprint clusters  (lenses_fired | em_kind) : count")
for (fired, kind), c in clusters.most_common():
    print(f"  {('+'.join(fired) or '(none)'):<30} {str(kind):<16} {c}")
print(f"\nTRIPWIRE LEADS (complex + unclassified): {len(leads)} / {len(ok)}")
for r in leads:
    print(f"  #{r['i']} em={r['em']} C={r['C']} struct={r['struct']} psi={r['psi']} "
          f"fired={'+'.join(r['fired']) or '(none)'}")

out = f"analysis/ring2/stage2_results_{FAMILY}.json"
with open(out, "w") as fh:
    json.dump(rows, fh, indent=1, default=str)
print(f"\nsaved {out}")
