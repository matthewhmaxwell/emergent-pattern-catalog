"""Round-2 #64: heavy-tail / power-law (self-organized criticality) channel.

Step 1 (diagnose): does the CURRENT instrument flag a power-law/SOC system, and is
there a clean power-law-vs-exponential signal to build a channel on? Probe = BTW
sandpile (avalanche sizes are power-law); null = random activity (exponential sizes).
"""
from __future__ import annotations

import numpy as np

import analysis.run_phase2a_panel as R
from epc.phase2a.emergence import generic_emergence


def plaw_llr(sizes, xmin=1.0):
    """Log-likelihood ratio power-law vs exponential (continuous approx, Clauset-style).
    Returns (alpha, LLR, decades). LLR>0 ⇒ power-law preferred over exponential."""
    x = np.asarray([s for s in sizes if s >= xmin], dtype=float)
    if x.size < 20:
        return None
    n = x.size
    alpha = 1.0 + n / np.sum(np.log(x / (xmin * 0.999)))
    lam = 1.0 / (x.mean() - xmin + 1e-9)
    ll_pl = np.sum(np.log((alpha - 1.0) / xmin) - alpha * np.log(x / xmin))
    ll_exp = np.sum(np.log(lam) - lam * (x - xmin))
    decades = np.log10(x.max() / max(x.min(), 1e-9))
    return float(alpha), float(ll_pl - ll_exp), float(decades)


def report(name, sizes, history=None):
    em = generic_emergence(history)["score"] if history is not None else float("nan")
    r = plaw_llr(sizes)
    if r is None:
        print(f"{name:<24} generic_em={em:.3f}  (too few events)")
        return
    a, llr, dec = r
    flag = "POWER-LAW" if (llr > 0 and dec >= 1.3) else "not-pl"
    print(f"{name:<24} generic_em={em:.3f}  alpha={a:.2f}  LLR(pl-exp)={llr:+.1f}  decades={dec:.2f}  -> {flag}")


# SOC probe: BTW sandpile (catalog P14 positive). Its history carries avalanche_sizes.
runs, meta = R.build_p14_positives(n_seeds=1)
hist = runs[0]
print("sandpile frame keys:", sorted(hist[-1].keys())[:10])
# pull avalanche sizes from the history
sizes = []
for f in hist:
    for k in ("avalanche_sizes", "avalanche_size", "size"):
        if k in f:
            v = f[k]
            sizes.extend(v if hasattr(v, "__len__") else [v])
            break
sizes = [s for s in sizes if s and s > 0]
report("sandpile (SOC)", sizes, hist)

# nulls: exponential and uniform size distributions (NOT power-law)
rng = np.random.default_rng(0)
report("exp-sizes (null)", rng.exponential(5.0, 4000) + 1)
report("uniform-sizes (null)", rng.uniform(1, 50, 4000))
print("\nWANT: sandpile -> POWER-LAW (LLR>0, ≥1.3 decades); nulls -> not-pl. "
      "And note whether generic_em already flags the sandpile.")
