"""Does multifractal earn its place (dense-fractal depth item)? 2D binomial multiplicative cascade
(strongly multifractal, wide Delta-alpha, spatial structure) vs uniform field (Delta-alpha~0) vs
iid lognormal (heavy-tailed values, NO spatial structure -> surrogate gate must zero it). Want:
cascade mf_width clearly >0; uniform ~0; iid-heavytail ~0 (shuffle leaves its Delta-alpha unchanged)."""
import numpy as np
from epc.metrics.multifractal import multifractal

G = 96


def H(field, n=4):
    return [{"field": field} for _ in range(n)]


def cascade(seed=0, levels=6, p=0.3):
    r = np.random.default_rng(seed)
    m = np.ones((1, 1))
    for _ in range(levels):
        w = r.choice([p, 1 - p], size=(m.shape[0] * 2, m.shape[1] * 2))
        m = np.kron(m, np.ones((2, 2))) * w
    # resize to GxG by cropping/padding
    s = m.shape[0]
    mm = m[:G, :G] if s >= G else np.pad(m, ((0, G - s), (0, G - s)), mode="wrap")
    return mm


def uniform(seed=0):
    r = np.random.default_rng(seed)
    return np.ones((G, G)) + 0.01 * r.standard_normal((G, G))


def iid_heavytail(seed=0):
    r = np.random.default_rng(seed)
    return r.lognormal(0, 2.0, (G, G))      # heavy-tailed VALUES, spatially iid (no structure)


print(f"{'system':<16}{'mf_width':>10}{'alpha_span':>12}")
rows = {}
for nm, fn in [("cascade", cascade), ("uniform", uniform), ("iid_heavytail", iid_heavytail)]:
    r = multifractal(H(fn(1)))
    rows[nm] = r
    print(f"{nm:<16}{r['mf_width']:>10.4f}{r['alpha_span']:>12.4f}")

casc = rows["cascade"]["mf_width"]
nulls = max(rows["uniform"]["mf_width"], rows["iid_heavytail"]["mf_width"])
print(f"\ncascade mf_width {casc:.3f} vs max null(uniform/iid-heavytail) {nulls:.3f}")
print(f"VERDICT: {'ADMIT (surrogate gate zeros iid heavy-tail)' if casc > 2 * nulls + 0.1 and casc > 0.2 else 'review'}")
