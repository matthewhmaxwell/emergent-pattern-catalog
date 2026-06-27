"""Does hyperuniformity earn its place? Point patterns with KNOWN long-wavelength scaling:
poisson (alpha~2), jittered lattice (hyperuniform, alpha~1), clustered/Matern clumps (alpha>2 =
giant fluctuations). Want number_variance_exponent to order them hyperuniform < poisson < clustered."""
import numpy as np
from epc.metrics.hyperuniformity import hyperuniformity

L = 100.0


def H(P):
    return [{"positions": P}]


def poisson(seed=0, n=3000):
    return np.random.default_rng(seed).uniform(0, L, (n, 2))


def lattice(seed=0, m=55):
    g = np.linspace(2, L - 2, m)
    X, Y = np.meshgrid(g, g)
    P = np.column_stack([X.ravel(), Y.ravel()])
    r = np.random.default_rng(seed)
    return P + 0.25 * (g[1] - g[0]) * r.standard_normal(P.shape)   # small jitter (stays hyperuniform)


def clustered(seed=0, n_clumps=60, per=50):
    r = np.random.default_rng(seed)
    centers = r.uniform(0, L, (n_clumps, 2))
    pts = np.vstack([c + 2.0 * r.standard_normal((per, 2)) for c in centers])
    return np.clip(pts, 0, L)


print(f"{'system':<14}{'number_variance_exponent':>26}")
rows = {}
for nm, fn in [("poisson", poisson), ("lattice_hyperu", lattice), ("clustered", clustered)]:
    r = hyperuniformity(H(fn(1)))
    rows[nm] = r
    print(f"{nm:<14}{r['number_variance_exponent']:>26.3f}")

hu = rows["lattice_hyperu"]["number_variance_exponent"]
po = rows["poisson"]["number_variance_exponent"]
cl = rows["clustered"]["number_variance_exponent"]
print(f"\nhyperuniform {hu:.2f} < poisson {po:.2f} < clustered {cl:.2f}")
print(f"VERDICT: {'ADMIT' if hu < po - 0.3 and cl > po + 0.3 else 'review'}")
