"""Does cross_substrate earn its place? Field pairs (A,B):
- coupled_white: B = shifted A -> coherent at all k, both autos white.
- buried_coupled: A=C+3n, B=C+3n with C a hidden pattern -> EACH auto noise-dominated (flat),
  but jointly coherent at C's wavelength (the reciprocal-cross-diffusion 'absent from marginals' case).
- independent_white: A,B independent -> coherence ~0.
- independent_patterned: A,B each individually patterned but UNCORRELATED -> coherence ~0 despite
  strong marginal structure (the lens must not be fooled by marginal peaks).
Want: coupled/buried coherence HIGH; both independent cases LOW."""
import numpy as np
from epc.metrics.cross_substrate import cross_substrate

G = 64
yy, xx = np.indices((G, G)) * (2 * np.pi / G)


def H(fn, n=8):
    return [{"field_a": fn(t)[0], "field_b": fn(t)[1]} for t in range(n)]


def coupled_white(t):
    r = np.random.default_rng(t)
    A = r.standard_normal((G, G))
    B = np.roll(A, 5, axis=0) + 0.1 * r.standard_normal((G, G))
    return A, B


def buried_coupled(t):
    r = np.random.default_rng(t)
    C = np.cos(6 * xx + 6 * yy)
    A = C + 3.0 * r.standard_normal((G, G))
    B = C + 3.0 * r.standard_normal((G, G))
    return A, B


def independent_white(t):
    r = np.random.default_rng(t); r2 = np.random.default_rng(1000 + t)
    return r.standard_normal((G, G)), r2.standard_normal((G, G))


def independent_patterned(t):
    r = np.random.default_rng(t); r2 = np.random.default_rng(1000 + t)
    A = np.cos(6 * xx) + 0.3 * r.standard_normal((G, G))
    B = np.cos(6 * yy) + 0.3 * r2.standard_normal((G, G))     # orthogonal orientation, independent
    return A, B


print(f"{'system':<22}{'cross_coherence':>16}")
rows = {}
for nm, fn in [("coupled_white", coupled_white), ("buried_coupled", buried_coupled),
               ("independent_white", independent_white), ("independent_patterned", independent_patterned)]:
    r = cross_substrate(H(fn))
    rows[nm] = r
    print(f"{nm:<22}{r['cross_coherence']:>16.4f}")

coupled = min(rows["coupled_white"]["cross_coherence"], rows["buried_coupled"]["cross_coherence"])
indep = max(rows["independent_white"]["cross_coherence"], rows["independent_patterned"]["cross_coherence"])
print(f"\ncoupled(min) {coupled:.3f} vs independent(max) {indep:.3f}")
print(f"VERDICT: {'ADMIT' if coupled > 2 * indep + 0.1 else 'review'}")
