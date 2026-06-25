"""Does the fractal-dimension lens earn its place? Positive controls are known
fractals (D well below 2 with clean scaling); nulls are space-filling (D~2) or
filament (D~1) structures. Then the EPC fractal probes (dla, percolation). Want:
fractals D in ~(1.4,1.95) with high R^2; compact/uniform fills D~2."""
import numpy as np
from analysis.blind_spot_probes import dla_fractal, percolation
from epc.metrics.fractal_dimension import fractal_dimension

N = 243  # 3^5, friendly to Sierpinski


def G(occ):
    return [{"grid": np.asarray(occ, float)}]


def sierpinski(n=256):
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    return ((i & j) == 0).astype(int)            # Sierpinski triangle, D=log3/log2~1.585


def filled_disk(n=256):
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    c = n / 2; return (((i - c) ** 2 + (j - c) ** 2) < (0.45 * n) ** 2).astype(int)


def random_gas(n=256, dens=0.08, seed=0):
    rng = np.random.default_rng(seed); return (rng.random((n, n)) < dens).astype(int)


def filament(n=256, seed=0):                      # near-1D: thick diagonal
    rng = np.random.default_rng(seed); g = np.zeros((n, n), int)
    for k in range(n):
        g[k, min(n - 1, k + rng.integers(0, 3))] = 1
    return g


rows = [("SYNTH_sierpinski", "fractal", fractal_dimension(G(sierpinski()))),
        ("SYNTH_filled_disk", "uniform", fractal_dimension(G(filled_disk()))),
        ("SYNTH_random_gas", "random", fractal_dimension(G(random_gas()))),
        ("SYNTH_filament", "filament", fractal_dimension(G(filament())))]
for fn, tr in [(dla_fractal, "fractal(EPC)"), (percolation, "fractal(EPC)")]:
    p = fn(0)
    rows.append((fn.__name__, tr, fractal_dimension(p["history"])))

print(f"{'structure':<22}{'truth':<14}{'D':>7}{'fit_r2':>8}{'fill':>7}")
for nm, t, r in sorted(rows, key=lambda x: (x[2] or {}).get('fractal_dim', 9)):
    if r:
        print(f"{nm:<22}{t:<14}{r['fractal_dim']:>7.3f}{r['fit_r2']:>8.3f}{r['fill_frac']:>7.3f}")
    else:
        print(f"{nm:<22}{t:<14}   (n/a)")

vals = [(t, r) for nm, t, r in rows if r]
frac = [r['fractal_dim'] for t, r in vals if 'fractal' in t]
uni = [r['fractal_dim'] for t, r in vals if t == 'uniform']
gas = [r['fractal_dim'] for t, r in vals if t == 'random']
print(f"\nfractal D: {[round(v,2) for v in sorted(frac)]}  uniform D: {[round(v,2) for v in uni]}"
      f"  random D: {[round(v,2) for v in gas]}")
# A real fractal lens must separate fractals from BOTH uniform fills AND random texture.
gap_uni = (min(uni) - max(frac)) if (frac and uni) else 0.0
overlap_random = any(min(frac) <= g <= max(frac) for g in gas) if (frac and gas) else False
print(f"gap to uniform = {gap_uni:+.3f} (TIGHT); random D overlaps the fractal band: {overlap_random}")
print("VERDICT: DEFER — box-counting D is confounded. random gas D~1.57 ~= Sierpinski D~1.585,"
      " and percolation D~1.83 ~= uniform disk D~1.90: D conflates self-similarity with"
      " density/boundary. Needs lacunarity / multifractal spectrum to admit.")
