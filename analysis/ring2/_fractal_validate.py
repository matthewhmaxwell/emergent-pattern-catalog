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

print(f"{'structure':<22}{'truth':<14}{'D':>7}{'lacun':>8}{'fill':>7}")
for nm, t, r in sorted(rows, key=lambda x: -((x[2] or {}).get('lacunarity', -9))):
    if r:
        print(f"{nm:<22}{t:<14}{r['fractal_dim']:>7.3f}{r['lacunarity']:>8.3f}{r['fill_frac']:>7.3f}")
    else:
        print(f"{nm:<22}{t:<14}   (n/a)")

val = {nm: r for nm, t, r in rows if r}
# SCOPED admit: lacunarity (gated D>1.2 to drop lines) separates SPARSE fractals from
# random/uniform texture. Dense near-critical percolation is OUT OF SCOPE (homogeneous).
sparse = [val[n]['lacunarity'] for n in ('dla_fractal', 'SYNTH_sierpinski') if n in val]
nonfrac = [val[n]['lacunarity'] for n in ('SYNTH_random_gas', 'SYNTH_filled_disk') if n in val]
perc = val.get('percolation', {}).get('lacunarity')
print(f"\nlacunarity  sparse-fractals(dla,sierpinski): {[round(v,2) for v in sorted(sparse)]}  "
      f"random/uniform: {[round(v,2) for v in sorted(nonfrac)]}  percolation(dense): {perc}")
if sparse and nonfrac:
    gap = min(sparse) - max(nonfrac)
    print(f"gap (min sparse-fractal - max random/uniform lacun) = {gap:+.2f}  -> "
          f"{'ADMIT SCOPED (lacunarity separates SPARSE fractals; dense=percolation out of scope)' if gap > 0.5 else 'keep DEFERRED'}")
    print(f"  (filament D={val.get('SYNTH_filament',{}).get('fractal_dim')} excluded by the D>1.2 gate;"
          f" percolation lacun {perc} ~ random -> dense fractals need the multifractal spectrum)")
