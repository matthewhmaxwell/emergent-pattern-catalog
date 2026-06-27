"""Does coarsening earn its place? Allen-Cahn quench (non-conserved coarsening, domains grow,
n~1/2 = positive) vs a STATIC stripe pattern (fixed wavelength, n~0) vs white noise (no domains).
Want: Allen-Cahn growth_exponent clearly >0 with good fit; static ~0; noise ~0/poor fit."""
import numpy as np
from epc.metrics.coarsening import coarsening

G = 96


def _step(phi, dt):
    lap = (np.roll(phi, 1, 0) + np.roll(phi, -1, 0) +
           np.roll(phi, 1, 1) + np.roll(phi, -1, 1) - 4 * phi)
    return phi + dt * (lap + phi - phi ** 3)


def allen_cahn(seed=0, warmup=20, steps_per=8, n=18, dt=0.15):
    r = np.random.default_rng(seed)
    phi = 0.1 * r.standard_normal((G, G))
    for _ in range(warmup):                                   # skip nucleation transient
        phi = _step(phi, dt)
    out = []
    for _ in range(n):
        for _ in range(steps_per):
            phi = _step(phi, dt)
        out.append({"field": phi.copy()})
    return out


def static_stripes(seed=0, n=14):
    y, x = np.indices((G, G)) * (2 * np.pi / G)
    base = np.cos(8 * x)
    r = np.random.default_rng(seed)
    return [{"field": base + 0.02 * r.standard_normal((G, G))} for _ in range(n)]


def noise(seed=0, n=14):
    r = np.random.default_rng(seed)
    return [{"field": r.standard_normal((G, G))} for _ in range(n)]


print(f"{'system':<16}{'growth_exponent':>16}{'fit_quality':>13}")
rows = {}
for nm, fn in [("allen_cahn", allen_cahn), ("static_stripes", static_stripes), ("noise", noise)]:
    r = coarsening(fn(1))
    rows[nm] = r
    print(f"{nm:<16}{r['growth_exponent']:>16.3f}{r['fit_quality']:>13.3f}")

ac = rows["allen_cahn"]; st = rows["static_stripes"]; ns = rows["noise"]
# Honest criterion: coarsening = clean monotonic scale-growth POWER LAW (positive n, high R2);
# static = no growth (n~0); noise = no clean power law (low R2, gated out).
coarsen_ok = ac["growth_exponent"] > 0.05 and ac["fit_quality"] > 0.8
static_ok = abs(st["growth_exponent"]) < 0.05
noise_gated = ns["fit_quality"] < 0.3
print(f"\nAllen-Cahn n={ac['growth_exponent']:.3f} R2 {ac['fit_quality']:.2f} (clean growth) | static n={st['growth_exponent']:.3f} | noise R2 {ns['fit_quality']:.2f} (gated)")
print(f"VERDICT: {'ADMIT (clean-growth vs static vs noise; abs exponent needs full scaling window)' if coarsen_ok and static_ok and noise_gated else 'review'}")
