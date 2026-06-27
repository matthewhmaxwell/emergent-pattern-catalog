"""Does chaos_dimension earn its place where recurrence failed? Canonical scalar series:
logistic chaos (r=4), Lorenz-x, sine (periodic), white noise, AR(1) colored noise. Want:
zero_one_K ~1 for chaos/noise & ~0 for sine; determinism HIGH for chaos+sine (deterministic),
LOW for white/AR1 noise (D2 keeps growing). The determinism axis is what separates chaos from
noise -- the thing a broadband spectrum cannot."""
import numpy as np
from epc.metrics.chaos_dimension import chaos_dimension

T = 2500


def S(x):
    return [{"scalar": float(v)} for v in x]


def logistic(seed=0):
    x = 0.4; out = []
    for _ in range(T + 100):
        x = 4.0 * x * (1 - x); out.append(x)
    return np.array(out[100:])


def lorenz(seed=0, sub=10):
    dt = 0.01; x, y, z = 1.0, 1.0, 1.0; out = []
    for _ in range(T * sub + 3000):
        dx = 10 * (y - x); dy = x * (28 - z) - y; dz = x * y - 8 / 3 * z
        x += dt * dx; y += dt * dy; z += dt * dz; out.append(x)
    return np.array(out[3000::sub])[:T]                      # subsample to decorrelate (0-1 test)


def sine(seed=0):
    return np.sin(np.linspace(0, 120 * np.pi, T))


def white(seed=0):
    return np.random.default_rng(seed).standard_normal(T)


def ar1(seed=0):
    r = np.random.default_rng(seed); x = 0.0; out = []
    for _ in range(T):
        x = 0.85 * x + r.standard_normal(); out.append(x)
    return np.array(out)


print(f"{'system':<12}{'zero_one_K':>12}{'determinism':>13}")
rows = {}
for nm, fn in [("logistic", logistic), ("lorenz", lorenz), ("sine", sine),
               ("white_noise", white), ("ar1_noise", ar1)]:
    r = chaos_dimension(S(fn(1)))
    rows[nm] = r
    print(f"{nm:<12}{r['zero_one_K']:>12.3f}{r['determinism']:>13.3f}")

chaos_det = min(rows["logistic"]["determinism"], rows["lorenz"]["determinism"])
noise_det = max(rows["white_noise"]["determinism"], rows["ar1_noise"]["determinism"])
K_sine = rows["sine"]["zero_one_K"]; K_chaos = min(rows["logistic"]["zero_one_K"], rows["lorenz"]["zero_one_K"])
print(f"\ndeterminism: chaos(min) {chaos_det:.2f} vs noise(max) {noise_det:.2f} | K: chaos {K_chaos:.2f} vs sine {K_sine:.2f}")
print(f"VERDICT: {'ADMIT' if chaos_det > noise_det + 0.2 and K_chaos > 0.6 and K_sine < 0.5 else 'review'}")
