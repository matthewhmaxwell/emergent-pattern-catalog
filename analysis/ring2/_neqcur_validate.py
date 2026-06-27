"""Does nonequilibrium_current earn its place? Two-channel series: directed_cycle (A=cos,B=sin =
circulating nonequilibrium flux), standing_oscillation (A,B IN PHASE = oscillates but reversible,
no net rotation -- the case spectral cannot tell from a directed cycle), equilibrium_noise
(independent OU). Want: directed circulation HIGH; standing oscillation + noise ~0."""
import numpy as np
from epc.metrics.nonequilibrium_current import nonequilibrium_current

T = 400


def C(A, B):
    return [{"chan_a": float(a), "chan_b": float(b)} for a, b in zip(A, B)]


def directed_cycle(seed=0):
    r = np.random.default_rng(seed); t = np.linspace(0, 40 * np.pi, T)
    return np.cos(t) + 0.1 * r.standard_normal(T), np.sin(t) + 0.1 * r.standard_normal(T)


def standing_oscillation(seed=0):
    r = np.random.default_rng(seed); t = np.linspace(0, 40 * np.pi, T)
    return np.cos(t) + 0.1 * r.standard_normal(T), np.cos(t) + 0.1 * r.standard_normal(T)   # in phase


def equilibrium_noise(seed=0):
    r1 = np.random.default_rng(seed); r2 = np.random.default_rng(1000 + seed)
    a = b = 0.0; A, B = [], []
    for _ in range(T):
        a = 0.9 * a + r1.standard_normal(); b = 0.9 * b + r2.standard_normal()
        A.append(a); B.append(b)
    return np.array(A), np.array(B)


print(f"{'system':<22}{'circulation':>12}")
rows = {}
for nm, fn in [("directed_cycle", directed_cycle), ("standing_oscillation", standing_oscillation),
               ("equilibrium_noise", equilibrium_noise)]:
    r = nonequilibrium_current(C(*fn(1)))
    rows[nm] = r
    print(f"{nm:<22}{r['circulation']:>12.4f}")

dc = rows["directed_cycle"]["circulation"]
nulls = max(rows["standing_oscillation"]["circulation"], rows["equilibrium_noise"]["circulation"])
print(f"\ndirected_cycle {dc:.3f} vs max null (standing-osc/noise) {nulls:.3f}")
print(f"VERDICT: {'ADMIT (directed flux vs reversible oscillation+noise)' if dc > 5 * nulls + 0.1 else 'review'}")
