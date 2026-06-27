"""Does sync_texture earn its place? A real Abrams-Strogatz nonlocal-Kuramoto chimera (coherent
+ incoherent coexist) vs uniform sync vs full incoherence. Want: chimera local_R_spread HIGH;
uniform sync + incoherent LOW -- even though sync has global R~1 and incoherent global R~0 (so
the GLOBAL channel cannot tell chimera from one of them, but the spread can)."""
import numpy as np
from epc.metrics.sync_texture import sync_texture

G = 64


def H(theta_fn, n=6):
    return [{"phases": theta_fn(t)} for t in range(n)]


def uniform_sync(t):
    r = np.random.default_rng(t)
    return 0.4 + 0.05 * r.standard_normal((G, G))            # R~1 everywhere


def incoherent(t):
    return np.random.default_rng(t).uniform(0, 2 * np.pi, (G, G))   # R~0 everywhere


def chimera(t):
    # left half: smooth coherent (phase ramp); right half: incoherent random
    r = np.random.default_rng(100 + t)
    th = np.zeros((G, G))
    yy, xx = np.indices((G, G))
    th[:, :G // 2] = 0.3 * yy[:, :G // 2]                    # coherent: smooth gradient
    th[:, G // 2:] = r.uniform(0, 2 * np.pi, (G, G // 2))    # incoherent
    return th


# a genuine dynamical chimera: Kuramoto-Battogtokh ring with exponential nonlocal kernel.
# dtheta/dt = omega - Im(e^{i(theta+alpha)} conj(m)),  m_i = sum_j G_ij e^{i theta_j}
def as_chimera(t, N=256, alpha=1.457, kappa=4.0, steps=2500, dt=0.05, seed=7):
    r = np.random.default_rng(seed)
    x = np.linspace(0, 2 * np.pi, N, endpoint=False)
    d = np.abs(x[:, None] - x[None, :]); d = np.minimum(d, 2 * np.pi - d)
    Gk = np.exp(-kappa * d); Gk /= Gk.sum(1, keepdims=True)
    # seed: coherent bump + random elsewhere (the standard chimera initial condition)
    theta = 6.0 * np.exp(-0.76 * ((x - np.pi) ** 2)) * r.uniform(-1, 1, N)
    theta = theta % (2 * np.pi)
    for _ in range(steps):
        m = Gk @ np.exp(1j * theta)
        drive = -np.imag(np.exp(1j * (theta + alpha)) * np.conj(m))
        theta = (theta + dt * drive) % (2 * np.pi)
    return theta                                              # 1D ring -> lens uses 1D local R


print(f"{'system':<16}{'local_R_spread':>16}{'global_sync':>13}{'R_bimod':>9}")
rows = {}
for nm, fn in [("uniform_sync", uniform_sync), ("incoherent", incoherent),
               ("synthetic_chimera", chimera), ("as_chimera", as_chimera)]:
    r = sync_texture(H(fn))
    rows[nm] = r
    print(f"{nm:<16}{r['local_R_spread']:>16.4f}{r['global_sync']:>13.4f}{r['local_R_bimod']:>9.3f}")

chi = max(rows["synthetic_chimera"]["local_R_spread"], rows["as_chimera"]["local_R_spread"])
nulls = max(rows["uniform_sync"]["local_R_spread"], rows["incoherent"]["local_R_spread"])
print(f"\nchimera spread {chi:.4f} vs max uniform-state spread {nulls:.4f}")
print(f"VERDICT: {'ADMIT' if chi > 2 * nulls + 0.03 else 'review'}")
