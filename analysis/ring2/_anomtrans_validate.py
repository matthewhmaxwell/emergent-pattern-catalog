"""Does anomalous_transport earn its place? Synthetic trajectory ensembles with KNOWN classes:
brownian (beta~1, EB~0), ballistic (beta~2), CTRW with heavy-tailed waits (beta<1 AND EB>0 =
ergodicity breaking). Want: msd_exponent recovers 1 / 2 / <1; ergodicity_break separates CTRW
from brownian."""
import numpy as np
from epc.metrics.anomalous_transport import anomalous_transport

T, N = 60, 300


def frames(X):
    return [{"positions": X[t]} for t in range(X.shape[0])]


def brownian(seed=0):
    r = np.random.default_rng(seed)
    steps = r.standard_normal((T, N, 2))
    X = np.cumsum(steps, axis=0)
    return frames(X)


def ballistic(seed=0):
    r = np.random.default_rng(seed)
    v = r.standard_normal((N, 2))
    t = np.arange(T)[:, None, None]
    X = t * v[None] + 0.05 * r.standard_normal((T, N, 2))
    return frames(X)


def fbm_sub(seed=0, H=0.3):
    # fractional Brownian motion via spectral synthesis: ergodic, MSD ~ t^(2H) -> beta=2H<1
    r = np.random.default_rng(seed)
    f = np.fft.rfftfreq(T).copy(); f[0] = f[1]
    psd = f ** (-(2 * H + 1))
    X = np.zeros((T, N, 2))
    for i in range(N):
        for d in range(2):
            ph = r.uniform(0, 2 * np.pi, len(f))
            x = np.fft.irfft(np.sqrt(psd) * np.exp(1j * ph), T)
            X[:, i, d] = x - x[0]
    return frames(X)


def ctrw(seed=0):
    r = np.random.default_rng(seed)
    X = np.zeros((T, N, 2))
    for i in range(N):
        t, pos = 0, np.zeros(2)
        while t < T:
            wait = 1 + int(r.pareto(0.7))            # heavy-tailed waiting -> subdiffusion + EB
            jump = r.standard_normal(2)
            te = min(T, t + wait)
            X[t:te, i] = pos
            pos = pos + jump
            t = te
    return frames(X)


print(f"{'system':<14}{'msd_exponent':>14}{'ergodicity_break':>18}")
rows = {}
for nm, fn in [("brownian", brownian), ("ballistic", ballistic), ("fbm_sub", fbm_sub), ("ctrw", ctrw)]:
    r = anomalous_transport(fn(1))
    rows[nm] = r
    print(f"{nm:<14}{r['msd_exponent']:>14.3f}{r['ergodicity_break']:>18.3f}")

b = rows["brownian"]; ba = rows["ballistic"]; fb = rows["fbm_sub"]; c = rows["ctrw"]
beta_ok = abs(b["msd_exponent"] - 1) < 0.25 and ba["msd_exponent"] > 1.6 and fb["msd_exponent"] < 0.9
eb_ok = c["ergodicity_break"] > 2 * b["ergodicity_break"] + 0.05
print(f"\nbeta: brownian~1 {b['msd_exponent']:.2f}, ballistic~2 {ba['msd_exponent']:.2f}, fbm-sub<1 {fb['msd_exponent']:.2f} | EB: ctrw {c['ergodicity_break']:.2f} vs brownian {b['ergodicity_break']:.2f}")
print(f"VERDICT: {'ADMIT' if beta_ok and eb_ok else 'review'}")
