"""Does velocity_order earn its place? The novel coordinate is velcorr_length xi. Synthetic
BARE-POSITION streams: correlated_flock (smooth velocity field -> high xi), drifting_gas (high
polarization but INDEPENDENT fluctuations -> xi~0, the confound polarization can't catch),
still_gas (independent jitter), milling (rotation). Want: xi HIGH for correlated motion, ~0 for
both gases; polarization HIGH for flock+drifting_gas (so xi is what tells them apart); mill HIGH
for milling."""
import numpy as np

from epc.metrics.velocity_order import velocity_order

N = 400
rng = np.random.default_rng(0)
P0 = rng.uniform(0, 40, (N, 2))


def stream(step_fn, T=10, seed=0):
    r = np.random.default_rng(seed)
    pos = P0.copy()
    out = [{"positions": pos.copy()}]
    for t in range(T):
        pos = pos + step_fn(pos, r)
        out.append({"positions": pos.copy()})
    return out


def correlated_flock(pos, r):
    # smooth spatial velocity field (correlated fluctuations) + mean drift
    k = 2 * np.pi / 40
    vx = 0.5 + 0.4 * np.sin(k * pos[:, 1])
    vy = 0.4 * np.cos(k * pos[:, 0])
    return np.stack([vx, vy], 1) + 0.02 * r.standard_normal((N, 2))


def drifting_gas(pos, r):
    # SAME mean drift as a flock, but INDEPENDENT fluctuations (no spatial correlation)
    return np.array([0.5, 0.0]) + 0.4 * r.standard_normal((N, 2))


def still_gas(pos, r):
    return 0.4 * r.standard_normal((N, 2))


def milling(pos, r):
    c = pos - pos.mean(0)
    perp = np.stack([-c[:, 1], c[:, 0]], 1)
    perp /= (np.hypot(c[:, 0], c[:, 1])[:, None] + 1e-9)
    return 0.5 * perp + 0.02 * r.standard_normal((N, 2))


print(f"{'system':<18}{'velcorr_xi':>12}{'polarization':>14}{'mill':>8}")
rows = {}
for nm, fn in [("correlated_flock", correlated_flock), ("drifting_gas", drifting_gas),
               ("still_gas", still_gas), ("milling", milling)]:
    r = velocity_order(stream(fn, seed=1))
    rows[nm] = r
    print(f"{nm:<18}{r['velcorr_length']:>12.4f}{r['polarization']:>14.4f}{r['mill_strength']:>8.4f}")

xi_flock = rows["correlated_flock"]["velcorr_length"]
xi_gas = max(rows["drifting_gas"]["velcorr_length"], rows["still_gas"]["velcorr_length"])
drift_pol = rows["drifting_gas"]["polarization"]
print(f"\nxi: correlated_flock {xi_flock:.4f} vs max gas {xi_gas:.4f}  (drifting_gas polarization {drift_pol:.3f} = high yet xi low -> confound caught)")
print(f"VERDICT: {'ADMIT' if xi_flock > 2 * xi_gas + 0.05 else 'review'}")
