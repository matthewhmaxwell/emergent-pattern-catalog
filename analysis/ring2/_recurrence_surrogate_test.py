"""Can the deferred recurrence lens be rescued? The confound: a smooth random walk scores
high DET (null_walk 0.54) without genuinely revisiting states. Test 3 candidate
discriminators that should separate genuine recurrence (limit cycle / chaos / patterned —
revisit FAR-in-time states) from a random walk (smooth but drifts, no long-range revisits):
  (1) FT-surrogate excess: DET_obs - mean(DET of phase-randomized surrogates).
  (2) shuffle-surrogate excess: DET_obs - mean(DET of time-shuffled surrogates).
  (3) long-range recurrence fraction: fraction of recurrence points with |i-j| > T/3
      (periodic revisits span the series; a random walk's recurrences are short-range).
Want a measure where null_walk is LOW and limit_cycle/chaos/patterned are HIGH."""
import numpy as np
from analysis.blind_spot_probes import PROBES
from epc.metrics.recurrence import _invariant_trajectory, _line_fraction

NAMES = {"null_walk", "null_noise", "limit_cycle", "traveling_wave",
         "spatiotemporal_chaos", "flocking", "vortex_milling", "aggregation"}


def build_R(traj, target_rr=0.10, theiler=None):
    T = traj.shape[0]
    if theiler is None:
        theiler = max(2, T // 12)
    D = np.sqrt(((traj[:, None, :] - traj[None, :, :]) ** 2).sum(-1))
    off = D[~np.eye(T, dtype=bool)]
    eps = float(np.quantile(off, target_rr)) if off.size else 0.0
    R = (D <= eps).astype(np.int8)
    ii, jj = np.indices((T, T))
    R[np.abs(ii - jj) <= theiler] = 0
    return R, np.abs(ii - jj)


def det_of(traj):
    R, _ = build_R(traj)
    return _line_fraction(R, 2, vertical=False)


def ft_surrogate(traj, rng):
    T, d = traj.shape
    out = np.empty_like(traj)
    for j in range(d):
        Xf = np.fft.rfft(traj[:, j]); mag = np.abs(Xf)
        ph = np.exp(1j * rng.uniform(0, 2 * np.pi, size=Xf.shape)); ph[0] = 1.0
        if T % 2 == 0:
            ph[-1] = 1.0
        out[:, j] = np.fft.irfft(mag * ph, n=T)
    return out


def longrange_frac(traj):
    R, sep = build_R(traj)
    s = sep[R > 0]
    return float((s > traj.shape[0] / 3).mean()) if s.size else 0.0


rng = np.random.default_rng(0)
rows = []
for p in PROBES:
    if p.__name__ not in NAMES:
        continue
    pr = p(0)
    traj = _invariant_trajectory(pr["history"])
    if traj is None:
        rows.append((p.__name__, pr.get("truth"), None)); continue
    det = det_of(traj)
    ft = np.mean([det_of(ft_surrogate(traj, rng)) for _ in range(12)])
    sh = []
    for _ in range(12):
        t2 = traj.copy(); rng.shuffle(t2); sh.append(det_of(t2))
    sh = np.mean(sh)
    rows.append((p.__name__, pr.get("truth"),
                 {"DET": det, "ft_exc": det - ft, "sh_exc": det - sh, "lr": longrange_frac(traj)}))

print(f"{'system':<22}{'truth':<10}{'DET':>7}{'ft_exc':>8}{'sh_exc':>8}{'longrange':>10}")
for nm, t, r in sorted(rows, key=lambda x: -((x[2] or {}).get('lr', -1))):
    if r:
        print(f"{nm:<22}{str(t):<10}{r['DET']:>7.3f}{r['ft_exc']:>8.3f}{r['sh_exc']:>8.3f}{r['lr']:>10.3f}")
    else:
        print(f"{nm:<22}{str(t):<10}   (n/a)")

for key in ("DET", "ft_exc", "sh_exc", "lr"):
    nul = [r[key] for nm, t, r in rows if r and t == "null"]
    emg = [r[key] for nm, t, r in rows if r and t == "emergent"]
    if nul and emg:
        print(f"{key:>9}:  null max {max(nul):+.3f}   emergent min {min(emg):+.3f}   "
              f"{'SEPARATES' if min(emg) > max(nul) else 'overlaps'}")
