"""Does the vorticity / rotational-order lens earn its place? The corpus has a ready
positive (vortex_milling) and a ready null (flocking = translation). Plus synthetic
controls: a rotating swarm (high |L|), a translating swarm (|L|~0), a random-velocity
swarm (|L|~0). Want: rotation HIGH ang_mom; translation/random/other LOW."""
import numpy as np
from analysis.blind_spot_probes import PROBES
from epc.metrics.vorticity import rotational_order


def H(pos_fn, vel_fn, n=12, N=200):
    out = []
    for t in range(n):
        out.append({"positions": pos_fn(t, N), "velocities": vel_fn(t, N)})
    return out


def rot(t, N):
    rng = np.random.default_rng(N); th = rng.uniform(0, 2 * np.pi, N); r = rng.uniform(2, 8, N)
    a = t * 0.1
    return np.c_[r * np.cos(th + a), r * np.sin(th + a)]


def rot_v(t, N):
    p = rot(t, N); return np.c_[-p[:, 1], p[:, 0]]                      # tangential -> rotation


def trans(t, N):
    rng = np.random.default_rng(N); return rng.uniform(-8, 8, (N, 2)) + np.array([t, 0.0])


def trans_v(t, N):
    return np.tile([1.0, 0.0], (N, 1))                                 # all same direction


def randv(t, N):
    rng = np.random.default_rng(1000 + t); return rng.standard_normal((N, 2))


rows = [("SYNTH_rotating", "rotation", rotational_order(H(rot, rot_v))),
        ("SYNTH_translating", "translation", rotational_order(H(trans, trans_v))),
        ("SYNTH_random_vel", "random", rotational_order(H(trans, randv)))]
for p in PROBES:
    try:
        pr = p(0)
        rows.append((p.__name__, pr.get("truth"), rotational_order(pr["history"])))
    except Exception as e:
        rows.append((p.__name__, "err", None))

print(f"{'system':<22}{'truth':<12}{'ang_mom':>9}")
for nm, t, r in sorted(rows, key=lambda x: -((x[2] or {}).get('ang_mom', -1))):
    print(f"{nm:<22}{str(t):<12}{r['ang_mom']:>9.3f}" if r else f"{nm:<22}{str(t):<12}   (n/a)")

val = {nm: r for nm, t, r in rows if r}
rotors = [val[n]['ang_mom'] for n in ('SYNTH_rotating', 'vortex_milling') if n in val]
others = [r['ang_mom'] for nm, t, r in rows if r and nm not in ('SYNTH_rotating', 'vortex_milling')]
if rotors and others:
    gap = min(rotors) - max(others)
    print(f"\nrotation(synth+vortex_milling): {[round(v,2) for v in sorted(rotors)]}  "
          f"max non-rotor: {round(max(others),2)}  gap {gap:+.3f}")
    print(f"VERDICT: {'ADMIT (rotation separates)' if gap > 0.1 else 'review'}")
