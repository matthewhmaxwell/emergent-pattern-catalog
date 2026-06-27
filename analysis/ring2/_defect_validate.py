"""Does the topological-defect-census lens earn its place? The active_nematic family gives a
ready positive (a +-1/2 defect SEA) and contrasts (polar flock = aligned/0, milling vortex =
one core, isotropic noise = random). Plus synthetic polar fields: a vortex pair, uniform, iid
noise. Want: defect SEA HIGH density; uniform/aligned ~0; noise ~0 (annihilates under
smoothing); single vortex tiny."""
import numpy as np
from epc.metrics.defect_census import defect_census
from epc.models.active_nematic import active_nematic_field, neg_polar_flock, neg_milling, neg_isotropic

N = 64


def H(theta_fn, n=10, polar=True):
    out = []
    for t in range(n):
        th = theta_fn(t)
        out.append({"phases": th} if polar else {"theta_field": th})
    return out


def vortex_pair(t):
    y, x = np.indices((N, N))
    return (np.arctan2(y - N * 0.35, x - N * 0.35) - np.arctan2(y - N * 0.65, x - N * 0.65)) % (2 * np.pi)


def uniform(t):
    return np.full((N, N), 0.7)


def noise(t):
    return np.random.default_rng(t).uniform(0, 2 * np.pi, (N, N))


rows = []
for nm, fn in [("active_nematic", active_nematic_field), ("polar_flock", neg_polar_flock),
               ("milling_vortex", neg_milling), ("isotropic_noise", neg_isotropic)]:
    h, _ = fn(0, G=64, n_frames=8)
    rows.append((nm, defect_census(h)))
rows.append(("SYNTH_vortex_pair", defect_census(H(vortex_pair))))
rows.append(("SYNTH_uniform", defect_census(H(uniform))))
rows.append(("SYNTH_noise", defect_census(H(noise))))

print(f"{'system':<20}{'defect_density':>15}{'coherence':>11}{'net_charge':>12}")
for nm, r in sorted(rows, key=lambda x: -((x[1] or {}).get('defect_density', -1))):
    print(f"{nm:<20}{r['defect_density']:>15.4f}{r['coherence']:>11.3f}{r['net_charge']:>12.3f}" if r else f"{nm:<20}      (n/a)")

val = {nm: r for nm, r in rows if r}
sea = val.get("active_nematic", {}).get("defect_density", 0)
nulls = [val[n]["defect_density"] for n in ("polar_flock", "SYNTH_uniform", "SYNTH_noise", "isotropic_noise") if n in val]
if sea and nulls:
    gap = sea - max(nulls)
    print(f"\ndefect SEA (active_nematic) {sea:.4f}  vs  max null(flock/uniform/noise) {max(nulls):.4f}  gap {gap:+.4f}")
    print(f"VERDICT: {'ADMIT (defect sea separates from ordered/noise)' if gap > 0.01 else 'review'}")
