"""Does the persistent-homology lens earn its place? Positive control: a noisy RING
should show one big persistent loop; a Gaussian BLOB should show ~0. Then the probe
corpus, leads, nulls. DISCRIMINATOR is h1_max (the single dominant loop) NOT h1_total
(sum is noise-confounded: dense random clouds have many tiny spurious holes). Want:
genuine loops/voids HIGH h1_max; blobs / gas / nulls / lattices LOW h1_max."""
import numpy as np
from analysis.blind_spot_probes import PROBES
from epc.metrics.persistent_homology import persistent_homology
from epc.models.particle_life import particle_life


def _hist(pts_fn, frames=18):
    out = []
    for t in range(frames):
        p = pts_fn(t)
        out.append({"positions": p, "velocities": np.zeros_like(p),
                    "types": np.zeros(len(p)), "step": t})
    return out


def ring(seed=0, n=200, R=4.0, noise=0.3):
    rng = np.random.default_rng(seed)
    return _hist(lambda t: (lambda th, r: np.c_[5 + r * np.cos(th), 5 + r * np.sin(th)])(
        rng.uniform(0, 2 * np.pi, n), R + rng.normal(0, noise, n)))


def blob(seed=0, n=200):
    rng = np.random.default_rng(seed)
    return _hist(lambda t: rng.normal(5, 1.5, (n, 2)))


rows = [("SYNTH_ring", "ctrl-loop", persistent_homology(ring())),
        ("SYNTH_blob", "ctrl-blob", persistent_homology(blob()))]
for probe in PROBES:
    try:
        p = probe(0)
        rows.append((probe.__name__, p.get("truth"), persistent_homology(p["history"])))
    except Exception:
        rows.append((probe.__name__, "err", None))
for i in [180, 270, 288, 298, 205]:
    rng = np.random.default_rng(i); K = int(rng.integers(3, 7)); F = rng.uniform(-1.0, 1.0, (K, K))
    h, _ = particle_life(i, K=K, F=F, N=200, steps=250)
    rows.append((f"pl_{i}", ("lead" if i != 205 else "pl-class"), persistent_homology(h)))

print(f"{'system':<22}{'truth':<12}{'h1_max':>8}{'h1_total':>9}{'ncomp':>7}  <- sorted by h1_max")
for nm, t, r in sorted(rows, key=lambda x: -((x[2] or {}).get('h1_max', -1))):
    if r:
        print(f"{nm:<22}{str(t):<12}{r['h1_max']:>8.3f}{r['h1_total']:>9.3f}{r['n_components']:>7.1f}")
    else:
        print(f"{nm:<22}{str(t):<12}   (n/a)")

# verdict: do genuine loops (ring ctrl, vortex) separate from everything non-loop?
applic = [(nm, t, r['h1_max']) for nm, t, r in rows if r]
loops = [v for nm, t, v in applic if nm in ("SYNTH_ring", "vortex_milling")]
nonloop = [v for nm, t, v in applic if nm not in ("SYNTH_ring", "vortex_milling")]
if loops and nonloop:
    gap = min(loops) - max(nonloop)
    print(f"\nh1_max  loops(ring,vortex): {[round(v,3) for v in sorted(loops)]}")
    print(f"h1_max  max non-loop:      {round(max(nonloop),3)}  (gap = {gap:+.3f})")
    print(f"VERDICT: {'ADMIT (clean gap)' if gap > 0.15 else 'DEFER (no clean gap)'}"
          f"  — threshold ~{(min(loops)+max(nonloop))/2:.2f}")
