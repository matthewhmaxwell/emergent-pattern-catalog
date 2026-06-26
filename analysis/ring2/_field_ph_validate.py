"""Does the field-topology PH path earn its place? Synthetic controls with KNOWN hole
counts (Swiss-cheese = many, ring = 1, blob = 0, noise = spurious) + the real field
substrates (RD worms/spots/stripes, Lenia). Want: holey fields HIGH field_loops;
blobs/spots/smooth/noise LOW (size filter suppresses noise pockets)."""
import numpy as np
from epc.metrics.persistent_homology import persistent_homology
from epc.models.reaction_diffusion import reaction_diffusion
from epc.models.lenia import lenia

N = 96


def H(field, n=4):
    return [{"field": np.asarray(field, float)} for _ in range(n)]


def swiss_cheese(holes=8, r=6, seed=0):
    rng = np.random.default_rng(seed); a = np.ones((N, N))
    y, x = np.indices((N, N))
    for _ in range(holes):
        cy, cx = rng.integers(r + 2, N - r - 2, size=2)
        a[(y - cy) ** 2 + (x - cx) ** 2 < r ** 2] = 0.0
    return a


def ring(R=30, w=6):
    y, x = np.indices((N, N)); d = np.hypot(y - N / 2, x - N / 2)
    return ((d > R - w) & (d < R + w)).astype(float)


def blob(R=30):
    y, x = np.indices((N, N)); return (np.hypot(y - N / 2, x - N / 2) < R).astype(float)


def noise(seed=0):
    return np.random.default_rng(seed).random((N, N))


rows = [
    ("SYNTH_swiss_cheese", "holey", persistent_homology(H(swiss_cheese()))),
    ("SYNTH_ring", "holey", persistent_homology(H(ring()))),
    ("SYNTH_blob", "solid", persistent_homology(H(blob()))),
    ("SYNTH_noise", "noise", persistent_homology(H(noise()))),
    ("RD_worms", "holey", persistent_homology(reaction_diffusion(0, F=0.039, k=0.058, N=96, steps=8000, record=12))),
    ("RD_spots", "solid", persistent_homology(reaction_diffusion(0, F=0.034, k=0.063, N=96, steps=8000, record=12))),
    ("RD_stripes", "?", persistent_homology(reaction_diffusion(0, F=0.030, k=0.057, N=96, steps=8000, record=12))),
    ("lenia_glider", "?", persistent_homology(lenia(0, N=64, steps=400, record=12, n_creatures=1))),
]

print(f"{'field':<22}{'truth':<8}{'field_loops':>12}{'loop_area':>11}")
for nm, t, r in sorted(rows, key=lambda x: -((x[2] or {}).get('field_loops', -1))):
    if r and r.get("kind") == "field":
        print(f"{nm:<22}{t:<8}{r['field_loops']:>12.2f}{r['field_loop_area']:>11.4f}")
    else:
        print(f"{nm:<22}{t:<8}   (n/a: {r})")

# loop_AREA is the discriminator (count is noise-confounded; real voids are LARGE)
holey = [r['field_loop_area'] for nm, t, r in rows if r and t == 'holey']
nonh = [r['field_loop_area'] for nm, t, r in rows if r and t in ('solid', 'noise')]
if holey and nonh:
    gap = min(holey) - max(nonh)
    print(f"\nholey loop_area: {[round(x,3) for x in sorted(holey)]}  "
          f"solid/noise: {[round(x,3) for x in sorted(nonh)]}")
    print(f"gap (min holey - max non-holey) = {gap:+.3f}  -> "
          f"{'ADMIT field path (area separates)' if gap > 0.02 else 'review'}")
