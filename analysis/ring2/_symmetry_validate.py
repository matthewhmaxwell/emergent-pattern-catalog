"""Does the pattern-symmetry lens earn its place? Synthetic fields that share a |k| but differ
in azimuthal symmetry: stripes (m=2), square (m=4), hexagonal (m=6), labyrinth (isotropic ring,
low strength), gas (no ring). Want: lattice_symmetry recovers 2/4/6; symmetry_strength HIGH for
lattices, LOW for labyrinth; ring_contrast HIGH for patterned, ~1 for gas."""
import numpy as np
from epc.metrics.pattern_symmetry import pattern_symmetry

G = 96
k = 8.0
y, x = np.indices((G, G)) * (2 * np.pi / G)


def H(fld_fn, n=6):
    return [{"field": fld_fn(t)} for t in range(n)]


def stripes(t):
    return np.cos(k * x)


def square(t):
    return np.cos(k * x) + np.cos(k * y)


def hexagonal(t):
    a1 = np.cos(k * x)
    a2 = np.cos(k * (0.5 * x + 0.8660254 * y))
    a3 = np.cos(k * (0.5 * x - 0.8660254 * y))
    return a1 + a2 + a3


def labyrinth(t):
    rng = np.random.default_rng(100 + t)
    f = np.zeros((G, G), complex)
    KX = np.fft.fftshift(np.fft.fftfreq(G) * G)
    KXX, KYY = np.meshgrid(KX, KX)
    rr = np.hypot(KXX, KYY)
    band = (rr > k - 1) & (rr < k + 1)
    ph = rng.uniform(0, 2 * np.pi, (G, G))
    f[band] = np.exp(1j * ph[band])
    return np.real(np.fft.ifft2(np.fft.ifftshift(f)))


def gas(t):
    return np.random.default_rng(t).standard_normal((G, G))


print(f"{'system':<14}{'symmetry(m)':>12}{'concentration':>14}{'ring_contrast':>15}")
rows = {}
for nm, fn in [("stripes", stripes), ("square", square), ("hexagonal", hexagonal),
               ("labyrinth", labyrinth), ("gas", gas)]:
    r = pattern_symmetry(H(fn))
    rows[nm] = r
    print(f"{nm:<14}{r['lattice_symmetry']:>12}{r['azimuthal_concentration']:>14.3f}{r['ring_contrast']:>15.2f}")

ok_sym = rows["stripes"]["lattice_symmetry"] == 2 and rows["square"]["lattice_symmetry"] == 4 and rows["hexagonal"]["lattice_symmetry"] == 6
lat_c = min(rows[n]["azimuthal_concentration"] for n in ("stripes", "square", "hexagonal"))
lab_c = rows["labyrinth"]["azimuthal_concentration"]
gas_ring = rows["gas"]["ring_contrast"]
print(f"\nsymmetry recovered 2/4/6: {ok_sym} | min lattice concentration {lat_c:.3f} vs labyrinth {lab_c:.3f} | gas ring_contrast {gas_ring:.2f}")
print(f"VERDICT: {'ADMIT' if ok_sym and lat_c > 2 * lab_c and gas_ring < 3 else 'review'}")
