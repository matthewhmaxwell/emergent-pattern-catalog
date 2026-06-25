"""Diagnose two RD hardening findings:
 (A) why structure_factor reads low on RD fields -> prototype a resolution-robust
     RADIAL S(k) peak (azimuthal average, the standard condensed-matter measure);
 (B) did 'death/null' actually die? print field stats."""
import numpy as np
from epc.models.reaction_diffusion import reaction_diffusion
from epc.metrics.structure_factor import _sk_peak

REGIMES = [("solitons",0.022,0.051),("mitosis",0.026,0.051),("stripes",0.030,0.057),
           ("spots",0.034,0.063),("worms",0.039,0.058),("death/null",0.060,0.062),
           ("low-F-death",0.010,0.050)]


def radial_sk_peak(field):
    a = field - field.mean()
    if a.std() < 1e-12:
        return 1.0
    S = np.abs(np.fft.fftshift(np.fft.fft2(a)))**2
    ny, nx = S.shape
    cy, cx = ny//2, nx//2
    S[cy, cx] = 0.0
    y, x = np.indices(S.shape)
    r = np.sqrt((y-cy)**2 + (x-cx)**2).astype(int)
    rmax = min(cy, cx)
    prof = np.array([S[r==rr].mean() for rr in range(1, rmax)])
    if prof.size < 3 or np.median(prof) <= 0:
        return 1.0
    return float(prof.max() / np.median(prof)), int(np.argmax(prof)+1)


print(f"{'regime':<12}{'F':>6}{'k':>7}{'v.mean':>8}{'v.std':>7}{'v.max':>7}{'old_sk':>8}{'radSk':>8}{'k0':>4}")
for name, F, k in REGIMES:
    h = reaction_diffusion(0, F=F, k=k, N=96, steps=8000, record=24)
    v = h[-1]["field"]
    old = _sk_peak(v)
    rad, k0 = radial_sk_peak(v)
    print(f"{name:<12}{F:>6.3f}{k:>7.3f}{v.mean():>8.3f}{v.std():>7.3f}{v.max():>7.3f}{old:>8.1f}{rad:>8.1f}{k0:>4d}")
print("\n(A) radSk = radial S(k) peak/baseline; k0 = peak radius. Turing patterns should show"
      " a sharp radial peak (radSk high) at a characteristic k0>0. (B) v.std~0 => died (true null);"
      " v.std>0 => still structured (em fire is correct, label was wrong).")
