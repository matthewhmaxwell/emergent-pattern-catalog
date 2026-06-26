"""Does the correlation-length / criticality lens earn its place? Controls with KNOWN
correlation structure: white noise (xi~0, no power law), exponential-correlated (finite xi,
exponential decay), power-law/critical (large xi, power-law decay -> plaw_gain>0), smooth
low-frequency order (xi large). Then the RD probe (Turing characteristic scale) + a couple
corpus systems. Want: plaw_gain separates power-law/critical fields from exponential+white;
xi_norm separates long-range from white noise."""
import numpy as np
from epc.metrics.correlation_length import correlation_length
from epc.models.reaction_diffusion import reaction_diffusion

N = 96


def H(field, n=4):
    return [{"field": np.asarray(field, float)} for _ in range(n)]


def white(seed=0):
    return np.random.default_rng(seed).standard_normal((N, N))


def exp_corr(seed=0, sigma=4.0):
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(np.random.default_rng(seed).standard_normal((N, N)), sigma)


def power_law(seed=0, beta=3.0):
    # field with power-law power spectrum P(k) ~ k^-beta -> power-law (scale-free) correlations
    rng = np.random.default_rng(seed)
    ky, kx = np.indices((N, N)) - N // 2
    k = np.sqrt(ky ** 2 + kx ** 2); k[N // 2, N // 2] = 1.0
    amp = k ** (-beta / 2.0)
    ph = np.exp(2j * np.pi * rng.random((N, N)))
    fld = np.fft.ifft2(np.fft.ifftshift(amp * ph)).real
    return fld


def smooth(seed=0):
    y, x = np.indices((N, N))
    return np.sin(2 * np.pi * x / N) + np.cos(2 * np.pi * y / N)


rows = [("white_noise", "disordered", correlation_length(H(white()))),
        ("exp_corr(sig4)", "short-range", correlation_length(H(exp_corr()))),
        ("power_law(b3)", "critical", correlation_length(H(power_law(0, 3.0)))),
        ("power_law(b2.5)", "critical", correlation_length(H(power_law(0, 2.5)))),
        ("smooth_order", "long-range", correlation_length(H(smooth()))),
        ("RD_stripes", "turing", correlation_length(reaction_diffusion(0, F=0.030, k=0.057, N=96, steps=8000, record=12)))]

print(f"{'system':<18}{'truth':<13}{'xi_norm':>9}{'plaw_gain':>11}")
for nm, t, r in sorted(rows, key=lambda x: -((x[2] or {}).get('plaw_gain', -9))):
    if r:
        print(f"{nm:<18}{t:<13}{r['xi_norm']:>9.3f}{r['plaw_gain']:>11.3f}")
    else:
        print(f"{nm:<18}{t:<13}   (n/a)")

val = {nm: r for nm, t, r in rows if r}
crit = [val[n]['plaw_gain'] for n in ('power_law(b3)', 'power_law(b2.5)') if n in val]
noncrit = [val[n]['plaw_gain'] for n in ('white_noise', 'exp_corr(sig4)') if n in val]
if crit and noncrit:
    gap = min(crit) - max(noncrit)
    print(f"\nplaw_gain  critical(power-law): {[round(v,2) for v in sorted(crit)]}  "
          f"exp/white: {[round(v,2) for v in sorted(noncrit)]}  gap {gap:+.3f}")
    print(f"xi_norm    white {val.get('white_noise',{}).get('xi_norm')} vs long-range "
          f"{val.get('smooth_order',{}).get('xi_norm')}")
    print(f"VERDICT: {'ADMIT (plaw_gain separates critical from exp/white)' if gap > 0.1 else 'review/defer'}")
