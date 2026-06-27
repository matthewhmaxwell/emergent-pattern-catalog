"""Finite-size-scaling criticality HARNESS (deferred-item resolution, proper method).

Why this is a harness, not a battery lens: criticality cannot be certified from a single
observation at one system size -- an ordered phase and a critical point both show long
correlations; only how observables SCALE WITH L distinguishes them. So this runs a model at
several sizes L and applies the gold-standard FSS criticality tests:

  - Binder cumulant U4(T,L) = 1 - <m^4>/(3 <m^2>^2): curves for different L CROSS at Tc
    (L-independent crossing = the criticality certificate).
  - susceptibility chi(T,L) = N (<m^2> - <|m|>^2): peak near Tc whose HEIGHT grows with L.

Demonstrated on the 2D Ising model (exact Tc = 2/ln(1+sqrt2) ~ 2.269), vectorized checkerboard
Metropolis. This is the depth-correct closure for the correlation-length/criticality axis: the
per-observation battery keeps spatial-correlation coordinates (Moran's I, structure factor); the
criticality QUESTION is answered here, across sizes.
"""
import numpy as np

TC = 2.0 / np.log(1.0 + np.sqrt(2.0))


def _ising(L, T, sweeps=240, burn=120, seed=0):
    r = np.random.default_rng(seed)
    s = r.choice([-1, 1], size=(L, L)).astype(np.int8)
    idx = np.indices((L, L)).sum(0) % 2
    masks = [idx == 0, idx == 1]
    mags = []
    for t in range(sweeps + burn):
        for mk in masks:
            nb = (np.roll(s, 1, 0) + np.roll(s, -1, 0) + np.roll(s, 1, 1) + np.roll(s, -1, 1))
            dE = 2 * s * nb
            flip = mk & (r.random((L, L)) < np.exp(-dE / T))
            s[flip] *= -1
        if t >= burn:
            mags.append(abs(s.mean()))
    m = np.asarray(mags)
    m2 = np.mean(m ** 2); m4 = np.mean(m ** 4)
    U4 = 1.0 - m4 / (3.0 * m2 ** 2 + 1e-12)
    chi = (L * L) * (m2 - np.mean(m) ** 2)
    return float(U4), float(chi)


def main():
    Ls = [12, 20, 32]
    Ts = [2.0, 2.15, 2.27, 2.4, 2.6]
    print(f"exact Tc = {TC:.4f}\n")
    print(f"{'T':>6} | " + " ".join(f"U4(L={L})" for L in Ls) + " | " + " ".join(f"chi(L={L})" for L in Ls))
    U = {L: [] for L in Ls}; C = {L: [] for L in Ls}
    for T in Ts:
        row = []
        for L in Ls:
            u, c = _ising(L, T, seed=1)
            U[L].append(u); C[L].append(c); row.append((u, c))
        print(f"{T:>6.2f} | " + " ".join(f"{u:7.3f}" for u, _ in row) + " | " + " ".join(f"{c:8.2f}" for _, c in row))
    # criticality certificate: chi peak grows with L, and peak T near Tc
    print("\nchi-peak location + height vs L (FSS criticality signature):")
    for L in Ls:
        i = int(np.argmax(C[L]))
        print(f"  L={L:>3}: chi_max={C[L][i]:8.2f} at T={Ts[i]:.2f}")
    grows = C[Ls[-1]][int(np.argmax(C[Ls[-1]]))] > C[Ls[0]][int(np.argmax(C[Ls[0]]))]
    peak_T = Ts[int(np.argmax(C[Ls[-1]]))]
    near_tc = abs(peak_T - TC) <= 0.2
    print(f"\nchi_max grows with L: {grows} | largest-L peak T={peak_T:.2f} near exact Tc {TC:.3f}: {near_tc}")
    print(f"VERDICT: {'FSS criticality CERTIFIED (chi-peak scales with L, peaks at Tc)' if grows and near_tc else 'review'}")


if __name__ == "__main__":
    main()
