"""Exact infinite-lattice interface density rho(t) of the 2D voter model, random m=0 start.

Theory (no voter-model simulation): Frachebourg & Krapivsky (1996) PRE 53:R3009 Eqs. (9)-(11) /
coalescing-random-walk duality (Dornic et al. 2001 PRL 87:045701, text before Eq. 1):
  rho(t) = (1/2) * S_e(t), S_e(t) = probability that the difference of two lineages, started at
  bond vector e, has not hit the origin by time t. Each lineage jumps at rate 1 to a uniformly chosen
  neighbour, so the difference walk has rate 2. S_e is computed exactly for the discrete-step walk by
  propagating the distribution with an absorbing origin, then Poissonised with mean 2t.
Time unit: each site updated once per unit time on average = one Monte Carlo sweep.
Moore value averages the n.n. and diagonal bond directions (matches voter.py _wall_density).
Runtime ~80 s on one core.
"""
import numpy as np
from scipy.stats import poisson

def surv_discrete(offsets, start, nmax, R):
    n = 2*R+1; p = np.zeros((n, n)); p[R+start[0], R+start[1]] = 1.0
    S = np.empty(nmax+1); S[0] = 1.0; k = len(offsets)
    for s in range(1, nmax+1):
        q = np.zeros_like(p)
        for dr, dc in offsets:
            q += np.roll(np.roll(p, dr, 0), dc, 1)
        p = q/k; p[R, R] = 0.0
        S[s] = p.sum()
    return S

VN = [(-1, 0), (1, 0), (0, -1), (0, 1)]
MO = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]
ts = [100, 300, 500, 1000, 2000, 3000]
tmax = max(ts); nmax = int(2*tmax + 8*np.sqrt(2*tmax)); R = int(6*np.sqrt(nmax)) + 5
res = {}
for name, offs, starts in [('vonNeumann', VN, [(1, 0)]), ('Moore', MO, [(1, 0), (1, 1)])]:
    Ss = [surv_discrete(offs, s, nmax, R) for s in starts]
    n = np.arange(nmax+1)
    for t in ts:
        w = poisson.pmf(n, 2*t)
        res[(name, t)] = 0.5*np.mean([float((w*S).sum()) for S in Ss])
g = 0.5772156649
for t in ts:
    fk = np.pi/(2*np.log(t) + np.log(256)); fkg = np.pi/(2*np.log(t) + np.log(256) + 2*g)
    print("t=%5d exact_vN=%.5f FK_Eq23=%.5f (%+.1f%%) FK+2gamma=%.5f exact_Moore=%.5f"
          % (t, res[('vonNeumann', t)], fk, 100*(fk/res[('vonNeumann', t)]-1), fkg, res[('Moore', t)]))
fit_t = [300, 500, 1000, 2000, 3000]
for name, asym in [('vonNeumann', 2/np.pi), ('Moore', 4/(3*np.pi))]:
    sl = np.polyfit(np.log(fit_t), [1/res[(name, t)] for t in fit_t], 1)[0]
    print("%s slope d(1/rho)/dln t over t=300..3000: %.4f; asymptotic 1/(2 pi D) = %.4f" % (name, sl, asym))
