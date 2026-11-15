"""Replica-symmetric (RS) predictions for the S=2 minority game, alpha > alpha_c.

Sources of the equations (no simulation involved):
  Challet, Marsili & Zecchina (2000) PRL 84:1824, Eq. (9): P(m), z = sqrt(alpha/(1+Q)),
      per-sign frozen weight (1-erf(z/sqrt2))/2  -> total frozen fraction phi = 1 - erf(z/sqrt2).
  Marsili, Challet & Zecchina (2000) Physica A 280:522, App. Eqs. (A.12)-(A.14):
      chi(y) = -zeta/(sqrt(alpha G)+zeta), alpha y^2 = G (alpha>alpha_c branch), H/N = G/(1+chi)^2.
      For S=2: G = (1+Q)/2, z = 1/(sqrt2 y), zeta = -erf(z/sqrt2)/(sqrt2 z)  =>  chi = erf/(alpha - erf).
  Approximate sigma^2/N = H/N + (1-Q)/2 (CMZ 2000 Eq. 5 with independent s_i) is NON-gating.
Units: a_i = +-1, A = sum a_i; random agents give sigma^2/N = 1, H/N = 1/2.
Check: the chi -> infinity condition reproduces the published alpha_c = 0.33740.
"""
import numpy as np
from scipy.special import erf
from scipy.optimize import brentq

def Q(z): return 1 - np.sqrt(2/np.pi)*np.exp(-z*z/2)/z - (1 - 1/z**2)*erf(z/np.sqrt(2))
def alpha(z): return z*z*(1+Q(z))
zc = brentq(lambda z: alpha(z) - erf(z/np.sqrt(2)), 0.05, 5)

def rs(a):
    z = brentq(lambda z: alpha(z)-a, zc*(1+1e-9), 50); q = Q(z); e = erf(z/np.sqrt(2)); c = e/(a-e)
    H = (1+q)/(2*(1+c)**2)
    return dict(z=z, Q=q, chi=c, phi=1-e, H=H, sigma2_approx=H+(1-q)/2)

if __name__ == "__main__":
    print("alpha_c = %.6f  (published 0.33740...)  Q_c=%.4f  phi(alpha_c+)=%.4f" % (alpha(zc), Q(zc), 1-erf(zc/np.sqrt(2))))
    for P, N in [(512, 511), (512, 255), (256, 255), (256, 127)]:
        r = rs(P/N)
        print("P=%d N=%d alpha=%.5f phi=%.4f H/N=%.4f Q=%.4f sigma2/N_approx=%.4f" % (P, N, P/N, r['phi'], r['H'], r['Q'], r['sigma2_approx']))
