# Solves D'Orsogna et al. 2006 PRL Eq. (4) (as printed, sum to N/2) and the exact equidistant-ring force balance (sum to N-1).
# Region II ring parameters of Fig. 2(c,d): m=1, alpha=1, beta=0.5, Ca=la=1, Cr=lr=0.5. Not a simulation.

# D'Orsogna et al. 2006 PRL Eq. (4): m*alpha/(2 R beta) = sum_{n=1}^{N/2} sin(pi n/N) [ (Ca/la) e^{-2R sin/la} - (Cr/lr) e^{-2R sin/lr} ]
import numpy as np
from scipy.optimize import brentq
m,al,be,Ca,la,Cr,lr = 1.0,1.0,0.5,1.0,1.0,0.5,0.5
def rhs(R,N,upper):
    n=np.arange(1,upper+1); s=np.sin(np.pi*n/N)
    return np.sum(s*((Ca/la)*np.exp(-2*R*s/la)-(Cr/lr)*np.exp(-2*R*s/lr)))
def solve(N, variant):
    if variant=='eq4':   f=lambda R: m*al/(2*R*be) - rhs(R,N,N//2)
    else:                f=lambda R: m*al/(R*be) - rhs(R,N,N-1)   # exact force balance, all N-1 partners
    Rs=np.linspace(1e-3,3,30000); v=[f(r) for r in Rs]
    roots=[brentq(f,Rs[i],Rs[i+1]) for i in range(len(Rs)-1) if np.sign(v[i])!=np.sign(v[i+1])]
    return roots
print(" N   R_eq4(roots)            R_exact(roots)   rel.diff  sqrt(2/N)")
for N in [50,60,80,100,150,200,250,300,350,400,450,460,500]:
    a=solve(N,'eq4'); b=solve(N,'exact')
    print(f"{N:4d} {str([round(x,5) for x in a]):24s} {str([round(x,5) for x in b]):18s} {((a[0]-b[0])/b[0] if a and b else float('nan')):+.4f}  {np.sqrt(2/N):.4f}")
