# Solves Abrams & Strogatz 2004 PRL Eqs. (12)-(13) (cosine kernel) and continues stable/unstable chimera branches in beta at fixed A.
# Usage: python3 p10_solve_as2004_selfconsistency_branches.py 0.995 . Not a simulation.

import numpy as np, sys, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import fsolve
def make(M):
    x=np.linspace(-np.pi,np.pi,M,endpoint=False)+np.pi/M; return np.cos(x)
cx=make(200000)
def resid(v,A,beta):
    c,ar,ai,D=v; a=ar+1j*ai
    R2=c*c+2*c*ar*cx+(ar*ar+ai*ai)*cx*cx
    sq=np.sqrt((D*D-R2).astype(complex)); h=(D-sq)/(c+np.conj(a)*cx); e=np.exp(1j*beta)
    r1=e*np.mean(h)-c; r2=A*e*np.mean(h*cx)-a
    return [r1.real,r1.imag,r2.real,r2.imag]
def props(v):
    c,ar,ai,D=v; R=np.abs(c+(ar+1j*ai)*cx); return c,D,R.min(),R.max(),np.mean(R<D)
A=float(sys.argv[1])
for label,start in [('stable',[0.62867,-0.16984,0.11559,0.73201]),('unstable',[0.89207,-0.0463,0.19712,0.88324])]:
    print(label)
    # continue downward and upward in beta from 0.18
    for direction in [np.arange(0.18,0.0,-0.01), np.r_[np.arange(0.18,0.255,0.01),np.arange(0.255,0.275,0.001)]]:
        v=np.array(start)
        for b in direction:
            s,info,ier,msg=fsolve(resid,v,args=(A,b),full_output=True,xtol=1e-12)
            ok= ier==1 and max(abs(np.array(resid(s,A,b))))<1e-8
            if not ok: print(f"   beta={b:.3f}  no convergence -> branch lost"); break
            c,D,Rmin,Rmax,fd=props(s)
            if not (0<fd<1): print(f"   beta={b:.3f} left chimera family fd={fd:.3f}"); break
            print(f"   beta={b:.3f} c(=r_global)={c:.4f} Delta={D:.4f} Rmin={Rmin:.4f} Rmax={Rmax:.4f} fdrift={fd:.4f}")
            v=s
