# Locates the saddle-node (fold) of the stable chimera branch of AS2004 Eqs. (12)-(13) for A in {0.995, 0.95, 0.90}. Not a simulation.

import numpy as np, sys, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import fsolve
x=np.linspace(-np.pi,np.pi,200000,endpoint=False)+np.pi/200000; cx=np.cos(x)
def resid(v,A,beta):
    c,ar,ai,D=v; a=ar+1j*ai
    R2=c*c+2*c*ar*cx+(ar*ar+ai*ai)*cx*cx
    sq=np.sqrt((D*D-R2).astype(complex)); h=(D-sq)/(c+np.conj(a)*cx); e=np.exp(1j*beta)
    r1=e*np.mean(h)-c; r2=A*e*np.mean(h*cx)-a
    return [r1.real,r1.imag,r2.real,r2.imag]
def track(A, v0, b0, b1, step):
    v=np.array(v0); last=None
    for b in np.arange(b0,b1,step):
        s,info,ier,msg=fsolve(resid,v,args=(A,b),full_output=True,xtol=1e-13)
        if ier!=1 or max(abs(np.array(resid(s,A,b))))>1e-9: return last, s
        last=(b,s); v=s
    return last, None
for A in [0.995, 0.95, 0.90]:
    # get stable chimera at beta=0.18 (A-continuation from the A=0.995 solution)
    v=np.array([0.62867,-0.16984,0.11559,0.73201])
    for a in np.linspace(0.995,A,20):
        v=fsolve(resid,v,args=(a,0.18),xtol=1e-13)
    last,_=track(A,v,0.18,0.30,0.002)
    b,s=last
    last2,_=track(A,s,b,b+0.004,0.00005)
    print(f"A={A}: stable branch last converged beta = {last2[0]:.5f} (fold just above);  c={last2[1][0]:.4f} Delta={last2[1][3]:.4f}")
