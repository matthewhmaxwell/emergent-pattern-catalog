# Provenance: digitizes Abrams & Strogatz 2004 Fig. 2 (existence boundary) from arXiv:nlin/0407045 page 3 vector paths.
# Usage: python3 p10_digitize_as2004_fig2.py <path to nlin_0407045.pdf>. Calibration check: dashed perturbation line fits A = 4.535*beta (paper: 4.535).

import pdfplumber, numpy as np, json
pdf=pdfplumber.open(__import__('sys').argv[1]); p=pdf.pages[2]
bx0,bx1=96.45,268.37; ay0,ay1=163.74,55.13
pts=[]
dashed=[]
for c in p.curves:
    if not (90<c['x0']<270 and 50<c['top']<170): continue
    w=c['x1']-c['x0']; h=c['bottom']-c['top']
    if len(c['pts'])==5 and w<1.2 and h<1.6:   # tiny dash pieces = perturbation-theory dashed line
        dashed+= [((x-bx0)/(bx1-bx0)*0.3,(ay0-y)/(ay0-ay1)) for x,y in c['pts']]; continue
    pts += [((x-bx0)/(bx1-bx0)*0.3,(ay0-y)/(ay0-ay1)) for x,y in c['pts']]
pts=np.array(sorted(set((round(b,5),round(a,5)) for b,a in pts), key=lambda t:t[1]))
print("solid boundary npts",len(pts), "A range",pts[:,1].min(),pts[:,1].max(), "beta range",pts[:,0].min(),pts[:,0].max())
# beta at given A by linear interpolation (boundary is monotone in A)
order=np.argsort(pts[:,1]); A=pts[order,1]; B=pts[order,0]
for a in [0.2,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.97,0.99,0.995,1.0]:
    print(f"A={a:5.3f}  beta_boundary={np.interp(a,A,B):.4f}")
d=np.array(dashed)
if len(d):
    # slope of dashed line A/beta
    sl=np.polyfit(d[:,0],d[:,1],1); print("dashed fit A = %.3f*beta + %.4f"%(sl[0],sl[1]), "n",len(d))
json.dump({'solid':pts.tolist(),'dashed':d.tolist()},open('as2004_fig2_boundary.json','w'))
