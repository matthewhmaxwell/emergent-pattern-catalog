# Provenance: digitizes Couzin et al. 2005 Nature 433:513 Fig. 1a from the publisher PDF (page 2) vector polylines.
# Usage: python3 p19_digitize_couzin2005_fig1a.py <path to nature03236.pdf>  (main-panel axes: p=0 at x=94.11, p=1 at 263.26; acc=0 at y=517.19, acc=1 at 396.92 pt)

import sys, pdfplumber, json
pdf=pdfplumber.open(sys.argv[1]); p=pdf.pages[1]
cols={(0.14,0.12,0.12):'N=10',(0.0,0.68,0.94):'N=30',(0.0,0.67,0.31):'N=50',(0.94,0.64,0.17):'N=100',(0.92,0.18,0.18):'N=200'}
X0,X1=94.11,263.26; Y0,Y1=517.19,396.92
H=p.height
out={}
for o in p.curves + p.lines:
    sc=o.get('stroking_color')
    if not isinstance(sc,(list,tuple)): continue
    col=tuple(round(x,2) for x in sc)
    if col not in cols: continue
    if not (90<o['x0']<270 and 390<o['top']<522): continue
    pts=o.get('pts',[])
    w=o['x1']-o['x0']
    if w<5: continue  # skip markers / short
    # pdfplumber pts are (x, top) in newer versions
    vs=[((x-X0)/(X1-X0),(Y0-y)/(Y0-Y1)) for x,y in pts]
    out.setdefault(cols[col],[]).append(vs)
for k in ['N=10','N=30','N=50','N=100','N=200']:
    for i,vs in enumerate(out.get(k,[])):
        # skip inset polylines: inset x-range maps to p in ~[0.37,0.96] and y in ~[0.17,0.66] with tiny span
        print(f"\n{k} polyline {i} nv={len(vs)} prange=({min(v[0] for v in vs):.3f},{max(v[0] for v in vs):.3f})")
        print("  ".join(f"{a:.3f}:{b:.3f}" for a,b in vs))
json.dump(out,open('couzin_fig1a_polylines.json','w'))
