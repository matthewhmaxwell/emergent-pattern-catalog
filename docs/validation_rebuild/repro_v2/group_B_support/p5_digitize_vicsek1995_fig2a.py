# Provenance: digitizes Vicsek et al. 1995 Fig. 2a from arXiv:cond-mat/0611743 (page 20), glyph positions via pdfplumber.
# Usage: python3 p5_digitize_vicsek1995_fig2a.py <path to cond-mat_0611743.pdf>  -> writes vicsek_fig2a_digitized.json
# Axis ticks were located on a 300-dpi raster: eta=0 at x=164.40 pt, eta=5 at 398.16 pt; v_a=1 at y=96.96 pt, v_a=0 at 312.72 pt.
# Glyph bitmap centres come from the Type-3 CharProcs; each series is anchored so v_a(eta=0.1)=0.9995 (shift <=0.004).

import sys, pdfplumber, json
from collections import defaultdict
pdf = pdfplumber.open(sys.argv[1]); page = pdf.pages[19]
U = 0.24  # pt per glyph unit
# data-font glyph bitmap centers (glyph units, y measured upward from origin), from CharProcs
centers = {'2':(37, 18.6,17.6,'N=40  (square)'), '+':(38,18.6,12.1,'N=100 (+)'),
           '(cid:2)':(39,19.1,12.1,'N=400 (x)'), '4':(45,22.1,18.1,'N=4000 (triangle)'),
           '3':(39,19.6,17.1,'N=10000 (diamond)')}
X0, XPU = 164.40, (398.16-164.40)/5.0     # eta=0 tick, pt per eta unit (from rendered ticks)
Y0, YPU = 312.72, (312.72-96.96)/1.0      # v_a=0 tick, pt per v_a unit
series = defaultdict(list)
for c in page.chars:
    t=c['text']
    if t not in centers or c['top']>320 or c['top']<95: continue
    adv,cx,cy,name = centers[t]
    if abs(c['width'] - adv*U) > 0.05: continue   # skip legend glyphs (larger font) and axis digits
    xc = c['x0'] + cx*U; yc = c['top'] - cy*U
    series[name].append(((xc-X0)/XPU, (Y0-yc)/YPU))
out={}
for name in sorted(series):
    pts = sorted(set((round(e,3),round(v,4)) for e,v in series[name]))
    # anchor correction: force the eta~0.1 point to 0.9995 (ordered-phase value), report shift
    first = [p for p in pts if abs(p[0]-0.1)<0.03]
    shift = (0.9995-first[0][1]) if first else 0.0
    out[name]={'anchor_shift':round(shift,4),'pts':[(e,round(v+shift,4)) for e,v in pts]}
    print(f"\n{name}: n={len(pts)} anchor_shift={shift:+.4f}")
    print("  ".join(f"{e:.2f}:{v+shift:.3f}" for e,v in pts))
json.dump(out, open('vicsek_fig2a_digitized.json','w'), indent=1)
