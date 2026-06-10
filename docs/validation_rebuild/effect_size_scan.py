import json, os, math
D = "/home/matthewhmaxwell/emergent-pattern-catalog/analysis/outputs"
inf = nopanel = 0
verdicts = {}
rows = []
for n in range(1, 33):
    f = os.path.join(D, "p%d_phase2a_panel.json" % n)
    if not os.path.exists(f):
        rows.append(("p%d" % n, "NO PANEL", "-", "-")); nopanel += 1; continue
    s = json.load(open(f)).get("summary", {})
    d = s.get("cohens_d_positive_vs_panel"); v = s.get("verdict", "?"); tnr = s.get("overall_tnr")
    verdicts[v] = verdicts.get(v, 0) + 1
    if isinstance(d, float) and math.isinf(d):
        inf += 1; disp = "Infinity"
    elif isinstance(d, (int, float)):
        disp = round(d, 3)
    else:
        disp = str(d)
    tnrd = round(tnr, 3) if isinstance(tnr, (int, float)) else tnr
    rows.append(("p%d" % n, v, disp, tnrd))
print("%-5s %-20s %-12s %s" % ("pat", "verdict", "cohens_d", "overall_tnr"))
for r in rows:
    print("%-5s %-20s %-12s %s" % r)
print("")
print("d=Infinity count: %d / 32" % inf)
print("verdict distribution: %s" % verdicts)
print("no-panel: %d" % nopanel)
