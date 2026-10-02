import json, collections, sys
from census.namer import FAMILY
rows = []
for t in sys.argv[1:]: rows += json.load(open(f"census/validation/{t}/placebo.json"))
agg = collections.defaultdict(collections.Counter)
for r in rows:
    hdr = r["prog"].split(" | ")[0]; null = hdr + " | " + (r["null"].split(" | ", 1)[1] if " | " in r["null"] else "-")
    for v, e in r["views"].items():
        a = agg[(FAMILY.get(v, v), null, v)]
        for c in ("runs", "generic", "complex", "gain"): a[c] += e[c]
        a["views"] += 1; a["flag"] += e["flagged"]; a["em2"] += e["em_seeds"] >= 2
        if e["flagged"]: print("FLAGGED:", r["prog"], "| view", v, "|", {k: e[k] for k in ("em_seeds", "generic", "complex", "gain", "runs", "kinds")})
for fam in ("agents", "network", "field", "grid", "lattice-phase", "network-phase"):
    ks = [(k, a) for k, a in agg.items() if k[0] == fam]
    print(f"\n=== {fam}: {len(ks)} distinct null programs, {sum(a['views'] for _, a in ks)} views")
    for ch in ("generic", "complex", "gain"):
        mid = [(k, a) for k, a in ks if 0.1 <= a[ch] / a["runs"] <= 0.9]
        if not mid: continue
        print(f"  {ch}: in-between fire rate (10-90% of seed-runs): {len(mid)} null programs, {sum(a['views'] for _, a in mid)} views, flags {sum(a['flag'] for _, a in mid)}")
        for k, a in sorted(mid, key=lambda x: -x[1]["views"])[:7]:
            print(f"     {a[ch]:4d}/{a['runs']:4d} runs  views={a['views']:3d} em>=2:{a['em2']:3d} flags={a['flag']}  {k[2]:5s} {k[1][:80]}")
