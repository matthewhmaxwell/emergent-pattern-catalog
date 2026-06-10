import json, os, time
from collections import Counter
D = "/home/matthewhmaxwell/emergent-pattern-catalog/analysis/outputs"
rows = []
guard_reasons = Counter()
now = time.time()
for n in range(1, 33):
    pid = "p%d" % n
    f = os.path.join(D, pid + "_phase2a_panel.json")
    if not os.path.exists(f):
        rows.append((pid, None)); continue
    try:
        d = json.load(open(f))
    except Exception as e:
        rows.append((pid, {"err": str(e)[:40]})); continue
    cat = d.get("catalog", []) or []
    c = Counter(); reasons = []
    for r in cat:
        rs = str(r.get("rejection_stage", "MISSING"))
        if rs.startswith("GUARD"):
            reason = rs.split(":", 1)[1] if ":" in rs else "?"
            c["GUARD"] += 1; reasons.append(reason); guard_reasons[reason] += 1
        else:
            c[rs] += 1
    rows.append((pid, {"n": len(cat), "FIRED": c.get("FIRED", 0), "GUARD": c.get("GUARD", 0),
                       "METRIC": c.get("METRIC", 0), "MISSING": c.get("MISSING", 0),
                       "reasons": reasons, "stale": (now - os.path.getmtime(f)) > 36000}))

print("## Phase 1 guard-dependence table (catalog-mate negatives)")
print("")
print("| pattern | mates | FIRED(fp) | GUARD-gated | METRIC | verdict |")
print("|--|--|--|--|--|--|")
flagged = []
for pid, s in rows:
    if s is None:
        print("| %s | - | - | - | - | **NO PANEL** |" % pid); flagged.append((pid, "no panel")); continue
    if "err" in s:
        print("| %s | - | - | - | - | ERR %s |" % (pid, s["err"])); flagged.append((pid, "err")); continue
    n, fired, guard, metric, miss = s["n"], s["FIRED"], s["GUARD"], s["METRIC"], s["MISSING"]
    tag = " (STALE)" if s["stale"] else ""
    if miss > 0:
        verdict = "? no rejection_stage field"; flagged.append((pid, "missing-field"))
    elif n == 0:
        verdict = "no catalog mates"
    elif fired > 0:
        verdict = "FALSE-POSITIVE"; flagged.append((pid, "%d false-pos" % fired))
    elif guard > 0 and guard >= metric:
        verdict = "GUARD-GATED"; flagged.append((pid, "%d/%d guard-gated" % (guard, n)))
    else:
        verdict = "ok (metric does the work)"
    print("| %s | %d | %d | %d | %d | %s%s |" % (pid, n, fired, guard, metric, verdict, tag))

print("")
print("## Distinct guard reasons doing the gating (across all patterns)")
print("")
for reason, ct in guard_reasons.most_common():
    print("- `%s` x%d" % (reason, ct))
print("")
print("## Flagged patterns (%d)" % len(flagged))
print("")
for pid, why in flagged:
    print("- %s: %s" % (pid, why))
