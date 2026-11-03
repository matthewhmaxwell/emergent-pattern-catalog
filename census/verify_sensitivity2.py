"""Post-hoc sensitivity analyses requested in pre-submission review (2026-11-02). Frozen classes in verify.json are not
changed; results go to census/runs/verify/SENSITIVITY2.md and sensitivity2.json.

1. All card features: re-classify every phenomenon with all its card features tested (no named-feature narrowing).
2. Stricter reproduction: a sub-program reproduces a robust feature only if it is within 3 s of the full program AND
   closer to the full program than to the knock-out (|R - P| < |P - K|).
3. Trivial knock-outs: the share of representatives whose interaction knock-out keeps no rule, so that it is the
   frozen random start.

  python -m census.verify_sensitivity2
"""
import collections, json
import numpy as np
from census import verify as V, knockout as KO


def strict_reproducer(it, rec, robust):
    cen = rec["conditions"]["census"]; view = it["view"]; best = None
    for sp in rec["subprograms"]:
        if "fps" not in sp: continue
        ok = True
        for f in robust:
            x, y, k = V._col(cen["real"], f, view), V._col(sp["fps"], f, view), V._col(cen["knockout"], f, view)
            if y is None or abs(np.median(x) - np.median(y)) > 3 * V._scale(x, y): ok = False; break
            if k is not None and not abs(np.median(x) - np.median(y)) < abs(np.median(y) - np.median(k)): ok = False; break
        if ok and (best is None or len(sp["rules"]) < len(best)): best = sp["rules"]
    return best


def main():
    inv = {it["pid"]: it for it in json.load(open(f"{V.OUT}/inventory.json"))}
    res = {r["pid"]: r for r in json.load(open(f"{V.OUT}/verify.json"))}
    allf, strict, trivial = {}, {}, {}
    for pid, it in inv.items():
        rec = json.load(open(f"{V.OUT}/sims/{pid.replace(':', '__')}.json"))
        it2 = dict(it, discrepant_features=None)
        allf[pid] = V.classify_one(it2, rec)["class"]
        r = res[pid]
        if r["class"] in ("DECOMPOSES", "COMBINATION EFFECT"):
            b = strict_reproducer(it, rec, r["robust"])
            strict[pid] = "DECOMPOSES" if b is not None else "COMBINATION EFFECT"
        p = V.grammar(it["grammar"]).decode(it["bits"]); kp = KO.knockout_program(p)
        trivial[pid] = kp is None or len(getattr(kp, "rules", [])) == 0
    named = [p for p in inv if inv[p]["flagged"] and inv[p]["discrepant_features"]]
    out = {
        "all_features_classes": collections.Counter(allf.values()),
        "frozen_classes": collections.Counter(r["class"] for r in res.values()),
        "named_60_frozen_artifact": sum(res[p]["class"] == "ARTIFACT" for p in named),
        "named_60_allfeatures_artifact": sum(allf[p] == "ARTIFACT" for p in named),
        "unnamed_154_frozen_artifact": sum(res[p]["class"] == "ARTIFACT" for p in inv if p not in named),
        "class_changes_all_features": collections.Counter((res[p]["class"], allf[p]) for p in inv if res[p]["class"] != allf[p]),
        "strict_reproduction": collections.Counter((res[p]["class"], strict[p]) for p in strict),
        "trivial_knockout_all": sum(trivial.values()), "n": len(inv),
        "trivial_knockout_stage3": sum(trivial[p] for p in res if res[p].get("stage3")),
        "trivial_knockout_by_class": collections.Counter(res[p]["class"] for p in inv if trivial[p]),
    }
    out = {k: (dict((" -> ".join(kk) if isinstance(kk, tuple) else kk, vv) for kk, vv in v.items()) if isinstance(v, collections.Counter) else v)
           for k, v in out.items()}
    json.dump(out, open(f"{V.OUT}/sensitivity2.json", "w"), indent=1)
    L = ["# Post-hoc sensitivity analyses (pre-submission review, 2026-11-02; frozen classes unchanged)", ""] + \
        [f"- {k}: {v}" for k, v in out.items()]
    open(f"{V.OUT}/SENSITIVITY2.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
