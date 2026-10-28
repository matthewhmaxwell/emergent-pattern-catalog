"""Post-hoc sensitivity analysis for VERIFICATION_PLAN.md (added after classification, found by an independent audit;
the frozen classes in verify.json are not changed).

The frozen mismatch flag searches each verdict's stored full text, which for strict-pass replies includes the JSON
block with the key "trivial_or_artifact"; the pattern "artifact" matches that key. This script recomputes the flag
with the key name removed, and re-classifies, with all card features, every phenomenon whose tested-feature set the
false flag narrowed. Writes census/runs/verify/SENSITIVITY.md and sensitivity.json.

  python -m census.verify_sensitivity
"""
import json, re
from census import verify as V


def main():
    inv = {it["pid"]: it for it in json.load(open(f"{V.OUT}/inventory.json"))}
    res = {r["pid"]: r for r in json.load(open(f"{V.OUT}/verify.json"))}
    strip = lambda t: re.sub(r'"?trivial_or_artifact"?', " ", t or "")
    fl = {p: bool(V.FLAG.search(strip(it["text"]))) for p, it in inv.items()}
    only_key = [p for p, it in inv.items() if it["flagged"] and not fl[p]]
    named = {p: [f for f in inv[p]["card_features"] if f in (inv[p]["text"] or "")] for p in inv}
    nm = [p for p in inv if fl[p] and named[p]]
    rows = []
    for p in only_key:
        it = dict(inv[p]); narrowed = bool(it["discrepant_features"])
        new = None
        if narrowed:
            it["discrepant_features"] = None; it["flagged"] = False
            new = V.classify_one(it, json.load(open(f"{V.OUT}/sims/{p.replace(':', '__')}.json")))
        rows.append({"pid": p, "frozen_class": res[p]["class"], "frozen_tested": len(res[p]["features"]),
                     "narrowed": narrowed, "all_features_class": new and new["class"],
                     "all_features_robust": new and new["robust"], "frozen_stage3": bool(res[p].get("stage3")),
                     "all_features_stage3": new and bool(new.get("stage3"))})
    out = {"flagged_frozen": sum(it["flagged"] for it in inv.values()), "flagged_without_key": sum(fl.values()),
           "census_flagged_without_key": sum(1 for p in fl if fl[p] and inv[p]["batch"] == "census"),
           "named_without_key": len(nm), "named_artifact_without_key": sum(1 for p in nm if res[p]["class"] == "ARTIFACT"),
           "false_flags": rows}
    json.dump(out, open(f"{V.OUT}/sensitivity.json", "w"), indent=1)
    L = ["# Sensitivity: mismatch flag without the JSON key name (post hoc; frozen classes unchanged)", "",
         f"Flagged: {out['flagged_frozen']} frozen, {out['flagged_without_key']} with the key name removed "
         f"(census batch {out['census_flagged_without_key']}). Flagged and naming card features: {out['named_without_key']}, "
         f"of which {out['named_artifact_without_key']} artifact under the frozen classes.", "",
         "| phenomenon | frozen class | features tested | narrowed by false flag | class with all card features | Stage 3 frozen / all features |",
         "|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['pid']} | {r['frozen_class']} | {r['frozen_tested']} | {r['narrowed']} | {r['all_features_class'] or ''} | "
                 f"{r['frozen_stage3']} / {'' if r['all_features_stage3'] is None else r['all_features_stage3']} |")
    open(f"{V.OUT}/SENSITIVITY.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
