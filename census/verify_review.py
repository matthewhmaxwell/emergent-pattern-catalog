"""Stage 3 (b) review aggregation for VERIFICATION_PLAN.md (deviation V4: the review was done by Claude research agents
at the owner's direction). Reads census/runs/verify/stage3_review/outputs/*.json, checks that every Stage 3 phenomenon
was reviewed exactly once, and writes stage3_review/final_labels.json and stage3_review/REVIEW.md.

  python -m census.verify_review
"""
import collections, glob, json, os

OUT = "census/runs/verify"   # same as census.verify.OUT; standard library only so it runs anywhere
R = os.path.join(OUT, "stage3_review")
LABELS = ("DOCUMENTED", "CLOSELY RELATED", "UNDOCUMENTED CANDIDATE")
MODEL = {"KNOWN": "DOCUMENTED", "KNOWN-VARIANT": "CLOSELY RELATED", "LITERATURE-SILENT": "UNDOCUMENTED CANDIDATE"}


def main():
    rows = []
    for f in sorted(glob.glob(f"{R}/outputs/*.json")):
        for r in json.load(open(f)):
            r["group"] = os.path.basename(f)[:-5]; rows.append(r)
    s3 = {}
    for f in os.listdir(os.path.join(OUT, "stage3_blind", "verdicts")):
        x = json.load(open(os.path.join(OUT, "stage3_blind", "verdicts", f))); s3[x["pid"]] = x
    res = {r["pid"]: r for r in json.load(open(f"{OUT}/verify.json"))}
    c = collections.Counter(r["pid"] for r in rows)
    problems = [f"reviewed {n} times: {p}" for p, n in c.items() if n != 1] + \
               [f"not reviewed: {p}" for p in s3 if p not in c] + \
               [f"not a Stage 3 phenomenon: {p}" for p in c if p not in s3] + \
               [f"{r['pid']}: label {r.get('final_label')!r}" for r in rows if r.get("final_label") not in LABELS]
    final = {}
    for r in rows:
        m = MODEL.get(s3.get(r["pid"], {}).get("verdict"))
        final[r["pid"]] = {"final_label": r.get("final_label"), "model_label": m, "changed": r.get("final_label") != m,
                           "class": res[r["pid"]]["class"], "group": r["group"], "confidence": r.get("confidence"),
                           "basis": r.get("basis"), "difference": r.get("difference"),
                           "citation_problems": r.get("citation_problems") or ""}
    json.dump(final, open(f"{R}/final_labels.json", "w"), indent=1)
    cites = [x for r in rows for x in (r.get("citations_checked") or [])]
    cc = lambda k, v: sum(1 for x in cites if str(x.get(k)).lower() == v)
    L = ["# Stage 3 (b) review (deviation V4: carried out by Claude research agents at the owner's direction)", "",
         f"{len(final)} phenomena reviewed. Final labels: " +
         ", ".join(f"{n} {k.lower()}" for k, n in collections.Counter(v["final_label"] for v in final.values()).most_common()) + ".",
         f"Changed from the plan-conforming Stage 3 label: {sum(v['changed'] for v in final.values())}.",
         f"Citations checked: {len(cites)}; exist: yes {cc('exists', 'yes')}, no {cc('exists', 'no')}, unclear {cc('exists', 'unclear')}; "
         f"attribution correct: yes {cc('attribution_correct', 'yes')}, no {cc('attribution_correct', 'no')}; "
         f"support the claim: yes {cc('supports_claim', 'yes')}, partly {cc('supports_claim', 'partly')}, no {cc('supports_claim', 'no')}.", ""]
    if problems: L += ["## Integrity problems", ""] + [f"- {p}" for p in problems] + [""]
    L += ["## Label changes", "", "| phenomenon | class | Stage 3 | final | why |", "|---|---|---|---|---|"]
    for p, v in sorted(final.items()):
        if v["changed"]:
            L.append(f"| {p} | {v['class']} | {v['model_label']} | {v['final_label']} | {(v['difference'] or v['basis'] or '')[:300]} |")
    L += ["", "## Citation problems", ""] + [f"- {p}: {v['citation_problems']}" for p, v in sorted(final.items()) if v["citation_problems"].strip()]
    L += ["", "## All phenomena", "", "| phenomenon | class | final label | confidence | basis |", "|---|---|---|---|---|"]
    for p, v in sorted(final.items()):
        L.append(f"| {p} | {v['class']} | {v['final_label']} | {v['confidence']} | {(v['basis'] or '')[:240]} |")
    open(f"{R}/REVIEW.md", "w").write("\n".join(L) + "\n")
    print("\n".join(L[:8]))
    if problems: print("PROBLEMS:", problems)


if __name__ == "__main__":
    main()
