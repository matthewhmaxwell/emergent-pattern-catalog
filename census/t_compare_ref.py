"""Pilot check: do cluster-first labels agree with the per-program battery (reference mode) on shared programs?

  python -m census.t_compare_ref <cluster_run_dir> <reference_run_dir>

For each program present in both runs (same canonical bits), compare, per view: reference per-program verdict
(MATCH P / UNCLASSIFIED / NONE) vs cluster-first outcome (label of the program's class / unlabelled class / trivial /
not emergent). Reports the agreement table and every disagreement.
"""
import collections, json, sys
from census.triage import load

run, ref = sys.argv[1], sys.argv[2]
clusters = json.load(open(f"{run}/clusters.json"))
cls_of = {}
for c in clusters:
    for idx, v in c["members"]: cls_of[(idx, v)] = c
rows = {r["bits"]: r for r in load(run)}; refrows = {r["bits"]: r for r in load(ref)}
table = collections.Counter(); dis = []
for bits, rr in refrows.items():
    r = rows.get(bits)
    if r is None or "views" not in rr or "views" not in r: continue
    for v, e in rr["views"].items():
        ref_out = f"MATCH {e['pattern']}" if e["verdict"] == "MATCH" else e["verdict"]
        ce = r["views"].get(v, {})
        if ce.get("verdict") != "EMERGENT": cf = "not-emergent"
        else:
            c = cls_of.get((r["idx"], v))
            cf = "trivial(null region)" if c is None else (f"class {c['label']}" if c.get("label") else "unlabelled class")
        table[(ref_out, cf)] += 1
        agree = (ref_out.startswith("MATCH") and cf == f"class {e['pattern']}") or \
                (ref_out == "NONE" and cf in ("not-emergent", "trivial(null region)")) or \
                (ref_out == "UNCLASSIFIED" and cf in ("unlabelled class",))
        if not agree: dis.append((r["prog"], v, ref_out, cf))
print("reference verdict  ->  cluster-first outcome : count")
for (a, b), n in sorted(table.items(), key=lambda x: -x[1]): print(f"  {a:22s} -> {b:26s} : {n}")
print(f"\n{len(dis)} disagreements (of {sum(table.values())} shared views):")
for d in dis[:40]: print("  ", d)
