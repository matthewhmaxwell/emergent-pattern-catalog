"""Census digest (pilot draft v0): summarize a run directory into DIGEST.md (+ digest.json for the dashboard).

  python -m census.digest census/runs/<tag>

Tiers (DESIGN §8): 4 vetted findings, 3 novelty candidates (both filled only after triage/vetting exist), 2 catalog
expansions, 1 unclassified clusters (pilot: unclassified programs, shortest first), then progress and the
emergence-complexity table so far (shortest program per matched pattern).
"""
import collections, glob, json, os, sys, time


def load(run):
    rows = []
    for f in sorted(glob.glob(os.path.join(run, "results_*.jsonl"))):
        for line in open(f):
            try: rows.append(json.loads(line))
            except Exception: pass
    return rows


def build(run):
    rows = load(run); man = json.load(open(os.path.join(run, "manifest.json")))
    st = collections.Counter(r["status"] for r in rows)
    shortest = {}                                            # pattern -> (len, prog, view)
    for r in rows:
        for v, e in (r.get("views") or {}).items():
            if e["verdict"] == "MATCH":
                cur = shortest.get(e["pattern"])
                if cur is None or r["len"] < cur[0]: shortest[e["pattern"]] = (r["len"], r["prog"], v)
    unc = sorted([r for r in rows if r["status"] == "UNCLASSIFIED"], key=lambda r: (r["len"], r["idx"]))
    errs = [r for r in rows if r["status"] == "ERROR"]
    t_sim = [r["sim_s"] for r in rows if "sim_s" in r]; t_fil = [r["filter_s"] for r in rows if "filter_s" in r]
    by_layer = collections.defaultdict(collections.Counter)
    for r in rows: by_layer[r["layers"]][r["status"]] += 1
    d = {"run": run, "generated": time.strftime("%Y-%m-%d %H:%M"), "done": len(rows), "total": man["n_programs"],
         "status": dict(st), "by_layer": {k: dict(v) for k, v in by_layer.items()},
         "mean_sim_s": sum(t_sim) / max(len(t_sim), 1), "mean_filter_s": sum(t_fil) / max(len(t_fil), 1),
         "table": {k: {"len": v[0], "prog": v[1], "view": v[2]} for k, v in sorted(shortest.items(), key=lambda x: x[1][0])},
         "unclassified": [{"idx": r["idx"], "len": r["len"], "prog": r["prog"],
                           "views": {v: (e["verdict"], [s.get("em_kind") for s in e["seeds"]])
                                     for v, e in r["views"].items() if e["verdict"] == "UNCLASSIFIED"}} for r in unc],
         "errors": [{"idx": r["idx"], "prog": r["prog"], "error": r.get("error")} for r in errs[:20]]}
    L = [f"# Emergence Census digest — {d['generated']}", "",
         f"Run `{run}` {'(NULL: shuffled)' if man.get('null') else ''}: **{d['done']:,} / {d['total']:,}** programs "
         f"({100 * d['done'] / max(d['total'], 1):.1f}%). Mean sim {d['mean_sim_s']:.1f}s, filter {d['mean_filter_s']:.1f}s.", "",
         "## Tier 4 — vetted findings", "_none (vetting not yet run)_", "",
         "## Tier 3 — novelty candidates", "_none (triage not yet run)_", "",
         "## Tier 2 — catalog expansions", "_none yet_", "",
         f"## Tier 1 — unclassified ({len(unc)})"]
    for u in d["unclassified"][:25]:
        L.append(f"- [{u['len']} bits] `{u['prog']}` — {u['views']}")
    L += ["", "## Status", "", "| status | n |", "|---|---|"] + [f"| {k} | {v} |" for k, v in st.most_common()]
    L += ["", "## Emergence-complexity table so far (shortest program per matched pattern)", "",
          "| pattern | bits | program | view |", "|---|---|---|---|"]
    L += [f"| {k} | {v['len']} | `{v['prog']}` | {v['view']} |" for k, v in d["table"].items()]
    if errs: L += ["", f"## Errors ({len(errs)})"] + [f"- {e['idx']}: {e['error']}" for e in d["errors"]]
    open(os.path.join(run, "DIGEST.md"), "w").write("\n".join(L) + "\n")
    json.dump(d, open(os.path.join(run, "digest.json"), "w"), indent=1, default=str)
    return d


if __name__ == "__main__":
    d = build(sys.argv[1]); print(open(os.path.join(sys.argv[1], "DIGEST.md")).read())
