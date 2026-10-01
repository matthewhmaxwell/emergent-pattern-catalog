"""Census triage (validated pipeline): name flagged views, build the behaviour map with routes, group the rest.

  python -m census.triage <run_dir> --lib census/reflib/<name>/library.json [--threshold 12]

Input: runner rows (census.runner) — per program, per view: flagged?, per-seed fingerprint + textbook measures.
1. NAME every flagged view with census.namer.CensusNamer (the gate-validated nearest-example namer + name-then-
   verify), reported at BEHAVIOUR level. Avalanche views have no fingerprint; they are flagged by the catalog's SOC
   detector itself and are reported as "power-law avalanches" (catalog P14, detector-named).
2. BEHAVIOUR MAP: for each behaviour — how many programs show it, the SHORTEST program (the emergence-complexity
   table), and its ROUTES: every distinct mechanism (layer combination + set of rule templates) that produces it,
   with a count and the shortest example. Ingredients common to all routes are listed (necessary-ingredient hint).
3. REVIEW: views that pass a known behaviour's textbook measure but do not look like its library examples
   (KNOWN-BEHAVIOUR-ATYPICAL-LOOK) — a possible new twist on a known behaviour.
4. UNNAMED: views with no name at all, grouped by fingerprint (Ward clustering per substrate family, robust z,
   distance threshold T) into behaviour classes for the literature check (census/LITCHECK_PROTOCOL.md).
Outputs triage.json + TRIAGE.md in the run directory.
"""
import argparse, collections, glob, json, os, re, time
import numpy as np
from census.namer import FAMILY, CensusNamer, robust

NULL_TEMPLATES = {"C.FLIP", "N.FLIP", "F.DECAY", "F.FEED", "C.EMIT"}


def is_null_program(prog):                                   # kept for census.validate2 / validate3
    names = set(re.findall(r"([CANF]\.[A-Z]+)\(", prog))
    return bool(names) and names <= NULL_TEMPLATES


def load(run):
    rows = {}
    for f in sorted(glob.glob(os.path.join(run, "results_*.jsonl"))):
        for line in open(f):
            try: r = json.loads(line); rows[r.get("idx", len(rows))] = r          # de-duplicate by program index
            except Exception: pass
    return list(rows.values())


def route_key(r): return (r["layers"], tuple(r.get("templates") or sorted(set(re.findall(r"([CANF]\.[A-Z]+)\(", r["prog"])))))


def mean_fp(seeds):
    fps = [s["fp"] for s in seeds if s.get("emergent") and s.get("fp")]
    keys = sorted(set().union(*fps)) if fps else []
    return {k: float(np.mean([f.get(k, 0.0) for f in fps])) for k in keys}


def cluster(items, threshold):
    """items: [(row, view, fp)] of one family -> cluster labels."""
    from scipy.cluster.hierarchy import linkage, fcluster
    if len(items) == 1: return np.array([1])
    keys = sorted(set().union(*[f for *_, f in items]))
    X = np.array([[f.get(k, 0.0) for k in keys] for *_, f in items], float); med, sc = robust(X)
    return fcluster(linkage(np.clip((X - med) / sc, -8, 8), "ward"), t=threshold, criterion="distance")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--lib", required=True)
    ap.add_argument("--threshold", type=float, default=12.0); a = ap.parse_args()
    rows = load(a.run); namer = CensusNamer(a.lib)
    named, review, unnamed = [], [], collections.defaultdict(list)
    for r in rows:
        for v, e in (r.get("views") or {}).items():
            if not e.get("flagged"): continue
            if v == "C.aval":
                nm = {"primary_class": None, "behaviour": "power-law avalanches", "catalog": "P14", "also": [],
                      "shown": ["power-law avalanches"], "status": "NAMED (catalog detector)"}
            else:
                nm = namer.name_view(v, e["seeds"])
            rec = {"idx": r["idx"], "len": r["len"], "prog": r["prog"], "view": v, "route": route_key(r), **nm}
            if nm["status"].startswith("NAMED"): named.append(rec)
            elif nm["status"] == "KNOWN-BEHAVIOUR-ATYPICAL-LOOK": review.append(rec)
            else: unnamed[FAMILY.get(v, v)].append((r, v, mean_fp(e["seeds"])))
    # behaviour map: a program counts for its primary behaviour and for every also-shown behaviour
    beh = collections.defaultdict(list)
    for rec in named:
        for b in [rec["behaviour"]] + rec["also"]: beh[b].append(rec)
    for rec in review:
        for b in rec["shown"]: beh[b].append({**rec, "atypical": True})
    table = {}
    for b, recs in beh.items():
        progs = {}                                           # one entry per program
        for x in recs:
            if x["idx"] not in progs or x["len"] < progs[x["idx"]]["len"]: progs[x["idx"]] = x
        ps = sorted(progs.values(), key=lambda x: (x["len"], x["idx"])); routes = collections.defaultdict(list)
        for x in ps: routes[x["route"]].append(x)
        tmpl_sets = [set(k[1]) for k in routes]
        table[b] = {"catalog": next((x["catalog"] for x in ps if x.get("catalog")), None), "programs": len(ps),
                    "typical_look": sum(1 for x in ps if not x.get("atypical")),
                    "shortest_len": ps[0]["len"], "shortest": ps[0]["prog"],
                    "common_templates": sorted(set.intersection(*tmpl_sets)) if tmpl_sets else [],
                    "routes": sorted([{"layers": k[0], "templates": list(k[1]), "programs": len(v),
                                       "shortest_len": v[0]["len"], "shortest": v[0]["prog"]} for k, v in routes.items()],
                                     key=lambda x: (x["shortest_len"], -x["programs"]))}
    clusters = []
    for fam, items in sorted(unnamed.items()):
        lab = cluster(items, a.threshold)
        for c in sorted(set(lab)):
            mem = sorted([items[i] for i in np.flatnonzero(lab == c)], key=lambda x: (x[0]["len"], x[0]["idx"]))
            clusters.append({"id": f"{fam}-{c}", "family": fam, "size": len(mem), "shortest_len": mem[0][0]["len"],
                             "shortest": mem[0][0]["prog"], "examples": [m[0]["prog"] for m in mem[:5]],
                             "members": [(m[0]["idx"], m[1]) for m in mem]})
    clusters.sort(key=lambda c: (c["shortest_len"], -c["size"]))
    st = collections.Counter(r["status"] for r in rows)
    fperr = sum(r.get("fp_errors", 0) for r in rows)
    out = {"run": a.run, "library": a.lib, "generated": time.strftime("%Y-%m-%d %H:%M"), "programs": len(rows),
           "status": dict(st), "fingerprint_errors": fperr, "flagged_views": len(named) + len(review) + sum(len(x) for x in unnamed.values()),
           "named_views": len(named), "review_views": len(review), "unnamed_views": sum(len(x) for x in unnamed.values()),
           "behaviours": table, "review": review[:200], "unnamed_classes": clusters}
    json.dump(out, open(os.path.join(a.run, "triage.json"), "w"), indent=1, default=str)
    L = [f"# Census triage — {out['generated']}", "",
         f"{len(rows):,} programs; {st.get('FLAGGED', 0):,} with a flagged view; errors {st.get('ERROR', 0)}; "
         f"fingerprint errors {fperr} (must be 0).",
         f"Flagged views: {out['flagged_views']:,} = {len(named):,} named + {len(review):,} known-behaviour-atypical-look "
         f"(review) + {out['unnamed_views']:,} unnamed in {len(clusters)} classes (literature check).", "",
         "## Behaviour map (shortest program per behaviour = the emergence-complexity table)", "",
         "| behaviour | catalog | bits | shortest program | programs | routes |", "|---|---|---|---|---|---|"]
    for b, t in sorted(table.items(), key=lambda x: x[1]["shortest_len"]):
        L.append(f"| {b} | {t['catalog'] or '—'} | {t['shortest_len']} | `{t['shortest']}` | {t['programs']} | {len(t['routes'])} |")
    L += ["", "## Routes (every distinct mechanism that produces each behaviour)"]
    for b, t in sorted(table.items(), key=lambda x: x[1]["shortest_len"]):
        L += ["", f"### {b} — {len(t['routes'])} routes; ingredient common to all: {', '.join(t['common_templates']) or 'none'}", "",
              "| layers | rule templates | programs | bits | shortest example |", "|---|---|---|---|---|"]
        L += [f"| {x['layers']} | {', '.join(x['templates'])} | {x['programs']} | {x['shortest_len']} | `{x['shortest']}` |" for x in t["routes"][:15]]
        if len(t["routes"]) > 15: L.append(f"| … | {len(t['routes']) - 15} more routes in triage.json | | | |")
    L += ["", f"## Review — known behaviour, atypical look ({len(review)})", ""]
    L += [f"- [{x['len']} bits] `{x['prog']}` view {x['view']} — shows: {', '.join(x['shown'])}" for x in sorted(review, key=lambda x: x["len"])[:25]]
    L += ["", f"## Unnamed behaviour classes -> literature check ({len(clusters)})", "",
          "| class | bits | shortest program | size |", "|---|---|---|---|"]
    L += [f"| {c['id']} | {c['shortest_len']} | `{c['shortest']}` | {c['size']} |" for c in clusters[:40]]
    open(os.path.join(a.run, "TRIAGE.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
