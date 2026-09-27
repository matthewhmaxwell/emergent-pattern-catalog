"""Census triage (pilot draft v0): cluster emergent programs into behaviour classes, then label each class.

  python -m census.triage census/runs/<tag> [--threshold T] [--reps 3] [--workers 4] [--no-battery]

1. Every EMERGENT view of every program -> feature vector (mean fingerprint over its emergent seeds + em_score).
2. Views are grouped by family (grid, agents, network, lattice-phase, network-phase, field, avalanche); within a
   family, features are robust-z-scored (median / MAD) and clustered by Ward agglomeration at distance threshold T.
3. Each cluster's REPS shortest programs are re-simulated exactly (stable seeds) and the full 36-detector battery
   is run on that view with its VALIDATED settings on full-resolution frames (3 seeds). A representative matches
   pattern P if P matches at >= confirmation in >= 2 seeds; the cluster is labelled P if >= 2 representatives agree.
4. Output clusters.json + TRIAGE.md: the emergence-complexity table (shortest program per label) and the unlabelled
   clusters, which go to the literature check (DESIGN §6 step 4).
"""
import argparse, collections, json, os, sys
from multiprocessing import Pool
import numpy as np

FAMILY = {"C": "grid", "C.transient": "grid", "A": "agents", "N": "network", "C.phase": "lattice-phase",
          "N.phase": "network-phase", "F0": "field", "F1": "field", "C.aval": "avalanche"}


def load(run):
    import glob
    rows = []
    for f in sorted(glob.glob(os.path.join(run, "results_*.jsonl"))):
        for line in open(f):
            try: rows.append(json.loads(line))
            except Exception: pass
    return rows


def items(rows):
    """[(family, view, row, featdict)] for every EMERGENT view."""
    out = []
    for r in rows:
        for v, e in (r.get("views") or {}).items():
            if e["verdict"] not in ("EMERGENT", "UNCLASSIFIED", "MATCH"): continue
            seeds = [s for s in e["seeds"] if s.get("emergent")]
            fps = [s.get("fp") or {} for s in seeds]
            keys = sorted(set().union(*fps)) if fps else []
            feat = {k: float(np.mean([f.get(k, 0.0) for f in fps])) for k in keys}
            ems = [s.get("em_score") for s in seeds if s.get("em_score") is not None]
            feat["em_score"] = float(np.mean(ems)) if ems else 0.0
            out.append((FAMILY.get(v, v), v, r, feat))
    return out


def cluster_family(its, threshold):
    from scipy.cluster.hierarchy import linkage, fcluster
    keys = sorted(set().union(*[f for *_, f in its]))
    X = np.array([[f.get(k, 0.0) for k in keys] for *_, f in its], float)
    med = np.median(X, 0); mad = np.median(np.abs(X - med), 0) * 1.4826; mad[mad < 1e-9] = 1.0
    Z = np.clip((X - med) / mad, -8, 8)
    if len(its) == 1: return np.array([1]), keys, Z
    lab = fcluster(linkage(Z, "ward"), t=threshold, criterion="distance")
    return lab, keys, Z


def _battery_rep(args):
    bits, view = args
    from census import grammar_g1 as G, sim, filter as FL
    from census.runner import seed_of
    from epc.phase2a.panel import _detected, _verdict
    p = G.decode(bits); out = sim.run(p, seed0=seed_of(bits)); B, FULL = FL._load(); per_seed = []
    for s in range(out["meta"].get("S", 3) if "S" in out["meta"] else 3):
        V = {n: (h, md) for n, h, md in FL.views(out, p, s)}
        if view not in V: continue
        h, md = V[view]; h = FL._thin(h, FL.CONFIRM_FRAMES); hits = []
        for pid, fn in FULL.items():
            try:
                res = fn(h, md)
                if _detected(res) and FL._TIER.get(_verdict(res), 0) >= 2: hits.append((pid, _verdict(res)))
            except Exception:
                pass
        hits.sort(key=lambda x: FL._TIER.get(x[1], 0), reverse=True)
        per_seed.append(hits[0][0] if hits else None)
    c = collections.Counter(x for x in per_seed if x)
    top = c.most_common(1)[0] if c else (None, 0)
    return {"bits": bits, "view": view, "per_seed": per_seed, "match": top[0] if top[1] >= 2 else None}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--threshold", type=float, default=12.0)
    ap.add_argument("--reps", type=int, default=3); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--no-battery", action="store_true"); a = ap.parse_args()
    rows = load(a.run); its = items(rows); fams = collections.defaultdict(list)
    for it in its: fams[it[0]].append(it)
    clusters = []
    for fam, fit in sorted(fams.items()):
        lab, keys, Z = cluster_family(fit, a.threshold)
        for c in sorted(set(lab)):
            mem = [fit[i] for i in np.flatnonzero(lab == c)]
            mem.sort(key=lambda x: (x[2]["len"], x[2]["idx"]))
            idx = np.flatnonzero(lab == c)
            clusters.append({"id": f"{fam}-{c}", "family": fam, "size": len(mem),
                             "shortest_len": mem[0][2]["len"], "shortest": mem[0][2]["prog"],
                             "reps": [(m[2]["bits"], m[1], m[2]["prog"]) for m in mem[:a.reps]],
                             "members": [(m[2]["idx"], m[1]) for m in mem],
                             "centroid": {k: float(v) for k, v in zip(keys, Z[idx].mean(0))},
                             "em_kinds": dict(collections.Counter(s.get("em_kind") for m in mem for s in m[2]["views"][m[1]]["seeds"] if s.get("emergent")))})
    print(f"{len(its)} emergent views -> {len(clusters)} clusters "
          f"({', '.join(f'{f}:{sum(1 for c in clusters if c["family"] == f)}' for f in sorted(fams))})", flush=True)
    if not a.no_battery:
        jobs = sorted({(b, v) for c in clusters for b, v, _ in c["reps"]})
        print(f"battery on {len(jobs)} representative views ...", flush=True)
        with Pool(a.workers) as pool: res = {(r["bits"], r["view"]): r for r in pool.imap_unordered(_battery_rep, jobs)}
        for c in clusters:
            ms = [res[(b, v)]["match"] for b, v, _ in c["reps"]]; c["rep_matches"] = ms
            cnt = collections.Counter(m for m in ms if m); top = cnt.most_common(1)[0] if cnt else (None, 0)
            c["label"] = top[0] if top[1] >= min(2, len(ms)) else None
    json.dump(clusters, open(os.path.join(a.run, "clusters.json"), "w"), indent=1, default=str)
    table = {}
    for c in clusters:
        if c.get("label") and (c["label"] not in table or c["shortest_len"] < table[c["label"]]["shortest_len"]): table[c["label"]] = c
    L = [f"# Triage — {a.run}", "", f"{len(its)} emergent views -> **{len(clusters)} behaviour classes** "
         f"(Ward, threshold {a.threshold}).", "", "## Emergence-complexity table (labelled classes)", "",
         "| pattern | bits | shortest program | class | size |", "|---|---|---|---|---|"]
    L += [f"| {k} | {c['shortest_len']} | `{c['shortest']}` | {c['id']} | {c['size']} |" for k, c in sorted(table.items(), key=lambda x: x[1]["shortest_len"])]
    un = sorted([c for c in clusters if not c.get("label")], key=lambda c: (c["shortest_len"], -c["size"]))
    L += ["", f"## Unlabelled classes -> literature check ({len(un)})", "", "| class | bits | shortest program | size | rep matches | generic lens |", "|---|---|---|---|---|---|"]
    L += [f"| {c['id']} | {c['shortest_len']} | `{c['shortest']}` | {c['size']} | {c.get('rep_matches')} | {max(c['em_kinds'], key=c['em_kinds'].get) if c['em_kinds'] else ''} |" for c in un]
    open(os.path.join(a.run, "TRIAGE.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
