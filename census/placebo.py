"""Placebo test of the FLAG rule: interaction-free runs compared with interaction-free runs must never be flagged.

  python -m census.placebo --from-run <run_dir> [<run_dir> ...] --out census/validation/<tag> [--workers 5] [--salt x]

For every program p of the given runs, take its knock-out kp (interaction rules removed — interaction-free by
construction) and run it TWICE on independent seeds: once in the role of the "real" run, once in the role of the
knock-out. The two runs are two realizations of the same interaction-free process, so any view the pipeline flags
(screen emergent in >= 2 of 3 seeds AND "interaction-driven") is a false flag by construction. This is the generic
form of the failure the hand-built negatives caught twice by luck (2af6ba9: network noise; final3: uncoupled
oscillators): a screen channel that fires at random on unordered data makes a negligible interaction look decisive.
Unlike the hand-built negatives it covers every substrate, header and view type the census will meet, in the census's
own proportions, and it is cheap enough to run in the thousands.

Reported per view family: views tested, false flags, and how often each screen channel fires per seed-run on
interaction-free data (generic score >= 0.5, model-free complexity, consensus gain) — a channel that fires often on
interaction-free data is unreliable for that family.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, collections, json, time
from multiprocessing import Pool


def one(job):
    bits, gname, salt = job
    from census import sim
    from census.runner import grammar, seed_of
    from census.knockout import knockout_program, per_seed_view_results, interaction_driven
    M = grammar(gname); p = M.decode(bits); kp = knockout_program(p); row = {"bits": bits, "prog": M.describe(p), "null": M.describe(kp) if kp.rules else "(no rules)", "views": {}}
    try:
        a = per_seed_view_results(sim.run(kp, seed0=seed_of(bits + "|A" + salt)), kp)
        b = per_seed_view_results(sim.run(kp, seed0=seed_of(bits + "|B" + salt)), kp)
        drv = interaction_driven(a, b)
        for v, seeds in a.items():
            em = sum(1 for x in seeds if x["emergent"])
            row["views"][v] = {"flagged": bool(sum(1 for s, x in enumerate(seeds) if x["emergent"] and drv[v][s]) >= 2),
                               "em_seeds": em, "driven": bool(drv[v][0]) if seeds else False,
                               "generic": sum(1 for x in seeds + b.get(v, []) if (x.get("em_score") or 0) >= 0.5),
                               "complex": sum(1 for x in seeds + b.get(v, []) if x.get("is_complex")),
                               "gain": sum(1 for x in seeds + b.get(v, []) if (x.get("consensus_gain") or 0) >= 0.2),
                               "aval": sum(1 for x in seeds + b.get(v, []) if x["emergent"]) if v == "C.aval" else 0,
                               "runs": len(seeds) + len(b.get(v, [])),
                               "kinds": [x.get("em_kind") for x in seeds + b.get(v, []) if (x.get("em_score") or 0) >= 0.5]}
    except Exception as e:
        row["error"] = f"{type(e).__name__}: {e}"
    return row


def main():
    from census.namer import FAMILY
    from census.triage import load
    ap = argparse.ArgumentParser(); ap.add_argument("--from-run", nargs="+", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=5); ap.add_argument("--salt", default=""); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); jobs = []; seen = set()
    for run in a.from_run:
        for r in load(run):
            g = r.get("grammar", "g1")
            if (g, r["bits"]) not in seen: seen.add((g, r["bits"])); jobs.append((r["bits"], g, a.salt))      # same bit string, different grammar = different program
    t = time.time()
    with Pool(a.workers) as pool: rows = list(pool.imap_unordered(one, jobs, chunksize=4))
    json.dump(rows, open(os.path.join(a.out, "placebo.json"), "w"))
    fam = collections.defaultdict(lambda: collections.Counter()); flagged = []
    for r in rows:
        for v, e in r["views"].items():
            f = fam[FAMILY.get(v, v)]; f["views"] += 1; f["flagged"] += e["flagged"]; f["runs"] += e["runs"]
            f["generic"] += e["generic"]; f["complex"] += e["complex"]; f["gain"] += e["gain"]; f["aval"] += e["aval"]
            f["em2"] += e["em_seeds"] >= 2
            if e["flagged"]: flagged.append((r["prog"], v, r["null"]))
    nv = sum(f["views"] for f in fam.values()); nf = sum(f["flagged"] for f in fam.values())
    L = [f"# Placebo test — {time.strftime('%Y-%m-%d %H:%M')}", "",
         f"{len(rows)} programs ({sum('error' in r for r in rows)} errors), {nv} interaction-free views compared with an independent "
         f"interaction-free run of the same program. Any flag is false by construction.", "",
         f"**False flags: {nf} / {nv}**" + (f" (rate <= {3 / nv:.2%} with 95% confidence, rule of three)" if nf == 0 and nv else ""), "",
         "| view family | views | FALSE FLAGS | screen emergent in >= 2 of 3 seeds | per seed-run: generic >= 0.5 | model-free complexity | consensus gain | avalanche |",
         "|---|---|---|---|---|---|---|---|"]
    for k, f in sorted(fam.items()):
        pc = lambda x: f"{f[x]} ({f[x] / max(f['runs'], 1):.1%})"
        L.append(f"| {k} | {f['views']} | {f['flagged']} | {f['em2']} | {pc('generic')} | {pc('complex')} | {pc('gain')} | {pc('aval')} |")
    if flagged:
        L += ["", "## False flags", ""] + [f"- `{p}` view {v} — null run: `{n}`" for p, v, n in flagged[:60]]
    L += ["", f"({time.time() - t:.0f} s; runs: {', '.join(a.from_run)}; salt '{a.salt}')"]
    open(os.path.join(a.out, "PLACEBO.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
