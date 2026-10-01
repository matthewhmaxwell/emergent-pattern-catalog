"""Census digest: one plain summary of a run directory -> DIGEST.md (+ digest.json for the dashboard).

  python -m census.digest <run_dir>

Reads the runner rows for progress, speed and integrity, and triage.json (census.triage) when present for the
findings. Tiers (DESIGN §8): 4 vetted findings, 3 novelty candidates (both filled only by the vetting step),
2 catalog expansions, 1 unnamed behaviour classes awaiting a literature check; then the behaviour map so far.
"""
import collections, json, os, sys, time
from census.triage import load


def build(run):
    rows = load(run); man = json.load(open(os.path.join(run, "manifest.json")))
    st = collections.Counter(r["status"] for r in rows); t = [r["sim_s"] + r["post_s"] for r in rows if "sim_s" in r]
    fperr = sum(r.get("fp_errors", 0) for r in rows); errs = [r for r in rows if r["status"] == "ERROR"]
    tri = json.load(open(os.path.join(run, "triage.json"))) if os.path.exists(os.path.join(run, "triage.json")) else None
    d = {"run": run, "generated": time.strftime("%Y-%m-%d %H:%M"), "done": len(rows), "total": man["n_programs"],
         "status": dict(st), "mean_seconds": sum(t) / max(len(t), 1), "fingerprint_errors": fperr,
         "errors": [{"idx": r["idx"], "prog": r["prog"], "error": r.get("error")} for r in errs[:20]],
         "longest_done": max((r["len"] for r in rows), default=0)}
    L = [f"# Emergence Census digest — {d['generated']}", "",
         f"**{d['done']:,} / {d['total']:,}** programs done ({100 * d['done'] / max(d['total'], 1):.1f}%), "
         f"up to {d['longest_done']} bits. Mean {d['mean_seconds']:.1f} s per program. "
         f"Flagged: {st.get('FLAGGED', 0):,}. Errors: {len(errs)}. Fingerprint errors: {fperr} (must be 0).", "",
         "## Tier 4 — vetted findings", "_none (vetting has not run)_", "",
         "## Tier 3 — novelty candidates", "_none (no class has passed the literature check as literature-silent)_", "",
         "## Tier 2 — catalog expansions", "_none yet_", ""]
    if tri:
        cl = tri["unnamed_classes"]
        L += [f"## Tier 1 — unnamed behaviour classes awaiting a literature check ({len(cl)})"]
        L += [f"- [{c['shortest_len']} bits, {c['size']} programs] `{c['shortest']}`" for c in cl[:15]]
        L += ["", f"## Review — known behaviour, unusual look ({tri['review_views']})", "",
              "## Behaviour map so far", "", "| behaviour | catalog | bits | simplest program | programs | routes |", "|---|---|---|---|---|---|"]
        L += [f"| {b} | {x['catalog'] or '—'} | {x['shortest_len']} | `{x['shortest']}` | {x['programs']} | {len(x['routes'])} |"
              for b, x in sorted(tri["behaviours"].items(), key=lambda kv: kv[1]["shortest_len"])]
        d["triage_generated"] = tri["generated"]
    else:
        L += ["## Findings", "_triage has not been run on this directory yet (python -m census.triage)_"]
    if errs: L += ["", f"## Errors ({len(errs)})"] + [f"- {e['idx']}: {e['error']}" for e in d["errors"]]
    open(os.path.join(run, "DIGEST.md"), "w").write("\n".join(L) + "\n")
    json.dump(d, open(os.path.join(run, "digest.json"), "w"), indent=1, default=str)
    return "\n".join(L)


if __name__ == "__main__":
    print(build(sys.argv[1]))
