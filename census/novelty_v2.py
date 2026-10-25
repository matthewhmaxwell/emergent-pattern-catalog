"""Quality-diversity arm v2 — DEVIATION 2 (owner-approved 2026-10-24, post-registration, exploratory).

The registered arm (census/novelty.py, unchanged) sized its CVT archive from a 200-program initialization:
k = floor(flagged/2) cells per substrate combination -> 25 cells in total, all filled almost at once, so its stop
rule fired mechanically at 20,200 evaluations; cells+agents and cells+fields combinations got no cells at all.

Changes here, and ONLY these:
  1. targeted initialization: random coupled programs are generated per combination until it has >= MIN_FLAGGED
     flagged programs or MAX_TRIES evaluations;
  2. a fixed CELLS cells per combination (k = min(CELLS, flagged//2) only if a combination stays below 2*CELLS);
  3. every combination with >= 4 flagged programs gets a CVT;
  4. resumable: cvt.json + replay of results_0.jsonl.
Imported unchanged from census.novelty: the search space (coupled G1 programs), evaluate() (the validated census
pipeline), quality (max mean evidence over flagged views), descriptor(), mutate(), the elite-replacement rule.
Budget and stop rule as registered: 150,000 evaluations, or 20,000 consecutive evaluations with no new niche.

  python -m census.novelty_v2 search    --out census/runs/census_v1_novelty_v2 [--workers 5]
  python -m census.novelty_v2 downstream --out census/runs/census_v1_novelty_v2     # name, group, literature-check
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, collections, json, random, re, time
from multiprocessing import Pool
import numpy as np
from census import grammar_g1 as G
from census import novelty as NV

COMBOS = ["CA", "CF", "AF", "CAF", "C", "N"]      # registered space: MULTI + single-layer adaptive programs
CELLS, MIN_FLAGGED, MAX_TRIES = 64, 128, 2500
BUDGET, PATIENCE, BATCH = 150_000, 20_000, 20
INIT_BATCH = 40


def log(*a):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), *a, flush=True)


def rand_program_combo(rng, combo):
    """census.novelty.random_program, with the substrate combination fixed instead of drawn."""
    for _ in range(200):
        k = {L: rng.randint(1, 4) for L in "CAN" if L in combo}
        agent = ((rng.choice([v for v, _ in G.hier_codes(G.GRIDS["speed"])]),
                  rng.choice([v for v, _ in G.hier_codes(G.GRIDS["angle"])])) if "A" in combo else ())
        nf = rng.choice([1, 2]) if "F" in combo else 0
        D = tuple(rng.choice([v for v, _ in G.hier_codes(G.GRIDS["diff"])]) for _ in range(nf))
        q = NV._repair(G.make(combo, k, nf, D, agent), rng)
        if q is not None:
            return q
    return None


def read_rows(path):
    rows = []
    if os.path.exists(path):
        for line in open(path):
            try: rows.append(json.loads(line))
            except Exception: pass
    return rows


def build_cvt(rows, seed):
    from scipy.cluster.vq import kmeans2
    cvt, info = {}, {}
    for combo in COMBOS:
        rs = [r for r in rows if r["layers"] == combo and r["status"] == "FLAGGED" and r.get("fp_mean")]
        info[combo] = {"tried": sum(1 for r in rows if r["layers"] == combo), "flagged": len(rs)}
        keys = sorted({(v, f) for r in rs for v, d in r["fp_mean"].items() for f in d})
        if not keys or len(rs) < 4:
            info[combo]["cells"] = 0; continue
        X = np.array([NV.descriptor(r, keys) for r in rs]); mu, sd = X.mean(0), X.std(0) + 1e-9
        k = CELLS if len(rs) >= 2 * CELLS else max(2, len(rs) // 2)
        C, _ = kmeans2((X - mu) / sd, k, minit="++", seed=seed)
        cvt[combo] = {"keys": keys, "mu": mu, "sd": sd, "C": C}; info[combo]["cells"] = int(k)
    return cvt, info


def save_cvt(cvt, info, path):
    json.dump({"info": info, "cvt": {c: {"keys": v["keys"], "mu": v["mu"].tolist(), "sd": v["sd"].tolist(),
               "C": v["C"].tolist()} for c, v in cvt.items()}}, open(path, "w"))


def load_cvt(path):
    d = json.load(open(path))
    return ({c: {"keys": [tuple(k) for k in v["keys"]], "mu": np.array(v["mu"]), "sd": np.array(v["sd"]),
                 "C": np.array(v["C"])} for c, v in d["cvt"].items()}, d["info"])


def cell_of(r, cvt):
    c = cvt.get(r["layers"])
    if c is None or r["status"] != "FLAGGED" or r.get("quality", 0) <= 0: return None
    z = (NV.descriptor(r, c["keys"]) - c["mu"]) / c["sd"]
    return (r["layers"], int(((c["C"] - z) ** 2).sum(1).argmin()))


def step_search(out, workers, seed):
    os.makedirs(out, exist_ok=True)
    res_path, cvt_path = os.path.join(out, "results_0.jsonl"), os.path.join(out, "cvt.json")
    rows = read_rows(res_path)
    rng = random.Random(seed * 1_000_003 + len(rows))
    arch = open(res_path, "a")

    def write(r, phase):
        r["idx"] = len(rows); r["phase"] = phase; rows.append(r); arch.write(json.dumps(r, default=str) + "\n")

    json.dump({"grammar": "g1", "arm": "quality-diversity v2 (DEVIATION 2)", "cells_per_combo": CELLS,
               "min_flagged": MIN_FLAGGED, "max_tries": MAX_TRIES, "budget": BUDGET, "patience": PATIENCE,
               "started_or_resumed": time.strftime("%Y-%m-%d %H:%M:%S")}, open(os.path.join(out, "manifest.json"), "w"), indent=1)
    with Pool(workers) as pool:
        if not os.path.exists(cvt_path):
            # 1. targeted initialization, per combination (resumes from rows already written)
            for combo in COMBOS:
                while True:
                    tried = sum(1 for r in rows if r["layers"] == combo)
                    fl = sum(1 for r in rows if r["layers"] == combo and r["status"] == "FLAGGED")
                    if fl >= MIN_FLAGGED or tried >= MAX_TRIES: break
                    batch = []
                    while len(batch) < INIT_BATCH:
                        q = rand_program_combo(rng, combo)
                        if q is not None: batch.append(G.encode(q))
                    for r in pool.imap_unordered(NV.evaluate, batch): write(r, "init")
                    arch.flush()
                    log(f"init {combo}: tried {tried + len(batch)}, flagged {sum(1 for r in rows if r['layers'] == combo and r['status'] == 'FLAGGED')}")
            cvt, info = build_cvt(rows, seed); save_cvt(cvt, info, cvt_path)
            log(f"CVT built: {info}  total cells {sum(v.get('cells', 0) for v in info.values())}")
        cvt, info = load_cvt(cvt_path)
        # 2. replay every row through the archive (deterministic given the CVT) -> elites, patience counter
        elites, since_new = {}, 0
        for r in rows:
            cell = cell_of(r, cvt)
            new = False
            if cell is not None:
                cur = elites.get(cell)
                if cur is None or (r["quality"], -r["len"]) > (cur["quality"], -cur["len"]):
                    new = cur is None; elites[cell] = r
            if r["phase"] == "search": since_new = 0 if new else since_new + 1
        log(f"resume: {len(rows)} rows, {len(elites)} niches filled of {sum(v.get('cells', 0) for v in info.values())}, {since_new} since last new niche")
        # 3. search (registered loop: mutate random elites, place, stop at budget or patience)
        while len(rows) < BUDGET and since_new < PATIENCE:
            parents = [G.decode(rng.choice(list(elites.values()))["bits"]) for _ in range(BATCH)] if elites else []
            kids = [G.encode(NV.mutate(p, rng)) for p in parents] if parents else [G.encode(NV.random_program(rng)) for _ in range(BATCH)]
            for r in pool.imap_unordered(NV.evaluate, kids):
                cell = cell_of(r, cvt); r["new_niche"] = False
                if cell is not None:
                    cur = elites.get(cell)
                    if cur is None or (r["quality"], -r["len"]) > (cur["quality"], -cur["len"]):
                        r["new_niche"] = cur is None; elites[cell] = r
                write(r, "search")
                since_new = 0 if r["new_niche"] else since_new + 1
            arch.flush()
            json.dump({f"{k[0]}:{k[1]}": {"quality": v["quality"], "len": v["len"], "prog": v["prog"], "bits": v["bits"]}
                       for k, v in elites.items()}, open(os.path.join(out, "elites.json"), "w"), indent=1)
            if len(rows) % 500 < BATCH:
                log(f"{len(rows)} evals, {len(elites)} niches, {since_new} since last new niche")
    stop = "budget" if len(rows) >= BUDGET else "patience"
    json.dump({"evaluations": len(rows), "niches_filled": len(elites),
               "cells": sum(v.get("cells", 0) for v in info.values()), "cvt_info": info, "stopped_on": stop,
               "finished": time.strftime("%Y-%m-%d %H:%M:%S")}, open(os.path.join(out, "SEARCH_DONE.json"), "w"), indent=1)
    log(f"SEARCH DONE: {len(rows)} evaluations, {len(elites)} niches, stopped on {stop}")


def step_downstream(out):
    """Name every flagged view (census namer), group unnamed + atypical-look views (fixed-radius), and run the
    literature check per distinct phenomenon (sonnet), with the strict opus second pass on any literature-silent or
    exact-combination-not-found verdict. Writes NOVELTY_V2.md."""
    import anthropic
    from census.namer import CensusNamer, FAMILY
    from census.triage import mean_fp, cluster
    from census import postrun as PR
    namer = CensusNamer(PR.LIB); rows = read_rows(os.path.join(out, "results_0.jsonl"))
    labs, beh = [], collections.defaultdict(list)
    for r in rows:
        if r.get("status") != "FLAGGED": continue
        for v, e in (r.get("views") or {}).items():
            if not e.get("flagged"): continue
            nm = ({"status": "NAMED", "behaviour": PR.AVAL, "also": [], "shown": [PR.AVAL]} if v == "C.aval"
                  else namer.name_view(v, e["seeds"]))
            L = {"idx": r["idx"], "len": r["len"], "bits": r["bits"], "prog": r["prog"], "grammar": "g1", "view": v,
                 "family": FAMILY.get(v, v), "status": nm["status"], "behaviour": nm.get("behaviour"),
                 "also": list(nm.get("also") or []), "shown": sorted(nm.get("shown") or []),
                 "fp": mean_fp(e["seeds"]) if v != "C.aval" else {}, "stage": "novelty_v2"}
            labs.append(L)
            for b in PR.behaviours_of(L): beh[b].append(L)
    json.dump([{k: v for k, v in L.items() if k != "fp"} for L in labs], open(os.path.join(out, "labels.json"), "w"))
    groups = collections.defaultdict(list)
    for L in labs:
        if L["status"] in ("UNNAMED", "KNOWN-BEHAVIOUR-ATYPICAL-LOOK"):
            groups[(L["status"], L["family"], tuple(L["shown"]))].append(L)
    phen = {}
    for key, items in groups.items():
        items.sort(key=lambda x: (x["len"], x["idx"]))
        lab = cluster([({"len": x["len"], "idx": x["idx"]}, x["view"], x["fp"]) for x in items], 7.0)
        for c in sorted(set(lab)):
            mem = [items[i] for i in np.flatnonzero(lab == c)]
            sig = tuple(sorted(set(re.findall(r"([CANF]\.[A-Z]+)\(", mem[0]["prog"]))))
            pk = (key[0], key[1], key[2], sig)
            p = phen.setdefault(pk, {"status": key[0], "rep": mem[0], "size": 0, "examples": [m["prog"] for m in mem[1:5]]})
            p["size"] += len(mem)
    plist = sorted(phen.values(), key=lambda p: (p["rep"]["len"], p["rep"]["idx"]))
    os.makedirs(os.path.join(out, "litcheck"), exist_ok=True)
    client = anthropic.Anthropic()
    for i, p in enumerate(plist, 1):
        fn = os.path.join(out, "litcheck", f"{i:03d}.json")
        if os.path.exists(fn): continue
        L = p["rep"]
        obs = ", ".join(f"{k}={round(v, 3)}" for k, v in sorted(L["fp"].items(), key=lambda kv: -abs(kv[1]))[:6])
        ctx = (f"\n\nIMPORTANT CONTEXT: this run passes the standard quantitative test for the known behaviour(s): "
               f"{', '.join(L['shown'])} -- but its overall fingerprint does not look like textbook examples of it."
               if p["status"] == "KNOWN-BEHAVIOUR-ATYPICAL-LOOK" and L["shown"] else "")
        card = (f"THE RULE SYSTEM (plain words):\n{PR.plain_rule(L)}\n\nSAME CLASS ALSO INCLUDES:\n" +
                "\n".join(f"  - {e}" for e in p["examples"]) + f"\n\nWHAT WAS OBSERVED (standout features; higher = "
                f"stronger): {obs}{ctx}\n\nThis appears across {p['size']} rule programs.")
        t = time.time()
        try:
            r = PR.ask(client, "claude-sonnet-5-5", PR.INSTR + "\n\n---\n\n" + card)
            if r.get("verdict") == "LITERATURE-SILENT" or (r.get("verdict") == "KNOWN-VARIANT" and PR.GRAY.search(str(r.get("notes", "")) + r.get("full_text", ""))):
                r2 = PR.ask(client, "claude-opus-5-5", PR.STRICT + "\n\n---\n\n" + card +
                            f"\n\nFIRST PASS said: {r.get('verdict')} / {r.get('standard_name')} -- {str(r.get('notes'))[:300]}")
                r = {"first_pass": {k: v for k, v in r.items() if k != "full_text"}, **r2, "second_pass": True}
        except Exception as e:
            r = {"verdict": "ERROR", "error": f"{type(e).__name__}: {e}"}
        r.update({"status_in_census": p["status"], "rep": L["prog"], "len": L["len"], "shown": L["shown"],
                  "size": p["size"], "checked": time.strftime("%Y-%m-%d %H:%M:%S")})
        json.dump(r, open(fn, "w"), indent=1)
        log(f"litcheck [{i}/{len(plist)}] {L['len']}b {p['status'][:7]} -> {r.get('verdict')} ({time.time() - t:.0f}s)")
        time.sleep(2)
    # summary
    sd = read_json(os.path.join(out, "SEARCH_DONE.json"))
    vs = [json.load(open(os.path.join(out, "litcheck", f))) for f in sorted(os.listdir(os.path.join(out, "litcheck")))]
    c = collections.Counter(v.get("verdict") for v in vs)
    census_behs = set()
    for st in ("g1", "g2"):
        tp = f"{PR.RUNS}/census_v1_{st}/triage.json"
        if os.path.exists(tp): census_behs |= set(json.load(open(tp))["behaviours"])
    lines = ["# Quality-diversity arm v2 (DEVIATION 2, post-registration, exploratory)", "",
             f"Search: {json.dumps(sd)}", "",
             f"Flagged views: {len(labs)}; status counts: {dict(collections.Counter(L['status'] for L in labs))}", "",
             "## Behaviours reached (shortest program)", "", "| behaviour | in census map? | shortest (bits) | programs |", "|---|---|---|---|"]
    for b, Ls in sorted(beh.items(), key=lambda kv: min(L["len"] for L in kv[1])):
        lines.append(f"| {b} | {'yes' if b in census_behs else '**NEW vs census**'} | {min(L['len'] for L in Ls)} | {len({L['idx'] for L in Ls})} |")
    lines += ["", "## Literature check of unnamed + atypical-look phenomena", "", f"- phenomena checked: {len(vs)}; verdicts: {dict(c)}",
              f"- **literature-silent after the strict second pass: {c.get('LITERATURE-SILENT', 0)}**", ""]
    for v in vs:
        if v.get("verdict") == "LITERATURE-SILENT":
            lines.append(f"  - CANDIDATE `{v['rep']}` — {v.get('undocumented_phenomenology')}  (needs prereg S4.4 vetting + owner review)")
    open(os.path.join(out, "NOVELTY_V2.md"), "w").write("\n".join(lines) + "\n")
    log("DOWNSTREAM DONE -> NOVELTY_V2.md")


def read_json(p):
    try: return json.load(open(p))
    except Exception: return {}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("step", choices=["search", "downstream"])
    ap.add_argument("--out", default="census/runs/census_v1_novelty_v2"); ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--seed", type=int, default=2)
    ap.add_argument("--smoke", action="store_true", help="tiny settings for a test run")
    a = ap.parse_args()
    if a.smoke:
        global CELLS, MIN_FLAGGED, MAX_TRIES, BUDGET, PATIENCE, INIT_BATCH
        CELLS, MIN_FLAGGED, MAX_TRIES, BUDGET, PATIENCE, INIT_BATCH = 3, 6, 24, 150, 40, 8
    step_search(a.out, a.workers, a.seed) if a.step == "search" else step_downstream(a.out)


if __name__ == "__main__":
    main()
