"""Simulation verification of the census literature verdicts (VERIFICATION_PLAN.md; post-registration, exploratory).

  python -m census.verify inventory              -> census/runs/verify/inventory.json (214 phenomena)
  python -m census.verify run --workers 6        -> census/runs/verify/sims/<pid>.json (resumable)
  python -m census.verify classify               -> census/runs/verify/verify.json + VERIFY.md
Frozen census code is imported unchanged; only system size and run length are varied, as the plan states.
"""
import argparse, collections, glob, hashlib, itertools, json, os, re, time
import numpy as np

OUT = "census/runs/verify"
FLAG = re.compile(r"contradict|inconsistent|does not match|doesn't match|not expected|unexpected|implausib|impossible|"
                  r"artifact|artefact|anomal|cannot produce|should not|would not produce|at odds|mismatch|opposite", re.I)
FAM_VIEWS = {"grid": ["C", "C.transient"], "agents": ["A"], "network": ["N"], "lattice-phase": ["C.phase"],
             "network-phase": ["N.phase"], "field": ["F0", "F1"]}
S_SEEDS = 10
CONDITIONS = {"census": dict(W=64, NA=400, NN=200, WA_ALONE=20, steps=1000),
              "area4": dict(W=128, NA=1600, NN=800, WA_ALONE=40, steps=1000),
              "len3": dict(W=64, NA=400, NN=200, WA_ALONE=20, steps=3000)}


def grammar(name):
    from census import grammar_g1, grammar_g2
    return grammar_g2 if name == "g2" else grammar_g1


def fresh_seed(bits):
    return int(hashlib.sha256(("verify:" + bits).encode()).hexdigest()[:8], 16)


# ------------------------------------------------------------------ inventory
TEXT_FIELDS = ("full_text", "notes", "phenomenon", "model_family", "standard_name", "undocumented_phenomenology")


def _text(*vs):
    # verdict prose only; citations and search strings are left out (a paper title must not trip the pattern)
    return " ".join(str(v.get(k) or "") for v in vs if v for k in TEXT_FIELDS)


def _card_features(fp):
    # the card builders listed the six largest-magnitude features of the representative's fingerprint
    return [k for k, _ in sorted(fp.items(), key=lambda kv: -abs(kv[1]))[:6]]


def _review_plist():
    """Rebuild census.postrun.step_review's phenomenon list (same inputs, same order) without calling a model."""
    from census import postrun as PR
    from census.triage import cluster
    import numpy as np
    by = collections.defaultdict(list)
    for stage in PR.STAGE_GRAMMAR:
        p = f"{PR.OUT}/labels_{stage}.jsonl"
        if not os.path.exists(p): continue
        for line in open(p):
            if "KNOWN-BEHAVIOUR-ATYPICAL-LOOK" not in line: continue
            L = json.loads(line)
            if L["status"] == "KNOWN-BEHAVIOUR-ATYPICAL-LOOK":
                L["stage"] = stage; by[(L["family"], tuple(L["shown"]))].append(L)
    groups = []
    for key, items in by.items():
        items.sort(key=lambda x: (x["len"], x["idx"]))
        lab = cluster([({"len": x["len"], "idx": x["idx"]}, x["view"], x["fp"]) for x in items], 7.0)
        for c in sorted(set(lab)):
            mem = [items[i] for i in np.flatnonzero(lab == c)]
            sig = tuple(sorted(set(re.findall(r"([CANF]\.[A-Z]+|[A-Z_]+\(\) -> [A-Z_]+)", mem[0]["prog"]))))
            groups.append({"id": f"{key[0]}|{'+'.join(key[1])}|{c}", "family": key[0], "shown": list(key[1]), "rep": mem[0], "sig": list(sig)})
    phen = {}
    for g in sorted(groups, key=lambda g: (g["rep"]["len"], g["rep"]["idx"])):
        pk = (g["family"], tuple(g["shown"]), tuple(g["sig"]))
        phen.setdefault(pk, {"groups": [], "rep": g["rep"]})
        phen[pk]["groups"].append(g["id"])
    return sorted(phen.values(), key=lambda p: (p["rep"]["len"], p["rep"]["idx"]))


def _qd_plist():
    """Rebuild census.novelty_v2's phenomenon list from its stored labels and result rows, without calling a model."""
    from census.novelty_v2 import read_rows
    from census.triage import mean_fp, cluster
    import numpy as np
    out = "census/runs/census_v1_novelty_v2"
    labs_stored = iter(json.load(open(os.path.join(out, "labels.json"))))
    labs = []
    for r in read_rows(os.path.join(out, "results_0.jsonl")):
        if r.get("status") != "FLAGGED": continue
        for v, e in (r.get("views") or {}).items():
            if not e.get("flagged"): continue
            L = dict(next(labs_stored))
            assert (L["idx"], L["view"]) == (r["idx"], v), "labels.json out of step with results_0.jsonl"
            L["fp"] = mean_fp(e["seeds"]) if v != "C.aval" else {}
            labs.append(L)
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
            phen.setdefault((key[0], key[1], key[2], sig), {"status": key[0], "rep": mem[0]})
    return sorted(phen.values(), key=lambda p: (p["rep"]["len"], p["rep"]["idx"]))


def inventory():
    R = "census/runs"
    items, problems = [], []
    reps = json.load(open(f"{R}/litcheck/representatives.json"))
    tri = {st: {c["id"]: c for c in json.load(open(f"{R}/census_v1_{st}/triage.json"))["unnamed_classes"]} for st in sorted({c["stage"] for c in reps})}
    rowbits = collections.defaultdict(dict)
    for c in reps:
        key = c["class_id"].replace(":", "__")
        v1 = json.load(open(f"{R}/litcheck/verdicts/{key}.json"))
        v2 = json.load(open(f"{R}/litcheck/secondpass/{key}.json")) if os.path.exists(f"{R}/litcheck/secondpass/{key}.json") else None
        cl = tri[c["stage"]][c["class_id"].split(":", 1)[1]]
        assert cl["shortest"] == c["program_terse"], c["class_id"]
        idx, view = cl["members"][0]
        rowbits[c["stage"]][idx] = None
        items.append({"pid": "census:" + c["class_id"], "batch": "census", "stage": c["stage"], "family": c["family"],
                      "rep": c["program_terse"], "idx": idx, "view": view,
                      "card_features": [o["feature"] for o in c.get("observed", [])],
                      "verdict": (v2 or v1).get("verdict"), "first_verdict": v1.get("verdict"), "text": _text(v1, v2)})
    # census representatives: bits and grammar from the stored result row of that program
    for st, want in rowbits.items():
        for fn in sorted(glob.glob(f"{R}/census_v1_{st}/results_*.jsonl")):
            for line in open(fn):
                try: r = json.loads(line)
                except Exception: continue
                if r.get("idx") in want: want[r["idx"]] = (r["bits"], r.get("grammar", "g2" if st == "g2" else "g1"))
    for it in items:
        it["bits"], it["grammar"] = rowbits[it["stage"]].get(it["idx"]) or (None, None)

    def add(batch, pid, v, L):
        fp = v.get("first_pass") if isinstance(v.get("first_pass"), dict) else None
        if L["prog"] != v["rep"] or list(L["shown"]) != list(v.get("shown") or []):
            problems.append(f"{pid}: rebuilt rep {L['prog']!r} {L['shown']} != stored {v['rep']!r} {v.get('shown')}")
        items.append({"pid": pid, "batch": batch, "stage": L.get("stage", "g1"), "family": L["family"], "rep": v["rep"],
                      "idx": L["idx"], "view": L["view"], "bits": L["bits"], "grammar": L.get("grammar", "g1"),
                      "card_features": _card_features(L["fp"]), "verdict": v.get("verdict"),
                      "first_verdict": (fp or v).get("verdict"), "text": _text(v, fp), "shown": v.get("shown")})

    rp = _review_plist(); files = sorted(glob.glob(f"{R}/postrun/review/*.json"))
    if len(files) != len(rp): problems.append(f"review: {len(files)} files vs {len(rp)} rebuilt phenomena")
    for f in files:
        i = int(os.path.basename(f)[:-5]); add("review", f"review:{i:03d}", json.load(open(f)), rp[i - 1]["rep"])
    qp = _qd_plist(); files = sorted(glob.glob(f"{R}/census_v1_novelty_v2/litcheck/*.json"))
    if len(files) != len(qp): problems.append(f"qd: {len(files)} files vs {len(qp)} rebuilt phenomena")
    for f in files:
        i = int(os.path.basename(f)[:-5]); add("qd", f"qd:{i:03d}", json.load(open(f)), qp[i - 1]["rep"])
    for it in items:
        it["flagged"] = bool(FLAG.search(it["text"] or ""))
        allf = it["card_features"] or []
        named = [f for f in allf if f in (it["text"] or "")]
        it["discrepant_features"] = named if (it["flagged"] and named) else None   # None -> all card features
    for p in problems: print("PROBLEM:", p)
    os.makedirs(OUT, exist_ok=True)
    json.dump(items, open(f"{OUT}/inventory.json", "w"), indent=1)
    c = collections.Counter((it["batch"], it["bits"] is not None, it["flagged"]) for it in items)
    print(len(items), "phenomena;", dict(c))


# ------------------------------------------------------------------ simulation
def _set_size(cond):
    from census import sim
    c = CONDITIONS[cond]
    sim.W, sim.NA, sim.NN, sim.WA_ALONE = c["W"], c["NA"], c["NN"], c["WA_ALONE"]
    return c["steps"]


def _fps(out_views, views):
    res = {}
    for v in views:
        if v in out_views:
            res[v] = [s["fp"] for s in out_views[v]]
    return res


def simulate(it):
    from census import sim, knockout as KO
    G = grammar(it["grammar"]); p = G.decode(it["bits"]); seed = fresh_seed(it["bits"])
    views = sorted(set(FAM_VIEWS.get(it["family"], ["C"])) | {it["view"]})
    rec = {"pid": it["pid"], "n_rules": len(p.rules), "conditions": {}, "subprograms": [], "errors": []}
    kp = KO.knockout_program(p)
    for cond in CONDITIONS:
        steps = _set_size(cond); t = time.time()
        try:
            r = KO.per_seed_view_results(sim.run(p, seed0=seed, steps=steps, S=S_SEEDS), p)
            k = KO.per_seed_view_results(sim.run(kp, seed0=seed, steps=steps, S=S_SEEDS), kp)
            rec["conditions"][cond] = {"real": _fps(r, views), "knockout": _fps(k, views), "secs": round(time.time() - t, 1)}
        except Exception as e:
            rec["errors"].append(f"{cond}: {type(e).__name__}: {e}")
    _set_size("census")
    if len(p.rules) >= 2:
        for n in range(1, len(p.rules)):
            for sub in itertools.combinations(range(len(p.rules)), n):
                try:
                    q = G.make(p.layers, p.k, p.nf, p.D, p.agent, [p.rules[i] for i in sub])
                    o = KO.per_seed_view_results(sim.run(q, seed0=seed, steps=1000, S=S_SEEDS), q)
                    rec["subprograms"].append({"rules": list(sub), "fps": _fps(o, views)})
                except Exception as e:
                    rec["subprograms"].append({"rules": list(sub), "error": f"{type(e).__name__}: {e}"})
    return rec


def _worker(it):
    fn = f"{OUT}/sims/{it['pid'].replace(':', '__')}.json"
    if os.path.exists(fn): return it["pid"], "cached"
    t = time.time()
    try: rec = simulate(it)
    except Exception as e: rec = {"pid": it["pid"], "fatal": f"{type(e).__name__}: {e}"}
    rec["secs_total"] = round(time.time() - t, 1)
    tmp = fn + ".tmp"; json.dump(rec, open(tmp, "w")); os.replace(tmp, fn)
    return it["pid"], f"{rec.get('secs_total')}s {('ERR ' + rec['fatal']) if 'fatal' in rec else ''}"


def run(workers, only=None):
    import multiprocessing as mp
    items = [it for it in json.load(open(f"{OUT}/inventory.json")) if it.get("bits")]
    if only: items = [it for it in items if only in it["pid"]]
    os.makedirs(f"{OUT}/sims", exist_ok=True)
    with mp.get_context("fork").Pool(workers, maxtasksperchild=1) as pool:
        for i, (pid, msg) in enumerate(pool.imap_unordered(_worker, items), 1):
            print(f"[{i}/{len(items)}] {pid} {msg}", flush=True)
    print("VERIFY RUN DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


# ------------------------------------------------------------------ classification (frozen rules)
def _col(fps_by_view, f, view):
    # one value per seed from the phenomenon's own view; seeds where the feature is undefined are dropped
    vals = [s[f] for s in (fps_by_view or {}).get(view, []) if s and s.get(f) is not None and np.isfinite(s[f])]
    return np.array(vals, float) if vals else None


def _scale(a, b):
    pooled = np.concatenate([a - np.median(a), b - np.median(b)])
    s = 1.4826 * np.median(np.abs(pooled))
    return max(s, 0.02 * (abs(np.median(a)) + abs(np.median(b))), 0.001)


def classify_one(it, rec):
    if "fatal" in rec or not rec.get("conditions", {}).get("census"):
        return {"class": "NOT RUNNABLE", "reason": rec.get("fatal") or "; ".join(rec.get("errors", []))}
    cen = rec["conditions"]["census"]
    view = it["view"]
    feats = it.get("discrepant_features") or it.get("card_features")
    per, robust = {}, []
    for f in feats:
        x, k = _col(cen["real"], f, view), _col(cen["knockout"], f, view)
        if x is None or k is None: per[f] = {"missing": True}; continue
        s = _scale(x, k); e = float(np.median(x) - np.median(k)); driven = abs(e) >= 3 * s
        pers = True; ctrl = {}
        for cond in ("area4", "len3"):
            c = rec["conditions"].get(cond)
            xc = _col(c["real"], f, view) if c else None; kc = _col(c["knockout"], f, view) if c else None
            if xc is None or kc is None: pers = False; ctrl[cond] = None; continue
            ec = float(np.median(xc) - np.median(kc)); ctrl[cond] = round(ec, 4)
            if not (np.sign(ec) == np.sign(e) and abs(ec) >= 0.5 * abs(e)): pers = False
        per[f] = {"median_real": round(float(np.median(x)), 4), "median_ko": round(float(np.median(k)), 4), "effect": round(e, 4),
                  "scale": round(s, 4), "interaction_driven": bool(driven), "persistent": bool(pers), "control_effects": ctrl}
        if driven and pers: robust.append(f)
    out = {"features": per, "robust": robust}
    if not robust: out["class"] = "ARTIFACT"; return out
    if rec["n_rules"] >= 2:
        best = None
        for sp in rec["subprograms"]:
            if "fps" not in sp: continue
            ok = True
            for f in robust:
                x, y = _col(cen["real"], f, view), _col(sp["fps"], f, view)
                if y is None or abs(np.median(x) - np.median(y)) > 3 * _scale(x, y): ok = False; break
            if ok and (best is None or len(sp["rules"]) < len(best)): best = sp["rules"]
        if best is not None: out.update({"class": "DECOMPOSES", "reproduced_by_rules": best})
        else: out["class"] = "COMBINATION EFFECT"
    else:
        out["class"] = "SINGLE MECHANISM"
    out["stage3"] = out["class"] == "COMBINATION EFFECT" or (out["class"] == "SINGLE MECHANISM" and it["flagged"])
    return out


def classify():
    items = {it["pid"]: it for it in json.load(open(f"{OUT}/inventory.json"))}
    res = []
    for pid, it in items.items():
        fn = f"{OUT}/sims/{pid.replace(':', '__')}.json"
        if not it.get("bits"): res.append({"pid": pid, "batch": it["batch"], "class": "NOT RUNNABLE", "reason": "bits not found"}); continue
        if not os.path.exists(fn): continue
        c = classify_one(it, json.load(open(fn)))
        res.append({"pid": pid, "batch": it["batch"], "rep": it["rep"], "family": it["family"], "flagged": it["flagged"],
                    "verdict": it["verdict"], **c})
    json.dump(res, open(f"{OUT}/verify.json", "w"), indent=1)
    tab = collections.Counter((r["batch"], r["class"]) for r in res)
    L = ["# Simulation verification of literature verdicts (VERIFICATION_PLAN.md)", "",
         f"{len(res)} phenomena classified; Stage 3 candidates: {sum(1 for r in res if r.get('stage3'))}", "",
         "| batch | class | n |", "|---|---|---|"] + [f"| {b} | {c} | {n} |" for (b, c), n in sorted(tab.items())]
    L += ["", "| phenomenon | rep | flagged | class | robust features | reproduced by rules |", "|---|---|---|---|---|---|"]
    for r in res:
        L.append(f"| {r['pid']} | `{r.get('rep','')}` | {r.get('flagged')} | {r['class']} | {', '.join(r.get('robust', []))} | {r.get('reproduced_by_rules', '')} |")
    open(f"{OUT}/VERIFY.md", "w").write("\n".join(L) + "\n"); print("\n".join(L[:12]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd"); ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only", default=None); a = ap.parse_args()
    {"inventory": inventory, "classify": classify}.get(a.cmd, lambda: run(a.workers, a.only))()
