"""Planted behaviors through the whole census pipeline (PLANTED_PIPELINE_PLAN.md, commit 85acbe1; post-registration,
exploratory). Frozen census code is imported unchanged; the planted models come from census/planted_models/.

  python -m census.planted_pipeline flag [--only N02]     screen + interaction knock-out + textbook measures
  python -m census.planted_pipeline name                  census namer (library census_v3) on every flagged view
  python -m census.planted_pipeline cards                 census literature / review cards + leak check (no model call)
  python -m census.planted_pipeline gate --run r1         census protocol, target blocked (also r2, r3). Env: ANTHROPIC_API_KEY
  python -m census.planted_pipeline verify [--workers 6]  VERIFICATION_PLAN Stage 2 simulations (resumable)
  python -m census.planted_pipeline classify              frozen verify.classify_one
  python -m census.planted_pipeline stage3 prep|ask       blind Stage 3 (a) cards and the strict call, target blocked
  python -m census.planted_pipeline status                what exists so far

Readings of the plan made here, before any planted output existed:
- One program per behavior (the plan's headline setting); seeds: census seed_of(bits) for N17, base.seed_of(CID) for the
  plug-ins; verification seeds verify.fresh_seed(bits or CID).
- The knock-out is skipped when no view is emergent in >= 2 seeds, as census.knockout.run_with_knockout does.
- Sub-programs keep every proper non-empty subset of the INTERACTION rules (plan); "n_rules" given to the frozen
  classifier is therefore the number of interaction rules (N17: its one census rule), so a model with one interaction
  rule is SINGLE MECHANISM.
- Gate verdicts per card come from 3 runs. A card is mismatch-flagged (verify.FLAG on verify._text, as in the
  verification plan) when the majority of its valid runs match; discrepant features = card features named verbatim in
  the flagged runs' text, else all card features.
- Card id: pipeline:NXX:<view>. Family: census.namer.FAMILY, the view name itself where FAMILY has no entry (F2).
"""
import argparse, collections, glob, hashlib, importlib, itertools, json, os, re, sys, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np

OUT = "census/runs/planted_pipeline"
LIB = "census/reflib/census_v3/library.json"
IN_SCOPE = ["N01", "N02", "N03", "N04", "N07", "N08", "N11", "N12", "N13", "N14", "N15", "N17", "N19"]
OUT_OF_SCOPE = {"N05": "one-dimensional ring", "N06": "three-component spin vectors", "N09": "three-dimensional space",
                "N10": "aperiodic tiling", "N16": "three-dimensional space", "N18": "three-component spin vectors",
                "N20": "triangular lattice, open edges, imposed defect textures"}
MIRRORS = ["pith.science", "ar5iv.org", "papers.cool", "emergentmind.com", "chatpaper.com", "synthical.com",
           "scholar.archive.org", "paperdigest.org", "catalyzex.com", "arxiv-sanity-lite.com", "connectedpapers.com",
           "academia.edu"]
SHARED = {"github.com", "github.io", "gitlab.com", "sites.google.com", "google.com", "wordpress.com", "medium.com",
          "youtube.com", "zenodo.org", "osf.io", "figshare.com"}
RUNS = ("r1", "r2", "r3")
N_SEEDS, STEPS, VSEEDS = 3, 1000, 10


def J(path, default=None):
    try: return json.load(open(path))
    except Exception: return default


def dump(obj, path):
    os.makedirs(os.path.dirname(path), exist_ok=True); tmp = path + ".tmp"
    json.dump(obj, open(tmp, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)); os.replace(tmp, path)


def truth():
    return json.load(open("census/controls/planted/N_truth.json"))


def module(n):
    return importlib.import_module(f"census.planted_models.{n.lower()}")


def prog(mod, keep):
    if hasattr(mod, "prog"): return mod.prog(keep)                    # N17: the census Program itself
    from census.planted_models.base import PlugProg
    return PlugProg(mod, keep)


def seed(mod):
    if hasattr(mod, "BITS"):
        from census.runner import seed_of
        return seed_of(mod.BITS)
    from census.planted_models.base import seed_of
    return seed_of(mod.CID)


def interaction_rules(mod):
    return [i + 1 for i, (_, it) in enumerate(mod.RULES) if it]


# ---------------------------------------------------------------- 1 flag (census runner logic)
def flag_one(n):
    from census import knockout as KO
    from census.reflib import applicable_checks
    mod = module(n); s0 = seed(mod); p = prog(mod, None); t0 = time.time()
    out = mod.simulate(s0, steps=STEPS, S=N_SEEDS, size="census"); real = KO.per_seed_view_results(out, p); t1 = time.time()
    ko_ran = any(sum(1 for x in sd if x["emergent"]) >= 2 for sd in real.values())
    if ko_ran:
        kp = prog(mod, ()); ko = KO.per_seed_view_results(mod.simulate(s0, steps=STEPS, S=N_SEEDS, size="census", keep=()), kp)
        driven = KO.interaction_driven(real, ko)
    else:
        driven = {v: [False] * len(sd) for v, sd in real.items()}
    views = {}
    for v, sd in real.items():
        n_em = sum(1 for x in sd if x["emergent"])
        flagged = sum(1 for s, x in enumerate(sd) if x["emergent"] and driven[v][s]) >= 2
        row = {"flagged": bool(flagged), "em_seeds": n_em, "driven": bool(driven[v][0]) if sd else False, "seeds": []}
        checks = applicable_checks(p, v) if flagged else {}
        for s, x in enumerate(sd):
            d = {"emergent": x["emergent"], "evidence": round(float(x["evidence"]), 4)}
            if flagged and x["emergent"]:
                d["fp"] = {k: round(float(val), 5) for k, val in x["fp"].items()}; d["checks"] = {}
                for c, fn in checks.items():
                    try: d["checks"][c] = bool(fn(out, s)[0])
                    except Exception: d["checks"][c] = False
            row["seeds"].append(d)
        views[v] = row
    try: pres = mod.presence(out)
    except Exception as e: pres = {"error": f"{type(e).__name__}: {e}"}
    return {"card": n, "cid": mod.CID, "seed0": s0, "layers": mod.LAYERS if hasattr(mod, "LAYERS") else p.layers,
            "status": "FLAGGED" if any(x["flagged"] for x in views.values()) else "NONE", "views": views,
            "knockout_run": ko_ran, "presence": pres, "sim_s": round(t1 - t0, 1), "total_s": round(time.time() - t0, 1),
            "meta": {k: out["meta"].get(k) for k in ("steps_run", "frozen_at", "absorbed_at", "unstable")}}


def _flag_worker(n):
    fn = f"{OUT}/rows/{n}.json"
    if os.path.exists(fn): return n, "cached"
    try: r = flag_one(n)
    except Exception as e:
        import traceback
        r = {"card": n, "status": "ERROR", "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-1500:]}
    dump(r, fn); return n, f"{r['status']} {r.get('total_s', '')}s"


def flag(only=None, workers=4):
    import multiprocessing as mp
    items = [n for n in IN_SCOPE if not only or n == only]
    with mp.get_context("fork").Pool(workers, maxtasksperchild=1) as pool:
        for n, msg in pool.imap_unordered(_flag_worker, items): print("flag", n, msg, flush=True)


# ---------------------------------------------------------------- 2 name (census triage / post-run labels)
def name():
    from census.namer import CensusNamer, FAMILY
    from census.triage import mean_fp
    namer = CensusNamer(LIB); labels = []
    for n in IN_SCOPE:
        r = J(f"{OUT}/rows/{n}.json")
        if not r: continue
        for v, e in sorted((r.get("views") or {}).items()):
            if not e.get("flagged"): continue
            nm = namer.name_view(v, e["seeds"])
            labels.append({"card": n, "view": v, "family": FAMILY.get(v, v), "status": nm["status"],
                           "behaviour": nm.get("behaviour"), "also": list(nm.get("also") or []),
                           "shown": sorted(nm.get("shown") or []), "fp": mean_fp(e["seeds"])})
    dump(labels, f"{OUT}/labels.json")
    for L in labels: print(L["card"], L["view"], L["status"], L["behaviour"] or "", L["shown"])


# ---------------------------------------------------------------- 3 cards
COMMON = set("""a an and the of in on at to for with from by as or is are be its into over under two one three model
models system systems phase phases dynamics behavior behaviour emergence emergent pattern patterns active matter
particles particle random theory simulation simulations""".split())


def leak_terms(t):
    words = [w for w in re.findall(r"[a-z]+", t["title"].lower()) if len(w) >= 5 and w not in COMMON]
    surn = [re.split(r"[ ,]", a.strip())[-1].lower() for a in t.get("authors") or [] if isinstance(a, str)]
    return sorted(set(words)), sorted({s for s in surn if len(s) > 3}), [c.lower() for c in t.get("coined_terms") or []]


def leak_check(text, t):
    low = text.lower(); words, surn, coined = leak_terms(t); hits = []
    if t["arxiv_id"].lower() in low: hits.append("arxiv id")
    hits += [f"author:{s}" for s in surn if re.search(r"\b" + re.escape(s) + r"\b", low)]
    hits += [f"coined:{c}" for c in coined if c and c in low]
    hits += [f"title-word:{w}" for w in words if re.search(r"\b" + re.escape(w) + r"\b", low)]
    return hits


def review_card_text(c):
    # census.postrun.step_review's card, with the planted rule text and one member
    return (f"THE RULE SYSTEM (plain words):\n{c['program_plain']}\n\nSAME CLASS ALSO INCLUDES:\n" +
            "\n".join(f"  - {e}" for e in c.get("also_in_class_plain", [])) +
            f"\n\nWHAT WAS OBSERVED (standout features; higher = stronger): {c['obs_inline']}\n\nIMPORTANT CONTEXT: this run "
            f"passes the standard quantitative test for the known behaviour(s): {', '.join(c['shown'])} -- but its "
            f"overall fingerprint does not look like textbook examples of it. The question is whether THIS form / "
            f"route of that behaviour (or anything else it is doing) is documented.\n\nThis appears across "
            f"{c['size']} rule programs.")


def cards():
    from census import litcheck_runner as LR
    T = truth(); out = []
    for L in J(f"{OUT}/labels.json", []):
        if L["status"] == "NAMED": continue
        mod = module(L["card"]); fp = L["fp"]
        top = sorted(fp.items(), key=lambda kv: -abs(kv[1]))[:6]
        c = {"class_id": f"pipeline:{L['card']}:{L['view']}", "card": L["card"], "view": L["view"], "stage": "planted-pipeline",
             "family": L["family"], "size": 1, "program_plain": mod.RULE_TEXT, "also_in_class_plain": [],
             "observed": [{"feature": k, "value": round(v, 3)} for k, v in top], "kind": "review" if L["status"] != "UNNAMED" else "litcheck",
             "shown": L["shown"], "obs_inline": ", ".join(f"{k}={round(v, 3)}" for k, v in top)}
        c["text"] = LR.card_text(c) if c["kind"] == "litcheck" else review_card_text(c)
        c["leak_hits"] = leak_check(c["text"], T[f"planted:{L['card']}"])
        dump(c, f"{OUT}/cards/{c['class_id'].replace(':', '__')}.json"); out.append(c)
        print(c["class_id"], c["kind"], "LEAK HITS:" if c["leak_hits"] else "no leak hits", c["leak_hits"])
    print(len(out), "cards")


# ---------------------------------------------------------------- 4 gate
def blocked(t):
    from census import litcheck_planted as LP
    keep = [d for d in LP.blocked_for(t) if d not in SHARED and not any(d.endswith("." + s) for s in SHARED)]
    for m in MIRRORS:
        if m not in keep: keep.append(m)
    return keep[:64]


def _ask_blocked(client, model, text, blocked_hosts, max_uses=10, max_tokens=8000):
    # census.postrun.ask with the target blocked
    from census.litcheck_planted import _record
    msg = client.messages.create(model=model, max_tokens=max_tokens,
        tools=[{"type": "web_search_20250305", "name": "web_search", "max_uses": max_uses, "blocked_domains": blocked_hosts}],
        messages=[{"role": "user", "content": text}])
    return {"model": model, **_record(msg), "checked": time.strftime("%Y-%m-%d %H:%M:%S")}


def gate(run):
    import anthropic
    from census import litcheck_runner as LR, litcheck_planted as LP, postrun as PR
    client = anthropic.Anthropic(); T = truth(); d = f"{OUT}/gate/{run}"; os.makedirs(d, exist_ok=True)
    for fn in sorted(glob.glob(f"{OUT}/cards/*.json")):
        c = J(fn); key = c["class_id"].replace(":", "__"); t = T[f"planted:{c['card']}"]; b = blocked(t)
        print(f"[{run}] {c['class_id']} ({c['kind']})", flush=True)
        if c["kind"] == "litcheck":
            v = LP._save(f"{d}/{key}.first.json", lambda: LP.first_pass(client, c, b, LR.INSTR))
            if v.get("verdict") == "KNOWN-VARIANT":
                LP._save(f"{d}/{key}.strict.json", lambda: LP.strict_pass(client, c, v, b))
        else:
            v = LP._save(f"{d}/{key}.first.json", lambda: _ask_blocked(client, "claude-sonnet-5-5", PR.INSTR + "\n\n---\n\n" + c["text"], b))
            if v.get("verdict") == "LITERATURE-SILENT" or (v.get("verdict") == "KNOWN-VARIANT" and PR.GRAY.search(str(v.get("notes", "")) + v.get("full_text", ""))):
                LP._save(f"{d}/{key}.strict.json", lambda: _ask_blocked(client, "claude-opus-5-5", PR.STRICT + "\n\n---\n\n" + c["text"] +
                         f"\n\nFIRST PASS said: {v.get('verdict')} / {v.get('standard_name')} -- {str(v.get('notes'))[:300]}", b))
    print("GATE", run, "DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


def final_verdicts(key):
    """[(run, verdict dict, source, first dict)] for one card over the runs that exist."""
    res = []
    for run in RUNS:
        f, s = f"{OUT}/gate/{run}/{key}.first.json", f"{OUT}/gate/{run}/{key}.strict.json"
        first = J(f)
        if first is None: continue
        strict = J(s)
        res.append((run, strict or first, "strict" if strict else "first", first))
    return res


# ---------------------------------------------------------------- 5 verification (VERIFICATION_PLAN Stage 2)
def _fps(res, views):
    return {v: [s["fp"] for s in res[v]] for v in views if v in res}


def vsim(job):
    from census import verify as V, knockout as KO
    n, kind, arg = job; mod = module(n); fs = V.fresh_seed(getattr(mod, "BITS", mod.CID)); views = job_views(n)
    if kind == "cond":
        steps = V._set_size(arg); size = "area4" if arg == "area4" else "census"
        try:
            r = KO.per_seed_view_results(mod.simulate(fs, steps=steps, S=VSEEDS, size=size), prog(mod, None))
            k = KO.per_seed_view_results(mod.simulate(fs, steps=steps, S=VSEEDS, size=size, keep=()), prog(mod, ()))
            return {"real": _fps(r, views), "knockout": _fps(k, views)}
        finally:
            V._set_size("census")
    V._set_size("census")
    o = KO.per_seed_view_results(mod.simulate(fs, steps=STEPS, S=VSEEDS, size="census", keep=tuple(arg)), prog(mod, tuple(arg)))
    return {"rules": list(arg), "fps": _fps(o, views)}


def job_views(n):
    from census import verify as V
    from census.namer import FAMILY
    vs = set()
    for c in glob.glob(f"{OUT}/cards/pipeline__{n}__*.json"):
        cc = J(c); vs |= set(V.FAM_VIEWS.get(FAMILY.get(cc["view"], cc["view"]), [])) | {cc["view"]}
    return sorted(vs)


def _vworker(job):
    n, kind, arg = job; tag = arg if kind == "cond" else "sub_" + "-".join(map(str, arg))
    fn = f"{OUT}/verify/parts/{n}__{tag}.json"
    if os.path.exists(fn): return job, "cached"
    t = time.time()
    try: r = vsim(job)
    except Exception as e: r = {"error": f"{type(e).__name__}: {e}"}
    r["secs"] = round(time.time() - t, 1); dump(r, fn)
    return job, f"{r['secs']}s {'ERR ' + r['error'] if 'error' in r else ''}"


def verify(workers=6, only=None):
    import multiprocessing as mp
    jobs = []
    for n in IN_SCOPE:
        if only and n != only: continue
        if not glob.glob(f"{OUT}/cards/pipeline__{n}__*.json"): continue          # no gate verdict -> not verified
        if not any(final_verdicts(J(c)["class_id"].replace(":", "__")) for c in glob.glob(f"{OUT}/cards/pipeline__{n}__*.json")): continue
        mod = module(n); inter = interaction_rules(mod)
        jobs += [(n, "cond", c) for c in ("census", "area4", "len3")]
        if len(inter) >= 2:
            jobs += [(n, "sub", sub) for k in range(1, len(inter)) for sub in itertools.combinations(inter, k)]
    print(len(jobs), "verification jobs", flush=True)
    with mp.get_context("fork").Pool(workers, maxtasksperchild=1) as pool:
        for job, msg in pool.imap_unordered(_vworker, jobs): print("verify", job, msg, flush=True)


def assemble(n):
    mod = module(n); rec = {"pid": n, "n_rules": len(interaction_rules(mod)) if not hasattr(mod, "BITS") else len(mod.prog(None).rules),
                            "conditions": {}, "subprograms": [], "errors": []}
    for c in ("census", "area4", "len3"):
        p = J(f"{OUT}/verify/parts/{n}__{c}.json")
        if p is None: rec["errors"].append(f"{c}: not run")
        elif "error" in p: rec["errors"].append(f"{c}: {p['error']}")
        else: rec["conditions"][c] = {"real": p["real"], "knockout": p["knockout"], "secs": p.get("secs")}
    for f in sorted(glob.glob(f"{OUT}/verify/parts/{n}__sub_*.json")):
        p = J(f); rec["subprograms"].append({"rules": p.get("rules"), "fps": p["fps"]} if "fps" in p else {"rules": p.get("rules"), "error": p.get("error")})
    return rec


def mismatch(key, features):
    from census import verify as V
    fv = final_verdicts(key); flags = []
    for run, v, src, first in fv:
        if v.get("verdict") in ("ERROR", "PARSE_FAILED", None): continue
        txt = V._text(first, v if src == "strict" else None); flags.append((run, bool(V.FLAG.search(txt)), txt))
    fl = sum(f for _, f, _ in flags) * 2 > len(flags) if flags else False
    named = sorted({f for _, ok, txt in flags if ok for f in features if f in txt}) if fl else []
    return fl, named


def classify():
    from census import verify as V
    res = []
    for fn in sorted(glob.glob(f"{OUT}/cards/*.json")):
        c = J(fn); key = c["class_id"].replace(":", "__")
        if not final_verdicts(key): continue
        feats = [o["feature"] for o in c["observed"]]; fl, named = mismatch(key, feats)
        it = {"pid": c["class_id"], "view": c["view"], "card_features": feats, "discrepant_features": named or None, "flagged": fl}
        rec = assemble(c["card"]); cl = V.classify_one(it, rec)
        res.append({"pid": c["class_id"], "card": c["card"], "view": c["view"], "flagged": fl, "discrepant": named or feats, **cl})
        print(c["class_id"], cl["class"], "stage3" if cl.get("stage3") else "", cl.get("robust"))
    dump(res, f"{OUT}/verify/verify.json")


# ---------------------------------------------------------------- 6 Stage 3 (a)
def stage3_card(c, cl, rec, mod):
    # census.verify_stage3.card(..., BLIND) with the planted rule text and the plan's rule numbers
    from census import litcheck_runner as LR, verify as V
    view = c["view"]; lines = []
    for f in cl["robust"]:
        d = cl["features"][f]; g = LR.GLOSS.get(f); gl = f" ({g})" if g else ""
        if d.get("knockout_view_absent"):
            lines.append(f"  - {f}{gl}: median {d['median_real']} with interaction; with interaction knocked out the program "
                         f"never produces this view, at the original size, at 4x area or at 3x run length"); continue
        lines.append(f"  - {f}{gl}: median {d['median_real']} with interaction, {d['median_ko']} with interaction knocked out "
                     f"(difference {d['effect']}, robust spread {d['scale']}); difference at 4x area {d['control_effects'].get('area4')}, "
                     f"at 3x run length {d['control_effects'].get('len3')}")
    subs = []
    for sp in rec.get("subprograms", []):
        if "fps" not in sp: continue
        med = lambda f: (lambda a: None if a is None else round(float(np.median(a)), 4))(V._col(sp["fps"], f, view))
        subs.append(f"  - interaction rules {', '.join(map(str, sp['rules']))} only: " + "; ".join(f"{f} {med(f)}" for f in cl["robust"]))
    ev = ("SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, "
          "averaged over the late part of the run. The interaction knock-out is the same program with every rule that "
          "reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase "
          "coupling is set to zero). Unless stated otherwise below, each listed feature differs from its knock-out by at least 3 robust spreads at the "
          "original size, and keeps its sign and at least half its size at 4x area and at 3x run length.\n" + "\n".join(lines))
    if subs:
        names = "; ".join(f"{i + 1} {nm}" for i, (nm, it) in enumerate(mod.RULES) if it)
        ev += (f"\n\nSUB-PROGRAMS (same program with some interaction rules removed, original size, 10 seeds; median of each "
               f"feature above; interaction rules: {names}):\n" + "\n".join(subs))
        if cl["class"] == "COMBINATION EFFECT": ev += "\nNo sub-program came within 3 robust spreads of the full program on every feature above."
    return (f"RULE (plain words):\n{c['program_plain']}\n\nOBSERVED (features that passed the simulation checks):\n  "
            + ", ".join(cl["robust"]) + "\n\n" + ev)


def stage3(cmd):
    from census import litcheck_runner as LR, litcheck_secondpass as SP, litcheck_planted as LP
    d = f"{OUT}/stage3"; T = truth()
    if cmd == "prep":
        os.makedirs(f"{d}/cards", exist_ok=True); n = 0
        for r in J(f"{OUT}/verify/verify.json", []):
            if not r.get("stage3"): continue
            c = J(f"{OUT}/cards/{r['pid'].replace(':', '__')}.json"); mod = module(c["card"])
            txt = stage3_card(c, r, assemble(c["card"]), mod)
            hits = leak_check(txt, T[f"planted:{c['card']}"])
            open(f"{d}/cards/{r['pid'].replace(':', '__')}.txt", "w").write(txt); n += 1
            print(r["pid"], "LEAK HITS" if hits else "", hits)
        print(n, "stage 3 cards")
        return
    import anthropic
    client = anthropic.Anthropic(); os.makedirs(f"{d}/verdicts", exist_ok=True)
    LABEL = {"KNOWN": "DOCUMENTED", "KNOWN-VARIANT": "CLOSELY RELATED", "LITERATURE-SILENT": "UNDOCUMENTED CANDIDATE"}
    for fn in sorted(glob.glob(f"{d}/cards/*.txt")):
        key = os.path.basename(fn)[:-4]; card = key.split("__")[1]; b = blocked(T[f"planted:{card}"])
        txt = open(fn).read()
        r = LP._save(f"{d}/verdicts/{key}.json", lambda: _ask_blocked(client, "claude-opus-5-5", SP.INSTR + "\n\n---\n\n" + txt, b,
                                                                      max_uses=12, max_tokens=4500))
        r["provisional_label"] = LABEL.get(r.get("verdict"), "UNRESOLVED"); dump(r, f"{d}/verdicts/{key}.json")


# ---------------------------------------------------------------- report (plan "Scoring")
def gate_score(c, T, adj):
    """per run: (run, final verdict, correct or None, leak, how); card majority over valid runs."""
    from census.controls.score_planted import leak as leak_fn
    key = c["class_id"].replace(":", "__"); t = T[f"planted:{c['card']}"]; runs = []
    for run, v, src, first in final_verdicts(key):
        a = adj.get(f"{run}|{c['class_id']}", {}); vd = v.get("verdict")
        lk = leak_fn(v, t) or (leak_fn(first, t) if src == "strict" else "")
        if lk and "leak" in a: lk = a["leak"]
        if vd == "LITERATURE-SILENT": ok, how = True, "verdict"
        elif vd == "KNOWN": ok, how = False, "verdict"
        elif vd == "KNOWN-VARIANT": ok, how = (a["correct"], "read") if "correct" in a else (None, "needs reading")
        else: ok, how = None, f"invalid ({vd})"
        runs.append({"run": run, "verdict": vd, "source": src, "correct": ok, "how": how, "leak": lk,
                     "close": bool(a.get("close")), "silent": vd == "LITERATURE-SILENT"})
    valid = [r for r in runs if not r["leak"] and r["correct"] is not None]
    maj = (sum(r["correct"] for r in valid) * 2 > len(valid)) if valid else None
    silent_maj = (sum(r["silent"] for r in valid) * 2 > len(valid)) if valid else None
    return runs, maj, silent_maj


def cp(k, n):
    from census.controls.score_planted import cp as _cp
    return _cp(k, n)


def report():
    T = truth(); adj = J(f"{OUT}/reading/decisions.json", {}); ver = {r["pid"]: r for r in J(f"{OUT}/verify/verify.json", [])}
    s3 = {os.path.basename(f)[:-5].replace("__", ":"): J(f) for f in glob.glob(f"{OUT}/stage3/verdicts/*.json")}
    s3b = J(f"{OUT}/stage3_review/final_labels.json", {})
    fid = {n: J(f"{OUT}/fidelity/{n}.json") for n in IN_SCOPE}
    labels = J(f"{OUT}/labels.json", []); rows = []
    for i in range(1, 21):
        n = f"N{i:02d}"
        if n in OUT_OF_SCOPE:
            rows.append({"card": n, "scope": "out", "reason": OUT_OF_SCOPE[n], "lost_at": "engine scope"}); continue
        r = J(f"{OUT}/rows/{n}.json") or {}; L = [x for x in labels if x["card"] == n]
        cs = [J(f) for f in sorted(glob.glob(f"{OUT}/cards/pipeline__{n}__*.json"))]
        e = {"card": n, "scope": "in", "flagged": sorted(v for v, x in (r.get("views") or {}).items() if x.get("flagged")),
             "labels": [(x["view"], x["status"], x["behaviour"]) for x in L], "presence": r.get("presence"),
             "fidelity": (fid.get(n) or {}).get("pass"), "cards": []}
        surv_a = surv_b = surv_lenient = False
        for c in cs:
            runs, maj, smaj = gate_score(c, T, adj); v = ver.get(c["class_id"], {}); st = s3.get(c["class_id"], {})
            fin = s3b.get(c["class_id"], {}).get("final_label") if isinstance(s3b, dict) else None
            e["cards"].append({"id": c["class_id"], "kind": c["kind"], "features": [o["feature"] for o in c["observed"]],
                               "runs": runs, "majority_correct": maj, "majority_silent": smaj,
                               "verify_class": v.get("class"), "stage3": v.get("stage3"), "stage3a": st.get("provisional_label"),
                               "stage3b": fin})
            surv_a |= bool(smaj); surv_lenient |= bool(maj); surv_b |= fin == "UNDOCUMENTED CANDIDATE"
        if not r: lost = "not run"
        elif not e["flagged"]: lost = "flag"
        elif not cs: lost = "name"
        elif not surv_lenient: lost = "gate"
        else: lost = "survives (gate)"
        e.update({"survives_registered": surv_a, "survives_lenient": surv_lenient, "survives_verification": surv_b, "lost_at": lost})
        rows.append(e)
    dump(rows, f"{OUT}/report.json")
    ins = [x for x in rows if x["scope"] == "in"]
    cnt = lambda f: sum(1 for x in ins if f(x))
    L = ["# Planted behaviors through the whole census pipeline (PLANTED_PIPELINE_PLAN.md)", "",
         f"In engine scope: {len(ins)} of 20. Flagged: {cnt(lambda x: x['flagged'])}. Reached a card (not named as known): "
         f"{cnt(lambda x: x['cards'])}. Gate majority correct (plan's reading rule): {cnt(lambda x: x['survives_lenient'])}. "
         f"LITERATURE-SILENT majority (registered path): {cnt(lambda x: x['survives_registered'])}. "
         f"Undocumented candidate after Stage 3: {cnt(lambda x: x['survives_verification'])}.", ""]
    for lab, f in (("registered path", "survives_registered"), ("verification path", "survives_verification"), ("lenient gate rule", "survives_lenient")):
        k = cnt(lambda x: x[f]); a, b = cp(k, 20); c_, d_ = cp(k, len(ins))
        L.append(f"- End-to-end, {lab}: {k} of 20 (95% CI {a:.2f} to {b:.2f}); {k} of {len(ins)} in scope ({c_:.2f} to {d_:.2f})")
    L += ["", "| card | scope | fidelity | flagged views | names | cards: gate runs (final verdicts) | verification | Stage 3 | lost at |",
          "|---|---|---|---|---|---|---|---|---|"]
    for x in rows:
        if x["scope"] == "out":
            L.append(f"| {x['card']} | out ({x['reason']}) | | | | | | | engine scope |"); continue
        g = "; ".join(f"{c['id'].split(':')[-1]}: " + ", ".join(f"{r['verdict']}{'(leak)' if r['leak'] else ''}{'+' if r['correct'] else ''}" for r in c["runs"]) for c in x["cards"])
        vc = "; ".join(str(c["verify_class"]) for c in x["cards"]); s3t = "; ".join(f"{c['stage3a'] or ''}/{c['stage3b'] or ''}" for c in x["cards"] if c["stage3"])
        nm = "; ".join(v + ":" + s + ((" " + b) if b else "") for v, s, b in x["labels"])
        L.append(f"| {x['card']} | in | {x['fidelity']} | {', '.join(x['flagged'])} | {nm} | {g} | {vc} | {s3t} | {x['lost_at']} |")
    open(f"{OUT}/RESULTS.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


def status():
    for n in IN_SCOPE:
        r = J(f"{OUT}/rows/{n}.json")
        print(n, r.get("status") if r else "-", {v: x["flagged"] for v, x in (r or {}).get("views", {}).items()})


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd"); ap.add_argument("sub", nargs="?")
    ap.add_argument("--only"); ap.add_argument("--run", choices=RUNS); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.cmd == "flag": flag(a.only, a.workers)
    elif a.cmd == "name": name()
    elif a.cmd == "cards": cards()
    elif a.cmd == "gate": gate(a.run)
    elif a.cmd == "verify": verify(a.workers, a.only)
    elif a.cmd == "classify": classify()
    elif a.cmd == "stage3": stage3(a.sub)
    elif a.cmd == "status": status()
    elif a.cmd == "report": report()
