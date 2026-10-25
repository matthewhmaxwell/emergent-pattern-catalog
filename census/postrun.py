"""Post-run REGISTERED analyses for the EPC Emergence Census (prereg S3, S4.5, S5, S10), on the frozen census code.

These were promised in PREREG_census_v1.md but not run before the results report was filed (found in the
2026-10-24 review). Runs unattended on the VPS; every step writes its own output and a failure in one step does not
stop the others.

  python -m census.postrun <step> [--limit N]
  steps: label | ffaudit | minimal | robust | review | summary

  label    name every flagged view of all three arms once; cache labels_<stage>.jsonl (used by the others)
  ffaudit  S5 false-flag audit: placebo test on a fresh random sample of 3,000 census programs of all lengths
  minimal  S10 10-seed recheck of each behaviour's minimal program (falls back to the next-shortest if it fails),
           S5 rule knock-outs on each confirmed G1 minimal program, S3 5,000-step confirmation of the shortest
           power-law-avalanche programs at the detector's confirmation tier
  robust   S5 robustness: for each confirmed minimal program, the fraction of all census programs with the same rule
           structure (substrates, type counts, field count, rule templates) that still produce the behaviour
  review   S4.5 known-behaviour-atypical-look views reviewed like unnamed groups: grouped by fingerprint
           (fixed-radius, radius 7), one literature check per distinct phenomenon, strict second pass on any
           literature-silent or exact-combination-not-found verdict
  summary  POSTRUN.md

Outputs: census/runs/postrun/.  API key for 'review': ANTHROPIC_API_KEY (source ~/.config/epc-litcheck.env).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, collections, glob, json, random, re, time
from multiprocessing import Pool
import numpy as np

RUNS = "census/runs"
OUT = "census/runs/postrun"
LIB = "census/reflib/census_v3/library.json"
STAGE_GRAMMAR = {"g1": "g1", "g2": "g2", "novelty": "g1"}
AVAL = "power-law avalanches"


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def load_rows(stage):
    rows = {}
    for f in sorted(glob.glob(f"{RUNS}/census_v1_{stage}/results_*.jsonl")):   # results_fix.jsonl sorts last -> wins
        for line in open(f):
            try:
                r = json.loads(line); rows[r["idx"]] = r
            except Exception:
                pass
    return list(rows.values())


def labels(stage):
    p = f"{OUT}/labels_{stage}.jsonl"
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


def grammar(name):
    from census.runner import grammar as g
    return g(name)


# ---------------------------------------------------------------- label
def step_label(limit=0):
    from census.namer import CensusNamer, FAMILY
    from census.triage import mean_fp
    namer = CensusNamer(LIB)
    for stage, gname in STAGE_GRAMMAR.items():
        rows = load_rows(stage)
        if limit: rows = rows[:limit]
        n = 0
        with open(f"{OUT}/labels_{stage}.jsonl", "w") as fh:
            for r in rows:
                if r.get("status") != "FLAGGED": continue
                for v, e in (r.get("views") or {}).items():
                    if not e.get("flagged"): continue
                    if v == "C.aval":
                        nm = {"status": "NAMED", "behaviour": AVAL, "also": [], "shown": [AVAL]}
                    else:
                        nm = namer.name_view(v, e["seeds"])
                    shown = nm.get("shown") or []
                    fh.write(json.dumps({"idx": r["idx"], "len": r["len"], "bits": r["bits"], "prog": r["prog"],
                        "grammar": r.get("grammar", gname), "layers": r["layers"], "view": v,
                        "family": FAMILY.get(v, v), "status": nm["status"], "behaviour": nm.get("behaviour"),
                        "also": list(nm.get("also") or []), "shown": sorted(shown),
                        "fp": mean_fp(e["seeds"]) if v != "C.aval" else {}}) + "\n")
                    n += 1
        log(f"label {stage}: {len(rows)} programs, {n} flagged views labelled")


def behaviours_of(L):
    return ({L["behaviour"]} if L.get("behaviour") else set()) | set(L.get("also") or []) | set(L.get("shown") or [])


# ---------------------------------------------------------------- ffaudit
def _placebo_one(job):
    from census.placebo import one
    return one(job)


def step_ffaudit(limit=0, n=3000, workers=5):
    from census.runner import program_list
    from census.namer import FAMILY
    g1 = program_list(22, "census/programs", "g1"); g2 = program_list(21, "census/programs", "g2")
    pool = [(b, "g1") for b in g1] + [(b, "g2") for b in g2]
    sample = random.Random(20261024).sample(pool, limit or n)
    jobs = [(b, g, "ffaudit-2026-10-24") for b, g in sample]
    t = time.time()
    with Pool(workers) as P:
        rows = list(P.imap_unordered(_placebo_one, jobs, chunksize=4))
    json.dump(rows, open(f"{OUT}/ffaudit.json", "w"))
    fam = collections.defaultdict(collections.Counter); flagged = []
    for r in rows:
        for v, e in r.get("views", {}).items():
            f = fam[FAMILY.get(v, v)]; f["views"] += 1; f["flagged"] += bool(e["flagged"])
            if e["flagged"]: flagged.append((r["prog"], v))
    nv = sum(f["views"] for f in fam.values()); nf = sum(f["flagged"] for f in fam.values())
    errs = sum("error" in r for r in rows)
    bound = f"rate <= {3 / nv:.3%} (95%, rule of three)" if nf == 0 and nv else f"rate = {nf / max(nv, 1):.3%}"
    res = {"programs": len(rows), "errors": errs, "views": nv, "false_flags": nf, "bound": bound,
           "by_family": {k: dict(v) for k, v in fam.items()}, "flagged": flagged[:50], "seconds": round(time.time() - t)}
    json.dump(res, open(f"{OUT}/ffaudit_summary.json", "w"), indent=1)
    log(f"ffaudit: {len(rows)} programs ({errs} errors), {nv} views, false flags {nf} -> {bound}")


# ---------------------------------------------------------------- minimal
def eval_S(p, bits, S):
    """runner.evaluate with S seeds instead of 3 (same flag rule, same textbook checks)."""
    from census.knockout import run_with_knockout
    from census.reflib import applicable_checks
    from census.runner import seed_of
    out, real, driven, _ = run_with_knockout(p, seed_of(bits), S=S)
    views = {}
    for v, seeds in real.items():
        flagged = sum(1 for s, x in enumerate(seeds) if x["emergent"] and driven[v][s]) >= 2
        checks = applicable_checks(p, v) if flagged else {}
        sl = []
        for s, x in enumerate(seeds):
            sd = {"emergent": bool(x["emergent"]), "evidence": float(x["evidence"])}
            if flagged and x["emergent"]:
                sd["fp"] = {k: float(val) for k, val in x["fp"].items()}; sd["checks"] = {}
                for c, fn in checks.items():
                    try: sd["checks"][c] = bool(fn(out, s)[0])
                    except Exception: sd["checks"][c] = False
            sl.append(sd)
        views[v] = {"flagged": bool(flagged), "n_em": sum(1 for x in seeds if x["emergent"]), "S": len(seeds),
                    "driven": bool(driven[v][0]) if seeds else False, "seeds": sl}
    return views


def shows(views, beh, namer):
    for v, e in views.items():
        if not e["flagged"]: continue
        if v == "C.aval":
            got = {AVAL}
        else:
            nm = namer.name_view(v, e["seeds"])
            got = ({nm["behaviour"]} if nm.get("behaviour") else set()) | set(nm.get("also") or []) | set(nm.get("shown") or [])
        if beh in got:
            return {"shown": True, "view": v, "emergent_seeds": e["n_em"], "S": e["S"], "driven": e["driven"]}
    best = max(views.values(), key=lambda e: e["n_em"], default=None)
    return {"shown": False, "emergent_seeds": best["n_em"] if best else 0, "S": best["S"] if best else 0}


def candidates(stage_list, beh, k):
    seen = {}; out = []
    for stage in stage_list:
        for L in labels(stage):
            if beh in behaviours_of(L) and (L["grammar"], L["bits"]) not in seen:
                seen[(L["grammar"], L["bits"])] = 1; out.append(L)
    out.sort(key=lambda L: (L["len"], L["idx"]))
    return out[:k]


def aval_confirm(L):
    from census import sim, filter as FL
    from census.runner import seed_of
    from epc.phase2a.panel import _verdict
    _, FULL = FL._load()
    p = grammar(L["grammar"]).decode(L["bits"])
    out = sim.run(p, seed0=seed_of(L["bits"]), steps=5000, rec_scale=5)
    tiers = []
    for s in range(3):
        for name, hist, md in FL.views(out, p, s):
            if name != "C.aval": continue
            try: tiers.append(_verdict(FULL["P14"](hist, md)))
            except Exception as e: tiers.append(f"error:{type(e).__name__}")
    ok = sum(t in ("confirmation", "definitive") for t in tiers)
    return {"prog": L["prog"], "len": L["len"], "grammar": L["grammar"], "tiers": tiers, "confirmed": ok >= 2}


def step_minimal(limit=0, S=10, k=3):
    from census.namer import CensusNamer
    namer = CensusNamer(LIB)
    G1 = grammar("g1")
    res = {"criterion": f"PASS = flagged under the census rule with {S} seeds (>=2 seeds emergent AND interaction-"
                        f"driven by the paired knock-out) AND the behaviour named/shown by the namer, AND emergent in "
                        f">= half of the {S} seeds", "behaviours": {}}
    behs = sorted({b for st in ("g1", "g2") for L in labels(st) for b in behaviours_of(L)})
    if limit: behs = behs[:limit]
    for beh in behs:
        res["behaviours"][beh] = {}
        for gname, stages in (("g1", ["g1"]), ("g2", ["g2"])):
            cands = candidates(stages, beh, k); tried = []; confirmed = None
            for L in cands:
                t = time.time()
                try:
                    v = eval_S(grammar(gname).decode(L["bits"]), L["bits"], S)
                    sh = shows(v, beh, namer)
                    sh["pass"] = bool(sh["shown"] and sh["emergent_seeds"] * 2 >= sh["S"])
                except Exception as e:
                    sh = {"pass": False, "error": f"{type(e).__name__}: {e}"}
                sh.update({"prog": L["prog"], "len": L["len"], "seconds": round(time.time() - t)})
                tried.append(sh)
                log(f"minimal {beh} [{gname}] {L['len']}b pass={sh.get('pass')} em={sh.get('emergent_seeds')}/{sh.get('S')}")
                if sh.get("pass"):
                    confirmed = L; break
            entry = {"census_minimum_len": cands[0]["len"] if cands else None, "tried": tried,
                     "confirmed_minimum": {"prog": confirmed["prog"], "len": confirmed["len"], "bits": confirmed["bits"]} if confirmed else None}
            # S5 rule knock-outs on the confirmed G1 minimal program
            if gname == "g1" and confirmed:
                p = G1.decode(confirmed["bits"]); kos = []
                if len(p.rules) == 1:
                    kos.append({"removed": p.rules[0].tmpl, "result": "sole rule -> interaction-free -> behaviour absent by construction (necessary)"})
                else:
                    for i, r in enumerate(p.rules):
                        q = G1.make(p.layers, p.k, p.nf, p.D, p.agent, [x for j, x in enumerate(p.rules) if j != i])
                        try:
                            qb = G1.encode(q); sh = shows(eval_S(q, qb, 3), beh, namer)
                            kos.append({"removed": r.tmpl, "behaviour_survives": sh["shown"], "necessary": not sh["shown"]})
                        except Exception as e:
                            kos.append({"removed": r.tmpl, "error": f"{type(e).__name__}: {e}"})
                entry["rule_knockouts"] = kos
            res["behaviours"][beh][gname] = entry
            json.dump(res, open(f"{OUT}/minimal.json", "w"), indent=1)
    # S3 avalanche confirmation: the three shortest avalanche programs per grammar, 5,000 steps, confirmation tier
    res["avalanche_confirmation"] = {}
    for gname in (() if limit else ("g1", "g2")):          # skipped in smoke tests (slow 5,000-step runs)
        rows = []
        for L in candidates([gname], AVAL, 3):
            try: r = aval_confirm(L)
            except Exception as e: r = {"prog": L["prog"], "error": f"{type(e).__name__}: {e}"}
            rows.append(r); log(f"avalanche {gname} {L['len']}b -> {r.get('tiers')} confirmed={r.get('confirmed')}")
        res["avalanche_confirmation"][gname] = rows
    json.dump(res, open(f"{OUT}/minimal.json", "w"), indent=1)
    log("minimal: done")


# ---------------------------------------------------------------- robust
def struct_key(prog, gname):
    header, _, rules = prog.partition(" | ")
    hp = header.split()
    ks = tuple(t for t in hp[1:] if re.match(r"k[CAN]=", t))
    nf = next((t for t in hp[1:] if t.startswith("fields=")), "fields=0")
    if gname == "g1":
        tm = tuple(re.findall(r"([CANF]\.[A-Z]+)\(", rules))
    else:
        tm = tuple(re.sub(r"\s+p=[0-9.e-]+", "", re.sub(r"\(.*?\)", "()", x)).strip() for x in rules.split(";"))
    return (hp[0], ks, nf, tm)


def step_robust(limit=0):
    mn = json.load(open(f"{OUT}/minimal.json"))
    res = {"definition": "fraction of all census programs with the same rule structure as the confirmed minimal "
                         "program (substrates, type counts, field count, rule templates; all rule/agent/diffusion "
                         "parameters free) whose run is flagged and shows the behaviour", "behaviours": {}}
    for gname, stage in (("g1", "g1"), ("g2", "g2")):
        tot = collections.Counter(struct_key(r["prog"], gname) for r in load_rows(stage) if "prog" in r)
        showing = collections.defaultdict(set)
        for L in labels(stage):
            key = struct_key(L["prog"], gname)
            for b in behaviours_of(L): showing[(b, key)].add(L["idx"])
        for beh, by in mn["behaviours"].items():
            c = (by.get(gname) or {}).get("confirmed_minimum")
            if not c: continue
            key = struct_key(c["prog"], gname); den = tot.get(key, 0); num = len(showing.get((beh, key), ()))
            res["behaviours"].setdefault(beh, {})[gname] = {"minimal": c["prog"], "structure": list(map(str, key)),
                "programs_with_structure": den, "showing_behaviour": num, "fraction": round(num / den, 3) if den else None}
            log(f"robust {beh} [{gname}] {num}/{den}")
    json.dump(res, open(f"{OUT}/robust.json", "w"), indent=1)
    log("robust: done")


# ---------------------------------------------------------------- review (S4.5)
INSTR = (
 "You are checking whether a behaviour produced by a very simple simulated rule system is already described in the "
 "scientific literature. You have a web_search tool -- use it, and prefer sources you can actually read. Most such "
 "behaviours are known; \"known\" is the expected answer and a perfectly good one. Work independently and skeptically.\n\n"
 "DO THIS:\n1. Identify the model family.\n2. Search for the SPECIFIC behaviour described (not just the family): is "
 "this behaviour, in this or an equivalent system, reported? Give the earliest clear source and one review if one "
 "exists.\n3. Verify every citation against what the search returns. Do not cite from memory alone.\n4. Decide exactly "
 "one verdict: KNOWN / KNOWN-VARIANT (say precisely what differs and whether it likely matters) / LITERATURE-SILENT "
 "(nothing describes it after a genuine search -- list the searches).\n5. Say how confident you are and why.\n\n"
 "End your reply with a fenced json block, exactly:\n```json\n{\"verdict\": \"KNOWN|KNOWN-VARIANT|LITERATURE-SILENT\", "
 "\"model_family\": \"...\", \"standard_name\": \"...\", \"confidence\": \"low|medium|high\", \"citations\": [\"...\"], "
 "\"searches\": [\"...\"], \"notes\": \"...\"}\n```")
STRICT = (
 "STRICT ADVERSARIAL SECOND PASS. A first pass could not settle this cleanly. Do not assume it was right. WHAT MATTERS "
 "IS THE BEHAVIOUR, NOT THE RULE TUPLE: many rule combinations produce the same known behaviour, so an undocumented "
 "exact rule combination is NOT a discovery. Use web_search and verify. Return KNOWN / KNOWN-VARIANT (a twist in the "
 "PHENOMENOLOGY you could not find) / LITERATURE-SILENT (the observed behaviour itself appears undocumented -- state "
 "exactly what is unaccounted for). Also say whether the features look like a genuine emergent pattern or a "
 "trivial/artifact outcome.\n\nEnd with a fenced json block exactly:\n```json\n{\"verdict\": \"...\", "
 "\"behaviour_documented\": true, \"trivial_or_artifact\": false, \"standard_name\": \"...\", \"citations\": [\"...\"], "
 "\"undocumented_phenomenology\": \"\", \"confidence\": \"...\", \"notes\": \"...\"}\n```")
GRAY = re.compile(r"not found|no source|no paper|not matched|no exact|wasn.t found|not located|could not find", re.I)


def plain_rule(L):
    try:
        if L["grammar"] == "g1":
            from census import plain
            return plain.describe(grammar("g1").decode(L["bits"]))
    except Exception:
        pass
    return L["prog"] + "  (compact grammar notation)"


def ask(client, model, text):
    msg = client.messages.create(model=model, max_tokens=8000,
        tools=[{"type": "web_search_20250305", "name": "web_search", "max_uses": 10}],
        messages=[{"role": "user", "content": text}])
    body = "".join(b.text for b in msg.content if getattr(b, "type", "") == "text")
    m = re.search(r"```json\s*(\{.*?\})\s*```", body, re.S)
    parsed = json.loads(m.group(1)) if m else {"verdict": "PARSE_FAILED"}
    parsed["searches_used"] = sum(1 for b in msg.content if getattr(b, "type", "") == "server_tool_use")
    parsed["full_text"] = body
    return parsed


def step_review(limit=0, radius=7.0):
    from census.triage import cluster
    import anthropic
    os.makedirs(f"{OUT}/review", exist_ok=True)
    by = collections.defaultdict(list)
    for stage in STAGE_GRAMMAR:
        for L in labels(stage):
            if L["status"] == "KNOWN-BEHAVIOUR-ATYPICAL-LOOK":
                L["stage"] = stage; by[(L["family"], tuple(L["shown"]))].append(L)
    groups = []
    for key, items in by.items():
        items.sort(key=lambda x: (x["len"], x["idx"]))
        lab = cluster([({"len": x["len"], "idx": x["idx"]}, x["view"], x["fp"]) for x in items], radius)
        for c in sorted(set(lab)):
            mem = [items[i] for i in np.flatnonzero(lab == c)]
            sig = tuple(sorted(set(re.findall(r"([CANF]\.[A-Z]+|[A-Z_]+\(\) -> [A-Z_]+)", mem[0]["prog"]))))
            groups.append({"id": f"{key[0]}|{'+'.join(key[1])}|{c}", "family": key[0], "shown": list(key[1]),
                           "size": len(mem), "rep": mem[0], "sig": list(sig), "examples": [m["prog"] for m in mem[1:5]]})
    # one literature check per distinct phenomenon (family, shown behaviours, rule signature of the representative)
    phen = {}
    for g in sorted(groups, key=lambda g: (g["rep"]["len"], g["rep"]["idx"])):
        pk = (g["family"], tuple(g["shown"]), tuple(g["sig"]))
        phen.setdefault(pk, {"groups": [], "rep": g["rep"], "size": 0, "examples": g["examples"]})
        phen[pk]["groups"].append(g["id"]); phen[pk]["size"] += g["size"]
    plist = sorted(phen.values(), key=lambda p: (p["rep"]["len"], p["rep"]["idx"]))
    json.dump({"atypical_views": sum(len(v) for v in by.values()), "groups": len(groups), "phenomena": len(plist),
               "list": [{"rep": p["rep"]["prog"], "len": p["rep"]["len"], "shown": p["rep"]["shown"], "size": p["size"],
                         "groups": p["groups"]} for p in plist]}, open(f"{OUT}/review_groups.json", "w"), indent=1)
    log(f"review: {sum(len(v) for v in by.values())} atypical views -> {len(groups)} groups -> {len(plist)} phenomena")
    if limit: plist = plist[:limit]
    client = anthropic.Anthropic()
    for i, p in enumerate(plist, 1):
        L = p["rep"]; fn = f"{OUT}/review/{i:03d}.json"
        if os.path.exists(fn): continue
        obs = ", ".join(f"{k}={round(v, 3)}" for k, v in sorted(L["fp"].items(), key=lambda kv: -abs(kv[1]))[:6])
        card = (f"THE RULE SYSTEM (plain words):\n{plain_rule(L)}\n\nSAME CLASS ALSO INCLUDES:\n" +
                "\n".join(f"  - {e}" for e in p["examples"]) +
                f"\n\nWHAT WAS OBSERVED (standout features; higher = stronger): {obs}\n\nIMPORTANT CONTEXT: this run "
                f"passes the standard quantitative test for the known behaviour(s): {', '.join(L['shown'])} -- but its "
                f"overall fingerprint does not look like textbook examples of it. The question is whether THIS form / "
                f"route of that behaviour (or anything else it is doing) is documented.\n\nThis appears across "
                f"{p['size']} rule programs.")
        t = time.time()
        try:
            r = ask(client, "claude-sonnet-5-5", INSTR + "\n\n---\n\n" + card)
            if r.get("verdict") == "LITERATURE-SILENT" or (r.get("verdict") == "KNOWN-VARIANT" and GRAY.search(str(r.get("notes", "")) + r.get("full_text", ""))):
                r2 = ask(client, "claude-opus-5-5", STRICT + "\n\n---\n\n" + card +
                         f"\n\nFIRST PASS said: {r.get('verdict')} / {r.get('standard_name')} -- {str(r.get('notes'))[:300]}")
                r = {"first_pass": {k: v for k, v in r.items() if k != "full_text"}, **r2, "second_pass": True}
        except Exception as e:
            r = {"verdict": "ERROR", "error": f"{type(e).__name__}: {e}"}
        r.update({"rep": L["prog"], "len": L["len"], "stage": L["stage"], "shown": L["shown"], "size": p["size"],
                  "groups": p["groups"], "checked": time.strftime("%Y-%m-%d %H:%M:%S")})
        json.dump(r, open(fn, "w"), indent=1)
        log(f"review [{i}/{len(plist)}] {L['len']}b {','.join(L['shown'])[:40]} -> {r.get('verdict')} "
            f"{'(2nd pass)' if r.get('second_pass') else ''} ({time.time() - t:.0f}s)")
        time.sleep(2)
    log("review: done")


# ---------------------------------------------------------------- summary
def step_summary(limit=0):
    def J(p):
        try: return json.load(open(f"{OUT}/{p}"))
        except Exception: return None
    L = ["# EPC Emergence Census — post-run registered analyses", "",
         f"Generated {time.strftime('%Y-%m-%d %H:%M')} on the frozen census code. These are the analyses promised in "
         "PREREG_census_v1.md (S3, S4.5, S5, S10) that were not run before the first results report.", ""]
    ff = J("ffaudit_summary.json")
    L += ["## S5 False-flag audit (placebo on fresh census programs)", ""]
    L += [f"- {ff['programs']} programs ({ff['errors']} errors), {ff['views']} interaction-free views compared with an "
          f"independent interaction-free run: **{ff['false_flags']} false flags** -> {ff['bound']}"] if ff else ["- not run / failed"]
    mn = J("minimal.json")
    L += ["", "## S10 Minimal programs re-checked on 10 seeds; S5 rule knock-outs", ""]
    if mn:
        L += [f"Criterion: {mn['criterion']}", "", "| behaviour | grammar | census minimum (bits) | confirmed minimum (bits) | confirmed program | knock-outs |", "|---|---|---|---|---|---|"]
        for beh, by in sorted(mn["behaviours"].items()):
            for g in ("g1", "g2"):
                e = by.get(g)
                if not e: continue
                c = e.get("confirmed_minimum"); ko = e.get("rule_knockouts")
                kos = "; ".join((k.get("result") or f"{k['removed']}: {'necessary' if k.get('necessary') else 'not necessary'}") for k in ko) if ko else ""
                L.append(f"| {beh} | {g} | {e.get('census_minimum_len')} | {c['len'] if c else 'NOT CONFIRMED'} | `{c['prog'] if c else ''}` | {kos} |")
        L += ["", "## S3 Avalanche confirmation (5,000 steps, confirmation tier)", ""]
        for g, rows in (mn.get("avalanche_confirmation") or {}).items():
            for r in rows:
                L.append(f"- [{g}] {r.get('len')}b `{r.get('prog')}` -> tiers {r.get('tiers')} confirmed={r.get('confirmed')}")
    else:
        L += ["- not run / failed"]
    rb = J("robust.json")
    L += ["", "## S5 Robustness of each minimal program", ""]
    if rb:
        L += [rb["definition"], "", "| behaviour | grammar | minimal program | programs with that structure | showing the behaviour | fraction |", "|---|---|---|---|---|---|"]
        for beh, by in sorted(rb["behaviours"].items()):
            for g, e in by.items():
                L.append(f"| {beh} | {g} | `{e['minimal']}` | {e['programs_with_structure']} | {e['showing_behaviour']} | {e['fraction']} |")
    else:
        L += ["- not run / failed"]
    rg = J("review_groups.json")
    L += ["", "## S4.5 Known-behaviour-atypical-look views, reviewed like unnamed groups", ""]
    if rg:
        vs = [json.load(open(f)) for f in sorted(glob.glob(f"{OUT}/review/*.json"))]
        c = collections.Counter(v.get("verdict") for v in vs)
        L += [f"- {rg['atypical_views']} atypical-look views -> {rg['groups']} fingerprint groups -> {rg['phenomena']} distinct phenomena; {len(vs)} literature-checked",
              f"- verdicts: {dict(c)}",
              f"- **literature-silent: {c.get('LITERATURE-SILENT', 0)}**", ""]
        for v in vs:
            if v.get("verdict") == "LITERATURE-SILENT":
                L.append(f"  - CANDIDATE: `{v['rep']}` ({', '.join(v['shown'])}) — {v.get('undocumented_phenomenology')}")
    else:
        L += ["- not run / failed"]
    open(f"{OUT}/POSTRUN.md", "w").write("\n".join(L) + "\n")
    log("summary: wrote POSTRUN.md")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("step"); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args(); os.makedirs(OUT, exist_ok=True)
    {"label": step_label, "ffaudit": step_ffaudit, "minimal": step_minimal, "robust": step_robust,
     "review": step_review, "summary": step_summary}[a.step](limit=a.limit)


if __name__ == "__main__":
    main()
