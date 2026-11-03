"""Score the planted positive controls under the rules fixed in PLANTED_CONTROLS_PLAN.md.

Each run's final verdict follows the census protocol (strict verdict where a strict pass ran, else the first pass).
N cards: LITERATURE-SILENT is correct; KNOWN is a miss; KNOWN-VARIANT is correct only if the reply states that the
reported behaviour itself is not described. A strict reply's `behaviour_documented` field decides that when present;
otherwise the case is listed for reading and decided in adjudication.json (with the quoted sentence). D cards: KNOWN, or
KNOWN-VARIANT naming the phenomenon (decided by reading), is correct; whether the reply cites the source paper itself is
reported alongside. A blocked-arm reply that cites the target paper (arXiv id, title, or an author surname with 2026) is a
leak and leaves the denominator. The leak detector is a screen: a reading decision in adjudication.json ("leak": "" or a
reason) overrides it, e.g. when the surname belongs to the same authors' earlier paper.

  python -m census.controls.score_planted     # writes runs/litcheck_planted/SCORES.md and scores.json
"""
import json, math, os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
PL = os.path.join(HERE, "planted"); RUNS = os.path.join(HERE, "..", "runs", "litcheck_planted")
ADJ = os.path.join(RUNS, "adjudication.json")


def cp(k, n):
    """exact two-sided 95% Clopper-Pearson interval"""
    if n == 0: return (float("nan"), float("nan"))
    def bisect(f, lo, hi):
        for _ in range(100):
            mid = (lo + hi) / 2
            if f(mid): hi = mid
            else: lo = mid
        return (lo + hi) / 2
    tail = lambda p, kk: sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(kk, n + 1))
    lo = 0.0 if k == 0 else bisect(lambda p: tail(p, k) >= 0.025, 0, 1)
    hi = 1.0 if k == n else bisect(lambda p: 1 - tail(p, k + 1) < 0.025, 0, 1)   # P(X <= k) falls below 0.025
    return (lo, hi)


def leak(v, t):
    blob = " ".join([v.get("full_text", "") or "", json.dumps(v.get("citations", "")), " ".join(v.get("result_urls", []))]).lower()
    if t["arxiv_id"].lower() in blob: return "arxiv id"
    title = re.sub(r"[^a-z0-9 ]", " ", t["title"].lower()); title = " ".join(title.split())
    if len(title) > 20 and title in " ".join(re.sub(r"[^a-z0-9 ]", " ", blob).split()): return "title"
    for a in t.get("authors") or []:
        sur = re.split(r"[ ,]", a.strip())[-1].lower() if isinstance(a, str) else ""
        if len(sur) > 3 and re.search(r"\b" + re.escape(sur) + r"\b[^.\n]{0,60}2026", blob): return f"author {sur} 2026"
    return ""


def final(arm_dir, key):
    f, s = os.path.join(arm_dir, f"{key}.first.json"), os.path.join(arm_dir, f"{key}.strict.json")
    if os.path.exists(s): return json.load(open(s)), "strict"
    if os.path.exists(f): return json.load(open(f)), "first"
    return None, None


def judge_N(v, adj, akey):
    vd = v.get("verdict")
    if vd == "LITERATURE-SILENT": return True, "verdict"
    if vd == "KNOWN": return False, "verdict"
    if vd == "KNOWN-VARIANT":
        bd = v.get("behaviour_documented")
        if isinstance(bd, bool) or str(bd).lower() in ("true", "false"):
            return str(bd).lower() == "false", "behaviour_documented"
        if akey in adj: return adj[akey]["correct"], "adjudicated"
        return None, "needs reading"
    return None, f"invalid ({vd})"


def main():
    adj = json.load(open(ADJ)) if os.path.exists(ADJ) else {}
    out = {"arms": {}, "pending": []}
    for setname, arms in (("N", ["N_blocked_r1", "N_blocked_r2", "N_blocked_r3", "N_open", "N_neutral", "N_strictall"]),
                          ("D", ["D_open"])):
        if not os.path.exists(os.path.join(PL, f"{setname}_cards.json")): continue
        cards = json.load(open(os.path.join(PL, f"{setname}_cards.json"))); truth = json.load(open(os.path.join(PL, f"{setname}_truth.json")))
        for arm in arms:
            d = os.path.join(RUNS, arm)
            if not os.path.isdir(d): continue
            rows = []
            for c in cards:
                key = c["class_id"].replace(":", "__"); v, src = final(d, key)
                if v is None: continue
                t = truth[c["class_id"]]; akey = f"{arm}|{c['class_id']}"
                lk = leak(v, t) if "blocked" in arm or arm in ("N_neutral", "N_strictall") else ""
                if lk and "leak" in adj.get(akey, {}): lk = adj[akey]["leak"]
                cites_target = bool(leak(v, t))
                if setname == "N":
                    ok, how = judge_N(v, adj, akey)
                else:
                    vd = v.get("verdict")
                    if vd == "KNOWN": ok, how = True, "verdict"
                    elif vd == "KNOWN-VARIANT" and akey in adj: ok, how = adj[akey]["correct"], "adjudicated"
                    elif vd == "KNOWN-VARIANT": ok, how = None, "needs reading"
                    elif vd == "LITERATURE-SILENT": ok, how = False, "verdict"
                    else: ok, how = None, f"invalid ({vd})"
                if ok is None and not lk and not how.startswith("invalid"): out["pending"].append(akey)
                rows.append({"card": c["class_id"], "verdict": v.get("verdict"), "source": src, "correct": ok, "how": how,
                             "leak": lk, "cites_target": cites_target, "searches": v.get("web_searches_used"), "stop_reason": v.get("stop_reason"),
                             "standard_name": str(v.get("standard_name", ""))[:120]})
            valid = [r for r in rows if not r["leak"] and r["correct"] is not None]
            k = sum(r["correct"] for r in valid)
            out["arms"][arm] = {"n_rows": len(rows), "valid": len(valid), "correct": k, "ci95": cp(k, len(valid)),
                                "leaks": sum(bool(r["leak"]) for r in rows), "cites_target": sum(r["cites_target"] for r in rows),
                                "invalid": sum(r["how"].startswith("invalid") for r in rows), "verdicts": collections.Counter(r["verdict"] for r in rows),
                                "rows": rows}
    # majority over the three blocked runs, per card (primary)
    maj = []
    for cid in sorted({r["card"] for a in ("N_blocked_r1", "N_blocked_r2", "N_blocked_r3") for r in out["arms"].get(a, {}).get("rows", [])}):
        vs = [r["correct"] for a in ("N_blocked_r1", "N_blocked_r2", "N_blocked_r3") for r in out["arms"].get(a, {}).get("rows", [])
              if r["card"] == cid and not r["leak"] and r["correct"] is not None]
        if vs: maj.append({"card": cid, "runs": len(vs), "correct_runs": sum(vs), "majority_correct": sum(vs) * 2 > len(vs),
                           "unanimous": len(set(vs)) == 1})
    k = sum(m["majority_correct"] for m in maj)
    pooled_n = sum(m["runs"] for m in maj); pooled_k = sum(m["correct_runs"] for m in maj)
    out["primary_pooled_runs"] = {"valid_runs": pooled_n, "correct_runs": pooled_k}
    out["primary_majority"] = {"cards": len(maj), "correct": k, "ci95": cp(k, len(maj)),
                               "unanimous": sum(m["unanimous"] for m in maj), "per_card": maj}
    json.dump(out, open(os.path.join(RUNS, "scores.json"), "w"), indent=1, default=lambda o: dict(o))
    L = ["# Planted positive controls: scores (PLANTED_CONTROLS_PLAN.md)", "",
         "| arm | runs | valid | correct | 95% CI | leaks | invalid | cites target | verdicts |",
         "|---|---|---|---|---|---|---|---|---|"]
    for arm, a in out["arms"].items():
        L.append(f"| {arm} | {a['n_rows']} | {a['valid']} | {a['correct']} | {a['ci95'][0]:.2f} to {a['ci95'][1]:.2f} | {a['leaks']} "
                 f"| {a['invalid']} | {a['cites_target']} | {dict(a['verdicts'])} |")
    pm, pr = out["primary_majority"], out["primary_pooled_runs"]
    L += ["", f"Primary (N blocked): {pr['correct_runs']} of {pr['valid_runs']} valid runs correct; by card majority over the "
          f"3 runs, {pm['correct']} of {pm['cards']} cards correct (exact 95% CI {pm['ci95'][0]:.2f} to {pm['ci95'][1]:.2f}); "
          f"unanimous across runs: {pm['unanimous']}.",
          f"Pending reading: {len(out['pending'])} {out['pending'][:20]}"]
    open(os.path.join(RUNS, "SCORES.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
