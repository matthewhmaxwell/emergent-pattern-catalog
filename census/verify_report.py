"""Reporting for VERIFICATION_PLAN.md (no class is changed; reads census/runs/verify/verify.json).

Adds what the plan reports but does not classify on: whether each tested feature replicates its literature-card value
on the fresh seeds, and cross-tabulations of class against batch, flag and verdict. Writes
census/runs/verify/VERIFY_REPORT.md (summary) and VERIFY_TABLE.md (every phenomenon, every tested feature).

Card values: census batch from census/runs/litcheck/representatives.json ("observed"); post-run review and
quality-diversity batches from each rebuilt representative's stored fingerprint (verify._review_plist / _qd_plist,
no model call). Replication uses the plan's rule |median(R) - card| <= max(3 s, 0.25 |card|), with s from R's own ten
census-setting seeds (the card is one stored value, so there is no second seed set to pool), floored as in the plan.

  python -m census.verify_report
"""
import collections, glob, json, os
import numpy as np
from census import verify as V


def _card_values():
    R = "census/runs"; out = {}
    for c in json.load(open(f"{R}/litcheck/representatives.json")):
        out["census:" + c["class_id"]] = {o["feature"]: o["value"] for o in c.get("observed", [])}
    rp = V._review_plist()
    for f in sorted(glob.glob(f"{R}/postrun/review/*.json")):
        i = int(os.path.basename(f)[:-5]); out[f"review:{i:03d}"] = rp[i - 1]["rep"]["fp"]
    qp = V._qd_plist()
    for f in sorted(glob.glob(f"{R}/census_v1_novelty_v2/litcheck/*.json")):
        i = int(os.path.basename(f)[:-5]); out[f"qd:{i:03d}"] = qp[i - 1]["rep"]["fp"]
    return out


def _replicates(x, card):
    m = float(np.median(x)); s = 1.4826 * float(np.median(np.abs(x - m)))
    s = max(s, 0.02 * (abs(m) + abs(card)), 0.001)
    return abs(m - card) <= max(3 * s, 0.25 * abs(card))


def main():
    res = json.load(open(f"{V.OUT}/verify.json"))
    inv = {it["pid"]: it for it in json.load(open(f"{V.OUT}/inventory.json"))}
    cards = _card_values()
    rep_n = rep_ok = 0; rep_by_class = collections.Counter(); rep_tot_class = collections.Counter()
    T = ["# Simulation verification: every phenomenon and tested feature", "",
         "Medians over 10 fresh seeds at census settings. e = median(real) - median(knock-out); s = robust spread; "
         "e/s in brackets. 4x and 3x are the effects at 4x area and 3x run length. 'ko absent' = the knock-out never "
         "produces the view (deviation V2). card = value on the literature card; rep = replicates the card (reported "
         "only). R = robust (interaction-driven and persistent).", ""]
    for r in res:
        it = inv[r["pid"]]; cv = cards.get(r["pid"], {})
        rec = json.load(open(f"{V.OUT}/sims/{r['pid'].replace(':', '__')}.json"))
        T += [f"## {r['pid']}  ({r['class']}{', ' + r['reason'] if r.get('reason') else ''})", "",
              f"`{r.get('rep', '')}` · view {it['view']} · flagged {r.get('flagged')} · verdict {r.get('verdict')}"
              + (f" · reproduced by rules {[i + 1 for i in r['reproduced_by_rules']]}" if r.get("reproduced_by_rules") else ""), "",
              "| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |", "|---|---|---|---|---|---|---|---|---|"]
        for f, d in (r.get("features") or {}).items():
            x = V._col(rec["conditions"]["census"]["real"], f, it["view"])
            card = cv.get(f); rep = None
            if x is not None and card is not None and np.isfinite(card):
                rep = _replicates(x, float(card)); rep_n += 1; rep_ok += rep
                rep_tot_class[r["class"]] += 1; rep_by_class[r["class"]] += rep
            cs = "" if card is None else f"{card:.4g}"
            rs = "" if rep is None else ("yes" if rep else "no")
            R = "R" if f in r.get("robust", []) else ""
            if d.get("missing"):
                T.append(f"| {f} | {'view absent' if d.get('real_view_absent') else 'undefined'} | | | | | {cs} | {rs} | |")
            elif d.get("knockout_view_absent"):
                c = d["control_knockout_view_absent"]
                T.append(f"| {f} | {d['median_real']:.4g} | ko absent | | {'ko absent' if c['area4'] else 'no'} | "
                         f"{'ko absent' if c['len3'] else 'no'} | {cs} | {rs} | {R} |")
            else:
                ce = d["control_effects"]; z = d["effect"] / d["scale"] if d["scale"] else float("nan")
                fmt = lambda v: "n/a" if v is None else f"{v:.4g}"
                T.append(f"| {f} | {d['median_real']:.4g} | {d['median_ko']:.4g} | {d['effect']:.4g} [{z:.1f}] | "
                         f"{fmt(ce.get('area4'))} | {fmt(ce.get('len3'))} | {cs} | {rs} | {R} |")
        T.append("")
    open(f"{V.OUT}/VERIFY_TABLE.md", "w").write("\n".join(T) + "\n")

    classes = ["ARTIFACT", "DECOMPOSES", "COMBINATION EFFECT", "SINGLE MECHANISM", "NOT RUNNABLE"]
    batches = ["census", "review", "qd"]
    S = ["# Simulation verification: summary", "", f"{len(res)} phenomena. Stage 3 candidates: "
         f"{sum(1 for r in res if r.get('stage3'))}.", "",
         "| class | " + " | ".join(batches) + " | all |", "|---|" + "---|" * (len(batches) + 1)]
    for c in classes:
        n = [sum(1 for r in res if r["batch"] == b and r["class"] == c) for b in batches]
        S.append(f"| {c} | " + " | ".join(map(str, n)) + f" | {sum(n)} |")
    S += ["", "## Flagged (verdict text notes a mismatch or artifact) against class", "",
          "| class | flagged | not flagged |", "|---|---|---|"]
    for c in classes:
        S.append(f"| {c} | {sum(1 for r in res if r['class'] == c and r.get('flagged'))} | "
                 f"{sum(1 for r in res if r['class'] == c and not r.get('flagged'))} |")
    S += ["", "## Literature verdict against class", "", "| verdict | " + " | ".join(classes[:4]) + " |",
          "|---|" + "---|" * 4]
    for v in sorted({str(r.get("verdict")) for r in res}):
        S.append(f"| {v} | " + " | ".join(str(sum(1 for r in res if str(r.get('verdict')) == v and r['class'] == c))
                                          for c in classes[:4]) + " |")
    nf = sum(len(r.get("features") or {}) for r in res); nr = sum(len(r.get("robust", [])) for r in res)
    v2 = sum(1 for r in res for f in r.get("robust", []) if r["features"][f].get("knockout_view_absent"))
    S += ["", f"Tested features: {nf}; robust: {nr} (of which {v2} by deviation V2).",
          f"Card replication (reported only): {rep_ok} of {rep_n} tested features with a card value replicate it.", "",
          "| class | replicate / with card value |", "|---|---|"]
    for c in classes[:4]:
        S.append(f"| {c} | {rep_by_class[c]} / {rep_tot_class[c]} |")
    st = [r for r in res if r.get("stage3")]
    S += ["", "## Stage 3 candidates by batch and class", "", "| batch | COMBINATION EFFECT | SINGLE MECHANISM (flagged) |",
          "|---|---|---|"]
    for b in batches:
        S.append(f"| {b} | {sum(1 for r in st if r['batch'] == b and r['class'] == 'COMBINATION EFFECT')} | "
                 f"{sum(1 for r in st if r['batch'] == b and r['class'] == 'SINGLE MECHANISM')} |")
    open(f"{V.OUT}/VERIFY_REPORT.md", "w").write("\n".join(S) + "\n"); print("\n".join(S))


def packet():
    """Stage 3 (b) review packet for the author: every Stage 3 phenomenon with its card, verdict, citations and notes,
    closely related first, then documented non-trivial, then documented trivial. Writes stage3/REVIEW_PACKET.md."""
    S = os.path.join(V.OUT, "stage3")
    res = {r["pid"]: r for r in json.load(open(f"{V.OUT}/verify.json"))}
    rows = [json.load(open(f"{S}/verdicts/{f}")) for f in sorted(os.listdir(f"{S}/verdicts"))]
    order = lambda r: ({"CLOSELY RELATED": 0, "UNDOCUMENTED CANDIDATE": 0}.get(r.get("provisional_label"), 1),
                       bool(r.get("trivial_or_artifact")), r["pid"])
    L = ["# Stage 3 (b): author review packet", "",
         "Model labels are provisional (VERIFICATION_PLAN.md, Stage 3). For each phenomenon: check that the cited sources "
         "exist and describe the behaviour as stated, then record a final label (documented / closely related / "
         "undocumented candidate). Order: closely related first, then documented and judged non-trivial, then documented "
         "and judged trivial or artifactual.", "",
         f"{len(rows)} phenomena: " + ", ".join(f"{n} {k.lower()}" for k, n in
                                               collections.Counter(r.get('provisional_label') for r in rows).most_common()) + ".", ""]
    for r in sorted(rows, key=order):
        c = res[r["pid"]]
        card = open(f"{S}/cards/{r['pid'].replace(':', '__')}.txt").read()
        rule = card.split("\n\n")[0].replace("RULE (plain words):\n", "")
        L += [f"## {r['pid']}: {r.get('provisional_label')} ({r.get('verdict')})", "",
              f"- Class: {c['class']}; flagged: {c.get('flagged')}; robust features: {', '.join(c.get('robust', []))}",
              f"- Rule: {rule}",
              f"- Standard name given: {r.get('standard_name')}",
              f"- Trivial or artifact: {r.get('trivial_or_artifact')}; confidence: {r.get('confidence')}; "
              f"overturns earlier verdict: {r.get('overturns_first_pass')}; web searches: {r.get('web_searches_used')}",
              f"- Twist (if any): {r.get('undocumented_phenomenology') or ''}",
              f"- Notes: {r.get('notes')}", "- Citations:"] + [f"  - {x}" for x in (r.get("citations") or [])] + \
             ["- Author review: [ ] sources checked  [ ] behaviour as described  Final label: ________", ""]
    open(f"{S}/REVIEW_PACKET.md", "w").write("\n".join(L) + "\n"); print(f"{S}/REVIEW_PACKET.md", len(rows), "entries")


if __name__ == "__main__":
    import sys
    packet() if sys.argv[1:] == ["packet"] else main()
