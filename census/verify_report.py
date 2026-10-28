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


if __name__ == "__main__":
    main()
