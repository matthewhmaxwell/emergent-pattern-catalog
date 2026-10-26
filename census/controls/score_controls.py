"""Score the literature-gate controls against census/controls/SCORING.md (post-registration, exploratory).

  python census/controls/score_controls.py     -> census/runs/litcheck_controls/SCORES.md + scores.json

Verdict-level scoring follows SCORING.md exactly. For each fabricated card it also records whether the reply TEXT
flags the mismatch (contradiction words) and, for the strict pass, its behaviour_documented / trivial_or_artifact
fields. These extras are descriptive; they do not change the pre-specified score (KNOWN = false known).
F13 is reported separately as a defective control: the first pass found a documented phenomenon matching its
'fabricated' observation (flocking by turning away), so its truth label is contaminated.
"""
import json, os, re
from scipy.stats import beta

HERE = os.path.dirname(os.path.abspath(__file__)); RUNS = os.path.join(HERE, "..", "runs", "litcheck_controls")
truth = json.load(open(os.path.join(HERE, "truth.json")))
FLAG = re.compile(r"contradict|inconsistent|does not match|doesn't match|not expected|unexpected|implausib|impossible|"
                  r"artifact|anomal|cannot produce|should not|would not produce|at odds|mismatch|opposite", re.I)
DEFECTIVE = {"ctrl:F13"}


def cp(k, n):
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return round(float(lo), 3), round(float(hi), 3)


def load(arm):
    d = os.path.join(RUNS, arm)
    return {json.load(open(os.path.join(d, f)))["class_id"]: json.load(open(os.path.join(d, f)))
            for f in sorted(os.listdir(d)) if f.endswith(".json")} if os.path.isdir(d) else {}


def main():
    arms = {a: load(a) for a in ("primed", "neutral", "secondpass", "strictall")}
    rows, summary = [], {}
    for arm, V in arms.items():
        kn = [c for c in V if truth[c]["kind"] == "known"]; fb = [c for c in V if truth[c]["kind"] == "fabricated"]
        fb_clean = [c for c in fb if c not in DEFECTIVE]
        known_ident = sum(V[c].get("verdict") in ("KNOWN", "KNOWN-VARIANT") for c in kn)
        false_known = sum(V[c].get("verdict") == "KNOWN" for c in fb)
        false_known_clean = sum(V[c].get("verdict") == "KNOWN" for c in fb_clean)
        flagged_text = sum(bool(FLAG.search(V[c].get("full_text", "") or "")) for c in fb if V[c].get("verdict") == "KNOWN")
        s = {"n_known": len(kn), "known_identified": known_ident,
             "known_ci": cp(known_ident, len(kn)) if kn else None,
             "n_fabricated": len(fb), "false_known": false_known, "false_known_ci": cp(false_known, len(fb)) if fb else None,
             "n_fabricated_clean": len(fb_clean), "false_known_clean": false_known_clean,
             "false_known_clean_ci": cp(false_known_clean, len(fb_clean)) if fb_clean else None,
             "false_known_whose_text_flags_mismatch": flagged_text}
        if arm in ("secondpass", "strictall"):
            s["fabricated_behaviour_documented_false"] = sum(V[c].get("behaviour_documented") is False for c in fb)
            s["fabricated_trivial_or_artifact"] = sum(bool(V[c].get("trivial_or_artifact")) for c in fb)
        summary[arm] = s
        for c in sorted(V):
            v = V[c]
            rows.append({"arm": arm, "card": c, "kind": truth[c]["kind"], "intended": truth[c]["intended"],
                         "verdict": v.get("verdict"), "standard_name": v.get("standard_name"),
                         "behaviour_documented": v.get("behaviour_documented"), "trivial": v.get("trivial_or_artifact"),
                         "text_flags_mismatch": bool(FLAG.search(v.get("full_text", "") or ""))})
    json.dump({"summary": summary, "rows": rows}, open(os.path.join(RUNS, "scores.json"), "w"), indent=1)
    L = ["# Literature-gate controls: scores (rules in census/controls/SCORING.md)", ""]
    for arm, s in summary.items():
        L.append(f"- **{arm}**: " + json.dumps(s))
    L += ["", "| arm | card | kind | verdict | documented? | trivial? | text flags mismatch | intended | gate's name |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['arm']} | {r['card']} | {r['kind']} | {r['verdict']} | {r['behaviour_documented']} | {r['trivial']} | "
                 f"{r['text_flags_mismatch']} | {r['intended'][:70]} | {str(r['standard_name'])[:70]} |")
    open(os.path.join(RUNS, "SCORES.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:8]))


if __name__ == "__main__":
    main()
