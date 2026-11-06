"""Reading inputs for the planted-pipeline gate (PLANTED_PIPELINE_PLAN "Scoring"; reading/BRIEF.md).

  python -m census.planted_reading prep    -> reading/inputs/Cxx.md, reading/KEY.json (never shown to the reader),
                                              reading/leaks.json (leak-detector hits, for adjudication by hand)
  python -m census.planted_reading stage3b-prep      -> stage3_review/inputs/Sxx.md, stage3_review/KEY.json
  python -m census.planted_reading stage3b-collect   -> stage3_review/final_labels.json (from the agents' outputs)

Every final KNOWN-VARIANT reply (strict where run, else first pass) of runs r1-r3 that the leak detector does not flag
is put in front of the reader. Card files are named C01.. and replies R01.. in a shuffled order (fixed seed), so the
reader sees neither the card, the run nor the pass.
"""
import json, os, random, sys

from census import planted_pipeline as P
from census.controls.score_planted import leak as leak_fn

D = f"{P.OUT}/reading"
SKIP = ("full_text", "search_inputs", "result_urls", "checked", "class_id", "model", "stop_reason", "web_searches_used")


def prep():
    import glob
    T = P.truth(); items, leaks = [], []
    for fn in sorted(glob.glob(f"{P.OUT}/cards/*.json")):
        c = P.J(fn); key = c["class_id"].replace(":", "__"); t = T[f"planted:{c['card']}"]
        for run, v, src, first in P.final_verdicts(key):
            lk = leak_fn(v, t) or (leak_fn(first, t) if src == "strict" else "")
            if lk: leaks.append({"run": run, "card": c["class_id"], "source": src, "detector": lk, "verdict": v.get("verdict")})
            if v.get("verdict") == "KNOWN-VARIANT" and not lk: items.append((c, run, src, v))
    rng = random.Random(20261104); rng.shuffle(items)
    cards = sorted({c["class_id"] for c, *_ in items}); rng.shuffle(cards)
    cid = {k: f"C{i + 1:02d}" for i, k in enumerate(cards)}
    key, per = {}, {k: [] for k in cards}
    for i, (c, run, src, v) in enumerate(items):
        r = f"R{i + 1:02d}"; key[r] = {"run": run, "card": c["class_id"], "source": src, "file": cid[c["class_id"]]}
        per[c["class_id"]].append((r, v))
    os.makedirs(f"{D}/inputs", exist_ok=True)
    for k in cards:
        c = P.J(f"{P.OUT}/cards/{k.replace(':', '__')}.json")
        L = [f"# {cid[k]}", "", "**The card the gate saw:**", "", "```", c["text"], "```", ""]
        for r, v in sorted(per[k]):
            L += [f"## Reply {r}", "", "**JSON fields:** " + json.dumps({f: v[f] for f in v if f not in SKIP}, ensure_ascii=False), "",
                  "**Reply text:**", "", v.get("full_text", ""), ""]
        open(f"{D}/inputs/{cid[k]}.md", "w").write("\n".join(L))
    P.dump(key, f"{D}/KEY.json"); P.dump(leaks, f"{D}/leaks.json")
    print(len(items), "replies to read in", len(cards), "files;", len(leaks), "leak-detector hits")
    for x in leaks: print("  leak?", x)


def stage3b_prep():
    """stage3_review/inputs/Sxx.md (Stage 3 (a) card + reply, neutral id) and KEY.json (id -> class_id, cutoff date)."""
    import glob
    T = P.truth(); S = f"{P.OUT}/stage3_review"; os.makedirs(f"{S}/inputs", exist_ok=True); key = {}
    fs = sorted(glob.glob(f"{P.OUT}/stage3/verdicts/*.json")); random.Random(20261105).shuffle(fs)
    for i, fn in enumerate(fs):
        k = os.path.basename(fn)[:-5]; cid = k.replace("__", ":"); card = cid.split(":")[1]; sid = f"S{i + 1:02d}"
        v = P.J(fn); txt = open(f"{P.OUT}/stage3/cards/{k}.txt").read()
        key[sid] = {"class_id": cid, "cutoff": T[f"planted:{card}"]["v1_date"]}
        L = [f"# {sid}", "", "## Card (rule in plain words + simulation evidence)", "", "```", txt, "```", "", "## Stage 3 (a) reply", "",
             "**JSON fields:** " + json.dumps({f: v[f] for f in v if f not in SKIP + ("provisional_label",)}, ensure_ascii=False), "",
             "**Full text:**", "", v.get("full_text", ""), ""]
        open(f"{S}/inputs/{sid}.md", "w").write("\n".join(L))
    P.dump(key, f"{S}/KEY.json"); print(len(key), "stage 3 (b) inputs;", {s: x["cutoff"] for s, x in key.items()})


def stage3b_collect():
    S = f"{P.OUT}/stage3_review"; key = P.J(f"{S}/KEY.json", {}); out = {}
    for sid, x in key.items():
        o = P.J(f"{S}/outputs/{sid}.json")
        o = o[0] if isinstance(o, list) and o else o
        out[x["class_id"]] = {"final_label": (o or {}).get("final_label"), "agent_output": f"outputs/{sid}.json"}
    P.dump(out, f"{S}/final_labels.json"); print(out)


if __name__ == "__main__":
    cmds = {"prep": prep, "stage3b-prep": stage3b_prep, "stage3b-collect": stage3b_collect}
    cmds[sys.argv[1]]() if sys.argv[1:] and sys.argv[1] in cmds else sys.exit(__doc__)
