"""Stage 3 of VERIFICATION_PLAN.md: literature check of the phenomena that Stage 2 sends on (written before any class
was read; the frozen verify.py is not changed).

For each phenomenon with stage3 = true in census/runs/verify/verify.json, the strict literature prompt
(census.litcheck_secondpass.INSTR, unchanged) is given the program, the robust features with their values and the
simulation evidence (10 seeds, 4x area, 3x length, knock-out, sub-programs). Model verdicts map to the plan's labels:
KNOWN -> DOCUMENTED, KNOWN-VARIANT -> CLOSELY RELATED, LITERATURE-SILENT -> UNDOCUMENTED CANDIDATE. These are
provisional until the author reviews the evidence and the cited sources (plan, Stage 3 (b)).

  python -m census.verify_stage3 prep        # build the cards (no model call); writes stage3/cards/*.txt
  python -m census.verify_stage3 ask         # one strict call per card; resumable. Env: ANTHROPIC_API_KEY
  python -m census.verify_stage3 report      # stage3/STAGE3.md
"""
import json, os, re, sys, time
import numpy as np
from census import verify as V
from census import litcheck_runner as LR
from census import litcheck_secondpass as SP

OUT = os.path.join(V.OUT, "stage3")
MODEL = "claude-opus-5-5"
LABEL = {"KNOWN": "DOCUMENTED", "KNOWN-VARIANT": "CLOSELY RELATED", "LITERATURE-SILENT": "UNDOCUMENTED CANDIDATE"}


def _plain(it):
    p = V.grammar(it["grammar"]).decode(it["bits"])
    if it["grammar"] == "g1":
        from census import plain
        return p, plain.describe(p)
    return p, V.grammar("g2").describe(p) + "  (compact grammar notation)"


def _med(fps, f, view):
    a = V._col(fps, f, view)
    return None if a is None else round(float(np.median(a)), 4)


def card(it, cl, rec):
    p, rule = _plain(it)
    view = it["view"]; cen = rec["conditions"]["census"]
    lines = []
    for f in cl["robust"]:
        d = cl["features"][f]
        g = LR.GLOSS.get(f)
        lines.append(f"  - {f}{f' ({g})' if g else ''}: median {d['median_real']} with interaction, "
                     f"{d['median_ko']} with interaction knocked out (difference {d['effect']}, robust spread {d['scale']}); "
                     f"difference at 4x area {d['control_effects'].get('area4')}, at 3x run length {d['control_effects'].get('len3')}")
    subs = []
    for sp in rec.get("subprograms", []):
        if "fps" not in sp: continue
        kept = ", ".join(str(i + 1) for i in sp["rules"])
        vals = "; ".join(f"{f} {_med(sp['fps'], f, view)}" for f in cl["robust"])
        subs.append(f"  - rules {kept} only: {vals}")
    ev = ("SIMULATION EVIDENCE (added for this check). Each feature below was re-measured on 10 fresh random seeds, "
          "averaged over the late part of the run. The interaction knock-out is the same program with every rule that "
          "reads other entities removed (one-way rules such as spontaneous flips, decay and emission are kept; phase "
          "coupling is set to zero). Each listed feature differs from its knock-out by at least 3 robust spreads at the "
          "original size, and keeps its sign and at least half its size at 4x area and at 3x run length.\n"
          + "\n".join(lines))
    if subs:
        ev += ("\n\nSUB-PROGRAMS (same program with some rules removed, original size, 10 seeds; median of each feature "
               f"above; the full program has {len(p.rules)} rules):\n" + "\n".join(subs))
        if cl["class"] == "COMBINATION EFFECT":
            ev += "\nNo sub-program came within 3 robust spreads of the full program on every feature above."
    return (f"RULE (plain words):\n{rule}\n\n"
            f"OBSERVED (features that passed the simulation checks):\n  " + ", ".join(cl["robust"]) + "\n\n"
            f"FIRST PASS concluded: {it['verdict']}\n"
            f"First-pass text (excerpt): {(it.get('text') or '').strip()[:1500]}\n\n" + ev)


def prep():
    inv = {it["pid"]: it for it in json.load(open(f"{V.OUT}/inventory.json"))}
    res = [r for r in json.load(open(f"{V.OUT}/verify.json")) if r.get("stage3")]
    os.makedirs(f"{OUT}/cards", exist_ok=True)
    for r in res:
        it = inv[r["pid"]]; rec = json.load(open(f"{V.OUT}/sims/{r['pid'].replace(':', '__')}.json"))
        txt = card(it, r, rec)
        open(f"{OUT}/cards/{r['pid'].replace(':', '__')}.txt", "w").write(txt)
    print(len(res), "cards written to", f"{OUT}/cards")


def ask():
    import anthropic
    client = anthropic.Anthropic(); os.makedirs(f"{OUT}/verdicts", exist_ok=True)
    cards = sorted(os.listdir(f"{OUT}/cards"))
    for i, fn in enumerate(cards, 1):
        fp = f"{OUT}/verdicts/{fn[:-4]}.json"
        if os.path.exists(fp): continue
        t = time.time(); txt = open(f"{OUT}/cards/{fn}").read()
        try:
            msg = client.messages.create(model=MODEL, max_tokens=4500,
                tools=[{"type": "web_search_20250305", "name": "web_search", "max_uses": 12}],
                messages=[{"role": "user", "content": SP.INSTR + "\n\n---\n\n" + txt}])
            text = "".join(b.text for b in msg.content if getattr(b, "type", "") == "text")
            searches = sum(1 for b in msg.content if getattr(b, "type", "") == "server_tool_use")
            m = re.search(r"```json\s*(\{.*?\})\s*```", text, re.S)
            parsed = json.loads(m.group(1)) if m else {"verdict": "PARSE_FAILED"}
            r = {"pid": fn[:-4].replace("__", ":"), "model": MODEL, "web_searches_used": searches, **parsed,
                 "provisional_label": LABEL.get(parsed.get("verdict"), "UNRESOLVED"), "full_text": text,
                 "checked": time.strftime("%Y-%m-%d %H:%M:%S")}
        except Exception as e:
            r = {"pid": fn[:-4].replace("__", ":"), "verdict": "ERROR", "error": f"{type(e).__name__}: {e}"}
        json.dump(r, open(fp, "w"), indent=1)
        print(f"[{i}/{len(cards)}] {r['pid']} -> {r.get('verdict')} / {r.get('provisional_label')} "
              f"({r.get('web_searches_used', '?')} searches, {time.time() - t:.0f}s)", flush=True)
        time.sleep(2)
    print("STAGE3 DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


def report():
    rows = [json.load(open(f"{OUT}/verdicts/{f}")) for f in sorted(os.listdir(f"{OUT}/verdicts"))]
    L = ["# Stage 3 literature check (provisional until author review)", "",
         "| phenomenon | model verdict | provisional label | standard name | trivial or artifact | confidence |",
         "|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['pid']} | {r.get('verdict')} | {r.get('provisional_label')} | {str(r.get('standard_name', ''))[:80]} | "
                 f"{r.get('trivial_or_artifact')} | {r.get('confidence')} |")
    open(f"{OUT}/STAGE3.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    {"prep": prep, "ask": ask, "report": report}[sys.argv[1]]()
