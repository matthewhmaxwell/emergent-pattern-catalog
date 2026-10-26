"""Run the literature-gate controls (POST-REGISTRATION, EXPLORATORY). See census/controls/SCORING.md.

  python -m census.litcheck_controls primed     # 40 cards, frozen instructions, first-pass model
  python -m census.litcheck_controls neutral    # 20 fabricated cards, priming sentence removed
  python -m census.litcheck_controls second     # strict second pass on primed KNOWN-VARIANT verdicts
  python -m census.litcheck_controls strictall  # strict second pass on EVERY primed verdict (diagnostic, added later)
Resumable: existing verdict files are skipped. Env: ANTHROPIC_API_KEY.
"""
import json, os, re, sys, time
import anthropic
import census.litcheck_runner as LR
import census.litcheck_secondpass as SP

HERE = os.path.dirname(os.path.abspath(__file__))
CARDS = os.path.join(HERE, "controls", "controls.json")
OUT = os.path.join(HERE, "runs", "litcheck_controls")
PRIMING = "Most such behaviours are known; \"known\" is the expected answer and a perfectly good one. "


def run(mode):
    cards = json.load(open(CARDS)); client = anthropic.Anthropic()
    if mode == "neutral":
        assert PRIMING in LR.INSTR, "priming sentence not found verbatim"
        LR.INSTR = LR.INSTR.replace(PRIMING, "")
        cards = [c for c in cards if c["class_id"].startswith("ctrl:F")]
    out = os.path.join(OUT, mode); os.makedirs(out, exist_ok=True)
    for i, c in enumerate(cards, 1):
        fp = LR.verdict_path(out, c)
        if os.path.exists(fp): continue
        t = time.time()
        try: r = LR.run_one(client, "claude-sonnet-5-5", c)
        except Exception as e: r = {"class_id": c["class_id"], "verdict": "ERROR", "error": f"{type(e).__name__}: {e}"}
        r["instructions"] = mode
        json.dump(r, open(fp, "w"), indent=1)
        print(f"[{mode} {i}/{len(cards)}] {c['class_id']} -> {r.get('verdict')} ({time.time()-t:.0f}s) "
              f"{str(r.get('standard_name',''))[:70]}", flush=True)
        time.sleep(2)


def second(all_cards=False):
    """Strict second pass. Default: primed KNOWN-VARIANT verdicts only (as in the census). all_cards=True: every primed
    verdict (diagnostic arm added after seeing the primed results; see SCORING.md addendum)."""
    src = os.path.join(OUT, "primed"); out = os.path.join(OUT, "strictall" if all_cards else "secondpass")
    os.makedirs(out, exist_ok=True)
    cards = {c["class_id"]: c for c in json.load(open(CARDS))}; client = anthropic.Anthropic()
    for f in sorted(os.listdir(src)):
        v = json.load(open(os.path.join(src, f)))
        if not all_cards and v.get("verdict") != "KNOWN-VARIANT": continue
        fp = os.path.join(out, f)
        if os.path.exists(fp): continue
        c = cards[v["class_id"]]
        case = {"class_id": v["class_id"], "rule": c["program_plain"],
                "observed": [f"{o['feature']} = {o['value']} ({LR.GLOSS.get(o['feature'], '')})" for o in c["observed"]],
                "first_pass_name": v.get("standard_name"), "first_pass_family": v.get("model_family"),
                "first_pass_note": v.get("notes")}
        t = time.time()
        try:
            msg = client.messages.create(model="claude-opus-5-5", max_tokens=4500,
                tools=[{"type": "web_search_20250305", "name": "web_search", "max_uses": 12}],
                messages=[{"role": "user", "content": SP.INSTR + "\n\n---\n\n" + SP.card(case)}])
            text = "".join(b.text for b in msg.content if getattr(b, "type", "") == "text")
            m = re.search(r"```json\s*(\{.*?\})\s*```", text, re.S)
            r = {"class_id": v["class_id"], "model": "claude-opus-5-5",
                 **(json.loads(m.group(1)) if m else {"verdict": "PARSE_FAILED"}), "full_text": text}
        except Exception as e:
            r = {"class_id": v["class_id"], "verdict": "ERROR", "error": f"{type(e).__name__}: {e}"}
        json.dump(r, open(fp, "w"), indent=1)
        print(f"[second] {v['class_id']} -> {r.get('verdict')} ({time.time()-t:.0f}s)", flush=True)
        time.sleep(2)


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "second": second()
    elif m == "strictall": second(all_cards=True)
    else: run(m)
    print("CONTROLS", m, "DONE", time.strftime("%H:%M:%S"), flush=True)
