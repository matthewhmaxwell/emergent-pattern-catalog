"""Strict adversarial SECOND-PASS literature check on the gray-zone KNOWN-VARIANT phenomena.

Reads secondpass.json (class_id, rule, observed, first_pass_*), asks a harder question with a stronger model and
web search: is the observed BEHAVIOUR documented, or is it genuinely literature-silent? Rule-combination novelty
does NOT count -- only behaviour novelty. Writes census/runs/litcheck/secondpass/<id>.json.

  ANTHROPIC_API_KEY=... python litcheck_secondpass.py [--model claude-opus-5-5]
"""
import argparse, json, os, re, time
import anthropic

INSTR = (
 "You are doing a STRICT, adversarial second-pass literature check on a behaviour produced by a very simple "
 "simulated rule system. A first automated pass judged it a VARIANT of a known model but could not find the exact "
 "system in the literature. DO NOT assume the first pass is right -- be willing to overturn it in either direction.\n\n"
 "WHAT MATTERS IS THE BEHAVIOUR, NOT THE RULE TUPLE. Many different rule combinations produce the same known "
 "behaviour, so 'this exact combination of rules is not in any paper' does NOT make it a discovery. Judge whether "
 "the observed emergent BEHAVIOUR / phenomenology is documented.\n\n"
 "Use web_search; open and verify sources. Then decide exactly one verdict:\n"
 "- KNOWN: the observed behaviour is documented for this or an equivalent system (cite).\n"
 "- KNOWN-VARIANT: the behaviour is essentially a documented one, but with a specific twist in the PHENOMENOLOGY "
 "(not merely the rules) that you could not find described -- state the twist and whether it plausibly matters.\n"
 "- LITERATURE-SILENT: after a genuine search, the observed BEHAVIOUR itself appears undocumented. Only use this if "
 "the phenomenology (not the rule combination) is what is missing; state precisely what is unaccounted for and the "
 "searches you ran.\n\n"
 "Also judge whether the observed features look like a GENUINE emergent pattern or like a TRIVIAL/ARTIFACT outcome "
 "(e.g. a random colouring, an inert coupled field, a numerically diverged field, a frozen absorbing state).\n\n"
 "End with a fenced json block exactly:\n"
 "```json\n{\"verdict\": \"KNOWN|KNOWN-VARIANT|LITERATURE-SILENT\", \"behaviour_documented\": true, "
 "\"trivial_or_artifact\": false, \"standard_name\": \"...\", \"citations\": [\"...\"], \"undocumented_phenomenology\": "
 "\"(only if literature-silent)\", \"overturns_first_pass\": false, \"confidence\": \"low|medium|high\", \"notes\": \"...\"}\n```"
)


def card(c):
    return (f"RULE (plain words):\n{c['rule']}\n\n"
            f"OBSERVED (standout quantitative features; higher = stronger):\n  " + ", ".join(c["observed"]) + "\n\n"
            f"FIRST PASS concluded: {c['first_pass_name']}\n"
            f"First-pass family: {c.get('first_pass_family')}\n"
            f"First-pass note (what it could not find): {c.get('first_pass_note')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default="/private/tmp/claude-501/-Users-matthewmaxwell-Downloads/ed5d6e8e-67aa-43e5-b2b6-30638ce75e8d/scratchpad/secondpass.json")
    ap.add_argument("--out", default="census/runs/litcheck/secondpass")
    ap.add_argument("--model", default="claude-opus-5-5")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cases = json.load(open(a.src))
    client = anthropic.Anthropic()
    print(f"second pass: {len(cases)} gray-zone cases, model {a.model}", flush=True)
    for i, c in enumerate(cases, 1):
        fp = os.path.join(a.out, c["class_id"].replace(":", "__") + ".json")
        if os.path.exists(fp):
            print(f"[{i}] {c['class_id']} already done, skip", flush=True); continue
        t = time.time()
        try:
            msg = client.messages.create(model=a.model, max_tokens=4500,
                tools=[{"type": "web_search_20250305", "name": "web_search", "max_uses": 12}],
                messages=[{"role": "user", "content": INSTR + "\n\n---\n\n" + card(c)}])
            text = "".join(b.text for b in msg.content if getattr(b, "type", "") == "text")
            searches = sum(1 for b in msg.content if getattr(b, "type", "") == "server_tool_use")
            m = re.search(r"```json\s*(\{.*?\})\s*```", text, re.S)
            parsed = json.loads(m.group(1)) if m else {"verdict": "PARSE_FAILED"}
            r = {"class_id": c["class_id"], "first_pass_name": c["first_pass_name"], "web_searches_used": searches,
                 "model": a.model, **parsed, "full_text": text, "checked": time.strftime("%Y-%m-%d %H:%M:%S")}
        except Exception as e:
            r = {"class_id": c["class_id"], "verdict": "ERROR", "error": f"{type(e).__name__}: {e}"}
        json.dump(r, open(fp, "w"), indent=1)
        print(f"[{i}/{len(cases)}] {c['class_id']:22s} -> {r.get('verdict'):16s} "
              f"(behaviour_documented={r.get('behaviour_documented')}, trivial={r.get('trivial_or_artifact')}, "
              f"overturns={r.get('overturns_first_pass')}, {r.get('web_searches_used','?')} searches, {time.time()-t:.0f}s)", flush=True)
        time.sleep(2)
    print("SECONDPASS DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    main()
