"""Planted positive controls for the literature gate (PLANTED_CONTROLS_PLAN.md; post-registration, exploratory).

Replays the frozen census gate (census.litcheck_runner first pass, census.litcheck_secondpass strict pass, same models,
prompts, token limits and search caps) on the planted cards, adding only `blocked_domains` to the web-search tool in
the blocked arms. Every reply is stored with the URLs its searches returned, so leaks of the target paper can be found.

  python -m census.litcheck_planted N_blocked_r1   # also N_blocked_r2, N_blocked_r3 (census protocol, blocked)
  python -m census.litcheck_planted N_open         # census protocol, no blocking
  python -m census.litcheck_planted N_neutral      # first pass without the priming sentence, blocked
  python -m census.litcheck_planted N_strictall    # strict prompt on every card (first-pass info from N_blocked_r1)
  python -m census.litcheck_planted D_open         # documented obscure set, census protocol, no blocking
Resumable: existing reply files are skipped. Env: ANTHROPIC_API_KEY.
"""
import json, os, re, sys, time
import anthropic
import census.litcheck_runner as LR
import census.litcheck_secondpass as SP

HERE = os.path.dirname(os.path.abspath(__file__))
PL = os.path.join(HERE, "controls", "planted")
OUT = os.path.join(HERE, "runs", "litcheck_planted")
PRIMING = "Most such behaviours are known; \"known\" is the expected answer and a perfectly good one. "
BASE_BLOCK = ["arxiv.org", "alphaxiv.org", "semanticscholar.org", "researchgate.net", "scholar.google.com",
              "paperswithcode.com", "huggingface.co", "scirate.com", "x.com", "twitter.com", "bsky.app", "reddit.com",
              "ycombinator.com"]
FIRST_MODEL, STRICT_MODEL = "claude-sonnet-5-5", "claude-opus-5-5"


def load(setname):
    cards = json.load(open(os.path.join(PL, f"{setname}_cards.json")))
    truth = json.load(open(os.path.join(PL, f"{setname}_truth.json")))
    for c in cards:                                   # card glosses join the frozen gloss table for rendering only
        for o in c["observed"]: LR.GLOSS.setdefault(o["feature"], o.get("gloss", ""))
    return cards, truth


def _host(d):
    d = re.sub(r"^[a-z]+://", "", d.strip().lower()).split("/")[0]
    return d[4:] if d.startswith("www.") else d


def blocked_for(t):
    extra = [_host(d) for d in [t.get("journal_domain")] + list(t.get("author_domains") or []) if d]
    out = []
    for d in BASE_BLOCK + extra:
        if d and "." in d and d not in out: out.append(d)
    return out[:64]


def _search_tool(max_uses, blocked):
    tool = {"type": "web_search_20250305", "name": "web_search", "max_uses": max_uses}
    if blocked: tool["blocked_domains"] = blocked
    return tool


def _record(msg):
    text = "".join(b.text for b in msg.content if getattr(b, "type", "") == "text")
    searches = [getattr(b, "input", {}) for b in msg.content if getattr(b, "type", "") == "server_tool_use"]
    urls = []
    for b in msg.content:
        if getattr(b, "type", "") == "web_search_tool_result" and isinstance(getattr(b, "content", None), list):
            urls += [getattr(r, "url", "") for r in b.content]
    m = re.search(r"```json\s*(\{.*?\})\s*```", text, re.S)
    try: parsed = json.loads(m.group(1)) if m else {"verdict": "PARSE_FAILED"}
    except json.JSONDecodeError: parsed = {"verdict": "PARSE_FAILED"}
    return {**parsed, "stop_reason": msg.stop_reason, "web_searches_used": len(searches), "search_inputs": searches,
            "result_urls": urls, "full_text": text}


def first_pass(client, c, blocked, instr):
    # identical to census.litcheck_runner.run_one except for the optional blocked_domains on the search tool
    msg = client.messages.create(model=FIRST_MODEL, max_tokens=4000, tools=[_search_tool(10, blocked)],
                                 messages=[{"role": "user", "content": instr + "\n\n---\n\n" + LR.card_text(c)}])
    return {"class_id": c["class_id"], "model": FIRST_MODEL, **_record(msg), "checked": time.strftime("%Y-%m-%d %H:%M:%S")}


def strict_pass(client, c, v, blocked):
    # identical to census.litcheck_controls.second except for the optional blocked_domains on the search tool
    case = {"class_id": c["class_id"], "rule": c["program_plain"],
            "observed": [f"{o['feature']} = {o['value']} ({LR.GLOSS.get(o['feature'], '')})" for o in c["observed"]],
            "first_pass_name": v.get("standard_name"), "first_pass_family": v.get("model_family"),
            "first_pass_note": v.get("notes")}
    msg = client.messages.create(model=STRICT_MODEL, max_tokens=4500, tools=[_search_tool(12, blocked)],
                                 messages=[{"role": "user", "content": SP.INSTR + "\n\n---\n\n" + SP.card(case)}])
    return {"class_id": c["class_id"], "model": STRICT_MODEL, **_record(msg), "checked": time.strftime("%Y-%m-%d %H:%M:%S")}


def _save(path, fn):
    if os.path.exists(path): return json.load(open(path))
    t = time.time()
    try: r = fn()
    except Exception as e: r = {"verdict": "ERROR", "error": f"{type(e).__name__}: {e}"}
    json.dump(r, open(path, "w"), indent=1)
    print(f"  {os.path.relpath(path, OUT)} -> {r.get('verdict')} ({r.get('web_searches_used', '?')} searches, "
          f"{time.time() - t:.0f}s)", flush=True)
    time.sleep(2)
    return r


def arm(name):
    setname = name.split("_")[0]
    cards, truth = load(setname)
    blocked = "blocked" in name or name in ("N_neutral", "N_strictall")
    instr = LR.INSTR
    if name == "N_neutral":
        assert PRIMING in instr, "priming sentence not found verbatim"
        instr = instr.replace(PRIMING, "")
    client = anthropic.Anthropic()
    out = os.path.join(OUT, name); os.makedirs(out, exist_ok=True)
    for i, c in enumerate(cards, 1):
        key = c["class_id"].replace(":", "__"); b = blocked_for(truth[c["class_id"]]) if blocked else None
        print(f"[{name} {i}/{len(cards)}] {c['class_id']}", flush=True)
        if name == "N_strictall":
            src = os.path.join(OUT, "N_blocked_r1", f"{key}.first.json")
            if not os.path.exists(src): print("  needs N_blocked_r1 first"); continue
            _save(os.path.join(out, f"{key}.strict.json"), lambda: strict_pass(client, c, json.load(open(src)), b))
            continue
        v = _save(os.path.join(out, f"{key}.first.json"), lambda: first_pass(client, c, b, instr))
        if name != "N_neutral" and v.get("verdict") == "KNOWN-VARIANT":     # census protocol: strict pass on KV
            _save(os.path.join(out, f"{key}.strict.json"), lambda: strict_pass(client, c, v, b))
    print("PLANTED", name, "DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    for a in sys.argv[1:]: arm(a)
