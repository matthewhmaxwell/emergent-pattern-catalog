"""EPC census literature-check runner (durable VPS job).

Runs the frozen LITCHECK_PROTOCOL on each distinct unnamed phenomenon using the Anthropic API with its server-side
web_search tool (no browser needed). Independent per call: the model gets only a self-contained class card plus the
protocol's verbatim instructions -- no EPC context, no hoped-for answer. Resumable: a phenomenon whose verdict file
already exists is skipped. One verdict JSON per phenomenon under census/runs/litcheck/verdicts/.

  python -m census.litcheck_runner --reps census/runs/litcheck/representatives.json [--limit N] [--only substr] [--model ...]

Env: ANTHROPIC_API_KEY (from ~/.config/epc-litcheck.env). Run in tmux, e.g.
  set -a; . ~/.config/epc-litcheck.env; set +a; nohup .../python -m census.litcheck_runner ... &
"""
import argparse, json, os, re, time, glob
import anthropic

GLOSS = {
 "g_largest_domain":"one cell-type fills most of the lattice","g_n_domains":"many separate domains",
 "g_moran":"like-valued cells cluster spatially","g_moran_trend":"spatial clustering grows over time",
 "g_corrlen":"long spatial correlation length","g_interface":"long/rough domain boundaries",
 "g_activity":"many cells change per step","g_activity_trend":"activity rising","g_activity_cv":"bursty activity",
 "g_cyclic_changes":"cells cycle through states above chance","g_reversals":"cells flip back and forth",
 "g_type_entropy":"cell types evenly mixed","g_n_types":"number of surviving types","g_spec_peak":"a dominant temporal frequency",
 "p_r":"global phase order (full sync)","p_local_r":"local/partial phase order","p_r_std":"phase order fluctuates",
 "p_vortex_density":"phase vortices/spirals","p_freq_spread":"spread of oscillator natural frequencies",
 "a_polar":"agents align direction (polar order)","a_local_align":"local alignment","a_nn_ratio":"clumping (below 1 = clustered)",
 "a_type_segregation":"agents sort by type","a_speed":"mean speed","a_milling":"rotational/milling motion",
 "n_consensus_gain":"nodes converge to agreement","n_type_modularity":"network splits into like-type modules",
 "n_edge_turnover":"edges rewire","n_cyclic_changes":"node states cycle","n_oscillation":"node states oscillate over time",
 "n_components":"number of connected components","n_giant":"fraction in the giant component","n_type_entropy":"node types evenly mixed",
 "n_degree_cv":"spread of node degrees","f_moran":"field forms spatial structure","f_peak_sharp":"sharp field pattern",
}
INSTR = (
 "You are checking whether a behaviour produced by a very simple simulated rule system is already described in the "
 "scientific literature. You have a web_search tool -- use it, and prefer sources you can actually read. Most such "
 "behaviours are known; \"known\" is the expected answer and a perfectly good one. Work independently and skeptically.\n\n"
 "DO THIS:\n"
 "1. Identify the model family this rule system belongs to (statistical physics, cellular automata, active matter, "
 "opinion dynamics, evolutionary game theory, synchronization, reaction-diffusion, network science, ...).\n"
 "2. Search for the SPECIFIC behaviour described (not just the family): is this behaviour, in this or an equivalent "
 "system, reported? Give the earliest clear source and one review if one exists.\n"
 "3. Verify every citation (title, authors, year, venue) against what the search returns. Do not cite from memory alone.\n"
 "4. Decide exactly one verdict: KNOWN (described for this or an equivalent system) / KNOWN-VARIANT (a close relative is "
 "described -- say precisely what differs and whether it likely matters) / LITERATURE-SILENT (after a genuine search, "
 "nothing describes it -- list the searches you ran).\n"
 "5. Say how confident you are (low/medium/high) and why.\n\n"
 "End your reply with a fenced json block, exactly:\n"
 "```json\n{\"verdict\": \"KNOWN|KNOWN-VARIANT|LITERATURE-SILENT\", \"model_family\": \"...\", \"standard_name\": \"...\", "
 "\"confidence\": \"low|medium|high\", \"citations\": [\"author year, venue - what it shows\"], \"searches\": [\"...\"], "
 "\"notes\": \"any caveat\"}\n```"
)


def card_text(c):
    obs = "\n".join(f"  - {o['feature']} = {o['value']}  ({GLOSS.get(o['feature'],'')})" for o in c.get("observed",[]))
    variants = "\n".join(f"  - {v}" for v in c.get("also_in_class_plain",[]) if v) or "  (none distinct)"
    return (f"THE RULE SYSTEM (plain words):\n{c['program_plain']}\n\n"
            f"OTHER RULE SYSTEMS IN THE SAME BEHAVIOUR CLASS (parameter variants):\n{variants}\n\n"
            f"WHAT WAS OBSERVED IN THE SIMULATION (standout quantitative features; higher = stronger):\n{obs}\n\n"
            f"This behaviour appears across {c.get('size','?')} distinct rule programs of this structure.")


def has_plain(c):
    return c.get("program_plain") and not c["program_plain"].strip().startswith(("C ", "N ", "A ", "F ", "CA ", "CF ", "AF ", "CAF "))


def verdict_path(outdir, c):
    return os.path.join(outdir, c["class_id"].replace(":", "__") + ".json")


def run_one(client, model, c):
    msg = client.messages.create(
        model=model, max_tokens=4000,
        tools=[{"type": "web_search_20250305", "name": "web_search", "max_uses": 10}],
        messages=[{"role": "user", "content": INSTR + "\n\n---\n\n" + card_text(c)}],
    )
    text = "".join(b.text for b in msg.content if getattr(b, "type", "") == "text")
    searches = sum(1 for b in msg.content if getattr(b, "type", "") == "server_tool_use")
    m = re.search(r"```json\s*(\{.*?\})\s*```", text, re.S)
    parsed = json.loads(m.group(1)) if m else {"verdict": "PARSE_FAILED"}
    return {"class_id": c["class_id"], "stage": c["stage"], "family": c["family"],
            "program_terse": c.get("program_terse"), "program_plain": c.get("program_plain"),
            "web_searches_used": searches, "model": model, **parsed, "full_text": text,
            "checked": time.strftime("%Y-%m-%d %H:%M:%S")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", default="census/runs/litcheck/representatives.json")
    ap.add_argument("--out", default="census/runs/litcheck/verdicts")
    ap.add_argument("--model", default="claude-sonnet-5-5")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    reps = json.load(open(a.reps))
    done = {os.path.basename(f)[:-5] for f in glob.glob(a.out + "/*.json")}
    client = anthropic.Anthropic()
    todo = []
    for c in sorted(reps, key=lambda c: c["shortest_len"]):
        key = c["class_id"].replace(":", "__")
        if key in done: continue
        if a.only and a.only not in c["class_id"]: continue
        if not has_plain(c):   # G2-only cryptic twin of a G1 phenomenon -> record as covered, do not spend a call
            json.dump({"class_id": c["class_id"], "stage": c["stage"], "family": c["family"],
                       "program_terse": c.get("program_terse"), "verdict": "COVERED-BY-TWIN",
                       "notes": "Same behaviour as its G1-grammar twin (grammar-only difference); see the G1 verdict for this family/mechanism.",
                       "checked": time.strftime("%Y-%m-%d %H:%M:%S")},
                      open(verdict_path(a.out, c), "w"), indent=1)
            continue
        todo.append(c)
    if a.limit: todo = todo[:a.limit]
    print(f"{len(reps)} phenomena; {len(done)} already have verdicts; running {len(todo)} now (model {a.model})", flush=True)
    for i, c in enumerate(todo, 1):
        t = time.time()
        try:
            r = run_one(client, a.model, c)
        except Exception as e:
            r = {"class_id": c["class_id"], "stage": c["stage"], "verdict": "ERROR", "error": f"{type(e).__name__}: {e}",
                 "checked": time.strftime("%Y-%m-%d %H:%M:%S")}
        json.dump(r, open(verdict_path(a.out, c), "w"), indent=1)
        print(f"[{i}/{len(todo)}] {c['class_id']:28s} -> {r.get('verdict'):16s} "
              f"({r.get('web_searches_used','?')} searches, {time.time()-t:.0f}s) {str(r.get('standard_name',''))[:60]}", flush=True)
        time.sleep(2)
    print("LITCHECK DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    main()
