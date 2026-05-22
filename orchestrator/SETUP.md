# Orchestrator Setup — Instructions for Claude Code

You are Claude Code on the VPS. This bundle contains a sprint orchestrator
that will remove the user from the per-sprint handoff loop. Your job: install
it, configure it for this VPS environment, verify it on a dry-run, then start
it under systemd. After that, the user is only involved when something
actually requires their judgment.

The user is not the right person to fill in VPS-specific values — *you* are.
You have the filesystem, you have your own auth setup, you know what the
shell aliases are. Make the calls.

## Step 1 — Place files on the VPS

```bash
sudo mkdir -p /opt/epc-orchestrator
sudo chown $USER:$USER /opt/epc-orchestrator
cp epc_orchestrator.py config.json SETUP.md /opt/epc-orchestrator/
sudo cp epc-orchestrator.service /etc/systemd/system/
```

If `/opt/epc-orchestrator` is wrong for this VPS (e.g., the user prefers
`/home/USER/epc-orchestrator`, or `/opt` isn't writable), pick a path that
makes sense and update the `WorkingDirectory=` and `ExecStart=` lines in
the systemd unit accordingly.

## Step 2 — Fill in config.json

Open `/opt/epc-orchestrator/config.json` and replace every `TODO:` placeholder
with the real value for this VPS. The orchestrator refuses to start if any
TODO remains. You decide:

- **repo_path** — absolute path to the local clone of
  `emergent-pattern-catalog`. You already pull/push from this clone, so
  use the same path. Verify with `git -C <path> rev-parse --show-toplevel`.

- **working_dir** — where briefs and state live. Suggest
  `/var/lib/epc-orchestrator/` if you have write access there, otherwise
  `~/epc-orchestrator-state/`. Must be writable by the service user.

- **state_file** and **log_file** — usually `<working_dir>/state.json` and
  `<working_dir>/orchestrator.log`. Create the directory first.

- **claude_code_cmd** — how to invoke yourself non-interactively. You know
  this; the user does not. Common forms:
    `["claude", "-p", "{prompt}"]`  — current Claude Code CLI, prompt as arg
    `["claude", "--print"]`         — older form, reads from stdin (set
                                       `claude_code_uses_stdin: true`)
  If you use a wrapper or alias, expand it to the real binary path.

- **claude_code_cwd** — usually leave `null` (the orchestrator uses
  `repo_path`). Override only if your install needs a different working dir
  to find its config.

- **notify_method** — start with `"file"` (simplest, writes escalations to
  `notify_file`). If the user has a webhook (Slack/Discord) or wants
  email, you can switch later. For `"file"`, set `notify_file` to a path
  the user will see — e.g., add a line to their `.bashrc` to print the
  file contents on shell start, or symlink it to their home directory.

- **anthropic_api_key_env** — leave as `"ANTHROPIC_API_KEY"` unless the
  user has chosen a different env var name. The key itself goes in
  `/etc/epc-orchestrator/env` (Step 4 below).

## Step 3 — Verify config

```bash
cd /opt/epc-orchestrator
python3 -c "import json; c = json.load(open('config.json')); \
  todos = [k for k,v in c.items() if isinstance(v, str) and v.startswith('TODO:')]; \
  print('TODOs remaining:', todos)"
```

Expect `TODOs remaining: []`. If anything remains, fix it.

## Step 4 — Anthropic API key

The orchestrator needs an `ANTHROPIC_API_KEY` for brief drafting. This is a
**separate key from whatever your own Claude Code auth uses** — it's billed
to the user's API account and used only for the brief-drafting call (not
for the sprint execution, which goes through you).

If the user hasn't provided a key, escalate via `notify_method=file` —
write to the configured notify file with:

> Orchestrator needs an Anthropic API key. Create one at
> https://console.anthropic.com, then run:
>   sudo mkdir -p /etc/epc-orchestrator
>   echo 'ANTHROPIC_API_KEY=sk-ant-...' | sudo tee /etc/epc-orchestrator/env
>   sudo chmod 600 /etc/epc-orchestrator/env
>   sudo chown root:<service-user> /etc/epc-orchestrator/env

Don't proceed past this step without the key. (You can still do the
dry-run in Step 6 without it, but live dispatch needs it.)

## Step 5 — Bootstrap state

The orchestrator needs to know which sprint return doc was the last
"already processed" one so it doesn't try to dispatch Sprint 36 again when
that one lands. The first time you run it, bootstrap state to point at the
most recent return doc currently on `main`:

```bash
cd /opt/epc-orchestrator
python3 -c "
import json, subprocess, pathlib
repo = json.load(open('config.json'))['repo_path']
head = subprocess.check_output(['git', '-C', repo, 'rev-parse', 'HEAD'], text=True).strip()
import re
log = subprocess.check_output(['git', '-C', repo, 'log', '--name-only', '--pretty=format:%H', '--', 'docs/sprint_returns/'], text=True)
shas = re.findall(r'^[0-9a-f]{40}', log, re.MULTILINE)
nums = [int(m.group(1)) for m in re.finditer(r'sprint_(\d+)_return\.md', log)]
last_sprint = max(nums) if nums else None
state = {'last_processed_sha': head, 'last_processed_sprint': last_sprint,
         'chain_depth': 0, 'in_flight': None, 'last_escalation': None}
out = json.load(open('config.json'))['state_file']
pathlib.Path(out).parent.mkdir(parents=True, exist_ok=True)
json.dump(state, open(out, 'w'), indent=2)
print(f'bootstrapped state at {out}: head={head[:8]}, last sprint={last_sprint}')
"
```

This tells the orchestrator: "the world starts here; only watch for sprints
that land after this." If Sprint 36 has already been pushed by the time you
run this, the orchestrator will wait for Sprint 36's return doc and
dispatch Sprint 37 from there. (Recommended — Sprint 36 goes through the
existing manual handoff as a known-good baseline.)

## Step 6 — Dry-run

Verify everything works end-to-end *without* actually dispatching a sprint.
Manually rewind state to the previous sprint so the orchestrator sees a
"new" return doc on its next cycle:

```bash
# Temporarily set last_processed_sprint to N-1 where N is current.
# Then run --once --dry-run, which will draft a brief but not dispatch.
python3 epc_orchestrator.py --once --dry-run --verbose
```

Expect:
- The script reads CLAUDE.md, depth_gap.md, and the latest return doc.
- It calls the Anthropic API once.
- It writes a `sprint_NN_brief.DRAFT.md` to the working directory.
- It does NOT invoke Claude Code.
- It exits 0.

Read the drafted brief. If it looks reasonable (pre-flight, parts, acceptance
criteria, post-flight), the brief-drafting half works. Restore state from
backup, move on.

If the draft is structurally wrong (missing sections, hallucinated
methodology, picks the wrong batch), STOP and notify the user. The
brief-drafting prompt may need to be revised before letting it dispatch
real sprints.

## Step 7 — Live one-shot

Once the dry-run looks good and the user has approved the brief format:

```bash
python3 epc_orchestrator.py --once --verbose
```

This will draft *and dispatch* a brief if there's a pending un-dispatched
return doc with a GO decision. Otherwise it exits cleanly. Watch the log
file. When you (Claude Code) take over from the dispatch, you'll be invoked
with the brief in hand and the standard sprint instructions.

## Step 8 — Daemon

When you're confident, enable the systemd unit:

```bash
sudo systemctl daemon-reload
sudo systemctl enable --now epc-orchestrator
sudo systemctl status epc-orchestrator
journalctl -u epc-orchestrator -f   # tail the logs
```

The daemon now polls every `poll_interval_sec`. It will keep going until:
- A GO-LIMITED or NO-GO decision arrives (escalate, pause)
- A chat-led sprint is required (escalate, pause)
- `max_chained_sprints` is hit (escalate, pause — safety cap)
- The service is stopped

Stop the daemon any time with `sudo systemctl stop epc-orchestrator`.
Resume with `--resume` if you cleared an in-flight lock manually:
`python3 epc_orchestrator.py --resume`.

## Step 9 — Commit the configured orchestrator

Once it's working, commit the orchestrator into the repo so future
revisions ship as commits, not as out-of-band file transfers:

```bash
mkdir -p <repo_path>/orchestrator
cp /opt/epc-orchestrator/epc_orchestrator.py <repo_path>/orchestrator/
cp /opt/epc-orchestrator/config.json <repo_path>/orchestrator/config.example.json
# Strip secrets from the example config — leave structure, replace values
# with TODO placeholders for the public commit.
cp /etc/systemd/system/epc-orchestrator.service <repo_path>/orchestrator/
cp SETUP.md <repo_path>/orchestrator/
cd <repo_path>
git add orchestrator/
git commit -m "Add sprint orchestrator (Sprint 36 follow-up infrastructure)"
git push origin main
```

**Important**: do NOT commit the real `config.json` (it may contain VPS
paths or notification credentials). Commit `config.example.json` with
placeholder values. Add `/opt/epc-orchestrator/config.json` (the real one)
to a private location or a `.gitignore` entry covering the live path.

## What's running after this

- `epc-orchestrator.service` polls the repo every 5 min (configurable).
- When Sprint 36's return doc lands on main, the orchestrator sees it,
  parses the decision (GO per the Sprint 36 brief's expectations), drafts
  a Sprint 37 brief via the Anthropic API, writes it to working_dir, and
  invokes Claude Code (you) on it.
- You execute Sprint 37, commit, push, write the return doc, exit.
- Orchestrator detects Sprint 37's return doc on its next poll, drafts
  Sprint 38's brief, dispatches, and so on — up to `max_chained_sprints`.
- Anything other than a clean GO triggers an escalation and pauses the
  chain. The user reads the escalation, acts (or asks Claude Chat to act),
  and resumes the daemon.

## Limits and known gaps

- **The orchestrator assumes you can be invoked non-interactively** with
  the brief as a prompt and that you'll run to completion. If your install
  needs a different invocation pattern (an SDK call rather than CLI, a
  wrapper script, MCP, ...), edit `claude_code_cmd` accordingly.
- **The brief-drafting LLM call is the riskiest piece.** A bad brief
  produces a bad sprint. The dry-run mode (Step 6) is the safety net; use
  it any time you suspect the prompt has drifted. The brief-drafting
  system prompt is at the top of `epc_orchestrator.py` —
  `BRIEF_SYSTEM_PROMPT`. Adjust if necessary.
- **No multi-sprint planning.** The orchestrator drafts one sprint at a
  time, based on the most recent return doc. It does not look ahead.
  Multi-sprint design (e.g., "the next 4 sprints close the lattice_2d
  batch") still requires Claude Chat.
- **Escalation == pause.** The orchestrator never tries to be clever
  about GO-LIMITED or NO-GO outcomes. Those are always handed back to the
  user. Don't try to teach it to recover automatically — the whole point
  is that those decisions are humans-required.

## When something goes wrong

- Check the log: `journalctl -u epc-orchestrator -n 200`
  or `tail -200 /var/lib/epc-orchestrator/orchestrator.log`
- Check current state: `python3 epc_orchestrator.py --status`
- Stuck "in-flight" lock from a crash: `python3 epc_orchestrator.py --resume`
- Brief drafting looks wrong on a dry-run: edit `BRIEF_SYSTEM_PROMPT` in
  `epc_orchestrator.py`, re-run dry-run, commit when satisfied.
- Anything else: escalate to the user with the log excerpt.
