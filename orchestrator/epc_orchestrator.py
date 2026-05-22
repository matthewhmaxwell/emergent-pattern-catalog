#!/usr/bin/env python3
"""
EPC Sprint Orchestrator
=======================

Runs code-led sprints of the Emergent Pattern Catalog autonomously on the VPS.
The orchestrator is the bridge between Claude Chat (which designs briefs and
reviews escalations) and Claude Code (which executes briefs). Both sides
communicate via the git repo as a shared message bus; this script removes the
human-as-courier step for code-led-to-code-led transitions.

Lifecycle of one sprint
-----------------------
1. Orchestrator detects a new `docs/sprint_returns/sprint_NN_return.md` on
   main since the last processed SHA.
2. Parses the return doc for the GO / GO-LIMITED / NO-GO decision and the
   suggested next batch.
3. Branches on the decision:
     GO + code-led-eligible  -> draft brief via API -> dispatch Claude Code
     GO-LIMITED / NO-GO      -> notify the user, pause
     chat-led required       -> notify the user, pause
4. On dispatch, Claude Code consumes the brief, executes the sprint, commits
   and pushes a new return doc. Orchestrator detects it on the next poll
   and the cycle repeats.
5. A `max_chained_sprints` safety cap prevents runaway loops.

Run modes
---------
  --once       Run one polling cycle and exit (suitable for cron).
  --daemon     Long-running loop, polls every config['poll_interval_sec'].
  --dry-run    Draft the next brief but do NOT dispatch to Claude Code.
               Combine with --once for a single inspect-only test.
  --resume     Clear the in-flight lock and resume. Use after a manual fix
               or crash recovery.
  --status     Print current state (last processed SHA, chain depth, locks)
               and exit. Does not modify state.

Configuration
-------------
All environment-specific values live in `config.json` next to this script
(or at the path given by --config). The script ships with a placeholder
config.json that Claude Code is expected to populate during first-time
setup — see SETUP.md.

Required external state
-----------------------
- A working git clone of the repo, with push credentials configured.
- The `claude` CLI on PATH, authenticated against the repo's GitHub account.
- ANTHROPIC_API_KEY in the environment (or path-resolvable via the config).
  This key is used ONLY for brief drafting — Claude Code uses its own auth.

Exit codes
----------
  0  clean exit
  1  config error or unrecoverable runtime error
  2  escalation notified, paused (user action required)
  3  --status reported in-flight lock (sprint currently executing)
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import os
import re
import shutil
import smtplib
import ssl
import subprocess
import sys
import textwrap
import time
import urllib.error
import urllib.request
from email.message import EmailMessage
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = SCRIPT_DIR / "config.json"


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

REQUIRED_KEYS = {
    "repo_path",
    "main_branch",
    "poll_interval_sec",
    "working_dir",
    "state_file",
    "log_file",
    "anthropic_api_key_env",
    "anthropic_model",
    "anthropic_max_tokens",
    "claude_code_cmd",
    "notify_method",
    "max_chained_sprints",
}


def load_config(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(
            f"Config file not found at {path}.\n"
            "Run setup first (see SETUP.md). The default config.json shipped "
            "with this script has TODO placeholders that Claude Code is "
            "expected to fill in for the local VPS environment."
        )
    with path.open() as f:
        cfg = json.load(f)
    missing = REQUIRED_KEYS - cfg.keys()
    if missing:
        raise SystemExit(f"config.json missing required keys: {sorted(missing)}")
    todos = [k for k, v in cfg.items() if isinstance(v, str) and v.startswith("TODO:")]
    if todos:
        raise SystemExit(
            f"config.json still has unfilled TODOs for: {todos}\n"
            "Claude Code: please fill these in based on the VPS environment, "
            "then commit and restart the orchestrator."
        )
    return cfg


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def setup_logging(log_file: str, verbose: bool = False) -> logging.Logger:
    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("epc_orchestrator")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    logger.addHandler(fh)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    return logger


# ---------------------------------------------------------------------------
# State management
# ---------------------------------------------------------------------------

def load_state(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {
            "last_processed_sha": None,
            "last_processed_sprint": None,
            "chain_depth": 0,
            "in_flight": None,
            "last_escalation": None,
        }
    with path.open() as f:
        return json.load(f)


def _rotate_state_backups(path: Path, keep: int = 5) -> None:
    """Rotate state.json.bak.N → state.json.bak.(N+1); keep last N versions."""
    for i in range(keep - 1, 0, -1):
        src = Path(str(path) + f".bak.{i}")
        dst = Path(str(path) + f".bak.{i + 1}")
        if src.exists():
            src.replace(dst)
    if path.exists():
        # current state → .bak.1 (a copy, since save_state will replace `path` next)
        import shutil as _sh
        _sh.copy2(path, Path(str(path) + ".bak.1"))


def save_state(path: Path, state: dict[str, Any]) -> None:
    """Atomic write via tmpfile + rename; rotate 5 backups."""
    path.parent.mkdir(parents=True, exist_ok=True)
    _rotate_state_backups(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w") as f:
        json.dump(state, f, indent=2, sort_keys=True)
    tmp.replace(path)


# ---------------------------------------------------------------------------
# Git operations
# ---------------------------------------------------------------------------

def run_git(repo: str, *args: str, check: bool = True) -> subprocess.CompletedProcess:
    cmd = ["git", "-C", repo, *args]
    return subprocess.run(cmd, capture_output=True, text=True, check=check)


def fetch_main(repo: str, branch: str, logger: logging.Logger) -> str:
    """Fetch latest, fast-forward main, return current HEAD SHA on the branch."""
    run_git(repo, "fetch", "origin", branch)
    run_git(repo, "checkout", branch)
    run_git(repo, "merge", "--ff-only", f"origin/{branch}")
    head = run_git(repo, "rev-parse", "HEAD").stdout.strip()
    logger.debug(f"main HEAD = {head}")
    return head


def find_new_return_doc(repo: str, branch: str, since_sha: str | None,
                        logger: logging.Logger) -> tuple[str, str, int] | None:
    """
    Return (sha, return_doc_path, sprint_number) for the newest sprint return
    doc that landed since `since_sha`, or None if there isn't one. If
    `since_sha` is None, returns the most recent return doc on the branch
    (used for initial state bootstrap).
    """
    if since_sha:
        log_range = f"{since_sha}..origin/{branch}"
    else:
        log_range = f"origin/{branch}"
    proc = run_git(repo, "log", log_range, "--name-only", "--pretty=format:%H",
                   "--", "docs/sprint_returns/")
    lines = [l for l in proc.stdout.splitlines() if l]
    if not lines:
        return None
    # Pair each SHA with the files in that commit; pick the latest return doc.
    current_sha = None
    candidates: list[tuple[str, str, int]] = []
    pattern = re.compile(r"docs/sprint_returns/sprint_(\d+)_return\.md$")
    for line in lines:
        if re.fullmatch(r"[0-9a-f]{40}", line):
            current_sha = line
        else:
            m = pattern.search(line.strip())
            if m and current_sha is not None:
                candidates.append((current_sha, line.strip(), int(m.group(1))))
    if not candidates:
        return None
    # Highest sprint number wins; tie-break on SHA recency (first in log).
    candidates.sort(key=lambda t: -t[2])
    sha, path, sprint = candidates[0]
    logger.info(f"new return doc: sprint {sprint} at {path} (sha {sha[:8]})")
    return sha, path, sprint


def read_file_at_head(repo: str, path: str) -> str:
    proc = run_git(repo, "show", f"HEAD:{path}")
    return proc.stdout


# ---------------------------------------------------------------------------
# Return doc parsing
# ---------------------------------------------------------------------------

DECISION_PATTERN = re.compile(
    r"(?:Sprint\s+\d+\s+decision|decision\s*[:\-]|\*\*Decision[:\-])\s*"
    r"(GO-LIMITED|GO\s*LIMITED|NO-?GO|GO)\b",
    re.IGNORECASE,
)

CHAT_LED_MARKERS = [
    r"\bchat-led\b",
    r"\brequires?\s+chat\b",
    r"\bdesign\s+pass\b",
    r"\bmethodology\s+revision\b",
]


def parse_decision(return_doc_text: str) -> str:
    """Returns one of: 'GO', 'GO-LIMITED', 'NO-GO', or 'UNKNOWN'."""
    m = DECISION_PATTERN.search(return_doc_text)
    if not m:
        return "UNKNOWN"
    raw = m.group(1).upper().replace(" ", "-").replace("NOGO", "NO-GO")
    if raw in ("GO", "GO-LIMITED", "NO-GO"):
        return raw
    return "UNKNOWN"


def detect_chat_led_required(return_doc_text: str) -> bool:
    """True if the return doc indicates the next sprint must be chat-led."""
    lower = return_doc_text.lower()
    # Look in the section discussing the next sprint specifically.
    for marker in CHAT_LED_MARKERS:
        if re.search(marker, lower):
            return True
    return False


# ---------------------------------------------------------------------------
# Brief drafting via Anthropic API
# ---------------------------------------------------------------------------

API_URL = "https://api.anthropic.com/v1/messages"

# ---------------------------------------------------------------------------
# Brief queue (Sprint 36 patch)
# ---------------------------------------------------------------------------
# Briefs pre-written by chat are dropped at <working_dir>/briefs/sprint_NN_brief.md
# and picked up here instead of being auto-drafted. This eliminates
# brief-drafting LLM cost for sprints we've planned in advance, and makes the
# escalation rules cleaner (queued briefs are pre-approved; LLM-drafted briefs
# from an empty queue are subject to chat-led-escalation).

def _queue_dir(cfg: dict[str, Any]) -> Path:
    return Path(cfg["working_dir"]) / "briefs"


def read_queued_brief(cfg: dict[str, Any], sprint_num: int) -> str | None:
    """Return a pre-written brief for sprint_num if present in the queue dir."""
    path = _queue_dir(cfg) / f"sprint_{sprint_num}_brief.md"
    if path.exists():
        return path.read_text()
    return None


def _disk_free_fraction(path: str = "/") -> float:
    """Fraction of free disk space at the given mount."""
    st = os.statvfs(path)
    if st.f_blocks == 0:
        return 1.0
    return st.f_bavail / st.f_blocks


BRIEF_SYSTEM_PROMPT = textwrap.dedent("""\
    You are drafting a code-led sprint brief for the Emergent Pattern Catalog
    (EPC) project. The brief specifies, mechanically, what Claude Code must
    do to execute the next sprint.

    Use the conventions documented in CLAUDE.md: pre-flight checks, scope
    parts (A, B, C...), acceptance criteria, post-flight, tag bump. Mirror
    the structure of the prior return doc's referenced precedent sprints.
    Do not invent new methodology. Do not modify detectors or models — that
    requires a chat-led sprint and an explicit escalation.

    Output the brief as Markdown, ready to commit to the repo. Do not
    include preamble, commentary, or explanation of your choices outside
    the brief body itself.
""")


def draft_brief(return_doc: str, claude_md: str, depth_gap: str,
                api_key: str, model: str, max_tokens: int,
                logger: logging.Logger) -> str:
    user_msg = textwrap.dedent(f"""\
        Below are three documents from the EPC repository.

        === CLAUDE.md (project ground truth) ===
        {claude_md}

        === docs/depth_gap.md (current audit matrix) ===
        {depth_gap}

        === Latest sprint return doc ===
        {return_doc}

        ---

        Draft the next code-led sprint brief. Use the prior return doc's
        "Suggested next batch" recommendation as your scope unless it is
        clearly contradicted by depth_gap.md. Include all required sections:
        pre-flight, scope parts, acceptance criteria, post-flight, tag bump.
    """)
    body = json.dumps({
        "model": model,
        "max_tokens": max_tokens,
        "system": BRIEF_SYSTEM_PROMPT,
        "messages": [{"role": "user", "content": user_msg}],
    }).encode("utf-8")
    req = urllib.request.Request(
        API_URL,
        data=body,
        headers={
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        },
        method="POST",
    )
    logger.info(f"drafting brief via {model} ({max_tokens=})")
    try:
        with urllib.request.urlopen(req, timeout=180) as resp:
            payload = json.loads(resp.read())
    except urllib.error.HTTPError as e:
        err_body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"Anthropic API HTTP {e.code}: {err_body}") from e

    # Concatenate text blocks from the response.
    chunks = [
        block.get("text", "")
        for block in payload.get("content", [])
        if block.get("type") == "text"
    ]
    brief = "".join(chunks).strip()
    if not brief:
        raise RuntimeError(f"API returned no text content: {payload}")
    logger.info(f"brief drafted ({len(brief)} chars)")
    return brief


# ---------------------------------------------------------------------------
# Claude Code dispatch
# ---------------------------------------------------------------------------

DISPATCH_PROMPT_TEMPLATE = textwrap.dedent("""\
    You are the Claude Code thread for sprint {next_sprint}. The brief is
    in your working directory at {brief_path}. Execute it.

    Standing instructions:

    1. Read the brief in full. Execute every part (A, B, C, ...) in the
       order specified.
    2. Run pre-flight checks before touching anything. If they fail, STOP
       and write an aborted return doc to docs/sprint_returns/ explaining
       the failure. Do not attempt repair.
    3. Run the post-flight checks before committing. Same rule on failure.
    4. Commit, tag (v0.{next_sprint}.0 unless the brief specifies
       otherwise), push to origin/{branch}.
    5. Write docs/sprint_returns/sprint_{next_sprint}_return.md. Use the
       structural template from the most recent prior return doc.
    6. Do not modify detectors, models, or the orchestration registry
       unless the brief explicitly scopes that work. Out-of-scope findings
       are carry-forwards, not in-sprint fixes.
    7. When done, exit cleanly. The orchestrator polls for your return doc.
""")


def dispatch_to_claude_code(brief: str, next_sprint: int, branch: str,
                            cfg: dict[str, Any], logger: logging.Logger,
                            chain_depth: int = 0) -> int:
    working_dir = Path(cfg["working_dir"])
    working_dir.mkdir(parents=True, exist_ok=True)
    brief_path = working_dir / f"sprint_{next_sprint}_brief.md"
    brief_path.write_text(brief)
    logger.info(f"brief written to {brief_path}")

    prompt = DISPATCH_PROMPT_TEMPLATE.format(
        next_sprint=next_sprint, brief_path=brief_path, branch=branch,
    )

    cmd_template = cfg["claude_code_cmd"]
    cwd = cfg.get("claude_code_cwd") or cfg["repo_path"]

    # Chain-aware pre-flight skip: if we're already mid-chain (depth > 0), the
    # immediately-preceding sprint's POST-flight pytest just ran and passed at
    # this HEAD. Skipping pre-flight saves ~11 min per chained sprint. The
    # brief template reads this env var to decide whether to skip.
    env = os.environ.copy()
    env["EPC_CHAIN_DEPTH"] = str(chain_depth)
    if chain_depth > 0:
        env["EPC_SKIP_PREFLIGHT_PYTEST"] = "1"

    # Two invocation modes depending on whether the CLI takes the prompt
    # as an argument or via stdin.
    if cfg.get("claude_code_uses_stdin", False):
        cmd = list(cmd_template)
        logger.info(f"dispatching to Claude Code: {' '.join(cmd)} (stdin) "
                    f"chain_depth={chain_depth}")
        proc = subprocess.run(cmd, input=prompt, cwd=cwd, text=True, env=env)
    else:
        cmd = [arg.replace("{prompt}", prompt) for arg in cmd_template]
        logger.info(f"dispatching to Claude Code: {cmd_template[0]} ... "
                    f"chain_depth={chain_depth}")
        proc = subprocess.run(cmd, cwd=cwd, text=True, env=env)

    logger.info(f"Claude Code exited with code {proc.returncode}")
    return proc.returncode


# ---------------------------------------------------------------------------
# Notification
# ---------------------------------------------------------------------------

def notify_user(reason: str, details: str, cfg: dict[str, Any],
                logger: logging.Logger) -> None:
    method = cfg["notify_method"]
    timestamp = dt.datetime.now().isoformat(timespec="seconds")
    body = f"[{timestamp}] EPC orchestrator escalation: {reason}\n\n{details}\n"

    if method == "file":
        path = Path(cfg.get("notify_file") or
                    Path(cfg["working_dir"]) / "NOTIFY_USER.txt")
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a") as f:
            f.write(body + "\n" + ("-" * 60) + "\n")
        logger.warning(f"escalation written to {path}: {reason}")

    elif method == "webhook":
        url = cfg["notify_webhook_url"]
        payload = json.dumps({"text": body}).encode("utf-8")
        req = urllib.request.Request(url, data=payload,
                                     headers={"content-type": "application/json"})
        try:
            urllib.request.urlopen(req, timeout=10)
            logger.warning(f"escalation posted to webhook: {reason}")
        except Exception as e:
            logger.error(f"webhook notify failed: {e}")

    elif method == "email":
        msg = EmailMessage()
        msg["From"] = cfg["notify_email_from"]
        msg["To"] = cfg["notify_email_to"]
        msg["Subject"] = f"EPC orchestrator: {reason}"
        msg.set_content(body)
        host = cfg["notify_email_smtp_host"]
        port = int(cfg.get("notify_email_smtp_port", 587))
        user = cfg["notify_email_smtp_user"]
        pw_env = cfg["notify_email_smtp_password_env"]
        password = os.environ.get(pw_env) if pw_env else None
        try:
            with smtplib.SMTP(host, port) as s:
                s.starttls(context=ssl.create_default_context())
                if user and password:
                    s.login(user, password)
                s.send_message(msg)
            logger.warning(f"escalation emailed to {cfg['notify_email_to']}: {reason}")
        except Exception as e:
            logger.error(f"email notify failed: {e}")

    else:
        logger.error(f"unknown notify_method {method!r}; details: {body}")


# ---------------------------------------------------------------------------
# Main cycle
# ---------------------------------------------------------------------------

def run_cycle(cfg: dict[str, Any], state: dict[str, Any], dry_run: bool,
              logger: logging.Logger) -> bool:
    """One polling cycle. Returns True if a sprint was dispatched (or would
    have been in dry-run), False if nothing to do."""
    repo = cfg["repo_path"]
    branch = cfg["main_branch"]

    head = fetch_main(repo, branch, logger)
    if state["last_processed_sha"] == head:
        logger.debug("no new commits since last processed sha; nothing to do")
        return False

    found = find_new_return_doc(repo, branch, state["last_processed_sha"], logger)
    if not found:
        # No new return doc, but head moved (e.g., docs commit). Advance.
        state["last_processed_sha"] = head
        return False

    return_sha, return_path, sprint_num = found
    if state["last_processed_sprint"] is not None and \
            sprint_num <= state["last_processed_sprint"]:
        logger.info(f"sprint {sprint_num} already processed; skipping")
        state["last_processed_sha"] = head
        return False

    # Read the return doc, CLAUDE.md, depth_gap.md at HEAD.
    return_doc = read_file_at_head(repo, return_path)
    claude_md = read_file_at_head(repo, "CLAUDE.md")
    depth_gap = read_file_at_head(repo, "docs/depth_gap.md")

    decision = parse_decision(return_doc)
    next_sprint = sprint_num + 1

    # Queued briefs are pre-approved by chat — they trump escalate_if_chat_led
    # AND the chat-led-required text detector. They do NOT trump explicit
    # escalate_on decisions (GO-LIMITED / NO-GO still pause for review).
    queued_brief_present = read_queued_brief(cfg, next_sprint) is not None

    chat_led = (
        detect_chat_led_required(return_doc)
        and cfg.get("escalate_if_chat_led", True)
        and not queued_brief_present
    )

    logger.info(f"sprint {sprint_num} decision: {decision} "
                f"(chat-led required: {chat_led}; "
                f"queued brief for sprint {next_sprint}: {queued_brief_present})")

    # Branching rules.
    if decision in cfg.get("escalate_on", ["GO-LIMITED", "NO-GO"]) or chat_led \
            or decision == "UNKNOWN":
        reason = (
            f"sprint {sprint_num} return needs your attention "
            f"(decision={decision}, chat_led={chat_led})"
        )
        notify_user(reason, return_doc[:4000], cfg, logger)
        state["last_processed_sha"] = head
        state["last_processed_sprint"] = sprint_num
        state["last_escalation"] = {
            "sprint": sprint_num, "decision": decision, "chat_led": chat_led,
            "at": dt.datetime.now().isoformat(timespec="seconds"),
        }
        state["chain_depth"] = 0
        return False

    if decision not in cfg.get("auto_dispatch_on", ["GO"]):
        logger.warning(f"decision {decision} not in auto_dispatch_on; escalating")
        notify_user(
            f"unexpected decision {decision} for sprint {sprint_num}",
            return_doc[:4000], cfg, logger,
        )
        state["last_processed_sha"] = head
        state["last_processed_sprint"] = sprint_num
        return False

    if state["chain_depth"] >= cfg["max_chained_sprints"]:
        logger.warning(f"chain depth {state['chain_depth']} hit cap; escalating")
        notify_user(
            f"max chained sprints reached ({cfg['max_chained_sprints']}); "
            f"pausing for review before sprint {sprint_num + 1}",
            "Resume by clearing chain_depth in state.json (or run --resume).",
            cfg, logger,
        )
        state["last_processed_sha"] = head
        state["last_processed_sprint"] = sprint_num
        return False

    # Safety gates before dispatch
    free_frac = _disk_free_fraction("/")
    if free_frac < 0.10:
        notify_user(
            f"disk free below 10% ({free_frac:.1%}); refusing to start new chain link",
            "Free space, then run --resume.", cfg, logger,
        )
        return False

    # (next_sprint computed earlier for queue lookup)

    # Pre-chain rollback tag (only when starting a fresh chain, depth==0).
    if state.get("chain_depth", 0) == 0:
        tag = f"pre-chain-{dt.datetime.now().strftime('%Y%m%d-%H%M%S')}"
        try:
            run_git(repo, "tag", tag, head)
            logger.info(f"pre-chain rollback tag created: {tag} -> {head[:8]}")
            state["last_pre_chain_tag"] = tag
        except subprocess.CalledProcessError as e:
            logger.warning(f"could not create pre-chain tag (non-fatal): {e}")

    # QUEUE PRIORITY: pre-written brief beats LLM draft
    brief = read_queued_brief(cfg, next_sprint)
    brief_source = "queue"
    if brief is None:
        # Empty queue → fall back to LLM-drafted brief. Apply chat-led escalation
        # check ONLY in this path; queued briefs are pre-approved.
        if chat_led:
            reason = (
                f"sprint {sprint_num} return implies chat-led work AND queue "
                f"has no sprint_{next_sprint}_brief.md; pausing for chat input."
            )
            notify_user(reason, return_doc[:4000], cfg, logger)
            state["last_processed_sha"] = head
            state["last_processed_sprint"] = sprint_num
            state["last_escalation"] = {
                "sprint": sprint_num, "decision": decision, "chat_led": True,
                "queue_miss_for": next_sprint,
                "at": dt.datetime.now().isoformat(timespec="seconds"),
            }
            return False
        api_key = os.environ.get(cfg["anthropic_api_key_env"])
        if not api_key:
            raise SystemExit(
                f"environment variable {cfg['anthropic_api_key_env']} is not set; "
                "cannot draft brief"
            )
        brief = draft_brief(
            return_doc=return_doc, claude_md=claude_md, depth_gap=depth_gap,
            api_key=api_key, model=cfg["anthropic_model"],
            max_tokens=cfg["anthropic_max_tokens"], logger=logger,
        )
        brief_source = "llm"

    import hashlib as _h
    brief_sha = _h.sha256(brief.encode()).hexdigest()[:12]
    logger.info(f"brief for sprint {next_sprint}: source={brief_source} "
                f"len={len(brief)} sha256={brief_sha}")
    if dry_run:
        out = Path(cfg["working_dir"]) / f"sprint_{next_sprint}_brief.DRAFT.md"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(brief)
        logger.info(f"[dry-run] brief written to {out}; NOT dispatching")
        state["last_processed_sha"] = head
        state["last_processed_sprint"] = sprint_num
        return True

    state["in_flight"] = {
        "sprint": next_sprint,
        "started": dt.datetime.now().isoformat(timespec="seconds"),
    }
    save_state(Path(cfg["state_file"]), state)

    rc = dispatch_to_claude_code(brief, next_sprint, branch, cfg, logger,
                                 chain_depth=state.get("chain_depth", 0))
    if rc != 0:
        notify_user(
            f"Claude Code exited non-zero ({rc}) on sprint {next_sprint}",
            f"Brief at {Path(cfg['working_dir']) / f'sprint_{next_sprint}_brief.md'}",
            cfg, logger,
        )
        state["in_flight"] = None
        state["last_processed_sha"] = head
        state["last_processed_sprint"] = sprint_num
        return False

    # Code is expected to push a new return doc. Next poll cycle picks it up.
    state["in_flight"] = None
    state["last_processed_sha"] = head
    state["last_processed_sprint"] = sprint_num
    state["chain_depth"] = state["chain_depth"] + 1
    return True


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH),
                        help="path to config.json")
    parser.add_argument("--once", action="store_true",
                        help="run one polling cycle and exit")
    parser.add_argument("--daemon", action="store_true",
                        help="long-running daemon; polls at config interval")
    parser.add_argument("--dry-run", action="store_true",
                        help="draft the next brief but do not dispatch")
    parser.add_argument("--resume", action="store_true",
                        help="clear in-flight lock and reset chain_depth")
    parser.add_argument("--status", action="store_true",
                        help="print state.json and exit")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    cfg = load_config(Path(args.config))
    logger = setup_logging(cfg["log_file"], verbose=args.verbose)
    state_path = Path(cfg["state_file"])
    state = load_state(state_path)

    if args.status:
        print(json.dumps(state, indent=2, sort_keys=True))
        return 3 if state.get("in_flight") else 0

    if args.resume:
        state["in_flight"] = None
        state["chain_depth"] = 0
        save_state(state_path, state)
        logger.info("state reset: in_flight cleared, chain_depth=0")
        return 0

    if not (args.once or args.daemon):
        parser.error("specify one of --once, --daemon, --status, --resume")

    if state.get("in_flight"):
        logger.error(f"in-flight lock present: {state['in_flight']}. "
                     "Use --resume to clear, or wait for the in-flight sprint "
                     "to complete and push its return doc.")
        return 3

    try:
        if args.once:
            run_cycle(cfg, state, dry_run=args.dry_run, logger=logger)
            save_state(state_path, state)
            return 0

        # --daemon
        idle_poll = cfg["poll_interval_sec"]
        active_poll = int(cfg.get("poll_interval_sec_active", 60))
        logger.info(f"daemon starting; idle poll = {idle_poll}s, "
                    f"active-chain poll = {active_poll}s")
        while True:
            try:
                run_cycle(cfg, state, dry_run=args.dry_run, logger=logger)
                save_state(state_path, state)
            except Exception:
                logger.exception("cycle raised; sleeping and continuing")
            # Dynamic poll: if we're mid-chain (depth>0), poll fast so the
            # next sprint fires as soon as the last commit lands. Else idle.
            in_chain = state.get("chain_depth", 0) > 0
            sleep_s = active_poll if in_chain else idle_poll
            time.sleep(sleep_s)

    except KeyboardInterrupt:
        logger.info("interrupted by user")
        return 0


if __name__ == "__main__":
    sys.exit(main())
