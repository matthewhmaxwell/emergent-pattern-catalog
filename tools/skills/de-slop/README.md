# de-slop skill snapshot

A copy of Matt's personal `de-slop` writing skill as revised on 2026-09-02
(v1.1.0). The canonical copy lives in the claude.ai skill profile and is
synced from there; this directory is a durable snapshot committed alongside
the repository that was open when the revision was made. It is deliberately
not under `.claude/skills/`, so it does not auto-load in this project and
cannot shadow or duplicate the synced personal skill.

To install or refresh the personal skill from this snapshot, package it with
skill-creator's `package_skill.py` and save the resulting `de-slop.skill`, or
copy the directory into `~/.claude/skills/de-slop/`.

Contents: `SKILL.md`, `references/` (patterns, living corpus, protect list,
ingestion flow, current register), `scripts/slop_scan.py` (the scanner), and
`evals/` (the five test cases used to benchmark the revision).
