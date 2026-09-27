# Deviations log — EPC Ring-3 pre-registration v1 (osf.io/up7m9; PREREG_v1.md SHA-256 47c770b1…d47aaa75)

Every departure from, or clarification of, the frozen PREREG_v1.md is logged here with its date and whether it was
made before or after any data existed. Entries are append-only.

## D1 — Evaluation action sampling (2026-09-27, BEFORE any data; before environment code was run)
- **Frozen text (§2 Evaluation):** "greedy for movement and sampled for signals (as in `dial2_ppo.py`)".
- **Fact:** `dial2_ppo.py` at the frozen commit `ac27a8d` evaluates with `greedy=False` (`_hits` → `rollout(...,
  greedy=False)`), i.e. it **samples both** movement and signal actions. The parenthetical is an inaccurate
  description of the trainer the protocol anchors to.
- **Why it matters:** greedy action selection is degenerate for exchangeable agents with identical observations
  (P1): identical greedy actions always collide, making the registered chance level (0.50) unreachable by
  construction. The project's own prior finding (#13 fair alternation) is that symmetric agents can only break
  symmetry stochastically.
- **Resolution:** follow the frozen trainer — evaluation **samples both** movement and signal actions for all
  tasks. No other part of the protocol is affected.
