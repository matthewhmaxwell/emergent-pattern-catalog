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

## D2 — P1 audit criterion B3 was mis-specified as two-sided (2026-09-27, AFTER the audit ran — POST-HOC, disclosed)
- **Pre-committed criterion** (`p1_audit.py` docstring, commit `8741a4c`, committed before the audit ran):
  `max_k |B3_k − 0.5| ≤ 0.10`, where B3 forces both agents to the same symbol k.
- **Result as specified:** 3 of 6 models FAIL B3 (seeds 1, 3, 4; max |dev| 0.122–0.156) → formal verdict
  "ALL MODELS PASS = False". This is reported as the official audit outcome.
- **Why the criterion was wrong:** forcing identical symbols gives the agents identical information, so independent
  identical agents can reach at most P(differ) = 2p(1−p) ≤ 0.5 (a fair coin). Only a LEAK (extra asymmetric
  information) can push success ABOVE 0.5; success BELOW 0.5 means a biased tie-breaking coin — suboptimal, not a leak.
  The leak hypothesis therefore predicts a ONE-sided test.
- **Post-hoc one-sided analysis (labelled post-hoc):** the highest B3 across all 48 (model × symbol) cells is 0.527,
  below the 0.5 + 2·SE bound of 0.532 (SE 0.016, n = 1000 per cell). No cell indicates a leak. Every failing cell lies
  BELOW 0.5 (minimum 0.344) — biased tie-breaking on some symbols, consistent with S | tied = 0.46–0.49.
- **All other tests pass on all 6 models:** symbol entropy 2.99–3.00 bits (uniform); S | distinct 0.997–0.999;
  antisymmetry 0.998–0.999; greedy (non-random) talk 0.49–0.50; mute 0.48–0.50; static exchangeability and
  reward-leak audits pass. Replication seeds 4–6: S_full 0.937/0.937/0.938, mute and no-channel ≈ 0.50.
- **Conclusion:** mechanism = a jointly-controlled lottery (the channel carries endogenous randomness); no leak detected.
  The formal pre-committed audit verdict ("not all pass") stands in the record next to this correction.
