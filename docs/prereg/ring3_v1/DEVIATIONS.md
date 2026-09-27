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

## D3 — P3 conflict-test interpretation label (2026-09-27, AFTER the test ran — POST-HOC reading, disclosed)
- **Pre-committed** (`p3_conflict.py`, commit `6f1e910`, committed before the script touched any trained model):
  outcome classes for full rho=1 models and a three-row interpretation table.
- **Result as specified.** Controls valid (memoryless rho=1 follow the landmark 0.92; full rho=0 follow the cue 0.93).
  Full rho=1 models, original seeds 1–3 and replication seeds 4–6 alike: **PERCEPTION-FIRST** (conflict landmark share
  L 0.858–0.875, all 6 seeds) and memory content **NECESSARY** (nu_resample 1.13–1.23, all 6 seeds). The pre-committed
  interpretation for this combination is "integrates both cues; no artifact attribution". **This is the official outcome.**
  Replication (protocol step 1): seeds 4–6 reproduce the registered pattern (S 0.922–0.925, memwipe 0.52–0.61;
  memoryless fair baseline 0.914–0.926).
- **Why the label is imprecise.** The table anticipated (A) "follows landmark + memory dispensable" and (M) "follows cue
  + memory necessary". The observed combination dissociates goal CHOICE from policy EXECUTION:
  - choice is driven by the landmark: cue share under conflict 0.13 vs 0.08 for memoryless agents, which see the cue
    only at t=0 — the stored cue adds ~5 points of weight;
  - the recurrent state is needed to reach any goal: under resample, timeouts rise from 0.00 to 0.31–0.39 and wrong-goal
    arrivals from 0.07 to 0.21–0.25, although a memoryless-trained policy scores 0.92.
- **Post-hoc reading (labelled).** Ablation necessity (memwipe, and the marginally in-distribution resample) measures the
  policy's dependence on its recurrent pathway, not its use of the stored cue. Neither (A) nor (M) holds as stated.
  The formal P3 MISS stands: C2 operationalized as ablation necessity fails at rho=1; C2 read as "the cheapest
  sufficient information source drives the decision" is supported by the conflict test. No novelty claim: the gap
  between a lesion's effect and the information a pathway carries is a known distinction.
- **Exploratory (not scored).** Across rho the decision source switches abruptly: median L 0.86 at rho=1.0 but
  0.19 / 0.12 / 0.09 / 0.07 at rho=0.75 / 0.5 / 0.25 / 0.0 — once the landmark is sometimes absent, agents follow the
  always-available stored cue even when a valid landmark disagrees (with elevated conflict timeouts at rho=0.75).
