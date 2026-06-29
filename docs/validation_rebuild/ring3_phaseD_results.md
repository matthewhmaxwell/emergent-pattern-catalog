# Ring 3 Phase D — the richer world (multi-agent): results

Validated rationale (after retracting the "low odds" over-claim): the single-agent toy ceiling is
*local* and does not transfer; multi-agent worlds add competency CLASSES that cannot exist with one
agent, so this is the highest-expected-value direction, not a long shot. Same loop as before
(open-ended generator → capable learner → discriminating baseline/ablation → agent-observer → literature
gate). Honest bound unchanged: literature gate before any novelty claim.

## D1 — substrate validation via emergent communication (`commworld_ppo.py`)
A referential signaling game in a gridworld: 3 colored goals; one is the hidden TARGET. The SPEAKER
sees the target color and emits a symbol (K=3); the LISTENER sees the symbol + directions to each
colored goal and navigates; shared reward for reaching the target-color goal. No built-in symbol→color
map — speaker and listener must CO-ADAPT a convention (Lewis signaling). Joint PPO trains both MLP
policies. A single agent cannot do this.

**Result (800 iters):** success 0.05 → **0.63** (plateau ~0.6–0.66). **Channel-scramble debunk** (feed
the listener a random symbol): success → **0.31 = chance (1/3)**. The collapse under scramble proves the
symbol genuinely carried the target information — **genuine emergent communication**, not a shortcut.
Imperfect (0.63, i.e. the co-adapted convention is partial/noisy), but real.

**Classification:** emergent communication / Lewis signaling — **known to science** (large emergent-comms
literature) → **0 new-to-science**. The significance is that it is the FIRST genuinely multi-agent
competency the pipeline has reached: it empirically confirms the richer world produces qualitatively new
competency classes (communication) the single-agent regime structurally could not, validating the
substrate + the channel-ablation discriminator for multi-agent work.

## Next — D2: the open-ended multi-agent hunt
Generalize from the hand-built referential game to a GENERATOR over multi-agent task structures
(coordination on a shared goal; anti-coordination / role-division; conditional/division-of-labor;
adversarial/deception variants; with/without a channel), train joint PPO, and have the agent-observer
NAME whatever multi-agent competency emerges (debunk via channel-ablation / role-swap / partner-swap),
classify vs the map or flag OFF-MAP, literature-gate any survivor. This is where the qualitatively
larger competency space actually gets searched for something not nameable in the existing literature.
