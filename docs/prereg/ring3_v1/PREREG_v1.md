# EPC Ring-3 Pre-registration v1
## Adversarial, pre-registered tests of the communication-niche claim and the minimal-sufficient-mechanism law

**Author:** Matt Maxwell (Emergent Pattern Catalog, EPC). **Status:** FROZEN 2026-09-27, before any data were collected.
At freeze this document is hashed (SHA-256 in `PREREG_v1.sha256`), committed to the EPC repository, and
registered on OSF. Nothing below may change after freeze; changes go in a dated deviations log.

---

### 0. Status at registration
- **No data exist for any task below.** The environments for P1–P4 have not been written. The only code that
  exists is the trainer family they derive from (`analysis/ring3_competency/dial2_ppo.py`, built on
  `openhunt4_ppo.py`), frozen at repository commit `ac27a8dba88ca511059e177f52d7b0bbcdd7edef`
  (github.com/matthewhmaxwell/emergent-pattern-catalog; `dial2_ppo.py` SHA-256
  `95cb793e09d38ec38f63598edf40995d38b3a9dd9691a45421a528e3638b523b`).
- **Code-before-data rule:** environment + analysis code (`prereg_envs.py`, `prereg_analysis.py`) will be committed
  and SHA-256 hashed *before the first training run*; the hash is appended to the OSF record as an update.
- **Smoke-test policy:** one ≤50-iteration run per environment is permitted solely to verify the code executes
  and to measure the random-policy chance level. Smoke-test task performance is not analysed or reported as data.
- **No hyperparameter tuning.** All hyperparameters are inherited verbatim from `dial2_ppo.py` at the frozen
  commit (recurrent GRU policy, hidden size 96, batch 192, PPO settings as coded, 1500 training iterations).

### 1. Claims under test (stated before this registration, in dated repository documents)
- **C1 — Communication niche** (`docs/validation_rebuild/ring3_coordination_cascade.md`; `dial_results.md`):
  *"Symmetric coordination never forces communication. … It is forced only by information ASYMMETRY."*
  Emergent coordination uses the cheapest available medium in the order
  `environmental focal point < observation < coordinate focal point < communication`.
- **C2 — Minimal sufficient mechanism** (`docs/validation_rebuild/ring3_synthesis.md` §5):
  *"A learner realizes the cheapest mechanism the task's structure permits."*

Every prediction below is the prediction **C1/C2 make as stated**. Competing hypotheses and the investigator's
prior credence that the stated claim survives are recorded so that the tests are genuinely risky.

### 2. Common methods
- **Learner:** parameter-shared recurrent PPO (the `dial2_ppo.py` trainer). Single-agent tasks use the same
  trainer with one agent. **Seeds:** 3 independent training seeds per condition (1, 2, 3).
- **Evaluation:** ≥2000 evaluation episodes per trained model (fixed evaluation seeds, disjoint from training),
  greedy for movement and sampled for signals (as in `dial2_ppo.py`).
- **Ablations (evaluation-time):** *mute* = partner's signal replaced by a uniformly random symbol;
  *blind* = partner-position observation replaced by a uniformly random position; *memwipe* = recurrent state reset
  to zero at every step.
- **Fair baselines (training-time):** for every claim that a channel is *unused*, an identical model is trained from
  scratch **without** that channel (constant input). This is the project's own correction for eval-time ablation
  artifacts (`ring3_synthesis.md` §2).
- **Necessity** of a channel, per seed: `ν = (S_full − S_ablated) / (S_full − S_chance)`, where `S` is the task
  success metric and `S_chance` is the random-policy level measured at smoke test.
  **Collapse:** ν ≥ 0.50. **Dead:** |ν| ≤ 0.15 **and** the fair baseline reaches ≥ 0.90·S_full. **Partial:** otherwise.
- **Ceiling:** `S_ceiling` = score of a scripted oracle policy for each task (P1: "higher symbol takes L";
  P2: speaker walks to the target, listener walks to the speaker; P3: go to the cued goal; P4: both walk to a fixed
  cell), measured at smoke test before any training.
- **Learned criterion:** a seed is "learned" if `S_full − S_chance ≥ 0.20·(S_ceiling − S_chance)`. Unlearned seeds
  are replaced by the next seed (max 3 replacements, all logged). Exception: P1, where failure to learn is itself
  an outcome (see P1).
- **Decision rule (per prediction):** **HIT** if the predicted fingerprint holds in ≥2 of 3 seeds *and* at the seed
  median. **MISS** if the competing fingerprint holds in ≥2 of 3 seeds *and* at the seed median.
  **INCONCLUSIVE** otherwise. All four outcomes are reported regardless of direction.

---

### 3. Predictions

#### P1 (headline) — Symmetric anti-coordination: does a talk channel get recruited to break ties?
*Tests C1 at its most vulnerable point: identical agents, identical information, but the task needs them to differ.*
- **Environment (non-spatial, one-shot):** two exchangeable agents (shared parameters, **no agent-ID input**),
  identical constant observation at t=0. **Talk step:** each emits a symbol from K=8 (sampled from its policy).
  **Act step:** each observes its own and its partner's symbol, then chooses L or R. **Reward:** 1 to both if their
  choices differ, else 0.
- **Levels:** chance (no usable channel, best symmetric mixed strategy) = 0.50. Ceiling with one talk round
  and K=8 = 0.9375 (distinct symbols → deterministic split; equal symbols → coin flip).
- **Positive control:** a hand-coded "higher symbol takes L" protocol is run in the environment and must reach
  0.9375 ± 0.02; this confirms the task is solvable and the channel is capable of carrying the solution.
- **C1 prediction (symmetric environmental information ⇒ communication not recruited):**
  `S_full ≤ 0.60` and the fair no-channel baseline within 0.05 of `S_full`. (ν is undefined when `S_full` ≈ chance,
  so this fingerprint is defined directly by `S_full` and the fair baseline rather than by ν.)
- **Competing hypothesis (jointly-controlled lottery / cheap-talk symmetry breaking; Aumann, Maschler & Stearns
  1968):** `S_full ≥ 0.70`, mute **collapses** (ν ≥ 0.5), fair no-channel baseline ≤ 0.60.
- **If C1's outcome occurs:** reported with the positive control as *"the channel could solve the task but PPO did
  not discover it within budget"* — a learnability result (cf. cascade rung 6), **not** evidence that
  communication is useless under symmetry.
- **Investigator credence that C1 survives: 0.35.**
- **If MISS:** C1 is refined to *"communication is forced by information asymmetry **or** by a symmetry-breaking
  demand among exchangeable agents, where the channel carries endogenously generated randomness."*

#### P2 — Asymmetric information with observation available: is observation cheaper than communication?
*Tests the ordering `observation < communication` where the cascade never tested it (rung 3 removed observation).*
- **Environment:** the `dial2_ppo.py` world at α = 1 (only the speaker sees the target cell), modified so that
  **each agent also observes its partner's position**. N=7, T=24, K=8, hits per episode with respawn, as coded.
- **C1 prediction:** the listener follows the speaker: blind **collapses** (ν_blind ≥ 0.5) and mute **dead**.
- **Competing hypothesis (efficiency under time pressure):** signalling beats following, so mute is load-bearing
  (ν_mute ≥ 0.30) regardless of blind.
- **Investigator credence that C1 survives: 0.50.**

#### P3 — Memory–perception dial (single agent): does the learner use the cheapest sufficient mechanism?
*Tests C2 outside coordination, as a dose–response curve analogous to the §8 communication dial.*
- **Environment:** one agent, N=7 grid, two candidate goal cells visible every step. At t=0 only, a cue identifies
  the correct goal. With probability ρ (drawn once per episode) the correct goal also carries a visible landmark
  at every step; otherwise the two goals are indistinguishable after t=0. Reward +1 correct goal, −0.3 wrong
  goal, −0.01 per step; episode ends on reaching a goal or at T=24.
- **Conditions:** ρ ∈ {1.0, 0.75, 0.5, 0.25, 0.0}, 3 seeds each. Fair memoryless (feed-forward) baselines at
  ρ = 1.0 and ρ = 0.0. Chance (random choice between goals) = 0.50.
- **C2 prediction:** memwipe **dead** at ρ = 1.0 (perception suffices); memwipe **collapses** at ρ = 0.0; and the
  seed-median ν_memwipe is non-decreasing in (1 − ρ), with Spearman correlation ≥ 0.8 across the 5 levels.
- **Competing hypothesis:** the reliable t=0 cue is always sufficient, so the learner uses memory at every ρ
  (ν_memwipe ≥ 0.30 at ρ = 1.0) — "cheapest" is not perception-first when a stored cue is equally sufficient.
- **Investigator credence that C2 survives: 0.55.**

#### P4 — Cross-play of learned conventions: were the cascade's focal points just shared training?
*Addresses the obvious reviewer critique of C1: rung-2 "coordinate focal points" might be conventions that only
exist because both agents were trained together.*
- **Environment:** cascade rung-2 world — two agents, **blind** (self-position only), shared absolute frame,
  K=8 channel, reward when both occupy the same cell, random respawn after each meeting; N=7, T=24.
- **Evaluation:** self-play (both agents from the same run) and cross-play (agents from different runs, all 6
  ordered pairs of the 3 runs).
- **C1 prediction:** mute **dead** in self-play (replicates rung 2) **and** the channel does not rescue cross-play:
  mean cross-play success of channel-trained pairs ≤ mean cross-play success of no-channel-trained pairs + 0.05
  (symbol meanings are arbitrary across independent runs; this form is used because cross-play may sit near chance,
  where ν is undefined).
- **Exploratory (not scored):** whether independent runs converge on the same meeting cell (e.g., the grid centre —
  Schelling salience), and the cross-play/self-play success ratio.
- **Investigator credence that C1 survives: 0.80.**

---

### 4. What a MISS triggers (the novelty protocol)
A miss is the only principled place a genuinely new result could appear, so it is handled conservatively:
1. Replicate with 3 additional seeds.
2. Run the EPC debunking battery: observation-leak and reward-leak audits, provably-symmetric controls,
   fair-baseline re-check.
3. Independent literature pass before any novelty claim (P1's competing hypothesis is already known in game
   theory; the open question is whether it *emerges* and whether it falsifies C1 as stated).
4. Report in full regardless of outcome.

### 5. Deviations, exclusions, stopping
- Exclusions: unlearned seeds (Section 2) and crashed runs (rerun with the same seed). No other exclusions.
- Stopping: fixed 1500 iterations per run; no early stopping; no additional iterations after viewing results.
- Any bug fix to environment code after the code-hash commit is logged as a deviation with its diff and rationale.

### 6. Budget
≈ 45 training runs on CPU (VPS, niced), plus evaluation. No paid APIs.
