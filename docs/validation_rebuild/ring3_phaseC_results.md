# Ring 3 Phase C — enriched grammar (meta/conditional rules): results

**Code:** `analysis/ring3_competency/metaworld_ppo.py` (+ `metaworld_verify.py`). Ran durably under
`vps-job phasec-meta2`. **Bound unchanged: 0 competencies new to science.** The contribution is a
further **map extension** — a META-LEVEL family of competencies (inferring/selecting/tracking a latent
rule from feedback) the capable learner reaches but the object-level grammar never demanded.

## The probe
Past the count/order/collect grammar, three rule-types whose competency is a runtime/meta capability;
the observation adds a **reward-feedback channel** (last collection good/bad) that makes inference
possible. GRU-PPO vs a true memoryless MLP control; MCC by good-collection precision (chance = 1/3).

## Catch first (8th over-claim avoided): a reward-design artifact
The first run came back all-null with PPO precision 0.00 and netR ≈ −0.9 — because a −1 penalty for bad
collections made **non-participation optimal** (collect nothing, eat only step costs). That tests
nothing about meta-learning. Retuned the bad penalty to −0.3 (keeping the full ±1 feedback *signal*);
agents then participated (netR +0.6 → +3..+7) and the rules became learnable.

## Results (600 iters, retuned)
| rule | what it demands | PPO prec | memoryless | gap | curve |
|---|---|---|---|---|---|
| `cued` | conditional rule-selection (read start-cue, hold, act) | 0.85 | 0.33 | +0.52 | 0.34→0.85 |
| `infer` | in-episode rule inference (infer hidden good-type from feedback) | 0.71–0.73 | 0.38 | +0.34 | 0.34→0.49→0.67→0.71 |
| `switch` | continual track-and-adapt (good-type flips after each success) | 0.81 | 0.47 | +0.34 | 0.34→…→0.81 |

All three demanding (PPO ≥ 0.6, gap ≥ 0.2), memoryless control near chance.

## Interventional verification of `infer` (the headline) — `metaworld_verify.py`
- **Feedback ablation → 0.33 (exactly chance).** Zeroing the reward-feedback channel collapses the agent
  to chance: it GENUINELY uses feedback to identify the hidden rule. Not a shortcut — this is real
  in-episode rule inference (meta-learning / RL²-style).
- **Mid-episode rule flip → post-flip 0.15 (below chance).** When the good-type is resampled at the
  midpoint, the agent keeps collecting the old type: it does **one-shot** inference (infer once, commit,
  exploit) and does NOT continually re-infer. The flip debunk bounds the claim precisely and blocks a
  "continual meta-learner" over-statement.
- Clean dissociation: `switch` (rule changes every success) learned *continual* tracking (0.81), while
  `infer` (rule fixed) learned the cheaper *one-shot* inference — the learner fits the competency the
  environment actually demands.

## What this adds to the governing-dynamics map
Earlier the map was object-level: state-modes {commitment (navigation), storage (memory),
accumulation (counting)} + sequencing, vs reactive forms. Phase C reaches a **meta-level family**: the
state holds a *belief about a latent rule*, updated from feedback, and drives behavior —
- one-shot **rule inference** (`infer`, verified genuine + verified one-shot),
- continual **rule tracking** (`switch`),
- cue-conditioned **rule selection** (`cued`).
Whether this is best modeled as a new state-MODE (belief-over-latents) or a new LEVEL above the
object-level modes is a taxonomic question the data now forces — but it is unmistakably a competency
family the prior (evolutionary, object-level) regime never produced. All three are **known to science**
(meta-RL / contextual policies) → **0 new-to-science**; this is a Tier-2 map enrichment via the capable
learner, the richest face of the arc.

## Bottom line
With a capable learner + a feedback channel, the open-ended generator produces genuine **meta-level**
competencies — and the self-debunking discipline both rescued the probe (the reward-design artifact)
and precisely characterized the result (genuine inference, but one-shot, via ablation + flip). The
honest ceiling holds: the toy regime keeps yielding *known* competencies (now up to meta-learning), not
new-to-science ones. The remaining theoretical route to off-map would be a competency that is not
nameable in the existing literature at all — which toy single-agent gridworlds are unlikely to afford;
reaching it credibly would need a materially richer world (multi-agent, open-ended artifacts, or scale).
