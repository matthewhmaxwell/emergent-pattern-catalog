# Ring-3 demo bundle — six learned competencies for animation

Each folder holds the substrate `.py` (which also defines the policy/network class), trained checkpoint(s), and this
README's per-competency block tells you how to rebuild the net, load it, step one greedy episode, and which per-step
quantities to capture for a faithful 2-D animation. All policies are numpy + torch, CPU-only. No GPU needed.

**Common loading pattern.** Each `*world.py` defines a `Pol` class (except #30, whose policy is `TransformerPolicy`
in `adaptheteroTF_train.py`). Build the world once to read its observation dim, build the policy at that dim, then
`load_state_dict`:

```python
import allocworld as W                      # the substrate module
env = W.AllocWorld(1, 0, cfg)               # B=1 probe instance
net = W.Pol(env.odim if hasattr(env,'odim') else env.obs().shape[-1])
net.load_state_dict(torch.load("division_of_labor_s0.pt")); net.eval()
```

**Greedy rollout** = argmax the policy each step (no sampling), loop the world's horizon `T` (module-level constant
in each file), record state each step. Ablation kwargs on the world constructor are the diagnostics; leave them all
False for the natural behavior.

---

## #22 division-of-labor
- **files:** `allocworld.py` (env + `Pol`), `alloc_dynamics.py` (curriculum trainer used to make the checkpoint), `division_of_labor_s0.pt`
- **policy class:** `Pol(odim, n=None, id_skip=False)` — build with `Pol(AllocWorld(1,0,cfg).odim)`. GRU-recurrent (pass hidden state through the rollout).
- **run one rollout:** `cfg = dict(n_ag=2, has_obs=True, has_sig=False, hetero=False, reward_shared=True, asym_init=False)`; `env = AllocWorld(B=1, seed=1, cfg=cfg)`; grid `N=9`, horizon `T=30`, `NMOVE=5`. Step with the movement head argmax (`env.step(mv, sg)`; `sg` is the signal head, unused when `has_sig=False`). Metric: `env.alloc_score()`. NOTE: this is an unstable attractor — seed 1 is an expressing seed (alloc 0.845); seeds 0/2/3 collapse. Animate the expressing one.
- **record for animation:** agent positions `env.ap` shape `(B, n=2, 2)` on a 9×9 grid; target cells `env.tg` shape `(B, 2, 2)`; per-agent type `env.atype` `(B,2)` and target type `env.ttype` `(B,2)` (only matters if hetero=True). Spatially: two agents start together, then split so each settles on a different target — one-agent-per-target is the behavior; watch the assignment hold once formed.

## #24 stigmergic-coordination
- **files:** `stigworld.py` (env + `Pol`), `stigmergic_coordination_s0.pt`
- **policy class:** `Pol(odim)` — GRU-recurrent.
- **run one rollout:** `cfg = dict(n_ag=4, grid=11, n_tgt=8, decay=0.85)`; `env = StigWorld(B=1, seed=0, cfg=cfg)`; horizon `T=24`, `NMOVE=5`. Step `env.step(mv)`. Metric: `env.cover_score()` (retrained checkpoint reaches 0.983). Diagnostic: `nomark=True` (zero the field each step) collapses it.
- **record for animation:** agent positions `env.ap` `(B, 4, 2)` on the 11×11 torus; the **pheromone field** `env.field` `(B, 11, 11)` — draw this as a decaying heat-map underlay (the key visual); target cells `env.tgt` `(B,11,11)` binary mask; coverage `env.visited` `(B,11,11)`. Spatially: four blind agents lay pheromone trails and steer off each other's marks to divide the torus, covering nearly all target cells without ever seeing one another.

## #26 niche-construction
- **files:** `terraformworld.py` (env + `Pol`), `niche_construction_s0/s1/s2.pt` (3 seeds)
- **policy class:** `Pol(odim)` — GRU-recurrent.
- **run one rollout:** `env = TerraformWorld(B=1, seed=0, n_ag=2)`; grid `N` (module const), `WALL_ROW=N//2`, wall health `HARD=2`, horizon `T=120`. Step `env.step(act)`. Metric: `env.score()` (goals reached). Diagnostic: `revert=True` (undo edits each step) collapses to 0.
- **record for animation:** agent positions `env.ap` `(B, 2, 2)`; the **wall/terrain state** `env.wall` shape `(B, N)` = health per column of the wall row (draw as a breachable barrier that changes over time — this is the constructed artifact); the goal `env.goal` `(B,2)` and which side it's on `env.goal_far` `(B,)`; facing `env.face` `(B,2)`. Spatially: agents dig through a wall (health drops to 0 at breached columns) to open a durable passage, then use it to reach goals that alternate sides.

## #28 momentum-control
- **files:** `momentumworld.py` (env + `Pol`), `momentum_control_s0/s1/s2.pt` (3 seeds; s0 = the 8000-iter scale net)
- **policy class:** `Pol(odim)` — GRU-recurrent.
- **run one rollout:** `env = MomentumWorld(B=1, seed=0, n_ag=2)`; continuous arena `BOX=10.0`, `DT=0.2`, `DRAG=0.06`, `VMAX=4.0`, puck `PUCK_MASS=2.0`, `GOAL_R=2.0`, horizon `T=50`, `NMOVE=5` (thrust directions). Step `env.step(mv)`. Metric: `env.score()`. Diagnostic: `novel=True` (zero own-velocity) collapses it.
- **record for animation:** agent positions `env.ap` `(B, 2, 2)` **continuous** coords in [0,BOX]; agent velocities `env.av` `(B,2,2)`; **puck** position `env.pp` `(B,2)` and velocity `env.pv` `(B,2)`; goal `env.goal` `(B,2)` with radius GOAL_R. Spatially: agents thrust into a puck and, accounting for inertia and drag, drive it across a continuous arena into the goal circle — smooth momentum-carrying motion, not grid steps.

## #29 morphology-specialization
- **files:** `heteroworld.py` (env + `Pol`), `morphology_specialization_s0/s1/s2.pt` (3 seeds, 8000-iter)
- **policy class:** `Pol(odim)` — GRU-recurrent.
- **run one rollout:** `env = HeteroWorld(B=1, seed=0, n_ag=4)`; grid `N`, `NRED=4, NBLUE=4`, speeds `FAST=2.0/SLOW=0.5`, horizon `T=60`. Step `env.step(act)`. Metric: `env.score()`. Diagnostic: `nobody=True` (zero own body one-hot) collapses; `noid=True` does NOT (specializes by morphology, not identity).
- **record for animation:** agent positions `env.ap` `(B, 4, 2)`; **per-agent body type** `env.body` `(B,4)` (0=red, 1=blue — colour the agents by this); red food `env.red` `(B,NRED,2)` with quantities `env.rq` `(B,NRED)`, blue food `env.blue` `(B,NBLUE,2)` with `env.bq`. Spatially: red-bodied agents converge on red food and blue on blue (each body harvests its matching type efficiently), so the population visibly sorts by colour onto matching resources.

## #30 compositional-attention
- **files:** `adapthetero.py` (env), `adaptheteroTF_train.py` (defines `TransformerPolicy`), `compositional_attention_s0/s1/s2.pt` (3 seeds, 8000-iter)
- **policy class:** `TransformerPolicy(odim, d=64, nhead=4, layers=2, maxT=80)` — causal Transformer, NOT a GRU. Rebuild with `TransformerPolicy(AdaptHetero(1,0).odim)` then load_state_dict.
- **run one rollout:** `env = AdaptHetero(B=1, seed=0, n_ag=1)` (single-agent); grid `G=7`, horizon `T=80` (module-level `E.T`), `NACT=5`, opponent block period `PERIOD=20`, colour-switch period `CPERIOD=13`. IMPORTANT interface quirk: `step()` returns a dummy array, so call `env.obs()` AFTER each step; the Transformer needs the growing observation prefix each step (feed `obs[:, :t+1]` with a causal mask, argmax the last position). Metric: `env.score()`. Diagnostic: BOTH `nobody=True` AND `freezeblock=True` are individually necessary — this is the joint-conditioning signature.
- **record for animation:** agent position `env.ap` `(B,2)`; **own body** `env.body` `(B,)` (0 red/1 blue); the two sites `SITES = [[1,3],[5,3]]`; **which site is currently blocked** = call `env._cur_block()` → `(B,)` index (0 or 1), and the **colour each site serves** = call `env._site_color(0)` and `env._site_color(1)` → `(B,)` each (these are METHODS computed from the timestep, not attributes; both switch on their own clocks — PERIOD=20 for the block, CPERIOD=13 for colour); own energy `env.e` `(B,)`. Spatially: a single agent moves to whichever site is (a) currently unblocked by the opponent AND (b) serving its own body colour — tracking two independently-switching signals at once. That dual-tracking is the behavior to show.

---

## Provenance note
#22 and #24 checkpoints were regenerated this session from the archived trainers (no `.pt` had been saved as an
artifact previously); #22 reproduces the unstable-attractor result (seed 1 expresses at 0.845, matching the
finding's 4/6-express pattern), #24 reproduces cover 0.983 (reported 0.98). #26/#28/#29/#30 checkpoints are the
originals behind the reported results (#28/#29/#30 are the 8000-iter scale nets). All gate KNOWN to science.
