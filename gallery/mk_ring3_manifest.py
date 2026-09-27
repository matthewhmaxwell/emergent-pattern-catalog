"""Append ALL 30 Ring-3 competency demos to gallery/manifest.json as a separate 'competency' track,
each with a core-algorithm CODE PANEL. Idempotent (drops any prior C## entries) + backs up manifest.json.
Sprite metadata (frames/cols/rows) is computed from each MP4 via ffprobe, so it never drifts from the
rendered asset. Run from repo root on the VPS (epc-venv)."""
import json, os, shutil, math, subprocess
from PIL import Image
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
MAN = ROOT + "/gallery/manifest.json"
ASSETS = ROOT + "/gallery/assets"
CAT = json.load(open(ROOT + "/analysis/ring3_competency/canonical_catalog.json"))
cat = {c["id"]: c for c in CAT["competencies"]}


def sprite_meta(base):
    mp4 = f"{ASSETS}/ring3_{base}.mp4"
    n = subprocess.run(["ffprobe", "-v", "error", "-count_frames", "-select_streams", "v:0",
                        "-show_entries", "stream=nb_read_frames", "-of", "csv=p=0", mp4],
                       capture_output=True, text=True).stdout.strip()
    frames = int(n)
    cols = math.ceil(math.sqrt(frames))                    # save_animation's grid convention
    rows = math.ceil(frames / cols)
    # derive the TRUE tile size from the sprite sheet: save_animation caps very large sheets, so
    # tiles are not always 300px (e.g. an 8x7 sheet is downscaled). Reading W/cols keeps the JS in sync.
    w, h = Image.open(f"{ASSETS}/ring3_{base}_sprite.png").size
    return dict(frames=frames, cols=cols, rows=rows, fw=round(w / cols), fh=round(h / rows))


# base, summary (what it is), effect, watch (what to look for), module, where, code (core algorithm)
DEMOS = {
 1: dict(base="navigation", module="stateful_funnel.py", where="run_fsm()",
   summary="A trained finite-state controller reaches a goal on unseen mazes, including a trap whose only opening points away from the goal.",
   effect="The agent commits to a detour that temporarily increases distance to the goal — a purely reactive policy oscillates at the barrier.",
   watch="The agent moves AWAY from the goal to clear the wall, then curls back — holding the detour rather than bouncing against the barrier.",
   code='''# a tiny finite-state controller: (sensory_code, state) -> (action, next_state).
# sensory_code = coarse goal-direction x local wall-bits; the internal state s lets the
# agent COMMIT to a detour that points away from the goal (a reactive agent oscillates).
def run_fsm(rule, start, goal, walls, steps=160):
    p = start; s = 0
    for _ in range(steps):
        if p == goal: return True
        code = (dir_code(p, goal) * 16 + wall_bits(p, walls)) * NSTATES + s
        a, s = rule[code]                 # look up (action, next state)
        if free(p + ACTIONS[a], walls): p = p + ACTIONS[a]
    return p == goal
# DIAGNOSTIC: the away-from-goal trap (gap only on the far side). A memoryless reactive
# policy provably caps out; only the state variable s can hold the detour.'''),

 2: dict(base="memory", module="evolve_memory.py", where="episode()",
   summary="See a cue, survive a variable delay full of distractors, then report the cue — working memory.",
   effect="The cue is held across delays far longer than any seen in training; flip the cue and the reported answer flips.",
   watch="The cue appears, vanishes under a run of distractor inputs, then is correctly reported at the decision step.",
   code='''# cue at t0, then a VARIABLE delay of distractors (cue gone), then a decision step.
# a memoryless agent has NO information at decision time -> 50% BY CONSTRUCTION.
def episode(rule, cue, delay, r):
    seq = [cue] + [r.choice(DISTRACT) for _ in range(delay)] + [DECISION]
    s = 0; out = 0
    for inp in seq:
        out, s = rule[inp * NSTATES + s]  # (emit, next state) — cue must persist in s
    return out == cue
# held-out = delays 12-25 (never trained) + a cue-FLIP test (flip cue -> output flips),
# proving it reads stored memory, not a fixed policy.'''),

 3: dict(base="counting", module="counting_ppo.py", where="VecCount.step()",
   summary="Collect EXACTLY k targets from an over-supplied field, then reach the goal — with no oracle count.",
   effect="A running tally is internalized in the GRU hidden state; the agent stops at k and stays supply-invariant when supply doubles.",
   watch="The agent picks up targets one by one, then STOPS (leaving extras uncollected) and heads to the goal.",
   code='''# collect EXACTLY k targets (supply > k, so it must STOP), then reach the goal.
# no oracle count — the running tally lives in the GRU hidden state.
on = (self.obj == self.pos[:, None, :]).all(2) & self.alive
for b in np.where(on.any(1))[0]:                 # picked up a target
    self.count[b] += 1
    rew[b] += 0.5 if self.count[b] <= self.k else -1.0    # reward under k, punish over
atgoal = (self.pos == self.goal).all(1)
for b in np.where(atgoal)[0]:
    rew[b] += 2.0 if self.count[b] == self.k else -0.5    # bonus only for EXACTLY k
# TEST: double the supply -> a real counter still stops at k (supply-invariant).'''),

 4: dict(base="sequencing", module="rnn_es.py", where="episode()",
   summary="Touch station A before goal G — track a prerequisite order.",
   effect="Order is held in the RNN hidden state; once it leaves A a reactive agent oscillates, so memory is load-bearing.",
   watch="Even starting next to G, the agent first detours to station A (which turns green), THEN returns to G.",
   code='''# touch station A BEFORE goal G. Reward 2 only if A was visited first; once you leave A
# a reactive agent oscillates, so the ORDER must be held in the RNN hidden state h.
def episode(th, start, A, G, T=50):
    Wx, Wh, b, Wo, bo = unpack(th); p = start; h = np.zeros(H); done = 0
    for _ in range(T):
        x = np.concatenate([sgn(p, A), sgn(p, G)])   # signs toward A and toward G
        h = np.tanh(Wx @ x + Wh @ h + b)             # recurrent memory
        p = step(p, argmax(Wo @ h + bo))
        if done == 0 and p == A: done = 1            # reached A first
        if done == 1 and p == G: return 2.0          # then G -> full credit
    return float(done)
# memoryless baseline provably caps at 0.49; RNN + evolution strategies reaches ~0.86.'''),

 5: dict(base="delayed", module="evolve_delayed.py", where="episode()",
   summary="Cash out at the peak of a reward that rises then vanishes — with the peak height unknown each episode.",
   effect="No fixed threshold works; cashing at the peak requires noticing the value just turned down (optimal stopping).",
   watch="The tracked value climbs, and the agent cashes out right at the top — before it decays away.",
   code='''# a reward rises then vanishes; its HEIGHT is random per episode, so no fixed threshold
# works. Cashing at the peak needs noticing the value just turned DOWN.
def curve(r):
    p = r.randint(3, T - 3); h = r.uniform(0.4, 1.0); w = r.uniform(2, 4)
    return [max(0.0, h * (1 - abs(t - p) / w)) for t in range(T)], h
def episode(rule, c, h):
    s = 0
    for v in c:
        a, s = rule[bucket(v) * NSTATES + s]   # (cash-out?, next state)
        if a == 1: return v / h                 # 1.0 = took it exactly at the peak
    return 0.0'''),

 6: dict(base="regulation", module="evolve_regulation.py", where="episode()",
   summary="Hold a variable at a setpoint against random disturbances, in a system with inertia.",
   effect="Feedback control: the agent sees position (not velocity) yet damps the momentum-driven overshoot to stay near the setpoint.",
   watch="Disturbances kick the value off zero; the controller pushes back and settles it near the setpoint without wild oscillation.",
   code='''# hold a variable at setpoint 0 against random disturbances, WITH inertia (momentum).
# the agent sees position only, not velocity; a purely reactive push overshoots.
def episode(rule, r, T=40, damping=0.85, force=0.18, sigma=0.25):
    x = r.uniform(-1, 1); v = 0.0; s = 0; sc = 0.0
    for t in range(T):
        a, s = rule[quant(x) * NSTATES + s]      # choose thrust in {-1, 0, +1}
        d = r.gauss(0, sigma)                     # random disturbance
        v = damping * v + force * ACT[a] + d      # momentum integrates the thrust
        x = clip(x + v, -1.5, 1.5)
        sc += max(0.0, 1.0 - abs(x))              # 1 at setpoint, 0 at the rails
    return sc / T'''),

 7: dict(base="infer", module="metaworld_ppo.py", where="VecMeta.step()  (rule='infer')",
   summary="A hidden 'good type' pays off; the agent infers which type from reward feedback, then exploits it (one-shot).",
   effect="Belief over a latent rule: it probes, reads the reward sign, and commits — but does not re-adapt if the rule later flips.",
   watch="Early moves sample different object types; once it gets a hit it sticks to that type for the rest of the episode.",
   code='''# a hidden good TYPE g gives +1; other types give -0.3. g is unknown, so the agent must
# INFER it from reward feedback, then exploit. (One-shot: it does not re-adapt to a flip.)
on = (self.obj == self.pos[:, None, :]).all(2) & self.alive
for b in np.where(on.any(1))[0]:
    t = self.otype[b, j]; good = (t == self.g[b])   # g = hidden good type
    rew[b] += 1.0 if good else -0.3
    self.last_r[b] = 1.0 if good else -1.0          # feedback the agent reads next step
# DIAGNOSTIC: ablate feedback -> 0.33 (chance); mid-episode flip -> 0.15 (no re-adapt).'''),

 8: dict(base="switch", module="metaworld_ppo.py", where="VecMeta.step()  (rule='switch')",
   summary="Like rule-inference, but the good type flips after every success — continual re-tracking.",
   effect="The agent keeps re-inferring the current rule instead of locking one answer, tracking a changing latent.",
   watch="Each time the agent scores, the winning type changes, and it switches which type it pursues next.",
   code='''# like rule-inference, but the good type FLIPS after each success — the agent must keep
# re-inferring the current rule (continual tracking), not lock one answer.
good = (t == self.g[b]); rew[b] += 1.0 if good else -0.3
self.last_r[b] = 1.0 if good else -1.0
if self.rule == "switch" and good:
    self.g[b] = (self.g[b] + 1) % T        # the rule changes the instant you succeed
# feedback-ablation -> chance; the continual tracker reaches ~0.81.'''),

 9: dict(base="cued", module="metaworld_ppo.py", where="VecMeta.obs()  (rule='cued')",
   summary="A start-cue shown only at t=0 names the good type; the agent must hold it and select the matching rule.",
   effect="Conditional policy: the held cue switches which type is correct — remove the cue and a memoryless agent is at chance.",
   watch="A brief cue at the start determines which type the agent then commits to for the whole episode.",
   code='''# a start CUE (shown only at t=0) names the good type; the agent must HOLD it and
# condition its policy on it. After t=0 the cue is gone, so the GRU must store it.
if self.rule == "cued" and self.step_i == 0:
    for b in range(B):
        out[b, cue_off + self.g[b]] = 1.0    # cue: the good type, at t=0 ONLY
# reward +1 for the cued good type; a memoryless MLP without the stored cue is at chance.
# metric 0.85.'''),

 10: dict(base="comms", module="openhunt4_ppo.py", where="VecRef.step()",
   summary="A speaker that can see the target emits a signal; a listener that cannot must reach the target from that signal alone.",
   effect="Emergent referential communication (a Lewis game): the pair co-adapt a signal→location code. Scramble the channel and it collapses.",
   watch="One agent stays put and signals; the other moves to the hidden target guided only by that signal, and both meet there.",
   code='''# speaker sees the target and emits a signal sg; listener sees ONLY the signal (target
# zeroed in its obs). Reward when BOTH reach the target -> the pair must co-adapt a code.
def step(self, mv, sg):
    rew = np.full(B, -0.01, np.float32)
    for i in range(2): self.ap[:, i] = clip(self.ap[:, i] + DIRS[mv[:, i]])
    self.lsig = sg.copy()                        # signal broadcast to the listener's obs
    both = (self.ap[:, 0] == self.tgt).all(1) & (self.ap[:, 1] == self.tgt).all(1)
    rew[both] += 1.0                             # only paid when the listener arrives too
    return rew
# DIAGNOSTIC: MUTE (scramble the signal) -> collapse to chance = the channel is required.'''),

 11: dict(base="roldiv", module="commhunt_ppo.py", where="Vec2.terminal()  (task='role_div')",
   summary="Two agents are rewarded only if they end on DIFFERENT goals — they must anti-coordinate into complementary roles.",
   effect="Role division via mutual position-observation: symmetric agents break symmetry so each takes a distinct goal.",
   watch="The two agents start alike, then peel apart and settle onto different goals — never doubling up on one.",
   code='''# team reward ONLY if the two agents end on DIFFERENT goals (anti-coordination / roles).
def terminal(self):
    g0, g1 = self.on_goal(0), self.on_goal(1)   # which goal each agent is on (-1 = none)
    r = np.zeros((B, 2), np.float32)
    if self.task == "role_div":
        win = (g0 >= 0) & (g1 >= 0) & (g0 != g1)     # both on goals AND on different ones
        r[:, 0] += win; r[:, 1] += win
    return r
# DIAGNOSTIC: blind the partner -> 0.99 -> 0.44 (needs to observe the other to divide).'''),

 12: dict(base="contest", module="contest_ppo.py", where="VecContest (reward shaping + closeness)",
   summary="Two agents compete for one contested resource; they resolve it by who-is-closer instead of colliding.",
   effect="Symmetry-breaking on an uncorrelated asymmetry (the Bourgeois strategy): the nearer agent claims, the other yields.",
   watch="Both approach the contested cell, but one backs off while the nearer one takes it — collisions stay near zero.",
   code='''# two agents, one contested resource. Per-step shaping rewards being on a goal; the
# CONTEST is resolved by relative position (who-is-closer), not by colliding.
r = np.full((B, 2), -0.01, np.float32)
for i in range(2):
    r[:, i] += 0.02 * (self.on_goal(i) >= 0)    # small pull to sit on the resource
# the learned equilibrium: the closer agent claims, the farther yields (uncorrelated
# asymmetry / Bourgeois). DIAGNOSTIC: blind the partner -> 0.63 -> 0.14; collisions ~0.02.'''),

 13: dict(base="fairness", module="fairness_ppo.py", where="Pol + sampled rollout",
   summary="Two symmetric agents share a resource fairly by TAKING TURNS claiming it (alternation).",
   effect="Reactive anti-coordination: each flips its last action so exactly one claims per round — fair turn-taking, broken by sampling.",
   watch="The two agents alternate claim/yield in lock-step — exactly one claims each round, and who claims flips.",
   code='''# two SYMMETRIC agents on identical observations. The fair solution is turn-taking:
# exactly one claims each round, alternating. Symmetry can ONLY be broken stochastically,
# so the policy MUST sample (argmax makes both claim).
O, A, ... = rollout(net, seed, greedy=False)     # greedy=False is the mechanism
one_claim = ((A[0] + A[1]) == 1).mean()          # exactly one claims each round
flips     = (A[:, 1:] != A[:, :-1]).mean()        # agents alternate turns
# selected by one_claim (fairness 0.98) + flips (efficiency 0.93). Neither channel-
# scramble NOR memory-wipe collapses it -> reactive (flip your own last action).'''),

 14: dict(base="comp", module="comp_ppo.py", where="episode()  (speaker + listener)",
   summary="A speaker counts an item-stream and encodes the count in a 3-slot message; a listener decodes the count.",
   effect="Emergent COMPOSITIONAL communication: a 2-symbol alphabet can't name 6 counts in one slot, so a multi-slot numeral code emerges.",
   watch="Items stream past the speaker; it emits a short multi-slot message, and the listener reports the exact count.",
   code='''# speaker (GRU) counts a K-item stream -> emits an M-slot, S-symbol message.
# listener (MLP) sees ONLY the message -> must output the count. Shared reward = exact count.
def episode(spk, lis, B, seed, scramble=False, countablate=False, dropslot=None):
    count = rng.integers(0, K + 1, size=B)        # uniform 0..K forces a multi-slot code
    stin  = stream_of(count)
    if countablate: stin = np.zeros_like(stin)    # speaker sees nothing -> can't count
    msg = spk(stin)                               # the emergent numeral message
    if scramble:      msg = random_message()      # channel-scramble -> collapse
    if dropslot is not None: msg[:, dropslot] = 0 # per-slot drop tests compositionality
    pred = lis(msg)
    return (pred == count).astype(float)          # 0.84 (chance 0.17); 3/3 slots informative'''),

 15: dict(base="stig", module="stigmergy_ppo.py", where="VecStig.step()",
   summary="Two BLIND agents (no channel, no mutual view) forage by reading and writing a shared mark trail.",
   effect="Stigmergy: coordination is externalized into a trail field in the world; each reads the partner's marks to divide territory.",
   watch="Two agents that never see each other spread apart to cover the food, steering off the trails left behind.",
   code='''# two BLIND agents forage a grid; they read+write a shared MARK trail (state lives in the
# WORLD, not in memory). Reading the partner's trail is what lets them divide territory.
for i in range(2):
    self.p[:, i] = clip(self.p[:, i] + DIRS[acts[:, i]])
    x, y = self.p[:, i, 0], self.p[:, i, 1]
    got = self.food[arange(B), x, y]
    rew += got; self.food[arange(B), x, y] = False    # eat food here
    self.visited[arange(B), x, y] = True
    self.vis_by[arange(B), x, y] |= (i + 1)           # mark WHO visited (the trail)
return rew - 0.01
# DIAGNOSTIC: no-trail -> collapse; own-trail-only 5.7 vs combined 8.1 (partner's matters).'''),

 16: dict(base="tool", module="tool_ppo.py", where="VecTool.step()",
   summary="Push a block into an impassable gap to build a bridge, then cross it — the goal is unreachable otherwise.",
   effect="Instrumental construction: the agent performs a non-rewarding sub-goal (build) to make the rewarding goal reachable.",
   watch="The agent shoves a block into the gap (a bridge forms), then walks across it to the goal it couldn't reach before.",
   code='''# the goal sits across an impassable GAP. The agent must PUSH a block into the gap to
# build a bridge (a non-rewarding sub-goal), then cross. No-block control -> 0.00.
if (not self.used[b]) and (tx, ty) == tuple(self.bp[b]):     # step into the block = push
    bx, by = self.bp[b] + d
    if by == GAP and not self.filled[b, bx]:
        self.filled[b, bx] = True; self.used[b] = True
        rew[b] += 0.5                                        # bridge built (sub-goal)
        self.ap[b] = (tx, ty)
if tuple(self.ap[b]) == tuple(self.goal[b]):
    rew[b] += 1.0; self.done[b] = True                       # goal, reachable only via bridge'''),

 17: dict(base="reputation", module="reputation_ppo.py", where="rollout()",
   summary="Recurring partners each have a hidden type; help known cooperators and refuse known defectors.",
   effect="Indirect reciprocity via a reputation map: the agent tracks each partner-id's revealed history and discriminates on it.",
   watch="The agent helps partners it has learned are cooperators and passes on ones it has learned defect.",
   code='''# recurring partners, each a hidden cooperator/defector. Help a coop -> +1; help a
# defector -> -1; pass -> 0. After each round the partner's type is REVEALED.
a = net(obs, h)                                   # HELP (1) or PASS (0), given partner id
coop = ptype[arange(B), pid]                       # this partner's hidden type
r = np.where(a == 1, np.where(coop, 1.0, -1.0), 0.0)
last_pid, last_type = pid, np.where(coop, 1.0, -1.0)   # written into a reputation memory
# DIAGNOSTIC: anonymize partners -> 0.50 (needs WHO); memory-wipe -> 0.50 (needs history).'''),

 18: dict(base="abstraction", module="abstraction_ppo.py", where="VecAbs.step()",
   summary="Go to the goal whose colour-vector is most similar to a cue — a relation that partially transfers to unseen colours.",
   effect="Relational abstraction: a learned similarity(cue, goal) generalizes above chance to colours seen only as distractors in training.",
   watch="Given a cue colour, the agent walks to the goal whose colour best matches it, even for held-out colours.",
   code='''# go to the goal whose COLOUR-VECTOR is most similar to the cue's. The comparison must
# transfer to held-out colours (seen only as distractors in training).
for g in range(G):
    if tuple(self.ap[b]) == tuple(self.gpos[b, g]):
        win = (g == self.match[b])                # match = argmax similarity(cue, goal_g)
        rew[b] += 1.0 if win else -0.5
        self.done[b] = True; break
# GENERALIZATION GAP (not an ablation): train 0.88 -> held-out 0.59 (chance 0.33) = PARTIAL
# transfer. Pure memorization would leave held-out AT chance; it did not.'''),

 19: dict(base="imitation", module="imitation_ppo.py", where="VecImit.step()",
   summary="A scripted demonstrator walks toward the hidden-correct goal; the learner must identify and go to the goal it chose.",
   effect="Goal emulation: the correct goal is not in the learner's own observation — the other agent's behaviour is the only signal.",
   watch="The learner watches the demonstrator's path and follows it to the goal the demonstrator picked.",
   code='''# a scripted demonstrator walks toward the hidden-correct goal. The correct goal is NOT
# in the LEARNER's observation -> the demonstrator's choice is the only signal.
if self.mode != "frozen":
    self.dp[b] = self.dp[b] + greedy_step(self.dp[b], self.gpos[b, self.demo_goal[b]])
...
if tuple(self.ap[b]) == tuple(self.gpos[b, g]):
    win = (g == self.correct[b]); rew[b] += 1.0 if win else -0.3
# DIAGNOSTIC: FROZEN demonstrator -> 0.08 (unfindable); MISLEADING demo -> 0.08 (the
# learner follows it to the WRONG goal = it copies the demonstrator's choice).'''),

 20: dict(base="tom", module="tom_ppo.py", where="VecToM.step()",
   summary="Watch a mover travel and stop; attribute its goal from its heading — the frozen frame alone is symmetric.",
   effect="Predictive intention-reading: only the integrated MOTION identifies the goal (the three goals are equidistant from the stop).",
   watch="A mover slides toward one of three equidistant goals and stops; the agent picks the goal in its direction of travel.",
   code='''# a mover steps along its heading for K steps then FREEZES on a stop that is EQUIDISTANT
# from all 3 goals (a provably symmetric single frame). Only the MOTION identifies g*.
if self.mode != "frozen" and self.t < K:
    self.mpos = clip(self.mpos + self.sv)         # mover travels along its true heading
if self.t >= K:                                   # after it stops, the agent attributes:
    correct = (a == self.gstar); rew = correct.astype(float)   # goal in the heading dir
# DIAGNOSTIC: memoryless baseline provably at chance 0.32 (symmetric frame); recurrent
# 0.86. Freeze the mover (no motion) -> 0.32.'''),

 21: dict(base="culture", module="culture_scale_ppo.py", where="run()  (inherit vs no-inheritance)",
   summary="Each generation inherits the culture so far, copies it, and explores one new frontier symbol — the recipe is solved across generations.",
   effect="A cumulative-culture ratchet: discoveries accumulate in a buffer that outlives any individual; a single lifetime cannot find the recipe.",
   watch="The known-correct prefix grows generation by generation as each learner adds one symbol beyond the inherited buffer.",
   code='''# each generation INHERITS the culture prefix, copies it (imitation), and explores the
# FRONTIER symbol (innovation). Discoveries accumulate across lifetimes.
obs = make_obs(recipe, p if inherit else 0, L)    # p = inherited known-prefix length
sym = policy(obs); pl = prefix_len(sym, recipe, L)
if inherit: p = np.maximum(p, pl)                 # the buffer ratchets forward, never back
# DIAGNOSTIC (collapse-ablation at the CROSS-GENERATIONAL scale): remove inheritance ->
# 11.69/12 (97%) collapses to 0.38/12 (3%). A single generation can't find a length-12 recipe.'''),

 22: dict(base="division", module="allocworld.py", where="AllocWorld.step()",
   summary="n identical agents must each reach a DIFFERENT target — emergent task allocation (a perfect matching) with no assigned roles.",
   effect="Two symmetric agents self-organize into complementary roles — one per target — anchored on agent identity.",
   watch="Two agents start together, then split so each settles on a different target; the assignment holds once formed.",
   code='''# n homogeneous agents, n targets. Reward each agent for sitting ALONE on a target and
# penalize collisions -> the team is pushed to a perfect matching (division of labor).
for i in range(n):
    on_t = any(all(self.ap[:, i] == self.tg[:, k], 1) for k in range(n))   # covers a target
    same = any(all(self.ap[:, i] == self.ap[:, j], 1) for j in range(n) if j != i)
    rew[:, i] += 0.1 * on_t - 0.05 * same          # cover a target, avoid doubling up
if self.reward_shared: rew[:] = rew.mean(1, keepdims=True)
# DIAGNOSTIC: memwipe + blind + noid all collapse it; teleport one agent off its target
# mid-episode -> the allocation RE-FORMS (recovery 0.90) = actively maintained.'''),

 23: dict(base="group_reg", module="group_reg_build.py", where="rollout() reward",
   summary="n agents each choose active/inactive; the active COUNT must track a target K*(t) that switches every few rounds.",
   effect="Homeostatic set-point tracking: agents break symmetry by identity into a stable active subset and re-regulate when K* moves.",
   watch="The dashed target line jumps around; the active-count line steps to meet it each time it moves.",
   code='''# n agents each pick active/inactive. Shared reward peaks when the ACTIVE COUNT equals a
# moving target K*(t). Agents use identity to take stable complementary roles that sum to K*.
x = concat([onehot_id, K_star / N, last_count / N, bias])   # per-agent obs (identity + K*)
logits, h = net(x, h)
a = Categorical(logits).sample()                 # sampling breaks symmetry (who is active)
count = int(a.sum())
rew = 1.0 - abs(count - K_star) / N              # peaks at exact match
# DIAGNOSTIC: drop the identity channel (noid) -> agents identical -> can only hit p*N in
# expectation, never the exact moving target -> collapse (necessity 1.00).'''),

 24: dict(base="stigmergy", module="stigworld.py", where="StigWorld.step()",
   summary="n agents on a torus cover target cells; they cannot see each other and coordinate only through a decaying pheromone field.",
   effect="Environment-mediated coordination: coordination lives in the field, not in memory or identity (the profile inverts role-partitioning).",
   watch="Four blind agents lay pheromone trails and steer off each other's marks to divide the torus, covering nearly all targets.",
   code='''# n agents on a torus; they coordinate ONLY through a shared, decaying pheromone field.
self.field *= self.decay                          # pheromone evaporates
for i in range(n):
    self.field[arange(B), *self.ap[:, i]] += DEPOSIT   # deposit at the agent's cell
if self.nomark: self.field[:] = 0.0               # (ablation: destroy the field)
new = count_newly_covered_targets(self.ap)        # reward = NEW targets covered this step
rew = repeat(new / self.n_tgt, n)
# DIAGNOSTIC INVERSION: nomark load-bearing (0.54); noid ~0.00 and memwipe ~0.02 are DEAD
# (agents interchangeable, no memory) — the opposite of the identity-anchored competencies.'''),

 25: dict(base="metabolic", module="metabolic_build.py", where="Commons.step()",
   summary="Survive by reading your own internal energy and harvesting a regrowing commons sustainably.",
   effect="Homeostatic foraging: a cell driven to zero dies, so the agent must harvest below the regrowth rate — sustainability is discovered, not rewarded.",
   watch="The agent nibbles cells (they pale then regrow), the food row stays green, and the energy bar stays full.",
   code='''# survive T steps: metabolism drains energy each step; eating refills it. A food cell
# driven to 0 is DEAD forever (the commons pressure), so harvest must stay below regrowth.
def step(self, a):
    if   a == 0: self.p = max(0, self.p - 1)
    elif a == 1: self.p = min(L - 1, self.p + 1)
    elif self.f[self.p] > 0:                       # eat here
        self.f[self.p] -= 1; self.e = min(EMAX, self.e + EAT)
    self.e -= METAB
    if self.t % REGROW == 0:
        self.f = where(self.f > 0, minimum(FMAX, self.f + 1), 0.0)   # dead cells stay dead
    return self.e > 0
# DIAGNOSTIC: zero the own-energy channel (novel) -> can't time intake -> collapse.'''),

 26: dict(base="niche", module="terraformworld.py", where="TerraformWorld.step()",
   summary="Durably modify the environment — dig a breach through a wall — to make an otherwise-unreachable goal reachable.",
   effect="Niche construction: agents reshape the world (a breached wall) with no reward for building; undo the edits each step and it collapses.",
   watch="Agents dig through a wall (its columns lose health and breach), then use the passage to reach goals that alternate sides.",
   code='''# the goal alternates sides of a wall. Agents DIG (action 5) to reduce wall-cell health;
# a breach is a durable passage. There is NO reward for digging itself.
dig = (a == 5); tr, tc = cell_faced(self.ap[:, i], self.face[:, i])
hit = dig & (tr == WALL_ROW) & (self.wall[arange(B), tc] > 0)
self.wall[arange(B), tc] -= hit                   # wall health drops -> eventually breaches
if self.revert: self.wall = self.wall0.copy()     # (ablation: edits don't persist)
rew += on_goal                                    # reward only for REACHING the goal
# DIAGNOSTIC: revert (undo edits each step) -> 0.00; the durable world-change is essential.'''),

 27: dict(base="nonstat", module="nonstat_build.py", where="rollout()  (opponent switch)",
   summary="A co-player follows a periodic pattern, then switches to a new one mid-episode; the agent must detect it and re-adapt.",
   effect="Online opponent-modelling with switch-readiness: it re-infers the co-player's pattern after a regime change, unlike a fixed model.",
   watch="Counter-accuracy is high, dips sharply right at the switch, then climbs back as the agent re-locks onto the new pattern.",
   code='''# the co-player plays a fixed periodic pattern, then SWITCHES at mid-episode. The agent
# scores when its move counters the co-player's current move -> it must model + re-model.
seq = [p1[t % len(p1)] if t < SWITCH else p2[t % len(p2)] for t in range(T)]
x = concat([onehot(prev_opp_move), onehot(prev_own_move), bias])   # opponent history
a = Categorical(net(x, h)).sample()
hit = int(a == (seq[t] + 1) % A)                  # did the agent's move counter the co-player?
# DIAGNOSTIC: zero the opponent-history channel -> collapse; a STATIONARY-trained agent
# (never saw a switch) holds phase 1 (0.82) but fails phase 2 (0.34).'''),

 28: dict(base="momentum", module="momentumworld.py", where="MomentumWorld.step()",
   summary="Steer a puck that carries momentum into a goal, against inertia and drag — genuine continuous-physics control.",
   effect="Own-velocity conditioning under continuous physics: agents transfer momentum to the puck rather than stepping it on a grid.",
   watch="Agents thrust into the puck and, accounting for inertia and drag, drive it across the arena into the goal circle.",
   code='''# continuous physics: thrust updates velocity (with drag); contact transfers momentum to
# the puck. Reward = the puck's progress toward the goal (not the striking itself).
self.av = clip_v((self.av + thrust * DT) * (1 - DRAG))
self.ap = clip(self.ap + self.av * DT, 0, BOX)
# on contact, project agent velocity along the normal into the puck:
self.pv = where(hit, self.pv + max(vproj, 0) * nrm / PUCK_MASS, self.pv)
puck_prog = dot(self.pv, unit(goal - self.pp))    # puck speed toward the goal
rew = puck_prog * 0.3 + (dist(puck, goal) < GOAL_R) * 2.0
# DIAGNOSTIC: zero the own-velocity channel -> collapse (necessity 0.53-0.63).'''),

 29: dict(base="morphology", module="heteroworld.py", where="HeteroWorld.step()",
   summary="Agents of different body types each specialize on the resource their morphology harvests best.",
   effect="Specialization by BODY, not identity index: red-bodied agents converge on red food, blue on blue — the population sorts by morphology.",
   watch="Red-body agents converge on red food and blue-body on blue — the mixed population sorts by colour onto matching resources.",
   code='''# each agent has a fixed BODY type. It harvests its matching resource for full value and
# the mismatched resource for none -> specialize on what your morphology harvests best.
for i in range(n):
    body = self.body[:, i]
    on_red  = harvested(self.ap[:, i], self.red)
    on_blue = harvested(self.ap[:, i], self.blue)
    rew += on_red  * (body == 0)                   # body-0 credited for red food
    rew += on_blue * (body == 1)                   # body-1 credited for blue food
# DIAGNOSTIC: zero the own-body one-hot (nobody) -> collapse; noid does NOT collapse
# (it specializes by BODY, not by identity index).'''),

 30: dict(base="compositional", module="adapthetero.py", where="AdaptHetero.step()  (Transformer policy)",
   summary="A single policy conditions on TWO independently-switching signals at once — its own body colour and the opponent's block.",
   effect="Joint conditioning that emerges only under an attention (Transformer) policy: a recurrent bottleneck holds one channel, not both.",
   watch="The agent moves to the site that is both unblocked by the opponent and serving its own body colour — tracking two signals at once.",
   code='''# feed only at the site that is (a) NOT blocked this step AND (b) serving your OWN body
# colour -> the policy must condition on two independently-switching channels at once.
blk = self._cur_block(); openidx = 1 - blk
for si in range(2):
    on = at_site(self.ap, si)
    canfeed = on & (openidx == si) & (self._site_color(si) == self.body) & self.alive
    self.e = where(canfeed, min(self.e + INTAKE, EMAX), self.e)
self.alive &= (self.e > 0)
# emerges ONLY with a Transformer (attention over obs-history; a GRU holds one channel).
# DIAGNOSTIC: nobody -> partial collapse AND freezeblock -> partial collapse (both needed).'''),
}


def short_ref(lit):
    for sep in (" (", ";", " / "):
        if sep in lit: lit = lit.split(sep)[0]
    return lit.strip()[:52]


m = json.load(open(MAN))
m = [e for e in m if not str(e.get("id", "")).startswith("C")]      # idempotent: drop prior C##
added = 0
for cid, D in DEMOS.items():
    c = cat[cid]
    m.append({
        "id": f"C{cid}", "name": c["name"], "ref": short_ref(c["lit_name"]),
        "summary": D["summary"], "effect": D["effect"], "viz": "point_cloud", "watch": D["watch"],
        "metric": c["metric"], "mechanism": c["mechanism"], "where": "",
        "detector": {"note": f"gate: KNOWN ({short_ref(c['lit_name'])}). Diagnostic — {c['diagnostic']}"},
        "asset": f"ring3_{D['base']}_sprite.png", "asset_type": "sprite", "sprite": sprite_meta(D["base"]),
        "mp4": f"ring3_{D['base']}.mp4", "contact_sheet": f"ring3_{D['base']}_sprite.png", "track": "competency",
        "code": D["code"], "code_module": D["module"], "code_where": D["where"],
    })
    added += 1
shutil.copy(MAN, MAN + ".bak")
json.dump(m, open(MAN, "w"), indent=2)
print("manifest entries:", len(m), "| competency:", sum(1 for e in m if e.get("track") == "competency"),
      "| added:", added, "| with code:", sum(1 for e in m if e.get("track") == "competency" and e.get("code")))
