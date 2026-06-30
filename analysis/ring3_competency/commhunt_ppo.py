"""Phase D2: the open-ended MULTI-AGENT hunt. Generalizes D1's referential game to a GENERATOR over
2-mobile-agent task structures, to search the qualitatively-larger multi-agent competency space:
  - referential : agent0 privately sees a target color; team reward iff agent1 reaches the target goal
                  (demands COMMUNICATION).
  - coordination: team reward iff BOTH agents end on the SAME goal (demands CONSENSUS).
  - role_div    : team reward iff the agents end on DIFFERENT goals (demands ANTI-coordination / roles).
  - independent : each agent has its own private target; reward is per-agent (CONTROL: no multi-agent
                  competency demanded).
Two agents share a dual-head policy (move + symbol) + value; each observes its goal-directions, the other
agent's relative position, the other's last symbol, an agent-id, and (referential/independent) its private
target. Joint PPO. Discriminators (the agent-observer's debunks, applied to the trained policy):
  - channel-scramble: feed random symbols -> collapse => the task genuinely used COMMUNICATION.
  - blind-partner: zero the other-agent channels (rel-pos + symbol) -> collapse => uses the PARTNER.
A structure is "multi-agent-demanding" if normal success is high AND an ablation collapses it toward the
no-coordination floor. Demanding structures then go to the agent-observer to NAME + classify (vs map or
OFF-MAP), literature-gated. Honest bound: all of {communication, consensus, role-division} are KNOWN ->
Tier-2; the hunt is for a combination/structure that resists every name.

Run on VPS epc-venv: python commhunt_ppo.py [--iters 500] [--tasks referential,coordination,role_div,independent]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 5; G = 3; K = 4; MAXT = 18; H = 64
ODIM = G * 2 + 2 + K + 1 + G     # goal-dirs, other rel-dir, other last-symbol, agent-id, private-target
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Agent(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(ODIM, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.move = nn.Linear(H, 5); self.sym = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x):
        z = self.net(x); return self.move(z), self.sym(z), self.v(z).squeeze(-1)


class Vec2:
    def __init__(self, B, task, seed):
        self.B, self.task = B, task; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.goal = np.zeros((B, G, 2), int); self.p = np.zeros((B, 2, 2), int)
        for b in range(B):
            idx = self.rng.permutation(N * N)[:G + 2]
            self.goal[b, :, 0], self.goal[b, :, 1] = idx[:G] % N, idx[:G] // N
            self.p[b, 0] = (idx[G] % N, idx[G] // N); self.p[b, 1] = (idx[G + 1] % N, idx[G + 1] // N)
        self.tgt = self.rng.integers(0, G, size=B)        # referential target / agent0's own target
        self.tgt1 = self.rng.integers(0, G, size=B)       # independent: agent1's own target
        self.last_sym = np.zeros((B, 2), int)
        return None

    def on_goal(self, i):                                  # goal index each of B agents-i is on, else -1
        eq = (self.goal == self.p[:, i, None, :]).all(2)   # (B,G)
        idx = np.where(eq.any(1), eq.argmax(1), -1)
        return idx

    def obs(self, i, scramble=False, blind=False):
        B = self.B; out = np.zeros((B, ODIM), np.float32); o = 0
        for g in range(G):
            out[:, o:o + 2] = np.sign(self.goal[:, g, :] - self.p[:, i, :]); o += 2
        oth = 1 - i
        rel = np.sign(self.p[:, oth, :] - self.p[:, i, :])
        out[:, o:o + 2] = 0.0 if blind else rel; o += 2
        sym = self.rng.integers(0, K, size=B) if scramble else self.last_sym[:, oth]
        if not blind:
            out[np.arange(B), o + sym] = 1.0
        o += K
        out[:, o] = i; o += 1                              # agent id
        if self.task == "referential" and i == 0:
            out[np.arange(B), o + self.tgt] = 1.0          # only agent0 sees the target
        if self.task == "independent":
            out[np.arange(B), o + (self.tgt if i == 0 else self.tgt1)] = 1.0
        o += G
        return out

    def step(self, mv):                                    # mv: (B,2) moves
        B = self.B
        for i in range(2):
            np_ = self.p[:, i, :] + DIRS[mv[:, i]]
            inb = (np_[:, 0] >= 0) & (np_[:, 0] < N) & (np_[:, 1] >= 0) & (np_[:, 1] < N)
            self.p[:, i, :] = np.where(inb[:, None], np_, self.p[:, i, :])
        r = np.full((B, 2), -0.01, np.float32)
        for i in range(2):
            r[:, i] += 0.02 * (self.on_goal(i) >= 0)       # shaping: be on a goal
        return r

    def terminal(self):                                    # task reward at episode end, per agent (B,2)
        B = self.B; g0, g1 = self.on_goal(0), self.on_goal(1); r = np.zeros((B, 2), np.float32)
        if self.task == "referential":
            win = (g1 == self.tgt); r[:, 0] += win; r[:, 1] += win
        elif self.task == "coordination":
            win = (g0 >= 0) & (g0 == g1); r[:, 0] += win; r[:, 1] += win
        elif self.task == "role_div":
            win = (g0 >= 0) & (g1 >= 0) & (g0 != g1); r[:, 0] += win; r[:, 1] += win
        else:  # independent
            r[:, 0] += 0.5 * (g0 == self.tgt); r[:, 1] += 0.5 * (g1 == self.tgt1)
        return r

    def success(self):
        g0, g1 = self.on_goal(0), self.on_goal(1)
        if self.task == "referential": return float((g1 == self.tgt).mean())
        if self.task == "coordination": return float(((g0 >= 0) & (g0 == g1)).mean())
        if self.task == "role_div": return float(((g0 >= 0) & (g1 >= 0) & (g0 != g1)).mean())
        return float((0.5 * (g0 == self.tgt) + 0.5 * (g1 == self.tgt1)).mean())


def gae(R, V, gamma=0.99, lam=0.95):
    A, Tt = R.shape; adv = np.zeros((A, Tt), np.float32); last = np.zeros(A, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(A, np.float32)
        delta = R[:, t] + gamma * nextv - V[:, t]
        last = delta + gamma * lam * last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False, scramble=False, blind=False):
    B = env.B; env.reset()
    O = np.zeros((B, 2, MAXT, ODIM), np.float32); MV = np.zeros((B, 2, MAXT), int); SY = np.zeros((B, 2, MAXT), int)
    LP = np.zeros((B, 2, MAXT), np.float32); V = np.zeros((B, 2, MAXT), np.float32); R = np.zeros((B, 2, MAXT), np.float32)
    for t in range(MAXT):
        for i in range(2):
            ob = env.obs(i, scramble=scramble, blind=blind); O[:, i, t] = ob
            with torch.no_grad():
                ml, sl, v = net(torch.from_numpy(ob))
            md, sd = torch.distributions.Categorical(logits=ml), torch.distributions.Categorical(logits=sl)
            mv = ml.argmax(1) if greedy else md.sample()
            sy = sl.argmax(1) if greedy else sd.sample()
            MV[:, i, t] = mv.numpy(); SY[:, i, t] = sy.numpy()
            LP[:, i, t] = (md.log_prob(mv) + sd.log_prob(sy)).numpy(); V[:, i, t] = v.numpy()
        r = env.step(np.stack([MV[:, 0, t], MV[:, 1, t]], 1))
        env.last_sym = np.stack([SY[:, 0, t], SY[:, 1, t]], 1)
        R[:, :, t] = r
    R[:, :, -1] += env.terminal()
    return O, MV, SY, LP, V, R


def train(task, iters, seed=0, B=320):
    torch.manual_seed(seed); net = Agent(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = Vec2(B, task, 1000 + seed * 99999 + it)
        O, MV, SY, LP, V, R = rollout(net, env)
        # flatten agents -> (2B, MAXT)
        O2 = O.transpose(0, 1, 2, 3).reshape(2 * B, MAXT, ODIM)
        MV2, SY2, LP2 = MV.reshape(2 * B, MAXT), SY.reshape(2 * B, MAXT), LP.reshape(2 * B, MAXT)
        V2, R2 = V.reshape(2 * B, MAXT), R.reshape(2 * B, MAXT)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, MVt, SYt, LPt = torch.from_numpy(O2), torch.from_numpy(MV2), torch.from_numpy(SY2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            ml, sl, v = net(Ot.reshape(-1, ODIM))
            ml = ml.reshape(2 * B, MAXT, 5); sl = sl.reshape(2 * B, MAXT, K); v = v.reshape(2 * B, MAXT)
            md, sd = torch.distributions.Categorical(logits=ml), torch.distributions.Categorical(logits=sl)
            nlp = md.log_prob(MVt) + sd.log_prob(SYt)
            ratio = torch.exp(nlp - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((v - rett) ** 2).mean() - 0.01 * (md.entropy() + sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            ev = np.mean([_succ(net, task, 7000 + it + k) for k in range(2)])
            print(f"  [{task}] iter {it:>3}: success {ev:.2f}", flush=True)
    return net


def _succ(net, task, seed, scramble=False, blind=False, B=400):
    env = Vec2(B, task, seed); rollout(net, env, greedy=True, scramble=scramble, blind=blind); return env.success()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 500
    tasks = a[a.index("--tasks") + 1].split(",") if "--tasks" in a else ["referential", "coordination", "role_div", "independent"]
    here = __file__.rsplit("/", 1)[0]; out = f"{here}/commhunt_sweep.json"
    prior = {r["task"]: r for r in json.load(open(out))} if os.path.exists(out) else {}
    results = []
    print(f"PHASE D2 multi-agent hunt: tasks={tasks}, iters={iters} (resuming {sorted(prior)})", flush=True)
    for task in tasks:
        if task in prior:
            results.append(prior[task]); print(f"{task}: (cached)", flush=True); continue
        net = train(task, iters)
        norm = np.mean([_succ(net, task, 9000 + k) for k in range(3)])
        scr = np.mean([_succ(net, task, 9100 + k, scramble=True) for k in range(3)])
        bld = np.mean([_succ(net, task, 9200 + k, blind=True) for k in range(3)])
        comms = norm >= 0.6 and scr <= norm - 0.2
        partner = norm >= 0.6 and bld <= norm - 0.2
        torch.save(net.state_dict(), f"{here}/commhunt_{task}.pt")
        results.append({"task": task, "success": round(float(norm), 2), "scrambled": round(float(scr), 2),
                        "blind": round(float(bld), 2), "uses_comms": bool(comms), "uses_partner": bool(partner)})
        json.dump(results, open(out, "w"), indent=1)
        print(f"{task:13s}: success {norm:.2f} | chan-scramble {scr:.2f} | blind-partner {bld:.2f} "
              f"| comms={comms} partner={partner}", flush=True)
    json.dump(results, open(out, "w"), indent=1); print("saved commhunt_sweep.json", flush=True)
