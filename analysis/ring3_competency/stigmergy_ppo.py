"""Track-1 probe: STIGMERGY -- coordination mediated by the ENVIRONMENT, not a channel or mutual
observation. Two agents forage a grid for hidden food. They CANNOT see each other and have NO direct
channel; they do NOT observe food directly. Every cell an agent steps on is MARKED (a shared trail), and
each agent observes only the mark-state of its current + 4 adjacent cells. To collect food efficiently the
agents must steer toward UNMARKED (unexplored) ground and away from already-covered territory -- including
the PARTNER's trail. Coordination therefore rides entirely on the shared mark field (the ant-pheromone
mechanism). Shared reward = food collected (each cell once).

Discriminator: HIDE the marks (zero the mark bits) -> the agents can no longer read the trail -> they
re-cover each other's territory -> joint food drops and overlap rises. If joint-food(marks) >> joint-food
(no-marks) AND overlap(marks) < overlap(no-marks), the coordination is genuinely stigmergic.

Run on VPS epc-venv: python stigmergy_ppo.py [--iters 700]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 7; F = 12; T = 18; H = 64; ODIM = 2 + 5; NACT = 5
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Agent(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(ODIM, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x):
        z = self.net(x); return self.pi(z), self.v(z).squeeze(-1)


class VecStig:
    def __init__(self, B, seed, hide=False, own_only=False):
        self.B, self.hide, self.own_only = B, hide, own_only; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.food = np.zeros((B, N, N), bool); self.visited = np.zeros((B, N, N), bool)
        self.vis_by = np.zeros((B, N, N), np.int8)            # which agent(s) visited (1,2,3=both)
        self.p = np.zeros((B, 2, 2), int)
        for b in range(B):
            cells = self.rng.permutation(N * N)[:F + 2]
            fx, fy = cells[:F] % N, cells[:F] // N; self.food[b, fx, fy] = True
            self.p[b, 0] = (cells[F] % N, cells[F] // N); self.p[b, 1] = (cells[F + 1] % N, cells[F + 1] // N)
        for i in range(2):
            self.visited[np.arange(B), self.p[:, i, 0], self.p[:, i, 1]] = True
            self.vis_by[np.arange(B), self.p[:, i, 0], self.p[:, i, 1]] |= (i + 1)
        self.collected = np.zeros(B, int)
        return None

    def obs(self, i):
        B = self.B; out = np.zeros((B, ODIM), np.float32)
        out[:, 0] = self.p[:, i, 0] / (N - 1); out[:, 1] = self.p[:, i, 1] / (N - 1)
        if not self.hide:
            for k, (dx, dy) in enumerate(DIRS):                # current + 4 neighbours
                cx, cy = self.p[:, i, 0] + dx, self.p[:, i, 1] + dy
                inb = (cx >= 0) & (cx < N) & (cy >= 0) & (cy < N)
                cxc, cyc = np.clip(cx, 0, N - 1), np.clip(cy, 0, N - 1)
                if self.own_only:
                    marked = ((self.vis_by[np.arange(B), cxc, cyc] & (i + 1)) > 0) & inb   # only THIS agent's trail
                else:
                    marked = self.visited[np.arange(B), cxc, cyc] & inb                     # combined trail
                out[:, 2 + k] = np.where(inb, marked, 1.0)      # OOB counts as marked (avoid walls)
        return out

    def step(self, acts):                                      # acts (B,2)
        B = self.B; rew = np.zeros(B, np.float32)
        for i in range(2):
            np_ = self.p[:, i, :] + DIRS[acts[:, i]]
            inb = (np_[:, 0] >= 0) & (np_[:, 0] < N) & (np_[:, 1] >= 0) & (np_[:, 1] < N)
            self.p[:, i, :] = np.where(inb[:, None], np_, self.p[:, i, :])
            x, y = self.p[:, i, 0], self.p[:, i, 1]
            got = self.food[np.arange(B), x, y]
            rew += got; self.food[np.arange(B), x, y] = False
            self.collected += got
            self.visited[np.arange(B), x, y] = True
            self.vis_by[np.arange(B), x, y] |= (i + 1)
        return rew - 0.01

    def overlap(self):
        return float((self.vis_by == 3).sum(axis=(1, 2)).mean())   # cells visited by BOTH agents


def gae(R, V, gamma=0.99, lam=0.95):
    A, Tt = R.shape; adv = np.zeros((A, Tt), np.float32); last = np.zeros(A, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(A, np.float32)
        delta = R[:, t] + gamma * nextv - V[:, t]
        last = delta + gamma * lam * last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B = env.B
    O = np.zeros((B, 2, T, ODIM), np.float32); A = np.zeros((B, 2, T), int)
    LP = np.zeros((B, 2, T), np.float32); V = np.zeros((B, 2, T), np.float32); R = np.zeros((B, 2, T), np.float32)
    for t in range(T):
        acts = np.zeros((B, 2), int)
        for i in range(2):
            ob = env.obs(i); O[:, i, t] = ob
            with torch.no_grad():
                lg, v = net(torch.from_numpy(ob))
            d = torch.distributions.Categorical(logits=lg); a = lg.argmax(1) if greedy else d.sample()
            A[:, i, t] = a.numpy(); LP[:, i, t] = d.log_prob(a).numpy(); V[:, i, t] = v.numpy(); acts[:, i] = a.numpy()
        r = env.step(acts); R[:, 0, t] = r; R[:, 1, t] = r
    return O, A, LP, V, R


def train(iters, seed=0, B=256):
    torch.manual_seed(seed); net = Agent(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecStig(B, 1000 + seed * 99999 + it)
        O, A, LP, V, R = rollout(net, env)
        O2, A2, LP2 = O.reshape(2 * B, T, ODIM).reshape(-1, ODIM), A.reshape(2 * B, T), LP.reshape(2 * B, T)
        V2, R2 = V.reshape(2 * B, T), R.reshape(2 * B, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot = torch.from_numpy(O.reshape(2 * B, T, ODIM)); At, LPt = torch.from_numpy(A2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            lg, v = net(Ot.reshape(-1, ODIM)); lg = lg.reshape(2 * B, T, NACT); v = v.reshape(2 * B, T)
            d = torch.distributions.Categorical(logits=lg)
            ratio = torch.exp(d.log_prob(At) - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((v - rett) ** 2).mean() - 0.01 * d.entropy().mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            jf, ov = _eval(net, 7000 + it)
            print(f"  iter {it:>3}: joint food {jf:.1f}/{F}  overlap {ov:.1f}", flush=True)
    return net


def _eval(net, seed, hide=False, own_only=False, B=400):
    env = VecStig(B, seed, hide=hide, own_only=own_only); rollout(net, env, greedy=True)
    return float(env.collected.mean()), env.overlap()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 700
    here = __file__.rsplit("/", 1)[0]
    print(f"STIGMERGY probe: 2 blind agents forage {F} food via a shared MARK trail (no channel, no mutual obs)", flush=True)
    net = train(iters)
    jf, ov = np.mean([_eval(net, 9000 + k) for k in range(3)], 0)                 # combined trail
    jo, oo = np.mean([_eval(net, 9050 + k, own_only=True) for k in range(3)], 0)  # OWN trail only
    jfh, ovh = np.mean([_eval(net, 9100 + k, hide=True) for k in range(3)], 0)    # no trail
    print(f"\nFINAL joint-food/{F}: combined-trail {jf:.1f} (overlap {ov:.1f}) | OWN-trail-only {jo:.1f} | NO-trail {jfh:.1f}", flush=True)
    uses_trail = jf - jfh >= 1.0                       # depends on the mark field at all (stigmergy, broad)
    inter_agent = jf - jo >= 1.0                       # depends on the PARTNER's trail specifically (division of labour)
    print(f"VERDICT: trail-mediated foraging={uses_trail} (gain {jf-jfh:+.1f}); inter-agent stigmergy={inter_agent} "
          f"(combined {jf:.1f} vs own-only {jo:.1f}) -> {'INTER-AGENT' if inter_agent else 'SELF-trail-following only'}", flush=True)
    torch.save(net.state_dict(), f"{here}/stigmergy_net.pt")
    json.dump({"food_combined": round(float(jf), 1), "food_own_only": round(float(jo), 1), "food_none": round(float(jfh), 1),
               "overlap_combined": round(float(ov), 1), "uses_trail": bool(uses_trail), "inter_agent": bool(inter_agent)},
              open(f"{here}/stigmergy_result.json", "w"))
