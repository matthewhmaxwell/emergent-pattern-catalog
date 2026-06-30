"""Track-1 probe: TOOL USE / instrumental construction. A single agent must reach a goal on the far side
of an impassable GAP (the middle row). It cannot cross directly; it must PUSH a movable BLOCK into a gap
cell to build a BRIDGE, then walk across. The goal is UNREACHABLE without modifying the environment, so the
task forces genuine instrumental construction (use an object as a means to an otherwise-impossible end).

Cleanest diagnostic in the whole arc: the NO-BLOCK control (no block present) makes the goal physically
unreachable -> success must be ~0. WITH a block, success should be high. The gap between them is the
tool-use competency. (Also: reaching the goal requires building the bridge FIRST, then crossing -- an
instrumental sub-goal that is not itself rewarding.)

Run on VPS epc-venv: python tool_ppo.py [--iters 800]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 5; GAP = 2; T = 26; H = 64; ODIM = 2 + 2 + 2 + N + 1; NACT = 5
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        out, h = self.gru(x, h); return self.pi(out), self.v(out).squeeze(-1), h


class VecTool:
    def __init__(self, B, seed, no_block=False):
        self.B, self.no_block = B, no_block; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.ap = np.zeros((B, 2), int); self.goal = np.zeros((B, 2), int); self.bp = np.zeros((B, 2), int)
        self.used = np.zeros(B, bool) | self.no_block; self.filled = np.zeros((B, N), bool); self.done = np.zeros(B, bool)
        for b in range(B):
            self.ap[b] = (self.rng.integers(0, N), self.rng.integers(0, GAP))        # start above the gap
            self.goal[b] = (self.rng.integers(0, N), self.rng.integers(GAP + 1, N))  # goal below the gap
            while True:
                bx, by = self.rng.integers(0, N), self.rng.integers(0, GAP)
                if (bx, by) != tuple(self.ap[b]): break
            self.bp[b] = (bx, by)
        return None

    def passable(self, b, x, y):
        if y != GAP: return True
        return self.filled[b, x]

    def obs(self):
        B = self.B; out = np.zeros((B, ODIM), np.float32); o = 0
        out[:, o:o + 2] = self.ap / (N - 1); o += 2
        out[:, o:o + 2] = np.sign(self.goal - self.ap); o += 2
        out[:, o:o + 2] = np.sign(self.bp - self.ap) * (~self.used)[:, None]; o += 2
        out[:, o:o + N] = self.filled.astype(np.float32); o += N
        out[:, o] = (~self.used).astype(np.float32)
        return out

    def step(self, a):
        B = self.B; rew = np.full(B, -0.02, np.float32)
        for b in range(B):
            if self.done[b]: rew[b] = 0.0; continue
            d = DIRS[a[b]]; tx, ty = self.ap[b, 0] + d[0], self.ap[b, 1] + d[1]
            if not (0 <= tx < N and 0 <= ty < N): continue
            if (not self.used[b]) and (tx, ty) == tuple(self.bp[b]):       # push the block
                bx, by = self.bp[b, 0] + d[0], self.bp[b, 1] + d[1]
                if 0 <= bx < N and 0 <= by < N:
                    if by == GAP and not self.filled[b, bx]:
                        self.filled[b, bx] = True; self.used[b] = True; rew[b] += 0.5  # bridge built
                        self.ap[b] = (tx, ty)
                    elif by != GAP:
                        self.bp[b] = (bx, by); self.ap[b] = (tx, ty)                   # block slides
            elif self.passable(b, tx, ty):
                self.ap[b] = (tx, ty)
            if tuple(self.ap[b]) == tuple(self.goal[b]):
                rew[b] += 1.0; self.done[b] = True
        return rew

    def success(self): return float(self.done.mean())


def gae(R, V, TERM, ACT, gamma=0.98, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(B, np.float32)
        nonterm = 1.0 - TERM[:, t]
        delta = R[:, t] + gamma * nextv * nonterm - V[:, t]
        last = delta + gamma * lam * nonterm * last; adv[:, t] = last * ACT[:, t]
    return adv, adv + V


def rollout(net, env, greedy=False):
    B = env.B; obs = env.obs()
    O = np.zeros((B, T, ODIM), np.float32); A = np.zeros((B, T), int); LP = np.zeros((B, T), np.float32)
    V = np.zeros((B, T), np.float32); R = np.zeros((B, T), np.float32); TERM = np.zeros((B, T), np.float32); ACT = np.zeros((B, T), np.float32)
    h = None
    for t in range(T):
        with torch.no_grad():
            lg, v, h = net(torch.from_numpy(obs)[:, None, :], h)
        d = torch.distributions.Categorical(logits=lg[:, 0]); a = lg[:, 0].argmax(1) if greedy else d.sample()
        O[:, t] = obs; A[:, t] = a.numpy(); LP[:, t] = d.log_prob(a).numpy(); V[:, t] = v[:, 0].numpy()
        ACT[:, t] = (~env.done).astype(np.float32)
        before = env.done.copy(); r = env.step(a.numpy()); R[:, t] = r
        TERM[:, t] = (env.done & ~before).astype(np.float32); obs = env.obs()
    return O, A, LP, V, R, TERM, ACT


def train(iters, seed=0, B=192):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecTool(B, 1000 + seed * 99999 + it)
        O, A, LP, V, R, TERM, ACT = rollout(net, env)
        adv, ret = gae(R, V, TERM, ACT); m = ACT.sum()
        adv = (adv - (adv * ACT).sum() / m) / (np.sqrt(((adv - (adv * ACT).sum() / m) ** 2 * ACT).sum() / m) + 1e-8)
        Ot, At, LPt = torch.from_numpy(O), torch.from_numpy(A), torch.from_numpy(LP)
        advt, rett, actt = torch.from_numpy(adv), torch.from_numpy(ret), torch.from_numpy(ACT)
        for _ in range(4):
            lg, v, _ = net(Ot); d = torch.distributions.Categorical(logits=lg)
            ratio = torch.exp(d.log_prob(At) - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -(torch.min(s1, s2) * actt).sum() / actt.sum() + 0.5 * (((v - rett) ** 2) * actt).sum() / actt.sum() - 0.01 * (d.entropy() * actt).sum() / actt.sum()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: reach-goal {_succ(net, 7000 + it):.2f}", flush=True)
    return net


def _succ(net, seed, no_block=False, B=400):
    env = VecTool(B, seed, no_block=no_block); rollout(net, env, greedy=True); return env.success()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 800
    here = __file__.rsplit("/", 1)[0]
    print("TOOL-USE probe: reach a goal across an impassable gap by pushing a block to build a bridge", flush=True)
    net = train(iters)
    wb = np.mean([_succ(net, 9000 + k) for k in range(3)])
    nb = np.mean([_succ(net, 9100 + k, no_block=True) for k in range(3)])
    print(f"\nFINAL reach-goal: WITH-block {wb:.2f} | NO-block (goal physically unreachable) {nb:.2f}", flush=True)
    tool = wb >= 0.6 and nb <= 0.1
    print(f"VERDICT: {'GENUINE TOOL USE' if tool else 'partial/none'} (with-block {wb:.2f}, no-block {nb:.2f} -> tool gain {wb-nb:+.2f})", flush=True)
    torch.save(net.state_dict(), f"{here}/tool_net.pt")
    json.dump({"with_block": round(float(wb), 2), "no_block": round(float(nb), 2), "tool_use": bool(tool)},
              open(f"{here}/tool_result.json", "w"))
