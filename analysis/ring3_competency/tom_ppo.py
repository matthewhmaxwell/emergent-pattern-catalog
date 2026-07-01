"""SCALE step 1 -- close the emulation shortcut from #19 to force PREDICTIVE INTENTION-READING (a component of
theory-of-mind). A scripted MOVER starts at a RANDOM location and walks toward its private goal g* for K steps
then FREEZES partway. The geometry is rigged so the goal NEAREST the mover's frozen position is a DECOY (with
a RANDOM offset direction, carrying no heading info), while g* is FAR along the mover's heading and a FILLER is
equally far in a DIFFERENT direction. So neither proximity (->decoy) nor distance-from-anywhere separates g*
from the filler -- ONLY the mover's observed VELOCITY does, and (start being random) velocity is recoverable
ONLY by integrating the mover's motion over time.

To isolate the intention-reading competency from navigation, the helper is a DISEMBODIED PREDICTOR: it watches
the mover for T steps and outputs a goal CHOICE (which goal is g*); it is scored on its final-step choice.

Discriminators:
  - RECURRENT predictor (GRU, integrates the mover's motion=velocity) vs MEMORYLESS (MLP, sees only the current
    frame -> after the mover freezes it has no motion, and with random start a single frame gives no heading).
    Recurrent >> memoryless == the competency genuinely requires integrating the trajectory (velocity/intention).
  - mover-FROZEN (never moves) -> no velocity to read -> chance.
  - proximity-baseline error: nearest-goal-to-stop != g* (high by construction) -- proximity/emulation fails.

Run on VPS epc-venv: python tom_ppo.py [--iters 800]
"""
import numpy as np, sys, json, torch, torch.nn as nn

N = 17; G = 3; T = 7; K = 3; RHO = 5; H = 96        # RHO = goal radius from the mover's stop; NACT = G (goal choice)
NACT = G
ODIM = 2 + G * 2                                     # mover xy(2), G goal xy(G*2)
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
AX = np.array([(1, 0), (-1, 0), (0, 1), (0, -1)])   # axis headings -> exact integer motion, mover lands ON the stop
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self, recurrent):
        super().__init__(); self.recurrent = recurrent
        if recurrent: self.gru = nn.GRU(ODIM, H, batch_first=True)
        else: self.mlp = nn.Sequential(nn.Linear(ODIM, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        if self.recurrent:
            z, h = self.gru(x, h); return self.pi(z), self.v(z).squeeze(-1), h
        z = self.mlp(x); return self.pi(z), self.v(z).squeeze(-1), None


def clip(p): return np.clip(p, 0, N - 1)


class VecToM:
    def __init__(self, B, seed, mode="normal"):
        self.B, self.mode = B, mode; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.gpos = np.zeros((B, G, 2), int); self.stop = np.zeros((B, 2), int)
        self.mpos = np.zeros((B, 2), int); self.gstar = np.zeros(B, int); self.sv = np.zeros((B, 2), int)
        for b in range(B):
            stop = self.rng.integers(RHO, N - RHO, size=2)                 # anchor: ALL goals sit at distance RHO from here
            dirs = self.rng.permutation(4)[:G]                             # G DISTINCT axis directions
            sv = AX[dirs[0]]                                               # heading = g* direction
            goals = [(stop + RHO * AX[i]).astype(int) for i in dirs]        # EQUIDISTANT (exact integer) -> single frame is SYMMETRIC
            order = self.rng.permutation(G); goals = [goals[i] for i in order]
            self.gpos[b] = goals; self.gstar[b] = int(np.where(order == 0)[0][0])
            self.sv[b] = sv; self.stop[b] = stop
            self.mpos[b] = stop if self.mode == "frozen" else (stop - K * sv)   # frozen: sits AT the stop (same final frame, NO motion)
        self.t = 0; self.done = np.zeros(B, bool); self.win = np.zeros(B, bool)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, ODIM), np.float32); o = 0
        out[:, o:o+2] = self.mpos/(N-1); o += 2                             # frozen vs normal already baked into mpos
        for g in range(G): out[:, o:o+2] = self.gpos[:, g]/(N-1); o += 2
        return out

    def step(self, a):
        B = self.B; rew = np.zeros(B, np.float32)
        # mover steps EXACTLY along its axis heading for K steps -> lands ON the stop (symmetric final frame), then freezes
        if self.mode != "frozen" and self.t < K:
            self.mpos = clip(self.mpos + self.sv).astype(int)
        if self.t >= K:                                                    # DENSE: score every post-motion step
            correct = (a == self.gstar); rew = correct.astype(np.float32)
            if self.t == T - 1: self.win = correct; self.done[:] = True     # final-step choice = the reported accuracy
        self.t += 1
        return rew

    def acc(self): return float(self.win.mean())

    def proximity_err(self):
        errs = 0
        for b in range(self.B):
            s = self.mpos[b]; dists = [abs(s[0]-self.gpos[b,g,0])+abs(s[1]-self.gpos[b,g,1]) for g in range(G)]
            if int(np.argmin(dists)) != self.gstar[b]: errs += 1
        return errs / self.B


def gae(R, V, TERM, gamma=0.98, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32); nt = 1.0 - TERM[:, t]
        delta = R[:, t] + gamma*nextv*nt - V[:, t]; last = delta + gamma*lam*nt*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B = env.B; obs = env.obs(); h = None
    O = np.zeros((B, T, ODIM), np.float32); A = np.zeros((B, T), int); LP = np.zeros((B, T), np.float32)
    V = np.zeros((B, T), np.float32); Rw = np.zeros((B, T), np.float32); TERM = np.zeros((B, T), np.float32)
    for t in range(T):
        with torch.no_grad():
            if net.recurrent: lg, v, h = net(torch.from_numpy(obs)[:, None, :], h); lg, v = lg[:, 0], v[:, 0]
            else: lg, v, _ = net(torch.from_numpy(obs))
        d = torch.distributions.Categorical(logits=lg); a = lg.argmax(1) if greedy else d.sample()
        O[:, t] = obs; A[:, t] = a.numpy(); LP[:, t] = d.log_prob(a).numpy(); V[:, t] = v.numpy()
        Rw[:, t] = env.step(a.numpy()); TERM[:, t] = env.done.astype(np.float32); obs = env.obs()
    return O, A, LP, V, Rw, TERM


def train(iters, recurrent, seed=0, B=256):
    torch.manual_seed(seed); net = Pol(recurrent); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    tag = "recurrent" if recurrent else "memoryless"
    for it in range(iters):
        env = VecToM(B, 1000 + it)
        O, A, LP, V, Rw, TERM = rollout(net, env)
        adv, ret = gae(Rw, V, TERM); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, At, LPt = torch.from_numpy(O), torch.from_numpy(A), torch.from_numpy(LP)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            if recurrent: lg, v, _ = net(Ot)
            else: lg, v, _ = net(Ot.reshape(-1, ODIM)); lg = lg.reshape(B, T, NACT); v = v.reshape(B, T)
            dd = torch.distributions.Categorical(logits=lg)
            ratio = torch.exp(dd.log_prob(At) - LPt); s1 = ratio*advt; s2 = torch.clamp(ratio, 0.8, 1.2)*advt
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.01*dd.entropy().mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  [{tag}] iter {it:>3}: acc {_acc(net, 7000+it,'normal'):.2f}", flush=True)
    return net


def _acc(net, seed, mode, B=500):
    env = VecToM(B, seed, mode); rollout(net, env, greedy=True); return env.acc()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 800
    here = __file__.rsplit("/", 1)[0]; chance = 1/G
    print(f"ToM probe (choice): 3 EQUIDISTANT goals -> the single frame is SYMMETRIC; ONLY the mover's observed motion picks g* (chance {chance:.2f})", flush=True)
    print("== training RECURRENT predictor ==", flush=True); rec = train(iters, True)
    print("== training MEMORYLESS predictor ==", flush=True); mem = train(iters, False)
    rn = np.mean([_acc(rec, 9000+k, "normal") for k in range(3)])
    rf = np.mean([_acc(rec, 9100+k, "frozen") for k in range(3)])
    mn = np.mean([_acc(mem, 9000+k, "normal") for k in range(3)])
    print(f"\nFINAL: recurrent {rn:.2f} | recurrent mover-frozen {rf:.2f} | memoryless {mn:.2f} (chance {chance:.2f})", flush=True)
    genuine = rn >= 0.7 and (rn - mn) >= 0.15 and rf <= chance + 0.12
    print(f"VERDICT: {'GENUINE PREDICTIVE INTENTION-READING (velocity/ToM)' if genuine else 'partial/none'} "
          f"(needs-memory={rn-mn:+.2f}; needs-mover={rn-rf:+.2f})", flush=True)
    torch.save(rec.state_dict(), f"{here}/tom_net.pt")
    json.dump({"recurrent": round(float(rn),2), "recurrent_frozen": round(float(rf),2), "memoryless": round(float(mn),2),
               "chance": round(float(chance),2), "genuine": bool(genuine)}, open(f"{here}/tom_result.json","w"))
