"""Track-1 probe: ANTICIPATION / forward-modeling. A single agent must INTERCEPT a target that moves on a
PREDICTABLE trajectory (constant velocity, bouncing off walls). Both move 1 cell/step, so a reactive chaser
(always head to the target's CURRENT cell) stays a constant distance behind and never catches it. Success
requires inferring the target's velocity from its motion and moving to where it WILL be -- an internal
forward model. The agent observes only positions (its own + the target's current cell), NOT the velocity,
so it must infer velocity over time (recurrent state) and extrapolate.

Diagnostics: (a) UNPREDICTABLE target (velocity re-randomized each step) -> anticipation is impossible, so
genuine forward-modeling collapses to reactive-chase level; (b) a MEMORYLESS baseline cannot infer velocity
-> reactive only. anticipation = predictable-success - unpredictable-success, and recurrent - memoryless.
New state-content type for the ladder: a *predicted future* (extrapolation), distinct from a held cue,
count, or belief-over-rule.

Run on VPS epc-venv: python anticipate_ppo.py [--iters 600]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 9; T = 22; H = 64; OBS = 4; NACT = 5
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
TVEL = np.array([(1, 0), (-1, 0), (0, 1), (0, -1)])       # speed 1; RENDEZVOUS demands predicting the exact final cell
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self, recurrent=True):
        super().__init__(); self.recurrent = recurrent
        if recurrent: self.gru = nn.GRU(OBS, H, batch_first=True)
        else: self.mlp = nn.Sequential(nn.Linear(OBS, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        if self.recurrent:
            out, h = self.gru(x, h); z = out
        else:
            z = self.mlp(x); h = None
        return self.pi(z), self.v(z).squeeze(-1), h


class VecChase:
    def __init__(self, B, seed, predictable=True):
        self.B, self.predictable = B, predictable; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.ap = self.rng.integers(0, N, size=(B, 2)); self.tp = self.rng.integers(0, N, size=(B, 2))
        self.tv = TVEL[self.rng.integers(0, 4, size=B)].copy()
        self.done = np.zeros(B, bool); self.hit = np.zeros(B, bool)
        return self.obs()

    def obs(self):
        return np.concatenate([self.ap / (N - 1), self.tp / (N - 1)], 1).astype(np.float32)

    def step(self, a):
        B = self.B
        npos = self.ap + DIRS[a]
        inb = (npos[:, 0] >= 0) & (npos[:, 0] < N) & (npos[:, 1] >= 0) & (npos[:, 1] < N)
        self.ap = np.where(inb[:, None], npos, self.ap)
        if not self.predictable:
            self.tv = TVEL[self.rng.integers(0, 4, size=B)].copy()
        for ax in (0, 1):
            out = (self.tp[:, ax] + self.tv[:, ax] < 0) | (self.tp[:, ax] + self.tv[:, ax] >= N)
            self.tv[out, ax] *= -1
        self.tp = np.clip(self.tp + self.tv, 0, N - 1)
        return np.full(B, -0.01, np.float32), np.zeros(B, bool), np.ones(B, bool)

    def terminal(self):                                    # rendezvous reward at the final step only
        d = np.abs(self.ap - self.tp).sum(1)
        return (1.0 * (d == 0) + 0.5 * np.maximum(0.0, 1 - np.minimum(d, 4) / 4) * (d > 0)).astype(np.float32)

    def success(self): return float((np.abs(self.ap - self.tp).sum(1) == 0).mean())


def gae(R, V, TERM, ACT, gamma=0.97, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(B, np.float32)
        nonterm = 1.0 - TERM[:, t]
        delta = R[:, t] + gamma * nextv * nonterm - V[:, t]
        last = delta + gamma * lam * nonterm * last; adv[:, t] = last * ACT[:, t]
    return adv, adv + V


def rollout(net, env, greedy=False):
    B = env.B; obs = env.reset()
    O = np.zeros((B, T, OBS), np.float32); A = np.zeros((B, T), int); LP = np.zeros((B, T), np.float32)
    V = np.zeros((B, T), np.float32); R = np.zeros((B, T), np.float32); TERM = np.zeros((B, T), np.float32); ACT = np.zeros((B, T), np.float32)
    h = None
    for t in range(T):
        with torch.no_grad():
            lg, v, h = net(torch.from_numpy(obs)[:, None, :], h)
        d = torch.distributions.Categorical(logits=lg[:, 0]); a = lg[:, 0].argmax(1) if greedy else d.sample()
        O[:, t] = obs; A[:, t] = a.numpy(); LP[:, t] = d.log_prob(a).numpy(); V[:, t] = v[:, 0].numpy()
        ACT[:, t] = 1.0
        r, _, _ = env.step(a.numpy()); R[:, t] = r; obs = env.obs()
    R[:, -1] += env.terminal()                             # rendezvous reward at the final step
    return O, A, LP, V, R, TERM, ACT


def train(iters, recurrent=True, predictable=True, seed=0, B=320):
    torch.manual_seed(seed); net = Pol(recurrent); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecChase(B, 1000 + seed * 99999 + it, predictable)
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
            ev = np.mean([_succ(net, 7000 + it + k, predictable) for k in range(2)])
            print(f"  iter {it:>3}: intercept {ev:.2f}", flush=True)
    return net


def _succ(net, seed, predictable, B=500):
    env = VecChase(B, seed, predictable); rollout(net, env, greedy=True); return env.success()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 600
    here = __file__.rsplit("/", 1)[0]
    print("ANTICIPATION probe: intercept a constant-velocity bouncing target (positions only -> infer velocity)", flush=True)
    rec = train(iters, recurrent=True, predictable=True)
    mem = train(iters, recurrent=False, predictable=True)
    rp = np.mean([_succ(rec, 9000 + k, True) for k in range(3)])      # recurrent on predictable
    ru = np.mean([_succ(rec, 9100 + k, False) for k in range(3)])     # recurrent on UNPREDICTABLE
    mp = np.mean([_succ(mem, 9200 + k, True) for k in range(3)])      # memoryless on predictable
    print(f"\nFINAL intercept: recurrent/predictable {rp:.2f} | recurrent/UNPREDICTABLE {ru:.2f} | memoryless/predictable {mp:.2f}", flush=True)
    antic = (rp - ru >= 0.2) and (rp - mp >= 0.2)
    print(f"VERDICT: {'GENUINE ANTICIPATION' if rp >= 0.6 and antic else 'partial/none'} "
          f"(predictability gain {rp-ru:+.2f}, memory gain {rp-mp:+.2f})", flush=True)
    torch.save(rec.state_dict(), f"{here}/anticipate_net.pt")
    json.dump({"rec_predictable": round(float(rp), 2), "rec_unpredictable": round(float(ru), 2),
               "mem_predictable": round(float(mp), 2), "anticipation": bool(antic)}, open(f"{here}/anticipate_result.json", "w"))
