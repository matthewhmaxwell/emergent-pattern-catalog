"""OPEN-ENDED N-agent hunt -- RUNG 4 (climb past single-medium: DISTRIBUTED-INFORMATION INTEGRATION, 3 agents).

The cascade's communication niche (rung 3) was ONE speaker who knew the whole target. Rung 4 breaks the target
across TWO speakers so no single signal suffices -- the minimal sufficient mechanism is forced to be a
COMPOSITION (the multi-agent analog of the sufficiency-pressure / compositionality finding, #14).

World: target cell (r, c) on an N x N grid. Agent 0 (ROW-knower) sees ONLY r; agent 1 (COL-knower) sees ONLY
c; agent 2 (SEEKER) sees neither -- only its own position + both informants' last signals. The seeker must
reach (r, c). To succeed, agent 0 must encode r, agent 1 must encode c, and the seeker must INTEGRATE the two
partial signals into a cell and navigate there. Informants are stationary; only the seeker moves.

Fingerprint (integration is genuine iff BOTH partial channels are required):
  - MUTE-0 (scramble the row-knower's signal) -> collapse (seeker loses the row).
  - MUTE-1 (scramble the col-knower's signal) -> collapse (seeker loses the column).
If EITHER-mute collapses it => distributed two-source integration (neither speaker alone suffices). Expected
KNOWN (distributed source coding / multi-party emergent comms); the climb asks whether a composed multi-party
protocol produces anything off-catalog, and hardens the cost-hierarchy law with a >2-agent, composed rung.

Run on VPS epc-venv: python openhunt5_ppo.py [--iters 1500]
"""
import numpy as np, sys, json, torch, torch.nn as nn

N = 7; N_AG = 3; K = 8; T = 22; H = 96; NMOVE = 5
ODIM = 3 + 2 + 1 + 2 * K            # role one-hot(3), own pos(2), my partial info(1), the 2 others' last signals(2K)
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.mv = nn.Linear(H, NMOVE); self.sg = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.mv(z), self.sg(z), self.v(z).squeeze(-1), h


class VecDist:
    def __init__(self, B, seed, mute0=False, mute1=False):
        self.B, self.mute0, self.mute1 = B, mute0, mute1
        self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.seek = self.rng.integers(0, N, size=(B, 2))          # seeker position
        self.tgt = self.rng.integers(0, N, size=(B, 2))           # (row, col) target
        self.lsig = np.zeros((B, N_AG), int); self.hit = np.zeros(B, int)
        return self.obs()

    def heard(self, i):
        B = self.B; out = np.zeros((B, 2 * K), np.float32); slot = 0
        for j in range(N_AG):
            if j == i: continue
            s = self.lsig[:, j].copy()
            if j == 0 and self.mute0: s = self.rng.integers(0, K, size=B)
            if j == 1 and self.mute1: s = self.rng.integers(0, K, size=B)
            out[np.arange(B), slot * K + s] = 1.0; slot += 1
        return out

    def obs(self):
        B = self.B; out = np.zeros((B, N_AG, ODIM), np.float32)
        for i in range(N_AG):
            o = 0; out[:, i, o + i] = 1.0; o += 3                 # role one-hot
            if i == 2: out[:, i, o:o+2] = self.seek / (N - 1)     # only the seeker's position matters
            o += 2
            if i == 0: out[:, i, o] = self.tgt[:, 0] / (N - 1)    # row-knower sees r
            elif i == 1: out[:, i, o] = self.tgt[:, 1] / (N - 1)  # col-knower sees c
            o += 1
            out[:, i, o:o + 2*K] = self.heard(i)
        return out

    def step(self, mv, sg):
        B = self.B; rew = np.full(B, -0.01, np.float32)
        self.seek = np.clip(self.seek + DIRS[mv[:, 2]], 0, N - 1)  # only the seeker moves
        self.lsig = sg.copy()
        at = np.all(self.seek == self.tgt, axis=1)
        for b in range(B):
            if at[b]:
                rew[b] += 1.0; self.hit[b] += 1
                self.tgt[b] = self.rng.integers(0, N, size=2); self.seek[b] = self.rng.integers(0, N, size=2)
        return rew

    def hits(self): return float(self.hit.mean())


def gae(R, V, gamma=0.97, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma*nextv - V[:, t]; last = delta + gamma*lam*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B, A = env.B, N_AG; obs = env.obs(); h = [None]*A
    O = np.zeros((B, A, T, ODIM), np.float32); MV = np.zeros((B, A, T), int); SG = np.zeros((B, A, T), int)
    LP = np.zeros((B, A, T), np.float32); V = np.zeros((B, A, T), np.float32); R = np.zeros((B, T), np.float32)
    for t in range(T):
        mv = np.zeros((B, A), int); sg = np.zeros((B, A), int)
        for i in range(A):
            with torch.no_grad():
                ml, sl, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
            md = torch.distributions.Categorical(logits=ml[:, 0]); sd = torch.distributions.Categorical(logits=sl[:, 0])
            m = ml[:, 0].argmax(1) if greedy else md.sample(); s = sl[:, 0].argmax(1) if greedy else sd.sample()
            O[:, i, t] = obs[:, i]; MV[:, i, t] = m.numpy(); SG[:, i, t] = s.numpy()
            LP[:, i, t] = (md.log_prob(m) + sd.log_prob(s)).numpy(); V[:, i, t] = v[:, 0].numpy()
            mv[:, i] = m.numpy(); sg[:, i] = s.numpy()
        R[:, t] = env.step(mv, sg); obs = env.obs()
    return O, MV, SG, LP, V, R


def train(iters, seed=0, B=192):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecDist(B, 1000 + it)
        O, MV, SG, LP, V, R = rollout(net, env)
        Rag = np.repeat(R[:, None, :], N_AG, axis=1)
        O2 = O.reshape(B*N_AG, T, ODIM); MV2, SG2 = MV.reshape(B*N_AG, T), SG.reshape(B*N_AG, T)
        LP2, V2, R2 = LP.reshape(B*N_AG, T), V.reshape(B*N_AG, T), Rag.reshape(B*N_AG, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, MVt, SGt, LPt = torch.from_numpy(O2), torch.from_numpy(MV2), torch.from_numpy(SG2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            ml, sl, v, _ = net(Ot)
            md = torch.distributions.Categorical(logits=ml); sd = torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(md.log_prob(MVt) + sd.log_prob(SGt) - LPt)
            s1 = ratio*advt; s2 = torch.clamp(ratio, 0.8, 1.2)*advt
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.02*(md.entropy()+sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: hits/episode {_hits(net, 7000+it):.2f}", flush=True)
    return net


def _hits(net, seed, B=400, **kw):
    env = VecDist(B, seed, **kw); rollout(net, env, greedy=False); return env.hits()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 1500
    here = __file__.rsplit("/", 1)[0]
    print(f"OPEN-ENDED HUNT rung 4: 3-agent DISTRIBUTED-INFO integration (row-knower + col-knower -> seeker)", flush=True)
    net = train(iters)
    norm = np.mean([_hits(net, 9000+k) for k in range(3)])
    m0 = np.mean([_hits(net, 9300+k, mute0=True) for k in range(3)])
    m1 = np.mean([_hits(net, 9400+k, mute1=True) for k in range(3)])
    print(f"\nFINGERPRINT (hits/episode): normal {norm:.2f} | mute-ROW {m0:.2f} | mute-COL {m1:.2f}", flush=True)
    needs_row = norm - m0; needs_col = norm - m1
    integrates = norm >= 1.0 and needs_row >= 0.4 and needs_col >= 0.4
    print(f"  needs-ROW-signal = {needs_row:+.2f}; needs-COL-signal = {needs_col:+.2f}", flush=True)
    print(f"RUNG-4 MECHANISM: {'DISTRIBUTED-INFO INTEGRATION (composed 2-source protocol; BOTH partial channels required)' if integrates else 'degenerate/single-source or weak'}", flush=True)
    torch.save(net.state_dict(), f"{here}/openhunt5_net.pt")
    json.dump({"normal": round(float(norm),2), "mute_row": round(float(m0),2), "mute_col": round(float(m1),2),
               "integrates": bool(integrates)}, open(f"{here}/openhunt5_result.json","w"))
