"""Cost-ladder DIAL (paper §8): the communication niche as ONE continuous causal dial.

Built on the PROVEN rung-3 trainer openhunt4_ppo.py (recurrent GRU; target given as absolute (r,c)
coords; comms verified to emerge -> normal 3.02 / mute 0.30). The only change is a continuous
ASYMMETRY DIAL alpha:

  a private TARGET cell is always shown to the SPEAKER; the LISTENER also sees it with prob (1 - alpha).
    alpha = 0  -> both agents always see the target (SYMMETRIC; a shared referent) -> each just navigates
                  there; the channel is present but should be DECLINED (mute-invariant).
    alpha = 1  -> only the speaker sees it (the original rung-3 asymmetry) -> the speaker must SIGNAL the
                  location; muting the channel COLLAPSES success (communication forced).

Sweep alpha and measure rides-CHANNEL = hits(normal) - hits(mute) per alpha. Prediction (communication
niche, as a causal dial, not 8 separate worlds): rides-CHANNEL ~ 0 at alpha=0 (available-but-declined under
symmetry) and rises monotonically toward alpha=1; MI(speaker symbol; target) rises in lockstep on the
asymmetric episodes.

Run (VPS epc-venv; polite): DIAL_THREADS=4 PYTHONDONTWRITEBYTECODE=1 nice -n 10 \
  /home/matthewhmaxwell/epc-venv/bin/python dial2_ppo.py --iters 1500 --alphas 0,0.5,1 --seeds 1
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 7; K = 8; T = 24; H = 96; NMOVE = 5
ODIM = 2 + 2 + 1 + K              # own pos(2), target(2; zeroed unless this agent sees it), am-speaker(1), other signal(K)
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(int(os.environ.get("DIAL_THREADS", "4")))


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.mv = nn.Linear(H, NMOVE); self.sg = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.mv(z), self.sg(z), self.v(z).squeeze(-1), h


class VecDial:
    def __init__(self, B, alpha, seed, mute=False):
        self.B, self.alpha, self.mute = B, alpha, mute; self.rng = np.random.default_rng(seed); self.reset()

    def _newmask(self, idx):
        self.lis_sees[idx] = self.rng.random(len(idx)) >= self.alpha    # listener also sees target iff True

    def reset(self):
        B = self.B
        self.ap = self.rng.integers(0, N, size=(B, 2, 2))
        self.tgt = self.rng.integers(0, N, size=(B, 2))
        self.speaker = self.rng.integers(0, 2, size=B)
        self.lis_sees = np.zeros(B, bool); self._newmask(np.arange(B))
        self.lsig = np.zeros((B, 2), int); self.hit = np.zeros(B, int)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, 2, ODIM), np.float32)
        for i in range(2):
            j = 1 - i; o = 0
            out[:, i, o:o + 2] = self.ap[:, i] / (N - 1); o += 2
            is_sp = (self.speaker == i)
            sees = is_sp | (~is_sp & self.lis_sees)                      # speaker always; listener iff lis_sees
            out[sees, i, o:o + 2] = self.tgt[sees] / (N - 1); o += 2
            out[:, i, o] = is_sp.astype(np.float32); o += 1             # designated-speaker role flag
            sj = self.rng.integers(0, K, size=B) if self.mute else self.lsig[:, j]
            out[np.arange(B), i, o + sj] = 1.0
        return out

    def step(self, mv, sg):
        B = self.B; rew = np.full(B, -0.01, np.float32)
        for i in range(2): self.ap[:, i] = np.clip(self.ap[:, i] + DIRS[mv[:, i]], 0, N - 1)
        self.lsig = sg.copy()
        both = np.all(self.ap[:, 0] == self.tgt, axis=1) & np.all(self.ap[:, 1] == self.tgt, axis=1)
        hb = np.where(both)[0]
        if len(hb):
            rew[hb] += 1.0; self.hit[hb] += 1
            self.tgt[hb] = self.rng.integers(0, N, size=(len(hb), 2))
            self.speaker[hb] = self.rng.integers(0, 2, size=len(hb))
            self._newmask(hb)
            self.ap[hb] = self.rng.integers(0, N, size=(len(hb), 2, 2))
        return rew

    def hits(self): return float(self.hit.mean())


def gae(R, V, gamma=0.97, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma * nextv - V[:, t]; last = delta + gamma * lam * last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False, log=False):
    B = env.B; obs = env.obs(); h = [None, None]
    O = np.zeros((B, 2, T, ODIM), np.float32); MV = np.zeros((B, 2, T), int); SG = np.zeros((B, 2, T), int)
    LP = np.zeros((B, 2, T), np.float32); V = np.zeros((B, 2, T), np.float32); R = np.zeros((B, T), np.float32)
    spk_sym = np.zeros((B, T), int) if log else None
    tgt0 = env.tgt.copy() if log else None; lis0 = env.lis_sees.copy() if log else None; spk0 = env.speaker.copy() if log else None
    for t in range(T):
        mv = np.zeros((B, 2), int); sg = np.zeros((B, 2), int)
        for i in range(2):
            with torch.no_grad():
                ml, sl, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
            md = torch.distributions.Categorical(logits=ml[:, 0]); sd = torch.distributions.Categorical(logits=sl[:, 0])
            m = ml[:, 0].argmax(1) if greedy else md.sample(); s = sl[:, 0].argmax(1) if greedy else sd.sample()
            O[:, i, t] = obs[:, i]; MV[:, i, t] = m.numpy(); SG[:, i, t] = s.numpy()
            LP[:, i, t] = (md.log_prob(m) + sd.log_prob(s)).numpy(); V[:, i, t] = v[:, 0].numpy()
            mv[:, i] = m.numpy(); sg[:, i] = s.numpy()
        if log:                                             # record the SPEAKER's symbol this step
            spk_sym[:, t] = sg[np.arange(B), spk0]
        R[:, t] = env.step(mv, sg); obs = env.obs()
    if log:
        return spk_sym, tgt0, lis0, spk0
    return O, MV, SG, LP, V, R


def train(alpha, iters, seed=0, B=192):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecDial(B, alpha, 1000 + seed * 99991 + it)
        O, MV, SG, LP, V, R = rollout(net, env)
        Rag = np.repeat(R[:, None, :], 2, axis=1)
        O2 = O.reshape(B * 2, T, ODIM); MV2, SG2 = MV.reshape(B * 2, T), SG.reshape(B * 2, T)
        LP2, V2, R2 = LP.reshape(B * 2, T), V.reshape(B * 2, T), Rag.reshape(B * 2, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, MVt, SGt, LPt = torch.from_numpy(O2), torch.from_numpy(MV2), torch.from_numpy(SG2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            ml, sl, v, _ = net(Ot)
            md = torch.distributions.Categorical(logits=ml); sd = torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(md.log_prob(MVt) + sd.log_prob(SGt) - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((v - rett) ** 2).mean() - 0.02 * (md.entropy() + sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 150 == 0 or it == iters - 1:
            print(f"  [alpha={alpha:.2f} s{seed}] iter {it:>4}: hits/ep {_hits(net, alpha, 7000 + it):.2f}", flush=True)
    return net


def _hits(net, alpha, seed, B=400, mute=False):
    env = VecDial(B, alpha, seed, mute=mute); rollout(net, env, greedy=False); return env.hits()


def _mi(net, alpha, seed, B=1500):
    """MI(speaker modal symbol ; target cell) in bits, on ASYMMETRIC episodes (listener does NOT see target)."""
    env = VecDial(B, alpha, seed); spk_sym, tgt0, lis0, spk0 = rollout(net, env, greedy=True, log=True)
    m = ~lis0
    if m.sum() < 60: return float("nan")
    modal = np.array([np.bincount(spk_sym[b], minlength=K).argmax() for b in np.where(m)[0]])
    cell = (tgt0[m, 0] * N + tgt0[m, 1]).astype(int)
    joint = np.zeros((K, N * N))
    for a, c in zip(modal, cell): joint[a, c] += 1
    joint /= joint.sum(); ps = joint.sum(1, keepdims=True); pg = joint.sum(0, keepdims=True)
    nz = joint > 0
    return float((joint[nz] * np.log2(joint[nz] / (ps @ pg)[nz])).sum())


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 1500
    alphas = [float(x) for x in (a[a.index("--alphas") + 1].split(",") if "--alphas" in a else ["0", "0.5", "1"])]
    seeds = int(a[a.index("--seeds") + 1]) if "--seeds" in a else 1
    here = __file__.rsplit("/", 1)[0]; out = f"{here}/dial2_sweep.json"
    prior = {f'{r["alpha"]:.2f}': r for r in json.load(open(out))} if os.path.exists(out) else {}
    print(f"COMM-NICHE DIAL (recurrent): alphas={alphas} iters={iters} seeds={seeds} (resuming {sorted(prior)})", flush=True)
    results = []
    for alpha in alphas:
        key = f"{alpha:.2f}"
        if key in prior:
            results.append(prior[key]); print(f"alpha={key}: (cached)", flush=True); continue
        norm = mute = mi = 0.0
        for s in range(seeds):
            net = train(alpha, iters, seed=s)
            norm += np.mean([_hits(net, alpha, 9000 + k) for k in range(3)]) / seeds
            mute += np.mean([_hits(net, alpha, 9300 + k, mute=True) for k in range(3)]) / seeds
            mi += _mi(net, alpha, 9400) / seeds
            torch.save(net.state_dict(), f"{here}/dial2_a{key}_s{s}.pt")
        rides = norm - mute
        row = {"alpha": alpha, "hits_normal": round(norm, 3), "hits_mute": round(mute, 3),
               "rides_channel": round(rides, 3), "mi_sym_tgt_bits": round(mi, 3)}
        results.append(row); json.dump(results, open(out, "w"), indent=1)
        print(f"alpha={key}: normal {norm:.2f} | mute {mute:.2f} | RIDES-CHANNEL {rides:+.2f} | MI {mi:.2f} bits", flush=True)
    json.dump(results, open(out, "w"), indent=1); print("DONE saved dial2_sweep.json", flush=True)
