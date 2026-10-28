"""Asymmetry dial, revised for review (paper §8.1). Supersedes dial2_ppo.py for the paper; dial2 is kept unchanged.

What changed relative to dial2 (each answers a referee point):
  1. FAIR BASELINE. Reliance on the channel is measured against a policy TRAINED FROM SCRATCH WITHOUT the channel
     (signal input absent during training and test), not by scrambling a channel-trained policy at test time.
     The test-time scramble is still reported, as a secondary number, for comparison with dial2.
  2. SEEDS. 10 seeds per (alpha, condition); every seed's result is saved as its own JSON file, so CIs are possible.
  3. FINER ALPHA near 0 (0.1 added).
  4. PER-TARGET-TYPE RELIANCE. The listener's view is redrawn for every target, so an episode mixes shared-target
     and private-target trials. Hits are split by target type and normalized by exposure (steps spent on targets of
     that type), so reliance can be read separately for shared and private targets.
  5. MI FIXED. dial2 paired the speaker's modal symbol over all 24 steps with the FIRST target only, although targets
     are redrawn after every hit. Here every step contributes one (current speaker's symbol, current target cell)
     pair, restricted to steps where the listener cannot see the target; MI is Miller-Madow bias-corrected and
     compared with a shuffle null (symbols permuted across pairs).

Environment, network, PPO settings and training length are identical to dial2 (N=7, K=8, T=24, GRU H=96, 1500 iters,
batch 192, lr 3e-3, clip 0.2, 4 epochs, gamma 0.97, lambda 0.95, entropy 0.02).

Usage (VPS epc-venv, as matthewhmaxwell):
  python dial3_ppo.py run --alpha 0.5 --seed 3 --cond channel     # one training run -> dial3_runs/*.json
  python dial3_ppo.py summarize                                    # aggregate -> dial3_summary.json
"""
import json
import os
import sys

import numpy as np
import torch
import torch.nn as nn

N = 7; K = 8; T = 24; H = 96; NMOVE = 5
ODIM = 2 + 2 + 1 + K
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(int(os.environ.get("DIAL_THREADS", "1")))
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.environ.get("DIAL3_RUNS", os.path.join(HERE, "dial3_runs"))


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.mv = nn.Linear(H, NMOVE); self.sg = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.mv(z), self.sg(z), self.v(z).squeeze(-1), h


class VecDial:
    """mode: 'intact' (partner's last symbol shown), 'scramble' (random symbol shown), 'absent' (no symbol input)."""

    def __init__(self, B, alpha, seed, mode="intact"):
        self.B, self.alpha, self.mode = B, alpha, mode; self.rng = np.random.default_rng(seed); self.reset()

    def _newmask(self, idx):
        self.lis_sees[idx] = self.rng.random(len(idx)) >= self.alpha

    def reset(self):
        B = self.B
        self.ap = self.rng.integers(0, N, size=(B, 2, 2))
        self.tgt = self.rng.integers(0, N, size=(B, 2))
        self.speaker = self.rng.integers(0, 2, size=B)
        self.lis_sees = np.zeros(B, bool); self._newmask(np.arange(B))
        self.lsig = np.zeros((B, 2), int); self.hit = np.zeros(B, int)
        self.hit_by = np.zeros(2); self.steps_by = np.zeros(2)   # index 0 = private target, 1 = shared target
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, 2, ODIM), np.float32)
        for i in range(2):
            j = 1 - i; o = 0
            out[:, i, o:o + 2] = self.ap[:, i] / (N - 1); o += 2
            is_sp = (self.speaker == i)
            sees = is_sp | (~is_sp & self.lis_sees)
            out[sees, i, o:o + 2] = self.tgt[sees] / (N - 1); o += 2
            out[:, i, o] = is_sp.astype(np.float32); o += 1
            if self.mode == "absent":
                continue
            sj = self.rng.integers(0, K, size=B) if self.mode == "scramble" else self.lsig[:, j]
            out[np.arange(B), i, o + sj] = 1.0
        return out

    def step(self, mv, sg):
        B = self.B; rew = np.full(B, -0.01, np.float32)
        typ = self.lis_sees.astype(int)
        np.add.at(self.steps_by, typ, 1)
        for i in range(2): self.ap[:, i] = np.clip(self.ap[:, i] + DIRS[mv[:, i]], 0, N - 1)
        self.lsig = sg.copy()
        both = np.all(self.ap[:, 0] == self.tgt, axis=1) & np.all(self.ap[:, 1] == self.tgt, axis=1)
        hb = np.where(both)[0]
        if len(hb):
            np.add.at(self.hit_by, typ[hb], 1)
            rew[hb] += 1.0; self.hit[hb] += 1
            self.tgt[hb] = self.rng.integers(0, N, size=(len(hb), 2))
            self.speaker[hb] = self.rng.integers(0, 2, size=len(hb))
            self._newmask(hb)
            self.ap[hb] = self.rng.integers(0, N, size=(len(hb), 2, 2))
        return rew

    def hits(self): return float(self.hit.mean())

    def hits_by_type(self):
        """Hits per T steps of exposure to each target type (comparable in scale to hits per episode)."""
        return [float(self.hit_by[k] / self.steps_by[k] * T) if self.steps_by[k] > 0 else float("nan") for k in (0, 1)]


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
    pairs = []  # (speaker symbol, target cell) on steps where the listener cannot see the target
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
        if log:
            priv = ~env.lis_sees
            sym = sg[np.arange(B), env.speaker]
            cell = env.tgt[:, 0] * N + env.tgt[:, 1]
            pairs.append(np.stack([sym[priv], cell[priv]], 1))
        R[:, t] = env.step(mv, sg); obs = env.obs()
    if log:
        return np.concatenate(pairs, 0) if pairs else np.zeros((0, 2), int)
    return O, MV, SG, LP, V, R


def train(alpha, iters, seed, mode, B=192):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecDial(B, alpha, 1000 + seed * 99991 + it, mode=mode)
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
        if it % 300 == 0 or it == iters - 1:
            print(f"  [a={alpha:.2f} s{seed} {mode}] iter {it:>4}", flush=True)
    return net


def evaluate(net, alpha, mode, seeds=(9000, 9001, 9002), B=400):
    tot, by = [], []
    for s in seeds:
        env = VecDial(B, alpha, s, mode=mode); rollout(net, env); tot.append(env.hits()); by.append(env.hits_by_type())
    by = np.array(by, float)
    return float(np.mean(tot)), [float(np.nanmean(by[:, 0])) if not np.all(np.isnan(by[:, 0])) else None,
                                 float(np.nanmean(by[:, 1])) if not np.all(np.isnan(by[:, 1])) else None]


def mi_bits(pairs, nshuf=200, seed=0):
    """Miller-Madow-corrected MI(symbol; cell) in bits, the mean of a shuffle null, and p = P(null >= observed)."""
    if len(pairs) < 200: return None
    def plugin(a, c):
        j = np.zeros((K, N * N)); np.add.at(j, (a, c), 1); n = j.sum()
        p = j / n; pa = p.sum(1, keepdims=True); pc = p.sum(0, keepdims=True); nz = p > 0
        mi = (p[nz] * np.log2(p[nz] / (pa @ pc)[nz])).sum()
        kj, ka, kc = nz.sum(), (pa > 0).sum(), (pc > 0).sum()
        return mi - (kj - ka - kc + 1) / (2 * n * np.log(2))   # Miller-Madow bias correction
    a, c = pairs[:, 0].astype(int), pairs[:, 1].astype(int)
    obs = plugin(a, c); rng = np.random.default_rng(seed)
    null = np.array([plugin(rng.permutation(a), c) for _ in range(nshuf)])
    return {"mi": float(obs), "null_mean": float(null.mean()), "p": float((null >= obs).mean()), "n_pairs": int(len(pairs))}


def run(alpha, seed, cond, iters):
    os.makedirs(RUNS, exist_ok=True)
    out = os.path.join(RUNS, f"{cond}_a{alpha:.2f}_s{seed}.json")
    if os.path.exists(out):
        print("exists:", out); return
    mode = "intact" if cond == "channel" else "absent"
    net = train(alpha, iters, seed, mode)
    row = {"alpha": alpha, "seed": seed, "cond": cond, "iters": iters}
    row["hits"], row["hits_by_type"] = evaluate(net, alpha, mode)           # by_type = [private, shared]
    if cond == "channel":
        row["hits_scrambled"], row["hits_scrambled_by_type"] = evaluate(net, alpha, "scramble")
        env = VecDial(1500, alpha, 9400, mode="intact")
        row["mi"] = mi_bits(rollout(net, env, greedy=True, log=True))
    tmp = out + ".tmp"; json.dump(row, open(tmp, "w"), indent=1); os.replace(tmp, out)
    torch.save(net.state_dict(), out.replace(".json", ".pt"))
    print("saved", out, json.dumps({k: row[k] for k in ("hits", "hits_by_type")}), flush=True)


def ci95(x):
    x = np.array([v for v in x if v is not None and not np.isnan(v)], float)
    if len(x) < 2: return None
    from scipy import stats
    m = x.mean(); h = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return {"mean": float(m), "lo": float(m - h), "hi": float(m + h), "n": int(len(x))}


def summarize():
    rows = [json.load(open(os.path.join(RUNS, f))) for f in sorted(os.listdir(RUNS)) if f.endswith(".json")]
    alphas = sorted({r["alpha"] for r in rows}); summary = []
    for a in alphas:
        val = lambda v: np.nan if v is None else v
        ch = {r["seed"]: r for r in rows if r["alpha"] == a and r["cond"] == "channel"}
        nc = {r["seed"]: r for r in rows if r["alpha"] == a and r["cond"] == "nochannel"}
        nc_mean = np.mean([r["hits"] for r in nc.values()]) if nc else np.nan
        val = lambda v: np.nan if v is None else v
        nc_priv = np.nanmean([val(r["hits_by_type"][0]) for r in nc.values()]) if nc and any(r["hits_by_type"][0] is not None for r in nc.values()) else np.nan
        nc_shared = np.nanmean([val(r["hits_by_type"][1]) for r in nc.values()]) if nc and any(r["hits_by_type"][1] is not None for r in nc.values()) else np.nan
        s = {"alpha": a, "n_channel": len(ch), "n_nochannel": len(nc),
             "hits_channel": ci95([r["hits"] for r in ch.values()]),
             "hits_nochannel": ci95([r["hits"] for r in nc.values()]),
             # fair reliance: each channel seed minus the mean of the independently trained no-channel seeds
             "reliance_fair": ci95([r["hits"] - nc_mean for r in ch.values()]),
             "reliance_fair_private": ci95([val(r["hits_by_type"][0]) - nc_priv for r in ch.values()]),
             "reliance_fair_shared": ci95([val(r["hits_by_type"][1]) - nc_shared for r in ch.values()]),
             "reliance_scramble": ci95([r["hits"] - r["hits_scrambled"] for r in ch.values()]),
             "mi_minus_null": ci95([r["mi"]["mi"] - r["mi"]["null_mean"] for r in ch.values() if r.get("mi")]),
             "mi_p_max": max([r["mi"]["p"] for r in ch.values() if r.get("mi")], default=None)}
        summary.append(s)
    json.dump(summary, open(os.path.join(HERE, "dial3_summary.json"), "w"), indent=1)
    for s in summary:
        f = lambda d: "NA" if not d else f"{d['mean']:+.2f} [{d['lo']:+.2f},{d['hi']:+.2f}] n={d['n']}"
        print(f"a={s['alpha']:.2f} fair {f(s['reliance_fair'])} | private {f(s['reliance_fair_private'])} | "
              f"shared {f(s['reliance_fair_shared'])} | scramble {f(s['reliance_scramble'])} | MI-null {f(s['mi_minus_null'])}")


if __name__ == "__main__":
    a = sys.argv
    if len(a) > 1 and a[1] == "summarize":
        summarize()
    else:
        alpha = float(a[a.index("--alpha") + 1]); seed = int(a[a.index("--seed") + 1]); cond = a[a.index("--cond") + 1]
        iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 1500
        assert cond in ("channel", "nochannel")
        run(alpha, seed, cond, iters)
