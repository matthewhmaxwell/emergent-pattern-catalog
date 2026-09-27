"""P1 MISS -> novelty protocol step 2: the EPC debunking battery for P1 (osf.io/up7m9, PREREG §4).
Criteria below are fixed BEFORE this script is run on any trained model.

Claim to debunk: exchangeable agents with zero environmental information asymmetry recruit the talk channel to
break symmetry, by transmitting ENDOGENOUS RANDOMNESS (a jointly-controlled lottery), not world information.

Static audits (code-level, no trained model):
  S1 exchangeability: agent 0's observation under symbols (a,b) == agent 1's observation under (b,a); t=0 obs identical.
  S2 reward leak: reward depends only on whether the two actions differ (never on symbols).
Behavioral audits (per trained model, 4000 episodes each):
  B1 symbol randomization: talk-symbol entropy (bits; uniform = 3.0).
  B2 conditional success: S | symbols distinct  (lottery => ~1.0)   and   S | symbols tied (=> ~0.5).
  B3 provably-symmetric control (LEAK TEST): force BOTH agents to the same symbol k, for every k => S must be ~0.5.
  B4 antisymmetry: force every ordered pair (a,b), a!=b => P(actions differ) (lottery => ~1.0).
  B5 randomness is the resource: greedy (argmax) talk keeps the channel but removes randomness => ~0.5.
  B6 mute: partner symbol scrambled => ~0.5.
PASS (mechanism = jointly-controlled lottery, no leak) iff, per model: B2 distinct >= 0.90 and |B2 tied - 0.5| <= 0.10;
max_k |B3_k - 0.5| <= 0.10; B4 mean >= 0.90; B5 <= 0.60; B6 <= 0.60. Report every model regardless.
"""
import os, sys, json, glob, numpy as np, torch
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prereg_envs as P
K = P.K


def static_audits():
    B = 64; rng = np.random.default_rng(0); a = rng.integers(0, K, B); b = rng.integers(0, K, B)
    e1 = P.P1Env(B, 1); e1.t = 1; e1.sym = np.stack([a, b], 1)
    e2 = P.P1Env(B, 1); e2.t = 1; e2.sym = np.stack([b, a], 1)
    s1 = bool(np.array_equal(e1.obs()[:, 0], e2.obs()[:, 1]) and np.array_equal(P.P1Env(B, 2).obs()[:, 0], P.P1Env(B, 2).obs()[:, 1]))
    e = P.P1Env(4, 3); e.step(np.zeros((4, 2), int), np.array([[0, 0], [0, 7], [3, 3], [5, 1]]))
    r = e.step(np.array([[0, 1], [0, 0], [1, 0], [1, 1]]), np.zeros((4, 2), int))
    s2 = r.tolist() == [1.0, 0.0, 1.0, 0.0]
    return {"S1_exchangeable_obs": s1, "S2_reward_depends_only_on_actions_differing": s2}


def episode(net, B, seed, force=None, greedy_talk=False, mute=False):
    env = P.P1Env(B, seed, mute=mute); obs = env.obs(); h = [None, None]; sg = np.zeros((B, 2), int)
    for i in range(2):
        with torch.no_grad(): _, sl, _, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
        sg[:, i] = (sl[:, 0].argmax(1) if greedy_talk else torch.distributions.Categorical(logits=sl[:, 0]).sample()).numpy()
    if force is not None: sg[:, 0], sg[:, 1] = force
    env.step(np.zeros((B, 2), int), sg); obs = env.obs(); mv = np.zeros((B, 2), int)
    for i in range(2):
        with torch.no_grad(): ml, _, _, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
        mv[:, i] = torch.distributions.Categorical(logits=ml[:, 0]).sample().numpy()
    env.step(mv, np.zeros((B, 2), int)); return sg, env.succ


def audit(path, B=4000):
    net = P.Pol(P.P1Env.odim, P.P1Env.nmove, P.P1Env.k); net.load_state_dict(torch.load(path)); net.eval()
    torch.manual_seed(7)
    sg, succ = episode(net, B, 11)
    p = np.bincount(sg.ravel(), minlength=K) / sg.size; H = float(-(p[p > 0] * np.log2(p[p > 0])).sum())
    tie = sg[:, 0] == sg[:, 1]
    b3 = {k: float(episode(net, 1000, 20 + k, force=(k, k))[1].mean()) for k in range(K)}
    b4 = [float(episode(net, 300, 100 + a * K + b, force=(a, b))[1].mean()) for a in range(K) for b in range(K) if a != b]
    r = {"B1_symbol_entropy_bits": round(H, 3), "S_normal": round(float(succ.mean()), 3),
         "B2_S_given_distinct": round(float(succ[~tie].mean()), 3), "B2_S_given_tied": round(float(succ[tie].mean()), 3) if tie.any() else None,
         "B2_tie_rate": round(float(tie.mean()), 3),
         "B3_forced_identical_by_k": {k: round(v, 3) for k, v in b3.items()},
         "B4_antisymmetry_mean": round(float(np.mean(b4)), 3), "B4_min": round(float(np.min(b4)), 3),
         "B5_greedy_talk": round(float(episode(net, B, 12, greedy_talk=True)[1].mean()), 3),
         "B6_mute": round(float(episode(net, B, 13, mute=True)[1].mean()), 3)}
    r["PASS"] = bool(r["B2_S_given_distinct"] >= 0.90 and r["B2_S_given_tied"] is not None and abs(r["B2_S_given_tied"] - 0.5) <= 0.10
                     and max(abs(v - 0.5) for v in b3.values()) <= 0.10 and r["B4_antisymmetry_mean"] >= 0.90
                     and r["B5_greedy_talk"] <= 0.60 and r["B6_mute"] <= 0.60)
    return r


if __name__ == "__main__":
    out = {"static": static_audits(), "models": {}}
    paths = sorted(glob.glob(os.path.join(HERE, "prereg_runs", "P1_full_s*.pt")) +
                   glob.glob(os.path.join(HERE, "replication_P1", "prereg_runs", "P1_full_s*.pt")))
    for pth in paths:
        tag = os.path.basename(pth)[:-3]; out["models"][tag] = audit(pth); print(tag, json.dumps(out["models"][tag]), flush=True)
    out["all_models_PASS"] = all(m["PASS"] for m in out["models"].values())
    json.dump(out, open(os.path.join(HERE, "P1_audit.json"), "w"), indent=1)
    print("STATIC:", out["static"]); print("ALL MODELS PASS:", out["all_models_PASS"])
