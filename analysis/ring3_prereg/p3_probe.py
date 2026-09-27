"""P3 follow-up (post-hoc, exploratory; osf.io/up7m9, DEVIATIONS D3): does the rho=1 agent STORE the t=0 cue at all?
Criteria below are fixed BEFORE this script is run on any trained model (committed to git first).

D3 found: rho=1 agents choose by the landmark (conflict landmark share 0.86-0.88) yet need their recurrent state to
navigate. Open question: is the cue absent from memory, or stored and then overridden by the landmark?

Design (per trained model, 2000 episodes, episode seeds 960000+, disjoint from all earlier evaluations):
  - cue index c and landmark index l are drawn INDEPENDENTLY (half congruent, half conflict); landmark always shown.
  - the agent is held in place ("stay" action) for steps t = 0..3, so its trajectory is identical whatever c and l
    are; the hidden state after each step (h_0..h_3) is recorded. Only h_1..h_3 are scored: the cue is visible at t=0
    only, so from t=1 on any cue information in h must come from memory.
  - linear probe (L2 logistic regression, C=1.0, standardized h, 5-fold stratified CV accuracy) decodes c and l.
Models: full rho=1 seeds 1-6 (test), full rho=0.5 seeds 1-3 (positive control: cue-first agents, D3 L 0.10-0.15).
Controls (INVALID if either fails, seed median): rho=0.5 cue decodability at t=3 >= 0.90; rho=1 landmark decodability
at t=3 >= 0.90 (the landmark is in the input, so this checks the probe works on these hidden states).
Scored (full rho=1 models, cue decodability at t=3; >= 2 of 3 original seeds AND the seed median; seeds 4-6 judged
separately with the same rule):
  CUE NOT STORED if <= 0.60;  CUE STORED (and overridden by the landmark) if >= 0.90;  WEAKLY STORED otherwise.
Reported for every model and t = 1..3: cue and landmark decodability.
"""
import os, re, sys, json, glob, numpy as np, torch
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prereg_envs as P

NEP, SEED0, HOLD = 2000, 960000, 4
PAT = re.compile(r"P3_full_rho([0-9.]+)_s(\d+)\.pt$")


def hidden_states(net, seed=SEED0):
    env = P.P3Env(NEP, seed, rho=1.0); rng = np.random.default_rng(seed)
    cue = env.correct.copy(); lmk = rng.integers(0, 2, size=NEP); r = np.arange(NEP); h = None; H = []
    for t in range(HOLD):
        o = env.obs(); o[:, 0, 8:10] = 0.0; o[r, 0, 8 + lmk] = 1.0          # independent landmark, always shown
        with torch.no_grad(): _, _, _, h = net(torch.from_numpy(o[:, 0])[:, None, :], h)
        H.append(h[0].numpy().copy()); env.step(np.zeros((NEP, 1), int), np.zeros((NEP, 1), int))   # stay
    return H, cue, lmk


def decode(X, y):
    clf = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=2000))
    return float(cross_val_score(clf, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=0)).mean())


def rule(vals, pred):
    return len(vals) >= 2 and sum(pred(v) for v in vals) >= 2 and pred(float(np.median(vals)))


def verdict(vals):
    return ("CUE NOT STORED" if rule(vals, lambda x: x <= 0.60) else
            "CUE STORED (overridden by landmark)" if rule(vals, lambda x: x >= 0.90) else "WEAKLY STORED")


if __name__ == "__main__":
    if "--smoke" in sys.argv:                                  # untrained network: execution check only
        torch.manual_seed(0); net = P.Pol(P.P3Env.odim, P.P3Env.nmove, P.P3Env.k); net.eval()
        H, c, l = hidden_states(net); print({t: (decode(H[t], c), decode(H[t], l)) for t in (1, 3)}); print("SMOKE OK"); sys.exit(0)
    paths = sorted(glob.glob(os.path.join(HERE, "prereg_runs", "P3_full_rho*_s*.pt")) +
                   glob.glob(os.path.join(HERE, "replication_P3", "prereg_runs", "P3_full_rho*_s*.pt")))
    out = {"models": {}}; res = {}
    for pth in paths:
        rho, s = PAT.search(pth).groups(); rho, s = float(rho), int(s)
        if rho not in (1.0, 0.5): continue
        net = P.Pol(P.P3Env.odim, P.P3Env.nmove, P.P3Env.k); net.load_state_dict(torch.load(pth)); net.eval()
        H, c, l = hidden_states(net)
        m = {"rho": rho, "seed": s, "cue_decode": {t: round(decode(H[t], c), 4) for t in (1, 2, 3)},
             "landmark_decode": {t: round(decode(H[t], l), 4) for t in (1, 2, 3)}}
        res[(rho, s)] = m; out["models"][os.path.relpath(pth, HERE)[:-3]] = m; print(json.dumps(m), flush=True)
    c05 = [m["cue_decode"][3] for (r, s), m in res.items() if r == 0.5 and s in (1, 2, 3)]
    l1 = [m["landmark_decode"][3] for (r, s), m in res.items() if r == 1.0 and s in (1, 2, 3)]
    out["controls"] = {"rho0.5_cue_decode_t3_median": float(np.median(c05)), "rho1_landmark_decode_t3_median": float(np.median(l1))}
    out["controls"]["valid"] = bool(np.median(c05) >= 0.90 and np.median(l1) >= 0.90)
    for name, seeds in (("original", (1, 2, 3)), ("replication", (4, 5, 6))):
        vals = [res[(1.0, s)]["cue_decode"][3] for s in seeds if (1.0, s) in res]
        out["verdict_" + name] = {"rho1_cue_decode_t3": vals, "verdict": verdict(vals)}
    json.dump(out, open(os.path.join(HERE, "P3_probe.json"), "w"), indent=1)
    print("CONTROLS:", json.dumps(out["controls"])); print("ORIGINAL:", json.dumps(out["verdict_original"]))
    print("REPLICATION:", json.dumps(out["verdict_replication"]))
