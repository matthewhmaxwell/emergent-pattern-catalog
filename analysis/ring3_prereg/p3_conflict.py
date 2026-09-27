"""P3 MISS -> novelty protocol (osf.io/up7m9, PREREG §4), step 2: cue-landmark CONFLICT test + resample ablation.
Criteria below are fixed BEFORE this script is run on any trained model (committed to git first).

Registered result under test: P3 verdict MISS. The registered memwipe ablation (recurrent state zeroed at every step)
collapses performance even at rho=1 (S 0.93 -> 0.55-0.61, nu 0.75-0.89), although the memoryless-trained fair
baseline scores 0.92 there, and memwipe scores fall BELOW chance at low rho (0.27-0.41). Two readings:
  (M) memory-first: the rho=1 agent navigates by its stored t=0 cue, not the always-present landmark -> genuine MISS.
  (A) artifact: the agent follows the landmark; zeroing the hidden state is out-of-distribution and breaks the policy.

Conditions (per trained model; 4000 episodes = 10 x 400; actions sampled as in DEVIATIONS D1; episode seeds 950000+
are disjoint from every registered evaluation seed; the same episode seeds are used in every condition -> paired):
  normal   the training distribution; outcome split correct / wrong / timeout (sanity check vs the registered S).
  memwipe  the registered ablation (h = 0 at every step), re-run with the outcome split.
  conflict the t=0 cue marks goal c; the landmark marks the OTHER goal (1-c) at every step. Record the goal reached
           first: cue goal / landmark goal / timeout. For rho>0 models each input is in-distribution on its own; only
           their disagreement is new (the classic cue-conflict paradigm).
  resample in-distribution ablation of memory CONTENT: at every step t>=1 each episode's hidden state is replaced by
           the hidden state of a randomly drawn other still-running episode at the same step (fresh draw every step).
           Hidden-state statistics stay on-distribution; the episode's own memory is destroyed.
Controls (the test is INVALID if either fails; seed median over the original seeds 1-3):
  C-landmark  memoryless rho=1 models (only the landmark is visible after t=0) under conflict: landmark share >= 0.80.
  C-cue       full rho=0 models (never trained with a landmark) under conflict: cue share >= 0.80.
Scored quantities for full rho=1 models, per seed:
  L     = landmark share of decided conflict episodes = n_landmark / (n_landmark + n_cue)
  nu_rs = (S_normal - S_resample) / (S_normal - 0.50)
Decision rule (the registered rule: holds in >= 2 of 3 seeds AND at the seed median; applied to the original seeds
1-3; replication seeds 4-6 are judged separately with the same rule):
  PERCEPTION-FIRST if L >= 0.80;  MEMORY-FIRST if L <= 0.20;  MIXED otherwise.
  memory content DISPENSABLE if nu_rs <= 0.15;  NECESSARY if nu_rs >= 0.50;  PARTIAL otherwise.
Interpretation (fixed now):
  PERCEPTION-FIRST + DISPENSABLE -> (A): the registered memwipe collapse at rho=1 is an ablation artifact; behaviour
      matches C2's prediction. The formal P3 verdict (MISS) stays in the record with this attribution.
  MEMORY-FIRST + NECESSARY      -> (M): the competing hypothesis holds; genuine MISS -> step 3 literature pass.
  any other combination          -> the agent integrates both cues; reported as such, no artifact attribution.
Exploratory (not scored): L and nu_rs as functions of rho across all full models; timeout share under memwipe.

  python p3_conflict.py --smoke     # execution check on an UNTRAINED network only (no trained model is loaded)
  python p3_conflict.py             # the test
"""
import os, re, sys, json, glob, numpy as np, torch
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import prereg_envs as P

B, REPS, SEED0 = 400, 10, 950000
PAT = re.compile(r"P3_(full|memoryless)_rho([0-9.]+)_s(\d+)\.pt$")


class ConflictEnv(P.P3Env):
    """t=0 cue on goal `correct` (= the cue goal); landmark on the OTHER goal at every step."""
    def obs(self):
        o = super().obs(); r = np.arange(self.B)
        o[:, 0, 8:10] = 0.0; o[r, 0, 8 + (1 - self.correct)] = 1.0
        return o


def rollout(net, env, hmode, rng):
    """hmode: 'keep' (recurrent), 'zero' (h=0 every step: memwipe / memoryless), 'resample' (donor hidden state)."""
    obs = env.obs(); h = None
    for t in range(env.T):
        hin = None if hmode == "zero" else h
        if hmode == "resample" and t >= 1:
            pool = np.where(~env.done)[0]
            if len(pool) < 2: pool = np.arange(env.B)
            idx = rng.choice(pool, size=env.B)
            for _ in range(10):
                same = idx == np.arange(env.B)
                if not same.any(): break
                idx[same] = rng.choice(pool, size=int(same.sum()))
            hin = h[:, torch.from_numpy(idx)].contiguous()
        with torch.no_grad(): ml, _, _, h = net(torch.from_numpy(obs[:, 0])[:, None, :], hin)
        m = torch.distributions.Categorical(logits=ml[:, 0]).sample().numpy()
        env.step(m[:, None], np.zeros((env.B, 1), int)); obs = env.obs()
    return env


def outcomes(net, rho, hmode, conflict=False, reps=REPS):
    torch.manual_seed(SEED0); rng = np.random.default_rng(SEED0); c = w = to = 0
    for r in range(reps):
        env = (ConflictEnv if conflict else P.P3Env)(B, SEED0 + r, rho=rho)
        rollout(net, env, hmode, rng)
        c += int(env.right.sum()); w += int((env.done & ~env.right).sum()); to += int((~env.done).sum())
    n = reps * B
    k = ("cue", "landmark") if conflict else ("correct", "wrong")
    return {k[0]: round(c / n, 4), k[1]: round(w / n, 4), "timeout": round(to / n, 4)}


def audit(net, rho, memoryless, reps=REPS):
    base = "zero" if memoryless else "keep"
    r = {"rho": rho, "memoryless": memoryless, "normal": outcomes(net, rho, base, reps=reps)}
    cf = outcomes(net, rho, base, conflict=True, reps=reps); r["conflict"] = cf
    dec = cf["cue"] + cf["landmark"]; r["L_landmark_share"] = round(cf["landmark"] / dec, 4) if dec > 0 else None
    if not memoryless:
        r["memwipe"] = outcomes(net, rho, "zero", reps=reps)
        r["resample"] = rs = outcomes(net, rho, "resample", reps=reps)
        S = r["normal"]["correct"]
        r["nu_resample"] = round((S - rs["correct"]) / (S - 0.5), 4) if S > 0.5 else None
    return r


def rule(vals, pred):
    vals = [v for v in vals if v is not None]
    return len(vals) >= 2 and sum(pred(v) for v in vals) >= 2 and pred(float(np.median(vals)))


def verdict(models, seeds):
    full1 = {s: m for (v, rho, s), m in models.items() if v == "full" and rho == 1.0 and s in seeds}
    if len(full1) < 2: return None
    Ls = [m["L_landmark_share"] for m in full1.values()]; nus = [m["nu_resample"] for m in full1.values()]
    beh = ("PERCEPTION-FIRST" if rule(Ls, lambda x: x >= 0.80) else
           "MEMORY-FIRST" if rule(Ls, lambda x: x <= 0.20) else "MIXED")
    mem = ("DISPENSABLE" if rule(nus, lambda x: x <= 0.15) else
           "NECESSARY" if rule(nus, lambda x: x >= 0.50) else "PARTIAL")
    interp = ("(A) ablation artifact: agent follows the landmark; registered memwipe collapse at rho=1 is OOD"
              if (beh, mem) == ("PERCEPTION-FIRST", "DISPENSABLE") else
              "(M) memory-first: competing hypothesis holds; genuine MISS -> literature pass"
              if (beh, mem) == ("MEMORY-FIRST", "NECESSARY") else
              "integrates both cues; no artifact attribution")
    return {"seeds": sorted(full1), "L_by_seed": dict(zip(sorted(full1), [full1[s]["L_landmark_share"] for s in sorted(full1)])),
            "nu_resample_by_seed": dict(zip(sorted(full1), [full1[s]["nu_resample"] for s in sorted(full1)])),
            "behaviour": beh, "memory_content": mem, "interpretation": interp}


if __name__ == "__main__":
    if "--smoke" in sys.argv:                     # untrained network: execution check only
        torch.manual_seed(0); net = P.Pol(P.P3Env.odim, P.P3Env.nmove, P.P3Env.k); net.eval()
        print(json.dumps(audit(net, 1.0, False, reps=1))); print(json.dumps(audit(net, 1.0, True, reps=1)))
        print("SMOKE OK"); sys.exit(0)
    paths = sorted(glob.glob(os.path.join(HERE, "prereg_runs", "P3_*.pt")) +
                   glob.glob(os.path.join(HERE, "replication_P3", "prereg_runs", "P3_*.pt")))
    models, out = {}, {"models": {}}
    for pth in paths:
        v, rho, s = PAT.search(pth).groups(); rho, s = float(rho), int(s)
        net = P.Pol(P.P3Env.odim, P.P3Env.nmove, P.P3Env.k); net.load_state_dict(torch.load(pth)); net.eval()
        m = audit(net, rho, v == "memoryless")
        reg = os.path.join(os.path.dirname(pth), os.path.basename(pth)[:-3] + ".json")
        if os.path.exists(reg): m["registered"] = {k: x for k, x in json.load(open(reg)).items() if k.startswith("S")}
        models[(v, rho, s)] = m; tag = os.path.relpath(pth, HERE)[:-3]; out["models"][tag] = m
        print(tag, json.dumps(m), flush=True)
    orig = (1, 2, 3)
    cl = [m["L_landmark_share"] for (v, rho, s), m in models.items() if v == "memoryless" and rho == 1.0 and s in orig]
    cc = [1 - m["L_landmark_share"] for (v, rho, s), m in models.items() if v == "full" and rho == 0.0 and s in orig
          and m["L_landmark_share"] is not None]
    out["controls"] = {"C_landmark_median_landmark_share": float(np.median(cl)) if cl else None,
                       "C_cue_median_cue_share": float(np.median(cc)) if cc else None}
    out["controls"]["valid"] = bool(cl and cc and np.median(cl) >= 0.80 and np.median(cc) >= 0.80)
    out["verdict_original_seeds"] = verdict(models, orig)
    out["verdict_replication_seeds"] = verdict(models, (4, 5, 6))
    expl = {}
    for rho in sorted({rho for (v, rho, s) in models if v == "full"}, reverse=True):
        ms = [m for (v, r, s), m in models.items() if v == "full" and r == rho and s in orig]
        expl[str(rho)] = {"median_L": float(np.median([m["L_landmark_share"] for m in ms if m["L_landmark_share"] is not None])),
                          "median_nu_resample": float(np.median([m["nu_resample"] for m in ms if m["nu_resample"] is not None])),
                          "median_memwipe_timeout": float(np.median([m["memwipe"]["timeout"] for m in ms]))}
    out["exploratory_by_rho"] = expl
    json.dump(out, open(os.path.join(HERE, "P3_conflict.json"), "w"), indent=1)
    print("CONTROLS:", json.dumps(out["controls"])); print("ORIGINAL:", json.dumps(out["verdict_original_seeds"]))
    print("REPLICATION:", json.dumps(out["verdict_replication_seeds"])); print("BY RHO:", json.dumps(expl))
