"""Within-episode dynamics of division of labor — the biologically-honest measure.

Reframe (user): demanding a behavior be maintained across ALL training seeds is an inherited detector-calibration
criterion, not a property of the behavior. In biology, task allocation is dynamic — reliably ARISES, is
functional when present, and is abandoned/reassigned over time. That non-permanence is normal, not degenerate.
So instead of scoring the final policy across seeds, measure what happens WITHIN a single episode on policies
that express the behavior:

  1. FORMATION   — does a valid allocation (each target covered by exactly one agent) arise? time-to-form.
  2. OCCUPANCY   — once formed, what fraction of remaining steps is coverage held? (duty cycle; "present but not permanent")
  3. RECOVERY    — displace one agent off its target mid-episode; does the pair re-establish coverage? rate + time.
                   (the sharp test: a rigid one-shot solution cannot recover; a dynamically-maintained one can.)
  4. FLUIDITY    — how often does WHICH agent covers WHICH target switch during the episode? ("more alive": reassignment, not a frozen lock.)

Policies are trained fresh via curriculum (high hit-rate for expressing DoL), kept in memory, probed on an
EXTENDED horizon (T_eval > T_train) so a form->perturb->recover cycle fits in one episode. Eval always at full
difficulty (random targets). No probe changes the reward or task.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, torch, torch.nn as nn
torch.set_num_threads(6)
import allocworld as A

N = A.N
HERE = os.path.dirname(os.path.abspath(__file__))
RESULT = os.path.join(HERE, "alloc_dynamics_result.json")

def train_curriculum(cfg, iters, seed, B=192, curr_frac=0.6):
    torch.manual_seed(seed); n = cfg["n_ag"]
    base = dict(cfg); base["difficulty"] = 1.0
    od = A.AllocWorld(1, 0, base).odim
    net = A.Pol(od, n=n); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    ramp = int(iters * curr_frac)
    for it in range(iters):
        c = dict(cfg); c["difficulty"] = min(1.0, it / ramp) if ramp > 0 else 1.0
        env = A.AllocWorld(B, 1000 + it, c)
        O, MV, SG, LP, V, R = A.rollout(net, env)
        O2 = O.reshape(B*n, A.T, od); MV2, SG2 = MV.reshape(B*n, A.T), SG.reshape(B*n, A.T)
        LP2, V2, R2 = LP.reshape(B*n, A.T), V.reshape(B*n, A.T), R.reshape(B*n, A.T)
        adv, ret = A.gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, MVt, SGt, LPt = torch.from_numpy(O2), torch.from_numpy(MV2), torch.from_numpy(SG2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            ml, sl, v, _ = net(Ot)
            md = torch.distributions.Categorical(logits=ml); sd = torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(md.log_prob(MVt) + sd.log_prob(SGt) - LPt)
            s1 = ratio*advt; s2 = torch.clamp(ratio, 0.8, 1.2)*advt
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.01*(md.entropy()+sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
    return net, base

def cover_match(env):
    """Per-batch: covered (every target covered by exactly one agent) + who[b,k]=agent id on target k."""
    B, n = env.B, env.n
    cov = np.zeros((B, n), int); who = np.full((B, n), -1, int)
    for k in range(n):
        for i in range(n):
            on = np.all(env.ap[:, i] == env.tg[:, k], axis=1)
            if env.hetero: on = on & (env.atype[:, i] == env.ttype[:, k])
            cov[:, k] += on.astype(int); who[on, k] = i
    covered = (cov == 1).all(axis=1)
    return covered, who

def probe(net, cfg, seed, B=400, T_eval=60, perturb_at=None):
    """Instrumented greedy rollout. Returns cov_ts (T_eval x B bool), who_ts (T_eval x B x n)."""
    n = cfg["n_ag"]; env = A.AllocWorld(B, seed, {**cfg, "difficulty": 1.0})
    obs = env.obs(); h = [None]*n; cov_ts = []; who_ts = []
    for t in range(T_eval):
        if perturb_at is not None and t == perturb_at:
            env.ap[:, 0] = env.rng.integers(0, N, size=(B, 2))   # displace agent 0 to a random cell
            obs = env.obs()
        mv = np.zeros((B, n), int); sg = np.zeros((B, n), int)
        for i in range(n):
            with torch.no_grad():
                ml, sl, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
            mv[:, i] = ml[:, 0].argmax(1).numpy(); sg[:, i] = sl[:, 0].argmax(1).numpy()
        env.step(mv, sg); obs = env.obs()
        covered, who = cover_match(env); cov_ts.append(covered); who_ts.append(who)
    return np.array(cov_ts), np.array(who_ts)

def analyze(cov, who, perturb_at=None):
    """cov: T x B bool; who: T x B x n. Returns dynamics metrics."""
    T, B = cov.shape; n = who.shape[2]
    out = {}
    # FORMATION
    ever = cov.any(axis=0)                               # B: formed at least once
    out["formation_rate"] = float(ever.mean())
    firsts = np.array([np.argmax(cov[:, b]) if ever[b] else -1 for b in range(B)])
    out["time_to_form"] = float(np.mean([firsts[b] for b in range(B) if ever[b]])) if ever.any() else None
    # OCCUPANCY: fraction of steps covered after first formation
    occ = []
    for b in range(B):
        if ever[b]:
            occ.append(cov[firsts[b]:, b].mean())
    out["occupancy"] = float(np.mean(occ)) if occ else None
    # FLUIDITY: among covered steps, does the matching signature switch? (n=2: who[:,0])
    switchers = 0; nformed = 0
    for b in range(B):
        if not ever[b]: continue
        nformed += 1
        sig = who[cov[:, b], b, 0]                       # agent covering target0 at each covered step
        if len(sig) > 1 and (np.diff(sig) != 0).any(): switchers += 1
    out["fluidity_rate"] = float(switchers / nformed) if nformed else None
    # RECOVERY (perturbed run): among b covered just before perturbation, do they re-cover after?
    if perturb_at is not None and perturb_at >= 1:
        pre = cov[perturb_at - 1]                        # covered right before the displacement
        rec = 0; rectimes = []
        idx = np.where(pre)[0]
        for b in idx:
            post = cov[perturb_at:, b]
            if post.any():
                rec += 1; rectimes.append(int(np.argmax(post)))
        out["n_perturbed"] = int(len(idx))
        out["recovery_rate"] = float(rec / len(idx)) if len(idx) else None
        out["recovery_time"] = float(np.mean(rectimes)) if rectimes else None
    return out

if __name__ == "__main__":
    a = sys.argv
    nseeds = int(a[a.index("--seeds")+1]) if "--seeds" in a else 4
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 1200
    T_eval = int(a[a.index("--teval")+1]) if "--teval" in a else 60
    perturb_at = int(a[a.index("--perturb")+1]) if "--perturb" in a else 30
    cfg = dict(n_ag=2, has_obs=True, has_sig=False, hetero=False, reward_shared=True, asym_init=False)
    results = []
    for s in range(nseeds):
        t0 = time.time(); net, base = train_curriculum(cfg, iters, s)
        cov_u, who_u = probe(net, base, 7000, T_eval=T_eval, perturb_at=None)
        cov_p, who_p = probe(net, base, 7000, T_eval=T_eval, perturb_at=perturb_at)
        m = analyze(cov_u, who_u); mp = analyze(cov_p, who_p, perturb_at=perturb_at)
        rec = dict(seed=s, sec=round(time.time()-t0, 1), **m,
                   recovery_rate=mp.get("recovery_rate"), recovery_time=mp.get("recovery_time"),
                   n_perturbed=mp.get("n_perturbed"),
                   coverage_curve=[round(float(x), 3) for x in cov_u.mean(axis=1)],
                   coverage_curve_perturbed=[round(float(x), 3) for x in cov_p.mean(axis=1)])
        results.append(rec)
        print(f"s{s}: form={rec['formation_rate']:.2f} t2form={rec['time_to_form']:.1f} occ={rec['occupancy']:.2f} "
              f"fluid={rec['fluidity_rate']:.2f} recov={rec['recovery_rate']} ({rec['sec']}s)", flush=True)
        json.dump({"perturb_at": perturb_at, "T_eval": T_eval, "seeds": results}, open(RESULT, "w"), indent=1)
    print("DONE", flush=True)
