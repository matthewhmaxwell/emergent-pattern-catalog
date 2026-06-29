"""Phase C verification (agent-observer debunk) of the `infer` agent's claimed in-episode rule
inference. Two interventions on the trained metaworld_infer.pt:
  - feedback ablation: zero the reward-feedback channel. A GENUINE inferer NEEDS feedback to identify
    the hidden good type -> precision must collapse toward chance (0.33). If it stays high, the "0.71"
    was a shortcut, not inference.
  - mid-episode flip: resample the good type at the midpoint. A genuine *continual* inferer notices its
    collections turning bad and re-infers -> post-flip precision recovers above chance.
Chance precision = 1/T = 0.33.
"""
import sys, torch, numpy as np
import metaworld_ppo as M

net = M.Policy(); net.load_state_dict(torch.load(M.__file__.rsplit("/", 1)[0] + "/metaworld_infer.pt")); net.eval()


def run(ablate=False, flip=None, seed=900, B=500):
    env = M.VecMeta(B, "infer", seed, ablate_fb=ablate, flip_at=flip)
    obs = env.reset(); h = None; snap = None
    for t in range(M.MAXT):
        if flip is not None and t == flip:
            snap = (env.good.copy(), env.bad.copy())
        with torch.no_grad():
            logit, _, h = net(torch.from_numpy(obs)[:, None, :], h)
        a = logit[:, 0].argmax(1).numpy(); obs, r, _, _ = env.step(a)
    coll = env.good + env.bad
    out = {"overall_prec": round(float((env.good[coll > 0] / coll[coll > 0]).mean()), 2)}
    if snap is not None:
        pg = env.good - snap[0]; pb = env.bad - snap[1]; pc = pg + pb
        out["post_flip_prec"] = round(float((pg[pc > 0] / pc[pc > 0]).mean()), 2)
    return out


print("infer-agent verification (chance precision = 0.33):")
print("  normal            :", run())
print("  feedback ABLATED  :", run(ablate=True), "<- genuine inference => collapses toward 0.33")
print("  mid-episode FLIP  :", run(flip=M.MAXT // 2), "<- genuine continual inference => post_flip recovers > 0.33")
