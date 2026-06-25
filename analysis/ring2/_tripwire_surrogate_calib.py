"""Calibrate the surrogate structure-score that fixes the tripwire false-positive.
structure_score = mean(H_shuffle) - H_obs  (permutation-entropy deficit vs time-shuffled
surrogates of the mean-macro): ~0 for iid/exchangeable series, >0 for genuine temporal
structure; bias-cancels (shuffle shares the finite-size sparsity). Goal: a STRUCT_THR
that (a) kills the uniform-noise-field false positive, (b) keeps genuinely-complex
systems above it, (c) doesn't regress the corpus nulls."""
import numpy as np
from analysis.blind_spot_probes import PROBES
from epc.phase2a.info_channels import micro_macro, mpr_complexity, psi_ce_best
from epc.models.reaction_diffusion import reaction_diffusion

C_THR, PSI_THR = 0.16, 0.05


def structure_score(series, n_surr=24, seed=0):
    s = np.asarray(series, float).ravel()
    if s.size < 8:
        return 0.0
    H_obs = mpr_complexity(s)["H"]
    rng = np.random.default_rng(seed)
    Hs = []
    for _ in range(n_surr):
        sh = s.copy(); rng.shuffle(sh); Hs.append(mpr_complexity(sh)["H"])
    return float(np.mean(Hs) - H_obs)


def row(name, truth, history):
    micro, cands = micro_macro(history)
    if cands is None:
        return (name, truth, None)
    C = float(mpr_complexity(cands.get("mean"))["C"])
    st = structure_score(cands.get("mean"))
    psi, _ = psi_ce_best(micro, cands)
    psi = float(psi) if psi == psi else 0.0
    return (name, truth, {"C": C, "struct": st, "psi": psi})


def uniform_null(seed=0, std=0.001):
    rng = np.random.default_rng(seed)
    return [{"field": 0.5 + std * rng.standard_normal((96, 96)), "step": t} for t in range(24)]


rows = [row(p.__name__, p(0).get("truth"), p(0)["history"]) for p in PROBES]
rows.append(row("UNIFORM_field_null", "null(OOD)", uniform_null()))
rows.append(row("RD_stripes", "emergent(OOD)", reaction_diffusion(0, F=0.030, k=0.057, N=96, steps=8000, record=24)))
rows.append(row("RD_spots", "emergent(OOD)", reaction_diffusion(0, F=0.034, k=0.063, N=96, steps=8000, record=24)))

print(f"{'system':<24}{'truth':<14}{'C':>7}{'struct':>8}{'psi':>8}{'C>thr':>6}")
for nm, t, r in sorted(rows, key=lambda x: -((x[2] or {}).get('struct', -9))):
    if r:
        print(f"{nm:<24}{str(t):<14}{r['C']:>7.3f}{r['struct']:>8.3f}{r['psi']:>8.2f}{('Y' if r['C']>C_THR else '-'):>6}")
    else:
        print(f"{nm:<24}{str(t):<14}   (n/a)")

# focus: among systems with C>C_THR (the ones the buggy C-path would fire on), how does
# struct separate genuine emergent from null?
hot = [(t, r) for nm, t, r in rows if r and r['C'] > C_THR]
null_hot = [r['struct'] for t, r in hot if 'null' in str(t)]
emg_hot = [r['struct'] for t, r in hot if 'emergent' in str(t)]
print(f"\nAmong C>{C_THR}:  null struct: {[round(x,3) for x in sorted(null_hot)]}")
print(f"               emergent struct: {[round(x,3) for x in sorted(emg_hot)]}")
if null_hot and emg_hot:
    print(f"  -> STRUCT_THR between {max(null_hot):.3f} (null max) and {min(emg_hot):.3f} (emergent min)")
