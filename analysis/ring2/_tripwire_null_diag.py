"""Why did the model-free tripwire flag a UNIFORM-noise field as complex? Print the
internals (C, psi, macro_feat, thresholds) for the uniform null vs a real RD pattern
vs a couple of existing corpus nulls, to locate the false-positive mechanism."""
import numpy as np
from epc.models.reaction_diffusion import reaction_diffusion
from epc.phase2a.novelty_tripwire import model_free_complexity, C_THR, PSI_THR
from analysis.blind_spot_probes import null_noise, null_walk


def uniform_null(seed=0, std=0.001):
    rng = np.random.default_rng(seed)
    return [{"field": 0.5 + std * rng.standard_normal((96, 96)), "step": t} for t in range(24)]


def flat_zero_var(seed=0):
    return [{"field": np.full((96, 96), 0.5), "step": t} for t in range(24)]


cases = [
    ("uniform_noise(1e-3)", uniform_null(0, 0.001)),
    ("uniform_noise(1e-1)", uniform_null(0, 0.1)),
    ("flat_zero_var", flat_zero_var()),
    ("RD_stripes", reaction_diffusion(0, F=0.030, k=0.057, N=96, steps=8000, record=24)),
    ("probe_null_noise", null_noise(0)["history"]),
    ("probe_null_walk", null_walk(0)["history"]),
]
print(f"thresholds: C_THR={C_THR}  PSI_THR={PSI_THR}")
print(f"{'case':<22}{'C':>8}{'psi':>8}{'is_complex':>11}  macro_feat")
for name, h in cases:
    mf = model_free_complexity(h)
    C = mf['C'] if mf['C'] == mf['C'] else float('nan')
    print(f"{name:<22}{C:>8.4f}{mf['psi']:>8.4f}{str(mf['is_complex']):>11}  {mf['macro_feat']}")
