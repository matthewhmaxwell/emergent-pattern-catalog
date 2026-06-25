"""Exercise the broadened lens battery on a NEW substrate family (Gray-Scott RD), OOD
from the agent/network catalog corpus. Walk the (F,k) atlas and report what each field
lens + the generic emergence indicator + the model-free tripwire see per regime. This is
a HARDENING run: do the lenses behave sensibly (fire on patterned regimes, stay quiet on
uniform death) on a substrate they were never tuned against?"""
import numpy as np
from epc.models.reaction_diffusion import reaction_diffusion
from epc.metrics.structure_factor import structure_factor_peak
from epc.metrics.persistent_homology import persistent_homology
from epc.metrics.fractal_dimension import fractal_dimension
from epc.phase2a.novelty_tripwire import novelty_tripwire

REGIMES = [
    ("solitons",   0.022, 0.051),
    ("mitosis",    0.026, 0.051),
    ("stripes",    0.030, 0.057),
    ("spots",      0.034, 0.063),
    ("worms",      0.039, 0.058),
    ("high-F",     0.060, 0.062),   # NB: patterned (v.std high), NOT a null
]


def _uniform_null(seed=0):           # genuine null: flat field + tiny noise
    rng = np.random.default_rng(seed)
    return [{"field": 0.5 + 0.001 * rng.standard_normal((96, 96)), "step": t} for t in range(24)]


print(f"{'regime':<14}{'sk_peak':>9}{'h1_max':>8}{'fracD':>7}{'em':>6} {'em_kind':<18}{'cplx':>5}{'trip':>5}")
WORK = [(name, reaction_diffusion(0, F=F, k=k, N=96, steps=8000, record=24)) for name, F, k in REGIMES]
WORK.append(("UNIFORM-null", _uniform_null()))
for name, h in WORK:
    sf = structure_factor_peak(h)
    ph = persistent_homology(h)
    fd = fractal_dimension(h)
    tw = novelty_tripwire(h)
    sk = sf.get("sk_peak") if sf else None
    h1 = ph.get("h1_max") if ph else None
    fD = fd.get("fractal_dim") if fd else None
    print(f"{name:<14}"
          f"{(f'{sk:9.1f}' if sk is not None else '      n/a')}"
          f"{(f'{h1:8.3f}' if h1 is not None else '     n/a')}"
          f"{(f'{fD:7.2f}' if fD is not None else '    n/a')}"
          f"{tw['em_score']:>6.2f} {str(tw['em_kind']):<18}"
          f"{('Y' if tw['is_complex'] else '-'):>5}{('Y' if tw['tripped'] else '-'):>5}")
print("\nSanity: patterned RD regimes -> sk_peak HIGH (field characteristic scale) + emergence"
      " fires (correct, RD is Turing-class). UNIFORM-null -> sk_peak ~1, em low, not complex.")
