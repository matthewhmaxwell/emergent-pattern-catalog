"""Does metastability earn its place? Direct scalar series: telegraph (two-state Poisson
switching = the positive), unimodal_stationary (AR(1) noise = null), slow_ramp (drift, must NOT
read as switching), oscillation (sine = multimodal but REGULAR dwells -> dwell_cv must separate
it from telegraph). Want: telegraph n>=2 & dwell_cv high; unimodal n=1; ramp n=1; oscillation
n>=2 but dwell_cv LOW."""
import numpy as np
from epc.metrics.metastability import metastability

T = 400


def S(series):
    return [{"order_parameter": float(v)} for v in series]


def telegraph(seed=0):
    r = np.random.default_rng(seed)
    state, out = 1.0, []
    for _ in range(T):
        if r.random() < 0.03:           # rare, random switches -> irregular dwells
            state = -state
        out.append(state + 0.25 * r.standard_normal())
    return np.array(out)


def unimodal(seed=0):
    r = np.random.default_rng(seed)
    x, out = 0.0, []
    for _ in range(T):
        x = 0.9 * x + r.standard_normal()
        out.append(x)
    return np.array(out)


def ramp(seed=0):
    r = np.random.default_rng(seed)
    return np.linspace(-2, 2, T) + 0.1 * r.standard_normal(T)


def oscillation(seed=0):
    r = np.random.default_rng(seed)
    return np.sin(np.linspace(0, 30 * np.pi, T)) + 0.05 * r.standard_normal(T)


print(f"{'system':<20}{'n_macrostates':>14}{'dwell_cv':>10}{'bimodality':>12}")
rows = {}
for nm, fn in [("telegraph", telegraph), ("unimodal_stationary", unimodal),
               ("slow_ramp", ramp), ("oscillation", oscillation)]:
    r = metastability(S(fn(1)))
    rows[nm] = r
    print(f"{nm:<20}{r['n_macrostates']:>14}{r['dwell_cv']:>10.3f}{r['bimodality']:>12.3f}")

tg = rows["telegraph"]; uni = rows["unimodal_stationary"]; osc = rows["oscillation"]; rmp = rows["slow_ramp"]
switch_ok = tg["n_macrostates"] >= 2 and tg["dwell_cv"] > 0.5
null_ok = uni["n_macrostates"] == 1 and rmp["n_macrostates"] == 1
osc_sep = tg["dwell_cv"] > 2 * osc["dwell_cv"]   # telegraph dwells far more irregular than periodic
print(f"\ntelegraph switching {switch_ok} | unimodal+ramp single-mode {null_ok} | dwell_cv telegraph {tg['dwell_cv']:.2f} vs oscillation {osc['dwell_cv']:.2f} -> {osc_sep}")
print(f"VERDICT: {'ADMIT' if switch_ok and null_ok and osc_sep else 'review'}")
