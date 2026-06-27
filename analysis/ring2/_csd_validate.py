"""Does critical_slowing earn its place? OU process with VANISHING restoring rate lambda(t):1->0.05
(canonical CSD: AR(1)->1, variance rises = approaching tipping) vs stationary OU (flat) vs a
drift_confound (stationary OU + linear mean ramp -- must NOT trip after detrending). Want:
approaching ews_ar1_trend strongly positive; stationary ~0; drift ~0 (confound defused)."""
import numpy as np
from epc.metrics.critical_slowing import critical_slowing

T = 600


def S(x):
    return [{"scalar": float(v)} for v in x]


BURN = 300


def approaching(seed=0):
    r = np.random.default_rng(seed)
    x = 0.0
    for _ in range(BURN):                       # equilibrate at lambda=1 before recording
        x = x - 1.0 * x * 0.1 + 0.1 * r.standard_normal()
    lam = np.linspace(1.0, 0.02, T)            # restoring rate -> ~0 (recovery slows, strong CSD)
    out = []
    for t in range(T):
        x = x - lam[t] * x * 0.1 + 0.1 * r.standard_normal()
        out.append(x)
    return np.array(out)


def stationary(seed=0):
    r = np.random.default_rng(seed)
    x = 0.0
    for _ in range(BURN):
        x = x - 0.5 * x * 0.1 + 0.1 * r.standard_normal()
    out = []
    for _ in range(T):
        x = x - 0.5 * x * 0.1 + 0.1 * r.standard_normal()
        out.append(x)
    return np.array(out)


def drift_confound(seed=0):
    return stationary(seed) + np.linspace(0, 3, T)   # slow mean drift, constant dynamics


print(f"{'system':<16}{'ews_ar1_rise':>14}{'ews_var_rise':>14}")
rows = {}
for nm, fn in [("approaching", approaching), ("stationary", stationary), ("drift_confound", drift_confound)]:
    # average over seeds to reduce single-realization noise
    rs = [critical_slowing(S(fn(sd))) for sd in range(1, 6)]
    r = {"ews_ar1_rise": float(np.mean([z["ews_ar1_rise"] for z in rs])),
         "ews_var_rise": float(np.mean([z["ews_var_rise"] for z in rs]))}
    rows[nm] = r
    print(f"{nm:<16}{r['ews_ar1_rise']:>14.3f}{r['ews_var_rise']:>14.3f}")

apv = rows["approaching"]["ews_var_rise"]; nullv = max(abs(rows["stationary"]["ews_var_rise"]), abs(rows["drift_confound"]["ews_var_rise"]))
apa = rows["approaching"]["ews_ar1_rise"]; nulla = max(rows["stationary"]["ews_ar1_rise"], rows["drift_confound"]["ews_ar1_rise"])
print(f"\nvar_rise: approaching {apv:.3f} vs null {nullv:.3f} (primary EWS) | ar1_rise: approaching {apa:.3f} vs null {nulla:.3f} (corroborating, right sign)")
print(f"VERDICT: {'ADMIT (CSD via rising variance EWS; ar1 corroborates; drift confound defused by per-segment linear detrend)' if apv > 0.3 and apv > 5 * nullv and apa > nulla else 'review'}")
