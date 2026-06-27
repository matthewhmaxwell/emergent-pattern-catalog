"""Does extreme_value earn its place? Same heavy tail, different temporal clustering of extremes:
iid heavy-tail (theta~1, extremes independent) vs persistent burst process (theta<1, extremes
cluster). The heavy-tail lens sees these as identical; theta must separate them."""
import numpy as np
from epc.metrics.extreme_value import extreme_value

T = 1500


def S(x):
    return [{"scalar": float(v)} for v in x]


def iid_heavytail(seed=0):
    return np.random.default_rng(seed).standard_t(2.5, T)        # heavy-tailed, independent in time


def clustered_extremes(seed=0):
    # persistent AR(1) driven by heavy-tailed innovations -> excursions ABOVE stay high (bursts)
    r = np.random.default_rng(seed); x = 0.0; out = []
    for _ in range(T):
        x = 0.85 * x + r.standard_t(2.5)
        out.append(x)
    return np.array(out)


print(f"{'system':<20}{'extremal_index':>16}")
rows = {}
for nm, fn in [("iid_heavytail", iid_heavytail), ("clustered_extremes", clustered_extremes)]:
    r = extreme_value(S(fn(1)))
    rows[nm] = r
    print(f"{nm:<20}{r['extremal_index']:>16.3f}")

iid = rows["iid_heavytail"]["extremal_index"]; cl = rows["clustered_extremes"]["extremal_index"]
print(f"\niid theta {iid:.2f} (~1, independent) vs clustered theta {cl:.2f} (<1, bursts)")
print(f"VERDICT: {'ADMIT (extremal index separates clustering at equal tail weight)' if iid > 0.8 and cl < iid - 0.2 else 'review'}")
