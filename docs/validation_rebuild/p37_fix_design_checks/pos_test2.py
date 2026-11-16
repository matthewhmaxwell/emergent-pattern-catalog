# P37 v2 design check, run on the five canonical POSITIVES only (2026-11-15).
# Pooled band statistic vs permutation surrogates, and lag-0 abundance-resource
# correlations. Run: PYTHONPATH=. python docs/validation_rebuild/p37_fix_design_checks/pos_test2.py
import numpy as np
from epc.models.resource_competition import resource_competition
from epc.metrics.competition_dynamics import _detrend, _band_bins, band_peak_fraction
np.set_printoptions(precision=2, suppress=True)
for s in range(5):
    fr = resource_competition(s)[0]
    A = np.array([f["abundances"] for f in fr]); R = np.array([f["resources"] for f in fr])
    T = len(A); n = T - int(round(T*0.5)); W = A[n:]; RW = R[n:]
    Ar = np.column_stack([_detrend(W[:, i]) for i in range(5)]); Rr = np.column_stack([_detrend(RW[:, j]) for j in range(5)])
    C0 = np.corrcoef(Ar.T, Rr.T)[:5, 5:]
    print(f"seed {s} lag-0 corr(abund_i, res_j):\n{C0}\n  resource mean {RW.mean(0)} cv {Rr.std(0)/RW.mean(0)}")
    valid = _band_bins(len(W), 3, 0.45)
    rng = np.random.default_rng(42)
    S = np.array([band_peak_fraction(Ar[:, i], valid) for i in range(5)])
    null = np.array([np.mean([band_peak_fraction(rng.permutation(Ar[:, i]), valid) for i in range(5)]) for _ in range(199)])
    z = (S.mean() - null.mean()) / null.std()
    print(f"  pooled S={S.mean():.3f} null mean={null.mean():.3f} sd={null.std():.3f} z={z:.2f} p={(1+np.sum(null>=S.mean()))/200:.3f}")
