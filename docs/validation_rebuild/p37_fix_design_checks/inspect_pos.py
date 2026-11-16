# P37 v2 design check, run on the five canonical POSITIVES only (2026-11-15).
# Autocorrelations and periodogram peaks of the late-window abundance series.
# Run: PYTHONPATH=. python docs/validation_rebuild/p37_fix_design_checks/inspect_pos.py
import numpy as np
from epc.models.resource_competition import resource_competition
np.set_printoptions(precision=2, suppress=True, linewidth=160)
for s in range(5):
    fr, _ = resource_competition(s)
    A = np.array([f["abundances"] for f in fr]); R = np.array([f["resources"] for f in fr])
    T = len(A); late = A[T//2:]; N = len(late)
    print(f"seed {s} T={T} late N={N} step spacing={fr[1]['step']-fr[0]['step']} (dt=0.02 -> {0.02*(fr[1]['step']-fr[0]['step'])} time units)")
    for i in range(A.shape[1]):
        x = late[:, i]; t = np.arange(N)
        r = x - np.polyval(np.polyfit(t, x, 1), t)
        ac = np.array([np.corrcoef(r[:-k], r[k:])[0,1] for k in range(1, N//2)])
        P = np.abs(np.fft.rfft(r))**2
        print(f"  sp{i} mean={x.mean():.3f} cv_detr={r.std()/x.mean():.3f} ACF1..8={ac[:8]} argmaxP(k)={P[1:].argmax()+1} topk={np.argsort(P[1:])[::-1][:4]+1}")
    print("  full-series sp0 first 40:", A[:40,0])
