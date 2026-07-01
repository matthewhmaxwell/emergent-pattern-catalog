"""Lens-space outlier analysis -- the discovery mechanism the comprehensive 21-lens battery
enables. The em/MPR-C tripwire flags complex+unclassified configs; but a genuinely novel structure
might NOT trip em/C yet still be a strong OUTLIER in the lens fingerprint space (unusual along one
of the 15 new axes the prior battery could not measure). This finds those.

Per family (substrate -> a stable set of lenses fire, so fingerprints are comparable): take the
numeric coordinates present in >=80% of configs, robust-z-score (median/MAD), and score each
config by its mean distance to its 3 nearest neighbours in z-space (kNN outlier score). Report the
top outliers and the coordinates that make each one extreme -> the candidate-novel frontier to vet.
Honest: an outlier is a CANDIDATE, not novelty; it still needs the catalogued-detector + literature
vet (a lattice or a known phase is an outlier too).
"""
import json
import sys
import glob

import numpy as np

NUMERIC = [
    "em_score", "mf_C", "mf_struct", "mf_psi", "sk_peak", "h1_max", "ph_components",
    "field_loop_area", "field_loops", "fractal_dim", "lacunarity", "ang_mom",
    "degree_cv", "modularity", "clustering", "dte_mean_te", "dte_directionality",
    "defect_density", "defect_coherence", "lattice_symmetry", "azimuthal_concentration",
    "velcorr_length", "polarization", "mill_strength", "n_macrostates", "dwell_cv",
    "growth_exponent", "coarsen_fit", "msd_exponent", "ergodicity_break", "local_R_spread",
    "global_sync", "zero_one_K", "determinism", "lz_ratio", "number_variance_exponent",
    "cross_coherence", "mf_width", "ews_ar1_rise", "ews_var_rise", "circulation", "extremal_index",
]


def _num(v):
    try:
        f = float(v)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


def analyse(path, topk=6):
    rows = [r for r in json.load(open(path)) if "error" not in r]
    if len(rows) < 8:
        print(f"{path}: too few configs ({len(rows)})"); return
    cols = [c for c in NUMERIC if sum(_num(r.get(c)) is not None for r in rows) >= 0.8 * len(rows)]
    X = np.array([[(_num(r.get(c)) if _num(r.get(c)) is not None else np.nan) for c in cols] for r in rows])
    med = np.nanmedian(X, axis=0)
    X = np.where(np.isnan(X), med, X)
    # RANK-transform each coord to [0,1] (robust to heavy tails/scale: no single coord dominates;
    # a config extreme on MANY axes scores high, a scale artifact on ONE axis does not).
    n = len(rows)
    R = np.zeros_like(X)
    for j in range(X.shape[1]):
        order_j = np.argsort(np.argsort(X[:, j], kind="mergesort"))
        R[:, j] = order_j / max(1, n - 1)
    keep = R.std(0) > 1e-9
    R = R[:, keep]; colk = [c for c, k in zip(cols, keep) if k]
    D = np.sqrt(((R[:, None, :] - R[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)
    k = min(3, n - 1)
    knn = np.sort(D, axis=1)[:, :k].mean(1)                  # kNN outlier score in rank space
    mu, sd = knn.mean(), knn.std() + 1e-9
    order = np.argsort(-knn)
    fam = path.split("stage2_results_")[-1].replace(".json", "")
    print(f"\n=== {fam}: {n} configs, {len(colk)} active lens coords ===")
    print(f"{'rank':>4} {'cfg':>4} {'score':>6} {'z':>5} {'em_kind':<22} {'trip':>5}  most-extreme rank-coords")
    for rank, idx in enumerate(order[:topk]):
        r = rows[idx]
        zz = sorted(zip(colk, R[idx]), key=lambda t: -abs(t[1] - 0.5))[:4]
        extremes = ", ".join(f"{c}={z:.2f}" for c, z in zz)
        z = (knn[idx] - mu) / sd
        print(f"{rank+1:>4} {r.get('i'):>4} {knn[idx]:>6.3f} {z:>5.1f} {str(r.get('em_kind')):<22} {str(r.get('tripped'))[:5]:>5}  {extremes}")


paths = sys.argv[1:] or sorted(glob.glob("analysis/ring2/stage2_results_*.json"))
for p in paths:
    analyse(p)
