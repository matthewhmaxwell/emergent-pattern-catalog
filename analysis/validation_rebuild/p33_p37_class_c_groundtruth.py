"""Ground truth for the P33-P37 Class C (failed-regime) choices.

Class C must be the target's OWN model at parameters where the pattern does NOT
form (or with its mechanism disabled). This script establishes "does not form"
from each model's own physical observables, NOT from the detector verdict, so the
Class C regimes in epc/phase2a/failed_regimes/p3[3-7]_*.py are chosen without
looking at what the (frozen) detectors do with them.

  P33  lattice nematic: local nematic order S_loc and +-1/2 defect density, short
       (panel protocol, 400 steps) and long (4000 steps) runs, vs activity eta
       and initial condition. Pattern = sustained (steady-state) defect gas in a
       locally ordered field.
  P34  adaptive voter: final giant-component fraction, number of large
       components, opinion purity of components, vs rewiring probability phi.
       Pattern = the network FISSIONS into same-opinion components.
  P35  soft-disk compression: crystalline-particle fraction (local |psi6| > 0.7)
       of the final configuration vs packing fraction eta_end; the literature
       freezing range for hard disks is eta ~ 0.70-0.72 (Bernard & Krauth 2011).
  P37  chemostat: late-window coefficient of variation of abundances over a
       LONG run (fixed point vs sustained oscillation) and number of persisting
       species, for the standalone negatives and the added regimes.

Run:  PYTHONPATH=. python analysis/validation_rebuild/p33_p37_class_c_groundtruth.py [p33|p34|p35|p37|all]
Writes analysis/outputs/p33_p37_class_c_groundtruth.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from typing import Any, Dict

import numpy as np

OUT = "analysis/outputs/p33_p37_class_c_groundtruth.json"


def p33() -> Dict[str, Any]:
    from epc.models.active_nematic import active_nematic_field
    from epc.metrics.nematic_order import local_nematic_order, half_integer_defect_density
    rows = []
    for init in ("random", "aligned"):
        for eta in (0.0, 0.05, 0.10, 0.18, 0.30, 0.45, 0.60, 0.80, 1.20):
            for n_steps, tag in ((400, "panel_protocol"), (4000, "long")):
                t0 = time.time()
                fr, _ = active_nematic_field(0, G=64, eta=eta, n_steps=n_steps,
                                             n_frames=60, init_mode=init)
                th = [f["theta_field"] for f in fr]
                h = len(th)
                sl = [local_nematic_order(t) for t in th]
                hd = [half_integer_defect_density(t)[0] for t in th]
                rows.append({
                    "init": init, "eta": eta, "n_steps": n_steps, "run": tag,
                    "S_loc_first": sl[0], "S_loc_late": float(np.mean(sl[h // 2:])),
                    "hdef_first": hd[0],
                    "hdef_late_mean": float(np.mean(hd[h // 2:])),
                    "hdef_last": hd[-1],
                    "hdef_late_trend": float(hd[-1] - hd[h // 2]),
                    "secs": round(time.time() - t0, 1),
                })
                r = rows[-1]
                print(f"P33 init={init:7s} eta={eta:4.2f} {tag:14s} S_loc {r['S_loc_first']:.2f}->"
                      f"{r['S_loc_late']:.2f}  hdef first={r['hdef_first']:.4f} late={r['hdef_late_mean']:.4f} "
                      f"last={r['hdef_last']:.4f}", flush=True)
    return {"rows": rows}


def p34() -> Dict[str, Any]:
    from epc.models.adaptive_network import adaptive_voter
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    rows = []
    for phi in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99):
        for seed in (0, 1):
            t0 = time.time()
            fr, _ = adaptive_voter(seed, N=200, mean_deg=4, phi=phi, max_steps=400000, n_frames=80)
            A = fr[-1]["adjacency"] > 0
            op = fr[-1]["opinions"]
            n, lab = connected_components(csr_matrix(A), directed=False)
            sizes = np.bincount(lab)
            big = np.where(sizes >= 20)[0]
            purity = [float(max(np.mean(op[lab == c]), 1 - np.mean(op[lab == c]))) for c in big]
            rows.append({"phi": phi, "seed": seed,
                         "giant_fraction": float(sizes.max() / 200),
                         "n_components_ge_10pct": int(len(big)),
                         "large_component_opinion_purity": purity,
                         "final_minority_fraction": float(min(op.mean(), 1 - op.mean())),
                         "secs": round(time.time() - t0, 1)})
            r = rows[-1]
            print(f"P34 phi={phi:4.2f} s{seed} giant={r['giant_fraction']:.2f} "
                  f"big_comps={r['n_components_ge_10pct']} purity={['%.2f' % p for p in purity]} "
                  f"minority={r['final_minority_fraction']:.2f} ({r['secs']}s)", flush=True)
    return {"rows": rows}


def p35() -> Dict[str, Any]:
    from epc.models.entropy_crystal import hard_disk_crystallization

    def psi6_each(pos, L):
        rij = pos[:, None, :] - pos[None, :, :]
        rij -= L * np.round(rij / L)
        d = np.sqrt((rij ** 2).sum(-1)); np.fill_diagonal(d, 1e9)
        out = np.empty(len(pos))
        for i in range(len(pos)):
            nn = np.argsort(d[i])[:6]
            out[i] = abs(np.mean(np.exp(6j * np.arctan2(rij[i, nn, 1], rij[i, nn, 0]))))
        return out

    rows = []
    for eta_end in (0.20, 0.40, 0.55, 0.65, 0.70, 0.74):
        for seed in (0, 1):
            t0 = time.time()
            fr, _ = hard_disk_crystallization(seed, N=196, L=14.0, eta_end=eta_end,
                                              n_steps=2500, n_frames=80, noise0=0.16)
            p0 = psi6_each(fr[2]["positions"], 14.0)
            p1 = psi6_each(fr[-1]["positions"], 14.0)
            rows.append({"eta_end": eta_end, "seed": seed,
                         "crystalline_fraction_early": float(np.mean(p0 > 0.7)),
                         "crystalline_fraction_final": float(np.mean(p1 > 0.7)),
                         "mean_local_psi6_final": float(p1.mean()),
                         "secs": round(time.time() - t0, 1)})
            r = rows[-1]
            print(f"P35 eta_end={eta_end:.2f} s{seed} cryst_frac {r['crystalline_fraction_early']:.2f}->"
                  f"{r['crystalline_fraction_final']:.2f}  <|psi6|>={r['mean_local_psi6_final']:.3f} "
                  f"({r['secs']}s)", flush=True)
    return {"rows": rows}


def p37() -> Dict[str, Any]:
    from epc.models import resource_competition as rc
    rows = []

    def summarize(name, frames):
        A = np.array([f["abundances"] for f in frames])
        late = A[int(len(A) * 2 / 3):]
        peak = float(A.max()) or 1.0
        persist = late.mean(0) > 0.01 * peak
        cvs = [float(late[:, i].std() / late[:, i].mean()) for i in np.where(persist)[0]
               if late[:, i].mean() > 1e-9]
        rows.append({"regime": name, "n_persist": int(persist.sum()),
                     "late_cv_median": float(np.median(cvs)) if cvs else 0.0,
                     "late_cv_max": float(max(cvs)) if cvs else 0.0})
        r = rows[-1]
        print(f"P37 {name:40s} persist={r['n_persist']} late CV median={r['late_cv_median']:.2e} "
              f"max={r['late_cv_max']:.2e}", flush=True)

    for s in (0, 1):
        summarize(f"positive_HW_s{s}", rc.resource_competition(s)[0])
        summarize(f"competitive_exclusion_s{s}", rc.neg_competitive_exclusion(s)[0])
        summarize(f"stable_coexistence_s{s}", rc.neg_stable_coexistence(s)[0])
        summarize(f"single_resource_s{s}", rc.neg_single_resource(s)[0])
    # Added candidate regimes on the canonical Huisman-Weissing parameter set.
    K, C, S = rc._HW_K.copy(), rc._HW_C.copy(), rc._HW_S.copy()
    for s in (0, 1):
        rng = np.random.default_rng(s * 7919 + 1)
        N0 = 0.1 + 0.02 * rng.random(5)
        # (a) cross-consumption structure removed: every species consumes every
        #     resource at the same (mean) rate -> no intransitive consumption loop.
        Cu = np.full_like(C, float(C.mean()))
        summarize(f"HW_uniform_consumption_s{s}",
                  rc._integrate(N0.copy(), S.copy(), K, Cu, S, np.ones(5), np.full(5, 0.25), 0.25,
                                3000.0, 0.02, 80))
        # (b) high mortality m = 0.6 (growth cannot outpace loss for most species)
        summarize(f"HW_high_mortality_m0.6_s{s}",
                  rc._integrate(N0.copy(), S.copy(), K, C, S, np.ones(5), np.full(5, 0.6), 0.25,
                                3000.0, 0.02, 80))
        # (c) low supply (all S_j = 1): resources too scarce to sustain the loop
        summarize(f"HW_low_supply_S1_s{s}",
                  rc._integrate(N0.copy(), np.ones(5), K, C, np.ones(5), np.ones(5), np.full(5, 0.25), 0.25,
                                3000.0, 0.02, 80))
    return {"rows": rows}


def main(argv):
    which = argv[0] if argv else "all"
    res: Dict[str, Any] = {}
    for name, fn in (("p33", p33), ("p34", p34), ("p35", p35), ("p37", p37)):
        if which in (name, "all"):
            res[name] = fn()
    # merge into the existing file at write time (per-pattern runs may be separate)
    out: Dict[str, Any] = json.load(open(OUT)) if os.path.exists(OUT) else {}
    out.update(res)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=2)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
