"""P23 minority game: literature reproduction v2 (spec: group_A.md, section P23).

T1 (primary): location alpha* of the sigma^2/N minimum, sign payoff, tie_break='random',
    S = 2, P = 64 and 128, seeds 0-19, burn-in 200P, measurement 1000P rounds.
    Pass: alpha*(P = 128) in [0.29, 0.39].
T1 sensitivity (reported): same grid with tie_break='first', seeds 0-9.
T2 (gating for the linear-payoff variant; reported in the plan): payoff='linear',
    tie_break='random', (P, N) = (512, 511), (512, 255) gating, (256, 255), (256, 127)
    finite-size check, seeds 0-9, burn-in 500P, measurement 2000P; sign-payoff companion.

Subcommands
  t1 --tie {random,first} --seeds A-B [--burn-mult 200] [--points P:N,...] --out FILE
  t2 --payoff {linear,sign} --seeds A-B --out FILE
  stationarity --tie {random,first} [--rerun]   pre-registered check; reruns failing points
                                                with burn-in 400P if --rerun
  analyze                                       writes analysis/outputs/repro_v2/P23.json

Units: A = 2 n1 - N; sigma^2/N = <A^2>/N (literature units, random agents give 1).
"""
from __future__ import annotations

import argparse
import json
import runpy
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
import mg_fast  # noqa: E402

from epc.models.minority_game import MinorityGame, MinorityGameParams  # noqa: E402

GRID = {
    64: [427, 355, 305, 267, 247, 229, 213, 207, 201, 193, 189, 183, 177, 169, 161, 143, 121, 101, 81, 65],
    128: [853, 711, 609, 533, 493, 457, 427, 413, 401, 387, 377, 365, 355, 337, 321, 285, 241, 203, 161, 129],
}
T2_POINTS = [(512, 511), (512, 255), (256, 255), (256, 127)]
T2_GATING = [(512, 511), (512, 255)]
ALPHA_C = 0.3374
T1_WINDOW = (0.29, 0.39)
STAT_BAND = 0.03          # stationarity check applies where |alpha - 0.3374| < 0.03
STAT_TOL = 0.05           # halves differing by more than 5% -> double the burn-in and rerun
N_BOOT = 2000
BOOT_SEED = 20261116
T1_PS = (64, 128)
PRIMARY_P = 128
BURN_MULT = 200           # T1 burn-in, in units of P rounds
MEAS_MULT = 1000          # T1 measurement, in units of P rounds
T2_BURN_MULT = 500
T2_MEAS_MULT = 2000
RAW = common.RAW_DIR


def raw_t1(tie: str, burn_mult: int = 200) -> Path:
    return RAW / f"P23_T1_{tie}_burn{burn_mult}P.jsonl"


def raw_t2(payoff: str) -> Path:
    return RAW / f"P23_T2_{payoff}.jsonl"


def _m(P: int) -> int:
    m = int(round(np.log2(P)))
    assert 2 ** m == P
    return m


def run_t1_point(P: int, N: int, seed: int, tie: str, burn_mult: int | None = None) -> dict:
    burn_mult = BURN_MULT if burn_mult is None else burn_mult
    meas_mult = MEAS_MULT
    t0 = time.time()
    model = MinorityGame(MinorityGameParams(n_agents=N, memory=_m(P), n_strategies=2, n_rounds=P,
                                            seed=seed, tie_break=tie, payoff="sign"))
    for _ in range(burn_mult):
        mg_fast.advance(model, P)
    a2 = np.empty(meas_mult)
    a1 = np.empty(meas_mult)
    for k in range(meas_mult):
        att, _, _ = mg_fast.advance(model, P)
        A = 2 * att - N
        a2[k] = float((A * A).sum())
        a1[k] = float(A.sum())
    T = meas_mult * P
    half = meas_mult // 2
    sigma2 = a2.sum() / T / N
    mean_A = a1.sum() / T
    return {
        "P": P, "N": N, "alpha": P / N, "seed": seed, "tie_break": tie, "payoff": "sign",
        "burn_rounds": burn_mult * P, "meas_rounds": T,
        "sigma2_over_N": sigma2,
        "sigma2_over_N_first_half": a2[:half].sum() / (half * P) / N,
        "sigma2_over_N_second_half": a2[half:].sum() / ((meas_mult - half) * P) / N,
        "sigma2_over_N_deciles": (a2.reshape(10, -1).sum(1) / (T / 10) / N).tolist(),
        "mean_A": mean_A,
        "var_n1_over_N_detector_units": (sigma2 - mean_A ** 2 / N) / 4.0,
        "runtime_s": time.time() - t0,
    }


def run_t2_point(P: int, N: int, seed: int, payoff: str) -> dict:
    burn_mult, meas_mult = T2_BURN_MULT, T2_MEAS_MULT
    t0 = time.time()
    model = MinorityGame(MinorityGameParams(n_agents=N, memory=_m(P), n_strategies=2, n_rounds=P,
                                            seed=seed, tie_break="random", payoff=payoff))
    for _ in range(burn_mult):
        mg_fast.advance(model, P)
    d0 = model.scores[:, 0] - model.scores[:, 1]
    sign0 = np.sign(d0)
    same_all = sign0 != 0
    same_first = same_all.copy()
    d_mid = None
    sign_mid = None
    same_second = None
    count = np.zeros(P)
    sumA = np.zeros(P)
    sumA2 = 0.0
    mid = meas_mult // 2
    for k in range(1, meas_mult + 1):
        att, _, mu = mg_fast.advance(model, P)
        A = 2 * att - N
        count += np.bincount(mu, minlength=P)
        sumA += np.bincount(mu, weights=A.astype(np.float64), minlength=P)
        sumA2 += float((A * A).sum())
        d = model.scores[:, 0] - model.scores[:, 1]
        sg = np.sign(d)
        same_all &= sg == sign0
        if k <= mid:
            same_first &= sg == sign0
        if k == mid:
            d_mid = d.copy()
            sign_mid = sg.copy()
            same_second = sign_mid != 0
        if k > mid:
            same_second &= sg == sign_mid
    d_end = d
    T = meas_mult * P
    f = count / T
    cond = np.divide(sumA, count, out=np.zeros(P), where=count > 0)
    return {
        "P": P, "N": N, "alpha": P / N, "seed": seed, "tie_break": "random", "payoff": payoff,
        "burn_rounds": burn_mult * P, "meas_rounds": T, "n_snapshots": meas_mult + 1,
        "phi_hat": float(np.mean(same_all & (np.abs(d_end) > np.abs(d0)))),
        "phi_hat_first_half": float(np.mean(same_first & (np.abs(d_mid) > np.abs(d0)))),
        "phi_hat_second_half": float(np.mean(same_second & (np.abs(d_end) > np.abs(d_mid)))),
        "H_hat_over_N": float(np.sum(f * cond ** 2) / N),
        "sigma2_over_N": sumA2 / T / N,
        "n_mu_visited": int(np.count_nonzero(count)),
        "runtime_s": time.time() - t0,
    }


def _parse_seeds(s: str) -> list[int]:
    a, b = s.split("-")
    return list(range(int(a), int(b) + 1))


def _parse_points(s: str | None) -> list[tuple[int, int]]:
    if not s:
        return [(P, N) for P in T1_PS for N in GRID[P]]
    return [tuple(int(x) for x in p.split(":")) for p in s.split(",")]


def _done(path: Path, keys: tuple[str, ...]) -> set:
    return {tuple(r[k] for k in keys) for r in common.read_jsonl(path)}


def cmd_t1(tie: str, seeds: list[int], burn_mult: int, points: list[tuple[int, int]], out: Path) -> None:
    git = common.git_info()
    done = _done(out, ("P", "N", "seed"))
    for (P, N) in points:
        for s in seeds:
            if (P, N, s) in done:
                continue
            rec = run_t1_point(P, N, s, tie, burn_mult=burn_mult)
            rec["git_commit"] = git["commit"]
            rec["tracked_tree_clean"] = git["tracked_tree_clean"]
            rec["timestamp"] = common.now_iso()
            common.append_jsonl(out, rec)
            print(f"T1 {tie} burn{burn_mult}P P={P} N={N} seed={s} {rec['runtime_s']:.1f}s", flush=True)


def cmd_t2(payoff: str, seeds: list[int], out: Path) -> None:
    git = common.git_info()
    done = _done(out, ("P", "N", "seed"))
    for (P, N) in T2_POINTS:
        for s in seeds:
            if (P, N, s) in done:
                continue
            rec = run_t2_point(P, N, s, payoff)
            rec["git_commit"] = git["commit"]
            rec["tracked_tree_clean"] = git["tracked_tree_clean"]
            rec["timestamp"] = common.now_iso()
            common.append_jsonl(out, rec)
            print(f"T2 {payoff} P={P} N={N} seed={s} {rec['runtime_s']:.1f}s", flush=True)


# ----------------------------------------------------------------------------- analysis
def stationarity(records: list[dict]) -> list[dict]:
    """Pre-registered check: at |alpha - 0.3374| < 0.03, seed-mean of the first half vs
    the second half of the measurement window; fails if |h2 - h1| / min(h1, h2) > 0.05."""
    out = []
    pts = sorted({(r["P"], r["N"]) for r in records})
    for (P, N) in pts:
        alpha = P / N
        if abs(alpha - ALPHA_C) >= STAT_BAND:
            continue
        rs = [r for r in records if r["P"] == P and r["N"] == N]
        h1 = float(np.mean([r["sigma2_over_N_first_half"] for r in rs]))
        h2 = float(np.mean([r["sigma2_over_N_second_half"] for r in rs]))
        rel = abs(h2 - h1) / min(h1, h2)
        out.append({"P": P, "N": N, "alpha": alpha, "n_seeds": len(rs), "burn_rounds": rs[0]["burn_rounds"],
                    "seedmean_first_half": h1, "seedmean_second_half": h2, "rel_diff": rel,
                    "fails_5pct": bool(rel > STAT_TOL)})
    return out


def cmd_stationarity(tie: str, rerun: bool, seeds: list[int]) -> None:
    recs = common.read_jsonl(raw_t1(tie, BURN_MULT))
    chk = stationarity(recs)
    common.write_json(RAW / f"P23_T1_{tie}_stationarity_burn{BURN_MULT}P.json", chk)
    fails = [(c["P"], c["N"]) for c in chk if c["fails_5pct"]]
    for c in chk:
        print(f"stationarity {tie} P={c['P']} N={c['N']} alpha={c['alpha']:.4f} rel={c['rel_diff']:.4f} "
              f"{'FAIL' if c['fails_5pct'] else 'ok'}")
    if fails and rerun:
        print(f"rerunning {len(fails)} points with burn-in {2 * BURN_MULT}P", flush=True)
        cmd_t1(tie, seeds, 2 * BURN_MULT, fails, raw_t1(tie, 2 * BURN_MULT))
        chk2 = stationarity(common.read_jsonl(raw_t1(tie, 2 * BURN_MULT)))
        common.write_json(RAW / f"P23_T1_{tie}_stationarity_burn{2 * BURN_MULT}P.json", chk2)
        for c in chk2:
            print(f"stationarity(after rerun) {tie} P={c['P']} N={c['N']} rel={c['rel_diff']:.4f} "
                  f"{'FAIL' if c['fails_5pct'] else 'ok'}")


def _t1_records(tie: str) -> tuple[list[dict], list[tuple[int, int]]]:
    """Burn-in-200P records, with points rerun at 400P (stationarity) replaced."""
    base = common.read_jsonl(raw_t1(tie, BURN_MULT))
    rerun = common.read_jsonl(raw_t1(tie, 2 * BURN_MULT))
    rerun_pts = sorted({(r["P"], r["N"]) for r in rerun})
    recs = [r for r in base if (r["P"], r["N"]) not in rerun_pts] + rerun
    return recs, rerun_pts


def alpha_star(records: list[dict], P: int, seeds: list[int] | None = None, boot: bool = True) -> dict:
    Ns = sorted(GRID[P], key=lambda n: P / n)          # ascending alpha
    alphas = np.array([P / n for n in Ns])
    seeds_present = sorted({r["seed"] for r in records if r["P"] == P})
    if seeds is not None:
        seeds_present = [s for s in seeds_present if s in seeds]
    lut = {(r["N"], r["seed"]): r["sigma2_over_N"] for r in records if r["P"] == P}
    S = np.array([[lut[(n, s)] for s in seeds_present] for n in Ns])   # (n_alpha, n_seeds)
    curve = S.mean(1)
    se = S.std(1, ddof=1) / np.sqrt(S.shape[1])
    i = int(np.argmin(curve))
    res = {"P": P, "n_seeds": len(seeds_present), "seeds": seeds_present,
           "alpha_grid": alphas.tolist(), "N_grid": Ns,
           "sigma2_over_N_seed_mean": curve.tolist(), "sigma2_over_N_seed_se": se.tolist(),
           "alpha_star": float(alphas[i]), "N_at_alpha_star": int(Ns[i]),
           "min_sigma2_over_N": float(curve[i]), "min_sigma2_over_N_se": float(se[i])}
    if boot:
        rng = np.random.default_rng(BOOT_SEED)
        bs = np.empty(N_BOOT)
        ns = S.shape[1]
        for b in range(N_BOOT):
            idx = rng.integers(0, ns, ns)
            bs[b] = alphas[int(np.argmin(S[:, idx].mean(1)))]
        res["bootstrap"] = {"n_resamples": N_BOOT, "rng_seed": BOOT_SEED,
                            "resampling": "seed indices resampled with replacement, jointly across grid points",
                            "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                            "frac_in_pass_window": float(np.mean((bs >= T1_WINDOW[0]) & (bs <= T1_WINDOW[1]))),
                            "alpha_star_counts": {f"{a:.5f}": int(np.sum(bs == a)) for a in np.unique(bs)}}
    return res


def cmd_analyze() -> None:
    rs = runpy.run_path(str(common.SUPPORT_DIR / "rs_minority_game.py"))
    zc = rs["zc"]
    alpha_c_rs = float(rs["alpha"](zc))
    Qc = float(rs["Q"](zc))
    rs_min = (1 - Qc) / 2
    theory = {"script": "docs/validation_rebuild/repro_v2/group_A_support/rs_minority_game.py",
              "alpha_c": alpha_c_rs, "Q_c": Qc, "sigma2_over_N_at_alpha_c_RS": rs_min,
              "T2_points": {f"{P}:{N}": {"alpha": P / N, **{k: float(v) for k, v in rs["rs"](P / N).items()}}
                            for (P, N) in T2_POINTS}}

    out: dict = {"pattern": "P23", "title": "Minority game, anti-coordination",
                 "spec": "docs/validation_rebuild/repro_v2/group_A.md#p23", "theory": theory}

    # ---- T1 primary (random tie-break)
    recs_r, rerun_r = _t1_records("random")
    stat200 = RAW / f"P23_T1_random_stationarity_burn{BURN_MULT}P.json"
    stat400 = RAW / f"P23_T1_random_stationarity_burn{2 * BURN_MULT}P.json"
    stat_r = {"burn200P": json.loads(stat200.read_text()) if stat200.exists() else None,
              "burn400P_rerun": json.loads(stat400.read_text()) if stat400.exists() else None,
              "points_rerun_with_burn400P": rerun_r}
    t1 = {P: alpha_star(recs_r, P) for P in T1_PS}
    PO = [p for p in T1_PS if p != PRIMARY_P][0]
    a128 = t1[PRIMARY_P]["alpha_star"]
    t1_pass = T1_WINDOW[0] <= a128 <= T1_WINDOW[1]
    out["T1"] = {
        "rule": "alpha*(P=128) in [0.29, 0.39], i.e. |alpha* - 0.3374| <= 0.05 (primary)",
        "settings": {"payoff": "sign", "tie_break": "random", "S": 2, "seeds": "0-19",
                     "burn_in_rounds": "200P (400P at points rerun by the stationarity check)",
                     "measurement_rounds": "1000P", "chunking": "P rounds per call"},
        f"P{PRIMARY_P}": t1[PRIMARY_P], f"P{PO}_reported": t1[PO],
        f"alpha_star_P{PRIMARY_P}": a128, "pass": bool(t1_pass),
        "min_sigma2_vs_RS": {f"P{PRIMARY_P}_min": t1[PRIMARY_P]["min_sigma2_over_N"],
                             f"P{PO}_min": t1[PO]["min_sigma2_over_N"],
                             "RS_value_at_alpha_c": rs_min,
                             "note": "reported only; expected higher at finite N"},
        "stationarity_check": stat_r,
        "total_runtime_s": float(sum(r["runtime_s"] for r in recs_r)),
        "n_runs": len(recs_r),
    }

    # ---- T1 sensitivity: tie_break='first', seeds 0-9
    recs_f, rerun_f = _t1_records("first")
    if recs_f:
        t1f = {P: alpha_star(recs_f, P) for P in T1_PS}
        r_0_9 = {P: alpha_star(recs_r, P, seeds=list(range(10)), boot=False) for P in T1_PS}
        diffs = {P: abs(t1f[P]["alpha_star"] - t1[P]["alpha_star"]) for P in T1_PS}
        stat200f = RAW / f"P23_T1_first_stationarity_burn{BURN_MULT}P.json"
        stat400f = RAW / f"P23_T1_first_stationarity_burn{2 * BURN_MULT}P.json"
        out["T1_tie_first_sensitivity"] = {
            "reported_only": True,
            "flag_rule": "flag if |alpha*(first) - alpha*(random)| > 0.05 at P = 128 or P = 64",
            f"P{PRIMARY_P}": t1f[PRIMARY_P], f"P{PO}": t1f[PO],
            "random_seeds_0_9_for_reference": {P: {"alpha_star": r_0_9[P]["alpha_star"],
                                                   "min_sigma2_over_N": r_0_9[P]["min_sigma2_over_N"]}
                                               for P in T1_PS},
            "abs_diff_vs_primary": {str(P): diffs[P] for P in T1_PS},
            "flag_raised": bool(any(d > 0.05 + 1e-12 for d in diffs.values())),
            "stationarity_check": {
                "burn200P": json.loads(stat200f.read_text()) if stat200f.exists() else None,
                "burn400P_rerun": json.loads(stat400f.read_text()) if stat400f.exists() else None,
                "points_rerun_with_burn400P": rerun_f},
            "total_runtime_s": float(sum(r["runtime_s"] for r in recs_f)),
            "n_runs": len(recs_f),
        }

    # ---- T2
    t2 = {}
    for payoff in ("linear", "sign"):
        recs = common.read_jsonl(raw_t2(payoff))
        if not recs:
            continue
        pts = {}
        for (P, N) in T2_POINTS:
            rr = [r for r in recs if r["P"] == P and r["N"] == N]
            if not rr:
                continue
            th = theory["T2_points"][f"{P}:{N}"]
            agg = {}
            for key in ("phi_hat", "phi_hat_first_half", "phi_hat_second_half", "H_hat_over_N", "sigma2_over_N"):
                v = np.array([r[key] for r in rr])
                agg[key] = {"seed_mean": float(v.mean()), "seed_sd": float(v.std(ddof=1)),
                            "se": float(v.std(ddof=1) / np.sqrt(len(v))), "per_seed": v.tolist()}
            phi, H = agg["phi_hat"]["seed_mean"], agg["H_hat_over_N"]["seed_mean"]
            pts[f"{P}:{N}"] = {
                "P": P, "N": N, "alpha": P / N, "n_seeds": len(rr), "gating": (P, N) in T2_GATING and payoff == "linear",
                "RS": {"phi": th["phi"], "H_over_N": th["H"], "Q": th["Q"], "sigma2_over_N_approx": th["sigma2_approx"]},
                "measured": agg,
                "phi_abs_diff": abs(phi - th["phi"]), "phi_within_0.03": bool(abs(phi - th["phi"]) <= 0.03),
                "H_rel_diff": abs(H - th["H"]) / th["H"], "H_within_10pct": bool(abs(H - th["H"]) <= 0.10 * th["H"]),
                "runtime_s": float(sum(r["runtime_s"] for r in rr)),
            }
        t2[payoff] = pts
    if "linear" in t2:
        g = [t2["linear"][f"{P}:{N}"] for (P, N) in T2_GATING if f"{P}:{N}" in t2["linear"]]
        t2_pass = len(g) == 2 and all(p["phi_within_0.03"] and p["H_within_10pct"] for p in g)
        out["T2"] = {
            "rule": "payoff='linear', tie_break='random': at P = 512 for alpha ~ 1.002 and ~ 2.008, "
                    "|phi_hat - phi_RS| <= 0.03 and |H_hat/N - H_RS| <= 0.10 H_RS "
                    "(gating for the linear-payoff variant; reported, not gating, in the plan)",
            "estimators": {"phi_hat": "fraction of agents whose Delta_i = U_i0 - U_i1 has the same nonzero sign "
                                      "at all 2001 snapshots (every P rounds) and |Delta_i(end)| > |Delta_i(start)|; "
                                      "seed mean",
                           "H_hat_over_N": "(1/N) sum_mu f_mu <A|mu>^2 over the measurement window; seed mean",
                           "half_windows": "same rule on snapshots 0-1000 and 1000-2000"},
            "linear": t2["linear"], "sign_companion_reported_only": t2.get("sign"),
            "pass": bool(t2_pass),
        }

    out["deviations"] = common.load_deviations("P23")
    out["seeds"] = {"T1": "0-19", "T1_tie_first": "0-9", "T2": "0-9"}
    out["code_commits_of_runs"] = sorted({r["git_commit"] for r in recs_r + recs_f +
                                          common.read_jsonl(raw_t2("linear")) + common.read_jsonl(raw_t2("sign"))})
    out["all_runs_from_clean_tree"] = all(r["tracked_tree_clean"] for r in recs_r + recs_f +
                                          common.read_jsonl(raw_t2("linear")) + common.read_jsonl(raw_t2("sign")))
    out["analysis_git"] = common.git_info()
    out["bit_identity_checks"] = "analysis/outputs/repro_v2/bit_identity_checks.json"
    out["total_compute_runtime_s"] = float(sum(r["runtime_s"] for r in recs_r + recs_f +
                                               common.read_jsonl(raw_t2("linear")) + common.read_jsonl(raw_t2("sign"))))
    out["primary_rule"] = "T1"
    out["primary_pass"] = out["T1"]["pass"]
    out["timestamp"] = common.now_iso()
    common.write_json(common.OUT_DIR / "P23.json", out)
    print(f"T1 alpha*(128) = {a128:.4f} -> {'PASS' if t1_pass else 'FAIL'}")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("t1"); a.add_argument("--tie", required=True); a.add_argument("--seeds", required=True)
    a.add_argument("--burn-mult", type=int, default=None); a.add_argument("--points")
    b = sub.add_parser("t2"); b.add_argument("--payoff", required=True); b.add_argument("--seeds", required=True)
    c = sub.add_parser("stationarity"); c.add_argument("--tie", required=True); c.add_argument("--seeds", required=True)
    c.add_argument("--rerun", action="store_true")
    sub.add_parser("analyze")
    args = ap.parse_args()
    if args.cmd == "t1":
        bm = args.burn_mult or BURN_MULT
        cmd_t1(args.tie, _parse_seeds(args.seeds), bm, _parse_points(args.points), raw_t1(args.tie, bm))
    elif args.cmd == "t2":
        cmd_t2(args.payoff, _parse_seeds(args.seeds), raw_t2(args.payoff))
    elif args.cmd == "stationarity":
        cmd_stationarity(args.tie, args.rerun, _parse_seeds(args.seeds))
    elif args.cmd == "analyze":
        cmd_analyze()


if __name__ == "__main__":
    main()
