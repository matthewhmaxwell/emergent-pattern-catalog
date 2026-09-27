"""EPC Ring-3 pre-registration v1 (osf.io/up7m9) — the registered decision rules, frozen BEFORE any data.

Reads prereg_runs/{smoke.json, <task>_<variant>[_rho<r>]_s<seed>.json, P4_crossplay.json} and prints a verdict per
prediction: HIT (the stated claim's fingerprint holds in >=2 of 3 seeds AND at the seed median), MISS (the competing
fingerprint holds in >=2 of 3 seeds AND at the seed median), else INCONCLUSIVE. Thresholds are copied from PREREG_v1.md
§2-§3; nothing here may be changed after the first training run (changes go in DEVIATIONS.md).
"""
import json, os, glob, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__)); RUNS = os.path.join(HERE, "prereg_runs")
COLLAPSE, DEAD, FAIR_FRAC = 0.50, 0.15, 0.90


def nu(s_full, s_abl, chance):
    den = s_full - chance
    return None if abs(den) < 1e-9 else (s_full - s_abl) / den


def learned(s_full, chance, ceil):
    return (s_full - chance) >= 0.20 * (ceil - chance)


def verdict(per_seed_c, med_c, per_seed_m, med_m):
    if sum(per_seed_c) >= 2 and med_c: return "HIT"
    if sum(per_seed_m) >= 2 and med_m: return "MISS"
    return "INCONCLUSIVE"


def load():
    rows = {}
    for f in glob.glob(os.path.join(RUNS, "P*_s*.json")):
        r = json.load(open(f)); rows[(r["task"], r["variant"], r.get("rho"), r["seed"])] = r
    return rows


def main():
    smoke = json.load(open(os.path.join(RUNS, "smoke.json"))); rows = load(); report = {}
    med = lambda xs: st.median(xs) if xs else float("nan")

    # ---- P1 (learned criterion exempt; C1 fingerprint defined by S_full + fair baseline, not by nu)
    ch = smoke["P1"]["chance_random_policy"]; ceil = smoke["P1"]["ceiling_oracle"]
    seeds = sorted(s for (t, v, _, s) in rows if t == "P1" and v == "full" and ("P1", "nochannel", None, s) in rows)
    F = {s: rows[("P1", "full", None, s)] for s in seeds}; Fr = {s: rows[("P1", "nochannel", None, s)]["S"] for s in seeds}
    c = [F[s]["S"] <= 0.60 and abs(Fr[s] - F[s]["S"]) <= 0.05 for s in seeds]
    n1 = {s: nu(F[s]["S"], F[s]["S_mute"], ch) for s in seeds}
    m = [F[s]["S"] >= 0.70 and (n1[s] or 0) >= COLLAPSE and Fr[s] <= 0.60 for s in seeds]
    mS, mF = med([F[s]["S"] for s in seeds]), med([Fr[s] for s in seeds]); mN = med([n1[s] or 0 for s in seeds])
    report["P1"] = {"verdict": verdict(c, mS <= 0.60 and abs(mF - mS) <= 0.05, m, mS >= 0.70 and mN >= COLLAPSE and mF <= 0.60),
                    "positive_control_ok": abs(ceil - 0.9375) <= 0.02, "ceiling": ceil, "chance": ch,
                    "per_seed": {s: {"S_full": F[s]["S"], "S_mute": F[s]["S_mute"], "S_fair_nochannel": Fr[s], "nu_mute": n1[s]} for s in seeds},
                    "note_if_HIT": "channel could solve the task (positive control) but PPO did not discover it within budget: "
                                   "a learnability result, not evidence that communication is useless under symmetry"}

    # ---- P2
    ch = smoke["P2"]["chance_random_policy"]; ceil = smoke["P2"]["ceiling_oracle"]
    seeds = sorted(s for (t, v, _, s) in rows if t == "P2" and v == "full" and ("P2", "nochannel", None, s) in rows)
    seeds = [s for s in seeds if learned(rows[("P2", "full", None, s)]["S"], ch, ceil)]
    P = {s: rows[("P2", "full", None, s)] for s in seeds}; fair = {s: rows[("P2", "nochannel", None, s)]["S"] for s in seeds}
    nb = {s: nu(P[s]["S"], P[s]["S_blind"], ch) for s in seeds}; nm = {s: nu(P[s]["S"], P[s]["S_mute"], ch) for s in seeds}
    c = [nb[s] >= COLLAPSE and abs(nm[s]) <= DEAD and fair[s] >= FAIR_FRAC * P[s]["S"] for s in seeds]
    m = [nm[s] >= 0.30 for s in seeds]
    mb, mm = med(list(nb.values())), med(list(nm.values())); mS, mF = med([P[s]["S"] for s in seeds]), med(list(fair.values()))
    report["P2"] = {"verdict": verdict(c, mb >= COLLAPSE and abs(mm) <= DEAD and mF >= FAIR_FRAC * mS, m, mm >= 0.30),
                    "learned_seeds": seeds, "chance": ch, "ceiling": ceil,
                    "per_seed": {s: {"S_full": P[s]["S"], "nu_blind": nb[s], "nu_mute": nm[s], "S_fair_nochannel": fair[s],
                                     "secondary_S_noobs_trained": rows.get(("P2", "noobs", None, s), {}).get("S")} for s in seeds}}

    # ---- P3 (registered chance = 0.50: random choice between goals)
    ch, ceil = smoke["P3"]["chance_registered"], smoke["P3"]["ceiling_oracle"]
    levels = [1.0, 0.75, 0.5, 0.25, 0.0]; nuL = {}
    for rho in levels:
        for (t, v, r, s), row in rows.items():
            if t == "P3" and v == "full" and r == rho and learned(row["S"], ch, ceil):
                nuL.setdefault(rho, {})[s] = nu(row["S"], row["S_memwipe"], ch)
    def fairmem(rho, s): return rows.get(("P3", "memoryless", rho, s), {}).get("S")
    s1 = nuL.get(1.0, {}); s0 = nuL.get(0.0, {})
    dead1 = [abs(n) <= DEAD and (fairmem(1.0, s) or 0) >= FAIR_FRAC * rows[("P3", "full", 1.0, s)]["S"] for s, n in s1.items()]
    coll0 = [n >= COLLAPSE for n in s0.values()]
    medL = {rho: med(list(nuL.get(rho, {}).values())) for rho in levels}
    xs = [1 - r for r in levels]; ys = [medL[r] for r in levels]
    def ranks(v): o = sorted(range(len(v)), key=lambda i: v[i]); rk = [0] * len(v); [rk.__setitem__(i, j) for j, i in enumerate(o)]; return rk
    rx, ry = ranks(xs), ranks(ys); n_ = len(xs)
    spearman = 1 - 6 * sum((a - b) ** 2 for a, b in zip(rx, ry)) / (n_ * (n_ ** 2 - 1)) if all(y == y for y in ys) else float("nan")
    med_dead1 = abs(medL[1.0]) <= DEAD and med([fairmem(1.0, s) or 0 for s in s1]) >= FAIR_FRAC * med([rows[("P3", "full", 1.0, s)]["S"] for s in s1])
    c2 = (sum(dead1) >= 2 and med_dead1) and (sum(coll0) >= 2 and medL[0.0] >= COLLAPSE) and spearman >= 0.8
    miss = sum(n >= 0.30 for n in s1.values()) >= 2 and medL[1.0] >= 0.30
    report["P3"] = {"verdict": "HIT" if c2 else ("MISS" if miss else "INCONCLUSIVE"), "chance": ch, "ceiling": ceil,
                    "median_nu_memwipe_by_rho": medL, "spearman_(1-rho)_vs_median_nu": spearman,
                    "fair_memoryless_S": {str(r): {s: fairmem(r, s) for s in (1, 2, 3)} for r in (1.0, 0.0)},
                    "per_seed_nu": {str(k): v for k, v in nuL.items()}}

    # ---- P4 (no separate competing fingerprint was registered: MISS = negation of the C1 criteria)
    ch = smoke["P4"]["chance_random_policy"]; ceil = smoke["P4"]["ceiling_oracle"]
    seeds = sorted(s for (t, v, _, s) in rows if t == "P4" and v == "full" and ("P4", "nochannel", None, s) in rows)
    seeds = [s for s in seeds if learned(rows[("P4", "full", None, s)]["S"], ch, ceil)]
    Q = {s: rows[("P4", "full", None, s)] for s in seeds}; fair = {s: rows[("P4", "nochannel", None, s)]["S"] for s in seeds}
    nm = {s: nu(Q[s]["S"], Q[s]["S_mute"], ch) for s in seeds}
    dead = [abs(nm[s]) <= DEAD and fair[s] >= FAIR_FRAC * Q[s]["S"] for s in seeds]
    xp = json.load(open(os.path.join(RUNS, "P4_crossplay.json"))) if os.path.exists(os.path.join(RUNS, "P4_crossplay.json")) else None
    cross_ok = xp is not None and xp["full"]["cross_play_mean"] <= xp["nochannel"]["cross_play_mean"] + 0.05
    med_dead = abs(med(list(nm.values()))) <= DEAD and med(list(fair.values())) >= FAIR_FRAC * med([Q[s]["S"] for s in seeds])
    hit = sum(dead) >= 2 and med_dead and cross_ok
    miss = (sum((nm[s] or 0) >= COLLAPSE for s in seeds) >= 2 and med(list(nm.values())) >= COLLAPSE) or \
           (xp is not None and not cross_ok)
    report["P4"] = {"verdict": "HIT" if hit else ("MISS" if miss else "INCONCLUSIVE"), "learned_seeds": seeds,
                    "chance": ch, "ceiling": ceil, "crossplay": xp,
                    "per_seed": {s: {"S_full": Q[s]["S"], "nu_mute": nm[s], "S_fair_nochannel": fair[s]} for s in seeds}}

    json.dump(report, open(os.path.join(RUNS, "VERDICTS.json"), "w"), indent=1, default=str)
    for k, v in report.items(): print(f"{k}: {v['verdict']}")
    print("full report: prereg_runs/VERDICTS.json")


if __name__ == "__main__":
    main()
